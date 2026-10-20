"""Example figure for the liquid-cloud calibration (AMT paper), improving the operational
dashboard diagnostic. One representative day (Palaiseau CL31, 2026-06-12). Panels:
  (a) pcolor of log10(beta') time-height, ALL cloud base (black) + calibration profiles (red);
  (b) one representative profile at its acquisition time, integration window;
  (c) water-vapour two-way transmission profile T^2_WV(z) (the 910 nm correction), median + daily
      min-max spread;
  (d) apparent lidar ratio S time series, valid (green) vs rejected (black);
  (e) O'Connor (2004) Fig.5-style: apparent S vs cloud height (2D density), with the theoretical
      eta(z)*S curve from the Cloudnet-measured droplet effective radius THAT DAY;
  (f) calibration-coefficient C histogram (C = 1/C_L, with C_L the paper's lidar constant).

Everything shares one averaged grid so the axes align: pcolor and the S time series share the
hours-UTC x-axis; the profile, the WV panel and the pcolor share the range-AGL y-axis.
"""
from __future__ import annotations

import copy
import sys
import warnings
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
warnings.filterwarnings("ignore")

from calibration.cloud.calibration import (  # noqa: E402
    CloudCalConfig, set_defaults, compute_wv_transmission, apply_instrument_filters,
    calculate_lidar_ratio, apply_cloud_filters, apply_temporal_consistency_filter)
from tests._ceilo_reader import read_ceilometer_data, average_ceilo_data  # noqa: E402
from validation.multiple_scattering_eta import eta_for_cloud_base, eta_table, CBH_GRID_KM  # noqa: E402

L1 = "D:/E-PROFILE_L1_2026/0-250-1001-07151/2026/06/L1_0-250-1001-07151_A20260612.nc"
CN = "D:/Cloudnet_Droplet_radius/20260612_palaiseau_der.nc"   # moved from A: (2026-10)
CAMS = "D:/CAMS"
OUT = Path("C:/DATA/Projects/202606_E-PROFILE_calibration/old/ms_fov")   # moved from .../ms_fov (2026-10)
REPO_FIGS = Path(__file__).resolve().parents[1] / "doc" / "reports" / "figs_multiple_scattering"
DAY = "2026-06-12"
LAT, LON = 48.72, 2.21
S_LIQ = 18.8
FOV_CL31, DIV_CL31, WL_CL31 = 0.83e-3, 0.57e-3, 910e-9


def cloudnet_re_profiles(times_epoch, cbh_m, salt):
    """Cloudnet-measured droplet effective radius [um] at each calibration profile (median over
    cloud base .. +300 m), plus the day median."""
    from netCDF4 import Dataset
    with Dataset(CN) as nc:
        cnt = np.asarray(nc.variables["time"][:], float)     # hours since midnight
        cnh = np.asarray(nc.variables["height"][:], float)   # MSL
        der = np.asarray(nc.variables["der"][:], float)      # (time,height) m
    der[~np.isfinite(der) | (der <= 0) | (der > 40e-6)] = np.nan
    day0 = np.datetime64(f"{DAY}T00:00:00").astype("datetime64[s]").astype(np.int64)
    re_prof = []
    for tsec, cb in zip(times_epoch, cbh_m):
        hh = (tsec - day0) / 3600.0
        ti = int(np.argmin(np.abs(cnt - hh)))
        if abs(cnt[ti] - hh) > 0.2:
            re_prof.append(np.nan); continue
        cloud_msl = cb + salt
        hm = (cnh >= cloud_msl - 30) & (cnh <= cloud_msl + 300)
        col = der[ti, hm]
        re_prof.append(np.nanmedian(col) * 1e6 if np.isfinite(col).any() else np.nan)
    re_prof = np.array(re_prof)
    return re_prof, float(np.nanmedian(re_prof[np.isfinite(re_prof)]))


def main():
    REPO_FIGS.mkdir(parents=True, exist_ok=True)
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.gridspec import GridSpec

    cfg = set_defaults(CloudCalConfig(
        nc_file=L1, instrument="CL31", apply_wv_correction=True,
        apply_transmission_correction=True, aerosol_lidar_ratio=50.0,
        cams_folder=CAMS, abs_cs_lookup_table="",
        station_latitude=LAT, station_longitude=LON,
        average_time_s=30.0, average_range_m=10.0))
    data0, status = read_ceilometer_data(cfg.nc_file, cfg)
    salt = float(getattr(data0, "station_altitude", 157.0) or 157.0)

    # ---- replicate the calibration pipeline on ONE averaged grid, so every panel shares it ----
    data = average_ceilo_data(copy.deepcopy(data0), cfg)
    trans2 = np.asarray(compute_wv_transmission(data, cfg), float)   # (n_range, n_time) two-way T^2
    data.beta = data.beta / trans2                                   # WV-corrected calibrated beta
    beta_wv = np.asarray(data.beta, float).copy()                   # for the pcolor + profile
    rng = np.asarray(data.range, float)

    # apparent S per profile WITHOUT the eta correction (flat), so the O'Connor height dependence is
    # preserved in the data; the cloud + temporal filters then flag which profiles are valid.
    beta_filt, _ = apply_instrument_filters(beta_wv.copy(), data, cfg)
    S_app_all, _ = calculate_lidar_ratio(beta_filt, data.range, cfg)      # (n_time,) sr
    S_cloud, _ = apply_cloud_filters(S_app_all.copy(), beta_filt, data, cfg)
    S_cons, _ = apply_temporal_consistency_filter(S_cloud, cfg)
    valid = np.isfinite(S_cons)                                          # passed ALL five tests

    cbh_all = np.asarray(data.cbh, float)
    cloud_present = np.isfinite(cbh_all) & (cbh_all >= 300) & (cbh_all <= 2600) & np.isfinite(S_app_all)
    tsec_all = data.time.astype("datetime64[s]").astype(np.int64)
    day0 = data.time.astype("datetime64[D]")[0].astype("datetime64[s]").astype(np.int64)
    hours = (tsec_all - day0) / 3600.0

    # --- Cloudnet-measured droplet radius at the valid calibration scenes ---
    re_prof, re_day = cloudnet_re_profiles(tsec_all[valid], cbh_all[valid], salt)
    print(f"Cloudnet measured eff. radius (day median at cal scenes): {re_day:.1f} um "
          f"(=> {2*re_day:.0f} um diameter)")

    # --- PVC multiple-scattering eta at the MEASURED radius, per-profile O'Connor coefficient ---
    et = np.column_stack([CBH_GRID_KM, eta_table("CL31", re_day, alpha=10e-3)])
    eta_valid = np.interp(cbh_all[valid] / 1000.0, et[:, 0], et[:, 1])
    C_prof = (S_app_all[valid] / eta_valid) / S_LIQ
    C_med, C_std = float(np.median(C_prof)), float(np.std(C_prof, ddof=1))
    CL_med = 1.0 / C_med                       # paper convention: C = 1 / C_L
    # apparent S on the physical scale (single day-constant): ~= eta(z)*S_liq, rising to 18.8 sr
    # where multiple scattering is weakest (low clouds).
    S_phys = S_app_all / C_med

    # theoretical eta(z)*S for the measured radius (+/-2 um band for context)
    zc = np.linspace(0.3, 2.6, 80)
    def theo(re):
        th = WL_CL31 / (np.pi * re * 1e-6)
        return np.array([eta_for_cloud_base(z*1000, FOV_CL31, DIV_CL31, th, alpha=10e-3, dr=5.0)
                         for z in zc]) * S_LIQ
    theo_meas, theo_lo, theo_hi = theo(re_day), theo(re_day + 2), theo(re_day - 2)

    # representative profile: the valid profile whose cloud base is closest to the daily median
    vi = np.where(valid)[0]
    rep = vi[int(np.argmin(np.abs(cbh_all[vi] - np.median(cbh_all[vi]))))]
    rep_t = data.time[rep].astype("datetime64[s]")
    rep_hhmm = str(rep_t).split("T")[1][:5]

    # ---------------- figure ----------------
    from mpl_toolkits.axes_grid1 import make_axes_locatable
    fig = plt.figure(figsize=(17, 9.4))
    gs = GridSpec(2, 3, figure=fig, width_ratios=[2.4, 1.0, 1.3], height_ratios=[1, 1],
                  hspace=0.34, wspace=0.32)
    axPC = fig.add_subplot(gs[0, 0])                    # (a) pcolor
    axPr = fig.add_subplot(gs[0, 1], sharey=axPC)       # (b) profile        (shares altitude y)
    axWV = fig.add_subplot(gs[0, 2], sharey=axPC)       # (c) WV correction  (shares altitude y)
    axTS = fig.add_subplot(gs[1, 0], sharex=axPC)       # (d) S time series  (shares hours x)
    axOC = fig.add_subplot(gs[1, 1])                    # (e) apparent S vs height
    axCH = fig.add_subplot(gs[1, 2])                    # (f) C histogram

    # (a) pcolor + all cloud base (black) + calibration profiles (red)
    with np.errstate(all="ignore"):
        lb = np.log10(np.where(beta_wv > 0, beta_wv, np.nan))
    pm = axPC.pcolormesh(hours, rng / 1000.0, lb, cmap="viridis", vmin=-6.7, vmax=-4.5, shading="auto")
    axPC.scatter(hours[cloud_present], cbh_all[cloud_present] / 1000.0, s=9, c="k", marker=".",
                 edgecolor="none", alpha=0.55, label="cloud base", zorder=3)
    axPC.scatter(hours[valid], cbh_all[valid] / 1000.0, s=20, c="red", marker=".",
                 edgecolor="none", alpha=0.95, label="calibration profiles", zorder=4)
    axPC.axvline(hours[rep], color="red", ls=":", lw=1.0, alpha=0.7)
    axPC.set_ylim(0, 3.3); axPC.set_ylabel("range AGL [km]"); axPC.set_xlabel("hours (UTC)")
    axPC.set_title("(a) attenuated backscatter $\\beta'$ (log$_{10}$) and cloud base")
    # colorbar via a divider so it does NOT shrink the pcolor relative to the time series below
    # (an equal-width, invisible divider is added to the time series -> x-axes stay aligned).
    plt.colorbar(pm, cax=make_axes_locatable(axPC).append_axes("right", size="2.5%", pad=0.08),
                 label="log$_{10}\\beta'$")
    axPC.legend(loc="upper right", fontsize=8, framealpha=0.9)

    # (b) representative profile at its acquisition time
    axPr.plot(beta_wv[:, rep], rng / 1000.0, "-", color="C0", lw=1.3)
    axPr.axhspan(0.1, 2.4, color="gold", alpha=0.12, label="integration 0.1-2.4 km")
    axPr.axhline(cbh_all[rep] / 1000.0, color="firebrick", ls="--", lw=1.2, label="cloud base")
    axPr.set_xlabel("$\\beta'$")
    axPr.set_title(f"(b) profile at {rep_hhmm} UTC")
    axPr.ticklabel_format(axis="x", style="sci", scilimits=(0, 0))
    axPr.legend(fontsize=7.5, loc="upper right"); axPr.grid(alpha=0.3)

    # (c) water-vapour two-way transmission profile (the 910 nm correction), daily median + spread
    twv = np.nanmedian(trans2, axis=1)
    twv_lo = np.nanmin(trans2, axis=1)
    twv_hi = np.nanmax(trans2, axis=1)
    win = (rng >= 100) & (rng <= 2400)
    t_min, t_max = float(np.nanmin(trans2[win, :])), float(np.nanmax(trans2[win, :]))
    axWV.fill_betweenx(rng / 1000.0, twv_lo, twv_hi, color="C4", alpha=0.20, label="daily min-max")
    axWV.plot(twv, rng / 1000.0, "-", color="C4", lw=1.7, label="daily median")
    axWV.axhline(cbh_all[rep] / 1000.0, color="firebrick", ls="--", lw=1.0, alpha=0.7)
    axWV.set_xlabel("$T^2_{WV}$ (two-way)")
    axWV.set_title("(c) water-vapour correction")
    axWV.set_xlim(max(0.4, t_min - 0.04), 1.01); axWV.grid(alpha=0.3)
    axWV.text(0.04, 0.045, "$\\beta' \\leftarrow \\beta'/T^2_{WV}$  (CAMS+HITRAN)\n"
              f"0.1-2.4 km: $T^2$ = {t_min:.2f} - {t_max:.2f}",
              transform=axWV.transAxes, fontsize=7.6, va="bottom",
              bbox=dict(boxstyle="round", fc="white", ec="0.7", alpha=0.85))
    axWV.legend(fontsize=7.2, loc="upper left")

    # (d) apparent lidar ratio time series: valid (green) vs rejected (black) -- shares x with (a)
    inval = cloud_present & ~valid
    axTS.scatter(hours[inval], S_phys[inval], s=9, c="k", alpha=0.45, label="rejected profiles")
    axTS.scatter(hours[valid], S_phys[valid], s=12, c="C2", label="valid (calibration) profiles")
    axTS.axhline(S_LIQ, color="k", ls="--", lw=1.2, label="single scatter 18.8 sr")
    axTS.axhline(np.median(S_phys[valid]), color="C2", lw=1.2, alpha=0.7,
                 label=f"median valid {np.median(S_phys[valid]):.1f} sr")
    axTS.axvline(hours[rep], color="red", ls=":", lw=1.0, alpha=0.7)
    axTS.set_xlabel("hours (UTC)"); axTS.set_ylabel("apparent S [sr]"); axTS.set_ylim(8, 26)
    axTS.set_title("(d) apparent lidar ratio time series"); axTS.grid(alpha=0.3)
    axTS.legend(fontsize=7.5, loc="upper right", ncol=2)
    # invisible divider of the same width as the pcolor colorbar -> keeps (a) and (d) x-aligned
    make_axes_locatable(axTS).append_axes("right", size="2.5%", pad=0.08).axis("off")

    # (e) O'Connor Fig.5-style apparent S vs cloud height + theoretical eta*S(measured radius)
    hb = axOC.hist2d(S_phys[valid], cbh_all[valid] / 1000.0,
                     bins=[np.linspace(8, 24, 45), np.linspace(0.3, 2.6, 40)], cmap="turbo", cmin=1)
    axOC.plot(theo_meas, zc, "k-", lw=2.4, label=f"$\\eta$S theory, $r_e$={re_day:.1f} $\\mu$m (Cloudnet)")
    axOC.fill_betweenx(zc, theo_lo, theo_hi, color="k", alpha=0.12, label="$r_e\\pm2$ $\\mu$m")
    axOC.axvline(S_LIQ, color="0.25", ls="--", lw=1.4)
    axOC.text(S_LIQ + 0.2, 2.5, "single\nscatter\n18.8 sr", color="0.25", fontsize=7.5, va="top")
    axOC.set_xlim(8, 24); axOC.set_ylim(0, 3.3)
    axOC.set_xlabel("apparent lidar ratio S [sr]"); axOC.set_ylabel("cloud base height [km]")
    axOC.set_title("(e) apparent S vs height + $\\eta$S theory")
    axOC.legend(fontsize=7.5, loc="upper right")
    plt.colorbar(hb[3], cax=make_axes_locatable(axOC).append_axes("right", size="3%", pad=0.08),
                 label="count")

    # (f) C histogram (C = 1/C_L)
    axCH.hist(C_prof, bins=np.linspace(C_med - 4*C_std, C_med + 4*C_std, 40), color="C2", alpha=0.75)
    axCH.axvline(C_med, color="firebrick", lw=1.8, label=f"median C = {C_med:.2f}  ($C_L$ = {CL_med:.2f})")
    axCH.axvspan(C_med - C_std, C_med + C_std, color="firebrick", alpha=0.12, label="$\\pm1\\sigma$")
    axCH.set_xlabel("calibration coefficient C  ($C = 1/C_L$)"); axCH.set_ylabel("count")
    axCH.set_title("(f) calibration coefficient"); axCH.legend(fontsize=8, loc="upper right")

    fig.suptitle(f"Liquid-cloud calibration — Palaiseau CL31, {DAY}   "
                 f"(n={int(valid.sum())} profiles,  C = {C_med:.2f},  $r_e$ = {re_day:.1f} $\\mu$m)",
                 fontweight="bold", fontsize=14, y=0.98)
    for d in (OUT, REPO_FIGS):
        fig.savefig(d / "fig_cloud_calibration_example.png", dpi=160, bbox_inches="tight")
    plt.close(fig)
    print(f"C = {C_med:.3f} +- {C_std:.3f}  (C_L = {CL_med:.3f}),  n={int(valid.sum())}, "
          f"rejected shown={int(inval.sum())}")
    print("  fig_cloud_calibration_example.png")


if __name__ == "__main__":
    main()
