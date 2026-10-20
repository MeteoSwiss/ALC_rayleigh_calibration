"""Figure 4 of the paper (v0.9): representative daily liquid-cloud calibration, built from the
OPERATIONAL cloud pipeline (calibration.cloud.calibration, the code that produces the dashboard
constants) on one day: Palaiseau CL31, 2026-06-12.

Differences from validation/_showcase_cloud_figure.py (the v0.8 figure):
  * the accepted profiles are those of the operational pipeline, i.e. the cloud tests are applied to
    the profile AFTER the gate-wise multiple-scattering correction (the v0.8 script tested the
    uncorrected profile);
  * the per-profile coefficient C includes the sub-cloud aerosol two-way transmission correction
    (S_aer = 50 sr), as operationally (the v0.8 script omitted it);
  * the headline n and median C are asserted equal to liquid_cloud_calibration() for that day.
Panels (d) and (e) show the apparent lidar ratio WITHOUT the multiple-scattering correction (so the
height dependence is visible), scaled by the day's median coefficient before the aerosol
correction, with the theoretical eta(z)*S_liq of the operational table (a_G = 5.5 um, alpha =
10 km-1). Panel (f) is the operational per-profile C.
"""
from __future__ import annotations

import copy
import sys
import warnings
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
warnings.filterwarnings("ignore")

from calibration.cloud.calibration import (  # noqa: E402
    CloudCalConfig, set_defaults, liquid_cloud_calibration, _detect_cloud_input, build_cloud_input,
    compute_wv_transmission, apply_multiple_scattering_correction, apply_instrument_filters,
    calculate_lidar_ratio, apply_cloud_filters, apply_temporal_consistency_filter,
    apply_transmission_correction, S_THEORETICAL)
from calibration.io.data_loader import load_data, average_ceilometer_data  # noqa: E402
from validation.multiple_scattering_eta import eta_table, CBH_GRID_KM  # noqa: E402
import validation._showcase_cloud_figure as _SHOW  # noqa: E402
from validation._showcase_cloud_figure import cloudnet_re_profiles  # noqa: E402
_SHOW.CN = "D:/Cloudnet_Droplet_radius/20260612_palaiseau_der.nc"   # Cloudnet droplet radius archive

L1 = "D:/E-PROFILE_L1_2026/0-250-1001-07151/2026/06/L1_0-250-1001-07151_A20260612.nc"
DAY = "2026-06-12"
LAT, LON = 48.72, 2.21
OUT = Path("C:/DATA/Projects/202606_E-PROFILE_calibration/figs_paper_validation/paper_v09")
A_G_UM, ALPHA = 5.5, 10e-3            # operational multiple-scattering table parameters


def _cams_folder():
    for c in ("A:/CAMS_Monthly_04/", "D:/CAMS/"):
        if Path(c).exists() and any(Path(c).glob("CAMS_Beta_202606*.nc")):
            return c
    raise FileNotFoundError("no CAMS file for June 2026")


def operational_intermediates(cfg):
    """Re-run the steps of liquid_cloud_calibration in the same order, keeping the intermediates."""
    itype, level, units = _detect_cloud_input(cfg.nc_file, copy.deepcopy(cfg))
    cd = load_data([Path(cfg.nc_file)], itype, level)
    cd = average_ceilometer_data(cd, cfg.average_time_s, cfg.average_range_m)
    data = build_cloud_input(cd, copy.deepcopy(cfg), units, cd.calibration_constant_applied)
    trans2 = np.asarray(compute_wv_transmission(data, cfg), float)
    data.beta = data.beta / trans2
    data.trans2_wv = trans2
    beta_wv = np.asarray(data.beta, float).copy()
    beta_eta = apply_multiple_scattering_correction(data.beta, data.range, cfg)
    beta_f, _ = apply_instrument_filters(beta_eta, data, cfg)
    s_eta, _ = calculate_lidar_ratio(beta_f, data.range, cfg)
    s_cloud, _ = apply_cloud_filters(s_eta.copy(), beta_f, data, cfg)
    s_cons, _ = apply_temporal_consistency_filter(s_cloud, cfg)
    valid = np.isfinite(s_cons)
    c_aer, _, _ = apply_transmission_correction(beta_f, data, s_cons, cfg)
    # apparent S without the multiple-scattering correction (for the height dependence)
    beta_raw_f, _ = apply_instrument_filters(beta_wv.copy(), data, cfg)
    s_raw, _ = calculate_lidar_ratio(beta_raw_f, data.range, cfg)
    return dict(data=data, trans2=trans2, beta_wv=beta_wv, s_raw=np.asarray(s_raw, float),
                s_cons=np.asarray(s_cons, float), valid=valid, c_aer=np.asarray(c_aer, float))


def main():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.gridspec import GridSpec
    from mpl_toolkits.axes_grid1 import make_axes_locatable

    OUT.mkdir(parents=True, exist_ok=True)
    cfg = set_defaults(CloudCalConfig(
        nc_file=L1, instrument="CL31", apply_wv_correction=True,
        apply_transmission_correction=True, aerosol_lidar_ratio=50.0,
        cams_folder=_cams_folder(), abs_cs_lookup_table="",
        station_latitude=LAT, station_longitude=LON,
        average_time_s=30.0, average_range_m=10.0))
    res = liquid_cloud_calibration(copy.deepcopy(cfg))          # the operational result
    it = operational_intermediates(copy.deepcopy(cfg))
    data, valid = it["data"], it["valid"]
    c_prof = it["c_aer"][valid]
    assert c_prof.size == res.n_profiles, (c_prof.size, res.n_profiles)
    assert abs(np.median(c_prof) - res.cal_median) < 1e-9, (np.median(c_prof), res.cal_median)
    c_med, c_std = float(np.median(c_prof)), float(np.std(c_prof, ddof=1))
    cl_med = float(res.lidar_constant)                         # Wiegner C_L = 1e8 / C (CL31 L1)
    c_noaer_med = float(np.median(it["s_cons"][valid] / S_THEORETICAL))
    s_phys = it["s_raw"] / c_noaer_med                         # apparent S on the physical scale
    aer_factor = it["c_aer"][valid] / (it["s_cons"][valid] / S_THEORETICAL)

    rng = np.asarray(data.range, float)
    beta_wv, trans2 = it["beta_wv"], it["trans2"]
    cbh_all = np.asarray(data.cbh, float)
    if cbh_all.ndim == 2:
        cbh_all = cbh_all[:, 0]
    t = np.asarray(data.time)
    if np.issubdtype(t.dtype, np.datetime64):
        tsec = t.astype("datetime64[s]").astype(np.int64)
    else:                                                      # days since 1970-01-01
        tsec = np.round(np.asarray(t, float) * 86400.0).astype(np.int64)
    day0 = np.datetime64(f"{DAY}T00:00:00").astype("datetime64[s]").astype(np.int64)
    hours = (tsec - day0) / 3600.0
    cloud_present = np.isfinite(cbh_all) & (cbh_all >= 300) & (cbh_all <= 2600) & np.isfinite(s_phys)
    salt = float(getattr(data, "station_altitude", 157.0) or 157.0)
    _, re_day = cloudnet_re_profiles(tsec[valid], cbh_all[valid], salt)

    zc = np.asarray(CBH_GRID_KM, float)
    theo = eta_table("CL31", A_G_UM, alpha=ALPHA) * S_THEORETICAL
    theo_lo = eta_table("CL31", A_G_UM + 2.0, alpha=ALPHA) * S_THEORETICAL
    theo_hi = eta_table("CL31", A_G_UM - 2.0, alpha=ALPHA) * S_THEORETICAL

    vi = np.where(valid)[0]
    rep = vi[int(np.argmin(np.abs(cbh_all[vi] - np.median(cbh_all[vi]))))]
    rep_hhmm = f"{int(hours[rep]):02d}:{int(round((hours[rep] % 1) * 60)) % 60:02d}"

    fig = plt.figure(figsize=(17, 9.4))
    gs = GridSpec(2, 3, figure=fig, width_ratios=[2.4, 1.0, 1.3], height_ratios=[1, 1],
                  hspace=0.34, wspace=0.32)
    axPC = fig.add_subplot(gs[0, 0])
    axPr = fig.add_subplot(gs[0, 1], sharey=axPC)
    axWV = fig.add_subplot(gs[0, 2], sharey=axPC)
    axTS = fig.add_subplot(gs[1, 0], sharex=axPC)
    axOC = fig.add_subplot(gs[1, 1])
    axCH = fig.add_subplot(gs[1, 2])

    # (a) attenuated backscatter, cloud base, calibration profiles
    with np.errstate(all="ignore"):
        lb = np.log10(np.where(beta_wv > 0, beta_wv, np.nan))
    pm = axPC.pcolormesh(hours, rng / 1000.0, lb, cmap="viridis", vmin=-6.7, vmax=-4.5, shading="auto")
    axPC.scatter(hours[cloud_present], cbh_all[cloud_present] / 1000.0, s=9, c="k", marker=".",
                 edgecolor="none", alpha=0.55, label="cloud base", zorder=3)
    axPC.scatter(hours[valid], cbh_all[valid] / 1000.0, s=20, c="red", marker=".",
                 edgecolor="none", alpha=0.95, label="calibration profiles", zorder=4)
    axPC.axvline(hours[rep], color="red", ls=":", lw=1.0, alpha=0.7)
    axPC.set_ylim(0, 3.3)
    axPC.set_ylabel("range AGL [km]")
    axPC.set_xlabel("hours (UTC)")
    axPC.set_title("(a) attenuated backscatter $\\beta'$ (log$_{10}$) and cloud base")
    plt.colorbar(pm, cax=make_axes_locatable(axPC).append_axes("right", size="2.5%", pad=0.08),
                 label="log$_{10}\\beta'$")
    axPC.legend(loc="upper right", fontsize=8, framealpha=0.9)

    # (b) representative profile
    axPr.plot(beta_wv[:, rep], rng / 1000.0, "-", color="C0", lw=1.3)
    axPr.axhspan(0.1, 2.4, color="gold", alpha=0.12, label="integration 0.1-2.4 km")
    axPr.axhline(cbh_all[rep] / 1000.0, color="firebrick", ls="--", lw=1.2, label="cloud base")
    axPr.set_xlabel("$\\beta'$")
    axPr.set_title(f"(b) profile at {rep_hhmm} UTC")
    axPr.ticklabel_format(axis="x", style="sci", scilimits=(0, 0))
    axPr.legend(fontsize=7.5, loc="upper right")
    axPr.grid(alpha=0.3)

    # (c) water-vapour two-way transmission
    twv = np.nanmedian(trans2, axis=1)
    win = (rng >= 100) & (rng <= 2400)
    t_min, t_max = float(np.nanmin(trans2[win, :])), float(np.nanmax(trans2[win, :]))
    axWV.fill_betweenx(rng / 1000.0, np.nanmin(trans2, axis=1), np.nanmax(trans2, axis=1),
                       color="C4", alpha=0.20, label="daily min-max")
    axWV.plot(twv, rng / 1000.0, "-", color="C4", lw=1.7, label="daily median")
    axWV.axhline(cbh_all[rep] / 1000.0, color="firebrick", ls="--", lw=1.0, alpha=0.7)
    axWV.set_xlabel("$T^2_{WV}$ (two-way)")
    axWV.set_title("(c) water-vapour correction")
    axWV.set_xlim(max(0.4, t_min - 0.04), 1.01)
    axWV.grid(alpha=0.3)
    axWV.text(0.04, 0.045, "$\\beta' \\leftarrow \\beta'/T^2_{WV}$  (CAMS+HITRAN)\n"
              f"0.1-2.4 km: $T^2$ = {t_min:.2f} - {t_max:.2f}",
              transform=axWV.transAxes, fontsize=7.6, va="bottom",
              bbox=dict(boxstyle="round", fc="white", ec="0.7", alpha=0.85))
    axWV.legend(fontsize=7.2, loc="upper left")

    # (d) apparent lidar ratio time series (no multiple-scattering correction)
    inval = cloud_present & ~valid
    axTS.scatter(hours[inval], s_phys[inval], s=9, c="k", alpha=0.45, label="rejected profiles")
    axTS.scatter(hours[valid], s_phys[valid], s=12, c="C2", label="calibration profiles")
    axTS.axhline(S_THEORETICAL, color="k", ls="--", lw=1.2, label="single scatter 18.8 sr")
    axTS.axhline(np.median(s_phys[valid]), color="C2", lw=1.2, alpha=0.7,
                 label=f"median of calibration profiles {np.median(s_phys[valid]):.1f} sr")
    axTS.axvline(hours[rep], color="red", ls=":", lw=1.0, alpha=0.7)
    axTS.set_xlabel("hours (UTC)")
    axTS.set_ylabel("apparent S [sr]")
    axTS.set_ylim(8, 26)
    axTS.set_title("(d) apparent lidar ratio (no multiple-scattering correction)")
    axTS.grid(alpha=0.3)
    axTS.legend(fontsize=7.5, loc="upper right", ncol=2)
    make_axes_locatable(axTS).append_axes("right", size="2.5%", pad=0.08).axis("off")

    # (e) apparent S vs cloud-base height + eta*S of the operational table
    hb = axOC.hist2d(s_phys[valid], cbh_all[valid] / 1000.0,
                     bins=[np.linspace(8, 24, 45), np.linspace(0.3, 2.6, 40)], cmap="turbo", cmin=1)
    axOC.plot(theo, zc, "k-", lw=2.4, label=f"$\\eta S$, table $a_G$ = {A_G_UM} $\\mu$m")
    axOC.fill_betweenx(zc, theo_lo, theo_hi, color="k", alpha=0.12, label="$a_G \\pm 2$ $\\mu$m")
    axOC.axvline(S_THEORETICAL, color="0.25", ls="--", lw=1.4)
    axOC.text(S_THEORETICAL + 0.2, 2.5, "single\nscatter\n18.8 sr", color="0.25", fontsize=7.5, va="top")
    axOC.text(0.03, 0.03, f"Cloudnet $r_e$ at the\ncalibration scenes: {re_day:.1f} $\\mu$m",
              transform=axOC.transAxes, fontsize=7.5, va="bottom",
              bbox=dict(boxstyle="round", fc="white", ec="0.7", alpha=0.85))
    axOC.set_xlim(8, 24)
    axOC.set_ylim(0, 3.3)
    axOC.set_xlabel("apparent lidar ratio S [sr]")
    axOC.set_ylabel("cloud base height [km]")
    axOC.set_title("(e) apparent S vs height and $\\eta S$")
    axOC.legend(fontsize=7.5, loc="upper right")
    plt.colorbar(hb[3], cax=make_axes_locatable(axOC).append_axes("right", size="3%", pad=0.08),
                 label="count")

    # (f) operational per-profile calibration coefficient
    axCH.hist(c_prof, bins=np.linspace(c_med - 4 * c_std, c_med + 4 * c_std, 40), color="C2", alpha=0.75)
    axCH.axvline(c_med, color="firebrick", lw=1.8, label=f"median C = {c_med:.2f}")
    axCH.axvspan(c_med - c_std, c_med + c_std, color="firebrick", alpha=0.12, label="$\\pm1\\sigma$")
    axCH.set_xlabel("calibration coefficient C")
    axCH.set_ylabel("count")
    axCH.set_title("(f) calibration coefficient")
    axCH.legend(fontsize=8, loc="upper right")

    fig.suptitle(f"Liquid-cloud calibration — Palaiseau CL31, {DAY}   (n = {c_prof.size} profiles,  "
                 f"C = {c_med:.2f},  $C_L$ = {cl_med:.2e})", fontweight="bold", fontsize=14, y=0.98)
    out = OUT / "fig_v09_cloud_example.png"
    fig.savefig(out, dpi=160, bbox_inches="tight")
    plt.close(fig)
    print(f"saved {out}")
    print(f"n = {c_prof.size}, C median = {c_med:.4f} (operational {res.cal_median:.4f}), std = {c_std:.4f}, "
          f"C_L = {cl_med:.4e}, C without aerosol correction = {c_noaer_med:.4f}, "
          f"aerosol factor median {np.median(aer_factor):.3f} ({aer_factor.min():.3f}-{aer_factor.max():.3f}), "
          f"Cloudnet r_e {re_day:.2f} um, rejected shown {int(inval.sum())}")


if __name__ == "__main__":
    main()
