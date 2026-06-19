"""
compare_lindenberg_chm_cl61.py — §7.6 Lindenberg: E-PROFILE CHM15k vs ACTRIS-Cloudnet CL61.

Compares calibrated β_att profiles at Lindenberg:
  - CHM15k (1064 nm): E-PROFILE L2 monthly, Rayleigh calibration from Python CSV + Kalman.
  - CL61  ( 910 nm): ACTRIS-Cloudnet daily files, corrected for factory overcalibration
    via the Rayleigh Kalman estimate (run_lindenberg_cl61_cal.py), then WV-corrected (CAMS).

Corrections applied to CL61 before comparison:
  0. Rayleigh Kalman calibration correction: β_true = β_Cloudnet / CL_kalman(t).
     Cloudnet's O'Connor factory calibration is overcalibrated by ~2× at Lindenberg
     (CL_kalman ≈ 1.93 from Rayleigh against molecular background).
  1. Two-way water-vapour absorption at 910.74 nm (CAMS monthly means).
  2. Wavelength scaling 910 → 1064 nm: two-component aerosol (Ångström α=1.5) + molecular
     (King-factor ratio).

Output: comparison figure (bias/scatter/profile overlay) and updated §7.6 stats.
"""
from __future__ import annotations
import os, sys, warnings
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[_v] = "1"
try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass

from pathlib import Path
from datetime import date, timedelta
from collections import defaultdict
import numpy as np
import pandas as pd
import netCDF4 as nc4
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.dates as mdates

# ----- paths -----------------------------------------------------------------
L2_ROOT  = Path("A:/E-PROFILE_L2_monthly")
CL61_ROOT = Path("A:/CL61_Cloudnet/Lindenberg")
CL61_KALMAN_CSV = Path(
    "C:/DATA/Projects/202606_E-PROFILE_calibration"
    "/figs_paper_validation/lindenberg_cl61_cal/rayleigh_lindenberg_cl61_kalman.csv"
)
CAL_CSV  = Path("C:/DATA/Projects/202606_E-PROFILE_calibration/E-PROFILE_calibration_rayleigh/fullcal_all/0-20000-0-10393_0/0-20000-0-10393_0_cl.csv")
CAMS_DIR = Path("D:/CAMS")
OUT_DIR  = Path("C:/DATA/Projects/202606_E-PROFILE_calibration/figs_paper_validation")

WMO = "0-20000-0-10393"
IDENT = "0"

ANG_EXP  = 1.5           # aerosol Ångström exponent 910→1064 nm
LAMBDA_CL = 910.74        # nm CL61
LAMBDA_CH = 1064.47       # nm CHM15k
Z_MIN, Z_MAX = 300.0, 8000.0   # comparison range AGL (m)
DT_BIN   = 30             # minutes: time bin for aggregation


# =============================================================================
#  Helpers
# =============================================================================
def _king_ratio(lam_from_nm, lam_to_nm):
    """Ratio of King-factor (F≈1 for air) × (lam⁻⁴) scattering between two wavelengths.
    beta_mol(lam_to) / beta_mol(lam_from) ≈ (lam_from/lam_to)^4 (near λ⁻⁴ law)."""
    return (lam_from_nm / lam_to_nm) ** 4


def _wavelength_correct(beta_910, beta_mol_910, alpha_ang=ANG_EXP):
    """Convert CL61 β_att(910 nm) → β_att(1064 nm) profile-wise.

    Two components:
      β_mol: scale by (910/1064)⁴.
      β_aer: subtract molecular, scale by (910/1064)^α, re-add molecular_1064.
    """
    r_mol = _king_ratio(LAMBDA_CL, LAMBDA_CH)
    r_aer = (LAMBDA_CL / LAMBDA_CH) ** alpha_ang
    beta_mol_1064 = beta_mol_910 * r_mol
    beta_aer_910  = np.maximum(beta_910 - beta_mol_910, 0.0)
    beta_aer_1064 = beta_aer_910 * r_aer
    return beta_aer_1064 + beta_mol_1064


# =============================================================================
#  Step 1: load CHM15k E-PROFILE L2 + apply Rayleigh calibration
# =============================================================================
def load_chm15k_calibrated(start_ym="202401", end_ym="202603"):
    """Read L2 monthly files and apply the Kalman-smoothed Rayleigh correction.

    Returns (time_dt, height_m, beta_att) with beta in m⁻¹ sr⁻¹.
    """
    # load calibration CSV
    cal = pd.read_csv(CAL_CSV, parse_dates=["date"])
    cal = cal[cal["flag"].isin([1, 1.0, 0.5]) & (cal["lidar_constant"] > 0)].copy()
    cal["date"] = pd.to_datetime(cal["date"])
    cal = cal.sort_values("date")

    # Reject outliers via log-MAD>4 (same as load_rayleigh_python_kalman.m — CSV has
    # garbage CL up to 1e13+ that corrupts a linear rolling mean).
    from scipy.ndimage import uniform_filter1d
    cl_vals = cal["lidar_constant"].values.astype(float)
    log_cl = np.log(cl_vals)
    log_med = np.nanmedian(log_cl)
    log_mad = np.nanmedian(np.abs(log_cl - log_med)) * 1.4826
    good = np.abs(log_cl - log_med) < 4 * log_mad
    cl_clean = np.where(good, cl_vals, np.nan)
    # interpolate over rejected gaps
    fi = np.where(np.isfinite(cl_clean))[0]
    if fi.size > 1:
        cl_clean = np.interp(np.arange(len(cl_clean)), fi, cl_clean[fi])
    else:
        cl_clean[:] = np.exp(log_med)
    # smooth in log space to stay robust against residual scatter
    cl_smooth = np.exp(uniform_filter1d(np.log(cl_clean), size=min(7, len(cl_clean)), mode="nearest"))
    cal["cl_smooth"] = cl_smooth
    print(f"  Lindenberg CHM15k CSV: {len(cal)} nights, cl_med={np.exp(log_med):.3e}, "
          f"outliers_rejected={int((~good).sum())}", flush=True)

    # collect all L2 monthly files
    ym_start = int(start_ym)
    ym_end   = int(end_ym)
    times_all, betas_all, heights_ref = [], [], None
    for yr in range(ym_start // 100, ym_end // 100 + 1):
        yr_dir = L2_ROOT / WMO / str(yr)
        if not yr_dir.is_dir():
            continue
        for nc_f in sorted(yr_dir.glob(f"L2_{WMO}_{IDENT}*.nc")):
            ym = int(nc_f.stem[-6:])
            if ym < ym_start or ym > ym_end:
                continue
            try:
                _read_l2_month(nc_f, cal, times_all, betas_all, heights_ref)
                if heights_ref is None and betas_all:
                    heights_ref = _last_heights
            except Exception as e:
                print(f"  CHM15k skip {nc_f.name}: {e}")

    if not betas_all:
        return None, None, None
    time_arr = np.concatenate(times_all)
    beta_arr = np.concatenate(betas_all, axis=0)
    return time_arr, heights_ref, beta_arr


_last_heights = None


def _read_l2_month(nc_f, cal_df, times_out, betas_out, heights_ref):
    """Append one L2 monthly file's data to time and beta lists."""
    global _last_heights
    with nc4.Dataset(str(nc_f), "r") as ds:
        t_units = ds.variables["time"].units
        t_raw = np.asarray(ds.variables["time"][:], dtype="f8")
        # parse time
        import cftime
        try:
            t_dt = nc4.num2date(t_raw, t_units, only_use_cftime_datetimes=False,
                                only_use_python_datetimes=True)
        except Exception:
            t_dt = nc4.num2date(t_raw, t_units)
        t_dt = np.asarray([np.datetime64(t, "ns") for t in t_dt])

        # height — try range/height (AGL) then altitude (ASL → subtract station_altitude)
        if "range" in ds.variables:
            hgt = np.asarray(ds.variables["range"][:], dtype="f8")
        elif "height" in ds.variables:
            hgt = np.asarray(ds.variables["height"][:], dtype="f8")
        elif "altitude" in ds.variables:
            alt_asl = np.asarray(ds.variables["altitude"][:], dtype="f8")
            sa = 0.0
            if "station_altitude" in ds.variables:
                sa = float(np.asarray(ds.variables["station_altitude"][:]).ravel()[0])
            hgt = alt_asl - sa
        else:
            raise ValueError("no height variable")
        _last_heights = hgt

        # calibration_constant_0 (per-profile or scalar) — take median for the month
        cl_l2 = None
        if "calibration_constant_0" in ds.variables:
            cl_raw = np.asarray(ds.variables["calibration_constant_0"][:], dtype="f8")
            cl_l2 = float(np.nanmedian(cl_raw[cl_raw > 0])) if cl_raw.size > 1 else float(cl_raw.ravel()[0])

        # attenuated_backscatter_0 [time, height], units ~µm⁻¹ sr⁻¹
        bv = np.asarray(ds.variables["attenuated_backscatter_0"][:], dtype="f8")
        if bv.ndim == 2 and bv.shape[0] != len(t_dt):
            bv = bv.T  # ensure (time, height)

        # apply L2 scale to get β_att in m⁻¹ sr⁻¹
        beta = bv * 1e-6   # stored in 10⁶ m⁻¹ sr⁻¹ (micro-scaled)

        # apply Rayleigh correction: scale by CL_L2 / CL_kalman
        if cl_l2 is not None and cl_l2 > 0:
            # find closest calibration to the mid-month date
            mid = pd.Timestamp(t_dt[len(t_dt)//2])
            idx = (cal_df["date"] - mid).abs().idxmin()
            cl_kal = float(cal_df.loc[idx, "cl_smooth"])
            if cl_kal > 0:
                correction = cl_l2 / cl_kal
                beta = beta * correction

    times_out.append(t_dt)
    betas_out.append(beta)


# =============================================================================
#  Step 2: load CL61 Cloudnet daily files + WV correction
# =============================================================================
def load_cl61_factory(start_ym="202401", end_ym="202603"):
    """Read CL61 daily files, restrict height, pre-bin to DT_BIN-min medians per day.

    Returning already-binned data avoids loading the full ~74 GB into RAM simultaneously.
    No WV correction here — applied later per month using CAMS.
    """
    ym_start, ym_end = int(start_ym), int(end_ym)
    times_all, betas_all, height_ref = [], [], None
    for nc_f in sorted(CL61_ROOT.glob("????????.nc")):
        ymd = nc_f.stem
        ym = int(ymd[:6])
        if ym < ym_start or ym > ym_end:
            continue
        try:
            t_dt, hgt, beta = _read_cl61_daily(nc_f)
            if height_ref is None:
                height_ref = hgt
            times_all.append(t_dt)
            betas_all.append(beta)
        except Exception as e:
            print(f"  CL61 skip {nc_f.name}: {e}")
    if not betas_all:
        return None, None, None
    return np.concatenate(times_all), height_ref, np.concatenate(betas_all, axis=0)


def _read_cl61_daily(nc_f):
    """Load one CL61 day, apply height mask Z_MIN/Z_MAX, return DT_BIN-min medians (float32)."""
    import cftime
    with nc4.Dataset(str(nc_f), "r") as ds:
        t_raw = np.asarray(ds.variables["time"][:], dtype="f8")
        t_units = ds.variables["time"].units
        try:
            t_dt = nc4.num2date(t_raw, t_units, only_use_cftime_datetimes=False,
                                only_use_python_datetimes=True)
        except Exception:
            t_dt = nc4.num2date(t_raw, t_units)
        t_dt = np.asarray([np.datetime64(t, "ns") for t in t_dt])

        hgt_full = np.asarray(ds.variables["range"][:], dtype="f4")
        h0 = int(np.searchsorted(hgt_full, Z_MIN))
        h1 = int(np.searchsorted(hgt_full, Z_MAX, side="right"))
        hgt = hgt_full[h0:h1]

        # Slice read avoids loading unused height gates (halves disk IO vs boolean mask)
        beta = np.asarray(ds.variables["beta_att"][:, h0:h1], dtype="f4")
        beta[~np.isfinite(beta) | (beta < 0)] = np.nan

    # Pre-bin to DT_BIN-min medians so concatenation across 330 days stays small
    dt_ns = DT_BIN * 60 * 1_000_000_000
    t_ns = t_dt.astype("int64")
    t_bin_ns = (t_ns // dt_ns) * dt_ns
    unique_bins = np.unique(t_bin_ns)
    out = np.full((len(unique_bins), len(hgt)), np.nan, dtype="f4")
    for i, tb in enumerate(unique_bins):
        idx = (t_bin_ns == tb)
        if idx.sum() > 0:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", RuntimeWarning)
                out[i] = np.nanmedian(beta[idx], axis=0)
    return unique_bins.astype("datetime64[ns]"), hgt, out


# =============================================================================
#  Step 3: WV correction for CL61 (CAMS monthly means)
# =============================================================================
def apply_wv_correction_monthly(beta_cl61, time_cl61, height_cl61,
                                 lat=52.21, lon=14.12, alt=123.0):
    """Apply monthly-mean WV correction from CAMS to CL61 β_att.

    Uses the cloud_calibration WV chain (identical to the MATLAB pipeline).
    Returns beta_corrected and trans2 (for diagnostics).
    """
    import json as _json
    from calibration.cloud.calibration import (
        CloudCalConfig, set_defaults, compute_wv_transmission, CeiloData, _matlab_datenum)

    _opt = _json.loads(Path("options.json").read_text())
    # group by month
    t_pd = pd.to_datetime(time_cl61)
    months = sorted(set(t_pd.to_period("M")))
    beta_out = beta_cl61.copy()

    for ym in months:
        mask = np.asarray(t_pd.to_period("M") == ym)
        if not mask.any():
            continue
        # build minimal CeiloData
        t_sub = time_cl61[mask]
        b_sub = beta_cl61[mask]  # (n_time, n_range)
        cfg = CloudCalConfig(instrument="CL61",
                             apply_wv_correction=True,
                             cams_folder=str(CAMS_DIR),
                             abs_cs_lookup_table=_opt.get("abs_cs_lookup_table", ""),
                             station_latitude=lat, station_longitude=lon)
        cfg = set_defaults(cfg)
        data = CeiloData(
            time=t_sub,
            time_num=_matlab_datenum(t_sub),
            station_altitude=alt,
            station_latitude=lat,
            station_longitude=lon,
            range=height_cl61,
            range_resol=float(height_cl61[1] - height_cl61[0]) if len(height_cl61) > 1 else 10.0,
            beta=b_sub.T,   # CeiloData expects (n_range, n_time)
            cbh=np.full(len(t_sub), np.nan),
            quality_flag=None,
            window_transmission=None,
            laser_energy=None,
        )
        try:
            trans2 = compute_wv_transmission(data, cfg)
            trans2 = np.asarray(trans2, dtype="f8")
            if trans2.shape == b_sub.T.shape:
                beta_out[mask] = b_sub / trans2.T
        except Exception as e:
            print(f"  WV correction skip {ym}: {e}")

    return beta_out


# =============================================================================
#  Step 3b: apply Rayleigh Kalman CL61 calibration correction
# =============================================================================
def load_cl61_rayleigh_kalman(csv_path: Path):
    """Load the E-PROFILE Kalman CL_L (Wiegner 2014) series from run_lindenberg_cl61_cal output.

    Returns (dates_np, kalman_wv) where dates_np are numpy datetime64[ns] timestamps at noon
    on each calendar day and kalman_wv is the dimensionless C_L series (>1 → Cloudnet overcal.).
    """
    if not csv_path.is_file():
        raise FileNotFoundError(f"CL61 Rayleigh Kalman CSV not found: {csv_path}\n"
                                "Run run_lindenberg_cl61_cal.py first.")
    df = pd.read_csv(csv_path, parse_dates=["date"])
    df = df[df["kalman_wv"].notna() & (df["kalman_wv"] > 0)].copy()
    df = df.sort_values("date")
    dates_np = np.array([np.datetime64(d, "ns") for d in df["date"]])
    return dates_np, df["kalman_wv"].values.astype("f8")


def apply_cl61_calib_correction(beta_factory, time_cl61_ns, cl_dates_ns, cl_kalman):
    """Divide Cloudnet factory β_att by the interpolated Rayleigh Kalman C_L.

    beta_true(t) = beta_Cloudnet(t) / CL_kalman(t)

    CL_kalman > 1 ⟹ Cloudnet is overcalibrated; dividing brings it in line with
    the Rayleigh-derived scale.  The daily Kalman values are linearly interpolated
    to the 30-min β_att time stamps.
    """
    # convert everything to seconds for interpolation
    t_ns  = time_cl61_ns.astype("f8")
    cl_t  = cl_dates_ns.astype("f8")
    cl_v  = cl_kalman.astype("f8")

    if cl_t.size < 2:
        print("  WARNING: fewer than 2 Kalman points — correction not applied.")
        return beta_factory

    cl_interp = np.interp(t_ns, cl_t, cl_v, left=cl_v[0], right=cl_v[-1])
    # cl_interp shape: (n_time,) — broadcast over height
    beta_corrected = beta_factory / cl_interp[:, None]
    med_cl = float(np.nanmedian(cl_interp))
    print(f"  CL61 Rayleigh correction: median CL_kalman = {med_cl:.3f}  "
          f"→ β reduced by factor {med_cl:.2f}×", flush=True)
    return beta_corrected


# =============================================================================
#  Step 4: aggregate to common time-height grid
# =============================================================================
def bin_to_grid(time_dt, height_m, beta, dt_min=DT_BIN):
    """Bin (time, height) β_att to dt_min-minute medians.

    Returns (time_bins_ns, height_m, beta_median).
    """
    t_ns = time_dt.astype("int64")
    dt_ns = dt_min * 60 * 1_000_000_000
    t_bin_ns = (t_ns // dt_ns) * dt_ns

    unique_bins = np.unique(t_bin_ns)
    out = np.full((len(unique_bins), len(height_m)), np.nan)
    for i, tb in enumerate(unique_bins):
        idx = (t_bin_ns == tb)
        if idx.sum() > 0:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", RuntimeWarning)
                out[i] = np.nanmedian(beta[idx], axis=0)
    return unique_bins.astype("datetime64[ns]"), height_m, out


# =============================================================================
#  Step 5: synchronize + interpolate to common altitude grid
# =============================================================================
def sync_and_interp(t_chm, z_chm, beta_chm, t_cl61, z_cl61, beta_cl61):
    """Match CHM15k and CL61 to the same times and heights.

    Returns (t_common, z_common, beta_chm_sync, beta_cl61_sync).
    """
    # common height grid = CHM15k grid (keep native), interpolate CL61 onto it
    z_common = z_chm[(z_chm >= Z_MIN) & (z_chm <= Z_MAX)]

    # interpolate CL61 to CHM15k height grid
    def interp_heights(z_src, beta_src, z_dst):
        out = np.full((beta_src.shape[0], len(z_dst)), np.nan)
        for i in range(len(z_dst)):
            if z_dst[i] < z_src[0] or z_dst[i] > z_src[-1]:
                continue
            idx = np.searchsorted(z_src, z_dst[i])
            if idx == 0:
                out[:, i] = beta_src[:, 0]
            elif idx >= len(z_src):
                out[:, i] = beta_src[:, -1]
            else:
                w = (z_dst[i] - z_src[idx-1]) / (z_src[idx] - z_src[idx-1])
                out[:, i] = (1-w) * beta_src[:, idx-1] + w * beta_src[:, idx]
        return out

    z_mask_chm = (z_chm >= Z_MIN) & (z_chm <= Z_MAX)
    z_mask_cl61 = (z_cl61 >= Z_MIN) & (z_cl61 <= Z_MAX)

    # match time axis: find common time stamps
    t_chm_set = {int(t): i for i, t in enumerate(t_chm.astype("int64"))}
    common_idx_chm, common_idx_cl61 = [], []
    for j, tc in enumerate(t_cl61.astype("int64")):
        if tc in t_chm_set:
            common_idx_chm.append(t_chm_set[tc])
            common_idx_cl61.append(j)

    if not common_idx_chm:
        return None, None, None, None

    t_common = t_chm[common_idx_chm]
    b_chm_sub = beta_chm[common_idx_chm][:, z_mask_chm]
    b_cl61_sub_raw = beta_cl61[common_idx_cl61][:, z_mask_cl61]

    # interpolate CL61 onto CHM15k z grid
    z_cl61_sub = z_cl61[z_mask_cl61]
    b_cl61_sync = interp_heights(z_cl61_sub, b_cl61_sub_raw, z_common)

    return t_common, z_common, b_chm_sub, b_cl61_sync


# =============================================================================
#  Step 6: statistics + figure
# =============================================================================
def compute_stats(beta_ref, beta_comp, label="CL61"):
    """Compute bias, RMSE, r on matched (time, height) log-β pairs."""
    mask = np.isfinite(beta_ref) & np.isfinite(beta_comp) & (beta_ref > 0) & (beta_comp > 0)
    if mask.sum() < 100:
        return {}
    lr = np.log10(beta_ref[mask])
    lc = np.log10(beta_comp[mask])
    bias_log = float(np.mean(lc - lr))
    rmse_log  = float(np.sqrt(np.mean((lc - lr)**2)))
    r  = float(np.corrcoef(lr, lc)[0, 1])
    bias_pct = float((10**bias_log - 1) * 100)
    n = int(mask.sum())
    return dict(bias_log=bias_log, bias_pct=bias_pct, rmse_log=rmse_log, r=r, n=n, label=label)


def make_figure(t_common, z_common, b_chm, b_cl61, stats, out_path):
    """Generate the comparison figure."""
    fig = plt.figure(figsize=(18, 7))
    fig.patch.set_facecolor("white")
    gs = fig.add_gridspec(1, 4, wspace=0.35)
    ax_prof = fig.add_subplot(gs[0])
    ax_chm  = fig.add_subplot(gs[1])
    ax_cl61 = fig.add_subplot(gs[2])
    ax_sc   = fig.add_subplot(gs[3])

    CLIM = (-7.5, -4.5)
    CMAP = "plasma"

    # median profiles
    med_chm  = np.nanmedian(b_chm,  axis=0)
    med_cl61 = np.nanmedian(b_cl61, axis=0)
    ax_prof.semilogx(med_chm  * 1e6, z_common / 1e3, "b-",  lw=1.8, label="CHM15k 1064 nm")
    ax_prof.semilogx(med_cl61 * 1e6, z_common / 1e3, "r--", lw=1.8,
                     label="CL61 910 nm (cal.+WV+λ corr.)")
    ax_prof.set_xlabel("β_att (Mm⁻¹ sr⁻¹)", fontsize=9)
    ax_prof.set_ylabel("Height AGL (km)", fontsize=9)
    ax_prof.set_title("Median profiles", fontsize=9)
    ax_prof.legend(fontsize=7)
    ax_prof.grid(True, alpha=0.3)
    ax_prof.set_xlim(1e-3, 1e2)
    ax_prof.set_ylim(0, Z_MAX / 1e3)

    # time-height pcolor (log10 β, Mm⁻¹ sr⁻¹)
    t_num = t_common.astype("float64") * 1e-9  # seconds since epoch
    for ax, data, title in [(ax_chm, b_chm, "CHM15k (1064 nm)"),
                             (ax_cl61, b_cl61, "CL61 (910 nm, cal.+WV+λ corr.)")]:
        ldata = np.log10(data * 1e6 + 1e-10)
        pc = ax.pcolormesh(t_num, z_common / 1e3, ldata.T,
                           cmap=CMAP, vmin=CLIM[0], vmax=CLIM[1], shading="auto")
        ax.set_ylabel("Height AGL (km)", fontsize=8)
        ax.set_title(title, fontsize=8)
        plt.colorbar(pc, ax=ax, label="log₁₀ β_att [Mm⁻¹ sr⁻¹]", shrink=0.8)
        # x-axis: datetime labels
        t_dts = pd.to_datetime(t_common).to_list()
        tick_every = max(1, len(t_common) // 8)
        ax.set_xticks(t_num[::tick_every])
        ax.set_xticklabels([pd.Timestamp(t).strftime("%m-%d") for t in t_common[::tick_every]],
                            rotation=30, ha="right", fontsize=6)
        ax.set_ylim(0, Z_MAX / 1e3)

    # scatter
    mask = np.isfinite(b_chm) & np.isfinite(b_cl61) & (b_chm > 0) & (b_cl61 > 0)
    if mask.sum() > 0:
        lr = np.log10(b_chm[mask] * 1e6)
        lc = np.log10(b_cl61[mask] * 1e6)
        # hexbin
        hb = ax_sc.hexbin(lr, lc, gridsize=60, cmap="YlOrRd", mincnt=5,
                          norm=matplotlib.colors.LogNorm())
        plt.colorbar(hb, ax=ax_sc, label="counts", shrink=0.8)
        lim = (max(np.nanpercentile(lr, 1), -7), min(np.nanpercentile(lr, 99), -4))
        ax_sc.plot(lim, lim, "k--", lw=1, label="1:1")
        ax_sc.set_xlim(lim); ax_sc.set_ylim(lim)
        bias_s = stats.get("bias_pct", float("nan"))
        r_s    = stats.get("r", float("nan"))
        n_s    = stats.get("n", 0)
        ax_sc.text(0.05, 0.93, f"bias={bias_s:+.0f}%  r={r_s:.2f}  n={n_s/1e3:.0f}k",
                   transform=ax_sc.transAxes, fontsize=8,
                   bbox=dict(boxstyle="round,pad=0.3", fc="white", alpha=0.8))
        ax_sc.set_xlabel("log₁₀ β_att CHM15k", fontsize=9)
        ax_sc.set_ylabel("log₁₀ β_att CL61", fontsize=9)
        ax_sc.set_title("Scatter (all times/heights)", fontsize=9)
        ax_sc.legend(fontsize=8)
        ax_sc.grid(True, alpha=0.3)

    fig.suptitle(
        "Lindenberg — E-PROFILE CHM15k 1064 nm vs ACTRIS-Cloudnet CL61 910 nm\n"
        f"(Rayleigh cal. corr. + WV + wavelength-scaled, {DT_BIN}-min medians)",
        fontsize=11)
    fig.tight_layout()
    fig.savefig(out_path, dpi=140, bbox_inches="tight")
    plt.close(fig)
    print(f"saved {out_path}")


# =============================================================================
#  Main
# =============================================================================
def main():
    warnings.filterwarnings("ignore")
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    print("Loading CHM15k E-PROFILE L2 + Rayleigh calibration...", flush=True)
    t_chm, z_chm, beta_chm_raw = load_chm15k_calibrated()
    if t_chm is None:
        print("No CHM15k data found — check L2_ROOT path."); return

    print(f"  {len(t_chm)} CHM15k profiles over {len(z_chm)} gates", flush=True)

    print("Loading CL61 Cloudnet daily files...", flush=True)
    t_cl61, z_cl61, beta_cl61_factory = load_cl61_factory()
    if t_cl61 is None:
        print("No CL61 Cloudnet files found — check CL61_ROOT path."); return

    print(f"  {len(t_cl61)} CL61 profiles over {len(z_cl61)} gates", flush=True)

    print("Loading CL61 Rayleigh Kalman correction (Wiegner 2014 C_L)...", flush=True)
    try:
        cl_dates, cl_kalman = load_cl61_rayleigh_kalman(CL61_KALMAN_CSV)
        print(f"  Loaded {len(cl_dates)} daily C_L values, "
              f"median CL = {np.nanmedian(cl_kalman):.3f}", flush=True)
        beta_cl61_cal = apply_cl61_calib_correction(
            beta_cl61_factory, t_cl61, cl_dates, cl_kalman)
    except FileNotFoundError as e:
        print(f"  WARNING: {e}\n  Using uncorrected Cloudnet β_att (overcalibrated ~2×).")
        beta_cl61_cal = beta_cl61_factory

    print("Applying WV correction to CL61 (CAMS 910 nm)...", flush=True)
    beta_cl61_wv = apply_wv_correction_monthly(
        beta_cl61_cal, t_cl61, z_cl61, lat=52.21, lon=14.12, alt=123.0)

    print("Computing molecular β at 910 nm for wavelength correction...", flush=True)
    try:
        from calibration.rayleigh.atmosphere import load_standard_atmosphere, calculate_molecular_properties
        from pathlib import Path as _Path
        std_atm_file = None  # use the US Standard Atmosphere shipped as package data
        atm = load_standard_atmosphere(std_atm_file, z_cl61)
        mol = calculate_molecular_properties(atm.temperature, atm.pressure, z_cl61,
                                             wavelength_m=LAMBDA_CL * 1e-9)
        beta_mol_cl61 = mol.beta_mol
    except Exception as e:
        print(f"  standard atm fallback ({e})")
        # exponential air density proxy
        beta_mol_cl61 = 1.2e-6 * _king_ratio(LAMBDA_CH, LAMBDA_CL) * np.exp(-z_cl61 / 8500)

    print("Applying wavelength correction 910 → 1064 nm...", flush=True)
    beta_mol_910_grid = np.broadcast_to(beta_mol_cl61[None, :], beta_cl61_wv.shape).copy()
    beta_cl61_1064 = _wavelength_correct(beta_cl61_wv, beta_mol_910_grid)

    # CL61 is already pre-binned to DT_BIN-min medians per day; CHM15k is 5-min L2 → bin it
    print(f"Binning CHM15k to {DT_BIN}-min medians...", flush=True)
    t_chm_bin, z_chm_b, beta_chm_bin = bin_to_grid(t_chm, z_chm, beta_chm_raw)
    t_cl61_bin, z_cl61_b, beta_cl61_bin = t_cl61, z_cl61, beta_cl61_1064

    print("Synchronising time+height grids...", flush=True)
    t_comm, z_comm, b_chm_sync, b_cl61_sync = sync_and_interp(
        t_chm_bin, z_chm_b, beta_chm_bin, t_cl61_bin, z_cl61_b, beta_cl61_bin)

    if t_comm is None or len(t_comm) == 0:
        print("No common time steps — check that the overlap period has data."); return

    print(f"  {len(t_comm)} common 30-min bins over {len(z_comm)} heights", flush=True)

    stats = compute_stats(b_chm_sync, b_cl61_sync, label="CL61 910 nm (cal.+WV+λ corr.)")
    print(f"  bias={stats.get('bias_pct',float('nan')):+.0f}%  "
          f"r={stats.get('r',float('nan')):.3f}  n={stats.get('n',0)/1e3:.0f}k", flush=True)

    out_fig = OUT_DIR / "lindenberg_chm_cl61_comparison.png"
    make_figure(t_comm, z_comm, b_chm_sync, b_cl61_sync, stats, out_fig)
    print("LINDENBERG_COMPARE_DONE")


if __name__ == "__main__":
    main()
