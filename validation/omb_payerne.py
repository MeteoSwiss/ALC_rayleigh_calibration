"""OmB spot-check on Payerne L1 vs CAMS (both operational and our-calibrated).

Validates the calibration/omb path before the network run:
* loads the operational L1 we already have (Payerne 0-20000-0-06610, May 2026),
* converts raw rcs_0 -> attenuated backscatter two ways:
    - "op"   : / operational constant (L2 calibration_constant_0)
    - "ours" : / our Kalman best-estimate C_L (fullcal_l1_2026/<key>_kalman.csv),
* runs compute_omb (3 h x 150 m averaging, cloud screen, 910 nm Angstrom +
  water-vapour correction reusing compute_wv_transmission),
* writes per-stream figures + a report.

Run:  python validation/omb_payerne.py
Outputs: doc/reports/omb_payerne/figs/*.png and doc/reports/omb_payerne.md
"""

from __future__ import annotations

import glob
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import matplotlib

matplotlib.use("Agg")
import matplotlib.dates as mdates  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
import netCDF4  # noqa: E402

from calibration.cloud.calibration import INSTRUMENT_CAL_DEFAULT  # noqa: E402
from calibration.omb.omb import compute_omb  # noqa: E402

# ----------------------------------------------------------------------------
# Configuration
# ----------------------------------------------------------------------------
STATION = "0-20000-0-06610"
STATION_NAME = "Payerne"
YEAR, MONTH = 2026, 5
CAMS_FILE = f"D:/CAMS/CAMS_Beta_{YEAR}{MONTH:02d}.nc"
L1_DIR = Path(f"D:/E-PROFILE_L1_2026/{STATION}/{YEAR}/{MONTH:02d}")
L2_DIR = Path(f"D:/E-PROFILE_L2_2026/{STATION}/{YEAR}/{MONTH:02d}")
FULLCAL = Path("C:/DATA/Projects/202606_E-PROFILE_calibration/fullcal_l1_2026")

# stream ident -> (instrument, Kalman method used as the authoritative constant)
STREAMS = {"A": ("CHM15k", "rayleigh"), "B": ("CL31", "cloud"), "C": ("CL61", "rayleigh")}

UNIT = 1e6  # m^-1 sr^-1 -> Mm^-1 sr^-1 for display
OUT_DIR = Path("doc/reports/omb_payerne")
FIG_DIR = OUT_DIR / "figs"
REPORT = Path("doc/reports/omb_payerne.md")

COL = {"ours": (0.85, 0.10, 0.10), "op": (0.25, 0.25, 0.25),
       "ours_wv": (0.10, 0.30, 0.85), "cams": (0.0, 0.0, 0.0)}


# ----------------------------------------------------------------------------
# Loading
# ----------------------------------------------------------------------------
def load_l1_month(ident: str):
    """Concatenate the month's L1 daily files for one stream."""
    files = sorted(glob.glob(str(L1_DIR / f"L1_{STATION}_{ident}{YEAR}{MONTH:02d}*.nc")))
    if not files:
        return None
    times, rcs, cbhs = [], [], []
    rng = wl = lat = lon = alt = None
    for fp in files:
        ds = netCDF4.Dataset(fp)
        try:
            days = np.asarray(ds.variables["time"][:], dtype="float64")
            t = np.datetime64("1970-01-01") + (days * 86400 * 1e9).astype("timedelta64[ns]")
            r = np.asarray(ds.variables["rcs_0"][:], dtype="float32")  # time x range
            cb = np.asarray(ds.variables["cloud_base_height"][:], dtype="float64")
            cb = np.where(cb > 0, cb, np.nan)
            cbh = np.nanmin(cb, axis=1) if cb.ndim == 2 else cb
            if rng is None:
                rng = np.asarray(ds.variables["range"][:], dtype="float64")
                wl = float(ds.variables["l0_wavelength"][...])
                lat = float(ds.variables["station_latitude"][...])
                lon = float(ds.variables["station_longitude"][...])
                alt = float(ds.variables["station_altitude"][...])
        finally:
            ds.close()
        if r.shape[1] != rng.size:
            continue  # skip a day with a different range grid
        times.append(t)
        rcs.append(r)
        cbhs.append(cbh)
    if not times:
        return None
    time = np.concatenate(times)
    rcs = np.concatenate(rcs, axis=0)
    cbh = np.concatenate(cbhs)
    order = np.argsort(time)
    return dict(time=time[order], rcs=rcs[order], cbh=cbh[order],
                range=rng, wl=wl, lat=lat, lon=lon, alt=alt)


def kalman_map(ident: str, method: str) -> dict:
    """date(YYYYMMDD str) -> Kalman C_L for the given method."""
    fp = FULLCAL / f"{STATION}_{ident}" / f"{STATION}_{ident}_kalman.csv"
    out = {}
    if not fp.is_file():
        return out
    import csv
    with open(fp, newline="") as f:
        for row in csv.DictReader(f):
            if row.get("method") == method:
                try:
                    out[row["date"]] = float(row["kalman"])
                except (ValueError, KeyError):
                    pass
    return out


def opconst_map(ident: str) -> dict:
    """date(YYYYMMDD) -> operational L2 calibration_constant_0 (median over day)."""
    out = {}
    for fp in sorted(glob.glob(str(L2_DIR / f"L2_{STATION}_{ident}{YEAR}{MONTH:02d}*.nc"))):
        date = fp[-11:-3]  # YYYYMMDD before .nc
        try:
            ds = netCDF4.Dataset(fp)
            try:
                v = np.asarray(ds.variables["calibration_constant_0"][:], dtype="float64")
            finally:
                ds.close()
            v = v[np.isfinite(v) & (v > 0)]
            if v.size:
                out[date] = float(np.median(v))
        except Exception:  # noqa: BLE001
            pass
    return out


def _const_per_profile(time, cmap: dict, fallback: float) -> np.ndarray:
    """Map each profile's date to a daily constant; fill gaps with the series
    median (the Kalman has all days, so gaps are rare), else the default."""
    med = np.median(list(cmap.values())) if cmap else fallback
    dates = np.datetime_as_string(time.astype("datetime64[D]")).astype("U10")
    dates = np.char.replace(dates, "-", "")
    return np.array([cmap.get(d, med if cmap else fallback) for d in dates], dtype="float64")


# ----------------------------------------------------------------------------
# Figure
# ----------------------------------------------------------------------------
def make_figure(instr: str, res, wl: float) -> str:
    have_wv = "ours_wv" in res.obs_interp
    fig = plt.figure(figsize=(16, 9))
    gs = fig.add_gridspec(2, 3)
    t = res.time_cams

    # (0,0:2) observation pcolor (our-calibrated), log10
    ax = fig.add_subplot(gs[0, 0:2])
    obs = res.obs_mean["ours"].T * UNIT  # (n_r, n_cams)
    obs_disp = np.log10(np.clip(obs, 1e-2, None))
    pcm = ax.pcolormesh(t, res.range_mean + 0, obs_disp, cmap="jet", vmin=-2, vmax=2,
                        shading="auto")
    ax.set_ylim(0, 15000); ax.set_ylabel("Range AGL [m]")
    ax.set_title(f"{instr} observation (our C_L), {wl:.0f} nm")
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%d"))
    fig.colorbar(pcm, ax=ax, label=r"log$_{10}\beta_{att}$ [Mm$^{-1}$sr$^{-1}$]")

    # (0,2) median profiles
    ax = fig.add_subplot(gs[0, 2])
    z = res.z_cams
    ax.fill_betweenx(z, res.prof["ours"]["obs_p25"] * UNIT,
                     res.prof["ours"]["obs_p75"] * UNIT, color=COL["ours"], alpha=0.15)
    ax.plot(res.prof["ours"]["obs_med"] * UNIT, z, "-", color=COL["ours"], label="obs (ours)")
    ax.plot(res.prof["op"]["obs_med"] * UNIT, z, "-", color=COL["op"], label="obs (op)")
    if have_wv:
        ax.plot(res.prof["ours_wv"]["obs_med"] * UNIT, z, "--", color=COL["ours_wv"],
                label="obs (ours, WV)")
    ax.plot(res.cams_med * UNIT, z, "-", color=COL["cams"], lw=2, label="CAMS")
    ax.set_ylim(0, 15000); ax.set_xlim(left=0)
    ax.set_xlabel(r"$\beta_{att}$ [Mm$^{-1}$sr$^{-1}$]"); ax.set_ylabel("Altitude ASL [m]")
    ax.set_title("Median profiles"); ax.grid(alpha=0.3); ax.legend(fontsize=8)

    # (1,0:2) CAMS pcolor
    ax = fig.add_subplot(gs[1, 0:2])
    cams_disp = np.log10(np.clip(res.cams_beta * UNIT, 1e-2, None))
    pcm = ax.pcolormesh(t, z, cams_disp, cmap="jet", vmin=-2, vmax=2, shading="auto")
    ax.set_ylim(0, 15000); ax.set_ylabel("Altitude ASL [m]"); ax.set_xlabel("Day of month")
    ax.set_title(f"CAMS aerosol backscatter at {wl:.0f} nm")
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%d"))
    fig.colorbar(pcm, ax=ax, label=r"log$_{10}\beta_{att}$ [Mm$^{-1}$sr$^{-1}$]")

    # (1,2) bias profile (O - B)
    ax = fig.add_subplot(gs[1, 2])
    for src, lab in (("op", "op - CAMS"), ("ours", "ours - CAMS"),
                     ("ours_wv", "ours WV - CAMS")):
        if src not in res.prof:
            continue
        ax.plot(res.prof[src]["bias_med"] * UNIT, z,
                color=COL[src], label=lab,
                ls="--" if src == "ours_wv" else "-")
    ax.fill_betweenx(z, res.prof["ours"]["bias_p25"] * UNIT,
                     res.prof["ours"]["bias_p75"] * UNIT, color=COL["ours"], alpha=0.12)
    ax.axvline(0, color="k", lw=0.8)
    ax.set_ylim(0, 15000); ax.set_xlabel(r"O - B [Mm$^{-1}$sr$^{-1}$]")
    ax.set_ylabel("Altitude ASL [m]"); ax.set_title("Bias profile"); ax.grid(alpha=0.3)
    ax.legend(fontsize=8)

    fig.suptitle(f"OmB - {STATION_NAME} {instr} - {YEAR}-{MONTH:02d}",
                 fontsize=14, fontweight="bold")
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    path = FIG_DIR / f"omb_{instr}.png"
    fig.savefig(path, dpi=150, facecolor="w")
    plt.close(fig)
    return f"omb_payerne/figs/omb_{instr}.png"


# ----------------------------------------------------------------------------
# Main
# ----------------------------------------------------------------------------
def run_stream(ident: str, instr: str, method: str):
    data = load_l1_month(ident)
    if data is None:
        print(f"{instr}: no L1 data"); return None
    c_ours = _const_per_profile(data["time"], kalman_map(ident, method),
                                INSTRUMENT_CAL_DEFAULT.get(instr, 1.0))
    c_op = _const_per_profile(data["time"], opconst_map(ident),
                              INSTRUMENT_CAL_DEFAULT.get(instr, 1.0))
    rcs = data["rcs"].astype("float64")
    beta_ours = rcs / c_ours[:, None]
    beta_op = rcs / c_op[:, None]
    del rcs
    res = compute_omb(
        time=data["time"], range_agl=data["range"],
        beta_sources={"op": beta_op, "ours": beta_ours},
        station_lat=data["lat"], station_lon=data["lon"], station_alt=data["alt"],
        wavelength=data["wl"], cams_file=CAMS_FILE, instrument=instr,
        cloud_base_height=data["cbh"],
    )
    fig = make_figure(instr, res, data["wl"])
    print(f"{instr}: wl={data['wl']:.0f}  "
          f"median bias ours={res.scalar['ours']['median_bias']*UNIT:.3f} "
          f"op={res.scalar['op']['median_bias']*UNIT:.3f} Mm-1sr-1  "
          f"(WV ours={res.scalar.get('ours_wv',{}).get('median_bias',float('nan'))*UNIT:.3f})")
    return dict(instr=instr, wl=data["wl"], fig=fig, res=res,
                c_ours=np.median(c_ours), c_op=np.median(c_op))


def main():
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    rows = []
    for ident, (instr, method) in STREAMS.items():
        r = run_stream(ident, instr, method)
        if r is not None:
            rows.append(r)

    lines = [f"# OmB spot-check - {STATION_NAME} L1 vs CAMS ({YEAR}-{MONTH:02d})\n"]
    lines.append(
        "Per-station Observation-minus-Background against the CAMS aerosol "
        "forecast, comparing the operational L2 constant and our Kalman "
        "best-estimate C_L. 910 nm streams (CL31, CL61) use the Angstrom "
        "interpolation (532/1064) plus the reused water-vapour correction; "
        "CHM15k (1064 nm) is native. Bias = observation - background, in "
        "Mm^-1 sr^-1.\n")
    lines.append(
        "> CAMS is on a 1 deg grid; the nearest point to Payerne sits in an "
        "Alpine cell whose surface is ~1400 m ASL, so the comparison starts "
        "~900 m above the Payerne plateau (490 m). This is inherent to CAMS "
        "resolution over complex terrain, as in the MATLAB reference. Because the "
        "near-surface (where 910 nm water-vapour absorption is largest) is below "
        "this floor, the reported WV-correction shift is a **lower bound** on the "
        "full-column effect.\n")
    lines.append("## Summary (median bias, Mm^-1 sr^-1)\n")
    hdr = ["Instrument", "lambda [nm]", "median C_L ours", "median C_op",
           "bias ours", "bias op", "bias ours (WV)", "RMS ours"]
    lines.append("| " + " | ".join(hdr) + " |")
    lines.append("| " + " | ".join("---" for _ in hdr) + " |")
    for r in rows:
        s = r["res"].scalar
        wvb = s.get("ours_wv", {}).get("median_bias", float("nan")) * UNIT
        lines.append("| " + " | ".join(str(x) for x in [
            r["instr"], f"{r['wl']:.0f}", f"{r['c_ours']:.3e}", f"{r['c_op']:.3e}",
            f"{s['ours']['median_bias']*UNIT:.4f}", f"{s['op']['median_bias']*UNIT:.4f}",
            f"{wvb:.4f}", f"{s['ours']['rms']*UNIT:.4f}"]) + " |")
    lines.append("")
    for r in rows:
        lines.append(f"## {r['instr']} ({r['wl']:.0f} nm)\n")
        lines.append(f"![omb {r['instr']}]({r['fig']})\n")
    REPORT.write_text("\n".join(lines), encoding="utf-8")
    print(f"Wrote report: {REPORT}")


if __name__ == "__main__":
    main()
