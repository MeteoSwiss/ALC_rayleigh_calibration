#!/usr/bin/env python3
"""
Compare the MATLAB Rayleigh calibration (A:/E-PROFILE_L2_Calibration/rayleigh_per_station,
one rayleigh_<WMO>_<id>.mat per station: daily_C, daily_C_std, daily_C_kalman, meta)
with the Python calibration (C:/DATA/Projects/202606_E-PROFILE_calibration/E-PROFILE_calibration_rayleigh/fullcal_all/<WMO>/<WMO>_cl.csv).

Both pipelines reconstruct rcs from L2-monthly the same way
(rcs = attenuated_backscatter_0 x 1e-6 x calibration_constant_0 in MATLAB;
Python uses an auto-detected unit factor), so the dimensionless
**calibration coefficient = re-derived C_L / operational C_op** (ideal = 1.0)
is computed for both and shown as ranked-bar figures in the style of
"All Stations Ranked by Calibration Coefficient".

Outputs (C:/DATA/Projects/202606_E-PROFILE_calibration/E-PROFILE_calibration_rayleigh/comparison_matlab_python/):
  ranked_coeff_matlab.png    ranked bars, MATLAB method
  ranked_coeff_python.png    ranked bars, Python method
  ranked_coeff_both.png      two-panel combined figure
  scatter_common.png         MATLAB vs Python coefficient, common stations
  comparison_table.csv       per-station numbers for both methods
"""
import csv
import glob
import json
import os
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import scipy.io as sio

MAT_DIR = Path("A:/E-PROFILE_L2_Calibration/rayleigh_per_station")
PY_DIR = Path("C:/DATA/Projects/202606_E-PROFILE_calibration/E-PROFILE_calibration_rayleigh/fullcal_all")
OUT = Path("C:/DATA/Projects/202606_E-PROFILE_calibration/E-PROFILE_calibration_rayleigh/comparison_matlab_python")
COP = json.load(open("C:/DATA/Projects/202606_E-PROFILE_calibration/cop_lookup.json"))
# Keyed per instrument '<WMO>_<identifier>' (a WMO may host several instruments).
MANIFEST = {f"{s['wmo']}_{s['identifier']}": s for s in json.load(open("C:/DATA/Projects/202606_E-PROFILE_calibration/stations_l2_manifest.json"))}

TYPE_COLORS = {  # keep the example figure's colors for the Vaisala family
    "CL31": "#1f77b4", "CL51": "#d2691e", "CL61": "#e6c229",
    "CHM15k": "#2ca02c", "Mini-MPL": "#9467bd",
}
COEF_CLIP = 12.0     # display clip: coefficients beyond this are off-scale/unit issues
SUCCESS = {"1", "0.5", "1.0"}


def cop_for(key):
    e = COP.get(key)
    return (e["cop_median"] if e else None)


# --------------------------------------------------------------------------- #
def _read_mat(f):
    """Read daily_C + instrument_type from a station .mat (v7 via scipy, v7.3 via h5py)."""
    try:
        m = sio.loadmat(f, squeeze_me=True, struct_as_record=False)
        if "daily_C" not in m:
            return None, "?"
        c = np.atleast_1d(np.asarray(m["daily_C"], dtype="f8"))
        itype = "?"
        try:
            itype = str(m["meta"].instrument_type)
        except Exception:
            pass
        return c, itype
    except NotImplementedError:
        # v7.3 (HDF5) fallback
        import h5py
        with h5py.File(f, "r") as h:
            if "daily_C" not in h:
                return None, "?"
            c = np.atleast_1d(np.asarray(h["daily_C"], dtype="f8")).ravel()
            itype = "?"
            try:
                ref = h["meta"]["instrument_type"]
                itype = "".join(chr(int(x)) for x in np.asarray(ref).ravel())
            except Exception:
                pass
            return c, itype
    except Exception:
        return None, "?"


def load_matlab():
    """Per-station: median coefficient + IQR from daily_C / C_op."""
    rows = []
    for f in sorted(MAT_DIR.glob("rayleigh_*.mat")):
        # stem is 'rayleigh_<WMO>_<id>' -> instrument key '<WMO>_<id>'
        key = f.stem.split("_", 1)[1]
        c, itype = _read_mat(f)
        if c is None:
            continue
        c = c[np.isfinite(c) & (c > 0)]
        if c.size < 3:
            continue
        cop = cop_for(key)
        if not cop:
            continue
        if itype in ("?", ""):
            itype = MANIFEST.get(key, {}).get("itype", "?")
        coeff = c / cop
        rows.append(dict(key=key, itype=itype, n=c.size,
                         med=float(np.median(coeff)),
                         q25=float(np.percentile(coeff, 25)),
                         q75=float(np.percentile(coeff, 75))))
    return rows


def load_python():
    """Per-station: median coefficient + IQR from successful nights / C_op."""
    rows = []
    for f in sorted(PY_DIR.glob("*/*_cl.csv")):
        key = f.parent.name  # '<WMO>_<identifier>'
        cls = []
        for r in csv.DictReader(open(f)):
            if r["flag"] in SUCCESS:
                try:
                    v = float(r["lidar_constant"])
                    if v > 0:
                        cls.append(v)
                except ValueError:
                    pass
        if len(cls) < 3:
            continue
        cop = cop_for(key)
        if not cop:
            continue
        coeff = np.asarray(cls) / cop
        itype = MANIFEST.get(key, {}).get("itype", "?")
        rows.append(dict(key=key, itype=itype, n=len(cls),
                         med=float(np.median(coeff)),
                         q25=float(np.percentile(coeff, 25)),
                         q75=float(np.percentile(coeff, 75))))
    return rows


# --------------------------------------------------------------------------- #
def ranked_panel(ax, rows, title):
    rows = [r for r in rows if np.isfinite(r["med"])]
    shown = [r for r in rows if r["med"] <= COEF_CLIP]
    n_clip = len(rows) - len(shown)
    shown.sort(key=lambda r: r["med"])
    x = np.arange(len(shown))
    meds = np.array([r["med"] for r in shown])
    lo = meds - np.array([r["q25"] for r in shown])
    hi = np.array([r["q75"] for r in shown]) - meds
    colors = [TYPE_COLORS.get(r["itype"], "0.5") for r in shown]
    ax.bar(x, meds, width=0.85, color=colors, edgecolor="none")
    # clip error bars to the display window so a few huge-IQR stations don't blow the axis
    lo = np.clip(lo, 0, meds)
    hi = np.clip(hi, 0, COEF_CLIP - meds)
    ax.errorbar(x, meds, yerr=[lo, hi],
                fmt="none", ecolor="k", elinewidth=0.6, capsize=1.5)
    ax.set_ylim(0, COEF_CLIP)
    net_med = float(np.median(meds)) if len(meds) else np.nan
    ax.axhline(1.0, color="red", ls="--", lw=1.6, label="Ideal (1.0)")
    ax.axhline(net_med, color="blue", ls="--", lw=1.6,
               label=f"Network Median ({net_med:.3f})")
    handles = [Patch(color=c, label=t) for t, c in TYPE_COLORS.items()
               if any(r["itype"] == t for r in shown)]
    handles += [plt.Line2D([], [], color="red", ls="--", label="Ideal (1.0)"),
                plt.Line2D([], [], color="blue", ls="--",
                           label=f"Network Median ({net_med:.3f})")]
    ax.legend(handles=handles, loc="upper left", fontsize=9)
    extra = f" — {n_clip} off-scale station(s) not shown" if n_clip else ""
    ax.set_title(f"{title}  ({len(shown)} stations{extra})")
    ax.set_xlabel("Station (sorted by calibration coefficient)")
    ax.set_ylabel("Calibration Coefficient")
    ax.set_xlim(-1, len(shown))
    ax.grid(True, axis="y", alpha=0.3)
    return net_med


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    mat = load_matlab()
    py = load_python()
    print(f"MATLAB stations: {len(mat)} | Python stations: {len(py)}")

    # individual figures
    for rows, name, label in [(mat, "ranked_coeff_matlab", "MATLAB method"),
                              (py, "ranked_coeff_python", "Python method")]:
        fig, ax = plt.subplots(figsize=(16, 6))
        ranked_panel(ax, rows, f"All Stations Ranked by Calibration Coefficient — {label}")
        fig.tight_layout()
        fig.savefig(OUT / f"{name}.png", dpi=140, bbox_inches="tight")
        plt.close(fig)

    # combined two-panel
    fig, axes = plt.subplots(2, 1, figsize=(16, 11))
    ranked_panel(axes[0], mat, "MATLAB method (rayleigh_per_station)")
    ranked_panel(axes[1], py, "Python method (fullcal_all, in progress)")
    fig.suptitle("All Stations Ranked by Calibration Coefficient — re-derived C$_L$ / operational C$_{op}$",
                 fontsize=14)
    fig.tight_layout(rect=[0, 0, 1, 0.97])
    fig.savefig(OUT / "ranked_coeff_both.png", dpi=140, bbox_inches="tight")
    plt.close(fig)

    # scatter on common instruments
    md = {r["key"]: r for r in mat}
    pd_ = {r["key"]: r for r in py}
    common = sorted(set(md) & set(pd_))
    fig, ax = plt.subplots(figsize=(7.5, 7))
    for w in common:
        a, b = md[w], pd_[w]
        if a["med"] <= COEF_CLIP and b["med"] <= COEF_CLIP:
            ax.plot(a["med"], b["med"], "o", ms=5,
                    color=TYPE_COLORS.get(b["itype"], "0.5"), alpha=0.8)
    lim = [0, COEF_CLIP / 2]
    ax.plot(lim, lim, "k--", lw=1)
    ax.set_xlim(lim); ax.set_ylim(lim)
    ax.set_xlabel("MATLAB coefficient (median)")
    ax.set_ylabel("Python coefficient (median)")
    ax.set_title(f"Common stations ({len(common)}) — MATLAB vs Python")
    ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(OUT / "scatter_common.png", dpi=140, bbox_inches="tight")
    plt.close(fig)

    # table + console summary
    with open(OUT / "comparison_table.csv", "w", newline="") as fcsv:
        w = csv.writer(fcsv)
        w.writerow(["key", "itype", "matlab_med", "matlab_n", "python_med", "python_n",
                    "rel_diff_pct"])
        for key in sorted(set(md) | set(pd_)):
            a, b = md.get(key), pd_.get(key)
            rd = ""
            if a and b and a["med"] > 0:
                rd = f"{(b['med'] - a['med']) / a['med'] * 100:.2f}"
            w.writerow([key,
                        (b or a)["itype"], f"{a['med']:.4f}" if a else "", a["n"] if a else "",
                        f"{b['med']:.4f}" if b else "", b["n"] if b else "", rd])
    if common:
        rel = [(pd_[w]["med"] - md[w]["med"]) / md[w]["med"] * 100 for w in common
               if md[w]["med"] > 0 and md[w]["med"] < COEF_CLIP and pd_[w]["med"] < COEF_CLIP]
        print(f"common stations: {len(common)} | median rel diff (Py-Mat): "
              f"{np.median(rel):+.2f}% | IQR [{np.percentile(rel,25):+.1f}, {np.percentile(rel,75):+.1f}]%")
    print(f"outputs -> {OUT}")


if __name__ == "__main__":
    main()
