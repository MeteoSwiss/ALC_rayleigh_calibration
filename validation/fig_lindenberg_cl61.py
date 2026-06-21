"""
fig_lindenberg_cl61.py — figure for attbsc_validation_technical.md sec 7.6 (Lindenberg cross-source
Cloudnet CL61). Shows the native-grid handicap and its fix: per-method nights-calibrated on the NATIVE
RAW grid vs binned to the L2 grid (30 m x 300 s), plus the binned per-method robust CV. Uses the
already-computed binned run (rayleigh_Lindenberg_CL61.json) and re-runs the same nights on the native
grid (l1_bin_to_l2_grid=False) for the "before".
"""
from __future__ import annotations
import os, sys, json, warnings, glob
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[_v] = "1"
try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass
import logging; logging.disable(logging.CRITICAL); warnings.filterwarnings("ignore")
from pathlib import Path
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO)); sys.path.insert(0, str(REPO / "validation"))
from calibration import calibrate_rayleigh, CalibrationOptions, InstrumentInfo, DataLevel
from calibration.config import InstrumentType
from compare_molecular_methods import run_methods, calibrates

OUT = Path("C:/DATA/Projects/202606_E-PROFILE_calibration/figs_paper_validation/cloudnet_cl61")
ROOT = Path("A:/CL61_Cloudnet")
METHODS = ["eprof_v1.1", "eprof_v1.2", "eprof_v0.25", "eprof_v2", "earlinet", "bellini"]
LAB = {"eprof_v1.1": "v1.1", "eprof_v1.2": "v1.2", "eprof_v0.25": "v0.25", "eprof_v2": "v2 (C8)",
       "earlinet": "EARLINET", "bellini": "Bellini"}
COL = {"eprof_v2": "#d62728"}


def robcv(x):
    x = np.asarray(x, float)
    if x.size < 2:
        return np.nan
    m = np.median(x)
    return 1.4826 * np.median(np.abs(x - m)) / abs(m) * 100 if m != 0 else np.nan


def base_options(bin_to_l2):
    o = CalibrationOptions.from_json(REPO / "options.json")
    o.folder_root = ROOT; o.data_level = DataLevel.RAW; o.molecular_source = "standard"
    o.apply_wv_correction = True; o.plot_main = False; o.plot_all = False
    o.l1_bin_to_l2_grid = bin_to_l2; o.folder_output = OUT
    return o


def main():
    binned = json.loads((OUT / "rayleigh_Lindenberg_CL61.json").read_text())
    nights = sorted(binned)
    info = InstrumentInfo(site_name="Lindenberg_CL61", wmo_id="Lindenberg", identifier="",
                          instrument_type=InstrumentType.CL61, latitude=52.21, longitude=14.12, altitude=123.0)
    # native re-run (same nights) for the "before"
    native = {}
    o = base_options(False)
    for k, ds in enumerate(nights):
        fin = {}
        try:
            calibrate_rayleigh(ds, info, o, fit_inputs_out=fin)
            if fin:
                res = run_methods(fin["signal"], fin["p_mol"], fin["range_alc"], fin["signal_stack"])
                native[ds] = {m: [bool(res[m].ok), float(res[m].cl)] for m in METHODS}
        except Exception:
            pass
        if (k + 1) % 10 == 0:
            print(f"  native re-run {k+1}/{len(nights)}", flush=True)
    (OUT / "rayleigh_Lindenberg_CL61_native.json").write_text(json.dumps(native), encoding="utf-8")

    n_native = {m: sum(1 for ds in native if native[ds][m][0]) for m in METHODS}
    n_binned = {m: sum(1 for ds in binned if binned[ds][m][0]) for m in METHODS}
    cv_binned = {m: robcv([binned[ds][m][1] for ds in binned
                           if binned[ds][m][0] and np.isfinite(binned[ds][m][1]) and binned[ds][m][1] > 0])
                 for m in METHODS}
    Nn, Nb = len(native), len(binned)

    fig, (axA, axB) = plt.subplots(1, 2, figsize=(17, 6))
    x = np.arange(len(METHODS)); w = 0.38
    axA.bar(x - w/2, [n_native[m] for m in METHODS], w, label=f"native RAW (4.8 m × ~60 s, n={Nn})", color="#aaaaaa")
    axA.bar(x + w/2, [n_binned[m] for m in METHODS], w, label=f"binned to L2 grid (30 m × 300 s, n={Nb})", color="#1f77b4")
    for i, m in enumerate(METHODS):
        axA.annotate(str(n_native[m]), (i - w/2, n_native[m]), ha="center", va="bottom", fontsize=8)
        axA.annotate(str(n_binned[m]), (i + w/2, n_binned[m]), ha="center", va="bottom", fontsize=8,
                     fontweight="bold" if m == "eprof_v2" else "normal",
                     color="#d62728" if m == "eprof_v2" else "black")
    axA.set_xticks(x); axA.set_xticklabels([LAB[m] for m in METHODS])
    axA.set_ylabel(f"nights calibrated (of {Nb})")
    axA.set_title("Native-grid handicap and its fix: nights calibrated per method\n"
                  f"(gated methods recover once native RAW is binned to the L2 grid: "
                  f"v2 {n_native['eprof_v2']}→{n_binned['eprof_v2']}, Bellini {n_native['bellini']}→{n_binned['bellini']})",
                  fontsize=11)
    axA.legend(fontsize=9); axA.grid(axis="y", alpha=0.3)

    bars = axB.bar(x, [cv_binned[m] for m in METHODS], color=["#d62728" if m == "eprof_v2" else "#1f77b4" for m in METHODS])
    for i, m in enumerate(METHODS):
        v = cv_binned[m]
        if np.isfinite(v):
            axB.annotate(f"{v:.0f}%", (i, v), ha="center", va="bottom", fontsize=8)
    axB.set_xticks(x); axB.set_xticklabels([LAB[m] for m in METHODS])
    axB.set_ylabel("robust CV of lidar constant C_L (%)")
    axB.set_title("Night-to-night stability (binned to the L2 grid; lower = better)\n"
                  "v2 = best precision–yield balance of the high-yield methods", fontsize=11)
    axB.grid(axis="y", alpha=0.3)

    fig.suptitle("Lindenberg cross-source Cloudnet CL61 — molecular methods after the native-grid fix", fontsize=14)
    fig.tight_layout(); fig.savefig(OUT / "lindenberg_cl61_methods.png", dpi=130); plt.close(fig)
    print("native n_cal:", n_native)
    print("binned n_cal:", n_binned)
    print("saved lindenberg_cl61_methods.png")


if __name__ == "__main__":
    main()
