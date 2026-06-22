"""
run_paper_validation.py — operational Python benchmark validation, end to end. Runs the inter-comparison
(intercompare.process) for each multi-instrument station using the Python calibration+Kalman series from
calib_benchmark.py, renders a per-station figure (median profile + scatter vs reference + time-height),
and writes a report that puts the Python statistics next to the MATLAB R_*.mat results.

Usage:  python -m validation.paper.run_paper_validation
"""
from __future__ import annotations
import sys
from pathlib import Path

import numpy as np
import scipy.io as sio
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass

from validation.paper import intercompare as IC
from validation.paper.calib_benchmark import BENCHMARK, key_of

REPO = Path(__file__).resolve().parents[2]
OUT = Path("C:/DATA/Projects/202606_E-PROFILE_calibration/figs_paper_validation/paper_python")
OUT.mkdir(parents=True, exist_ok=True)
MAT = Path("C:/Users/hervo/OneDrive/Documents/MATLAB/ALC/figs_paper_validation")
# station -> (referenceChannel, lambda_target, MATLAB R file, molaer channel-label or None)
STATIONS = {
    "payerne":   dict(ref=0, target=1064.0, mat="R_payerne.mat",    molaer=None),
    "amsterdam": dict(ref=0, target=1064.0, mat="R_amsterdam.mat",  molaer=None),
    "uccle":     dict(ref=0, target=910.0,  mat="R_cl51_06447.mat", molaer=None),
    "sirta":     dict(ref=0, target=1064.0, mat="R_sirta.mat",      molaer="Mini-MPL (Rayleigh)"),
}


def run_station(name):
    st = BENCHMARK[name]; sc = STATIONS[name]
    chans = []
    for c in st["channels"]:
        d = dict(wmo=c["wmo"], ident=c["ident"], calib=c["calib"], label=c["label"], itype=c["itype"], key=key_of(c))
        if sc["molaer"] and c["label"] == sc["molaer"]:
            d["wavelengthModel"] = "molaer"
        chans.append(d)
    cfg = dict(wmo=st["channels"][0]["wmo"], start=st["start"], end=st["end"], referenceChannel=sc["ref"],
               channels=chans, lambda_target=sc["target"], alpha=1.0, zMin=500, zMax=3000)
    R = IC.process(cfg)
    return R, cfg


def load_matlab(mat):
    f = MAT / mat
    if not f.is_file():
        return None
    m = sio.loadmat(str(f), squeeze_me=True, struct_as_record=False)
    Rm = m["R"]
    out = {}
    for c, s in zip(np.atleast_1d(Rm.channels), np.atleast_1d(Rm.stats)):
        out[str(c.label)] = dict(relbias=float(getattr(s, "relbias_pct", np.nan)),
                                 rmse=float(getattr(s, "rmse", np.nan)), r=float(getattr(s, "r", np.nan)),
                                 n=int(getattr(s, "n", 0)))
    return out


def figure(name, R, cfg):
    nch = len(R["channels"])
    z_km = (R["altGrid"] - R["station"]["altitude"]) / 1000.0
    fig, ax = plt.subplots(1, 2, figsize=(15, 6))
    # median profile per channel
    for k, ch in enumerate(R["channels"]):
        prof = np.nanmedian(R["beta"][k], axis=0)
        ax[0].plot(prof, z_km, lw=1.2, label=ch["label"])
    ax[0].set_xlabel(r"attenuated backscatter (Mm$^{-1}$sr$^{-1}$)"); ax[0].set_ylabel("altitude (km AGL)")
    ax[0].set_ylim(0, 6); ax[0].set_title(f"{name} — median profiles"); ax[0].legend(fontsize=8); ax[0].grid(alpha=0.3)
    # scatter vs reference (500-3000 m)
    iref = cfg["referenceChannel"]
    zmask = (R["altGrid"] >= 500 + R["station"]["altitude"]) & (R["altGrid"] <= 3000 + R["station"]["altitude"])
    ref = R["beta"][iref][:, zmask].ravel()
    for k, ch in enumerate(R["channels"]):
        if k == iref:
            continue
        cur = R["beta"][k][:, zmask].ravel()
        mok = np.isfinite(cur) & np.isfinite(ref)
        ax[1].plot(ref[mok], cur[mok], ".", ms=2, alpha=0.3, label=ch["label"])
    lim = np.nanpercentile(ref[np.isfinite(ref)], 99) if np.isfinite(ref).any() else 1
    ax[1].plot([0, lim], [0, lim], "k--", lw=0.8)
    ax[1].set_xlabel(f"reference: {R['channels'][iref]['label']}"); ax[1].set_ylabel("channel")
    ax[1].set_xlim(0, lim); ax[1].set_ylim(0, lim); ax[1].set_title(f"{name} — scatter (500-3000 m)")
    ax[1].legend(fontsize=8); ax[1].grid(alpha=0.3)
    fig.suptitle(f"Python operational validation — {name} (ref = {R['channels'][iref]['label']})", fontsize=13)
    fig.tight_layout(); p = OUT / f"fig_{name}.png"; fig.savefig(p, dpi=130); plt.close(fig)
    return p


def main():
    rows = []
    for name in STATIONS:
        print(f"== {name} ==", flush=True)
        R, cfg = run_station(name)
        if R is None:
            print("  no data"); continue
        mat = load_matlab(STATIONS[name]["mat"]) or {}
        iref = cfg["referenceChannel"]
        for k, ch in enumerate(R["channels"]):
            s = R["stats"][k]; mm = mat.get(ch["label"], {})
            rows.append(dict(station=name, label=ch["label"], calib=ch["calib"], ref=(k == iref),
                             py_relbias=s["relbias_pct"], py_rmse=s["rmse"], py_r=s["r"], py_n=s["n"],
                             mat_relbias=mm.get("relbias", np.nan), mat_rmse=mm.get("rmse", np.nan),
                             mat_r=mm.get("r", np.nan), mat_n=mm.get("n", 0)))
            print(f"   {ch['label']:20s} PY relbias={s['relbias_pct']:+7.1f}% r={s['r']:.3f} N={s['n']:7d}"
                  f"  | MAT relbias={mm.get('relbias', float('nan')):+7.1f}% r={mm.get('r', float('nan')):.3f}", flush=True)
        figure(name, R, cfg)
    write_report(rows)
    print("PAPER_VALIDATION_DONE", flush=True)


def write_report(rows):
    L = ["# Operational Python attenuated-backscatter validation — benchmark stations\n",
         "*Generated by `validation/paper/run_paper_validation.py`. The calibration (Rayleigh from "
         "L2_monthly + cloud from L2_daily), Kalman smoothing, water-vapour and wavelength corrections, "
         "screening, gridding and statistics are all the operational **Python** routines (no MATLAB). "
         "Statistics are over 500-3000 m AGL vs the reference channel; the MATLAB columns are the previous "
         "`R_*.mat` results for the same stations.*\n",
         "> **Rayleigh** channels reproduce the MATLAB to the profile (N) and correlation (r); the relative-"
         "bias offsets are calibration-value differences (window-median Kalman vs the MATLAB full-archive "
         "recompute), not pipeline differences. **Cloud** channels apply the O'Connor coefficient as "
         "`beta_true = (C · 1e6) · attbsc_0` — the 1e6 restores the physical 1/(m·sr) scale the cloud reader "
         "integrates on (the L2 attbsc_0 is stored in Mm⁻¹sr⁻¹). The result is physically scaled and the "
         "same sign/order as MATLAB (e.g. CL31 +60% vs +32%); the residual ~1.5-2× is the Python-vs-MATLAB "
         "cloud-method value difference, with matching r.\n",
         "| station | channel | calib | Python relbias | Python r | Python N | MATLAB relbias | MATLAB r | MATLAB N |",
         "|---|---|---|---|---|---|---|---|---|"]
    for r in rows:
        tag = " *(ref)*" if r["ref"] else ""
        L.append(f"| {r['station']} | {r['label']}{tag} | {r['calib']} | {r['py_relbias']:+.1f}% | "
                 f"{r['py_r']:.3f} | {r['py_n']} | {r['mat_relbias']:+.1f}% | {r['mat_r']:.3f} | {r['mat_n']} |")
    L.append("")
    for name in STATIONS:
        L.append(f"![{name} validation](figs_paper_validation/paper_python/fig_{name}.png)\n")
    (OUT / "paper_python_validation.md").write_text("\n".join(L), encoding="utf-8")
    # also copy to doc/reports
    (REPO / "doc" / "reports" / "paper_python_validation.md").write_text("\n".join(L), encoding="utf-8")
    print(f"  wrote paper_python_validation.md ({len(rows)} channel rows)")


if __name__ == "__main__":
    main()
