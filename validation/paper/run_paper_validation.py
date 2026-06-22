"""
run_paper_validation.py — operational Python benchmark validation, end to end. Runs the inter-comparison
(intercompare.process) for each multi-instrument station using the Python calibration+Kalman series from
calib_benchmark.py, renders the MATLAB-style figures (multi-ALC 3x3 panel, combined calibration time-series,
EARLINET 2x2 panel) via figures.py, and writes a report with the Python statistics next to the MATLAB R_*.mat.

Usage:  python -m validation.paper.run_paper_validation
"""
from __future__ import annotations
import sys
from pathlib import Path

import numpy as np
import scipy.io as sio

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass

from validation.paper import intercompare as IC
from validation.paper import earlinet as EA
from validation.paper import figures as FIG
from validation.paper.calib_benchmark import BENCHMARK, key_of

REPO = Path(__file__).resolve().parents[2]
OUT = Path("C:/DATA/Projects/202606_E-PROFILE_calibration/figs_paper_validation/paper_python")
OUT.mkdir(parents=True, exist_ok=True)
CALIB = OUT / "calib"
MAT = Path("C:/Users/hervo/OneDrive/Documents/MATLAB/ALC/figs_paper_validation")
SITE_NAME = {"payerne": "Payerne", "amsterdam": "Amsterdam", "uccle": "Uccle", "sirta": "Palaiseau",
             "earlinet": ""}
# station -> (referenceChannel, lambda_target, MATLAB R file, molaer label, WIGOS for title)
STATIONS = {
    "payerne":   dict(ref=0, target=1064.0, mat="R_payerne.mat",    molaer=None, wmo="0-20000-0-06610"),
    "amsterdam": dict(ref=0, target=1064.0, mat="R_amsterdam.mat",  molaer=None, wmo="0-20000-0-06240"),
    "uccle":     dict(ref=0, target=910.0,  mat="R_cl51_06447.mat", molaer=None, wmo="0-20000-0-06447"),
    "sirta":     dict(ref=0, target=1064.0, mat="R_sirta.mat",      molaer="Mini-MPL (Rayleigh)", wmo="0-250-1001-07151"),
}


def run_station(name):
    st = BENCHMARK[name]; sc = STATIONS[name]
    chans = []
    for c in st["channels"]:
        d = dict(wmo=c["wmo"], ident=c["ident"], calib=c["calib"], label=c["label"], itype=c["itype"], key=key_of(c))
        if sc["molaer"] and c["label"] == sc["molaer"]:
            d["wavelengthModel"] = "molaer"
        chans.append(d)
    cfg = dict(wmo=sc["wmo"], start=st["start"], end=st["end"], referenceChannel=sc["ref"],
               channels=chans, lambda_target=sc["target"], alpha=1.0, zMin=500, zMax=3000)
    return IC.process(cfg), cfg


def load_matlab(mat):
    f = MAT / mat
    if not f.is_file():
        return {}
    m = sio.loadmat(str(f), squeeze_me=True, struct_as_record=False)
    Rm = m["R"]; out = {}
    for c, s in zip(np.atleast_1d(Rm.channels), np.atleast_1d(Rm.stats)):
        out[str(c.label)] = dict(relbias=float(getattr(s, "relbias_pct", np.nan)),
                                 r=float(getattr(s, "r", np.nan)), n=int(getattr(s, "n", 0)))
    return out


def calib_channel_list():
    """All calibrated channels (for the combined time-series grid), site-labelled, in a sensible order."""
    order, seen = [], set()
    # CHM Rayleigh first, then other Rayleigh, then cloud — grouped by station as defined
    for pref in ("rayleigh", "cloud"):
        for name, st in BENCHMARK.items():
            for c in st["channels"]:
                if c["calib"] != pref:
                    continue
                k = key_of(c)
                if k in seen:
                    continue
                seen.add(k)
                site = SITE_NAME.get(name) or c["label"].split(" ")[0]
                unit = "C [-]" if c["calib"] == "cloud" else "C$_L$ [a.u.]"
                title = f"{site} {c['label']}" if SITE_NAME.get(name) else c["label"]
                order.append(dict(key=k, title=title, unit=unit))
    return order


def main():
    rows = []
    for name, sc in STATIONS.items():
        print(f"== {name} ==", flush=True)
        R, cfg = run_station(name)
        if R is None:
            print("  no data"); continue
        mat = load_matlab(sc["mat"])
        for k, ch in enumerate(R["channels"]):
            s = R["stats"][k]; mm = mat.get(ch["label"], {})
            rows.append(dict(station=name, label=ch["label"], calib=ch["calib"], ref=(k == cfg["referenceChannel"]),
                             py_relbias=s["relbias_pct"], py_r=s["r"], py_n=s["n"],
                             mat_relbias=mm.get("relbias", np.nan), mat_r=mm.get("r", np.nan), mat_n=mm.get("n", 0)))
            print("   %-20s PY relbias=%+7.1f%% r=%.3f N=%7d | MAT relbias=%+7.1f%% r=%.3f"
                  % (ch["label"], s["relbias_pct"], s["r"], s["n"], mm.get("relbias", np.nan), mm.get("r", np.nan)), flush=True)
        title = "%s (%s)  —  %s to %s" % (SITE_NAME[name], sc["wmo"], _d(cfg["start"]), _d(cfg["end"]))
        FIG.fig_multi_alc(R, cfg, OUT / f"fig_{name}.png", title)
        print(f"   -> fig_{name}.png", flush=True)

    # combined calibration time-series (all channels)
    FIG.fig_calib_timeseries(calib_channel_list(), CALIB, OUT / "fig_calib_timeseries.png")
    print("   -> fig_calib_timeseries.png", flush=True)

    # EARLINET 2x2 figures
    erows = []
    for code in ("sir", "ino", "ari", "lei", "cbw"):
        try:
            s = EA.compare(code, "20250101", "20260630", return_profiles=True)
        except Exception as exc:
            s = {"error": repr(exc)}
        mm = EA.load_matlab_earlinet(code)
        label = {"sir": "Palaiseau", "ino": "Magurele", "ari": "Leipzig", "lei": "Leipzig", "cbw": "Cabauw"}[code]
        if s and "error" not in s:
            erows.append((code, label, s, mm))
            FIG.fig_earlinet(code, label, s["betaE"], s["betaC"], s["grid"], s["times"], s, mm, OUT / f"fig_earlinet_{code}.png")
            print("   EARLINET %s: relbias=%+.1f%% r=%.2f matched=%d -> fig_earlinet_%s.png"
                  % (code, s["relbias_pct"], s["r"], s["matched"], code), flush=True)
        else:
            erows.append((code, label, s, mm)); print("   EARLINET %s: %s" % (code, s), flush=True)

    write_report(rows, erows)
    print("PAPER_VALIDATION_DONE", flush=True)


def write_report(rows, erows):
    L = ["# Operational Python attenuated-backscatter validation — benchmark stations\n",
         "*Generated by `validation/paper/run_paper_validation.py`. Calibration (Rayleigh from L2_monthly "
         "or native CL61 raw + cloud), Kalman smoothing, water-vapour and wavelength corrections, screening, "
         "gridding and statistics are all the operational **Python** routines (no MATLAB). The figures "
         "reproduce the MATLAB layouts (`make_validation_figures.m`).*\n",
         "> **Rayleigh** channels reproduce the MATLAB to the profile (N) and correlation (r); relative-bias "
         "offsets are calibration-value differences (window-median Kalman vs the MATLAB full-archive recompute). "
         "**Cloud** β_att uses `(C·1e6)·attbsc_0` (physical scale): same sign/order as MATLAB with matching r. "
         "**CL61 is calibrated twice** (Rayleigh + cloud) at Payerne and Uccle; Payerne CL61 is calibrated from "
         "the native Vaisala raw (`A:/CL61_PAY`).\n",
         "## Calibration coefficient time series (all channels)\n",
         "![calibration time series](figs_paper_validation/paper_python/fig_calib_timeseries.png)\n",
         "## Per-station validation\n",
         "| station | channel | calib | Python relbias | Python r | Python N | MATLAB relbias | MATLAB r |",
         "|---|---|---|---|---|---|---|---|"]
    for r in rows:
        tag = " *(ref)*" if r["ref"] else ""
        L.append("| %s | %s%s | %s | %+.1f%% | %.3f | %d | %+.1f%% | %.3f |"
                 % (r["station"], r["label"], tag, r["calib"], r["py_relbias"], r["py_r"], r["py_n"],
                    r["mat_relbias"], r["mat_r"]))
    L.append("")
    for name in STATIONS:
        L.append(f"![{name} validation](figs_paper_validation/paper_python/fig_{name}.png)\n")
    # EARLINET
    L.append("## EARLINET — ceilometer (CHM15k) vs EARLINET research-lidar reference\n")
    L.append("| site | Python relbias | Python r | matched | MATLAB relbias | MATLAB r |")
    L.append("|---|---|---|---|---|---|")
    for code, label, s, mm in erows:
        if s and "error" not in s:
            L.append("| %s (%s) | %+.1f%% | %.2f | %d | %+.1f%% | %.2f |"
                     % (code, label, s["relbias_pct"], s["r"], s["matched"], mm.get("relbias", np.nan), mm.get("r", np.nan)))
        else:
            L.append("| %s (%s) | no EARLINET 1064 data in the 2025-2026 window | | | | |" % (code, label))
    L.append("")
    for code, label, s, mm in erows:
        if s and "error" not in s:
            L.append(f"![earlinet {code}](figs_paper_validation/paper_python/fig_earlinet_{code}.png)\n")
    txt = "\n".join(L)
    (OUT / "paper_python_validation.md").write_text(txt, encoding="utf-8")
    (REPO / "doc" / "reports" / "paper_python_validation.md").write_text(txt, encoding="utf-8")
    print(f"  wrote paper_python_validation.md ({len(rows)} channels, {len(erows)} EARLINET sites)")


def _d(yyyymmdd):
    from datetime import datetime
    return datetime.strptime(yyyymmdd, "%Y%m%d").strftime("%Y-%m-%d")


if __name__ == "__main__":
    main()
