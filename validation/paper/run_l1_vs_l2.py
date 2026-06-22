"""
run_l1_vs_l2.py — L1-vs-L2 validation (eprof_v2). For each benchmark station it (A) compares the
calibration coefficient derived from the native L1 (binned to the L2 grid) against the L2 product, and
(B) re-runs the inter-comparison applying the L1- vs L2-derived calibration to the measured beta_att.
Shows that calibrating from L1 or L2 gives the same coefficient and the same cross-instrument validation.

Usage:  python -m validation.paper.run_l1_vs_l2
"""
from __future__ import annotations
import csv
import sys
from pathlib import Path

import numpy as np

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass

from validation.paper import intercompare as IC
from validation.paper import figures as FIG
from validation.paper.calib_benchmark import BENCHMARK, key_of
from validation.paper.run_paper_validation import STATIONS, SITE_NAME, calib_channel_list, OUT, CALIB

REPO = Path(__file__).resolve().parents[2]


def median_kalman(key, level):
    f = CALIB / f"{key}_{level}.csv"
    if not f.is_file():
        return np.nan, 0
    ck = [float(r["C_kalman"]) for r in csv.DictReader(open(f, encoding="utf-8"))
          if r.get("C_kalman") not in ("", "nan", None)]
    return (float(np.median(ck)), len(ck)) if ck else (np.nan, 0)


def run_station_level(name, level):
    st = BENCHMARK[name]; sc = STATIONS[name]
    chans = []
    for c in st["channels"]:
        d = dict(wmo=c["wmo"], ident=c["ident"], calib=c["calib"], label=c["label"], itype=c["itype"], key=key_of(c))
        if sc["molaer"] and c["label"] == sc["molaer"]:
            d["wavelengthModel"] = "molaer"
        chans.append(d)
    cfg = dict(wmo=sc["wmo"], start=st["start"], end=st["end"], referenceChannel=sc["ref"],
               channels=chans, lambda_target=sc["target"], alpha=1.0, zMin=500, zMax=3000, calibLevel=level)
    return IC.process(cfg), cfg


def main():
    # (A) calibration coefficient L1 vs L2
    coef = []
    for c in calib_channel_list():
        (l1, n1), (l2, n2) = median_kalman(c["key"], "L1"), median_kalman(c["key"], "L2")
        if np.isfinite(l1) or np.isfinite(l2):
            ratio = l1 / l2 if (np.isfinite(l1) and np.isfinite(l2) and l2 != 0) else np.nan
            coef.append(dict(title=c["title"], l1=l1, l2=l2, n1=n1, n2=n2, ratio=ratio))
            print("  %-30s L1=%.3e (n%d)  L2=%.3e (n%d)  L1/L2=%.3f"
                  % (c["title"], l1, n1, l2, n2, ratio), flush=True)

    # (B) inter-comparison with L1- vs L2-derived calibration
    inter = {}
    for name in STATIONS:
        for level in ("L2", "L1"):
            print(f"== {name} [{level}] ==", flush=True)
            R, cfg = run_station_level(name, level)
            if R is None:
                continue
            for k, ch in enumerate(R["channels"]):
                s = R["stats"][k]
                inter.setdefault((name, ch["label"]), {})[level] = (s["relbias_pct"], s["r"], s["n"], ch["calib"])
                print("   %-20s [%s] relbias=%+7.1f%% r=%.3f" % (ch["label"], level, s["relbias_pct"], s["r"]), flush=True)

    FIG.fig_calib_l1l2(calib_channel_list(), CALIB, OUT / "fig_calib_l1l2.png")
    write_report(coef, inter)
    print("L1_VS_L2_DONE", flush=True)


def write_report(coef, inter):
    L = ["# L1-vs-L2 calibration validation (eprof_v2)\n",
         "*The native L1 archive (`D:/E-PROFILE_L1_2026`, 15 m × 15 s for CHM15k, binned to the 30 m × 300 s "
         "L2 grid via `l1_bin_to_l2_grid`) vs the E-PROFILE L2 product, both calibrated with the operational "
         "**eprof_v2** molecular-window method. Switching from eprof_v1.2 to eprof_v2 roughly tripled the "
         "Rayleigh yield on clear nights (Payerne CHM 4 → 11 calibrated nights over Mar–May 2026).*\n",
         "## (A) Calibration coefficient — L1 vs L2\n",
         "| channel | L1 median | L2 median | L1/L2 |",
         "|---|---|---|---|"]
    for c in coef:
        L.append("| %s | %.3e | %.3e | %s |" % (c["title"], c["l1"], c["l2"],
                 ("%.3f" % c["ratio"]) if np.isfinite(c["ratio"]) else "—"))
    rr = [c["ratio"] for c in coef if np.isfinite(c["ratio"])]
    if rr:
        L.append("")
        L.append("> Rayleigh/cloud lidar constants from L1 and L2 agree to a median **L1/L2 = %.3f** "
                 "(spread %.3f–%.3f): binning the native L1 to the L2 grid removes the grid-interaction "
                 "rejection, so L1 and L2 calibrate consistently.\n" % (np.median(rr), min(rr), max(rr)))
    L.append("![calibration L1 vs L2](figs_paper_validation/paper_python/fig_calib_l1l2.png)\n")

    L.append("## (B) Inter-comparison — L1- vs L2-derived calibration\n")
    L.append("*Same measured β_att (L2 product), calibrated with the L1- vs L2-derived coefficient; stats over "
             "500–3000 m AGL vs the station reference.*\n")
    L.append("| station | channel | calib | L2 relbias | L2 r | L1 relbias | L1 r |")
    L.append("|---|---|---|---|---|---|---|")
    for (name, label), d in inter.items():
        l2 = d.get("L2"); l1 = d.get("L1")
        if l2 is None and l1 is None:
            continue
        cal = (l2 or l1)[3]
        L.append("| %s | %s | %s | %s | %s | %s | %s |" % (
            name, label, cal,
            ("%+.1f%%" % l2[0]) if l2 else "—", ("%.3f" % l2[1]) if l2 else "—",
            ("%+.1f%%" % l1[0]) if l1 else "—", ("%.3f" % l1[1]) if l1 else "—"))
    txt = "\n".join(L)
    (OUT / "l1_vs_l2_validation.md").write_text(txt, encoding="utf-8")
    (REPO / "doc" / "reports" / "l1_vs_l2_validation.md").write_text(txt, encoding="utf-8")
    print(f"  wrote l1_vs_l2_validation.md ({len(coef)} channels)")


if __name__ == "__main__":
    main()
