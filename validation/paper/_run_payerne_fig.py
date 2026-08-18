# -*- coding: utf-8 -*-
"""Regenerate ONLY the Payerne multi-instrument intercomparison figure, now including the CHM15k
(Rayleigh, offset-corr) channel alongside the CL61 (Rayleigh, offset-corr) one, to see whether the
hood-offset corrections improve the agreement vs the CHM15k (Rayleigh) reference."""
import sys, warnings
warnings.filterwarnings("ignore")
sys.path.insert(0, r"C:/Users/hervo/OneDrive/Documents/ALC_rayleigh_calibration")
try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass
from validation.paper import run_paper_validation as RPV
from validation.paper import figures as FIG
from validation.paper.calib_benchmark import OUT

R, cfg = RPV.run_station("payerne")
if R is None:
    print("no data"); sys.exit(1)
print("\nchannel                              relbias   med-rel     r    N")
for k, ch in enumerate(R["channels"]):
    s = R["stats"][k]
    ref = "  <-- reference" if k == cfg["referenceChannel"] else ""
    print(f"  {ch['label']:34s} {s['relbias_pct']:+6.1f}%  {s.get('medrelbias_pct', float('nan')):+6.1f}%  "
          f"{s['r']:.3f}  {s['n']:6d}{ref}")
title = "%s (%s)  —  %s to %s" % (RPV.SITE_NAME["payerne"], cfg["wmo"], RPV._d(cfg["start"]), RPV._d(cfg["end"]))
FIG.fig_multi_alc(R, cfg, OUT / "fig_payerne.png", title)
print("\n-> fig_payerne.png (", OUT / "fig_payerne.png", ")")
print("PAYERNE_FIG_DONE")
