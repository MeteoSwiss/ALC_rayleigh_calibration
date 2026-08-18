"""
_amsterdam_noB_figure.py — the Amsterdam station figure (report Figure 3) with the problematic
unit B excluded: the healthy trio A (reference), C, D.

The site is processed once through the paper pipeline (run_paper_validation.run_site) and unit B
is dropped from the result before rendering. The pairwise C/D statistics vs A are unchanged by
construction (they never involved B); what changes is the common-hours sample of the median
panel (no longer constrained by B's availability) and the 3-curtain layout — i.e. the figure the
station would produce if B were serviced or excluded.

Usage: python -m validation.paper._amsterdam_noB_figure
"""
from __future__ import annotations
import sys
import warnings
from pathlib import Path

OUT = Path("C:/DATA/Projects/202606_E-PROFILE_calibration/figs_paper_validation/paper_python")


def main():
    warnings.filterwarnings("ignore")
    from validation.paper import run_paper_validation as RPV
    from validation.paper import figures as FIG

    print("processing amsterdam through the paper pipeline ...", flush=True)
    R, cfg, _rows = RPV.run_site("amsterdam")
    labels = [c["label"] for c in R["channels"]]
    iB = labels.index("CHM15k B")

    # drop unit B from every per-channel list (reference A stays at index 0)
    for key in ("channels", "beta", "beta_disp", "beta_noovl", "cbh", "stats"):
        R[key] = [v for k, v in enumerate(R[key]) if k != iB]
    cfg = dict(cfg, channels=[c for k, c in enumerate(cfg["channels"]) if k != iB])

    png = OUT / "fig_amsterdam_noB.png"
    FIG.fig_multi_alc(R, cfg, png,
                      "Amsterdam (0-20000-0-06240)  —  2026-03-01 to 2026-05-31"
                      "  —  excluding unit B")
    for k, ch in enumerate(R["channels"]):
        s = R["stats"][k]
        print("   %-12s med %+6.1f%%  logr %.2f  relbias %+6.1f%%  r %.2f  N=%7d"
              % (ch["label"], s["medrelbias_pct"], s["r_log"], s["relbias_pct"], s["r"], s["n"]))
    print(f"-> {png.name}", flush=True)
    print("AMSTERDAM_NOB_DONE", flush=True)


if __name__ == "__main__":
    try:
        sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    except Exception:
        pass
    main()
