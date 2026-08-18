#!/usr/bin/env python3
"""
CL comparison plots (with vs without the Klett sign + aerosol-OD fix) for any site.

Reads C:/DATA/Projects/202606_E-PROFILE_calibration/E-PROFILE_calibration_rayleigh/<site>/<site>_cl_{withfix,nofix}.csv and writes
to <site>/comparison/:
  <site>_cl_timeseries.png   CL vs time, both configs, uncertainty bands
  <site>_cl_difference.png   per-night relative difference + summary stats
Prints a text summary.

Usage:  python plot_cl_timeseries_site.py --site granada
"""
from __future__ import annotations

import argparse
import csv
from datetime import datetime
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import dates as mdates

BASE = Path("C:/DATA/Projects/202606_E-PROFILE_calibration/E-PROFILE_calibration_rayleigh")
SUCCESS_FLAGS = {1.0, 0.5}


def _load(csv_path: Path) -> dict[str, dict]:
    out = {}
    if not csv_path.exists():
        return out
    for row in csv.DictReader(open(csv_path, newline="", encoding="utf-8")):
        try:
            if float(row["flag"]) not in SUCCESS_FLAGS:
                continue
            out[row["date"]] = {
                "dt": datetime.strptime(row["date"], "%Y%m%d"),
                "cl": float(row["lidar_constant"]),
                "unc": float(row["uncertainty"]),
            }
        except (ValueError, KeyError):
            pass
    return out


def _series(d):
    items = sorted(d.values(), key=lambda r: r["dt"])
    return (np.array([r["dt"] for r in items]),
            np.array([r["cl"] for r in items]),
            np.array([r["unc"] for r in items]))


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--site", required=True)
    args = ap.parse_args()
    site = args.site
    out_dir = BASE / site / "comparison"
    out_dir.mkdir(parents=True, exist_ok=True)

    wf = _load(BASE / site / f"{site}_cl_withfix.csv")
    nf = _load(BASE / site / f"{site}_cl_nofix.csv")
    if not wf and not nf:
        print(f"[{site}] no data — skipping"); return 1

    wf_dt, wf_cl, wf_unc = _series(wf)
    nf_dt, nf_cl, nf_unc = _series(nf)

    # ---- time series ----
    fig, ax = plt.subplots(figsize=(15, 6))
    for dts, cl, unc, color, lab in [
        (nf_dt, nf_cl, nf_unc, "tab:red", "without fix (HEAD~1)"),
        (wf_dt, wf_cl, wf_unc, "tab:blue", "with fix (HEAD)"),
    ]:
        if len(dts):
            ax.fill_between(dts, cl - unc, cl + unc, color=color, alpha=0.15)
            ax.plot(dts, cl, "o-", color=color, ms=3, lw=0.8, label=lab)
    ax.set_ylabel("Lidar constant  (counts/s·m³/sr)")
    ax.set_title(f"{site.capitalize()} — Rayleigh lidar constant\n"
                 f"with vs without the Klett sign + aerosol-OD fix")
    ax.xaxis.set_major_locator(mdates.YearLocator())
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%Y"))
    ax.grid(True, alpha=0.3); ax.legend()
    fig.tight_layout()
    fig.savefig(out_dir / f"{site}_cl_timeseries.png", dpi=140, bbox_inches="tight")
    plt.close(fig)

    # ---- difference ----
    common = sorted(set(wf) & set(nf))
    rel = np.array([(wf[d]["cl"] - nf[d]["cl"]) / nf[d]["cl"] * 100 for d in common])
    cdt = np.array([wf[d]["dt"] for d in common])
    fig2, ax2 = plt.subplots(figsize=(15, 5))
    ax2.axhline(0, color="k", lw=0.6)
    if len(rel):
        ax2.plot(cdt, rel, "o", color="tab:purple", ms=3)
        txt = (f"median Δ = {np.median(rel):+.2f}%\nmean Δ   = {np.mean(rel):+.2f}%\n"
               f"max |Δ|  = {np.max(np.abs(rel)):.2f}%")
        ax2.text(0.01, 0.97, txt, transform=ax2.transAxes, va="top", family="monospace",
                 bbox=dict(boxstyle="round", fc="white", alpha=0.8))
    ax2.set_ylabel("(with − without) / without  [%]")
    ax2.set_title(f"{site.capitalize()} — per-night CL difference on {len(common)} common nights")
    ax2.xaxis.set_major_locator(mdates.YearLocator())
    ax2.xaxis.set_major_formatter(mdates.DateFormatter("%Y"))
    ax2.grid(True, alpha=0.3)
    fig2.tight_layout()
    fig2.savefig(out_dir / f"{site}_cl_difference.png", dpi=140, bbox_inches="tight")
    plt.close(fig2)

    # ---- summary ----
    print(f"=== {site} ===")
    print(f"  successful: withfix={len(wf)} nofix={len(nf)} common={len(common)}")
    if len(wf_cl):
        print(f"  median CL withfix = {np.median(wf_cl):.4e}")
    if len(nf_cl):
        print(f"  median CL nofix   = {np.median(nf_cl):.4e}")
    if len(rel):
        print(f"  median rel diff = {np.median(rel):+.3f}%  max|diff| = {np.max(np.abs(rel)):.3f}%")
    print(f"  plots -> {out_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
