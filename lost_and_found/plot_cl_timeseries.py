#!/usr/bin/env python3
"""
Build the Payerne CHM15k lidar-constant comparison plots: with vs without the last
commit (Klett sign + aerosol-OD fix).

Reads payerne_cl_withfix.csv and payerne_cl_nofix.csv and writes, under comparison/:
  - payerne_cl_timeseries.png : CL vs time, both configs, with uncertainty bands
  - payerne_cl_difference.png : per-night relative difference (withfix - nofix) + stats
Also prints a text summary.
"""
from __future__ import annotations

import csv
from datetime import datetime
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import dates as mdates

BASE = Path("C:/DATA/Projects/202606_E-PROFILE_calibration/E-PROFILE_calibration_rayleigh")
OUT = BASE / "comparison"
SUCCESS_FLAGS = {1.0, 0.5}


def _load(csv_path: Path) -> dict[str, dict]:
    """date -> {dt, cl, unc, flag} for successful nights only."""
    out = {}
    with open(csv_path, newline="", encoding="utf-8") as f:
        for row in csv.DictReader(f):
            try:
                flag = float(row["flag"])
            except (ValueError, KeyError):
                continue
            if flag not in SUCCESS_FLAGS:
                continue
            out[row["date"]] = {
                "dt": datetime.strptime(row["date"], "%Y%m%d"),
                "cl": float(row["lidar_constant"]),
                "unc": float(row["uncertainty"]),
                "flag": flag,
            }
    return out


def _series(d: dict[str, dict]):
    items = sorted(d.values(), key=lambda r: r["dt"])
    return (np.array([r["dt"] for r in items]),
            np.array([r["cl"] for r in items]),
            np.array([r["unc"] for r in items]))


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    wf = _load(BASE / "payerne_cl_withfix.csv")
    nf = _load(BASE / "payerne_cl_nofix.csv")

    wf_dt, wf_cl, wf_unc = _series(wf)
    nf_dt, nf_cl, nf_unc = _series(nf)

    # ---- Plot 1: time series ------------------------------------------------
    fig, ax = plt.subplots(figsize=(15, 6))
    for dts, cl, unc, color, lab in [
        (nf_dt, nf_cl, nf_unc, "tab:red", "without fix (HEAD~1)"),
        (wf_dt, wf_cl, wf_unc, "tab:blue", "with fix (HEAD)"),
    ]:
        ax.fill_between(dts, cl - unc, cl + unc, color=color, alpha=0.15)
        ax.plot(dts, cl, "o-", color=color, ms=3, lw=0.8, label=lab)
    ax.set_ylabel("Lidar constant  (counts/s·m³/sr)")
    ax.set_title("Payerne CHM15k (0-20000-0-06610) — Rayleigh lidar constant\n"
                 "with vs without the Klett sign + aerosol-OD fix")
    ax.xaxis.set_major_locator(mdates.YearLocator())
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%Y"))
    ax.grid(True, alpha=0.3)
    ax.legend()
    fig.tight_layout()
    fig.savefig(OUT / "payerne_cl_timeseries.png", dpi=140, bbox_inches="tight")
    plt.close(fig)

    # ---- Plot 2: relative difference on common nights -----------------------
    common = sorted(set(wf) & set(nf))
    cdt = np.array([wf[d]["dt"] for d in common])
    rel = np.array([(wf[d]["cl"] - nf[d]["cl"]) / nf[d]["cl"] * 100 for d in common])

    fig2, ax2 = plt.subplots(figsize=(15, 5))
    ax2.axhline(0, color="k", lw=0.6)
    ax2.plot(cdt, rel, "o", color="tab:purple", ms=3)
    ax2.set_ylabel("(with − without) / without  [%]")
    ax2.set_title(f"Payerne CHM15k — per-night CL difference on {len(common)} common nights")
    ax2.xaxis.set_major_locator(mdates.YearLocator())
    ax2.xaxis.set_major_formatter(mdates.DateFormatter("%Y"))
    ax2.grid(True, alpha=0.3)
    if len(rel):
        txt = (f"median Δ = {np.median(rel):+.2f}%\n"
               f"mean Δ   = {np.mean(rel):+.2f}%\n"
               f"max |Δ|  = {np.max(np.abs(rel)):.2f}%")
        ax2.text(0.01, 0.97, txt, transform=ax2.transAxes, va="top", family="monospace",
                 bbox=dict(boxstyle="round", fc="white", alpha=0.8))
    fig2.tight_layout()
    fig2.savefig(OUT / "payerne_cl_difference.png", dpi=140, bbox_inches="tight")
    plt.close(fig2)

    # ---- Text summary -------------------------------------------------------
    print("=" * 60)
    print("Payerne CHM15k — CL comparison summary")
    print("=" * 60)
    print(f"successful nights  : with fix = {len(wf)}   without fix = {len(nf)}")
    print(f"common nights      : {len(common)}")
    if len(wf_cl):
        print(f"median CL with fix : {np.median(wf_cl):.4e}")
    if len(nf_cl):
        print(f"median CL no fix   : {np.median(nf_cl):.4e}")
    if len(rel):
        print(f"median rel. diff   : {np.median(rel):+.3f}%   max |diff| = {np.max(np.abs(rel)):.3f}%")
    only_wf = sorted(set(wf) - set(nf))
    only_nf = sorted(set(nf) - set(wf))
    print(f"successful ONLY with fix : {len(only_wf)}  {only_wf[:10]}")
    print(f"successful ONLY without  : {len(only_nf)}  {only_nf[:10]}")
    print(f"plots -> {OUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
