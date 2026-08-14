# -*- coding: utf-8 -*-
"""Calibration time series per configuration — what the gate change actually does to C_L(t).

One row per station, one panel per configuration, so the effect of loosening the noise tolerance
is visible directly: which nights appear, where they sit relative to the nights v2 already had,
and whether the added points scatter or form a coherent series.

Points are coloured by ORIGIN, which is the whole question:
  grey   nights v2 already calibrated (identical under v2.2 by construction)
  blue   nights only v2.2 recovers
and sized by nothing else, so density is readable. The v2 Kalman is drawn as the reference line,
and the per-station C_L-vs-window-height gradient is annotated because it predicts how far the
recovered points sit from the retained ones.

Run:  python rayleigh_availability/plot_timeseries.py [station-substring ...]
Out:  rayleigh_availability/figs/timeseries_<SITE>.png
"""
from __future__ import annotations
import json
import sys
from datetime import datetime
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.dates as mdates

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))                       # monitoring.kalman
sys.path.insert(0, str(REPO / "rayleigh_availability"))
import indicators as IND                                                       # noqa: E402
from monitoring.kalman import kalman_best_estimate                             # noqa: E402

DATA = Path("C:/DATA/Projects/202606_E-PROFILE_calibration/rayleigh_availability")
BASE, CAND = DATA / "baselines", DATA / "candidates"
MANIFEST = json.loads((REPO / "rayleigh_availability" / "scope_availability.json").read_text())
FIG = REPO / "rayleigh_availability" / "figs"
REF = "eprof_v2"
# Default stations: the three that carry the argument (aged/high-gradient, healthy, big-recovery).
DEFAULT = ["PAYERNE", "AMSTERDAM", "GOTTFRIEDING"]


def load(p):
    return json.loads(p.read_text()) if p.exists() else None


def configs_present():
    return sorted({p.name.split("_", 2)[1] for p in CAND.glob("cand_*.json")},
                  key=lambda s: float(s[1:].split("_")[0]) if s[1:2].isdigit() else 99)


def complete_configs():
    """Configs covering the SAME dates as the baseline, so panels are comparable.

    A config left over from an earlier, shorter corpus would otherwise be drawn next to an extended
    one and the difference in date coverage would read as a difference in method.
    """
    ref_n = {}
    for i in MANIFEST:
        r = load(BASE / f"base_{REF}_{i['label']}.json")
        if r:
            ref_n[i["label"]] = len(r)
    out = []
    for cfg in configs_present():
        ok = 0
        for i in MANIFEST:
            c = load(CAND / f"cand_{cfg}_{i['label']}.json")
            if c and i["label"] in ref_n and len(c) == ref_n[i["label"]]:
                ok += 1
        if ok == len(ref_n):
            out.append(cfg)
    return out


def dt(d):
    return datetime.strptime(d, "%Y%m%d")


def panel(ax, ref, new, title, grad=None):
    """One configuration for one station."""
    kept = {d: new[d][1] for d in new
            if IND.is_valid(new[d][0]) and d in ref and IND.is_valid(ref[d][0]) and new[d][1]}
    added = {d: new[d][1] for d in new
             if IND.is_valid(new[d][0]) and not (d in ref and IND.is_valid(ref[d][0]))
             and new[d][1]}
    if kept:
        ax.plot([dt(d) for d in sorted(kept)], [kept[d] for d in sorted(kept)],
                "o", ms=3.4, color="#777", label=f"kept by v2 ({len(kept)})", zorder=3)
    if added:
        ax.plot([dt(d) for d in sorted(added)], [added[d] for d in sorted(added)],
                "o", ms=3.4, color="#1f77b4", label=f"recovered ({len(added)})", zorder=4)
    # the v2 Kalman as the reference line the operational series would follow
    dts, c = IND.constants(ref)
    if len(c) >= 5:
        kt, ks, _ = kalman_best_estimate([dt(d) for d in dts], c)
        if len(kt):
            ax.plot(kt.astype("datetime64[s]").astype(datetime), ks, "-", color="#d62728",
                    lw=1.6, label="v2 Kalman", zorder=5)
    ax.set_title(title, fontsize=10)
    ax.grid(alpha=0.3)
    # Linear y: a log axis compresses exactly the differences these panels exist to show.
    # The range is set from the RETAINED nights so a few low outliers cannot squash the series.
    if kept:
        v = np.array(sorted(kept.values()), float)
        lo, hi = np.percentile(v, [2, 98])
        pad = 0.6 * (hi - lo) if hi > lo else 0.3 * hi
        ax.set_ylim(max(0.0, lo - pad), hi + pad)
    if grad is not None and np.isfinite(grad):
        ax.text(0.015, 0.05, f"dC$_L$/dz = {grad:+.1f} %/km", transform=ax.transAxes,
                fontsize=8, color="#444",
                bbox=dict(fc="white", ec="#ccc", alpha=0.85, pad=2))


def main():
    want = [a.upper() for a in sys.argv[1:]] or DEFAULT
    cfgs = complete_configs()
    if not cfgs:
        print("no configuration is complete yet"); return
    print(f"complete configs: {cfgs}")

    for pat in want:
        insts = [i for i in MANIFEST if pat in i["site"].upper() and i["group"] == "CHM15k"]
        if not insts:
            print(f"  no CHM15k stream matching {pat!r}"); continue
        inst = insts[0]
        ref = load(BASE / f"base_{REF}_{inst['label']}.json")
        if not ref:
            continue
        ncol = len(cfgs) + 1
        fig, axes = plt.subplots(1, ncol, figsize=(4.1 * ncol, 4.0), sharey=True)
        axes = np.atleast_1d(axes)

        # first panel = the v2 baseline alone (what operations produces today)
        panel(axes[0], ref, ref, f"v2 (C8) — baseline")
        for k, cfg in enumerate(cfgs, start=1):
            new = load(CAND / f"cand_{cfg}_{inst['label']}.json")
            if not new:
                continue
            g = IND.altitude_gradient(new)["slope_pct_per_km"]
            panel(axes[k], ref, new, f"v2.2 {cfg}", grad=g)
        for ax in axes:
            ax.xaxis.set_major_locator(mdates.MonthLocator(interval=3))
            ax.xaxis.set_major_formatter(mdates.DateFormatter("%b\n%Y"))
            ax.legend(fontsize=7.5, loc="upper left")
        axes[0].set_ylabel(r"$C_L$")
        # Keep the title inside the axes width: with two panels the long form was clipped at both
        # ends, hiding the station name itself.
        fig.suptitle(f"{inst['site'][:22]} ({inst['ident']}) — Rayleigh $C_L$ per configuration "
                     f"[{inst['split']}, "
                     f"{'noisy' if inst.get('v2_minus2_pct_clear', 0) > 50 else 'moderate'}]",
                     fontsize=11, fontweight="bold")
        fig.tight_layout(rect=(0, 0, 1, 0.94))
        out = FIG / f"timeseries_{inst['site'].split('_')[0][:14]}.png"
        FIG.mkdir(parents=True, exist_ok=True)
        fig.savefig(out, dpi=140)
        plt.close(fig)
        print(f"  -> {out}")


if __name__ == "__main__":
    main()
