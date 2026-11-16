# -*- coding: utf-8 -*-
"""Fig 25 — the Christmas-2025 clean-airmass experiment (operator's idea): on 2025-12-25/26 an
exceptionally clean, shallow airmass covered NL + DE; 77 CHM15k measured the SAME near-empty
free troposphere. The cross-station median of the relative residual

    delta_u(z) = (S - C M) / (C M)      (fraction of the molecular signal)

is the shared atmosphere + model error; each station's DEVIATION from it is its own
instrument systematic (dark + overlap error + near-range artefacts) measured on the cleanest
nights the archive offers — the distributed-co-location idea: the airmass is the common
reference that a single station never has.

Reading guide: below ~1.6 km a deviation can be multiplicative (overlap error — Schiphol B's
known -25 % deficit is the built-in validation); higher up the additive dark dominates
because M collapses. The ranking (right panel) is each unit's instrument systematic in % of
molecular over 1-3 km.
"""
from __future__ import annotations

import json
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from remote_dark.common import DATA, ensure_out

FIG = ensure_out("figs")
NPZ = DATA / "remote_dark" / "xmas_residuals.npz"
ZG = np.arange(300.0, 9001.0, 30.0)
JBAND = (1000.0, 3000.0)


def smooth(x, k=9):
    pad = np.pad(x, k // 2, mode="edge")
    return np.array([np.nanmedian(pad[i:i + k]) for i in range(x.size)])


def region(key):
    if key.startswith("0-20000-0-06") or key.startswith("0-528"):
        return "NL"
    if key.startswith("0-20000-0-10") or key.startswith("0-276"):
        return "DE"
    if key.startswith("0-20000-0-03"):
        return "UK"
    return "other"


def main():
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    sites = {f"{s['wmo']}_{s['ident']}": s["site"]
             for s in json.load(open("remote_dark/chm_streams.json"))}
    z = np.load(NPZ, allow_pickle=True)
    keys = sorted({k[len("rng_"):] for k in z.files if k.startswith("rng_")})
    prof = {}
    for key in keys:
        rng, delta = z[f"rng_{key}"], z[f"delta_{key}"]
        d = np.nanmean(np.vstack([smooth(row) for row in delta]), axis=0)
        prof[key] = np.interp(ZG, rng, d, left=np.nan, right=np.nan)
    A = np.vstack([prof[k] for k in keys])
    com = np.nanmedian(A, axis=0)
    dev = A - com[None, :]

    aer = np.nanmedian(np.vstack([
        np.interp(ZG, z[f"rng_{k}"], smooth(z[f"aer_{k}"]), left=np.nan, right=np.nan)
        for k in keys]), axis=0)

    jm = (ZG >= JBAND[0]) & (ZG <= JBAND[1])
    rank = np.nanmedian(dev[:, jm], axis=1)

    fig = plt.figure(figsize=(16.5, 7.0))
    gs = fig.add_gridspec(1, 3, width_ratios=[1.0, 1.0, 0.9], wspace=0.24)
    ax1, ax2, ax3 = [fig.add_subplot(gs[0, i]) for i in range(3)]
    schiphol = {f"0-20000-0-06240_{i}": c for i, c in
                zip("ABCD", ("crimson", "darkorange", "seagreen", "purple"))}

    for i, key in enumerate(keys):
        ax1.plot(A[i] * 100, ZG / 1e3, "-", color=schiphol.get(key, "0.75"),
                 lw=1.6 if key in schiphol else 0.6,
                 alpha=1.0 if key in schiphol else 0.5)
    ax1.plot(com * 100, ZG / 1e3, "k-", lw=2.4, label="common (77-station median)")
    ax1.plot(aer * 100, ZG / 1e3, "--", color="steelblue", lw=1.6,
             label="CAMS-expected aerosol / molecular")
    ax1.axvline(0, color="k", lw=0.5)
    ax1.set_xlim(-40, 80)
    ax1.set_xlabel("residual / molecular (%)")
    ax1.set_ylabel("altitude (km)")
    ax1.set_title("δ(z) per station — 2025-12-25/26\n(grey; Schiphol A–D coloured)",
                  fontsize=10)
    ax1.legend(fontsize=8, loc="upper right")
    ax1.grid(alpha=0.25)

    for i, key in enumerate(keys):
        ax2.plot(dev[i] * 100, ZG / 1e3, "-", color=schiphol.get(key, "0.75"),
                 lw=1.6 if key in schiphol else 0.6,
                 alpha=1.0 if key in schiphol else 0.5)
    ax2.axvline(0, color="k", lw=0.5)
    spread = np.nanpercentile(dev[:, jm], [25, 75])
    ax2.set_xlim(-40, 80)
    ax2.set_xlabel("deviation from common (%)")
    ax2.set_title(f"per-unit instrument systematic\nIQR at 1–3 km: "
                  f"{spread[0] * 100:+.1f}..{spread[1] * 100:+.1f} %", fontsize=10)
    ax2.grid(alpha=0.25)

    order = np.argsort(rank)
    ypos = np.arange(len(keys))
    for iy, idx in enumerate(order):
        key = keys[idx]
        colr = schiphol.get(key, {"NL": "steelblue", "DE": "0.4",
                                  "UK": "darkorange"}.get(region(key), "0.7"))
        ax3.plot(rank[idx] * 100, iy, "o", ms=4, color=colr)
    ax3.axvline(0, color="k", lw=0.6)
    ax3.set_yticks(ypos[::2])
    ax3.set_yticklabels([sites.get(keys[i], keys[i])[:14] for i in order[::2]], fontsize=5)
    ax3.set_xlabel("median deviation 1–3 km (% of molecular)")
    ax3.set_title("unit ranking (blue NL, grey DE, orange UK;\nSchiphol A–D coloured)",
                  fontsize=10)
    ax3.grid(alpha=0.25, axis="x")

    print("top |deviation| units (1-3 km, % of molecular):")
    for idx in np.argsort(-np.abs(rank))[:10]:
        print(f"  {sites.get(keys[idx], keys[idx]):24.24s} {rank[idx] * 100:+7.1f} %")
    for k in schiphol:
        if k in keys:
            i = keys.index(k)
            lo = np.nanmedian(dev[i][(ZG >= 500) & (ZG <= 1500)])
            print(f"  Schiphol {k[-1]}: 1-3 km {rank[i] * 100:+.1f} %   "
                  f"0.5-1.5 km {lo * 100:+.1f} %")
    fig.suptitle("The Christmas-2025 clean-airmass experiment — 77 CHM15k under one near-empty "
                 "sky: common part = atmosphere+model, deviations = the instruments",
                 fontsize=11)
    f = FIG / "fig25_xmas_experiment.png"
    fig.savefig(f, dpi=140, bbox_inches="tight")
    plt.close(fig)
    print(f)


if __name__ == "__main__":
    main()
