"""Mean normalised in-cloud backscatter vs depth around the peak, per instrument type, from the
extraction of validation/eta_study/extract_incloud_profiles.py (C:/DATA/Projects/202611_eta_study/incloud/).

Every accepted calibration profile is normalised by its own value at the peak gate (the peak the
calibration's cloud filters use, i.e. argmax of beta_after), separately for the beta BEFORE and
AFTER the WV + eta corrections, then averaged per depth:
  pooled mean / median / P25 / P75 (profile-weighted) and the unit-weighted mean (mean of the
  per-unit means, so one prolific unit cannot dominate a type).
Groups: CL31 all, CL31 excl. Payerne CL31 (distorted unit), CL51 all, CL61 all; CBH bins
0.5-1 / 1-1.5 / 1.5-2.4 km and all. Grids: working (calibration) gates for every type, and the
native 4.8 m gates for the CL61.

Run:  python validation/eta_study/incloud_mean_profiles.py
Outputs (C:/DATA/Projects/202611_eta_study/incloud_mean/): incloud_mean_profile_by_type.csv,
incloud_mean_profile_by_unit.csv, incloud_decay_depths.csv, fig_incloud_mean_profiles.png.
"""
from __future__ import annotations

import json
from collections import defaultdict
from pathlib import Path

import numpy as np

SRC = Path("C:/DATA/Projects/202611_eta_study/incloud")
OUT = Path("C:/DATA/Projects/202611_eta_study/incloud_mean")
CBH_BINS = [("all", 0.0, 1e9), ("0.5-1 km", 500.0, 1000.0), ("1-1.5 km", 1000.0, 1500.0),
            ("1.5-2.4 km", 1500.0, 2400.0)]
GROUPS = [("CL31 all", lambda u: u["itype"] == "CL31"),
          ("CL31 excl. Payerne CL31", lambda u: u["itype"] == "CL31" and not u["distorted"]),
          ("CL51 all", lambda u: u["itype"] == "CL51"),
          ("CL61 all", lambda u: u["itype"] == "CL61")]


def normalised(z, grid, variant):
    """(profiles normalised at the peak gate, depth axis) for one unit, grid, variant."""
    pre = "native_" if grid == "native" else ""
    arr = z[f"{pre}beta_{variant}"].astype("float64")
    depth = z[f"{pre}range_rel_m"].astype("float64")
    k0 = int(np.argmin(np.abs(depth)))
    with np.errstate(all="ignore"):
        pk = arr[:, k0:k0 + 1]
        out = np.where(pk > 0, arr / pk, np.nan)
    return out, depth


def decay_depth(depth, prof, level):
    """First depth above the peak where the profile drops below `level` (linear interpolation)."""
    above = depth > 0
    d, p = depth[above], prof[above]
    idx = np.where(p < level)[0]
    if idx.size == 0:
        return np.nan
    i = idx[0]
    if i == 0:
        return float(d[0] * (1 - level) / (1 - p[0])) if p[0] < 1 else float(d[0])
    return float(d[i - 1] + (d[i] - d[i - 1]) * (p[i - 1] - level) / (p[i - 1] - p[i]))


def aggregate(cache, gkeys, gname, grid, variant, lab, lo, hi, type_rows, decay_rows):
    """Pooled and unit-weighted mean profile of one group / grid / variant / CBH bin."""
    Ps, unit_means = [], []
    depth = None
    for k in gkeys:
        P, depth, cbh = cache[(k, grid, variant)]
        m = (cbh >= lo) & (cbh < hi)
        if m.sum() >= 20:
            Ps.append(P[m])
            with np.errstate(all="ignore"):
                unit_means.append(np.nanmean(P[m], axis=0))
    if not Ps:
        return
    Q = np.vstack(Ps)
    uw = np.nanmean(np.vstack(unit_means), axis=0)
    with np.errstate(all="ignore"):
        mean = np.nanmean(Q, axis=0)
    for d, mu, md, q1, q3, nn, w in zip(depth, mean, np.nanmedian(Q, axis=0),
                                        np.nanpercentile(Q, 25, axis=0), np.nanpercentile(Q, 75, axis=0),
                                        np.isfinite(Q).sum(axis=0), uw):
        type_rows.append(dict(group=gname, grid=grid, variant=variant, cbh_bin=lab,
                              depth_m=round(float(d), 2), n_profiles=int(nn), n_units=len(Ps),
                              mean=mu, median=md, p25=q1, p75=q3, unit_weighted_mean=w))
    decay_rows.append(dict(group=gname, grid=grid, variant=variant, cbh_bin=lab,
                           n_profiles=Q.shape[0], n_units=len(Ps),
                           gate_m=float(np.median(np.diff(depth))),
                           half_depth_m=decay_depth(depth, mean, 0.5),
                           e_fold_depth_m=decay_depth(depth, mean, np.exp(-1)),
                           half_depth_unitw_m=decay_depth(depth, uw, 0.5),
                           mean_at_minus_gate=float(mean[np.argmin(np.abs(depth)) - 1]),
                           mean_at_plus_50m=float(np.interp(50.0, depth, mean))))


def main():
    import pandas as pd
    OUT.mkdir(parents=True, exist_ok=True)
    units = json.load(open(SRC / "units.json", encoding="utf-8"))
    data = {}
    for key, u in units.items():
        fp = SRC / f"{key}.npz"
        if fp.exists():
            data[key] = (u, np.load(fp, allow_pickle=True))
    print(f"{len(data)} unit files", flush=True)

    unit_rows, type_rows, decay_rows = [], [], []
    cache = {}
    for key, (u, z) in data.items():
        cbh = z["cbh_m"].astype("float64")
        for grid in ("working", "native"):
            if grid == "native" and "native_beta_after" not in z.files:
                continue
            for variant in ("before", "after"):
                P, depth = normalised(z, grid, variant)
                cache[(key, grid, variant)] = (P, depth, cbh)
                for lab, lo, hi in CBH_BINS:
                    m = (cbh >= lo) & (cbh < hi)
                    if m.sum() < 20:
                        continue
                    Q = P[m]
                    with np.errstate(all="ignore"):
                        mean = np.nanmean(Q, axis=0)
                    for d, mu, md, q1, q3, nn in zip(depth, mean, np.nanmedian(Q, axis=0),
                                                     np.nanpercentile(Q, 25, axis=0),
                                                     np.nanpercentile(Q, 75, axis=0),
                                                     np.isfinite(Q).sum(axis=0)):
                        unit_rows.append(dict(key=key, itype=u["itype"], site=u["site"],
                                              distorted_unit=u["distorted"], grid=grid,
                                              variant=variant, cbh_bin=lab, depth_m=round(float(d), 2),
                                              n=int(nn), mean=mu, median=md, p25=q1, p75=q3))

    for gname0, sel in GROUPS:
        keys = [k for k, (u, _) in data.items() if sel(u)]
        for grid in ("working", "native"):
            gk_all = [k for k in keys if (k, grid, "after") in cache]
            # units configured at another gate spacing (CL31 at 20 m) form their own sub-group
            by_dr = defaultdict(list)
            for k in gk_all:
                by_dr[round(float(np.median(np.diff(cache[(k, grid, "after")][1]))), 1)].append(k)
            for dr, gkeys in sorted(by_dr.items()):
                gname = gname0 if len(by_dr) == 1 or dr == min(by_dr) else f"{gname0} [{dr:g} m gates]"
                for variant in ("before", "after"):
                    for lab, lo, hi in CBH_BINS:
                        aggregate(cache, gkeys, gname, grid, variant, lab, lo, hi, type_rows, decay_rows)

    tdf = pd.DataFrame(type_rows)
    pd.DataFrame(unit_rows).to_csv(OUT / "incloud_mean_profile_by_unit.csv", index=False, float_format="%.5g")
    tdf.to_csv(OUT / "incloud_mean_profile_by_type.csv", index=False, float_format="%.5g")
    ddf = pd.DataFrame(decay_rows)
    ddf.to_csv(OUT / "incloud_decay_depths.csv", index=False, float_format="%.4g")
    print(ddf[ddf.cbh_bin == "all"].to_string(index=False), flush=True)
    print(ddf[(ddf.variant == "after") & (ddf.grid == "working")].to_string(index=False), flush=True)

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(1, 4, figsize=(20, 6.2), sharey=True)
    panels = [("CL31 all", "working"), ("CL51 all", "working"), ("CL61 all", "working"), ("CL61 all", "native")]
    colors = {"0.5-1 km": "C0", "1-1.5 km": "C1", "1.5-2.4 km": "C3"}
    for ax, (g, grid) in zip(axes, panels):
        sub = tdf[(tdf.group == g) & (tdf.grid == grid)]
        if sub.empty:
            continue
        a = sub[(sub.variant == "after") & (sub.cbh_bin == "all")]
        b = sub[(sub.variant == "before") & (sub.cbh_bin == "all")]
        ax.fill_betweenx(a.depth_m, a.p25, a.p75, color="k", alpha=0.10, label="after: P25-P75")
        ax.plot(a["mean"], a.depth_m, "k-", lw=2.2, label=f"after WV+eta, all CBH (n={int(a.n_profiles.max())})")
        ax.plot(b["mean"], b.depth_m, "k--", lw=1.5, label="before corrections, all CBH")
        for lab, c in colors.items():
            s = sub[(sub.variant == "after") & (sub.cbh_bin == lab)]
            if len(s):
                ax.plot(s["mean"], s.depth_m, "-", color=c, lw=1.4, label=f"after, CBH {lab} (n={int(s.n_profiles.max())})")
        if g == "CL31 all":
            s = tdf[(tdf.group == "CL31 excl. Payerne CL31") & (tdf.grid == grid) & (tdf.variant == "after")
                    & (tdf.cbh_bin == "all")]
            ax.plot(s["mean"], s.depth_m, ":", color="C2", lw=2, label="after, excl. Payerne CL31")
        ax.axhline(0, color="grey", lw=0.8)
        ax.set_xlim(0, 1.15)
        ax.set_title(f"{g} -- {grid} gates ({np.median(np.diff(a.depth_m)):.1f} m)")
        ax.set_xlabel("beta / beta(peak gate)")
        ax.grid(alpha=0.3)
        ax.legend(fontsize=7)
    axes[0].set_ylabel("depth relative to the backscatter peak [m]")
    fig.suptitle("Mean normalised in-cloud attenuated backscatter around the peak at the accepted "
                 "cloud-calibration profiles (2025-01..2026-06)", fontweight="bold")
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    fig.savefig(OUT / "fig_incloud_mean_profiles.png", dpi=140)
    plt.close(fig)
    print("done", flush=True)


if __name__ == "__main__":
    main()
