"""
variability_metrics_l1_2026.py — per-instrument variability of the L1-2026 Rayleigh
calibration, using the drift-insensitive metrics defined for the long-run study
(precision_longrun.py): valid%, sigma_night, sigma_SD (successive-difference), sigma_detrend,
sigma_within-month, plus the robust night-to-night CV (1.4826*MAD/median). calipso dropped.

Reads results_<label>.json (run_l1_2026_variability.py output), writes:
  figs_paper_validation/l1_2026_variability/*.png          (figures)
  .../metrics_tables.md                                    (markdown tables)
  .../metrics_summary.json                                 (key numbers for the report)
"""
from __future__ import annotations
import json
import sys
from collections import defaultdict
from pathlib import Path

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from datetime import datetime

REPO = Path(__file__).resolve().parents[1]
MANIFEST = json.loads((REPO / "validation" / "scope_l1_2026.json").read_text())
OUT = Path("C:/DATA/Projects/202606_E-PROFILE_calibration/figs_paper_validation/l1_2026_variability")
OUT.mkdir(parents=True, exist_ok=True)

QC = 15.0
I_OK, I_CL, I_CLERR, I_REL = 0, 1, 2, 3

from calibration.rayleigh.molecular_methods import METHODS as _LIVE, METHOD_LABELS
METHODS = tuple(_LIVE)                    # 6 live selectable methods (E-PROF version keys)
ALL_METHODS = ("eprof_v1.0",) + METHODS   # + the E-PROF v1.0 sign-error baseline (comparison only)
HEADLINE = "eprof_v2"                     # production-recommended (most precise by sigma_SD)
METHOD_LABEL = dict(METHOD_LABELS)
METHOD_COLORS = {"eprof_v1.0": "#000000", "eprof_v1.1": "#7f7f7f", "eprof_v1.2": "#1f77b4",
                 "eprof_v0.25": "#2ca02c", "earlinet": "#9467bd", "eprof_v2": "#d62728",
                 "bellini": "#8c564b"}
GROUP_COLORS = {"CHM15k": "#1f77b4", "Mini-MPL": "#2ca02c", "CL61": "#d62728"}
GROUPS = ("CHM15k", "Mini-MPL", "CL61")


# ---------------------------------------------------------------------------
# Metric definitions (verbatim from precision_longrun.py — "the metrics we defined")
# ---------------------------------------------------------------------------
def _robust_rel(x, med):
    if len(x) < 3 or med == 0:
        return np.nan
    return float(1.4826 * np.median(np.abs(np.asarray(x) - np.median(x))) / abs(med) * 100)


def _rolling_median(y, k=9):
    n = len(y)
    h = k // 2
    return np.array([np.median(y[max(0, i - h):min(n, i + h + 1)]) for i in range(n)])


def precision(dslist, cls, clerrs):
    cls = np.asarray(cls, float)
    med = np.median(cls) if cls.size else np.nan
    out = dict(n=len(cls), cv=np.nan, rob_cv=np.nan, sigma_night=np.nan,
               sigma_sd=np.nan, sigma_dt=np.nan, sigma_im=np.nan)
    if len(cls) < 4 or med <= 0:
        if len(cls):
            rel = np.asarray(clerrs) / np.abs(cls)
            out["sigma_night"] = float(np.median(rel[np.isfinite(rel)]) * 100) if np.any(np.isfinite(rel)) else np.nan
        return out
    out["cv"] = float(np.std(cls) / abs(np.mean(cls)) * 100)
    out["rob_cv"] = float(1.4826 * np.median(np.abs(cls - med)) / abs(med) * 100)
    rel = np.asarray(clerrs) / np.abs(cls)
    out["sigma_night"] = float(np.median(rel[np.isfinite(rel)]) * 100)
    diffs = np.diff(cls)
    out["sigma_sd"] = float(1.4826 * np.median(np.abs(diffs)) / np.sqrt(2) / abs(med) * 100)
    out["sigma_dt"] = _robust_rel(cls - _rolling_median(cls, 9), med)
    bymonth = defaultdict(list)
    for ds, c in zip(dslist, cls):
        bymonth[ds[:6]].append(c)
    pooled = []
    for mlist in bymonth.values():
        if len(mlist) >= 2:
            pooled.extend(np.asarray(mlist) - np.median(mlist))
    out["sigma_im"] = (float(1.4826 * np.median(np.abs(pooled)) / abs(med) * 100)
                       if len(pooled) >= 3 else np.nan)
    return out


# ---------------------------------------------------------------------------
# Load + per (instrument, method) metrics
# ---------------------------------------------------------------------------
def load_all():
    """rows: per (inst, method); series: {(label, method): (dslist, cls)} for time plots."""
    rows = []
    series = {}
    meta = {m["label"]: m for m in MANIFEST}
    for m in MANIFEST:
        label = m["label"]
        fp = OUT / f"results_{label}.json"
        if not fp.exists():
            print(f"  (no results for {label})")
            continue
        data = json.loads(fp.read_text())
        n_fit = len(data)
        order = sorted(data)
        for meth in ALL_METHODS:                 # incl. eprof_v1.0 (may be absent on some nights)
            ds_ok, cls, errs = [], [], []
            for ds in order:
                w = data[ds].get(meth)
                if w is None:
                    continue
                if w[I_OK] and np.isfinite(w[I_REL]) and w[I_REL] <= QC and np.isfinite(w[I_CL]) and w[I_CL] > 0:
                    ds_ok.append(ds); cls.append(w[I_CL]); errs.append(w[I_CLERR])
            p = precision(ds_ok, cls, errs)
            p.update(inst=label, group=m["group"], site=m["site"], type=m["type"],
                     method=meth, n_fit=n_fit, n_valid=len(cls),
                     valid_pct=100.0 * len(cls) / n_fit if n_fit else np.nan)
            rows.append(p)
            series[(label, meth)] = (ds_ok, cls)
    return rows, series, meta


def agg_by(rows, method, group=None):
    """Mean metric over instruments for one method (optionally one group)."""
    mr = [r for r in rows if r["method"] == method and (group is None or r["group"] == group)]
    def mn(k):
        v = [r[k] for r in mr if np.isfinite(r[k])]
        return float(np.mean(v)) if v else np.nan
    return {k: mn(k) for k in ("valid_pct", "sigma_night", "sigma_sd", "sigma_dt", "sigma_im", "rob_cv", "cv")}


# Preference order for the per-instrument "recommended" estimate: most-precise quality method
# first, falling back to the next when a stricter one yields too few valid nights to be robust.
# (E-PROF v1.0 sign-error is excluded — it is a baseline, never a recommendation.)
REC_PREF = ("eprof_v2", "earlinet", "eprof_v1.2", "eprof_v0.25", "eprof_v1.1")
MIN_VALID_REC = 6


def recommended(rows, label):
    """Per-instrument: the most-precise quality method that has >= MIN_VALID_REC valid nights."""
    by_m = {r["method"]: r for r in rows if r["inst"] == label}
    for m in REC_PREF:
        r = by_m.get(m)
        if r and r["n_valid"] >= MIN_VALID_REC and np.isfinite(r["sigma_sd"]):
            return m, r
    cand = [(by_m[m]["n_valid"], m) for m in by_m if np.isfinite(by_m[m]["sigma_sd"])]
    if cand:
        m = max(cand)[1]
        return m, by_m[m]
    return None, None


# ---------------------------------------------------------------------------
# Figures
# ---------------------------------------------------------------------------
def fig_method_precision(rows):
    """Pooled over all instruments: valid% (left) + drift-insensitive sigma metrics (right)."""
    agg = {m: agg_by(rows, m) for m in ALL_METHODS}
    order = sorted(ALL_METHODS, key=lambda m: agg[m]["sigma_sd"] if np.isfinite(agg[m]["sigma_sd"]) else 1e9)
    fig, ax = plt.subplots(1, 2, figsize=(16, 6))
    x = np.arange(len(order))
    ax[0].bar(x, [agg[m]["valid_pct"] for m in order], color=[METHOD_COLORS[m] for m in order])
    ax[0].set_title("Yield — valid calibrations (% of fit-nights)")
    ax[0].set_ylabel("valid nights (%)")
    ax[0].set_xticks(x); ax[0].set_xticklabels([METHOD_LABEL[m] for m in order], rotation=25, ha="right", fontsize=8)
    ax[0].grid(True, axis="y", alpha=0.25)
    metrics = [("sigma_night", "σ_night\n(within-night)"), ("sigma_sd", "σ_SD\n(successive)"),
               ("sigma_dt", "σ_detrend"), ("sigma_im", "σ_within-month")]
    nm = len(order); bw = 0.8 / nm
    gx = np.arange(len(metrics))
    for j, m in enumerate(order):
        off = (j - (nm - 1) / 2) * bw
        ax[1].bar(gx + off, [agg[m][k] for k, _ in metrics], bw, color=METHOD_COLORS[m], label=METHOD_LABEL[m])
    ax[1].set_title("Precision (drift-insensitive), colour = method — lower = more precise")
    ax[1].set_ylabel("σ (% of median C_L)")
    ax[1].set_xticks(gx); ax[1].set_xticklabels([lab for _, lab in metrics], fontsize=9)
    ax[1].legend(fontsize=8, ncol=2, title="method")
    ax[1].grid(True, axis="y", alpha=0.25)
    fig.suptitle("Molecular-window methods on L1 2026 (24 instruments; calipso dropped) — precision", fontsize=13)
    fig.tight_layout()
    fig.savefig(OUT / "method_precision_l1_2026.png", dpi=140, bbox_inches="tight")
    plt.close(fig)
    return agg, order


def fig_precision_by_group(rows):
    """σ_SD per method, one panel per instrument group (CHM15k / Mini-MPL / CL61)."""
    fig, axes = plt.subplots(1, 3, figsize=(17, 5.2), squeeze=False)
    for ax, grp in zip(axes[0], GROUPS):
        agg = {m: agg_by(rows, m, grp) for m in ALL_METHODS}
        order = sorted(ALL_METHODS, key=lambda m: agg[m]["sigma_sd"] if np.isfinite(agg[m]["sigma_sd"]) else 1e9)
        x = np.arange(len(order))
        ax.bar(x, [agg[m]["sigma_sd"] for m in order], color=[METHOD_COLORS[m] for m in order])
        for xi, m in zip(x, order):
            v = agg[m]["valid_pct"]
            if np.isfinite(v):
                ax.text(xi, 0.5, f"{v:.0f}%", ha="center", va="bottom", fontsize=7, color="white", rotation=90)
        n_inst = len({r["inst"] for r in rows if r["group"] == grp})
        ax.set_title(f"{grp}  (n={n_inst})")
        ax.set_xticks(x); ax.set_xticklabels([METHOD_LABEL[m].split(" ")[0] for m in order], rotation=25, ha="right", fontsize=8)
        ax.grid(True, axis="y", alpha=0.25)
    axes[0][0].set_ylabel("σ_SD (% of median C_L)\nsuccessive-difference precision")
    fig.suptitle("Per-group method precision (σ_SD, lower = better; bar label = valid%) — L1 2026", fontsize=13)
    fig.tight_layout()
    fig.savefig(OUT / "method_precision_by_group.png", dpi=140, bbox_inches="tight")
    plt.close(fig)


def fig_per_instrument(rows):
    """Per-instrument σ_SD and yield, using the recommended (best-precise, well-covered)
    method per instrument; bar colour = type, bar label = method initial + # valid nights."""
    labels = [m["label"] for m in MANIFEST if (OUT / f"results_{m['label']}.json").exists()]
    recs = [(lab, *recommended(rows, lab)) for lab in labels]
    recs = [(lab, m, r) for lab, m, r in recs if r is not None and np.isfinite(r["sigma_sd"])]
    recs.sort(key=lambda t: (GROUPS.index(t[2]["group"]), t[2]["sigma_sd"]))
    x = np.arange(len(recs))
    cols = [GROUP_COLORS[r["group"]] for _, _, r in recs]
    fig, ax = plt.subplots(2, 1, figsize=(16, 9))
    ax[0].bar(x, [r["sigma_sd"] for _, _, r in recs], color=cols)
    ax[0].set_ylabel("σ_SD (% of median C_L)")
    ax[0].set_title("Per-instrument night-to-night precision σ_SD (recommended method per instrument) — L1 2026")
    ax[0].grid(True, axis="y", alpha=0.25)
    ax[0].set_xticks(x); ax[0].set_xticklabels([])
    for xi, (_, m, r) in zip(x, recs):
        ax[0].text(xi, r["sigma_sd"] + 0.2, f"{m[0]}·{r['n_valid']}", ha="center", va="bottom",
                   fontsize=6.2, color="0.25")
    ax[1].bar(x, [r["valid_pct"] for _, _, r in recs], color=cols)
    ax[1].set_ylabel("valid nights (% of fit-nights)")
    ax[1].set_title("Per-instrument yield of the recommended method")
    ax[1].grid(True, axis="y", alpha=0.25)
    ax[1].set_xticks(x); ax[1].set_xticklabels([lab for lab, _, _ in recs], rotation=55, ha="right", fontsize=7.5)
    handles = [plt.Line2D([], [], marker="s", ls="", color=GROUP_COLORS[g], label=g) for g in GROUPS]
    ax[0].legend(handles=handles, title="instrument type", fontsize=9, ncol=3, loc="upper left")
    fig.suptitle("Per-instrument Rayleigh-calibration variability — L1 2026 "
                 "(bar label = method initial · # valid nights; o=optimal, e=earlinet, i=improved)", fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.97])
    fig.savefig(OUT / "per_instrument_variability.png", dpi=140, bbox_inches="tight")
    plt.close(fig)


def fig_timeseries(series, meta):
    """Small-multiples: per-instrument C_L time series (headline + production methods)."""
    insts = [m["label"] for m in MANIFEST if (OUT / f"results_{m['label']}.json").exists()]
    n = len(insts); ncol = 6; nrow = int(np.ceil(n / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(3.6 * ncol, 2.9 * nrow), squeeze=False)
    axA = axes.ravel()
    for ax, label in zip(axA, insts):
        grp = meta[label]["group"]
        for meth, lw, ms in [("eprof_v1.2", 0.7, 2.0), (HEADLINE, 1.1, 3.0)]:
            ds_ok, cls = series.get((label, meth), ([], []))
            if cls:
                dts = [datetime.strptime(d, "%Y%m%d") for d in ds_ok]
                ax.plot(dts, cls, "-o", ms=ms, lw=lw, color=METHOD_COLORS[meth], label=meth)
        # robust y-range from headline+improved
        allc = [c for meth in ("eprof_v1.2", HEADLINE) for c in series.get((label, meth), ([], []))[1]]
        if len(allc) >= 3:
            lo, hi = np.nanpercentile(allc, [3, 97]); pad = 0.6 * (hi - lo) + 1e-30
            ax.set_ylim(max(0, lo - pad), hi + pad)
        ax.set_title(label, fontsize=8.5, color=GROUP_COLORS[grp])
        ax.grid(True, alpha=0.25)
        ax.xaxis.set_major_formatter(mdates.DateFormatter("%m"))
        ax.tick_params(labelsize=6)
    for ax in axA[n:]:
        ax.axis("off")
    handles = [plt.Line2D([], [], color=METHOD_COLORS[m], marker="o", lw=1.0,
                          label=METHOD_LABEL[m]) for m in ("eprof_v1.2", HEADLINE)]
    fig.legend(handles=handles, loc="upper center", ncol=2, fontsize=10, frameon=False, bbox_to_anchor=(0.5, 1.005))
    fig.suptitle("Calibration-constant time series per instrument — L1 2026 (title colour = type)", y=1.02, fontsize=13)
    fig.tight_layout(rect=[0, 0, 1, 0.97])
    fig.savefig(OUT / "timeseries_l1_2026.png", dpi=135, bbox_inches="tight")
    plt.close(fig)


# ---------------------------------------------------------------------------
# Tables + summary
# ---------------------------------------------------------------------------
def write_tables(rows):
    agg_all = {m: agg_by(rows, m) for m in ALL_METHODS}
    order = sorted(ALL_METHODS, key=lambda m: agg_all[m]["sigma_sd"] if np.isfinite(agg_all[m]["sigma_sd"]) else 1e9)
    L = []
    L.append("### Per-method precision (mean over all 24 instruments)\n\n")
    L.append("method | valid_% | σ_night% | **σ_SD%** | σ_detrend% | σ_within-month% | rob_CV% | CV%\n")
    L.append("---|--:|--:|--:|--:|--:|--:|--:\n")
    for m in order:
        v = agg_all[m]
        L.append(f"{METHOD_LABEL[m]} | {v['valid_pct']:.0f} | {v['sigma_night']:.1f} | "
                 f"**{v['sigma_sd']:.1f}** | {v['sigma_dt']:.1f} | {v['sigma_im']:.1f} | "
                 f"{v['rob_cv']:.1f} | {v['cv']:.0f}\n".replace("nan", "–"))
    L.append(f"\n**Most precise (lowest σ_SD): `{order[0]}`.**\n\n")

    for grp in GROUPS:
        ag = {m: agg_by(rows, method=m, group=grp) for m in ALL_METHODS}
        og = sorted(ALL_METHODS, key=lambda m: ag[m]["sigma_sd"] if np.isfinite(ag[m]["sigma_sd"]) else 1e9)
        n_inst = len({r["inst"] for r in rows if r["group"] == grp})
        L.append(f"### {grp} (n={n_inst}) — per-method precision\n\n")
        L.append("method | valid_% | σ_night% | **σ_SD%** | σ_detrend% | σ_within-month% | rob_CV%\n")
        L.append("---|--:|--:|--:|--:|--:|--:\n")
        for m in og:
            v = ag[m]
            L.append(f"{METHOD_LABEL[m]} | {v['valid_pct']:.0f} | {v['sigma_night']:.1f} | "
                     f"**{v['sigma_sd']:.1f}** | {v['sigma_dt']:.1f} | {v['sigma_im']:.1f} | "
                     f"{v['rob_cv']:.1f}\n".replace("nan", "–"))
        L.append("\n")

    # Per-instrument variability with the recommended (best-precise, well-covered) method.
    labels = sorted({r["inst"] for r in rows},
                    key=lambda lab: (GROUPS.index(next(r["group"] for r in rows if r["inst"] == lab)),
                                     next(r["site"] for r in rows if r["inst"] == lab)))
    L.append("### Per-instrument variability — recommended method per instrument\n\n")
    L.append("The recommended method is the most precise quality method (preference "
             "optimal → earlinet → improved) with ≥6 valid nights, so every instrument is "
             "characterised even where the strictest method yields too few nights.\n\n")
    L.append("instrument | type | fit-nights | method | valid | valid_% | σ_night% | **σ_SD%** | σ_detrend% | σ_within-month% | rob_CV%\n")
    L.append("---|---|--:|:--|--:|--:|--:|--:|--:|--:|--:\n")
    for lab in labels:
        m, r = recommended(rows, lab)
        if r is None:
            r0 = next(x for x in rows if x["inst"] == lab)
            L.append(f"{lab} | {r0['type']} | {r0['n_fit']} | – | 0 | 0 | – | **–** | – | – | –\n")
            continue
        L.append(f"{lab} | {r['type']} | {r['n_fit']} | {m} | {r['n_valid']} | {r['valid_pct']:.0f} | "
                 f"{r['sigma_night']:.1f} | **{r['sigma_sd']:.1f}** | {r['sigma_dt']:.1f} | "
                 f"{r['sigma_im']:.1f} | {r['rob_cv']:.1f}\n".replace("nan", "–"))

    # Compact σ_SD-by-method comparison (the three quality methods).
    L.append("\n### Per-instrument σ_SD by quality method (valid nights in parentheses)\n\n")
    L.append("instrument | type | E-PROF v2 | EARLINET | E-PROF v1.2\n---|---|--:|--:|--:\n")
    by_im = {(r["inst"], r["method"]): r for r in rows}
    def cell(lab, m):
        r = by_im.get((lab, m))
        if r is None or not np.isfinite(r["sigma_sd"]):
            return f"– ({r['n_valid'] if r else 0})"
        return f"{r['sigma_sd']:.1f} ({r['n_valid']})"
    for lab in labels:
        typ = next(r["type"] for r in rows if r["inst"] == lab)
        L.append(f"{lab} | {typ} | {cell(lab,'eprof_v2')} | {cell(lab,'earlinet')} | {cell(lab,'eprof_v1.2')}\n")
    (OUT / "metrics_tables.md").write_text("".join(L), encoding="utf-8")

    # summary json (key numbers for the prose)
    summary = dict(
        n_instruments=len({r["inst"] for r in rows}),
        n_by_group={g: len({r["inst"] for r in rows if r["group"] == g}) for g in GROUPS},
        ranking_sigma_sd=order,
        agg_all={m: agg_all[m] for m in ALL_METHODS},
        agg_by_group={g: {m: agg_by(rows, method=m, group=g) for m in ALL_METHODS} for g in GROUPS},
        headline=HEADLINE,
        per_instrument={r["inst"]: {k: r[k] for k in ("type", "group", "n_fit", "n_valid",
                        "valid_pct", "sigma_night", "sigma_sd", "sigma_dt", "sigma_im", "rob_cv", "cv")}
                        for r in rows if r["method"] == HEADLINE},
    )
    rec = {}
    for lab in {r["inst"] for r in rows}:
        m, r = recommended(rows, lab)
        if r is not None:
            rec[lab] = {"method": m, **{k: r[k] for k in ("type", "group", "n_fit", "n_valid",
                        "valid_pct", "sigma_night", "sigma_sd", "sigma_dt", "sigma_im", "rob_cv")}}
    summary["per_instrument_recommended"] = rec
    (OUT / "metrics_summary.json").write_text(json.dumps(summary, indent=1), encoding="utf-8")
    return order, agg_all


def main():
    rows, series, meta = load_all()
    if not rows:
        print("No results found — run run_l1_2026_variability.py first.")
        return
    agg, order = fig_method_precision(rows)
    fig_precision_by_group(rows)
    fig_per_instrument(rows)
    fig_timeseries(series, meta)
    order2, agg_all = write_tables(rows)
    print("Most precise (σ_SD), all instruments:", order2[0])
    for m in order2:
        v = agg_all[m]
        print(f"  {m:9s} valid={v['valid_pct']:.0f}% σ_SD={v['sigma_sd']:.1f}% "
              f"σ_dt={v['sigma_dt']:.1f}% σ_im={v['sigma_im']:.1f}% rob_CV={v['rob_cv']:.1f}%")
    print("Saved figures + metrics_tables.md + metrics_summary.json to", OUT)
    print("METRICS_DONE")


if __name__ == "__main__":
    main()
