"""
precision_longrun.py — drift-INSENSITIVE precision metrics for the molecular-window methods,
from the long-run per-instrument JSON checkpoints (no re-run).

The plain CV (std/mean) mixes measurement PRECISION with real seasonal + long-term (laser
ageing / window) drift, so it unfairly penalises a precise instrument that simply has a
seasonal cycle. We instead report metrics that remove slow drift:

  valid_pct   % of sampled nights that yield a valid calibration            (yield)
  sigma_night avg single-night spread = in-window std(signal/molecular)/CL  (within-night)
  sigma_SD    successive-difference precision = robust |ΔCL| between
              time-ordered consecutive calibrations / sqrt(2) / median      (night-to-night)
  sigma_dt    detrended residual = robust scatter of CL minus its ~2-month
              rolling median                                                (drift-removed)
  sigma_im    within-month scatter = robust scatter of CL pooled within
              each calendar month (season ~ const within a month)          (drift-removed)
  CV          classic std/mean, kept for reference (precision + drift)
All scatters are relative (% of the median CL). sigma_SD / sigma_dt / sigma_im are the
recommended precision metrics; CV is NOT.
"""
from __future__ import annotations
import json, glob, os
from collections import defaultdict
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from compare_molecular_methods import METHODS_DISPLAY, METHOD_COLORS, METHOD_LABEL

OUT = Path("C:/DATA/Projects/202606_E-PROFILE_calibration/figs_paper_validation/molecular_methods_longrun")
QC = 15.0
I_OK, I_CL, I_CLERR, I_REL = 0, 1, 2, 3


def _robust_rel(x, med):
    """1.4826*MAD(x) / |med| in %."""
    if len(x) < 3 or med == 0:
        return np.nan
    return float(1.4826 * np.median(np.abs(np.asarray(x) - np.median(x))) / abs(med) * 100)


def _rolling_median(y, k=9):
    n = len(y)
    h = k // 2
    out = np.empty(n)
    for i in range(n):
        out[i] = np.median(y[max(0, i - h):min(n, i + h + 1)])
    return out


def precision(dslist, cls, clerrs):
    """dslist time-ordered (YYYYMMDD), cls/clerrs aligned. Returns the metric dict."""
    cls = np.asarray(cls, float)
    med = np.median(cls)
    out = dict(n=len(cls), cv=np.nan, sigma_night=np.nan, sigma_sd=np.nan,
               sigma_dt=np.nan, sigma_im=np.nan)
    if len(cls) < 4 or med <= 0:
        if len(cls):
            rel = np.asarray(clerrs) / np.abs(cls)
            out["sigma_night"] = float(np.median(rel[np.isfinite(rel)]) * 100) if np.any(np.isfinite(rel)) else np.nan
        return out
    out["cv"] = float(np.std(cls) / abs(np.mean(cls)) * 100)
    rel = np.asarray(clerrs) / np.abs(cls)
    out["sigma_night"] = float(np.median(rel[np.isfinite(rel)]) * 100)
    # successive difference (von Neumann), robust: diff std of i.i.d. noise = sqrt(2)*sigma
    diffs = np.diff(cls)
    out["sigma_sd"] = float(1.4826 * np.median(np.abs(diffs)) / np.sqrt(2) / abs(med) * 100)
    # detrended residual (remove ~2-month rolling median)
    resid = cls - _rolling_median(cls, 9)
    out["sigma_dt"] = _robust_rel(resid, med)
    # within-month pooled residual
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


def main():
    files = sorted(glob.glob(str(OUT / "results_*.json")))
    rows = []
    for fp in files:
        label = os.path.basename(fp)[len("results_"):-len(".json")]
        itype = "Mini-MPL" if "MPL" in label else "CHM15k"
        data = json.load(open(fp))
        n_nights = len(data)
        order = sorted(data)
        for m in METHODS_DISPLAY:
            ds_ok, cls, errs = [], [], []
            for ds in order:
                w = data[ds][m]
                if w[I_OK] and np.isfinite(w[I_REL]) and w[I_REL] <= QC and np.isfinite(w[I_CL]) and w[I_CL] > 0:
                    ds_ok.append(ds); cls.append(w[I_CL]); errs.append(w[I_CLERR])
            p = precision(ds_ok, cls, errs)
            p.update(inst=label, itype=itype, method=m,
                     valid_pct=100.0 * len(cls) / n_nights if n_nights else np.nan)
            rows.append(p)

    # aggregate per method (mean over instruments)
    agg = {}
    for m in METHODS_DISPLAY:
        mr = [r for r in rows if r["method"] == m]
        def mn(k):
            v = [r[k] for r in mr if np.isfinite(r[k])]
            return float(np.mean(v)) if v else np.nan
        agg[m] = {k: mn(k) for k in ("valid_pct", "sigma_night", "sigma_sd", "sigma_dt", "sigma_im", "cv")}
    ranking = sorted(agg.items(), key=lambda kv: (kv[1]["sigma_sd"] if np.isfinite(kv[1]["sigma_sd"]) else 1e9))

    L = ["# Long-run precision — drift-insensitive metrics (full archive, 14 sites)\n\n",
         "CV mixes precision with seasonal + long-term drift; the metrics below remove slow "
         "drift. All scatters are % of the median CL. `sigma_SD` (successive-difference), "
         "`sigma_dt` (detrended), `sigma_im` (within-month) are the precision metrics; CV is "
         "kept only for reference.\n\n",
         "## Per method (mean over the 14 instruments)\n\n",
         "method | valid_% | sigma_night% | **sigma_SD%** | sigma_dt% | sigma_im% | CV% (ref)\n",
         "---|---|---|---|---|---|---\n"]
    for m, v in ranking:
        L.append(f"{m} | {v['valid_pct']:.0f} | {v['sigma_night']:.1f} | **{v['sigma_sd']:.1f}** | "
                 f"{v['sigma_dt']:.1f} | {v['sigma_im']:.1f} | {v['cv']:.0f}\n".replace("nan", "-"))
    L.append(f"\n**Most precise (lowest sigma_SD): `{ranking[0][0]}`.** "
             "Note CV ≫ sigma_SD for every method — most of the CV was real drift, not noise.\n\n")
    L.append("## Per instrument × method\n\n")
    L.append("inst | type | method | valid_% | sigma_night% | sigma_SD% | sigma_dt% | sigma_im% | CV%\n")
    L.append("---|---|---|---|---|---|---|---|---\n")
    for r in rows:
        L.append(f"{r['inst']} | {r['itype']} | {r['method']} | {r['valid_pct']:.0f} | "
                 f"{r['sigma_night']:.1f} | {r['sigma_sd']:.1f} | {r['sigma_dt']:.1f} | "
                 f"{r['sigma_im']:.1f} | {r['cv']:.0f}\n".replace("nan", "-"))
    (OUT / "precision_longrun.md").write_text("".join(L), encoding="utf-8")
    print("saved precision_longrun.md  -> most precise (sigma_SD):", ranking[0][0])
    for m, v in ranking:
        print(f"  {m:9s} valid={v['valid_pct']:.0f}% night={v['sigma_night']:.1f}% "
              f"SD={v['sigma_sd']:.1f}% detrend={v['sigma_dt']:.1f}% inmonth={v['sigma_im']:.1f}% CV={v['cv']:.0f}%")

    # figure: the four drift-insensitive metrics + CV per method
    fig, ax = plt.subplots(1, 2, figsize=(16, 6))
    mlist = [m for m, _ in ranking]
    x = np.arange(len(mlist))
    for metric, lab, a in [("valid_pct", "valid days (%)", ax[0])]:
        a.bar(x, [agg[m][metric] for m in mlist], color=[METHOD_COLORS[m] for m in mlist])
        a.set_title("Yield: % of nights with a valid calibration")
        a.set_ylabel("valid nights (%)")
        a.set_xticks(x); a.set_xticklabels([METHOD_LABEL[m] for m in mlist], rotation=30, ha="right", fontsize=8)
        a.grid(True, axis="y", alpha=0.25)
    # Right panel grouped by METRIC (not by method); one bar per calibration method,
    # colour-coded by method (same colour code as the left panel). The four drift-
    # insensitive metrics share the left axis; CV (much larger, reference only) sits on a
    # twin axis so it stays visible without squashing the precision metrics.
    metrics_main = [("sigma_night", "σ_night\n(within-night)"), ("sigma_sd", "σ_SD\n(successive)"),
                    ("sigma_dt", "σ_detrend"), ("sigma_im", "σ_within-month")]
    nm = len(mlist)
    bw = 0.8 / nm
    gx = np.arange(len(metrics_main) + 1)          # last group (CV) is on the twin axis
    ax2 = ax[1].twinx()
    for j, m in enumerate(mlist):
        off = (j - (nm - 1) / 2) * bw
        ax[1].bar(gx[:len(metrics_main)] + off, [agg[m][k] for k, _ in metrics_main], bw,
                  color=METHOD_COLORS[m], label=METHOD_LABEL[m])
        ax2.bar(gx[-1] + off, agg[m]["cv"], bw, color=METHOD_COLORS[m])
    ax[1].axvline(gx[-1] - 0.5, color="0.6", lw=1, ls="--")
    ax[1].set_title("Precision by metric (drift-removed), colour = method — lower = more precise")
    ax[1].set_ylabel("σ (% of median CL) — drift-insensitive")
    ax2.set_ylabel("CV (% of median CL) — reference (drift + noise)")
    ax[1].set_xticks(gx)
    ax[1].set_xticklabels([lab for _, lab in metrics_main] + ["CV\n(right axis)"], fontsize=8)
    ax[1].set_ylim(0, max(30, np.nanmax([agg[m][k] for m in mlist for k, _ in metrics_main]) * 1.25))
    ax2.set_ylim(0, np.nanmax([agg[m]["cv"] for m in mlist]) * 1.2)
    ax[1].legend(fontsize=8, ncol=2, title="calibration method")
    ax[1].grid(True, axis="y", alpha=0.25)
    fig.suptitle("Calibration precision metrics (drift-insensitive) — full archive", fontsize=13)
    fig.tight_layout()
    fig.savefig(OUT / "precision_longrun.png", dpi=140, bbox_inches="tight")
    plt.close(fig)
    print("saved precision_longrun.png")


if __name__ == "__main__":
    main()
