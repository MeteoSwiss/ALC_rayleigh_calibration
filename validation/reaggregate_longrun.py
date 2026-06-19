"""Re-aggregate the long-run per-instrument JSON checkpoints with a ROBUST night-to-night
scatter (MAD-based), since std/mean CV is dominated by rare outlier nights (which the Kalman
would smooth anyway). Writes a robust ranking table + a summary figure."""
from __future__ import annotations
import json, glob, os
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from compare_molecular_methods import METHODS_DISPLAY, METHOD_COLORS, METHOD_LABEL

OUT = Path("C:/DATA/Projects/202606_E-PROFILE_calibration/figs_paper_validation/molecular_methods_longrun")
QC = 15.0
# JSON layout per night/method: [ok, cl, cl_err, rel_error, r2, temporal_cv, scat, start, end, center]
I_OK, I_CL, I_REL, I_R2, I_TCV = 0, 1, 3, 4, 5


def robust_cv(cls):
    cls = np.asarray([c for c in cls if np.isfinite(c) and c > 0], float)
    if cls.size < 3:
        return np.nan, np.nan
    med = np.median(cls)
    mad = np.median(np.abs(cls - med)) * 1.4826
    std_cv = float(np.std(cls) / abs(np.mean(cls)) * 100)
    rob_cv = float(mad / abs(med) * 100) if med != 0 else np.nan
    return rob_cv, std_cv


def main():
    files = sorted(glob.glob(str(OUT / "results_*.json")))
    rows = []
    for fp in files:
        label = os.path.basename(fp)[len("results_"):-len(".json")]
        data = json.load(open(fp))
        itype = "Mini-MPL" if "MPL" in label else "CHM15k"
        for m in METHODS_DISPLAY:
            cls, r2s, tcvs = [], [], []
            n_cal = 0
            for ds, rec in data.items():
                w = rec[m]
                if w[I_OK] and np.isfinite(w[I_REL]) and w[I_REL] <= QC:
                    n_cal += 1
                    cls.append(w[I_CL]); r2s.append(w[I_R2]); tcvs.append(w[I_TCV])
            rob, std = robust_cv(cls)
            rows.append(dict(inst=label, itype=itype, method=m, n_nights=len(data), n_cal=n_cal,
                             rob_cv=rob, std_cv=std,
                             med_r2=float(np.median(r2s)) if r2s else np.nan,
                             med_tcv=float(np.nanmedian(tcvs)) if tcvs else np.nan))
    # ranking (robust): frac * stability(robust CV, ref 20%) * cleanliness
    score = {}
    for m in METHODS_DISPLAY:
        mr = [r for r in rows if r["method"] == m]
        fr = np.mean([r["n_cal"] / r["n_nights"] for r in mr if r["n_nights"]])
        rc = np.nanmean([r["rob_cv"] for r in mr if np.isfinite(r["rob_cv"])])
        tc = np.nanmean([r["med_tcv"] for r in mr if np.isfinite(r["med_tcv"])])
        stab = 20.0 / max(rc, 20.0)
        clean = 1.0 / (1.0 + (tc if np.isfinite(tc) else 1.0))
        score[m] = dict(frac=fr, rob_cv=rc, tcv=tc, s=fr * stab * clean)
    ranking = sorted(score.items(), key=lambda kv: kv[1]["s"], reverse=True)

    L = ["# Long-run (full archive) — molecular-window methods, ROBUST aggregation\n\n",
         "14 instruments (10 CHM15k + 4 Mini-MPL), 5 nights/month over the full E-PROFILE L2 "
         "archive (~80-113 months each). `rob_CV` = robust night-to-night scatter "
         "(1.4826·MAD/median, %); `std_CV` = classic std/mean (outlier-sensitive, shown for "
         "reference). n_cal = nights calibrating through the pipeline.\n\n",
         "inst | type | method | nights | n_cal | rob_CV% | std_CV% | med_R2 | med_tcv\n",
         "---|---|---|---|---|---|---|---|---\n"]
    for r in rows:
        L.append(f"{r['inst']} | {r['itype']} | {r['method']} | {r['n_nights']} | {r['n_cal']} | "
                 f"{r['rob_cv']:.1f} | {r['std_cv']:.0f} | {r['med_r2']:.3f} | {r['med_tcv']:.2f}\n"
                 .replace("nan", "-"))
    L.append("\n## Ranking (mean over instruments, ROBUST CV)\n\n")
    L.append("method | calibrated-fraction | mean rob_CV% | mean temporal_cv | score\n---|---|---|---|---\n")
    for m, v in ranking:
        L.append(f"{m} | {v['frac']:.2f} | {v['rob_cv']:.1f} | {v['tcv']:.2f} | {v['s']:.3f}\n".replace("nan", "-"))
    L.append(f"\n**Best overall (robust, full archive): `{ranking[0][0]}`.**\n")
    (OUT / "ranking_robust_longrun.md").write_text("".join(L), encoding="utf-8")
    print("saved ranking_robust_longrun.md  -> best =", ranking[0][0])
    for m, v in ranking:
        print(f"  {m:9s} frac={v['frac']:.2f} robCV={v['rob_cv']:.0f}% tcv={v['tcv']:.2f} score={v['s']:.3f}")

    # summary figure: robust CV per method per instrument
    insts = sorted({r["inst"] for r in rows})
    fig, ax = plt.subplots(1, 2, figsize=(17, 6))
    x = np.arange(len(insts))
    nm = len(METHODS_DISPLAY)
    w = 0.8 / nm
    for k, m in enumerate(METHODS_DISPLAY):
        fr = [next((r["n_cal"] / r["n_nights"] for r in rows if r["inst"] == i and r["method"] == m), 0) for i in insts]
        rc = [next((r["rob_cv"] for r in rows if r["inst"] == i and r["method"] == m), np.nan) for i in insts]
        off = (k - (nm - 1) / 2) * w
        ax[0].bar(x + off, fr, w, color=METHOD_COLORS[m], label=METHOD_LABEL[m])
        ax[1].bar(x + off, rc, w, color=METHOD_COLORS[m])
    ax[0].set_title("Usable-night fraction (full archive)")
    ax[0].set_ylabel("fraction")
    ax[1].set_title("Robust night-to-night CV (MAD-based, lower=better)")
    ax[1].set_ylabel("robust CV (%)")
    ax[1].set_ylim(0, 120)
    for a in ax:
        a.set_xticks(x); a.set_xticklabels(insts, rotation=40, ha="right", fontsize=7)
        a.grid(True, axis="y", alpha=0.25)
    ax[0].legend(ncol=4, fontsize=7.5)
    fig.suptitle("Molecular-window methods over the full archive (14 sites)", fontsize=13)
    fig.tight_layout()
    fig.savefig(OUT / "summary_robust_longrun.png", dpi=140, bbox_inches="tight")
    plt.close(fig)
    print("saved summary_robust_longrun.png")


if __name__ == "__main__":
    main()
