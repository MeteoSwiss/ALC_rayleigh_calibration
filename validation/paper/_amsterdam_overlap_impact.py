"""
_amsterdam_overlap_impact.py — what does the temperature-dependent overlap correction (Hervo et
al. 2016) actually do to the four co-located Amsterdam CHM15k units?

Amsterdam is the only site with four unit-specific overlap models side by side, so it isolates
the correction's unit-to-unit character. Four probes, one landscape figure:

  (a) the correction factor 1 + (a(z)*T + b(z))/100 per unit, evaluated at each unit's OBSERVED
      internal-temperature median (solid) and p10-p90 span (shading) — how different the four
      models are, and how much the temperature term moves them;
  (b) median beta_att profiles over common hours, corrected (solid) vs uncorrected (dashed),
      zoomed 0-1200 m — the correction in product space;
  (c) the realised impact per unit: median of 100*(corrected/uncorrected - 1) per gate — this is
      the correction actually applied to the science stream (multiplicative, so it equals the
      effective factor - 1);
  (d) does the correction tighten the four-way agreement in the NEAR RANGE? median relative bias
      of B/C/D vs unit A over 300-700 m (below the 500 m analysis floor where the correction
      lives), with and without the correction.

Everything uses the paper pipeline pieces (read_l1 -> channel_beta twice -> grid_and_stats), so
the matrices are exactly the validation's science stream. Numbers go to
discrepancy_analysis.json under "amsterdam_overlap".
Usage: python -m validation.paper._amsterdam_overlap_impact
"""
from __future__ import annotations
import json
import sys
import warnings
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

OUT = Path("C:/DATA/Projects/202606_E-PROFILE_calibration/figs_paper_validation/paper_python")
RESULTS = OUT / "discrepancy_analysis.json"
ZPLOT = 1200.0                 # the correction is a no-op above ~720 m; plot with margin
NEAR_BAND = (300.0, 700.0)     # near-range agreement band, below the 500 m analysis floor


def _medrel(X, A, gates):
    m = np.isfinite(X[:, gates]) & np.isfinite(A[:, gates]) & (A[:, gates] > 0)
    x, a = X[:, gates][m], A[:, gates][m]
    return float(100 * np.median((x - a) / a)) if x.size > 50 else np.nan


def main():
    warnings.filterwarnings("ignore")
    from validation.paper import run_paper_validation as RPV
    from validation.paper import intercompare as IC
    from validation.paper import overlap as OV
    from validation.paper import figures as FIG
    from validation.paper.calib_benchmark import BENCHMARK

    st = BENCHMARK["amsterdam"]
    items, temps, models = [], {}, {}
    for ch in st["channels"]:
        print(f"reading + processing unit {ch['ident']} ...", flush=True)
        l1 = IC.read_l1(ch["wmo"], ch["ident"], st["start"], st["end"])
        on = RPV.channel_beta(l1, ch, 1064.0, apply_wv=True)
        off = RPV.channel_beta(l1, ch, 1064.0, apply_wv=True, apply_overlap=False)
        items.append(dict(ch=ch, l1=l1, **on, beta_scr_noovl=off["beta_scr"]))
        temps[ch["ident"]] = np.asarray(l1["temp_int"], "f8")
        models[ch["ident"]] = OV.load_overlap_model(ch["wmo"], ch["ident"])
    R = RPV.grid_and_stats(items, iref=0)
    z = np.asarray(R["altGrid"]) - R["station"]["altitude"]
    zp = (z >= 0) & (z <= ZPLOT)
    cols = FIG.channel_colors(R["channels"])
    labels = [c["label"] for c in R["channels"]]

    res = {"near_band_m": list(NEAR_BAND)}
    fig, axes = plt.subplots(1, 4, figsize=(19, 6.5))

    # (a) model correction factors at the observed temperatures --------------------------------
    ax = axes[0]
    zg = np.arange(0.0, ZPLOT + 1, 15.0)
    for k, it in enumerate(items):
        ident = it["ch"]["ident"]
        m = models[ident]
        T = temps[ident]
        t50, t10, t90 = np.nanpercentile(T, [50, 10, 90])
        a = np.interp(zg, m["range"], m["a"], left=m["a"][0], right=0.0)
        b = np.interp(zg, m["range"], m["b"], left=m["b"][0], right=0.0)
        f50 = 1 + (a * t50 + b) / 100.0
        f10 = 1 + (a * t10 + b) / 100.0
        f90 = 1 + (a * t90 + b) / 100.0
        ax.plot(f50, zg, "-", color=cols[k], lw=1.8,
                label=f"{labels[k]}  (T$_{{med}}$={t50:.0f}°C)")
        ax.fill_betweenx(zg, np.minimum(f10, f90), np.maximum(f10, f90), color=cols[k], alpha=0.15)
        res.setdefault("factor_at_Tmed", {})[ident] = {
            "T_med_C": float(t50), "T_p10_C": float(t10), "T_p90_C": float(t90),
            "factor_200m": float(np.interp(200, zg, f50)),
            "factor_350m": float(np.interp(350, zg, f50)),
            "factor_500m": float(np.interp(500, zg, f50)),
            "factor_700m": float(np.interp(700, zg, f50)), "model_file": m["file"]}
    ax.axvline(1.0, color="k", lw=0.8)
    ax.set_xlabel("overlap correction factor  1 + Dif/100")
    ax.set_ylabel("altitude a.g.l. [m]"); ax.set_ylim(0, ZPLOT)
    ax.grid(alpha=0.3); ax.legend(fontsize=7, loc="upper right")
    ax.set_title("(a) Model factor at observed T (band: p10–p90 T)", fontsize=10)

    # (b) median profiles, corrected vs uncorrected, common hours ------------------------------
    ax = axes[1]
    have = [np.any(np.isfinite(R["beta"][k][:, zp]), axis=1) for k in range(len(items))]
    common = np.logical_and.reduce(have)
    for k in range(len(items)):
        B, Bn = R["beta"][k][common][:, zp], R["beta_noovl"][k][common][:, zp]
        with np.errstate(all="ignore"):
            ax.plot(np.nanmedian(B, axis=0), z[zp], "-", color=cols[k], lw=1.8, label=labels[k])
            ax.plot(np.nanmedian(Bn, axis=0), z[zp], "--", color=cols[k], lw=1.2)
    ax.set_xlabel(r"median $\beta_{att}$ [Mm$^{-1}$ sr$^{-1}$]  (solid corr / dashed uncorr)")
    ax.set_ylim(0, ZPLOT); ax.grid(alpha=0.3); ax.legend(fontsize=7)
    ax.set_title(f"(b) Median profile, N={int(common.sum())} common hours", fontsize=10)

    # (c) realised impact on the science stream -------------------------------------------------
    ax = axes[2]
    for k in range(len(items)):
        with np.errstate(all="ignore"):
            ratio = R["beta"][k][:, zp] / R["beta_noovl"][k][:, zp]
            imp = 100 * (np.nanmedian(ratio, axis=0) - 1)
        ax.plot(imp, z[zp], "-", color=cols[k], lw=1.8, label=labels[k])
        ident = items[k]["ch"]["ident"]
        res.setdefault("impact_pct", {})[ident] = {
            "at_200m": float(np.interp(200, z[zp], imp)), "at_350m": float(np.interp(350, z[zp], imp)),
            "at_500m": float(np.interp(500, z[zp], imp)), "at_700m": float(np.interp(700, z[zp], imp))}
    ax.axvline(0, color="k", lw=0.8)
    for zb in NEAR_BAND:
        ax.axhline(zb, ls=":", color="k", lw=0.8)
    ax.set_xlabel(r"impact on $\beta_{att}$: 100·(corr/uncorr − 1) [%]")
    ax.set_ylim(0, ZPLOT); ax.grid(alpha=0.3); ax.legend(fontsize=7)
    ax.set_title("(c) Realised impact per unit (dotted: 300–700 m band)", fontsize=10)

    # (d) four-way near-range agreement with / without the correction ---------------------------
    ax = axes[3]
    gates = (z >= NEAR_BAND[0]) & (z <= NEAR_BAND[1])
    A_c, A_n = R["beta"][0], R["beta_noovl"][0]
    names, v_corr, v_unc = [], [], []
    for k in range(1, len(items)):
        names.append(labels[k].replace("CHM15k ", ""))
        v_corr.append(_medrel(R["beta"][k], A_c, gates))
        v_unc.append(_medrel(R["beta_noovl"][k], A_n, gates))
    x = np.arange(len(names))
    ax.bar(x - 0.18, v_unc, 0.36, color="0.65", label="uncorrected")
    ax.bar(x + 0.18, v_corr, 0.36, color="#2ca02c", label="overlap-corrected")
    for xi, (vu, vc) in enumerate(zip(v_unc, v_corr)):
        ax.text(xi - 0.18, vu, f"{vu:+.1f}", ha="center", va="bottom", fontsize=8)
        ax.text(xi + 0.18, vc, f"{vc:+.1f}", ha="center", va="bottom", fontsize=8)
    ax.axhline(0, color="k", lw=0.8)
    ax.set_xticks(x); ax.set_xticklabels([f"unit {n}" for n in names])
    ax.set_ylabel(f"med relbias vs unit A, {int(NEAR_BAND[0])}–{int(NEAR_BAND[1])} m [%]")
    ax.grid(alpha=0.3, axis="y"); ax.legend(fontsize=8)
    ax.set_title("(d) Near-range four-way agreement", fontsize=10)
    res["near_range_agreement"] = {n: {"uncorrected": u, "corrected": c}
                                   for n, u, c in zip(names, v_unc, v_corr)}

    print("impact @350m:", {k: round(v["at_350m"], 1) for k, v in res["impact_pct"].items()}, flush=True)
    print("near-range agreement:", res["near_range_agreement"], flush=True)

    fig.suptitle("Amsterdam: impact of the temperature-dependent overlap correction on the four"
                 " co-located CHM15k units", fontweight="bold", fontsize=13)
    fig.tight_layout(rect=(0, 0, 1, 0.92))
    fig.savefig(OUT / "fig_amsterdam_overlap_impact.png", dpi=160)
    plt.close(fig)
    print("-> fig_amsterdam_overlap_impact.png", flush=True)

    J = {}
    if RESULTS.is_file():
        try:
            J = json.loads(RESULTS.read_text())
        except Exception:
            J = {}
    J["amsterdam_overlap"] = res
    RESULTS.write_text(json.dumps(J, indent=1))
    print("AMSTERDAM_OVERLAP_DONE", flush=True)


if __name__ == "__main__":
    try:
        sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    except Exception:
        pass
    main()
