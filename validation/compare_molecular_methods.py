"""
compare_molecular_methods.py

Compare the six molecular-window detection methods (main / improved / matlab / calipso /
earlinet / optimal) across several E-PROFILE sites and instruments, select the best, and
produce diagnostics:

  * DETAIL figures (a few representative instruments): vertical profiles with each method's
    window bracketed, and a time-height signal/molecular pcolor with the selected window
    centres + the cells the "optimal" method flagged as aerosol/cloud (and excluded);
  * STATS across all instruments/nights: per-method calibrated-night count, night-to-night
    lidar-constant CV, median R²/temporal-CV/rel-error -> a recommended method.

Figures + tables -> figs_paper_validation/molecular_methods/ in the MATLAB ALC repo.
"""
from __future__ import annotations
import os, sys, logging, warnings
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[_v] = "1"
try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass
from pathlib import Path
from datetime import datetime
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.dates as mdates

from rayleigh_calibration import calibrate_rayleigh, CalibrationOptions, InstrumentInfo, DataLevel
from rayleigh_calibration.config import InstrumentType
from rayleigh_calibration.rayleigh.molecular_methods import METHODS, compute_window_grid, select_molecular_window

OUT = Path("C:/DATA/Projects/202606_E-PROFILE_calibration/figs_paper_validation/molecular_methods")
OUT.mkdir(parents=True, exist_ok=True)
HALF = tuple(range(250, 2000, 240))
QC_THR = 15.0  # pipeline threshold_quality: rel_error > this -> flag -2

METHOD_COLORS = {"eprof_v10": "#000000", "main": "#7f7f7f", "improved": "#1f77b4",
                 "matlab": "#2ca02c", "calipso": "#ff7f0e", "earlinet": "#9467bd",
                 "optimal": "#d62728", "bellini": "#8c564b"}
# The three E-PROFILE production versions are renamed to their release numbers:
#   main -> v1.1 (sign-fixed legacy), improved -> v1.2 (production), optimal -> v2.
# v1.0 (eprof_v10) is the historical pre-a4e7140 sign-error calibration (run with the
# same 'main' window but options.sign_error_v10=True). matlab/calipso/earlinet/bellini
# are the external reference methods and keep their own names.
METHOD_LABEL = {"eprof_v10": "E-PROF v1.0", "main": "E-PROF v1.1",
                "improved": "E-PROF v1.2", "optimal": "E-PROF v2",
                "matlab": "MATLAB Auto_Calib_25", "calipso": "CALIPSO-type",
                "earlinet": "EARLINET/SCC-type", "bellini": "Bellini/ALICENET"}
# Ordered tuple for display/precision/longrun aggregation — includes eprof_v10 for historical
# comparison. METHODS (in molecular_methods.py) only lists live selectable methods.
METHODS_DISPLAY = ("eprof_v10", "main", "improved", "optimal",
                   "matlab", "calipso", "earlinet", "bellini")

T = InstrumentType
INSTRUMENTS = [
    dict(label="Payerne_CHM15k",   wmo="0-20000-0-06610", ident="A", itype=T.CHM15k, lat=46.8137, lon=6.9425, alt=505.0),
    dict(label="Payerne_CL31",     wmo="0-20000-0-06610", ident="B", itype=T.CL31,   lat=46.8137, lon=6.9425, alt=505.0),
    dict(label="Payerne_CL61",     wmo="0-20000-0-06610", ident="C", itype=T.CL61,   lat=46.8137, lon=6.9425, alt=505.0),
    dict(label="Amsterdam_CHM15k", wmo="0-20000-0-06240", ident="A", itype=T.CHM15k, lat=52.3170, lon=4.8037, alt=6.0),
    dict(label="EDT_CL51",         wmo="0-20008-0-EDT",   ident="A", itype=T.CL51,   lat=53.5500, lon=-114.1000, alt=776.0),
    dict(label="EDT_CL61",         wmo="0-20008-0-EDT",   ident="B", itype=T.CL61,   lat=53.5500, lon=-114.1000, alt=776.0),
]
DETAIL_LABELS = {"Payerne_CL61", "Amsterdam_CHM15k", "EDT_CL61"}
DETAIL_PREF = {"Payerne_CL61": ["20260312", "20260316", "20260328"]}  # the nights discussed
N_DETAIL = 3
# Sample every 2nd day of Feb+Mar 2026 (EDT has only Jan-Mar; bounded runtime on the slow A:
# drive — cloudy nights return fast before the fit, so cost ~ number of clear nights). The
# pinned Payerne detail nights are added so they are always processed.
_BASE = [f"2026{mm:02d}{d:02d}" for mm in (2, 3) for d in range(1, 32, 2)]
STATS_NIGHTS = sorted(set(_BASE) | {n for ns in DETAIL_PREF.values() for n in ns})


def base_options():
    o = CalibrationOptions.from_json(Path("options.json"))
    o.folder_root = Path("A:/E-PROFILE_L2_monthly")
    o.data_level = DataLevel.L2_MONTHLY
    o.molecular_source = "standard"
    o.plot_main = False
    o.plot_all = False
    o.folder_output = OUT
    return o


def run_methods(signal, p_mol, rng, stack):
    """Run all 6 methods on one prepared profile (optimal gets the stack for flagging)."""
    grid = compute_window_grid(signal, p_mol, rng, HALF, range_start_m=2000,
                               range_end_m=6000, increment_bins=8, signal_stack=stack)
    out = {}
    for m in METHODS:
        if m == "optimal":
            out[m] = select_molecular_window("optimal", signal, p_mol, rng, HALF, signal_stack=stack)
        else:
            out[m] = select_molecular_window(m, signal, p_mol, rng, HALF, grid=grid)
    return out


def calibrates(w):
    return bool(w.ok and np.isfinite(w.rel_error) and w.rel_error <= QC_THR)


def run_instrument(inst):
    o = base_options()
    info = InstrumentInfo(site_name=inst["label"], wmo_id=inst["wmo"], identifier=inst["ident"],
                          instrument_type=inst["itype"], latitude=inst["lat"],
                          longitude=inst["lon"], altitude=inst["alt"])
    per_night = {}        # ds -> {method: MethodWindow}
    detail = {}           # ds -> inputs (only for DETAIL instruments)
    is_detail = inst["label"] in DETAIL_LABELS
    pref = DETAIL_PREF.get(inst["label"])
    for ds in STATS_NIGHTS:
        fin = {}
        try:
            calibrate_rayleigh(ds, info, o, fit_inputs_out=fin)
        except Exception:
            continue
        if not fin:
            continue
        res = run_methods(fin["signal"], fin["p_mol"], fin["range_alc"], fin["signal_stack"])
        per_night[ds] = res
        if is_detail:
            if pref is not None:
                if ds in pref:
                    detail[ds] = fin
            elif len(detail) < N_DETAIL:
                detail[ds] = fin
    n_clear = len(per_night)
    print(f"  {inst['label']}: {n_clear} nights reached the fit (of {len(STATS_NIGHTS)} sampled)")
    return per_night, detail


# ----------------------------------------------------------------------------
# DETAIL figures
# ----------------------------------------------------------------------------
SHORT = {"eprof_v10": "v1.0", "main": "main", "improved": "impr", "matlab": "matl",
         "calipso": "cali", "earlinet": "earl", "optimal": "opt", "bellini": "bell"}


def _cl_panel(ax, res):
    """Bottom subplot: lidar calibration constant per method with error bars (only methods
    that actually calibrate through the pipeline)."""
    any_pt = False
    for k, m in enumerate(METHODS):
        w = res[m]
        if w.ok and np.isfinite(w.rel_error) and w.rel_error <= QC_THR and np.isfinite(w.cl):
            any_pt = True
            err = w.cl_err if np.isfinite(w.cl_err) else 0.0
            ax.errorbar([k], [w.cl], yerr=[err], fmt="o", ms=5, color=METHOD_COLORS[m], capsize=3, lw=1.3)
    ax.set_xticks(range(len(METHODS)))
    ax.set_xticklabels([SHORT[m] for m in METHODS], rotation=45, fontsize=6.5, ha="right")
    ax.set_xlim(-0.6, len(METHODS) - 0.4)
    ax.set_ylabel("C_L (proxy)", fontsize=8)
    ax.tick_params(axis="y", labelsize=7)
    ax.grid(True, axis="y", alpha=0.25)
    if not any_pt:
        ax.text(0.5, 0.5, "no calibration", ha="center", va="center",
                transform=ax.transAxes, color="0.5", fontsize=8)


def plot_profiles(label, detail, per_night):
    nights = list(detail.keys())
    if not nights:
        return
    n = len(nights)
    fig, axes = plt.subplots(2, n, figsize=(4.8 * n, 9.2), squeeze=False,
                             gridspec_kw={"height_ratios": [3, 1]})
    for col, ds in enumerate(nights):
        ax = axes[0][col]
        d = detail[ds]
        rng_km = d["range_alc"] / 1e3
        sig = d["signal"]
        pmol = d["p_mol"]
        res = per_night[ds]
        cls = [w.cl for w in res.values() if w.ok and np.isfinite(w.cl) and w.cl > 0]
        cl_ref = np.median(cls) if cls else np.nanmedian(sig / pmol)
        good = np.isfinite(sig) & (sig > 0) & (rng_km < 7)
        vis = sig[good]
        if vis.size:
            xlo = float(np.nanmin(vis))            # TRUE minimum positive signal -> never clipped
            xhi = float(np.nanpercentile(vis, 99.8))
        else:
            xlo, xhi = 1e-9, 1e-3
        ax.plot(sig[good], rng_km[good], color="0.2", lw=0.9)
        ax.plot(cl_ref * pmol, rng_km, color="#1f77b4", lw=1.3, ls="--")
        ax.set_xscale("log")
        ax.set_ylim(0, 7)
        ax.set_xlim(xlo * 0.5, xhi * 2.0)
        xs = np.geomspace(xhi * 1.5e-2, xhi * 0.8, len(METHODS))
        for k, m in enumerate(METHODS):
            w = res[m]
            if not w.ok:
                continue
            c = METHOD_COLORS[m]
            x = xs[k]
            ax.plot([x, x], [w.start_m / 1e3, w.end_m / 1e3], color=c, lw=2.8, solid_capstyle="butt")
            ax.plot([x], [w.center_m / 1e3], marker="o", color=c, ms=5)
            for yy in (w.start_m / 1e3, w.end_m / 1e3):
                ax.plot([x * 0.84, x * 1.19], [yy, yy], color=c, lw=1.2)
        ax.set_title(ds, fontsize=10)
        ax.set_xlabel("range-norm. signal (a.u.)")
        ax.grid(True, alpha=0.25, which="both")
        _cl_panel(axes[1][col], res)
    axes[0][0].set_ylabel("altitude AGL (km)")
    handles = [plt.Line2D([], [], color=METHOD_COLORS[m], lw=2.8, label=METHOD_LABEL[m]) for m in METHODS]
    handles += [plt.Line2D([], [], color="0.2", lw=0.9, label="signal (night mean)"),
                plt.Line2D([], [], color="#1f77b4", lw=1.3, ls="--", label="molecular x CL")]
    fig.legend(handles=handles, loc="upper center", ncol=5, fontsize=8.5, frameon=False,
               bbox_to_anchor=(0.5, 1.02))
    fig.suptitle(f"{label} — profiles + selected windows (top) and calibration C_L ± err per method (bottom)",
                 y=1.10, fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.88])
    p = OUT / f"profiles_{label}.png"
    fig.savefig(p, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"    saved {p.name}")


def plot_pcolor(label, detail, per_night):
    nights = list(detail.keys())
    if not nights:
        return
    n = len(nights)
    fig, axes = plt.subplots(2, n, figsize=(5.0 * n, 8.6), squeeze=False,
                             gridspec_kw={"height_ratios": [3, 1]})
    for col, ds in enumerate(nights):
        ax = axes[0][col]
        d = detail[ds]
        rng_km = d["range_alc"] / 1e3
        stack, pmol = d["signal_stack"], d["p_mol"]
        hours = d.get("hours")
        if hours is None or len(hours) != stack.shape[0]:
            hours = np.arange(stack.shape[0], dtype=float)
        with np.errstate(divide="ignore", invalid="ignore"):
            ratio = stack / pmol[None, :]
        ratio_pos = np.where(ratio > 0, ratio, np.nan)
        fin = ratio_pos[np.isfinite(ratio_pos)]
        vlo, vhi = (np.nanpercentile(fin, [5, 97]) if fin.size else (1e-3, 1.0))
        msh = ax.pcolormesh(hours, rng_km, ratio_pos.T, shading="auto", cmap="turbo",
                            norm=matplotlib.colors.LogNorm(vmin=max(vlo, vhi / 1e3),
                                                           vmax=max(vhi, vlo * 10)))
        ax.set_ylim(0, 7)
        res = per_night[ds]
        opt = res.get("optimal")
        if opt is not None and opt.cell_flag is not None and np.any(opt.cell_flag):
            ax.contourf(hours, rng_km, opt.cell_flag.T.astype(float), levels=[0.5, 1.5],
                        colors="none", hatches=["xxx"], zorder=4)
            ax.contour(hours, rng_km, opt.cell_flag.T.astype(float), levels=[0.5],
                       colors="magenta", linewidths=0.7, zorder=4)
        for m in METHODS:
            w = res[m]
            if w.ok:
                ax.axhline(w.center_m / 1e3, color=METHOD_COLORS[m], lw=1.5, ls="--")
        ax.set_title(ds, fontsize=10)
        ax.set_xlabel("hours since start")
        plt.colorbar(msh, ax=ax, label="signal / molecular", pad=0.02)
        _cl_panel(axes[1][col], res)
    axes[0][0].set_ylabel("altitude AGL (km)")
    handles = [plt.Line2D([], [], color=METHOD_COLORS[m], lw=1.5, ls="--", label=METHOD_LABEL[m]) for m in METHODS]
    handles += [plt.Line2D([], [], color="magenta", lw=1.0, label="optimal: flagged aerosol/cloud (excluded)")]
    fig.legend(handles=handles, loc="upper center", ncol=4, fontsize=8.3, frameon=False,
               bbox_to_anchor=(0.5, 1.02))
    fig.suptitle(f"{label} — signal/molecular ratio + window centres (top) and calibration C_L ± err (bottom)",
                 y=1.10, fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.88])
    p = OUT / f"pcolor_{label}.png"
    fig.savefig(p, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"    saved {p.name}")


# ----------------------------------------------------------------------------
# STATS aggregation
# ----------------------------------------------------------------------------
def aggregate(all_results):
    """all_results: {inst_label: {ds: {method: mw}}} -> per (inst, method) stats."""
    rows = []
    for inst in INSTRUMENTS:
        label = inst["label"]
        pn = all_results.get(label, {})
        n_nights = len(pn)
        for m in METHODS:
            cls, r2s, tcvs, rels = [], [], [], []
            n_cal = 0
            for ds, res in pn.items():
                w = res[m]
                if calibrates(w):
                    n_cal += 1
                    cls.append(w.cl)
                    r2s.append(w.r2)
                    if np.isfinite(w.temporal_cv):
                        tcvs.append(w.temporal_cv)
                    rels.append(w.rel_error)
            cls = np.array(cls, float)
            cl_cv = float(np.std(cls) / np.abs(np.mean(cls)) * 100) if cls.size >= 3 else np.nan
            rows.append(dict(inst=label, itype=inst["itype"].value, method=m,
                             n_nights=n_nights, n_cal=n_cal,
                             cl_cv=cl_cv,
                             med_r2=float(np.median(r2s)) if r2s else np.nan,
                             med_tcv=float(np.median(tcvs)) if tcvs else np.nan,
                             med_rel=float(np.median(rels)) if rels else np.nan))
    return rows


def select_best(rows):
    """Rank methods balancing usability (calibrated fraction), night-to-night STABILITY
    (low CL_CV — the dominant quality for a calibration constant) and cleanliness (low
    temporal_cv = aerosol-free). A method that calibrates many nights but with a huge CV
    (e.g. by selecting noisy/aerosol windows) is penalised, not rewarded."""
    score = {}
    for m in METHODS:
        mr = [r for r in rows if r["method"] == m]
        fracs = [r["n_cal"] / r["n_nights"] for r in mr if r["n_nights"] > 0]
        # only count CV/tcv where the method actually produced enough nights to be meaningful
        cvs = [r["cl_cv"] for r in mr if np.isfinite(r["cl_cv"])]
        tcvs = [r["med_tcv"] for r in mr if np.isfinite(r["med_tcv"])]
        score[m] = dict(frac=float(np.mean(fracs)) if fracs else 0.0,
                        cv=float(np.mean(cvs)) if cvs else np.nan,
                        tcv=float(np.mean(tcvs)) if tcvs else np.nan)
    cv_ref = 20.0   # %, a "good" night-to-night CV; CV at/below this scores ~1
    for m, v in score.items():
        cv = v["cv"] if np.isfinite(v["cv"]) else 999.0
        tcv = v["tcv"] if np.isfinite(v["tcv"]) else 1.0
        stab = cv_ref / max(cv, cv_ref)           # 1 if CV<=20%, ->0 as CV grows
        clean = 1.0 / (1.0 + tcv)                  # 1 if temporal_cv=0, 0.5 at tcv=1
        # geometric-ish balance: need usable nights AND stability AND cleanliness
        v["rank_score"] = v["frac"] * stab * clean
    ranking = sorted(score.items(), key=lambda kv: kv[1]["rank_score"], reverse=True)
    return score, ranking


def plot_summary(rows):
    insts = [i["label"] for i in INSTRUMENTS]
    fig, axes = plt.subplots(1, 2, figsize=(16, 6))
    x = np.arange(len(insts))
    w = 0.13
    for k, m in enumerate(METHODS):
        fr = [next((r["n_cal"] / r["n_nights"] if r["n_nights"] else 0
                    for r in rows if r["inst"] == lab and r["method"] == m), 0) for lab in insts]
        cv = [next((r["cl_cv"] for r in rows if r["inst"] == lab and r["method"] == m), np.nan) for lab in insts]
        off = (k - (len(METHODS) - 1) / 2.0) * w
        axes[0].bar(x + off, fr, w, color=METHOD_COLORS[m], label=METHOD_LABEL[m])
        axes[1].bar(x + off, cv, w, color=METHOD_COLORS[m])
    axes[0].set_ylabel("calibrated-night fraction")
    axes[0].set_title("Usable nights per method")
    axes[1].set_ylabel("lidar-constant CV (%)")
    axes[1].set_title("Night-to-night stability (lower = better)")
    for a in axes:
        a.set_xticks(x)
        a.set_xticklabels(insts, rotation=30, ha="right", fontsize=8)
        a.grid(True, axis="y", alpha=0.25)
    axes[0].legend(ncol=3, fontsize=8, loc="upper right")
    fig.suptitle("Molecular-window methods across sites/instruments", fontsize=13)
    fig.tight_layout()
    p = OUT / "summary_methods_multisite.png"
    fig.savefig(p, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  saved {p.name}")


def plot_timeseries(all_results):
    """Per-instrument time series of the calibration constant (CL proxy) for each method,
    using the shared method colour code."""
    fig, axes = plt.subplots(2, 3, figsize=(18, 9), squeeze=False)
    axes = axes.ravel()
    for ax, inst in zip(axes, INSTRUMENTS):
        label = inst["label"]
        pn = all_results.get(label, {})
        any_pt = False
        for m in METHODS:
            dates, cls = [], []
            for ds in sorted(pn):
                w = pn[ds][m]
                if calibrates(w) and np.isfinite(w.cl) and w.cl > 0:
                    dates.append(datetime.strptime(ds, "%Y%m%d"))
                    cls.append(w.cl)
            if dates:
                any_pt = True
                ax.plot(dates, cls, "-o", ms=4, lw=1.0, color=METHOD_COLORS[m], label=METHOD_LABEL[m])
        # Robust y-range from the stable (gated) methods so they stay readable; the wild
        # main/calipso excursions then clip at the top rather than compressing everything.
        stable = []
        for sm in ("improved", "optimal", "earlinet", "matlab", "bellini"):
            for ds in pn:
                w = pn[ds][sm]
                if calibrates(w) and np.isfinite(w.cl) and w.cl > 0:
                    stable.append(w.cl)
        if len(stable) >= 3:
            lo, hi = np.nanpercentile(stable, [5, 95])
            pad = 0.6 * (hi - lo) + 1e-30
            ax.set_ylim(max(0.0, lo - pad), hi + pad)
        ax.set_title(label, fontsize=10)
        ax.grid(True, alpha=0.25)
        ax.set_ylabel("lidar constant (proxy, a.u.)")
        ax.xaxis.set_major_formatter(mdates.DateFormatter("%m-%d"))
        for lab in ax.get_xticklabels():
            lab.set_rotation(30)
            lab.set_fontsize(7)
            lab.set_ha("right")
        if not any_pt:
            ax.text(0.5, 0.5, "no calibrations", ha="center", va="center",
                    transform=ax.transAxes, color="0.5")
    handles = [plt.Line2D([], [], color=METHOD_COLORS[m], marker="o", lw=1.0, label=METHOD_LABEL[m])
               for m in METHODS]
    fig.legend(handles=handles, loc="upper center", ncol=7, fontsize=8.5, frameon=False,
               bbox_to_anchor=(0.5, 1.02))
    fig.suptitle("Calibration-constant time series per method (per instrument)", y=1.06, fontsize=13)
    fig.tight_layout(rect=[0, 0, 1, 0.93])
    p = OUT / "timeseries_methods_multisite.png"
    fig.savefig(p, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  saved {p.name}")


def write_tables(rows, score, ranking):
    lines = ["# Molecular-window methods — multi-site comparison (Payerne, Amsterdam, EDT)\n\n",
             f"Sampled {len(STATS_NIGHTS)} nights (every 2nd day, Mar+Apr 2026) per instrument. "
             "`n_cal` = nights that calibrate through the full pipeline (rel_error<=15%). "
             "`CL_CV` = night-to-night lidar-constant scatter (lower = more stable).\n\n",
             "inst | type | method | nights | n_cal | CL_CV% | med_R2 | med_tcv | med_rel%\n",
             "---|---|---|---|---|---|---|---|---\n"]
    for r in rows:
        lines.append(f"{r['inst']} | {r['itype']} | {r['method']} | {r['n_nights']} | {r['n_cal']} | "
                     f"{r['cl_cv']:.1f} | {r['med_r2']:.3f} | {r['med_tcv']:.2f} | {r['med_rel']:.1f}\n"
                     .replace("nan", "-"))
    lines.append("\n## Ranking (mean over instruments)\n\n")
    lines.append("method | calibrated-fraction | mean CL_CV% | mean temporal_cv | score\n---|---|---|---|---\n")
    for m, v in ranking:
        lines.append(f"{m} | {v['frac']:.2f} | {v['cv']:.1f} | {v['tcv']:.2f} | {v['rank_score']:.3f}\n"
                     .replace("nan", "-"))
    best = ranking[0][0]
    lines.append(f"\n**Best overall: `{best}`** "
                 f"(highest usable-night fraction at competitive night-to-night stability).\n")
    p = OUT / "method_comparison_multisite.md"
    p.write_text("".join(lines), encoding="utf-8")
    print(f"  saved {p.name}  -> best = {best}")
    return best


def main():
    warnings.filterwarnings("ignore")
    logging.getLogger().setLevel(logging.ERROR)
    all_results = {}
    for inst in INSTRUMENTS:
        print(f"== {inst['label']} ==")
        per_night, detail = run_instrument(inst)
        all_results[inst["label"]] = per_night
        if detail:
            plot_profiles(inst["label"], detail, per_night)
            plot_pcolor(inst["label"], detail, per_night)
    rows = aggregate(all_results)
    score, ranking = select_best(rows)
    plot_summary(rows)
    plot_timeseries(all_results)
    write_tables(rows, score, ranking)
    print("COMPARE_DONE")


if __name__ == "__main__":
    main()
