"""
analyze_rayleigh_diag.py — full-network Rayleigh diagnosis. Combines:
  * net_L1_<label>.json  (network run): per-night v2 constant -> time series, valid%, sigma_SD, outlier%
  * diaglite_<label>.json (light diag): sampled instrument-health + FT-aerosol metrics
  * scope_network_2026.json: n_days -> clear-sky fraction

Maps the four candidate failure causes:
  not enough clear sky  -> clear_frac = clear-reaching nights / archive days     (low = bad)
  lots of FT aerosol    -> scat_med   = median window scattering ratio            (high = bad)
  low laser             -> sig_med    = near-range signal strength (+ laser age)  (low = bad)
  electronic background -> bg_noise / snr = far-range noise                       (high bg / low snr = bad)

A station is problematic if (within its type) it is in the worst quartile of valid yield OR worst
decile of sigma_SD / outlier rate; each is assigned the dominant cause. Produces paginated network
time-series grids, a cause-overview figure, a problematic-station table, and a JSON summary.
"""
from __future__ import annotations
import glob
import json
import sys
from collections import Counter
from datetime import datetime
from pathlib import Path

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

REPO = Path(__file__).resolve().parents[1]
DIR = Path("C:/DATA/Projects/202606_E-PROFILE_calibration/figs_paper_validation/rayleigh_diag")
NET = Path("C:/DATA/Projects/202606_E-PROFILE_calibration/figs_paper_validation/network_v2_v11")
META = {m["label"]: m for m in json.loads((REPO / "validation" / "scope_network_2026.json").read_text())}
QC = 15.0
CAUSE_COLOR = {"clear-sky": "#9467bd", "FT-aerosol": "#d62728", "low-laser": "#ff7f0e",
               "background": "#8c564b", "other": "#7f7f7f", "ok": "#2ca02c"}


def rolling_median(x, win=9):
    n = len(x); h = win // 2
    return np.array([np.median(x[max(0, i - h):min(n, i + h + 1)]) for i in range(n)])


def sigma_sd(c):
    c = np.asarray(c, float)
    if c.size < 4:
        return np.nan
    m = np.median(c)
    return float(1.4826 * np.median(np.abs(np.diff(c))) / np.sqrt(2) / abs(m) * 100) if m > 0 else np.nan


def outlier_pct(c):
    c = np.asarray(c, float)
    if c.size < 5:
        return np.nan
    r = c - rolling_median(c); r = r - np.median(r)
    s = 1.4826 * np.median(np.abs(r))
    return float(100 * np.mean(np.abs(r) > 3 * s)) if s > 0 else 0.0


def _net_metrics(level, label, n_days):
    """v2 metrics from a net_<level>_<label>.json: (cl_dts, cls, n_fit, clear%, valid-of-clear%,
    valid-of-archive%, sigma_SD, outlier%). Returns None if the file is missing."""
    netf = NET / f"net_{level}_{label}.json"
    if not netf.exists():
        return None
    nights = json.loads(netf.read_text())
    order = sorted(nights)
    cl_dts, cls = [], []
    for ds in order:
        ok, cl, rel = nights[ds]["v2"]
        if ok and np.isfinite(cl) and cl > 0 and np.isfinite(rel) and rel <= QC:
            cl_dts.append(datetime.strptime(ds, "%Y%m%d")); cls.append(cl)
    cls = np.asarray(cls, float)
    n_fit = len(order)
    return dict(cl_dts=cl_dts, cls=cls, n_fit=n_fit, n_valid=int(cls.size),
                clear_frac=100.0 * n_fit / n_days,
                valid_of_clear=100.0 * cls.size / n_fit if n_fit else np.nan,
                valid_of_data=100.0 * cls.size / n_days,
                sigma_sd=sigma_sd(cls), outlier=outlier_pct(cls))


def load_station(label):
    m = META[label]
    n_days = max(int(m.get("n_days", 1)), 1)
    L1 = _net_metrics("L1", label, n_days)
    if L1 is None:
        return None
    L2 = _net_metrics("L2", label, n_days)
    hk = json.loads((DIR / f"diaglite_{label}.json").read_text()) if (DIR / f"diaglite_{label}.json").exists() else {}
    r = dict(
        label=label, group=m["group"], site=m.get("site", label), n_days=n_days,
        n_fit=L1["n_fit"], n_valid=L1["n_valid"], cl_dts=L1["cl_dts"], cls=L1["cls"],
        clear_frac=L1["clear_frac"], valid_of_clear=L1["valid_of_clear"],
        valid_of_data=L1["valid_of_data"], sigma_sd=L1["sigma_sd"], outlier=L1["outlier"],
        scat_med=float(hk.get("scat", np.nan)), sig_med=float(hk.get("sig", np.nan)),
        bg_noise=float(hk.get("bg", np.nan)), laser_med=float(hk.get("laser", np.nan)),
        wt_med=float(hk.get("wtrans", np.nan)))
    # L2 comparison metrics (NaN if no L2 stream)
    for k in ("clear_frac", "valid_of_clear", "valid_of_data", "sigma_sd", "outlier"):
        r[f"L2_{k}"] = L2[k] if L2 else np.nan
    return r


def classify(rows):
    for grp in set(r["group"] for r in rows):
        sub = [r for r in rows if r["group"] == grp]

        def pct(key, q):
            v = [r[key] for r in sub if np.isfinite(r[key])]
            return float(np.percentile(v, q)) if v else np.nan
        thr = dict(clear=pct("clear_frac", 25), valid=pct("valid_of_data", 25),
                   scat=pct("scat_med", 75), sig=pct("sig_med", 25), bg=pct("bg_noise", 75),
                   sd=pct("sigma_sd", 90), out=pct("outlier", 90), laser=pct("laser_med", 85))
        snr = [r["sig_med"] / r["bg_noise"] for r in sub
               if np.isfinite(r["sig_med"]) and np.isfinite(r["bg_noise"]) and r["bg_noise"] > 0]
        snr_lo = float(np.percentile(snr, 25)) if snr else np.nan
        for r in sub:
            r["snr"] = (r["sig_med"] / r["bg_noise"]) if (np.isfinite(r["bg_noise"]) and r["bg_noise"] > 0) else np.nan
            r["problematic"] = bool(
                (np.isfinite(r["valid_of_data"]) and r["valid_of_data"] <= thr["valid"])
                or (np.isfinite(r["sigma_sd"]) and r["sigma_sd"] >= thr["sd"])
                or (np.isfinite(r["outlier"]) and r["outlier"] >= thr["out"]))
            sev = {}
            if np.isfinite(r["clear_frac"]) and np.isfinite(thr["clear"]) and r["clear_frac"] <= thr["clear"]:
                sev["clear-sky"] = (thr["clear"] - r["clear_frac"]) / (abs(thr["clear"]) + 1e-9)
            if np.isfinite(r["scat_med"]) and np.isfinite(thr["scat"]) and r["scat_med"] >= thr["scat"]:
                sev["FT-aerosol"] = (r["scat_med"] - thr["scat"]) / (abs(thr["scat"]) + 1e-9)
            if (np.isfinite(r["sig_med"]) and np.isfinite(thr["sig"]) and r["sig_med"] <= thr["sig"]) or \
               (np.isfinite(r["laser_med"]) and np.isfinite(thr["laser"]) and r["laser_med"] >= thr["laser"]):
                sev["low-laser"] = (thr["sig"] - r["sig_med"]) / (abs(thr["sig"]) + 1e-9) if np.isfinite(r["sig_med"]) else 0.5
            if (np.isfinite(r["bg_noise"]) and np.isfinite(thr["bg"]) and r["bg_noise"] >= thr["bg"]) or \
               (np.isfinite(r["snr"]) and np.isfinite(snr_lo) and r["snr"] <= snr_lo):
                sev["background"] = (r["bg_noise"] - thr["bg"]) / (abs(thr["bg"]) + 1e-9) if np.isfinite(r["bg_noise"]) else 0.5
            r["causes"] = sorted(sev, key=lambda k: -sev[k])
            r["cause"] = (r["causes"][0] if r["causes"] else "other") if r["problematic"] else "ok"
    return rows


def ts_panel(ax, dts, cls):
    cls = np.asarray(cls, float)
    if cls.size:
        med = np.median(cls)
        r = cls - rolling_median(cls); r = r - np.median(r)
        s = 1.4826 * np.median(np.abs(r))
        out = np.abs(r) > 3 * s if s > 0 else np.zeros(cls.size, bool)
        ax.axhline(med, color="green", lw=0.8)
        ax.plot(dts, cls, "o", ms=2.2, color="#1f77b4", alpha=0.7)
        if out.any():
            ax.plot(np.array(dts)[out], cls[out], "o", ms=4, color="red")
    ax.tick_params(labelsize=5)
    for lab in ax.get_xticklabels():
        lab.set_rotation(30); lab.set_ha("right")


def ts_pages(rows, group, prefix, ncol=5, nrow=5):
    sub = [r for r in rows if r["group"] == group and r["cls"].size >= 1]
    sub.sort(key=lambda r: (np.inf if not np.isfinite(r["valid_of_data"]) else r["valid_of_data"]))
    per = ncol * nrow
    npages = max(int(np.ceil(len(sub) / per)), 1)
    for p in range(npages):
        chunk = sub[p * per:(p + 1) * per]
        fig, axes = plt.subplots(nrow, ncol, figsize=(4.0 * ncol, 2.4 * nrow), squeeze=False)
        for ax in axes.ravel():
            ax.axis("off")
        for k, r in enumerate(chunk):
            ax = axes[k // ncol][k % ncol]; ax.axis("on")
            ax.set_facecolor("white" if r["cause"] == "ok" else "#fff5f5")
            ts_panel(ax, r["cl_dts"], r["cls"])
            ax.set_title(f"{r['site'][:16]}\nvalid {r['valid_of_data']:.0f}% σ{r['sigma_sd']:.0f} out{r['outlier']:.0f}%"
                         + ("" if r["cause"] == "ok" else f" [{r['cause']}]"), fontsize=7)
        fig.suptitle(f"{group} nightly lidar constant (v2, L1) — page {p+1}/{npages} "
                     f"(sorted by valid%; pink = problematic)", fontsize=12)
        fig.tight_layout(rect=[0, 0, 1, 0.97])
        fig.savefig(DIR / f"{prefix}_p{p+1}.png", dpi=110); plt.close(fig)
    return npages


def overview(rows, path):
    from matplotlib.patches import Patch
    fig, axes = plt.subplots(2, 2, figsize=(20, 12))
    chm = [r for r in rows if r["group"] == "CHM15k"]
    ax = axes[0][0]
    s = sorted(chm, key=lambda r: (np.inf if not np.isfinite(r["valid_of_data"]) else r["valid_of_data"]))
    ax.bar(range(len(s)), [r["valid_of_data"] for r in s], color=[CAUSE_COLOR[r["cause"]] for r in s])
    ax.set_title("CHM15k valid-calibration yield (sorted), coloured by diagnosed cause", fontsize=11)
    ax.set_xlabel("stream"); ax.set_ylabel("valid calibrations / archive days (%)")
    ax.legend(handles=[Patch(color=c, label=k) for k, c in CAUSE_COLOR.items()], fontsize=8, ncol=2)
    ax.grid(axis="y", alpha=0.3)
    for ax, xk, yk, xl, yl, ttl, logx in [
        (axes[0][1], "clear_frac", "valid_of_data", "clear-sky fraction (%)", "valid calibration %", "Clear-sky availability vs yield", False),
        (axes[1][0], "scat_med", "sigma_sd", "median window scattering ratio (FT aerosol)", "σ_SD (%)", "FT aerosol vs short-term variability", False),
        (axes[1][1], "snr", "outlier", "SNR proxy (near signal / far noise)", "outlier %", "Signal-to-background vs outlier rate", True)]:
        for r in chm:
            ax.scatter(r[xk], r[yk], s=25, color=CAUSE_COLOR[r["cause"]], alpha=0.7, edgecolor="k", linewidth=0.2)
        if logx:
            ax.set_xscale("log")
        ax.set_xlabel(xl); ax.set_ylabel(yl); ax.set_title(ttl, fontsize=11); ax.grid(alpha=0.3)
    fig.suptitle("CHM15k network Rayleigh diagnosis — cause attribution", fontsize=14)
    fig.tight_layout(rect=[0, 0, 1, 0.97]); fig.savefig(path, dpi=130); plt.close(fig)


def fig_l1_vs_l2(rows, path):
    """Paired per-stream L1 vs L2 (v2): clear-sky reach, yield-on-clear, sigma_SD, outlier%."""
    from matplotlib.patches import Patch
    panels = [("clear_frac", "clear-sky reach (% of archive days)", (0, 100), "above=L2 reaches more"),
              ("valid_of_clear", "valid / clear nights (%)", (0, 100), "above=L2 higher yield"),
              ("sigma_sd", "σ_SD (%)", (0, 25), "below=L2 more precise"),
              ("outlier", "outlier rate (%)", (0, 35), "below=L2 fewer outliers")]
    fig, axes = plt.subplots(1, 4, figsize=(24, 6))
    for ax, (key, lab, lim, hint) in zip(axes, panels):
        for r in rows:
            x, y = r.get(key), r.get(f"L2_{key}")
            if x is not None and y is not None and np.isfinite(x) and np.isfinite(y):
                ax.scatter(x, y, s=26, color=CAUSE_COLOR.get(r["cause"], "#7f7f7f"),
                           alpha=0.7, edgecolor="k", linewidth=0.2)
        ax.plot(lim, lim, "k--", lw=1, alpha=0.6)
        ax.set_xlim(lim); ax.set_ylim(lim)
        ax.set_xlabel(f"L1 — {lab}"); ax.set_ylabel(f"L2 — {lab}")
        ax.set_title(f"{lab}\n({hint})", fontsize=10); ax.grid(alpha=0.3)
    axes[0].legend(handles=[Patch(color=c, label=k) for k, c in CAUSE_COLOR.items()], fontsize=7, ncol=2)
    fig.suptitle("Rayleigh calibration L1 vs L2, same streams (v2; colour = L1-diagnosed cause)", fontsize=14)
    fig.tight_layout(rect=[0, 0, 1, 0.96]); fig.savefig(path, dpi=130); plt.close(fig)


def l1l2_aggregate(rows):
    out = {}
    for t in ["CHM15k", "Mini-MPL"]:
        sub = [r for r in rows if r["group"] == t]
        def med(key):
            v = [r[key] for r in sub if isinstance(r.get(key), (int, float)) and np.isfinite(r.get(key))]
            return float(np.median(v)) if v else float("nan")
        out[t] = {k: med(k) for k in
                  ("clear_frac", "valid_of_clear", "sigma_sd", "outlier",
                   "L2_clear_frac", "L2_valid_of_clear", "L2_sigma_sd", "L2_outlier")}
        out[t]["n"] = len(sub)
    return out


def main():
    labels = [Path(f).stem[len("net_L1_"):] for f in glob.glob(str(NET / "net_L1_*.json"))]
    rows = [r for r in (load_station(l) for l in labels
                        if META.get(l, {}).get("group") in ("CHM15k", "Mini-MPL")) if r]
    classify(rows)
    keep = ["label", "group", "site", "n_days", "n_fit", "n_valid", "clear_frac", "valid_of_clear",
            "valid_of_data", "scat_med", "sig_med", "bg_noise", "snr", "laser_med", "wt_med",
            "sigma_sd", "outlier", "L2_clear_frac", "L2_valid_of_clear", "L2_sigma_sd", "L2_outlier",
            "problematic", "cause", "causes"]
    (DIR / "rayleigh_diag_summary.json").write_text(
        json.dumps([{k: r.get(k) for k in keep} for r in rows], indent=1, default=float), encoding="utf-8")
    npages = ts_pages(rows, "CHM15k", "fig_ts_CHM15k")
    ts_pages(rows, "Mini-MPL", "fig_ts_MiniMPL", ncol=5, nrow=1)
    overview(rows, DIR / "fig_cause_overview.png")
    fig_l1_vs_l2(rows, DIR / "fig_l1_vs_l2.png")
    agg = l1l2_aggregate(rows)
    (DIR / "l1l2_aggregate.json").write_text(json.dumps(agg, indent=1), encoding="utf-8")

    prob = sorted([r for r in rows if r["problematic"]], key=lambda r: (r["group"], r["valid_of_data"]))
    L = ["# Problematic Rayleigh stations — diagnosed cause\n\n",
         f"{len(prob)} of {len(rows)} CHM15k+Mini-MPL streams flagged problematic "
         "(worst-quartile valid yield, or worst-decile σ_SD / outlier rate, within type).\n\n",
         "site | stream | type | valid% | clear% | σ_SD% | out% | scat | SNR | laser(h) | **cause** | also\n",
         "---|---|---|---|---|---|---|---|---|---|---|---\n"]
    for r in prob:
        L.append(f"{r['site'][:18]} | {r['label']} | {r['group']} | {r['valid_of_data']:.0f} | "
                 f"{r['clear_frac']:.0f} | {r['sigma_sd']:.1f} | {r['outlier']:.0f} | {r['scat_med']:.2f} | "
                 f"{r['snr']:.0f} | {r['laser_med']:.0f} | **{r['cause']}** | {', '.join(r['causes'][1:]) or '-'}\n"
                 .replace("nan", "-"))
    (DIR / "problematic_stations.md").write_text("".join(L), encoding="utf-8")
    print(f"{len(rows)} streams; {len(prob)} problematic; CHM15k TS pages={npages}")
    print("cause breakdown:", dict(Counter(r["cause"] for r in rows)))
    print("problematic by cause:", dict(Counter(r["cause"] for r in prob)))
    print("\n=== L1 vs L2 (median over streams) ===")
    for t, a in agg.items():
        print(f"{t:9s} n={a['n']:3d} | clear% L1 {a['clear_frac']:.0f}->L2 {a['L2_clear_frac']:.0f} | "
              f"valid/clear% L1 {a['valid_of_clear']:.0f}->L2 {a['L2_valid_of_clear']:.0f} | "
              f"σ_SD L1 {a['sigma_sd']:.1f}->L2 {a['L2_sigma_sd']:.1f} | "
              f"out% L1 {a['outlier']:.0f}->L2 {a['L2_outlier']:.0f}".replace("nan", "-"))


if __name__ == "__main__":
    main()
