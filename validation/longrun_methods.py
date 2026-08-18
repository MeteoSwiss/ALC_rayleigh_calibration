"""
longrun_methods.py — full-archive comparison of the 7 molecular-window methods across
10 CHM15k + 4 Mini-MPL sites (E-PROFILE L2 monthly), sampling 5 nights/month over each
site's entire available period. Saves per-instrument JSON incrementally (so partial results
survive), then aggregates into per-method stats, a summary bar chart, per-instrument
calibration-constant time series, and a markdown table.
"""
from __future__ import annotations
import os, sys, glob, json, logging, warnings
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

from calibration import calibrate_rayleigh, CalibrationOptions, InstrumentInfo, DataLevel
from calibration.config import InstrumentType
import compare_molecular_methods as cm
from compare_molecular_methods import METHODS, METHOD_COLORS, METHOD_LABEL, run_methods, calibrates

OUT = Path("C:/DATA/Projects/202606_E-PROFILE_calibration/figs_paper_validation/molecular_methods_longrun")
OUT.mkdir(parents=True, exist_ok=True)
cm.OUT = OUT                          # redirect cm's figure/table writers here
ROOT = Path("A:/E-PROFILE_L2_monthly")
SAMPLE_DAYS = [3, 9, 15, 21, 27]      # 5 nights/month

# (wmo, identifier, itype) — 10 CHM15k (6 named + 4 best-covered) + 4 Mini-MPL
SEL = [
    ("0-20000-0-06610", "A", "CHM15k"), ("0-20000-0-10393", "0", "CHM15k"),
    ("0-380-5-1", "0", "CHM15k"),       ("0-250-1001-07151", "B", "CHM15k"),
    ("0-20008-0-UGR", "A", "CHM15k"),   ("0-20008-0-INO", "A", "CHM15k"),
    ("0-20000-0-01311", "A", "CHM15k"), ("0-20000-0-01492", "A", "CHM15k"),
    ("0-20000-0-10140", "0", "CHM15k"), ("0-20000-0-10962", "0", "CHM15k"),
    ("0-20000-0-07110", "A", "Mini-MPL"), ("0-20000-0-07617", "A", "Mini-MPL"),
    ("0-20000-0-07774", "A", "Mini-MPL"), ("0-20000-0-07145", "A", "Mini-MPL"),
]
NAMES = {
    "0-20000-0-06610": "Payerne", "0-20000-0-10393": "Lindenberg", "0-380-5-1": "Aosta",
    "0-250-1001-07151": "Palaiseau", "0-20008-0-UGR": "Granada", "0-20008-0-INO": "Magurele",
    "0-20000-0-01311": "Bergen", "0-20000-0-01492": "Oslo", "0-20000-0-10140": "Hamburg",
    "0-20000-0-10962": "Hohenpeiss", "0-20000-0-07110": "Brest-MPL", "0-20000-0-07617": "Toulouse-MPL",
    "0-20000-0-07774": "Corsica-MPL", "0-20000-0-07145": "SIRTA-MPL",
}
ITYPE = {"CHM15k": InstrumentType.CHM15k, "Mini-MPL": InstrumentType.MINI_MPL}

_man = {m["wmo"] + "_" + m.get("identifier", "A"): m for m in json.load(open("C:/DATA/Projects/202606_E-PROFILE_calibration/stations_l2_manifest.json"))}


def build_instruments():
    insts = []
    for wmo, ident, itype in SEL:
        m = _man.get(f"{wmo}_{ident}")
        if m is None:
            print(f"  WARN: {wmo}_{ident} not in manifest; skipping")
            continue
        label = f"{NAMES.get(wmo, wmo)}_{itype.replace('Mini-MPL', 'MPL')}"
        insts.append(dict(label=label, wmo=wmo, ident=ident, itype=ITYPE[itype],
                          lat=m["lat"], lon=m["lon"], alt=m["alt"], n_months=m["n_months"]))
    return insts


def avail_months(wmo, ident):
    months = []
    pre = f"L2_{wmo}_{ident}"
    for f in glob.glob(str(ROOT / wmo / "*" / f"{pre}*.nc")):
        core = os.path.basename(f)[len(pre):-3]
        if len(core) >= 6 and core[:6].isdigit():
            months.append(core[:6])
    return sorted(set(months))


def base_options():
    o = CalibrationOptions.from_json(Path("options.json"))
    o.folder_root = ROOT
    o.data_level = DataLevel.L2_MONTHLY
    o.molecular_source = "standard"
    o.plot_main = False
    o.plot_all = False
    o.folder_output = OUT
    return o


def mw_to_list(w):
    return [bool(w.ok), float(w.cl), float(w.cl_err), float(w.rel_error), float(w.r2),
            float(w.temporal_cv), float(w.scattering_ratio), float(w.start_m), float(w.end_m),
            float(w.center_m)]


def run_instrument(inst):
    o = base_options()
    info = InstrumentInfo(site_name=inst["label"], wmo_id=inst["wmo"], identifier=inst["ident"],
                          instrument_type=inst["itype"], latitude=inst["lat"],
                          longitude=inst["lon"], altitude=inst["alt"])
    months = avail_months(inst["wmo"], inst["ident"])
    per_night = {}
    n_cal = 0
    for ym in months:
        for d in SAMPLE_DAYS:
            ds = f"{ym}{d:02d}"
            fin = {}
            try:
                calibrate_rayleigh(ds, info, o, fit_inputs_out=fin)
            except Exception:
                continue
            if not fin:
                continue
            res = run_methods(fin["signal"], fin["p_mol"], fin["range_alc"], fin["signal_stack"])
            per_night[ds] = res
            if any(calibrates(w) for w in res.values()):
                n_cal += 1
    # incremental JSON backup
    dump = {ds: {m: mw_to_list(res[m]) for m in METHODS} for ds, res in per_night.items()}
    (OUT / f"results_{inst['label']}.json").write_text(json.dumps(dump), encoding="utf-8")
    print(f"  {inst['label']}: {len(months)} months, {len(per_night)} fit-nights, "
          f"{n_cal} with >=1 calibration  (saved JSON)")
    return per_night


def plot_timeseries_grid(all_results, insts):
    n = len(insts)
    ncol = 4
    nrow = int(np.ceil(n / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(5.0 * ncol, 3.2 * nrow), squeeze=False)
    axA = axes.ravel()
    for ax, inst in zip(axA, insts):
        pn = all_results.get(inst["label"], {})
        for m in METHODS:
            dates, cls = [], []
            for ds in sorted(pn):
                w = pn[ds][m]
                if calibrates(w) and np.isfinite(w.cl) and w.cl > 0:
                    dates.append(datetime.strptime(ds, "%Y%m%d"))
                    cls.append(w.cl)
            if dates:
                ax.plot(dates, cls, "-o", ms=2.5, lw=0.7, color=METHOD_COLORS[m])
        # robust y-range from the stable methods so main spikes clip
        stable = [pn[ds][sm].cl for sm in ("eprof_v1.2", "eprof_v2", "earlinet", "eprof_v0.25", "bellini")
                  for ds in pn if calibrates(pn[ds][sm]) and np.isfinite(pn[ds][sm].cl) and pn[ds][sm].cl > 0]
        if len(stable) >= 3:
            lo, hi = np.nanpercentile(stable, [3, 97])
            pad = 0.6 * (hi - lo) + 1e-30
            ax.set_ylim(max(0.0, lo - pad), hi + pad)
        ax.set_title(inst["label"], fontsize=9)
        ax.grid(True, alpha=0.25)
        ax.xaxis.set_major_formatter(mdates.DateFormatter("%Y"))
        for lab in ax.get_xticklabels():
            lab.set_fontsize(6)
        ax.tick_params(axis="y", labelsize=6)
    for ax in axA[n:]:
        ax.axis("off")
    handles = [plt.Line2D([], [], color=METHOD_COLORS[m], marker="o", lw=0.8, label=METHOD_LABEL[m]) for m in METHODS]
    fig.legend(handles=handles, loc="upper center", ncol=7, fontsize=8, frameon=False, bbox_to_anchor=(0.5, 1.01))
    fig.suptitle("Calibration-constant time series per method — full archive (per site)", y=1.03, fontsize=13)
    fig.tight_layout(rect=[0, 0, 1, 0.97])
    fig.savefig(OUT / "timeseries_longrun.png", dpi=140, bbox_inches="tight")
    plt.close(fig)
    print("  saved timeseries_longrun.png")


def main():
    warnings.filterwarnings("ignore")
    logging.getLogger().setLevel(logging.ERROR)
    insts = build_instruments()
    cm.INSTRUMENTS = insts            # so cm.aggregate/select_best/plot_summary/write_tables use these
    all_results = {}
    for inst in insts:
        print(f"== {inst['label']}  ({inst['wmo']}_{inst['ident']}, {inst['n_months']} mo) ==")
        all_results[inst["label"]] = run_instrument(inst)
    rows = cm.aggregate(all_results)
    score, ranking = cm.select_best(rows)
    cm.plot_summary(rows)
    plot_timeseries_grid(all_results, insts)
    cm.write_tables(rows, score, ranking)
    print("LONGRUN_DONE")


if __name__ == "__main__":
    main()
