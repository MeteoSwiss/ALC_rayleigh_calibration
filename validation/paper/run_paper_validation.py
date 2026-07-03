"""
run_paper_validation.py — operational Python benchmark validation, end to end. Runs the inter-comparison
(intercompare.process) for each multi-instrument station using the Python calibration+Kalman series from
calib_benchmark.py, renders the MATLAB-style figures (multi-ALC 3x3 panel, combined calibration time-series,
EARLINET 2x2 panel) via figures.py, and writes a report with the Python statistics next to the MATLAB R_*.mat.

Usage:  python -m validation.paper.run_paper_validation
"""
from __future__ import annotations
import sys
from pathlib import Path

import numpy as np
import scipy.io as sio

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass

from validation.paper import intercompare as IC
from validation.paper import earlinet as EA
from validation.paper import figures as FIG
from validation.paper.calib_benchmark import BENCHMARK, key_of

REPO = Path(__file__).resolve().parents[2]
OUT = Path("C:/DATA/Projects/202606_E-PROFILE_calibration/figs_paper_validation/paper_python")
OUT.mkdir(parents=True, exist_ok=True)
CALIB = OUT / "calib"
MAT = Path("C:/Users/hervo/OneDrive/Documents/MATLAB/ALC/figs_paper_validation")
SITE_NAME = {"payerne": "Payerne", "amsterdam": "Amsterdam", "uccle": "Uccle", "sirta": "Palaiseau",
             "lindenberg": "Lindenberg", "aosta": "Aosta", "camborne": "Camborne", "earlinet": ""}
# station -> (referenceChannel, lambda_target, MATLAB R file, molaer label, WIGOS for title)
STATIONS = {
    "payerne":   dict(ref=0, target=1064.0, mat="R_payerne.mat",    molaer=None, wmo="0-20000-0-06610"),
    "amsterdam": dict(ref=0, target=1064.0, mat="R_amsterdam.mat",  molaer=None, wmo="0-20000-0-06240"),
    "uccle":     dict(ref=0, target=910.0,  mat="R_cl51_06447.mat", molaer=None, wmo="0-20000-0-06447"),
    "sirta":     dict(ref=0, target=1064.0, mat="R_sirta.mat",      molaer="Mini-MPL (Rayleigh)", wmo="0-250-1001-07151"),
    # new CL61+CHM15k pairs (no MATLAB reference -> mat file absent -> load_matlab returns {})
    # lindenberg: the CL61 L2 product is not in the local archive -> read native L1 (rcs_0 / C_L)
    "lindenberg": dict(ref=0, target=1064.0, mat="R_lindenberg.mat", molaer=None, wmo="0-20000-0-10393",
                       dataLevel="L1"),
    "aosta":      dict(ref=0, target=1064.0, mat="R_aosta.mat",      molaer=None, wmo="0-380-5-1"),
    "camborne":   dict(ref=0, target=1064.0, mat="R_camborne.mat",   molaer=None, wmo="0-20000-0-03808"),
}


def run_station(name):
    st = BENCHMARK[name]; sc = STATIONS[name]
    chans = []
    for c in st["channels"]:
        d = dict(wmo=c["wmo"], ident=c["ident"], calib=c["calib"], label=c["label"], itype=c["itype"], key=key_of(c))
        if sc["molaer"] and c["label"] == sc["molaer"]:
            d["wavelengthModel"] = "molaer"
        chans.append(d)
    if name == "payerne":
        # second CL61 Rayleigh entry: constants recalibrated after removing the hood-measured
        # dark offset (cl61_rayleigh_investigation.md 1c) — shown alongside the native series
        chans.append(dict(wmo="0-20000-0-06610", ident="C", calib="rayleigh",
                          label="CL61 (Rayleigh, offset-corr)", itype="CL61",
                          key="0-20000-0-06610_C_rayleigh_offsetcorr"))
        # CHM15k Rayleigh recalibrated after removing its (photon-counting over-subtraction) hood
        # offset — to test whether the correction improves the multi-instrument agreement
        # (cl61_chm15k_offset_correction.md). NB the CHM15k offset is within its calibration scatter.
        chans.append(dict(wmo="0-20000-0-06610", ident="A", calib="rayleigh",
                          label="CHM15k (Rayleigh, offset-corr)", itype="CHM15k",
                          key="0-20000-0-06610_A_rayleigh_offsetcorr"))
    cfg = dict(wmo=sc["wmo"], start=st["start"], end=st["end"], referenceChannel=sc["ref"],
               channels=chans, lambda_target=sc["target"], alpha=1.0, zMin=500, zMax=3000,
               calibLevel="L1")   # native L1 (binned, eprof_v2 + fixed cloud); falls back to L2 where no L1
    if sc.get("dataLevel"):
        cfg["dataLevel"] = sc["dataLevel"]   # read native L1 rcs_0 and apply the calout C_L directly
    return IC.process(cfg), cfg


def load_matlab(mat):
    f = MAT / mat
    if not f.is_file():
        return {}
    m = sio.loadmat(str(f), squeeze_me=True, struct_as_record=False)
    Rm = m["R"]; out = {}
    for c, s in zip(np.atleast_1d(Rm.channels), np.atleast_1d(Rm.stats)):
        out[str(c.label)] = dict(relbias=float(getattr(s, "relbias_pct", np.nan)),
                                 r=float(getattr(s, "r", np.nan)), n=int(getattr(s, "n", 0)))
    return out


def calib_channel_list():
    """One entry per INSTRUMENT for the combined C_L time-series grid: all its calibration
    methods (CL61: Rayleigh + cloud) are drawn in the same panel, site-labelled."""
    order, seen = [], {}
    for name, st in BENCHMARK.items():
        for c in st["channels"]:
            stream = f"{c['wmo']}_{c['ident']}"
            k = key_of(c)
            if stream in seen:
                if k not in [s["key"] for s in seen[stream]["series"]]:
                    seen[stream]["series"].append(dict(key=k, calib=c["calib"], itype=c["itype"]))
                continue
            site = SITE_NAME.get(name) or c["label"].split(" ")[0]
            base = c["label"].split(" (")[0]     # drop the "(cloud)/(Rayleigh)" method suffix
            title = f"{site} {base}" if SITE_NAME.get(name) else base
            entry = dict(title=title, unit="C$_L$", series=[dict(key=k, calib=c["calib"], itype=c["itype"])])
            seen[stream] = entry
            order.append(entry)
    return order


def main():
    rows = []
    for name, sc in STATIONS.items():
        print(f"== {name} ==", flush=True)
        R, cfg = run_station(name)
        if R is None:
            print("  no data"); continue
        mat = load_matlab(sc["mat"])
        for k, ch in enumerate(R["channels"]):
            s = R["stats"][k]; mm = mat.get(ch["label"], {})
            rows.append(dict(station=name, label=ch["label"], calib=ch["calib"], ref=(k == cfg["referenceChannel"]),
                             py_relbias=s["relbias_pct"], py_medrel=s.get("medrelbias_pct", np.nan),
                             py_r=s["r"], py_rlog=s.get("r_log", np.nan), py_n=s["n"],
                             mat_relbias=mm.get("relbias", np.nan), mat_r=mm.get("r", np.nan), mat_n=mm.get("n", 0)))
            print("   %-20s PY relbias=%+7.1f%% (med %+6.1f%%) r=%.3f (log %.3f) N=%7d | MAT relbias=%+7.1f%% r=%.3f"
                  % (ch["label"], s["relbias_pct"], s.get("medrelbias_pct", np.nan), s["r"],
                     s.get("r_log", np.nan), s["n"], mm.get("relbias", np.nan), mm.get("r", np.nan)), flush=True)
        title = "%s (%s)  —  %s to %s" % (SITE_NAME[name], sc["wmo"], _d(cfg["start"]), _d(cfg["end"]))
        FIG.fig_multi_alc(R, cfg, OUT / f"fig_{name}.png", title)
        print(f"   -> fig_{name}.png", flush=True)

    # combined calibration time-series (all channels) — the operational dashboard L1 2025-2026 series
    FIG.fig_calib_timeseries(calib_channel_list(), CALIB, OUT / "fig_calib_timeseries.png")
    print("   -> fig_calib_timeseries.png (dashboard L1 2025-2026, eprof_v2 + cloud)", flush=True)

    # EARLINET 2x2 figures
    erows = []
    for code in ("sir", "ino", "ari", "lei", "cbw", "sir_532"):
        try:
            s = EA.compare(code, "20250101", "20260630", return_profiles=True)
        except Exception as exc:
            s = {"error": repr(exc)}
        mm = EA.load_matlab_earlinet(code)
        label = {"sir": "Palaiseau", "ino": "Magurele", "ari": "Leipzig", "lei": "Leipzig",
                 "cbw": "Cabauw", "sir_532": "Palaiseau 532 nm"}[code]
        if s and "error" not in s:
            erows.append((code, label, s, mm))
            site = EA.SITES[code]
            FIG.fig_earlinet(code, label, s["betaE"], s["betaC"], s["grid"], s["times"], s, mm,
                             OUT / f"fig_earlinet_{code}.png", betaC_raw=s.get("betaC_raw"),
                             instr=site.get("instr", "CHM15k (Rayleigh)"),
                             itype=site.get("itype", "CHM15k"))
            print("   EARLINET %s: relbias=%+.1f%% (med %+.1f%%) r=%.2f (log %.2f) matched=%d -> fig_earlinet_%s.png"
                  % (code, s["relbias_pct"], s.get("medrelbias_pct", np.nan), s["r"],
                     s.get("r_log", np.nan), s["matched"], code), flush=True)
        else:
            erows.append((code, label, s, mm)); print("   EARLINET %s: %s" % (code, s), flush=True)

    write_report(rows, erows)
    print("PAPER_VALIDATION_DONE", flush=True)


def write_report(rows, erows):
    L = ["# Operational Python attenuated-backscatter validation — benchmark stations\n",
         "*Generated by `validation/paper/run_paper_validation.py`. The calibration applied here is the "
         "operational **dashboard** calibration: every channel is calibrated per night (Rayleigh, `eprof_v2`) "
         "or per day (liquid-cloud, O'Connor/Hopkin) from the native **Level-1** archive over **2025-2026** and "
         "smoothed with the E-PROFILE Kalman filter — the same `fullcal_l1_2026` series that feeds the "
         "monitoring dashboard (`scripts/run_all_l1_2026.py`). Water-vapour (910 nm) and wavelength corrections, "
         "screening, gridding and statistics are the operational **Python** routines. The figures reproduce the "
         "MATLAB layouts (`make_validation_figures.m`); the MATLAB columns are kept as a **legacy reference** "
         "(the earlier full-archive MATLAB calibration).*\n",
         "> The Python **relative-bias** vs MATLAB is now a calibration-method difference (operational L1 Kalman "
         "vs the legacy MATLAB recompute), not a port-fidelity check; the correlation **r** and profile count "
         "**N** still track the MATLAB closely. **Cloud** β_att uses the O'Connor multiplier on the physical "
         "attenuated backscatter (`C·attbsc_0`). **CL61 is calibrated twice** (Rayleigh + cloud) at Payerne and "
         "Uccle. Note the **Payerne CL61 Rayleigh** series is thin (only ~6 successful nights in the 2025-2026 "
         "L1 record — the native-signal molecular fit rarely finds an eligible window), so its constant is "
         "largely Kalman-predicted; the CL61 cloud calibration is robust.\n",
         "*Metrics: **relbias** = 100·mean(channel−ref)/mean(ref); **med relbias** = "
         "100·median((channel−ref)/ref) over ref>0; **r** = Pearson on linear β_att; **log r** = Pearson on "
         "log₁₀ β_att over positive pairs. The linear moments are dominated by the rare large aerosol/cloud "
         "values while the band also contains the molecular floor, so the log-space r and the median relative "
         "bias are the robust indicators. Panel (a) of each station figure shows medians restricted to the "
         "**common hours** where every channel reports (N in the panel title), so the profiles describe the "
         "same atmospheric sample.*\n",
         "## Calibration coefficient time series (all channels)\n",
         "![calibration time series](figs_paper_validation/paper_python/fig_calib_timeseries.png)\n",
         "## Per-station validation\n",
         "| station | channel | calib | Python relbias | med relbias | Python r | log r | Python N | MATLAB relbias | MATLAB r |",
         "|---|---|---|---|---|---|---|---|---|---|"]
    for r in rows:
        tag = " *(ref)*" if r["ref"] else ""
        L.append("| %s | %s%s | %s | %+.1f%% | %+.1f%% | %.3f | %.3f | %d | %+.1f%% | %.3f |"
                 % (r["station"], r["label"], tag, r["calib"], r["py_relbias"], r["py_medrel"],
                    r["py_r"], r["py_rlog"], r["py_n"], r["mat_relbias"], r["mat_r"]))
    L.append("")
    for name in STATIONS:
        L.append(f"![{name} validation](figs_paper_validation/paper_python/fig_{name}.png)\n")
    # Mini-MPL wavelength note
    L.append("## Mini-MPL 532 nm -> 1064 nm wavelength conversion (Palaiseau)\n")
    L.append("The Palaiseau Mini-MPL operates at **532 nm** and is compared to the **1064 nm** CHM15k "
             "reference, so its β_att is converted 532->1064 with the **molaer** model "
             "(`intercompare.wavelength_correct`): β_aer = β_total − β_mol·T²_mol(532) (ATTENUATED "
             "molecular), scaled by (532/1064)^(−α), recombined with the attenuated 1064 nm molecular. "
             "A single Ångström exponent on the whole signal would scale the molecular by λ⁻¹ instead of "
             "λ⁻⁴ (8× too large in clean air). CAVEAT: in clean air the extraction subtracts two "
             "nearly-equal numbers (aerosol ≈ 2 % of the 532 nm signal above the SNR gate), amplifying "
             "any residual 532 nm scale/model error ≈ 13× into the converted value — the ≈ −36 % vs the "
             "CHM15k measures this conditioning, NOT the instrument: at native 532 nm the Mini-MPL agrees "
             "with the EARLINET SIRTA lidar to **−1.9 %** (see the `sir_532` comparison).\n")
    # EARLINET
    L.append("## EARLINET — ceilometer (CHM15k) vs EARLINET research-lidar reference\n")
    L.append("*The CHM stream is screened like the station intercomparison (quality flag, clouds via CBH, "
             "fog, ±15 min expansion); EARLINET profiles are SCC cloud-screened. EARLINET gates below the "
             "instrument overlap are excluded (not filled); the transmission integral extends the lowest "
             "trusted extinction to the ground; the lidar ratio is the per-scene SCC assumption (fallback 50 sr).*\n")
    L.append("| site | Python relbias | med relbias | Python r | log r | matched | MATLAB relbias | MATLAB r |")
    L.append("|---|---|---|---|---|---|---|---|")
    for code, label, s, mm in erows:
        if s and "error" not in s:
            L.append("| %s (%s) | %+.1f%% | %+.1f%% | %.2f | %.2f | %d | %+.1f%% | %.2f |"
                     % (code, label, s["relbias_pct"], s.get("medrelbias_pct", np.nan), s["r"],
                        s.get("r_log", np.nan), s["matched"], mm.get("relbias", np.nan), mm.get("r", np.nan)))
        else:
            L.append("| %s (%s) | no EARLINET 1064 data in the 2025-2026 window | | | | | | |" % (code, label))
    L.append("")
    for code, label, s, mm in erows:
        if s and "error" not in s:
            L.append(f"![earlinet {code}](figs_paper_validation/paper_python/fig_earlinet_{code}.png)\n")
    txt = "\n".join(L)
    (OUT / "paper_python_validation.md").write_text(txt, encoding="utf-8")
    (REPO / "doc" / "reports" / "paper_python_validation.md").write_text(txt, encoding="utf-8")
    print(f"  wrote paper_python_validation.md ({len(rows)} channels, {len(erows)} EARLINET sites)")


def _d(yyyymmdd):
    from datetime import datetime
    return datetime.strptime(yyyymmdd, "%Y%m%d").strftime("%Y-%m-%d")


if __name__ == "__main__":
    main()
