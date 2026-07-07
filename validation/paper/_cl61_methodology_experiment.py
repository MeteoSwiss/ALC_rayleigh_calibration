"""
_cl61_methodology_experiment.py — how SHOULD the 910 nm CL61 attenuated backscatter be converted
to compare against the 1064 nm CHM15k reference?  (the paper's key methodology experiment)

For every station with a CHM15k reference AND a co-located CL61, read L1 once, build the raw
calibrated beta for both instruments, then sweep a matrix of corrections on the CL61 and score the
agreement with the CHM15k over 500-3000 m AGL and per altitude band:

    WV      in {off, on}                         two-way water-vapour transmission (910 nm only)
    lambda  in {none, flat-alpha, molecular}     910 -> 1064 nm conversion
              none      : compare the 910 nm beta directly to the 1064 nm CHM (no conversion)
              flat-alpha: beta * (910/1064)^alpha           <- the CURRENT pipeline for CL61
              molecular : beta_mol*T2_mol(1064) + [beta - beta_mol*T2_mol(910)]*(910/1064)^alpha
                          i.e. the analytic Rayleigh part (lambda^-4) is handled EXACTLY and only the
                          aerosol residual carries the Angstrom law   <- the molaer treatment the
                          pipeline already uses for the Mini-MPL, applied here to the CL61
    alpha   swept in ALPHAS                        aerosol backscatter-related Angstrom exponent

Why this matters: the flat-alpha model applies one Angstrom exponent to the WHOLE signal, but
molecular backscatter scales as lambda^-4 (ratio (910/1064)^4 = 0.535) while aerosol scales as
lambda^-alpha (~0.855 for alpha=1). Where the molecular share is large (clean air, high altitude,
winter, night) the flat model leaves the molecular part ~1.6x over-scaled, biasing the CL61 high by
an amount that GROWS with altitude. Separating the molecular part removes that structure.

The treatments are applied on the hourly-gridded science matrix. This is median-exact in TIME: flat-alpha
is multiplicative and the molecular transform is affine per gate (constant in time within a bin), and the
median commutes with both. (In altitude the regrid is a nanmean applied before the molecular subtraction;
because beta_mol varies ~5%/500 m the commutation is not exact, but the error is sub-0.01 Mm^-1 sr^-1 and
identical for flat and molecular, so it does not favour either.) WV is carried as a per-gate transmission
array gridded on the SAME bins, so the WV-on and WV-off variants share an identical detection/screening
mask (an apples-to-apples matrix).

Usage: python -m validation.paper._cl61_methodology_experiment [station ...]
       (default: payerne lindenberg aosta camborne, + uccle as a 910-vs-910 control)
Outputs (paper_python/): methodology_experiment.csv, methodology_experiment.json,
       fig_method_<station>.png (per station), fig_method_summary.png (cross-station).
"""
from __future__ import annotations
import csv
import json
import os
import sys
import warnings
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "4")
try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from validation.paper import intercompare as IC
from validation.paper import overlap as OV
from validation.paper.calib_benchmark import key_of
from validation.paper.run_paper_validation import BENCHMARK, SITE, SITE_NAME
from calibration.rayleigh.atmosphere import load_standard_atmosphere, calculate_molecular_properties

OUT = Path("C:/DATA/Projects/202606_E-PROFILE_calibration/figs_paper_validation/paper_python")
ZMIN, ZMAX = 500.0, 3000.0
BANDS = [(500.0, 1000.0), (1000.0, 2000.0), (2000.0, 3000.0)]
BANDS_FINE = [(z, z + 250.0) for z in range(500, 3000, 250)]   # 10 x 250 m bands, flatness robustness
ALPHAS = [0.0, 0.5, 1.0, 1.5, 2.0]


def _relstats(cur, ref, zmask):
    """Robust agreement of cur vs ref over an altitude band. Returns (signed median rel bias %,
    MARD = median absolute relative difference %, IQR of the relative difference %, n). MARD/IQR are
    the tightening check the signed median alone cannot give: if |median| falls but MARD/IQR do not,
    the correction merely translated a skewed distribution rather than improving per-gate agreement."""
    a = cur[:, zmask]; b = ref[:, zmask]
    m = np.isfinite(a) & np.isfinite(b) & (b > 0)
    if m.sum() < 3:
        return np.nan, np.nan, np.nan, int(m.sum())
    r = (a[m] - b[m]) / b[m]
    return (float(100 * np.median(r)), float(100 * np.median(np.abs(r))),
            float(100 * (np.percentile(r, 75) - np.percentile(r, 25))), int(r.size))
# stations with a CHM15k reference AND a co-located CL61 (the wavelength-conversion cases);
# uccle is a control: reference is a 910 nm CL51, so wavelength is a no-op. It is a WV-MISMATCH probe,
# not a clean WV-cancellation test — the CL51 (910.0/3.4 nm) and CL61 (910.74/1.0 nm) laser lines have
# different two-way WV transmissions, so WV does not fully cancel (the ~12 pt WV shift there is real).
STATIONS = ["payerne", "lindenberg", "aosta", "camborne", "uccle"]

TYPE_COLORS = {"CHM15k": "#d62728", "CL61": "#1f77b4", "CL31": "#ff7f0e",
               "CL51": "#9467bd", "Mini-MPL": "#2ca02c"}


def _mol_att_abs(z_agl, station_alt, wavelength_nm):
    """Molecular attenuated backscatter beta_mol*T2_mol [Mm^-1 sr^-1] at ABSOLUTE altitude
    (station_alt + z_agl), with the two-way transmission referenced to the instrument. Corrects the
    IC._molecular_beta approximation that evaluates the US std atmosphere at AGL as if the station
    were at sea level -- which overestimates the air density (hence beta_mol) at elevated sites and
    over-subtracts the molecular part aloft. T2_mol(station->z) = T2(0->z) / T2(0->station)."""
    grid = np.arange(0, 20001, 30.0)                      # ASL grid, high enough for salt+3 km
    atm = load_standard_atmosphere(IC.STD_ATM, grid)
    mol = calculate_molecular_properties(atm.temperature, atm.pressure, grid, wavelength_nm * 1e-9)
    z_abs = station_alt + z_agl
    beta = np.interp(z_abs, grid, mol.beta_mol, left=np.nan, right=np.nan)
    t2 = np.interp(z_abs, grid, mol.transmission, left=np.nan, right=np.nan)
    t2_st = float(np.interp(station_alt, grid, mol.transmission))
    return beta * (t2 / t2_st) * 1e6


# --------------------------------------------------------------------------- per-channel gridding
def _calibrated_beta0(l1, ch):
    """Raw calibrated beta = rcs_0 / C_L * 1e6 [Mm^-1 sr^-1], BEFORE WV and wavelength. CHM15k also
    gets the temperature-dependent overlap correction (as in the pipeline). None if no calib series."""
    beta = l1["beta"].copy()
    if ch["itype"] in ("CHM15k", "CHM8k"):
        model = OV.load_overlap_model(ch["wmo"], ch["ident"])
        if model is not None:
            beta = OV.correct_rcs(beta, l1["alt"] - l1["station_alt"], l1["temp_int"], model)
    cal = IC.load_calib_series(key_of(ch), "L1")
    if cal is None:
        return None
    ck = IC.interp_calib(cal[0], cal[1], l1["time"])
    return beta / ck[:, None] * 1e6


def _grid_reference(l1, ch):
    """Reference science stream: raw calibrated -> (WV if 910 nm) -> screened -> hourly grid
    (SNR-gated). Returns (grid_time, beta_grid, alt) or None. A 1064 nm CHM reference gets no WV;
    a 910 nm reference (Uccle CL51) is WV-corrected so a 910-vs-910 comparison cancels WV."""
    b0 = _calibrated_beta0(l1, ch)
    if b0 is None:
        return None
    if ch["itype"] in ("CL31", "CL51", "CL61"):
        lam0, fwhm = IC.WV_PARAMS.get(ch["itype"], (910.0, 3.4))
        b0, _ = IC.apply_wv(b0, l1, lam0, fwhm)
    scr = IC.screen(b0, l1)
    grid, (bg,) = IC.retime_hourly(l1["time"], [scr], min_cov_s=IC.MIN_AVG_S, snr_idx=(0,))
    return grid, bg, l1["alt"]


def _grid_cl61(l1, ch):
    """CL61 science stream gridded, carrying the WV transmission so WV-on/off share one mask.
    Returns dict(grid, M_wv, T2, alt) where M_wv = gridded WV-corrected beta (science, SNR-gated)
    and T2 = gridded two-way WV transmission on the same bins; M_raw (WV-off) = M_wv * T2."""
    b0 = _calibrated_beta0(l1, ch)
    if b0 is None:
        return None
    lam0, fwhm = IC.WV_PARAMS.get(ch["itype"], (910.0, 3.4))
    beta_wv, info = IC.apply_wv(b0, l1, lam0, fwhm)
    with np.errstate(all="ignore"):
        t2 = b0 / beta_wv                      # per-gate two-way WV transmission (NaN where WV-masked)
    scr = IC.screen(beta_wv, l1)               # screen AFTER WV, as in the pipeline
    grid, (M_wv, T2) = IC.retime_hourly(l1["time"], [scr, t2], min_cov_s=IC.MIN_AVG_S, snr_idx=(0,))
    return dict(grid=grid, M_wv=M_wv, T2=T2, alt=l1["alt"],
                months_excluded=info.get("months_excluded", []))


# --------------------------------------------------------------------------- treatment matrix
def _apply_treatment(base, wl, alpha, mol, lam, target):
    """Apply a wavelength/molecular treatment to a gridded (time x alt) matrix. mol maps a molecular
    treatment name to (bml, bmt) = molecular attenuated backscatter at (lam, target) on the alt grid.
      flat          : base * (lam/target)^alpha              -- single Angstrom on the TOTAL signal
      molecular     : analytic Rayleigh (US std atm at AGL) removed, Angstrom on the aerosol residual
      molecular_abs : same but the Rayleigh profile at the station's ABSOLUTE altitude (recommended)"""
    if wl == "none" or abs(lam - target) < 1.0:
        return base
    f = (lam / target) ** alpha                 # aerosol Angstrom factor (0.855 for a=1 at 910->1064)
    if wl == "flat":
        return base * f
    if wl in mol:
        bml, bmt = mol[wl]
        return bmt[None, :] + (base - bml[None, :]) * f
    raise ValueError(wl)


def run_station(name):
    st = BENCHMARK[name]
    sc = SITE[name]
    target = sc["target"]
    channels = st["channels"]
    ref_ch = channels[sc["ref"]]
    cl61_chs = [c for c in channels if c["itype"] == "CL61"]
    if not cl61_chs:
        print(f"  {name}: no CL61 channel"); return None

    l1_cache = {}
    def getl1(ch):
        k = (ch["wmo"], ch["ident"])
        if k not in l1_cache:
            l1_cache[k] = IC.read_l1(ch["wmo"], ch["ident"], st["start"], st["end"])
        return l1_cache[k]

    ref_l1 = getl1(ref_ch)
    if ref_l1 is None:
        print(f"  {name}: no L1 for reference {ref_ch['label']}"); return None
    refg = _grid_reference(ref_l1, ref_ch)
    if refg is None:
        print(f"  {name}: no calib for reference {ref_ch['label']}"); return None
    ref_grid, ref_M, ref_alt = refg

    # grid every CL61 method
    cl61_grids = []
    for ch in cl61_chs:
        l1 = getl1(ch)
        if l1 is None:
            print(f"    [skip] {ch['label']}: no L1"); continue
        g = _grid_cl61(l1, ch)
        if g is None:
            print(f"    [skip] {ch['label']}: no calib series ({key_of(ch)})"); continue
        g["ch"] = ch
        cl61_grids.append(g)
    if not cl61_grids:
        return None

    # common altitude grid across the reference and all CL61 streams
    altGrid = IC.build_common_grid([ref_alt] + [g["alt"] for g in cl61_grids])
    salt = ref_l1["station_alt"]
    z_agl = altGrid - salt
    zf = (z_agl >= ZMIN) & (z_agl <= ZMAX)
    band_masks = [(z_agl >= lo) & (z_agl <= hi) for lo, hi in BANDS]
    fine_masks = [(z_agl >= lo) & (z_agl <= hi) for lo, hi in BANDS_FINE]
    # molecular attenuated backscatter at the CL61 line and the target, on the common alt grid:
    # AGL-from-sea-level (as the pipeline's molaer does) and at the station's absolute altitude.
    lam_cl = IC.WV_PARAMS.get("CL61", (910.74, 1.0))[0]
    mol = {"molecular":     (IC._molecular_beta(z_agl, salt, lam_cl), IC._molecular_beta(z_agl, salt, target)),
           "molecular_abs": (_mol_att_abs(z_agl, salt, lam_cl),       _mol_att_abs(z_agl, salt, target))}

    rows = []
    for g in cl61_grids:
        ch = g["ch"]
        # union time grid of this CL61 with the reference; reindex both
        union = np.unique(np.concatenate([ref_grid, g["grid"]]))
        def onto(grid, M):
            idx = {t: i for i, t in enumerate(grid)}
            pos = np.array([idx.get(t, -1) for t in union])
            out = np.full((union.size, M.shape[1]), np.nan)
            ok = pos >= 0
            out[ok] = M[pos[ok]]
            return out
        refU = onto(ref_grid, ref_M)
        M_wvU = onto(g["grid"], g["M_wv"])
        T2U = onto(g["grid"], g["T2"])
        with np.errstate(all="ignore"):
            M_rawU = M_wvU * T2U                              # WV-off base
        refCg = IC.regrid(refU, ref_alt, altGrid)
        base_on = IC.regrid(M_wvU, g["alt"], altGrid)
        base_raw = IC.regrid(M_rawU, g["alt"], altGrid)

        for wv, base in (("off", base_raw), ("on", base_on)):
            wl_opts = (["none"] if abs(lam_cl - target) < 1.0
                       else ["none", "flat", "molecular", "molecular_abs"])
            for wl in wl_opts:
                alphas = [np.nan] if wl == "none" else ALPHAS
                for a in alphas:
                    T = _apply_treatment(base, wl, (0.0 if wl == "none" else a),
                                         mol, lam_cl, target)
                    sf = IC._stats(T, refCg, zf)
                    med, mard, iqr, nrel = _relstats(T, refCg, zf)
                    bands = [_relstats(T, refCg, bm)[0] for bm in band_masks]
                    fine = [_relstats(T, refCg, fm)[0] for fm in fine_masks]
                    rows.append(dict(
                        station=name, calib=ch["calib"], itype=ch["itype"],
                        wv=wv, wl=wl, alpha=(None if wl == "none" else a),
                        medrel=med, mard=mard, iqr=iqr, logr=sf["r_log"], n=sf["n"],
                        band_lo=bands[0], band_mid=bands[1], band_hi=bands[2],
                        band_spread=float(np.nanstd(bands)),
                        fine_spread=float(np.nanstd(fine))))
        print(f"    {ch['label']}: {sum(r['calib']==ch['calib'] for r in rows)} treatments"
              f"  (WV-excluded months: {len(g['months_excluded'])})", flush=True)
    return dict(name=name, rows=rows, target=target, lam_cl=lam_cl)


# --------------------------------------------------------------------------- figures
def _pick(rows, **kw):
    for r in rows:
        if all(r.get(k) == v for k, v in kw.items()):
            return r
    return None


def fig_station(res, out_png):
    name = res["name"]; rows = res["rows"]
    calibs = sorted({r["calib"] for r in rows})
    ncol = len(calibs)
    fig, axes = plt.subplots(1, 2 * ncol, figsize=(7.5 * ncol, 5.0), squeeze=False)
    for ci, calib in enumerate(calibs):
        rr = [r for r in rows if r["calib"] == calib]
        axL = axes[0][2 * ci]; axR = axes[0][2 * ci + 1]
        # LEFT: med relbias vs alpha, flat vs molecular vs molecular-abs (WV on); reference lines
        for wl, col, mk in (("flat", "#ff7f0e", "o"), ("molecular", "#9467bd", "^"),
                            ("molecular_abs", "#1f77b4", "s")):
            xs = [r for r in rr if r["wv"] == "on" and r["wl"] == wl]
            if xs:
                xs = sorted(xs, key=lambda r: r["alpha"])
                axL.plot([r["alpha"] for r in xs], [r["medrel"] for r in xs],
                         mk + "-", color=col, label=f"{wl} (WV on)")
        none_on = _pick(rr, wv="on", wl="none")
        none_off = _pick(rr, wv="off", wl="none")
        cur = _pick(rr, wv="on", wl="flat", alpha=1.0)      # current pipeline cell
        for r, lab, col, ls in ((none_off, "no conv, WV off (native 910)", "#7f7f7f", ":"),
                                (none_on, "no conv, WV on", "#111111", "--"),
                                (cur, "current pipeline (flat a=1, WV on)", "#d62728", "-.")):
            if r and np.isfinite(r["medrel"]):
                axL.axhline(r["medrel"], color=col, ls=ls, lw=1.2, label=f"{lab}: {r['medrel']:+.1f}%")
        axL.axhline(0, color="k", lw=0.8)
        axL.set_xlabel("aerosol Ångström exponent α")
        axL.set_ylabel("CL61 − CHM15k median rel. bias  [%]")
        axL.set_title(f"{SITE_NAME.get(name, name)} — CL61 ({calib})", fontsize=11, fontweight="bold")
        axL.grid(alpha=0.3); axL.legend(fontsize=7, loc="best")
        # RIGHT: altitude-band bias, current vs proposed molecular-abs alpha=1 (WV on)
        prop = _pick(rr, wv="on", wl="molecular_abs", alpha=1.0)
        band_c = [(lo + hi) / 2 / 1000 for lo, hi in BANDS]
        if cur:
            axR.plot([cur["band_lo"], cur["band_mid"], cur["band_hi"]], band_c,
                     "o-", color="#d62728", label=f"current flat α=1: {cur['medrel']:+.1f}%")
        if prop:
            axR.plot([prop["band_lo"], prop["band_mid"], prop["band_hi"]], band_c,
                     "s-", color="#1f77b4", label=f"molecular-abs α=1: {prop['medrel']:+.1f}%")
        axR.axvline(0, color="k", lw=0.8)
        axR.set_xlabel("median rel. bias in band  [%]")
        axR.set_ylabel("altitude AGL  [km]")
        axR.set_title("altitude structure: does it flatten?", fontsize=10)
        axR.grid(alpha=0.3); axR.legend(fontsize=8, loc="best")
    fig.suptitle(f"910→1064 nm comparison methodology — {SITE_NAME.get(name, name)}",
                 fontsize=13, fontweight="bold")
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    fig.savefig(out_png, dpi=170); plt.close(fig)
    print(f"   -> {out_png.name}", flush=True)


def fig_summary(all_res, out_png):
    """Cross-station: residual bias + band spread, current (flat α=1) vs proposed (molecular α=1),
    for the CL61 cloud method (the operational one)."""
    labels, cur_b, prop_b, cur_s, prop_s = [], [], [], [], []
    for res in all_res:
        rr = [r for r in res["rows"] if r["calib"] == "cloud"]
        if not rr:
            continue
        cur = _pick(rr, wv="on", wl="flat", alpha=1.0)
        prop = _pick(rr, wv="on", wl="molecular_abs", alpha=1.0)
        if not (cur and prop):
            continue
        labels.append(SITE_NAME.get(res["name"], res["name"]))
        cur_b.append(cur["medrel"]); prop_b.append(prop["medrel"])
        cur_s.append(cur["band_spread"]); prop_s.append(prop["band_spread"])
    x = np.arange(len(labels)); w = 0.38
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))
    ax1.bar(x - w / 2, cur_b, w, color="#d62728", label="current (flat α=1, WV on)")
    ax1.bar(x + w / 2, prop_b, w, color="#1f77b4", label="proposed (molecular α=1, WV on)")
    ax1.axhline(0, color="k", lw=0.8); ax1.set_xticks(x); ax1.set_xticklabels(labels)
    ax1.set_ylabel("CL61 − CHM15k median rel. bias  [%]")
    ax1.set_title("Residual bias per station", fontsize=11, fontweight="bold")
    ax1.grid(axis="y", alpha=0.3); ax1.legend(fontsize=9)
    ax2.bar(x - w / 2, cur_s, w, color="#d62728", label="current")
    ax2.bar(x + w / 2, prop_s, w, color="#1f77b4", label="proposed")
    ax2.set_xticks(x); ax2.set_xticklabels(labels)
    ax2.set_ylabel("across-band bias spread  [%]  (lower = flatter)")
    ax2.set_title("Altitude flatness per station", fontsize=11, fontweight="bold")
    ax2.grid(axis="y", alpha=0.3); ax2.legend(fontsize=9)
    fig.suptitle("Current flat-α vs proposed molecular-aware 910→1064 nm conversion",
                 fontsize=13, fontweight="bold")
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    fig.savefig(out_png, dpi=170); plt.close(fig)
    print(f"   -> {out_png.name}", flush=True)


# --------------------------------------------------------------------------- main
def _process_station(name):
    """ProcessPool worker: run one station's treatment matrix + its figure; return (name, res)."""
    warnings.filterwarnings("ignore")
    res = run_station(name)
    if res is None:
        return name, None
    fig_station(res, OUT / f"fig_method_{name}.png")
    return name, res


def main():
    warnings.filterwarnings("ignore")
    req = [a for a in sys.argv[1:] if a in BENCHMARK]
    stations = req or STATIONS
    all_rows, all_res = [], []
    workers = max(1, min(len(stations), int(os.environ.get("ALC_VAL_SITE_WORKERS", "5"))))
    results = {}
    if workers > 1 and len(stations) > 1:
        print(f"[parallel] {len(stations)} stations over {workers} workers", flush=True)
        with ProcessPoolExecutor(max_workers=workers) as ex:
            for name, res in ex.map(_process_station, stations):
                results[name] = res
    else:
        for name in stations:
            n, res = _process_station(name)
            results[n] = res
    for name in stations:
        res = results.get(name)
        if res is None:
            print(f"== {name} ==   no result", flush=True)
            continue
        all_res.append(res)
        all_rows.extend(res["rows"])
    if not all_rows:
        print("no rows"); return
    cols = ["station", "calib", "itype", "wv", "wl", "alpha", "medrel", "mard", "iqr", "logr", "n",
            "band_lo", "band_mid", "band_hi", "band_spread", "fine_spread"]
    with open(OUT / "methodology_experiment.csv", "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=cols); w.writeheader()
        for r in all_rows:
            w.writerow({c: r.get(c, "") for c in cols})
    with open(OUT / "methodology_experiment.json", "w", encoding="utf-8") as f:
        json.dump(all_rows, f, indent=1, default=lambda o: None if (isinstance(o, float) and not np.isfinite(o)) else o)
    fig_summary(all_res, OUT / "fig_method_summary.png")
    print(f"   -> methodology_experiment.csv ({len(all_rows)} rows)")
    print("METHODOLOGY_EXPERIMENT_DONE", flush=True)


if __name__ == "__main__":
    main()
