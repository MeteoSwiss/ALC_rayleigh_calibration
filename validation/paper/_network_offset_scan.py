# -*- coding: utf-8 -*-
"""Network-scale electronic-offset scan: extract the clear-night offset for all usable CL61 (10),
10 CL51 and 10 CL31 (+ the hood-anchored Payerne CL31), fit a compact ripple model per instrument,
test the internal-temperature dependence, and emit an operational coefficient table.

Per instrument we report (validated on Payerne CL31 in _cl31_clearnight_vs_hood.py — clear-night
recovers the PERSISTENT ripple, not the near-range damped ring):
  * digitizer ripple at f_sample/4 (period = 4 range gates): the Vaisala ADC-interleave fixed pattern
  * the dominant coherent autocorrelation period (if any) and the ripple RMS in 2-8 km
  * a temperature slope (ripple RMS per deg C of internal laser temperature)
  * b_phys(range): the empirical correction (ripple x z^2, in rcs_0 units), zeroed below the
    aerosol-contaminated near range -> subtract from L1 rcs_0 to de-ripple the aerosol backscatter.

Outputs: network_offset_coeffs.csv (+ .json), per-instrument npz, and a summary figure."""
import sys, json, csv, warnings
warnings.filterwarnings("ignore")
sys.path.insert(0, r"C:/Users/hervo/OneDrive/Documents/ALC_rayleigh_calibration")
try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from validation.paper._offset_lib import (load_clearnight, offset_stats, highpass, autocorr,
                                          dominant_period, fit_ripple, gate_fold, split_half_repro, C_LIGHT)

D = Path("C:/DATA/Projects/202606_E-PROFILE_calibration/figs_paper_validation/paper_python")
OUTDIR = D / "network_offset"; OUTDIR.mkdir(exist_ok=True)
STRIDE = 2                       # every 2nd file: ~60 clear-night days is plenty for a fixed pattern
CLEAN_LO = 1800.0                # below this, boundary-layer aerosol contaminates the clear-night median
RMS_LO, RMS_HI = 2000.0, 8000.0

SAMPLE = json.loads(Path(r"C:/Users/hervo/AppData/Local/Temp/claude/"
                         r"C--Users-hervo-OneDrive-Documents-ALC-rayleigh-calibration/"
                         r"7771e6db-207a-4e4d-a515-e8c3bb4cd7c9/scratchpad/network_sample.json").read_text())
# add the hood-anchored Payerne CL31 and the Uccle CL51 reference (the detailed 40 m-ripple case)
SAMPLE["CL31"].append(["0-20000-0-06610", "B", 0, "PAYERNE(hood)"])
SAMPLE["CL51"].append(["0-20000-0-06447", "A", 0, "UCCLE(ref)"])

TEMP_VAR = {"CL61": "temperature_laser", "CL51": "temperature_laser", "CL31": "temperature_laser"}


def analyse(itype, wmo, ident, site):
    cn = load_clearnight(wmo, ident, temp_var=TEMP_VAR[itype], file_stride=STRIDE)
    if cn is None or cn["X"].shape[0] < 200:
        return None
    rng = cn["rng"]; dr = float(np.median(np.diff(rng)))
    Pmed, MAD, SE, N = offset_stats(cn["X"], rng)
    hp = highpass(Pmed, dr)
    zkm = rng / 1000.0
    # dimensionless reference signal (mid-range, present on clear nights) so CL61 (physical beta ~1e-7)
    # and CL31/CL51 (volts ~1-50) are comparable; rel_ripple = ripple as a fraction of it
    sig_ref = float(np.nanmedian(np.abs(Pmed[(rng >= 800) & (rng <= 1500)])))
    # dominant coherent period (autocorr) + ripple RMS
    Lam_dom, ac_dom, rms = dominant_period(hp, rng, dr, RMS_LO, RMS_HI, pmin=2 * dr, pmax=1500.0)
    # (2b) SPLIT-HALF REPRODUCIBILITY — the robust 'is this a real fixed pattern?' test
    repro, hpa, hpb = split_half_repro(cn["X"], rng, RMS_LO, RMS_HI)
    rel = 100 * rms / sig_ref if (sig_ref and sig_ref > 0) else np.nan
    N = int(round(Lam_dom / dr)) if np.isfinite(Lam_dom) else 0
    is_gate_ripple = (2 <= N <= 12 and abs(Lam_dom / dr - N) < 0.2 and np.isfinite(repro) and repro > 0.6
                      and np.isfinite(rel) and rel > 1.0)
    reproducible = np.isfinite(repro) and repro > 0.5 and np.isfinite(rel) and rel > 1.0
    correctable = bool(is_gate_ripple or reproducible)
    # (3) temperature slope of the ripple RMS
    T = cn["temp"].copy()
    if np.isfinite(T).any() and np.nanmedian(T) > 200:
        T = T - 273.15
    okT = np.isfinite(T); tslope = np.nan; trange = (np.nan, np.nan); tamps = None
    if okT.sum() > 500 and np.nanstd(T[okT]) > 0.5:
        q = np.nanpercentile(T[okT], [25, 50, 75])
        masks = [T <= q[0], (T > q[0]) & (T <= q[2]), T > q[2]]
        pts = []
        for mk in masks:
            if mk.sum() > 100:
                Pm, *_ = offset_stats(cn["X"][mk], rng)
                rr = np.sqrt(np.nanmean(highpass(Pm, dr)[(rng >= RMS_LO) & (rng <= RMS_HI)] ** 2))
                pts.append((float(np.nanmedian(T[mk])), float(rr)))
        if len(pts) >= 2:
            tt, aa = np.array([p[0] for p in pts]), np.array([p[1] for p in pts])
            tslope = float(np.polyfit(tt, aa, 1)[0]); trange = (float(tt.min()), float(tt.max()))
            tamps = (tt, aa)
    # (4) correction b_phys (rcs_0 units), applied at ALL altitudes:
    #   - a phase-locked digitizer ripple is reconstructed everywhere by the N-gate fold;
    #   - a reproducible slow distortion is taken as the empirical high-pass pattern over the clean
    #     free troposphere (near range needs a hood, so it is left uncorrected there);
    #   - otherwise no correction.
    if is_gate_ripple:
        ripple = gate_fold(hp, rng, N, RMS_LO, min(rng.max(), 9000.0))
        b_phys = ripple * zkm ** 2                      # ALL ranges (phase-locked)
        mode = f"gate-fold N={N}"
    elif reproducible:
        b = np.zeros_like(rng); c = (rng >= RMS_LO) & (rng <= 9000.0); b[c] = np.nan_to_num(hp[c])
        b_phys = b * zkm ** 2                            # free-troposphere distortion (empirical)
        mode = "empirical FT"
    else:
        b_phys = np.zeros_like(rng); mode = "none"
    np.savez(OUTDIR / f"{wmo}_{ident}_{itype}.npz", rng=rng, Pmed=Pmed, MAD=MAD, hp=hp, hpa=hpa,
             hpb=hpb, b_phys=b_phys, Lam_dom=Lam_dom, N=N, repro=repro, rel=rel,
             correctable=correctable, is_gate_ripple=is_gate_ripple, sig_ref=sig_ref, site=site)
    sr = sig_ref if (sig_ref and np.isfinite(sig_ref) and sig_ref > 0) else np.nan
    return dict(itype=itype, wmo=wmo, ident=ident, site=site, n=cn["n_kept"], dr=dr,
                Lam_dom=Lam_dom, ac_dom=ac_dom, ripple_rms=rms, sig_ref=sr, repro=repro,
                N=N, rel_ripple_pct=rel, correctable=correctable, is_gate_ripple=is_gate_ripple,
                mode=mode, tslope=tslope, tlo=trange[0], thi=trange[1],
                _rng=rng, _hp=hp, _tamps=tamps)


def main():
  rows = []
  for itype in ("CL61", "CL51", "CL31"):
    for wmo, ident, _nd, site in SAMPLE[itype]:
        print(f"[{itype}] {wmo}_{ident} {site} ...", flush=True)
        try:
            r = analyse(itype, wmo, ident, site)
        except Exception as exc:
            print(f"   FAILED: {exc}"); continue
        if r is None:
            print("   skipped (too few clear nights)"); continue
        rows.append(r)
        vf = "CORRECTABLE" if r["correctable"] else "-"
        print(f"   n={r['n']} dr={r['dr']:.0f}m  Λ={r['Lam_dom']:.0f}m relrms={r['rel_ripple_pct']:.1f}%  "
              f"repro={r['repro']:+.2f}  [{r['mode']}]  {vf}", flush=True)

  # ---------- compact scan CSV (the authoritative coefficient table is built by _network_offset_coeffs.py) ----------
  cols = ["itype", "site", "wmo", "ident", "n_profiles", "range_gate_m", "dominant_period_m",
          "ripple_gates_N", "rel_ripple_pct", "split_half_repro", "correction_mode",
          "temp_slope_per_C", "correctable"]
  with open(D / "network_offset_scan.csv", "w", newline="", encoding="utf-8") as f:
    w = csv.writer(f); w.writerow(cols)
    for r in rows:
        w.writerow([r["itype"], r["site"], r["wmo"], r["ident"], r["n"], f"{r['dr']:.0f}",
                    f"{r['Lam_dom']:.0f}", r["N"], f"{r['rel_ripple_pct']:.2f}", f"{r['repro']:.3f}",
                    r["mode"], f"{r['tslope']:.4g}", int(r["correctable"])])
  nc = {t: (sum(1 for r in rows if r["itype"] == t and r["correctable"]),
            sum(1 for r in rows if r["itype"] == t)) for t in ("CL51", "CL31", "CL61")}
  print(f"\nwrote {len(rows)} instruments. correctable: "
        f"CL51 {nc['CL51'][0]}/{nc['CL51'][1]}, CL31 {nc['CL31'][0]}/{nc['CL31'][1]}, CL61 {nc['CL61'][0]}/{nc['CL61'][1]}")
  print("NETWORK_OFFSET_SCAN_DONE")


if __name__ == "__main__":
    main()
