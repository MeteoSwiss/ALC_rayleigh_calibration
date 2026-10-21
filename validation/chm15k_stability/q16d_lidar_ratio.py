# -*- coding: utf-8 -*-
"""q16d: lidar ratio and aerosol-attenuation correction of the CHM15k Rayleigh constant (report 19, q16, agent d).

Question: can a wrong aerosol lidar ratio (options.json LRaer = 52 sr) or an error in the aerosol two-way
transmission correction below the molecular window create (or hide) the solstice-phased ~+-9 % annual cycle of
the CHM15k constant C_L that survives the kappa x T_int normalisation?

What the operational chain does (calibration/rayleigh, read 2026-10-20, see q16d_summary.md for the equations):
  * Rayleigh fit of the night-mean signal against the US-1976 molecular signal -> window [zb, zt] and slope a.
  * X = rcs / a (attenuated-backscatter scale; C_L does NOT depend on a, the Klett solution is scale-free).
  * Backward Klett-Fernald inversion from the window MID-POINT (reference: beta_tot(ref) = beta_mol(ref) x
    mean_window(X/beta_mol) / (X/beta_mol)(ref), i.e. R ~ 1 on average in the window) DOWN TO THE FIRST RANGE GATE
    (20 m) and up to the window top, LR = 52 sr, aerosol extinction = LR x max(beta_aer, 0).
  * C_L(z) = rcs(z) / beta_tot(z) x exp(2 tau_tot(r0 -> z)) (trapezoid from the first gate, tau_tot = aerosol + 8pi/3
    beta_mol), median over the window; final constant = median of 5 LR (LRaer-20..+20) x 5 window shifts
    (-200..+200 m) x 5 time samples.

This script re-runs the operational Rayleigh calibration (native L1, options.json, monkeypatched only to CAPTURE
the inputs of the Klett step; nothing is written outside q16/) for ~22 AERONET-matched CHM15k units x all their
operational nights, then re-evaluates the constant for LR = 10..100 sr, without any aerosol extinction, and with
the near-range (incomplete-overlap) extinction replaced by a constant extrapolation from 300/600/1000 m.
It then (i) converts LR scenarios (seasonal LR, marine LR, random night-to-night LR) into a predicted C_L bias
series (annual harmonic amplitude/phase, 31-day anomalies) and (ii) compares the Klett aerosol optical depth with
AERONET AOD (1020 nm -> 1064 nm with the 440-870 Angstrom exponent) to estimate the missing extinction.

Run:   python validation/chm15k_stability/q16d_lidar_ratio.py extract   (4 workers, resumable)
       python validation/chm15k_stability/q16d_lidar_ratio.py analyse
Output: C:/DATA/Projects/202610_CHM15k_stability/q16/q16d_*
"""
from __future__ import annotations

import contextlib
import csv
import io
import json
import logging
import sys
import tempfile
import warnings
from datetime import datetime, timedelta
from multiprocessing import Pool
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import CENSUS, OUT, load_constants  # noqa: E402

warnings.filterwarnings("ignore")
Q16 = OUT / "q16"
NIGHT_DIR = Q16 / "q16d_nights"
SCRATCH = Path("C:/Users/hervo/AppData/Local/Temp/claude/C--Users-hervo-OneDrive-Documents-ALC-rayleigh-calibration/"
               "e6120125-e8c9-4e0b-8352-26777f3067de/scratchpad")
AERONET_DIR = SCRATCH / "aeronet_raw"
EPOCH = datetime(1970, 1, 1)
W_YEAR = 2 * np.pi / 365.25
N_WORKERS = 4

LR_OP = 52.0
LR_GRID = np.array([10, 20, 30, 40, 52, 60, 70, 80, 100], float)     # single-config (nominal) evaluation
LR_CENTRES = np.array([20, 30, 40, 52, 60, 70, 80], float)          # full operational ensemble around these
LR_DELTAS = (-20, -10, 0, 10, 20)
ALT_SHIFTS = (-200, -100, 0, 100, 200)
NEAR_FILL_M = (300.0, 600.0, 1000.0)                                 # constant-extrapolation heights
TAU_LAYERS = ((0, 200), (200, 500), (500, 1000), (1000, 1500), (1500, 1e9))
PROF_EDGES = np.arange(0, 4001, 100.0)                                # 100 m bins for the stored profiles

# AERONET-matched CHM15k units (q5_aeronet.csv) -> geographically spread, coastal + continental + alpine.
UNITS = ["0-100-20000-0000001_A", "0-20000-0-01492_A", "0-20000-0-03354_B", "0-20000-0-03808_A",
         "0-20000-0-03838_A", "0-20000-0-06348_A", "0-20000-0-06784_A", "0-20000-0-10393_0",
         "0-20000-0-10471_0", "0-20000-0-10865_0", "0-20000-0-10962_0", "0-20000-0-73009_A",
         "0-20008-0-CGR_A", "0-20008-0-POT_A", "0-20008-0-UEX_A", "0-20008-0-UGR_A", "0-20008-0-UVA_A",
         "0-250-1001-07151_B", "0-250-1001-75105006_A", "0-276-0-10508_A", "0-380-5-1_0", "0-380-6-1_A",
         "0-620-3704-5480_B", "0-724-60000-0000000CIEMAT001_A", "0-20000-0-06610_A"]
PAYERNE = "0-20000-0-06610_A"
# Sites where marine aerosol dominates the boundary layer (coast < ~5 km, onshore flow most of the year).
MARINE = {"0-20000-0-03808_A", "0-20000-0-03838_A", "0-20008-0-CGR_A"}


def split_key(key):
    wmo, ident = key.rsplit("_", 1)
    return wmo, ident


def doy_of(t):
    return np.array([(EPOCH + timedelta(days=float(x))).timetuple().tm_yday for x in np.atleast_1d(t)], float)


# ---------------------------------------------------------------------------------------------
# Re-evaluation of the constant (mirror of calibration.rayleigh.calibration._compute_cl_for_perturbation)
# ---------------------------------------------------------------------------------------------
def cl_variant(rcs_mean, rng, bmol, fr, lr, shift=0.0, mode="nominal", z_fill=None, diag=False):
    """C_L for one (LR, window shift) with the operational Klett + Wiegner-Geiss chain.

    mode: 'nominal' (operational), 'noaer' (aerosol extinction set to 0 -> no transmission correction),
          'fill' (aerosol extinction below z_fill replaced by its mean over z_fill..z_fill+150 m).
    Returns C (and a diagnostics dict when diag=True).
    """
    from calibration.rayleigh.atmosphere import MOLECULAR_LIDAR_RATIO, klett_inversion
    from calibration.rayleigh.rayleigh_fit import calculate_lidar_constant

    s0, s1 = fr.range_start_m + shift, fr.range_end_m + shift
    if s0 < rng[0] or s1 > rng[-1]:
        return (np.nan, None) if diag else np.nan
    mol_mask = (rng >= s0) & (rng <= s1) & ~np.isnan(rcs_mean)
    idx = np.where(mol_mask)[0]
    if len(idx) < 3:
        return (np.nan, None) if diag else np.nan
    beta_att = rcs_mean / fr.slope                       # = (rcs/r^2)/slope * r^2 (subtract_background = False)
    ref = int((idx[0] + idx[-1]) / 2)
    ref_val = np.nanmean(beta_att[mol_mask] / bmol[mol_mask])
    beta_aer, beta_tot, ext_aer = klett_inversion(beta_att=beta_att, beta_mol=bmol, range_alc=rng,
                                                  reference_index=ref, lidar_ratio_aerosol=lr,
                                                  reference_value=ref_val, i_start=0, i_end=idx[-1])
    if mode == "noaer":
        ext_aer = np.zeros_like(ext_aer)
    elif mode == "fill":
        k = (rng >= z_fill) & (rng < z_fill + 150.0)
        ext_aer = ext_aer.copy()
        ext_aer[rng < z_fill] = np.nanmean(ext_aer[k])
    ext_tot = ext_aer + bmol * MOLECULAR_LIDAR_RATIO
    try:
        res = calculate_lidar_constant(rcs_mean=rcs_mean, beta_tot=beta_tot, ext_tot=ext_tot, range_alc=rng,
                                       molecular_mask=mol_mask, fit_result=fr, subtract_background=False)
        c = float(res.lidar_constant)
    except ValueError:
        c = np.nan
    if not (np.isfinite(c) and c > 0):
        c = np.nan
    if not diag:
        return c
    i0 = idx[0]
    tau_layers = []
    for lo, hi in TAU_LAYERS:
        m = (rng >= lo) & (rng < min(hi, rng[i0]))
        tau_layers.append(float(np.trapezoid(ext_aer[m], rng[m])) if m.sum() >= 2 else 0.0)
    d = dict(tau_aer_below=float(np.trapezoid(ext_aer[:i0], rng[:i0])),
             tau_mol_below=float(np.trapezoid(bmol[:i0] * MOLECULAR_LIDAR_RATIO, rng[:i0])),
             tau_layers=np.array(tau_layers),
             ext_prof=bin_profile(rng, ext_aer), R_prof=bin_profile(rng, beta_tot / bmol))
    return c, d


def bin_profile(rng, x):
    out = np.full(len(PROF_EDGES) - 1, np.nan)
    k = np.digitize(rng, PROF_EDGES) - 1
    for j in range(len(out)):
        m = (k == j) & np.isfinite(x)
        if m.any():
            out[j] = np.mean(x[m])
    return out


# ---------------------------------------------------------------------------------------------
# Extraction (one worker per unit)
# ---------------------------------------------------------------------------------------------
def census_info(key):
    wmo, ident = split_key(key)
    for s in json.load(open(CENSUS, encoding="utf-8")):
        if s["wmo"] == wmo and s["ident"] == ident:
            return s
    return None


def extract_unit(key):
    import calibration.rayleigh.calibration as RC
    from calibration.config import CalibrationOptions, DataLevel, InstrumentInfo, InstrumentType

    logging.disable(logging.CRITICAL)
    out_f = NIGHT_DIR / f"{key}.npz"
    if out_f.exists():
        return key, "exists"
    wmo, ident = split_key(key)
    cons = load_constants(wmo, ident)
    s = census_info(key)
    if cons is None or s is None:
        return key, "no constants/census"

    cap = {"samples": []}
    orig = RC._compute_cl_for_perturbation

    def wrap(*a, **k):
        if k.get("return_diagnostics"):
            cap["nominal"] = k
        elif k.get("lidar_ratio_aerosol") == LR_OP and k.get("altitude_shift_m") == 0:
            cap["samples"].append(k["rcs_mean"])          # one per time sample (all profiles + 4 subsets)
        return orig(*a, **k)

    RC._compute_cl_for_perturbation = wrap
    o = CalibrationOptions.from_json(REPO / "options.json")
    o.folder_root = Path("D:/E-PROFILE_L1_2026")
    o.data_level = DataLevel.L1
    o.apply_wv_correction = False
    o.folder_output = Path(tempfile.mkdtemp(prefix="q16d_"))
    o.plot_main = o.plot_all = False
    info = InstrumentInfo(site_name=s["site"], wmo_id=wmo, identifier=ident, instrument_type=InstrumentType.CHM15k,
                          latitude=float(s["lat"]), longitude=float(s["lon"]), altitude=float(s["alt"]))
    rows = []
    nl, nc, nf = len(LR_GRID), len(LR_CENTRES), len(NEAR_FILL_M)
    for i in range(len(cons["t"])):
        if not (np.isfinite(cons["c"][i]) and cons["c"][i] > 0):
            continue
        ds = (EPOCH + timedelta(days=float(cons["t"][i]) + 0.01)).strftime("%Y%m%d")
        cap.clear()
        cap["samples"] = []
        try:
            with contextlib.redirect_stdout(io.StringIO()):
                r = RC.calibrate_rayleigh(ds, info, o)
        except Exception:
            continue
        if "nominal" not in cap or not (np.isfinite(r.lidar_constant) and r.lidar_constant > 0):
            continue
        k = cap["nominal"]
        rcs, rng, bmol, fr = k["rcs_mean"], np.asarray(k["range_alc"], float), k["beta_mol"], k["fit_result"]
        rec = dict(t=cons["t"][i], c_op=cons["c"][i], c_re=float(r.lidar_constant), flag=float(r.flag),
                   module=str(cons["module"][i]), zb=fr.range_start_m, zt=fr.range_end_m)
        lnc_nom = np.full(nl, np.nan)
        dg = None
        for j, lr in enumerate(LR_GRID):
            if lr == LR_OP:
                c, dg = cl_variant(rcs, rng, bmol, fr, lr, diag=True)
            else:
                c = cl_variant(rcs, rng, bmol, fr, lr)
            lnc_nom[j] = np.log(c)
        lnc_noaer = np.log(cl_variant(rcs, rng, bmol, fr, LR_OP, mode="noaer"))
        lnc_fill = np.array([np.log(cl_variant(rcs, rng, bmol, fr, LR_OP, mode="fill", z_fill=zf))
                             for zf in NEAR_FILL_M])
        # Full operational ensemble (5 LR x 5 shifts x 5 time samples) re-centred on each LR centre
        samples = cap["samples"] if cap["samples"] else [rcs]
        lnc_ens = np.full(nc, np.nan)
        for j, lc in enumerate(LR_CENTRES):
            vals = []
            for smp in samples:
                for dlr in LR_DELTAS:
                    if lc + dlr <= 0:
                        continue
                    for sh in ALT_SHIFTS:
                        vals.append(cl_variant(smp, rng, bmol, fr, lc + dlr, shift=sh))
            vals = np.array(vals)
            vals = vals[np.isfinite(vals)]
            if vals.size:
                lnc_ens[j] = np.log(np.median(vals))
        rec.update(lnc_nom=lnc_nom, lnc_noaer=lnc_noaer, lnc_fill=lnc_fill, lnc_ens=lnc_ens,
                   tau_below=dg["tau_aer_below"] if dg else np.nan, tau_mol=dg["tau_mol_below"] if dg else np.nan,
                   tau_layers=dg["tau_layers"] if dg else np.full(len(TAU_LAYERS), np.nan),
                   ext_prof=dg["ext_prof"] if dg else np.full(len(PROF_EDGES) - 1, np.nan),
                   R_prof=dg["R_prof"] if dg else np.full(len(PROF_EDGES) - 1, np.nan),
                   n_samples=len(samples))
        rows.append(rec)
    if not rows:
        return key, "no nights"
    arr = {kk: np.array([r[kk] for r in rows]) for kk in rows[0]}
    np.savez(out_f, **arr)
    return key, f"{len(rows)} nights"


def run_extract():
    NIGHT_DIR.mkdir(parents=True, exist_ok=True)
    with Pool(N_WORKERS) as pool:
        for key, msg in pool.imap_unordered(extract_unit, UNITS):
            print(f"{key}: {msg}", flush=True)


# ---------------------------------------------------------------------------------------------
# Analysis helpers
# ---------------------------------------------------------------------------------------------
def harmonic(t, y, trend=True):
    """Annual harmonic (+ linear trend) by least squares: (amp, peak DOY, in-phase amp on DOY 176, cos, sin)."""
    ok = np.isfinite(t) & np.isfinite(y)
    if ok.sum() < 30 or (t[ok].max() - t[ok].min()) < 250:
        return np.nan, np.nan, np.nan
    doy = doy_of(t[ok])
    cols = [np.ones(ok.sum()), np.cos(W_YEAR * doy), np.sin(W_YEAR * doy)]
    if trend:
        cols.append((t[ok] - t[ok].mean()) / 365.25)
    coef, *_ = np.linalg.lstsq(np.column_stack(cols), y[ok], rcond=None)
    a, b = coef[1], coef[2]
    amp = float(np.hypot(a, b))
    peak = float((np.arctan2(b, a) / W_YEAR) % 365.25)
    inph = float(a * np.cos(W_YEAR * 176) + b * np.sin(W_YEAR * 176))
    return amp, peak, inph


def running_median(t, x, window=31.0):
    out = np.full_like(x, np.nan, dtype=float)
    lo = np.searchsorted(t, t - window / 2)
    hi = np.searchsorted(t, t + window / 2, side="right")
    for i in range(len(t)):
        seg = x[lo[i]:hi[i]]
        seg = seg[np.isfinite(seg)]
        if seg.size >= 5:
            out[i] = np.median(seg)
    return out


def anomaly_std(t, y):
    a = y - running_median(t, y)
    a = a[np.isfinite(a)]
    return float(np.std(a)) if a.size > 10 else np.nan


def interp_lnc(lnc_grid, lr):
    """ln C at an arbitrary LR (per night, linear in LR on the stored grid)."""
    return np.array([np.interp(l, LR_GRID, g) if np.all(np.isfinite(g)) else np.nan for g, l in zip(lnc_grid, lr)])


def read_aeronet_daily(path):
    """{date: aod1064} from an AERONET v3 daily-average file (1020 nm -> 1064 nm with the 440-870 Angstrom)."""
    out = {}
    lines = open(path, encoding="utf-8", errors="replace").read().splitlines()
    hdr = next((i for i, l in enumerate(lines) if l.startswith("AERONET_Site,")), None)
    if hdr is None:
        return out
    cols = lines[hdr].split(",")
    ix = {c: k for k, c in enumerate(cols)}
    for l in lines[hdr + 1:]:
        p = l.split(",")
        if len(p) < len(cols) - 2:
            continue
        try:
            d = datetime.strptime(p[ix["Date(dd:mm:yyyy)"]], "%d:%m:%Y").date()
            a1020 = float(p[ix["AOD_1020nm"]])
            ang = float(p[ix["440-870_Angstrom_Exponent"]])
        except (ValueError, KeyError):
            continue
        if a1020 < -900 or ang < -900:
            continue
        out[d] = a1020 * (1020.0 / 1064.0) ** ang
    return out


def aeronet_for(key):
    match = json.load(open(SCRATCH / "aeronet_match.json", encoding="utf-8"))
    for k, site, ae_site, dist in match:
        if k == key:
            p = AERONET_DIR / f"AOD15_daily_{ae_site}.txt"
            return read_aeronet_daily(p) if p.exists() else {}, ae_site
    return {}, ""


def night_aod(aer, t):
    """Mean of the daily AODs of the evening (D-1) and morning (D) of the night ending on day t."""
    d1 = (EPOCH + timedelta(days=float(t) + 0.01)).date()
    d0 = d1 - timedelta(days=1)
    v = [aer[d] for d in (d0, d1) if d in aer]
    return float(np.mean(v)) if v else np.nan


def season_of(doy):
    m = np.array([(datetime(2025, 1, 1) + timedelta(days=int(d) - 1)).month for d in doy])
    return np.select([np.isin(m, (12, 1, 2)), np.isin(m, (3, 4, 5)), np.isin(m, (6, 7, 8))], ["DJF", "MAM", "JJA"], "SON")


# ---------------------------------------------------------------------------------------------
# Analysis
# ---------------------------------------------------------------------------------------------
def load_unit(key):
    f = NIGHT_DIR / f"{key}.npz"
    if not f.exists():
        return None
    with np.load(f, allow_pickle=False) as z:
        U = {k: z[k] for k in z.files}
    # main optical module only (era split) and a one-pass 4-MAD screen on ln C_op within it
    mods, cnt = np.unique(U["module"], return_counts=True)
    keep = U["module"] == mods[np.argmax(cnt)]
    lnc = np.log(U["c_op"])
    med = np.nanmedian(lnc[keep])
    mad = 1.4826 * np.nanmedian(np.abs(lnc[keep] - med))
    keep &= np.abs(lnc - med) <= 4 * mad
    order = np.argsort(U["t"])
    order = order[keep[order]]
    return {k: v[order] for k, v in U.items()}


def scenario_bias(U, lr_true):
    """ln(C_op / C_true) when the real LR is lr_true but the chain uses 52 sr (nominal configuration)."""
    i52 = int(np.where(LR_GRID == LR_OP)[0][0])
    return U["lnc_nom"][:, i52] - interp_lnc(U["lnc_nom"], lr_true)


def analyse():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    i52 = int(np.where(LR_GRID == LR_OP)[0][0])
    j52 = int(np.where(LR_CENTRES == LR_OP)[0][0])
    unit_rows, night_rows = [], []
    pooled = {s: {"ext": [], "R": [], "aod": [], "tau": []} for s in ("DJF", "MAM", "JJA", "SON")}
    pooled_anom = []
    for key in UNITS:
        U = load_unit(key)
        if U is None or len(U["t"]) < 30:
            continue
        t = U["t"]
        doy = doy_of(t)
        seas = season_of(doy)
        lnc_op = np.log(U["c_op"])
        repro = np.nanmedian(np.abs(np.log(U["c_re"] / U["c_op"])))
        L = U["lnc_nom"]
        E = U["lnc_ens"]
        two_tau_eff = L[:, i52] - U["lnc_noaer"]                       # total ln-correction applied by Klett (LR 52)
        # dlnC/dLR around 52 (nominal and full ensemble), in % per 10 sr
        dnom = 100 * (L[:, LR_GRID == 60][:, 0] - L[:, LR_GRID == 40][:, 0]) / 20 * 10
        dens = 100 * (E[:, LR_CENTRES == 60][:, 0] - E[:, LR_CENTRES == 40][:, 0]) / 20 * 10
        ens_vs_nom = 100 * (E[:, j52] - L[:, i52])
        # scenarios ----------------------------------------------------------------------------
        lr_s1 = 52 + 18 * np.cos(W_YEAR * (doy - 182))                   # summer smoke/pollution: 34 -> 70 sr
        lr_s3 = np.full_like(doy, 25.0)                                  # marine
        b_s1 = scenario_bias(U, lr_s1)
        b_s3 = scenario_bias(U, lr_s3)
        b_zero = L[:, i52] - U["lnc_noaer"]                              # bias if the true extinction were 0
        rng_ = np.random.default_rng(0)
        b_rand = scenario_bias(U, np.clip(52 + 15 * rng_.standard_normal(len(doy)), 10, 100))
        fill = {zf: L[:, i52] - U["lnc_fill"][:, k] for k, zf in enumerate(NEAR_FILL_M)}  # op minus near-fill
        # AERONET -----------------------------------------------------------------------------
        aer, ae_site = aeronet_for(key)
        aod = np.array([night_aod(aer, x) for x in t]) if aer else np.full(len(t), np.nan)
        tau_k = U["tau_below"]
        okA = np.isfinite(aod) & np.isfinite(tau_k)
        ratio_tk_aod = np.nanmedian(tau_k[okA] / aod[okA]) if okA.sum() > 10 else np.nan
        r_tk_aod = np.corrcoef(tau_k[okA], aod[okA])[0, 1] if okA.sum() > 10 else np.nan
        miss = {}
        for f in (0.6, 0.75, 0.9):
            miss[f] = np.where(okA, -2 * (f * aod - tau_k), np.nan)      # ln(C_op/C_true) from missing tau
        row = dict(stream=key, n=len(t), aeronet=ae_site, n_aod=int(okA.sum()), repro_pct=100 * repro,
                   tau_k_med=np.nanmedian(tau_k), two_tau_eff_med=100 * np.nanmedian(two_tau_eff),
                   frac_tau_lt200=np.nanmedian(U["tau_layers"][:, 0] / tau_k),
                   frac_tau_lt500=np.nanmedian((U["tau_layers"][:, 0] + U["tau_layers"][:, 1]) / tau_k),
                   frac_tau_lt1000=np.nanmedian(U["tau_layers"][:, :3].sum(1) / tau_k),
                   dlnC_per10sr_nom=np.nanmedian(dnom), dlnC_per10sr_ens=np.nanmedian(dens),
                   ens_minus_nom=np.nanmedian(ens_vs_nom), aod_med=np.nanmedian(aod), tauk_over_aod=ratio_tk_aod,
                   r_tauk_aod=r_tk_aod, marine=key in MARINE)
        for name, y in (("obs", lnc_op), ("twotau", two_tau_eff), ("s1", b_s1), ("s3", b_s3), ("zero", b_zero),
                        ("rand", b_rand), ("fill300", fill[300.0]), ("fill600", fill[600.0]),
                        ("fill1000", fill[1000.0]), ("miss75", miss[0.75]), ("aod", aod), ("tauk", tau_k)):
            amp, peak, inph = harmonic(t, y)
            row[f"amp_{name}"] = 100 * amp if name not in ("aod", "tauk") else amp
            row[f"peak_{name}"] = peak
            row[f"inph_{name}"] = 100 * inph if name not in ("aod", "tauk") else inph
            row[f"n2n_{name}"] = 100 * anomaly_std(t, y) if name not in ("aod", "tauk") else anomaly_std(t, y)
            row[f"mean_{name}"] = 100 * np.nanmean(y) if name not in ("aod", "tauk") else np.nanmean(y)
        # seasonal means of the missing term and of tau_K / AOD
        for s in ("DJF", "MAM", "JJA", "SON"):
            m = seas == s
            row[f"tauk_{s}"] = np.nanmedian(tau_k[m]) if m.any() else np.nan
            row[f"aod_{s}"] = np.nanmedian(aod[m & okA]) if (m & okA).any() else np.nan
            row[f"twotau_{s}"] = 100 * np.nanmedian(two_tau_eff[m]) if m.any() else np.nan
            row[f"miss75_{s}"] = 100 * np.nanmedian(miss[0.75][m & okA]) if (m & okA).sum() >= 5 else np.nan
            row[f"fill600_{s}"] = 100 * np.nanmedian(fill[600.0][m]) if m.any() else np.nan
            if key != PAYERNE and m.any():
                pooled[s]["ext"].append(np.nanmedian(U["ext_prof"][m], axis=0))
                pooled[s]["R"].append(np.nanmedian(U["R_prof"][m], axis=0))
        # data-driven LR: which LR minimises the night-to-night scatter / the in-phase amplitude of ln C?
        lr_scan = np.r_[0.0, LR_GRID]
        lnc_scan = np.column_stack([U["lnc_noaer"], L])
        n2n_scan = np.array([anomaly_std(t, lnc_scan[:, j]) for j in range(len(lr_scan))])
        inph_scan = np.array([harmonic(t, lnc_scan[:, j])[2] for j in range(len(lr_scan))])
        row["lr_min_n2n"] = lr_scan[int(np.nanargmin(n2n_scan))] if np.any(np.isfinite(n2n_scan)) else np.nan
        row["n2n_lr0"], row["n2n_lr52"], row["n2n_lr100"] = 100 * n2n_scan[0], 100 * n2n_scan[5], 100 * n2n_scan[-1]
        row["inph_lr0"], row["inph_lr52"], row["inph_lr100"] = 100 * inph_scan[0], 100 * inph_scan[5], 100 * inph_scan[-1]
        # anomaly regression: ln C_op anomaly on 2 tau_K anomaly. If tau_true = k tau_K, slope = 1 - k.
        a_y = lnc_op - running_median(t, lnc_op)
        a_x = two_tau_eff - running_median(t, two_tau_eff)
        ok = np.isfinite(a_x) & np.isfinite(a_y)
        if ok.sum() > 20:
            sl = np.polyfit(a_x[ok], a_y[ok], 1)[0]
            row["anom_slope"] = sl
            row["anom_r"] = np.corrcoef(a_x[ok], a_y[ok])[0, 1]
            pooled_anom.append(np.column_stack([a_x[ok], a_y[ok]]) if key != PAYERNE else np.zeros((0, 2)))
        else:
            row["anom_slope"] = row["anom_r"] = np.nan
        unit_rows.append(row)
        for i in range(len(t)):
            night_rows.append(dict(stream=key, t=t[i], doy=doy[i], lnc_op=lnc_op[i], lnc_re=np.log(U["c_re"][i]),
                                   zb=U["zb"][i], zt=U["zt"][i], tau_k=tau_k[i], two_tau_eff=two_tau_eff[i],
                                   **{f"lnc_lr{int(l)}": L[i, j] for j, l in enumerate(LR_GRID)},
                                   **{f"lnc_ens{int(l)}": E[i, j] for j, l in enumerate(LR_CENTRES)},
                                   lnc_noaer=U["lnc_noaer"][i],
                                   **{f"lnc_fill{int(z)}": U["lnc_fill"][i, k] for k, z in enumerate(NEAR_FILL_M)},
                                   **{f"tau_{lo}_{int(min(hi, 9999))}": U["tau_layers"][i, k]
                                      for k, (lo, hi) in enumerate(TAU_LAYERS)},
                                   aod1064=aod[i]))
        print(f"{key}: n={len(t)} repro={100*repro:.3f}% 2tau={row['two_tau_eff_med']:.2f}% "
              f"dC/10sr={row['dlnC_per10sr_ens']:.2f}% amp_s1={row['amp_s1']:.2f} peak={row['peak_s1']:.0f}", flush=True)

    # ---- write tables ------------------------------------------------------------------------
    with open(Q16 / "q16d_units.csv", "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(unit_rows[0].keys()))
        w.writeheader()
        w.writerows(unit_rows)
    with open(Q16 / "q16d_nights.csv", "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(night_rows[0].keys()))
        w.writeheader()
        w.writerows(night_rows)

    net = [r for r in unit_rows if r["stream"] != PAYERNE]

    def med(col, rows=net):
        v = np.array([r[col] for r in rows], float)
        v = v[np.isfinite(v)]
        if not v.size:
            return "nan"
        return f"{np.median(v):.2f} [{np.percentile(v, 25):.2f}, {np.percentile(v, 75):.2f}]"

    summary_cols = ["n", "repro_pct", "tau_k_med", "two_tau_eff_med", "frac_tau_lt200", "frac_tau_lt500",
                    "frac_tau_lt1000", "lr_min_n2n", "n2n_lr0", "n2n_lr52", "n2n_lr100", "inph_lr0", "inph_lr52",
                    "inph_lr100", "anom_slope", "anom_r", "dlnC_per10sr_nom", "dlnC_per10sr_ens", "ens_minus_nom", "aod_med",
                    "tauk_over_aod", "r_tauk_aod"]
    for name in ("obs", "twotau", "s1", "s3", "zero", "rand", "fill300", "fill600", "fill1000", "miss75"):
        summary_cols += [f"amp_{name}", f"peak_{name}", f"inph_{name}", f"n2n_{name}", f"mean_{name}"]
    for s in ("DJF", "MAM", "JJA", "SON"):
        summary_cols += [f"tauk_{s}", f"aod_{s}", f"twotau_{s}", f"miss75_{s}", f"fill600_{s}"]
    with open(Q16 / "q16d_network.txt", "w", encoding="utf-8") as f:
        f.write(f"network units (Payerne excluded): {len(net)}\n")
        for c in summary_cols:
            f.write(f"{c:22s} {med(c)}\n")
        PA = np.vstack(pooled_anom)
        sl, ic = np.polyfit(PA[:, 0], PA[:, 1], 1)
        f.write(f"pooled anomaly slope dlnC_op/d(2tau_K) = {sl:.3f} (r = {np.corrcoef(PA.T)[0, 1]:.3f}, n = {len(PA)})"
                f" -> implied tau_true/tau_K = {1 - sl:.2f}, LR_eff ~ {52 * (1 - sl):.0f} sr" + "\n")
        rngb = np.random.default_rng(1)
        bs = [np.polyfit(*PA[rngb.integers(0, len(PA), len(PA))].T, 1)[0] for _ in range(300)]
        f.write(f"bootstrap 95% CI of the slope: {np.percentile(bs, 2.5):.3f} .. {np.percentile(bs, 97.5):.3f}" + "\n")
        pay = [r for r in unit_rows if r["stream"] == PAYERNE]
        if pay:
            f.write("\nPayerne A:\n")
            for c in summary_cols:
                f.write(f"{c:22s} {pay[0][c]}\n")
    print(open(Q16 / "q16d_network.txt", encoding="utf-8").read())

    # ---- figures ------------------------------------------------------------------------------
    zmid = 0.5 * (PROF_EDGES[:-1] + PROF_EDGES[1:]) / 1000
    cols = {"DJF": "#3b6fb6", "MAM": "#57a773", "JJA": "#d1495b", "SON": "#e0a030"}
    fig, ax = plt.subplots(1, 3, figsize=(16, 5.5))
    for s, c in cols.items():
        if pooled[s]["ext"]:
            e = np.nanmedian(np.array(pooled[s]["ext"]), axis=0) * 1e3
            r = np.nanmedian(np.array(pooled[s]["R"]), axis=0)
            ax[0].plot(e, zmid, color=c, label=s, lw=2)
            ax[1].plot(r, zmid, color=c, label=s, lw=2)
    ax[0].set_xlabel("Klett aerosol extinction at LR 52 sr [km$^{-1}$]")
    ax[0].set_ylabel("range above the instrument [km]")
    ax[0].set_title("Network median Klett extinction (calibration nights)")
    ax[0].axhspan(0, 0.8, color="0.85", alpha=0.5, zorder=0)
    ax[0].text(0.98, 0.4, "incomplete overlap\n(firmware-corrected)", transform=ax[0].get_yaxis_transform(),
               ha="right", fontsize=8, color="0.3")
    ax[0].legend()
    ax[1].set_xlabel("backscatter ratio R = beta_tot / beta_mol")
    ax[1].set_xscale("log")
    ax[1].set_title("Retrieved backscatter ratio")
    ax[1].axhspan(0, 0.8, color="0.85", alpha=0.5, zorder=0)
    # right: dlnC vs LR, network median curve per season (relative to 52)
    for s, c in cols.items():
        curves = []
        for key in UNITS:
            if key == PAYERNE:
                continue
            U = load_unit(key)
            if U is None:
                continue
            m = season_of(doy_of(U["t"])) == s
            if m.sum() < 5:
                continue
            curves.append(np.nanmedian(100 * (U["lnc_nom"][m] - U["lnc_nom"][m][:, [i52]]), axis=0))
        if curves:
            ax[2].plot(LR_GRID, np.nanmedian(np.array(curves), axis=0), "o-", color=c, label=s)
    ax[2].axvline(52, color="k", ls=":")
    ax[2].axhline(0, color="k", lw=0.5)
    ax[2].set_xlabel("aerosol lidar ratio used in the Klett inversion [sr]")
    ax[2].set_ylabel("ln C(LR) - ln C(52) [%]")
    ax[2].set_title("Sensitivity of C_L to the lidar ratio (network median)")
    ax[2].legend()
    fig.tight_layout()
    fig.savefig(Q16 / "q16d_profiles_sensitivity.png", dpi=130)
    plt.close(fig)

    # monthly composites: observed ln C anomaly vs predicted biases
    fig, ax = plt.subplots(1, 2, figsize=(16, 5.5))
    month_terms = {"observed ln C_op (raw, kappa not removed)": "obs",
                   "2 tau_K applied (LR 52) = bias if true extinction = 0": "twotau",
                   "bias if true LR = 34->70 sr (S1)": "s1",
                   "near-range fill 600 m minus op": "fill600", "-2 x missing tau (0.75 AOD - tau_K)": "miss75"}
    acc = {k: [[] for _ in range(12)] for k in month_terms}
    for key in UNITS:
        if key == PAYERNE:
            continue
        U = load_unit(key)
        if U is None:
            continue
        t = U["t"]
        mon = np.array([(EPOCH + timedelta(days=float(x))).month for x in t])
        doy = doy_of(t)
        aer, _ = aeronet_for(key)
        aod = np.array([night_aod(aer, x) for x in t]) if aer else np.full(len(t), np.nan)
        series = {"obs": np.log(U["c_op"]), "twotau": U["lnc_nom"][:, i52] - U["lnc_noaer"],
                  "s1": scenario_bias(U, 52 + 18 * np.cos(W_YEAR * (doy - 182))),
                  "zero": U["lnc_nom"][:, i52] - U["lnc_noaer"],
                  "fill600": U["lnc_nom"][:, i52] - U["lnc_fill"][:, 1],
                  "miss75": -2 * (0.75 * aod - U["tau_below"])}
        for lab, k in month_terms.items():
            y = series[k]
            ym = np.nanmedian(y)
            for mm in range(1, 13):
                sel = (mon == mm) & np.isfinite(y)
                if sel.sum() >= 3:
                    acc[lab][mm - 1].append(100 * (np.nanmedian(y[sel]) - ym))
    for lab, v in acc.items():
        ax[0].plot(range(1, 13), [np.median(x) if x else np.nan for x in v], "o-", label=lab)
    ax[0].axhline(0, color="k", lw=0.5)
    ax[0].set_xlabel("month")
    ax[0].set_ylabel("monthly median minus unit median [%]")
    ax[0].set_title("Observed C_L cycle vs aerosol-correction terms (network median)")
    ax[0].legend(fontsize=8)
    # scatter tau_K vs AERONET AOD
    nn = [r for r in night_rows if r["stream"] != PAYERNE and np.isfinite(r["aod1064"]) and np.isfinite(r["tau_k"])]
    x = np.array([r["aod1064"] for r in nn])
    y = np.array([r["tau_k"] for r in nn])
    sc = np.array([season_of([r["doy"]])[0] for r in nn])
    for s, c in cols.items():
        m = sc == s
        ax[1].scatter(x[m], y[m], s=5, color=c, alpha=0.4, label=f"{s} (median tau_K/AOD {np.median(y[m] / x[m]):.2f})")
    lim = [0, 0.3]
    ax[1].plot(lim, lim, "k-", lw=0.8)
    ax[1].plot(lim, [0, 0.225], "k--", lw=0.8, label="0.75 x AOD")
    ax[1].set_xlim(lim)
    ax[1].set_ylim(lim)
    ax[1].set_xlabel("AERONET AOD at 1064 nm (daily mean, evening + morning)")
    ax[1].set_ylabel("Klett tau below the window (LR 52)")
    ax[1].set_title("Klett aerosol optical depth vs AERONET")
    ax[1].legend(fontsize=8, markerscale=3)
    fig.tight_layout()
    fig.savefig(Q16 / "q16d_monthly_aeronet.png", dpi=130)
    plt.close(fig)


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "extract":
        run_extract()
    elif len(sys.argv) > 1 and sys.argv[1] == "analyse":
        analyse()
    else:
        print(__doc__)
