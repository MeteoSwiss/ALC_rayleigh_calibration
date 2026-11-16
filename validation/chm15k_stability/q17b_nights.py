# -*- coding: utf-8 -*-
"""q17b (nights): per-night clear-sky molecular-normalised signals X(z) of the CHM15k and the CL61.

X(z) = rcs_0 / (beta_mol * T_mol^2)            (CHM15k, 1064 nm)
X(z) = rcs_0 / (beta_mol * T_mol^2 * T_wv^2)   (CL61, 910.55 nm, FWHM 0.188 nm)

* Same 5-min blocks for both instruments: night (solar zenith > 100 deg, the operational criterion),
  both instruments recording, NO cloud base (any layer) reported by EITHER instrument in the block
  or its two neighbours (+-5 min); >= 6 such blocks per night.
* Molecular profile from CAMS 0.4 deg L137 T/p (time-averaged over the night's clear span +-1.5 h),
  Bucholtz backscatter via calibration.rayleigh.atmosphere.calculate_molecular_properties - NOT US-1976.
  Below the lowest CAMS level (valley sites) p is extrapolated hydrostatically, T with -6.5 K/km.
* CL61 two-way water-vapour transmission from the repo's operational kernel
  (calibration.water_vapor_correction.water_vapor.two_way_wv_transmission, HITRAN LUT shipped with the
  package) with the OPERATIONAL CL61 spectrum 910.55 nm / FWHM 0.188 nm and CAMS humidity. Two fills
  below the CAMS surface: 'op' = constant number density (what the operational code does) and
  'q' = constant specific humidity (moister valley air) -> bounds the orography error.
* Measured dark baseline (rayleigh_availability/dark_profiles.py, Payerne hood campaign): subtracted
  for Payerne A (CHM15k) and Payerne C (CL61); the Payerne CL61 template is also applied to the other
  CL61 as a sensitivity variant ('dark'), since the remote-dark campaign closed the CL61 template at
  theta ~ 1.
* AERONET (L1.5 daily, co-located) 31-day running AOD at 1064 / 910 nm (from 1020 / 870 nm + AE
  440-870) for the optional below-band aerosol two-way transmission correction.
* Also a CHM15k-only night table (CHM15k clear blocks, whole 2025-26 record) for the full-year cycle.

Outputs: q17/q17b_nights.csv (common nights), q17/q17b_nights_chm.csv (CHM15k-only nights)
Run:     .venv/Scripts/python.exe validation/chm15k_stability/q17b_nights.py [--workers 4]
"""
from __future__ import annotations

import argparse
import csv
import sys
import warnings
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1]))
sys.path.insert(0, str(HERE))
from q17b_extract import BLK, NBLK, Q17, SITES, Z_EDGES  # noqa: E402

warnings.filterwarnings("ignore")
BANDS = ((2000, 3000), (3000, 4000), (4000, 5000), (5000, 6000), (6000, 8000), (8000, 10000))
BAND_NAMES = tuple(f"{a // 1000}-{b // 1000}" for a, b in BANDS)
WL_CHM, WL_CL61, FWHM_CL61 = 1064e-9, 910.55, 0.188
DARK_NPZ = Path("C:/DATA/Projects/202606_E-PROFILE_calibration/rayleigh_availability/dark_profiles_payerne.npz")
AERONET = Path("C:/Users/hervo/AppData/Local/Temp/claude/C--Users-hervo-OneDrive-Documents-ALC-rayleigh-calibration/"
               "e6120125-e8c9-4e0b-8352-26777f3067de/scratchpad/aeronet_raw")
AERONET_SITE = {"Aosta": "Aosta", "Lindenberg": "MetObs_Lindenberg", "Camborne": "Camborne_MO",
                "Payerne": "Payerne"}
MIN_BLOCKS = 6
SZA_NIGHT = 100.0
KB, EPS, G0, RD = 1.380649e-23, 0.621981, 9.80665, 287.058
ZF = np.arange(5.0, 12000.0, 10.0)        # fine vertical grid above the instrument [m]
ZMID = 0.5 * (Z_EDGES[:-1] + Z_EDGES[1:])


def sza_deg(t, lat, lon):
    """Solar zenith angle (deg) for UTC days since 1970 (NOAA low-precision, includes eq. of time)."""
    d = np.asarray(t, float) - 10957.5
    g = np.deg2rad(357.529 + 0.98560028 * d)
    q = 280.459 + 0.98564736 * d
    L = np.deg2rad(q + 1.915 * np.sin(g) + 0.020 * np.sin(2 * g))
    e = np.deg2rad(23.439 - 0.00000036 * d)
    ra = np.rad2deg(np.arctan2(np.cos(e) * np.sin(L), np.cos(L)))
    dec = np.arcsin(np.sin(e) * np.sin(L))
    gmst = np.mod(18.697374558 + 24.06570982441908 * d, 24.0)
    ha = np.deg2rad(gmst * 15.0 + lon - ra)
    la = np.deg2rad(lat)
    cz = np.sin(la) * np.sin(dec) + np.cos(la) * np.cos(dec) * np.cos(ha)
    return np.rad2deg(np.arccos(np.clip(cz, -1, 1)))


def load_blocks(site, instr):
    f = BLK / f"{site}_{instr}.npz"
    if not f.exists():
        return None
    with np.load(f, allow_pickle=False) as z:
        return {k: z[k] for k in z.files}


def timeline(b, d0, nd):
    """Scatter a per-day block file on a continuous (nd*288) timeline starting at day d0."""
    n = nd * NBLK
    rcs = np.full((n, len(ZMID)), np.nan, "f4")
    npf = np.zeros(n, "f4")
    ncl = np.zeros(n, "f4")
    hk = np.full((n, b["hk"].shape[-1]), np.nan, "f4")
    for k, day in enumerate(b["day"]):
        i = int(day - d0)
        if 0 <= i < nd:
            s = slice(i * NBLK, (i + 1) * NBLK)
            rcs[s] = b["rcs"][k]
            npf[s] = b["nprof"][k]
            ncl[s] = b["ncloud"][k]
            hk[s] = b["hk"][k]
    return rcs, npf, ncl, hk


def dilate(mask, k=1):
    out = mask.copy()
    for s in range(1, k + 1):
        out[s:] |= mask[:-s]
        out[:-s] |= mask[s:]
    return out


def binned_dark(key, tilt):
    """Hood-measured dark b(z) (rcs units) of Payerne unit `key` on the 100 m vertical bins."""
    with np.load(DARK_NPZ) as z:
        rng, b = z[f"{key}_range"], z[f"{key}_b_rcs"]
    zz = rng * np.cos(np.deg2rad(tilt))
    idx = np.digitize(zz, Z_EDGES) - 1
    out = np.full(len(ZMID), np.nan)
    for k in range(len(ZMID)):
        m = (idx == k) & np.isfinite(b)
        if m.any():
            out[k] = b[m].mean()
    return out


# ---------------------------------------------------------------------------------------------
# CAMS molecular + water vapour
# ---------------------------------------------------------------------------------------------
def load_cams(site):
    with np.load(Q17 / f"q17b_cams_{site}.npz") as z:
        return {k: z[k] for k in z.files}


def night_profiles(cams, t0, t1, alt_st):
    """CAMS T, p, n_wv (two fills) on alt_st + ZF averaged over the steps within [t0-1.5h, t1+1.5h]."""
    sel = np.where((cams["time"] >= t0 - 0.0625) & (cams["time"] <= t1 + 0.0625))[0]
    if sel.size == 0:
        sel = np.array([int(np.argmin(np.abs(cams["time"] - 0.5 * (t0 + t1))))])
    alt_f = alt_st + ZF
    acc = {k: np.zeros(len(ZF)) for k in ("T", "p", "nop", "nq")}
    pwv = 0.0
    for i in sel:
        za = cams["alt"][i].astype(float)
        o = np.argsort(za)
        za, T, p, q = za[o], cams["T"][i][o].astype(float), cams["p"][i][o].astype(float), cams["q"][i][o].astype(float)
        Tf = np.interp(alt_f, za, T)
        pf = np.exp(np.interp(alt_f, za, np.log(p)))
        qf = np.interp(alt_f, za, q)
        below = alt_f < za[0]
        if below.any():
            dz = za[0] - alt_f[below]
            Tf[below] = T[0] + 0.0065 * dz
            Tm = 0.5 * (T[0] + Tf[below])
            pf[below] = p[0] * np.exp(G0 * dz / (RD * Tm))
            qf[below] = q[0]
        e = qf * pf / (EPS + (1 - EPS) * qf)
        n_q = e / (KB * Tf)
        e0 = q[0] * p[0] / (EPS + (1 - EPS) * q[0])
        n_op = n_q.copy()
        n_op[below] = e0 / (KB * T[0])            # operational: constant number density below CAMS
        acc["T"] += Tf
        acc["p"] += pf
        acc["nop"] += n_op
        acc["nq"] += n_q
        # precipitable water above the station (op fill), kg/m2 = mm
        rho_v = n_op * 18.015e-3 / 6.02214e23
        pwv += np.trapezoid(rho_v, ZF) if hasattr(np, "trapezoid") else np.trapz(rho_v, ZF)
    n = len(sel)
    return {k: v / n for k, v in acc.items()}, pwv / n, alt_f[0] < cams["alt"][sel[0]].min()


def denominators(prof, alt_st, tilt_chm, tilt_cl):
    """Bin-mean beta_mol*T_mol^2 (CHM15k), and for the CL61 beta_mol*T_mol^2 and the WV T^2 (two fills)."""
    from calibration.rayleigh.atmosphere import calculate_molecular_properties
    from calibration.water_vapor_correction.water_vapor import two_way_wv_transmission
    idx = np.digitize(ZF, Z_EDGES) - 1

    def binmean(x):
        s = np.bincount(idx, weights=x, minlength=len(ZMID))[:len(ZMID)]
        c = np.bincount(idx, minlength=len(ZMID))[:len(ZMID)]
        return s / np.maximum(c, 1)

    mp_a = calculate_molecular_properties(prof["T"], prof["p"], ZF / np.cos(np.deg2rad(tilt_chm)), WL_CHM)
    mp_b = calculate_molecular_properties(prof["T"], prof["p"], ZF / np.cos(np.deg2rad(tilt_cl)), WL_CL61 * 1e-9)
    alt_f = alt_st + ZF
    t2op = two_way_wv_transmission(alt_f, alt_st, alt_f, prof["nop"], None, WL_CL61, FWHM_CL61)
    t2q = two_way_wv_transmission(alt_f, alt_st, alt_f, prof["nq"], None, WL_CL61, FWHM_CL61)
    return dict(chm=binmean(mp_a.beta_att_mol), cl=binmean(mp_b.beta_att_mol),
                cl_wv=binmean(mp_b.beta_att_mol * t2op), cl_wvq=binmean(mp_b.beta_att_mol * t2q),
                t2op=binmean(t2op), t2q=binmean(t2q))


def band_lnx(rcs_mean, den):
    out = []
    for a, b in BANDS:
        m = (ZMID >= a) & (ZMID < b) & np.isfinite(rcs_mean) & np.isfinite(den)
        v = np.mean(rcs_mean[m] / den[m]) if m.sum() >= 0.6 * ((b - a) / 100) else np.nan
        out.append(np.log(v) if np.isfinite(v) and v > 0 else np.nan)
    return np.array(out)


def band_sigma(rcs_blocks, den):
    """Relative noise of the night-mean band X from the block-to-block scatter (upper bound: includes
    real atmospheric variability within the night)."""
    out = []
    for a, b in BANDS:
        m = (ZMID >= a) & (ZMID < b) & np.isfinite(den)
        xb = np.nanmean(rcs_blocks[:, m] / den[m], axis=1)
        xb = xb[np.isfinite(xb)]
        out.append(np.std(xb, ddof=1) / abs(np.mean(xb)) / np.sqrt(len(xb)) if len(xb) > 2 and np.mean(xb) != 0 else np.nan)
    return np.array(out)


# ---------------------------------------------------------------------------------------------
# AERONET
# ---------------------------------------------------------------------------------------------
def aeronet_running(site):
    import pandas as pd
    f = AERONET / f"AOD15_daily_{AERONET_SITE[site]}.txt"
    if not f.exists():
        return None
    df = pd.read_csv(f, skiprows=5, na_values=[-999.0, "-999."], encoding="latin-1")
    dt = pd.to_datetime(df["Date(dd:mm:yyyy)"], format="%d:%m:%Y")
    t = ((dt - pd.Timestamp("1970-01-01")).dt.days + 0.5).to_numpy(float)
    ae = df["440-870_Angstrom_Exponent"].to_numpy(float)
    t1064 = df["AOD_1020nm"].to_numpy(float) * (1064.0 / 1020.0) ** (-ae)
    t910 = df["AOD_870nm"].to_numpy(float) * (910.55 / 870.0) ** (-ae)
    ok = np.isfinite(t1064) & np.isfinite(t910)
    return t[ok], t1064[ok], t910[ok], ae[ok]


def running_at(tq, t, y, half=15.0, min_n=5):
    out = np.full(len(tq), np.nan)
    for k, x in enumerate(tq):
        m = np.abs(t - x) <= half
        if m.sum() >= min_n:
            out[k] = np.median(y[m])
    # fill gaps by linear interpolation in time
    ok = np.isfinite(out)
    if ok.sum() >= 2:
        out[~ok] = np.interp(tq[~ok], tq[ok], out[ok])
    return out


# ---------------------------------------------------------------------------------------------
# Per-site processing
# ---------------------------------------------------------------------------------------------
def process_site(site):
    wmo, ich, icl, lat, lon = SITES[site]
    A = load_blocks(site, "CHM15k")
    B = load_blocks(site, "CL61")
    cams = load_cams(site)
    d0 = int(min(A["day"].min(), B["day"].min()))
    nd = int(max(A["day"].max(), B["day"].max()) - d0 + 2)
    ra, npa, cla, hka = timeline(A, d0, nd)
    rb, npb, clb, hkb = timeline(B, d0, nd)
    tb = d0 + (np.arange(nd * NBLK) + 0.5) / NBLK
    night = sza_deg(tb, lat, lon) > SZA_NIGHT
    alt_st = float(np.nanmedian(A["alt"]))
    tilt_a = float(np.nanmedian(A["tilt"]))
    tilt_b = float(np.nanmedian(B["tilt"]))
    exp_a = np.nanmedian(A["nprof"][A["nprof"] > 0])
    exp_b = np.nanmedian(B["nprof"][B["nprof"] > 0])
    oka = npa >= 0.5 * exp_a
    okb = npb >= 0.5 * exp_b
    cloud_a = dilate(cla > 0, 1)
    cloud_b = dilate(clb > 0, 1)
    clear_common = night & oka & okb & ~cloud_a & ~cloud_b
    clear_chm = night & oka & ~cloud_a
    nid = np.floor(tb + lon / 360.0 - 0.5).astype(int)

    dark_a = binned_dark("A", tilt_a) if site == "Payerne" else np.zeros(len(ZMID))
    dark_b = binned_dark("C", tilt_b)          # Payerne CL61 hood dark (template for other CL61)
    aer = aeronet_running(site)
    hkv_b = [str(x) for x in B["hk_vars"]]
    rows_c, rows_a = [], []
    for kind, mask in (("common", clear_common), ("chm", clear_chm)):
        nights = np.unique(nid[mask])
        for n_ in nights:
            sel = np.where(mask & (nid == n_))[0]
            if len(sel) < MIN_BLOCKS:
                continue
            t0, t1 = tb[sel[0]], tb[sel[-1]]
            prof, pwv, below = night_profiles(cams, t0, t1, alt_st)
            den = denominators(prof, alt_st, tilt_a, tilt_b)
            rma = np.nanmean(ra[sel], axis=0)
            row = dict(site=site, night=int(n_), t=float(np.mean(tb[sel])), n_blk=len(sel),
                       T_int=float(np.nanmean(hka[sel, 0]) - 273.15), pwv=pwv, cams_below=int(below))
            lx = band_lnx(rma, den["chm"])
            sg = band_sigma(ra[sel], den["chm"])
            for j, bn in enumerate(BAND_NAMES):
                row[f"chm_{bn}"] = lx[j]
                row[f"sig_chm_{bn}"] = sg[j]
            if site == "Payerne":
                lx = band_lnx(rma - dark_a, den["chm"])
                for j, bn in enumerate(BAND_NAMES):
                    row[f"chmdark_{bn}"] = lx[j]
            if kind == "common":
                rmb = np.nanmean(rb[sel], axis=0)
                variants = {"cl_nowv": (rmb, den["cl"]), "cl_wv": (rmb, den["cl_wv"]),
                            "cl_wvq": (rmb, den["cl_wvq"]), "cl_wvdark": (rmb - dark_b, den["cl_wv"])}
                for v, (num, dd) in variants.items():
                    lx = band_lnx(num, dd)
                    for j, bn in enumerate(BAND_NAMES):
                        row[f"{v}_{bn}"] = lx[j]
                sg = band_sigma(rb[sel], den["cl_wv"])
                for j, (a, b) in enumerate(BANDS):
                    bn = BAND_NAMES[j]
                    row[f"sig_cl_{bn}"] = sg[j]
                    m = (ZMID >= a) & (ZMID < b)
                    row[f"t2wv_{bn}"] = float(np.mean(den["t2op"][m]))
                    row[f"t2wvq_{bn}"] = float(np.mean(den["t2q"][m]))
                for k, v in enumerate(hkv_b):
                    row[f"cl_{v}"] = float(np.nanmean(hkb[sel, k]))
            if aer is not None:
                tq = np.array([row["t"]])
                row["tau1064"] = float(running_at(tq, aer[0], aer[1])[0]) if len(aer[0]) else np.nan
                row["tau910"] = float(running_at(tq, aer[0], aer[2])[0]) if len(aer[0]) else np.nan
            (rows_c if kind == "common" else rows_a).append(row)
    print(f"{site}: {len(rows_c)} common nights, {len(rows_a)} CHM15k-only nights", flush=True)
    return rows_c, rows_a


def write_csv(rows, path):
    keys = []
    for r in rows:
        for k in r:
            if k not in keys:
                keys.append(k)
    with open(path, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=keys)
        w.writeheader()
        w.writerows(rows)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=4)
    a = ap.parse_args()
    allc, alla = [], []
    with ProcessPoolExecutor(a.workers) as ex:
        for rc, ra in ex.map(process_site, list(SITES)):
            allc += rc
            alla += ra
    write_csv(allc, Q17 / "q17b_nights.csv")
    write_csv(alla, Q17 / "q17b_nights_chm.csv")


if __name__ == "__main__":
    main()
