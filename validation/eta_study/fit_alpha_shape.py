"""Step 2 of the eta study: effective in-cloud extinction alpha from the SHAPE of the in-cloud
backscatter at every accepted cloud-calibration profile, against the PVC forward table.

Inputs
  C:/DATA/Projects/202611_eta_study/incloud/<key>.npz   (extract_incloud_profiles.py)
  C:/DATA/Projects/202611_eta_study/forward_lut_ext.npz (forward_lut.py, window base-100..base+380 m)
  D:/Cloudnet_Droplet_radius                            (closure at Payerne / Palaiseau / Lindenberg)

Observable (identical for data and model)
  beta_before (attenuated, uncorrected) on the calibration gates (working grid: 10 m CL31/CL51,
  9.6 m CL61) and, for the CL61, also on the native 4.8 m gates. Window = gates whose centre lies in
  [peak-60, peak+300] m, peak = gate maximum (data: the peak the calibration's cloud filters use).
  Each profile is normalised by its sum over the window (= integral / gate). F50 / F100 = fraction of
  that integral on gates up to peak+50 / peak+100 m. Model: LUT gate profiles (pulse-convolved,
  gate-averaged; 5 gate phases) peak-aligned at THEIR gate maximum, normalised the same way, then
  averaged over the phases.

Fits per profile (least squares on the normalised shape, all window gates, unweighted)
  homogeneous LUT: alpha on a 121-point log grid 3..60 /km (LUT interpolated linearly in log alpha
  and in CBH, CBH node every 50 m), a_G = 3.3 um (sensitivity 2.5 and 5 um), each pulse width
  0/60/110/160 ns; edge flag when the best alpha is a grid end. alpha from F100 alone (same grid).
  adiabatic LUT: N on a 61-point log grid 50..400 cm-3 for f_ad 0.6 and 1.0, best (f_ad, N) per
  pulse; alpha_adi = mean extinction of the best cloud over the lowest 100 m above base.
  "chi2" = sum of squared residuals of the unit-sum shapes (no noise model: comparable between
  models, not an absolute goodness of fit).
Model-light quantities: k_app = slope of ln(beta / beta_peak) on gates +20..+150 m above the peak
  (gates with beta > 0.5 % of the peak, >= 3 gates); half-height distance above the peak.

Also: beta_before vs beta_after shape check (random profiles), CL61 native vs working,
  Payerne CL31/CL61 FOV contrast (pairs within +-1 min, |dCBH| <= 30 m), Cloudnet closure (+-15 min),
  summaries per unit / type / CBH / month, figures (a)-(d) and the step-3 eta preview.

Run:  python validation/eta_study/fit_alpha_shape.py [--workers 22]
Outputs: C:/DATA/Projects/202611_eta_study/step2/
"""
from __future__ import annotations

import os
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[_v] = "1"

import argparse
import json
import sys
import time as time_mod
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO))

# Paths: home-PC defaults; the CSCS array overrides ETA_BASE (and ETA_OUT for the extraction).
BASE = Path(os.environ.get("ETA_BASE", "C:/DATA/Projects/202611_eta_study"))
SRC = Path(os.environ.get("ETA_OUT", str(BASE / "incloud")))
LUT_FILE = Path(os.environ.get("ETA_LUT", str(BASE / "forward_lut_ext.npz")))
OUT = BASE / "step2"
CLOUDNET = Path(os.environ.get("ETA_CLOUDNET", "D:/Cloudnet_Droplet_radius"))

WIN_LO, WIN_HI = -60.0, 300.0
PULSES = np.array([0.0, 60.0, 110.0, 160.0])
NOMINAL = {"CL31": 110.0, "CL51": 110.0, "CL61": 160.0}
AG_FIT = (2.5, 3.3, 5.0)
ALPHA_FINE = np.logspace(np.log10(3.0), np.log10(60.0), 121)
N_FINE = np.logspace(np.log10(50.0), np.log10(400.0), 61)
CBH_NODES = np.arange(200.0, 3000.0 + 1, 50.0)
KAPP_LO, KAPP_HI, KAPP_MIN_FRAC = 20.0, 150.0, 0.005
MIN_PEAK_M = 160.0                     # peak - 60 m must stay above the 100 m lowest usable gate
SIGMA_LN = 0.35
CN_SITES = {"0-20000-0-06610_C": "payerne", "0-20000-0-06610_B": "payerne",
            "0-20000-0-10393_C": "lindenberg", "0-250-1001-07151_A": "palaiseau"}
CN_WINDOW_S = 900.0
BASE_SEARCH_M = 150.0

_LUT = None
_MODELS = {}


# =============================================================================
#  Model construction (cached per worker)
# =============================================================================
def lut():
    global _LUT
    if _LUT is None:
        z = np.load(LUT_FILE, allow_pickle=True)
        _LUT = {k: z[k] for k in z.files if k != "README"}
    return _LUT


def gate_offsets(dg):
    return np.arange(int(np.ceil(WIN_LO / dg - 1e-6)), int(np.floor(WIN_HI / dg + 1e-6)) + 1)


def peak_aligned(arr, ks):
    """arr (..., phase, gate) -> phase-averaged, peak-aligned, unit-sum shapes (..., K)."""
    sh = arr.shape
    a = np.nan_to_num(arr.reshape(-1, sh[-2], sh[-1]).astype("float64"))
    pk = a.argmax(axis=-1)
    idx = pk[..., None] + ks
    if idx.min() < 0 or idx.max() >= sh[-1]:
        raise RuntimeError("LUT gate window too short for the peak-aligned window")
    w = np.take_along_axis(a, idx, axis=-1)
    w /= w.sum(-1, keepdims=True)
    m = w.mean(axis=1)
    m /= m.sum(-1, keepdims=True)
    return m.reshape(sh[:-2] + (ks.size,))


def _interp_axis(arr, x_old, x_new, axis):
    """Linear interpolation of `arr` along `axis` from x_old to x_new (both ascending)."""
    i = np.clip(np.searchsorted(x_old, x_new) - 1, 0, x_old.size - 2)
    w = (x_new - x_old[i]) / (x_old[i + 1] - x_old[i])
    w = np.clip(w, 0.0, 1.0)
    a0 = np.take(arr, i, axis=axis)
    a1 = np.take(arr, i + 1, axis=axis)
    shape = [1] * arr.ndim
    shape[axis] = -1
    w = w.reshape(shape)
    return a0 * (1 - w) + a1 * w


PHASES = np.arange(5) / 5.0


def gate_sample_conv(conv, z, dg):
    """Gate means (5 phases) of 1 m pulse-convolved profiles conv (..., z): gate centres at
    (k + p) dg relative to cloud base, only gates lying fully inside the stored z range. Same
    definition as forward_lut.gate_sample; used when a unit's gate differs from the LUT's stored
    gates (CL31 units configured at 20 m)."""
    c0 = np.nan_to_num(conv.astype("float64"))
    cum = np.concatenate((np.zeros(c0.shape[:-1] + (1,)), np.cumsum(c0, axis=-1)), axis=-1)
    e0 = z[0] - 0.5                                                  # first edge, 1 m spacing

    def G(x):
        f = x - e0
        i = np.clip(np.floor(f).astype(int), 0, cum.shape[-1] - 2)
        w = f - i
        return cum[..., i] * (1 - w) + cum[..., i + 1] * w

    k = np.arange(int(np.ceil((z[0] + dg) / dg)), int(np.floor((z[-1] - dg) / dg)) + 1)
    out = [(G((k + p) * dg + dg / 2) - G((k + p) * dg - dg / 2)) / dg for p in PHASES]
    return np.stack(out, axis=-2)                                    # (..., phase, gate)


def models(itype, grid, dg_data=None):
    """Fine model tensors for one type/grid/gate spacing (cached)."""
    L = lut()
    g = "native" if grid == "native" else "work"
    dg_lut = float(L[f"{itype}_{g}_gate_m"])
    dg = dg_lut if dg_data is None else float(dg_data)
    key = (itype, grid, round(dg, 2))
    if key in _MODELS:
        return _MODELS[key]
    ks = gate_offsets(dg)
    if abs(dg - dg_lut) < 0.05:
        S = peak_aligned(L[f"{itype}_{g}_beta"], ks)             # (alpha16, aG8, cbh15, pulse4, K)
        Sa_raw = L[f"{itype}_adi_{g}_beta"]
    else:
        conv = L["beta_conv"][list(L["types"]).index(itype)]        # (alpha, aG, cbh, pulse, z)
        S = peak_aligned(gate_sample_conv(conv, L["z_rel_m"], dg), ks)
        Sa_raw = gate_sample_conv(L[f"{itype}_adi_conv"], L["z_rel_m"], dg)
    ag = L["aG_um"]
    sel = []
    for a in AG_FIT:
        sel.append(_interp_axis(S, ag, np.array([a]), axis=1)[:, 0])  # (alpha, cbh, pulse, K)
    S3 = np.stack(sel, axis=1)                                     # (alpha, aG3, cbh, pulse, K)
    S3 = _interp_axis(S3, np.log(L["alpha_km"]), np.log(ALPHA_FINE), axis=0)
    S3 = _interp_axis(S3, L["cbh_m"], CBH_NODES, axis=2)            # (alphaF, aG3, cbhN, pulse, K)
    Mh = np.transpose(S3, (2, 3, 1, 0, 4))                         # (cbhN, pulse, aG3, alphaF, K)
    Mh = Mh / Mh.sum(-1, keepdims=True)
    # adiabatic: (N4, fad2, cbh15, pulse4, phase, gate)
    Sa = peak_aligned(Sa_raw, ks)                                   # (N4, fad2, cbh15, pulse4, K)
    Sa = _interp_axis(Sa, np.log(L["adi_N_cm3"]), np.log(N_FINE), axis=0)
    Sa = _interp_axis(Sa, L["cbh_m"], CBH_NODES, axis=2)            # (NF, fad, cbhN, pulse, K)
    Ma = np.transpose(Sa, (2, 3, 1, 0, 4))                         # (cbhN, pulse, fad, NF, K)
    Ma = Ma / Ma.sum(-1, keepdims=True)
    z = L["z_rel_m"]
    al = L[f"{itype}_adi_alpha"]                                    # (N4, fad2, cbh15, z) m-1
    a100 = np.nanmean(np.where((z >= 0) & (z <= 100), al, np.nan), axis=-1) * 1e3   # /km
    a100 = np.exp(_interp_axis(np.log(a100), np.log(L["adi_N_cm3"]), np.log(N_FINE), axis=0))
    a100 = _interp_axis(a100, L["cbh_m"], CBH_NODES, axis=2)       # (NF, fad, cbhN)
    A100 = np.transpose(a100, (2, 1, 0))                           # (cbhN, fad, NF)
    zg = ks * dg
    F50 = Mh[..., zg <= 50 + 1e-6].sum(-1)
    F100 = Mh[..., zg <= 100 + 1e-6].sum(-1)
    # model tail slope k_app(alpha) at a_G 3.3 um, same estimator as the data (cbhN, pulse, alphaF)
    kapp_m = model_kapp(Mh[:, :, 1], zg)[0]
    _MODELS[key] = dict(ks=ks, dg=dg, zg=zg, Mh=Mh, Ma=Ma, A100=A100, F50=F50, F100=F100,
                        kapp_m=kapp_m, from_conv=abs(dg - dg_lut) >= 0.05)
    return _MODELS[key]


def alpha_from_kapp(kapp_obs, cbh, md, pulse_idx):
    """Tail-slope alpha: invert the model k_app(alpha) at each profile's CBH node. A Gaussian
    response preserves an exponential slope, so this estimator is (nearly) independent of the
    unknown pulse/receiver width -- unlike the full-shape fit. NaN outside the 3-60 /km grid."""
    node = np.clip(np.round((cbh - CBH_NODES[0]) / 50.0).astype(int), 0, CBH_NODES.size - 1)
    out = np.full(kapp_obs.size, np.nan)
    for c in np.unique(node):
        m = node == c
        kk = md["kapp_m"][c, pulse_idx]
        ok = np.isfinite(kk)
        order = np.argsort(kk[ok])
        out[m] = np.interp(kapp_obs[m], kk[ok][order], ALPHA_FINE[ok][order], left=np.nan, right=np.nan)
    return out


# =============================================================================
#  Per-profile observables
# =============================================================================
def observables(P, zg):
    """P (n, K) normalised (unit sum) -> F50, F100, k_app [1/km], n_k gates, half-height [m]."""
    F50 = P[:, zg <= 50 + 1e-6].sum(1)
    F100 = P[:, zg <= 100 + 1e-6].sum(1)
    k0 = int(np.argmin(np.abs(zg)))
    y = P / P[:, k0:k0 + 1]
    sel = (zg >= KAPP_LO - 1e-6) & (zg <= KAPP_HI + 1e-6)
    ys, zs = y[:, sel], zg[sel]
    ok = ys > KAPP_MIN_FRAC
    lny = np.where(ok, np.log(np.where(ok, ys, 1.0)), 0.0)
    n = ok.sum(1)
    sz = (ok * zs).sum(1)
    szz = (ok * zs ** 2).sum(1)
    sy = lny.sum(1)
    szy = (lny * zs).sum(1)
    with np.errstate(all="ignore"):
        slope = (n * szy - sz * sy) / (n * szz - sz ** 2)
    kapp = np.where(n >= 3, -slope * 1e3, np.nan)                  # 1/km
    above = zg > 0
    ya, za = y[:, above], zg[above]
    below = ya < 0.5
    first = np.where(below.any(1), below.argmax(1), -1)
    hh = np.full(P.shape[0], np.nan)
    for i in np.nonzero(first >= 0)[0]:
        j = first[i]
        y1, z1 = (1.0, 0.0) if j == 0 else (ya[i, j - 1], za[j - 1])
        y2, z2 = ya[i, j], za[j]
        hh[i] = z1 + (z2 - z1) * (y1 - 0.5) / (y1 - y2) if y1 != y2 else z2
    return F50, F100, kapp, n, hh


def model_kapp(M, zg):
    """k_app of model shapes (..., K) with the same definition as the data."""
    sh = M.shape
    F50, F100, kapp, n, hh = observables(M.reshape(-1, sh[-1]), zg)
    return kapp.reshape(sh[:-1]), hh.reshape(sh[:-1])


def fit_block(P, cbh, md):
    """Least-squares fits for profiles P (n, K) at CBHs `cbh` -> dict of arrays."""
    n = P.shape[0]
    node = np.clip(np.round((cbh - CBH_NODES[0]) / 50.0).astype(int), 0, CBH_NODES.size - 1)
    res = {f"alpha_p{int(p)}_ag{a}": np.full(n, np.nan, np.float32) for p in PULSES for a in AG_FIT}
    res.update({f"chi2_p{int(p)}_ag{a}": np.full(n, np.nan, np.float32) for p in PULSES for a in AG_FIT})
    res.update({f"edge_p{int(p)}_ag{a}": np.zeros(n, np.int8) for p in PULSES for a in AG_FIT})
    res.update({f"alphaF100_p{int(p)}": np.full(n, np.nan, np.float32) for p in PULSES})
    res.update({f"alpha_adi_p{int(p)}": np.full(n, np.nan, np.float32) for p in PULSES})
    res.update({f"chi2_adi_p{int(p)}": np.full(n, np.nan, np.float32) for p in PULSES})
    res.update({f"N_adi_p{int(p)}": np.full(n, np.nan, np.float32) for p in PULSES})
    res.update({f"fad_adi_p{int(p)}": np.full(n, np.nan, np.float32) for p in PULSES})
    zg = md["zg"]
    f100 = P[:, zg <= 100 + 1e-6].sum(1)
    pp = (P ** 2).sum(1)
    for c in np.unique(node):
        m = node == c
        Pc = P[m]
        Mh = md["Mh"][c]                                            # (pulse, aG3, alphaF, K)
        flat = Mh.reshape(-1, Mh.shape[-1])
        sse = pp[m][:, None] + (flat ** 2).sum(1)[None, :] - 2.0 * Pc @ flat.T
        sse = sse.reshape((Pc.shape[0],) + Mh.shape[:-1])           # (n, pulse, aG, alphaF)
        best = sse.argmin(-1)
        for ip, p in enumerate(PULSES):
            for ia, a in enumerate(AG_FIT):
                b = best[:, ip, ia]
                res[f"alpha_p{int(p)}_ag{a}"][m] = ALPHA_FINE[b]
                res[f"chi2_p{int(p)}_ag{a}"][m] = sse[np.arange(b.size), ip, ia, b]
                res[f"edge_p{int(p)}_ag{a}"][m] = ((b == 0) | (b == ALPHA_FINE.size - 1)).astype(np.int8)
            fm = md["F100"][c, ip, 1]                               # (alphaF,), a_G 3.3
            bF = np.abs(f100[m][:, None] - fm[None, :]).argmin(1)
            res[f"alphaF100_p{int(p)}"][m] = ALPHA_FINE[bF]
        Ma = md["Ma"][c]                                            # (pulse, fad, NF, K)
        flat = Ma.reshape(-1, Ma.shape[-1])
        sse = pp[m][:, None] + (flat ** 2).sum(1)[None, :] - 2.0 * Pc @ flat.T
        sse = sse.reshape((Pc.shape[0], Ma.shape[0], -1))           # (n, pulse, fad*NF)
        best = sse.argmin(-1)
        for ip, p in enumerate(PULSES):
            b = best[:, ip]
            fi, ni = np.divmod(b, N_FINE.size)
            res[f"chi2_adi_p{int(p)}"][m] = sse[np.arange(b.size), ip, b]
            res[f"N_adi_p{int(p)}"][m] = N_FINE[ni]
            res[f"fad_adi_p{int(p)}"][m] = np.array([0.6, 1.0])[fi]
            res[f"alpha_adi_p{int(p)}"][m] = md["A100"][c, fi, ni]
    return res


# =============================================================================
#  Worker 1: one unit / grid
# =============================================================================
def fit_unit(job):
    key, u, grid = job
    import pandas as pd
    z = np.load(SRC / f"{key}.npz", allow_pickle=True)
    pre = "native_" if grid == "native" else ""
    B = z[f"{pre}beta_before"].astype("float64")
    Ba = z[f"{pre}beta_after"].astype("float64")
    depth = z[f"{pre}range_rel_m"].astype("float64")
    md = models(u["itype"], grid, float(z[f"{pre}dr_m"]))
    dg = md["dg"]
    # data window: the stored offsets that match the model window
    k_data = np.round(depth / dg).astype(int)
    cols = np.array([int(np.nonzero(k_data == k)[0][0]) for k in md["ks"]])
    W, Wa = B[:, cols], Ba[:, cols]
    peak = (z["native_peak_range_m"] if grid == "native" else z["peak_range_m"]).astype("float64")
    cbh = z["cbh_m"].astype("float64")
    cbh_use = np.where(np.isfinite(cbh), cbh, peak)
    tot = W.sum(1)
    valid = np.isfinite(W).all(1) & (tot > 0) & (peak >= MIN_PEAK_M)
    P = np.where(valid[:, None], W / np.where(tot > 0, tot, 1.0)[:, None], np.nan)
    out = dict(key=np.full(W.shape[0], key), itype=np.full(W.shape[0], u["itype"]),
               site=np.full(W.shape[0], u["site"]), distorted=np.full(W.shape[0], bool(u["distorted"])),
               grid=np.full(W.shape[0], grid), time_ns=z["time_ns"], date=z["date"],
               cbh_m=cbh.astype(np.float32), peak_range_m=peak.astype(np.float32),
               valid=valid, c_profile=z["c_profile"])
    F50 = np.full(W.shape[0], np.nan)
    F100 = F50.copy()
    kapp = F50.copy()
    nk = np.zeros(W.shape[0], int)
    hh = F50.copy()
    iv = np.nonzero(valid)[0]
    fits = {}
    if iv.size:
        f50, f100, ka, n_, h_ = observables(P[iv], md["zg"])
        F50[iv], F100[iv], kapp[iv], nk[iv], hh[iv] = f50, f100, ka, n_, h_
        for s in range(0, iv.size, 20000):
            blk = iv[s:s + 20000]
            r = fit_block(P[blk], cbh_use[blk], md)
            for k, v in r.items():
                fits.setdefault(k, np.full(W.shape[0], np.nan if v.dtype.kind == "f" else 0, v.dtype))[blk] = v
    alpha_k = alpha_from_kapp(kapp, cbh_use, md, list(PULSES).index(NOMINAL[u["itype"]]))
    out.update(F50=F50.astype(np.float32), F100=F100.astype(np.float32), k_app=kapp.astype(np.float32),
               n_kapp_gates=nk.astype(np.int16), half_height_m=hh.astype(np.float32),
               alpha_k=alpha_k.astype(np.float32), **fits)
    df = pd.DataFrame(out)
    (OUT / "profiles").mkdir(parents=True, exist_ok=True)
    try:
        import pyarrow  # noqa: F401
        df.to_parquet(OUT / "profiles" / f"{key}_{grid}.parquet", compression="zstd", index=False)
    except ImportError:                                  # CSCS 'alc' env has no pyarrow
        np.savez_compressed(OUT / "profiles" / f"{key}_{grid}.npz",
                            **{c: df[c].to_numpy() for c in df.columns})
    # before-vs-after check on up to 40 random valid profiles of this unit
    chk = None
    if iv.size:
        rng = np.random.default_rng(abs(hash(key)) % 2 ** 32)
        pick = rng.choice(iv, size=min(40, iv.size), replace=False)
        Pa = Wa[pick] / Wa[pick].sum(1, keepdims=True)
        okp = np.isfinite(Pa).all(1)
        pick, Pa = pick[okp], Pa[okp]
        if pick.size:
            fa50, fa100, kaa, _, _ = observables(Pa, md["zg"])
            ra = fit_block(Pa, cbh_use[pick], md)
            p_nom = int(NOMINAL[u["itype"]])
            chk = dict(key=key, grid=grid, n=int(pick.size),
                       max_abs_shape_diff=float(np.nanmax(np.abs(Pa - P[pick]))),
                       dF50_med=float(np.nanmedian(fa50 - F50[pick])),
                       dF100_med=float(np.nanmedian(fa100 - F100[pick])),
                       alpha_ratio_after_before_med=float(np.nanmedian(
                           ra[f"alpha_p{p_nom}_ag3.3"] / fits[f"alpha_p{p_nom}_ag3.3"][pick])))
    return key, grid, int(valid.sum()), int(W.shape[0]), chk


# =============================================================================
#  Worker 2: Cloudnet closure, one site-day
# =============================================================================
def cloudnet_day(job):
    """job = (cn_site, ds, [(row, t_epoch_s, cbh_m)]) -> [(row, alpha_cn_km, aG_cn_um, n_cn)]."""
    from netCDF4 import Dataset
    cn_site, ds, rows = job
    fp = CLOUDNET / f"{ds}_{cn_site}_der.nc"
    if not fp.exists():
        return []
    day0 = np.datetime64(f"{ds[:4]}-{ds[4:6]}-{ds[6:]}T00:00:00").astype("datetime64[s]").astype(np.int64)
    ts = np.array([r[1] for r in rows], float)
    try:
        with Dataset(fp) as nc:
            cnt = np.asarray(nc.variables["time"][:], float) * 3600.0 + day0
            h = np.asarray(nc.variables["height"][:], float)
            salt = float(np.asarray(nc.variables["altitude"][:]).ravel()[0])
            i0 = int(max(np.searchsorted(cnt, ts.min() - CN_WINDOW_S) - 1, 0))
            i1 = int(min(np.searchsorted(cnt, ts.max() + CN_WINDOW_S) + 1, cnt.size))
            sl = (slice(i0, i1), slice(None))

            def rd(name, fill=np.nan, dtype=float):
                return np.ma.filled(np.ma.asarray(nc.variables[name][sl], dtype=dtype), fill)
            st = rd("der_retrieval_status", 0, int)
            rs = rd("der_scaled")
            nn = rd("N_scaled")
    except Exception:  # noqa: BLE001
        return []
    cnt = cnt[i0:i1]
    dz = float(np.median(np.diff(h)))
    liq = (st >= 1) & (st <= 3)
    rs = np.where((rs > 0) & (rs < 1e-3), rs, np.nan)
    nn = np.where(nn > 0, nn, np.nan)
    alpha_px = 2 * np.pi * nn * rs ** 2 * np.exp(-3 * SIGMA_LN ** 2)
    ag_px = rs * np.exp(-1.5 * SIGMA_LN ** 2)
    # per Cloudnet profile: contiguous liquid runs -> (bottom gate, top gate)
    runs = []
    for r in range(liq.shape[0]):
        x = np.diff(np.concatenate(([0], liq[r].astype(int), [0])))
        runs.append(list(zip(np.nonzero(x == 1)[0], np.nonzero(x == -1)[0] - 1)))
    memo = {}

    def run_value(r, jb, jt):
        k = (r, jb)
        if k not in memo:
            base_edge = h[jb] - dz / 2
            g = np.arange(jb, jt + 1)
            g = g[(h[g] - base_edge) < 100.0]
            g = g[st[r, g] == 1]
            a = alpha_px[r, g]
            a = a[np.isfinite(a)]
            ag = ag_px[r, g]
            ag = ag[np.isfinite(ag)]
            memo[k] = (np.median(a) if a.size else np.nan, np.median(ag) if ag.size else np.nan)
        return memo[k]

    out = []
    for (row, t, cbh) in rows:
        cloud_msl = cbh + salt
        lo, hi = np.searchsorted(cnt, t - CN_WINDOW_S), np.searchsorted(cnt, t + CN_WINDOW_S, side="right")
        vals, ags = [], []
        for r in range(lo, hi):
            for jb, jt in runs[r]:
                if h[jt] >= cloud_msl - BASE_SEARCH_M and h[jb] <= cloud_msl + BASE_SEARCH_M:
                    a, ag = run_value(r, jb, jt)
                    if np.isfinite(a):
                        vals.append(a)
                        ags.append(ag)
                    break                                           # lowest run near the CBH only
        if vals:
            out.append((row, float(np.median(vals)) * 1e3, float(np.nanmedian(ags)) * 1e6, len(vals)))
    return out


# =============================================================================
#  Driver + summaries
# =============================================================================
def q(x, p):
    x = np.asarray(x, float)
    x = x[np.isfinite(x)]
    return float(np.percentile(x, p)) if x.size else np.nan


def main():
    import pandas as pd
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=22)
    ap.add_argument("--skip-fit", action="store_true", help="reuse step2/profiles/*.parquet")
    ap.add_argument("--unit", help="CSCS array task: fit ONE unit (working + native for CL61), no summaries")
    a = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    t0 = time_mod.time()
    if a.unit:
        u = json.load(open(SRC / "units" / f"{a.unit}.json", encoding="utf-8"))
        if not (SRC / f"{a.unit}.npz").exists():
            print(f"{a.unit}: no extracted profiles - nothing to fit", flush=True)
            return
        for grid in (("working", "native") if u["itype"] == "CL61" else ("working",)):
            key, grid, nv, nt, chk = fit_unit((a.unit, u, grid))
            print(f"{key} {grid}: {nv}/{nt} valid ({time_mod.time() - t0:.0f} s) check={chk}", flush=True)
        return
    log = open(OUT / "run_summary.txt", "w", encoding="utf-8")

    def say(s=""):
        print(s, flush=True)
        log.write(s + "\n")
        log.flush()

    if (SRC / "units.json").exists():
        units = json.load(open(SRC / "units.json", encoding="utf-8"))
    else:                                                   # CSCS array: one json per unit
        units = {p.stem: json.load(open(p, encoding="utf-8")) for p in sorted((SRC / "units").glob("*.json"))}
    units = {k: u for k, u in units.items() if (SRC / f"{k}.npz").exists()}
    jobs = [(k, u, "working") for k, u in units.items()]
    jobs += [(k, u, "native") for k, u in units.items() if u["itype"] == "CL61"]
    jobs.sort(key=lambda j: -(SRC / f"{j[0]}.npz").stat().st_size)
    checks = []
    if not a.skip_fit:
        say(f"[fit] {len(jobs)} unit/grid jobs -> {a.workers} processes")
        with ProcessPoolExecutor(max_workers=a.workers) as ex:
            for key, grid, nv, nt, chk in ex.map(fit_unit, jobs, chunksize=1):
                say(f"  {key} {grid}: {nv}/{nt} valid ({time_mod.time() - t0:.0f} s)")
                if chk:
                    checks.append(chk)
        pd.DataFrame(checks).to_csv(OUT / "check_before_vs_after.csv", index=False, float_format="%.5g")
    tables = [pd.read_parquet(p) for p in sorted((OUT / "profiles").glob("*.parquet"))]
    tables += [pd.DataFrame(dict(np.load(p, allow_pickle=True))) for p in sorted((OUT / "profiles").glob("*.npz"))]
    df = pd.concat(tables, ignore_index=True)
    say(f"[fit] {len(df)} profiles, {int(df.valid.sum())} valid  ({time_mod.time() - t0:.0f} s)")

    # ---------------- Cloudnet closure (+-15 min) ----------------
    w = df[(df.grid == "working") & df.valid & df.key.isin(list(CN_SITES))]
    grp = defaultdict(list)
    for row, k, t, cb in zip(w.index, w.key, w.time_ns, w.cbh_m):
        tt = int(t // 1_000_000_000)
        ds = str(np.datetime64(tt, "s").astype("datetime64[D]")).replace("-", "")
        if np.isfinite(cb):
            grp[(CN_SITES[k], ds)].append((row, tt, float(cb)))
    cjobs = [(s, d, v) for (s, d), v in sorted(grp.items())]
    say(f"[closure] {len(w)} profiles on {len(cjobs)} site-days")
    df["alpha_cn_km"] = np.nan
    df["aG_cn_um"] = np.nan
    df["n_cn"] = 0
    rows, acn, agcn, ncn = [], [], [], []
    with ProcessPoolExecutor(max_workers=a.workers) as ex:
        for res in ex.map(cloudnet_day, cjobs, chunksize=2):
            for row, al, ag, n in res:
                rows.append(row)
                acn.append(al)
                agcn.append(ag)
                ncn.append(n)
    df.loc[rows, "alpha_cn_km"] = acn
    df.loc[rows, "aG_cn_um"] = agcn
    df.loc[rows, "n_cn"] = ncn
    say(f"[closure] {len(rows)} profiles with Cloudnet alpha ({time_mod.time() - t0:.0f} s)")
    # per-profile fits live in profiles/<key>_<grid>.parquet; only the closure columns are added here
    df.loc[rows, ["key", "grid", "time_ns", "alpha_cn_km", "aG_cn_um", "n_cn"]].to_parquet(
        OUT / "closure_profiles.parquet", compression="zstd", index=False)

    summarise(df, say)
    figures(df, say)
    eta_preview(df, say)
    say(f"total {time_mod.time() - t0:.0f} s")
    log.close()


def nominal_cols(itype):
    p = int(NOMINAL[itype])
    return f"alpha_p{p}_ag3.3", f"edge_p{p}_ag3.3", f"chi2_p{p}_ag3.3"


def add_nominal(df):
    for c in ("alpha_nom", "edge_nom", "chi2_nom", "alpha_adi_nom", "chi2_adi_nom", "alphaF100_nom"):
        df[c] = np.nan
    for it in NOMINAL:
        m = df.itype == it
        p = int(NOMINAL[it])
        df.loc[m, "alpha_nom"] = df.loc[m, f"alpha_p{p}_ag3.3"]
        df.loc[m, "edge_nom"] = df.loc[m, f"edge_p{p}_ag3.3"]
        df.loc[m, "chi2_nom"] = df.loc[m, f"chi2_p{p}_ag3.3"]
        df.loc[m, "alpha_adi_nom"] = df.loc[m, f"alpha_adi_p{p}"]
        df.loc[m, "chi2_adi_nom"] = df.loc[m, f"chi2_adi_p{p}"]
        df.loc[m, "alphaF100_nom"] = df.loc[m, f"alphaF100_p{p}"]
    return df


GROUPS = [("CL31 all", lambda d: d.itype == "CL31"),
          ("CL31 excl. Payerne CL31", lambda d: (d.itype == "CL31") & ~d.distorted.astype(bool)),
          ("CL51 all", lambda d: d.itype == "CL51"),
          ("CL61 all", lambda d: d.itype == "CL61"),
          ("Vaisala excl. Payerne CL31", lambda d: ~d.distorted.astype(bool))]


def summarise(df, say):
    import pandas as pd
    from scipy.stats import spearmanr
    df = add_nominal(df)
    v = df[df.valid]
    vw = v[v.grid == "working"]

    say("\n== alpha_eff (homogeneous PVC, a_G 3.3 um, nominal pulse 110/110/160 ns; working gates) ==")
    rows = []
    for (k, it, site, dist), g in vw.groupby(["key", "itype", "site", "distorted"]):
        ne = g[g.edge_nom == 0]
        rows.append(dict(level="unit", group=f"{it} {site} ({k})", distorted=dist, n=len(g),
                         frac_edge=float((g.edge_nom == 1).mean()), n_noedge=len(ne),
                         alpha_med=q(ne.alpha_nom, 50), alpha_p10=q(ne.alpha_nom, 10), alpha_p90=q(ne.alpha_nom, 90),
                         alpha_all_med=q(g.alpha_nom, 50), alphaF100_med=q(g.alphaF100_nom, 50),
                         alpha_adi_med=q(g.alpha_adi_nom, 50), F50_med=q(g.F50, 50), F100_med=q(g.F100, 50),
                         kapp_med=q(g.k_app, 50), half_height_med=q(g.half_height_m, 50),
                         chi2_hom_med=q(g.chi2_nom, 50), chi2_adi_med=q(g.chi2_adi_nom, 50),
                         alpha_k_med=q(g.alpha_k, 50), alpha_k_p10=q(g.alpha_k, 10), alpha_k_p90=q(g.alpha_k, 90),
                         frac_alpha_k_outside=float(g.alpha_k.isna().mean()), cbh_med=q(g.cbh_m, 50)))
    for gname, sel in GROUPS:
        g = vw[sel(vw)]
        ne = g[g.edge_nom == 0]
        um = ne.groupby("key").alpha_nom.median()
        rows.append(dict(level="type", group=gname, distorted=False, n=len(g), n_units=g.key.nunique(),
                         frac_edge=float((g.edge_nom == 1).mean()), n_noedge=len(ne),
                         alpha_med=q(ne.alpha_nom, 50), alpha_p10=q(ne.alpha_nom, 10), alpha_p90=q(ne.alpha_nom, 90),
                         alpha_unit_med_of_meds=q(um, 50), alpha_all_med=q(g.alpha_nom, 50),
                         alphaF100_med=q(g.alphaF100_nom, 50), alpha_adi_med=q(g.alpha_adi_nom, 50),
                         F50_med=q(g.F50, 50), F100_med=q(g.F100, 50), kapp_med=q(g.k_app, 50),
                         half_height_med=q(g.half_height_m, 50), chi2_hom_med=q(g.chi2_nom, 50),
                         chi2_adi_med=q(g.chi2_adi_nom, 50), frac_adi_better=float((g.chi2_adi_nom < g.chi2_nom).mean()),
                         alpha_k_med=q(g.alpha_k, 50), alpha_k_p10=q(g.alpha_k, 10), alpha_k_p90=q(g.alpha_k, 90),
                         alpha_k_unit_med_of_meds=q(g.groupby("key").alpha_k.median(), 50),
                         frac_alpha_k_outside=float(g.alpha_k.isna().mean()), cbh_med=q(g.cbh_m, 50)))
    summ = pd.DataFrame(rows)
    summ.to_csv(OUT / "summary_alpha_by_unit_type.csv", index=False, float_format="%.4g")
    say(summ[summ.level == "type"].drop(columns=["level", "distorted"]).to_string(index=False))

    say("\n== pulse and a_G sensitivity: median alpha (no-edge at each setting) / frac at grid edge ==")
    prow = []
    for gname, sel in GROUPS:
        g = vw[sel(vw)]
        r = dict(group=gname, n=len(g))
        for p in PULSES:
            for ag in AG_FIT:
                c, e = f"alpha_p{int(p)}_ag{ag}", f"edge_p{int(p)}_ag{ag}"
                r[f"p{int(p)}_ag{ag}"] = q(g.loc[g[e] == 0, c], 50)
                r[f"p{int(p)}_ag{ag}_edge"] = float((g[e] == 1).mean())
            r[f"p{int(p)}_chi2hom"] = q(g[f"chi2_p{int(p)}_ag3.3"], 50)
            r[f"p{int(p)}_chi2adi"] = q(g[f"chi2_adi_p{int(p)}"], 50)
            r[f"p{int(p)}_alpha_adi"] = q(g[f"alpha_adi_p{int(p)}"], 50)
        prow.append(r)
    ps = pd.DataFrame(prow)
    ps.to_csv(OUT / "summary_pulse_aG_sensitivity.csv", index=False, float_format="%.4g")
    show = ["group"] + [f"p{int(p)}_ag3.3" for p in PULSES] + [f"p{int(p)}_ag3.3_edge" for p in PULSES]
    say(ps[show].to_string(index=False))
    show = ["group"] + [f"p{int(p)}_ag{ag}" for p in (110, 160) for ag in AG_FIT]
    say(ps[show].to_string(index=False))
    show = ["group"] + [f"p{int(p)}_chi2hom" for p in PULSES] + [f"p{int(p)}_chi2adi" for p in PULSES]
    say(ps[show].to_string(index=False))

    say("\n== model-light: k_app (1/km, +20..+150 m) and half-height (m) vs the LUT ==")
    L = lut()
    lrows = []
    for it in ("CL31", "CL51", "CL61"):
        for grid in ("working", "native") if it == "CL61" else ("working",):
            md = models(it, grid)
            ka, hh = model_kapp(md["Mh"], md["zg"])                   # (cbhN, pulse, aG, alphaF)
            F50m, F100m = md["F50"], md["F100"]
            c = int(np.argmin(np.abs(CBH_NODES - 1000.0)))
            for p_i, p in enumerate(PULSES):
                for al in (5.0, 10.0, 20.0, 40.0):
                    ai = int(np.argmin(np.abs(ALPHA_FINE - al)))
                    lrows.append(dict(itype=it, grid=grid, pulse_ns=p, alpha_km=al, cbh_m=1000,
                                      kapp_model=ka[c, p_i, 1, ai], half_height_model=hh[c, p_i, 1, ai],
                                      F50_model=F50m[c, p_i, 1, ai], F100_model=F100m[c, p_i, 1, ai]))
    lk = pd.DataFrame(lrows)
    lk.to_csv(OUT / "lut_kapp_F_table.csv", index=False, float_format="%.4g")
    say(lk.pivot_table(index=["itype", "grid", "pulse_ns"], columns="alpha_km",
                       values=["kapp_model", "half_height_model", "F100_model"]).round(2).to_string())
    obs = vw.groupby("itype")[["k_app", "half_height_m", "F50", "F100"]].median()
    say("observed medians (working gates):\n" + obs.round(3).to_string())
    nat = v[v.grid == "native"][["k_app", "half_height_m", "F50", "F100"]].median()
    say("observed medians CL61 native 4.8 m:\n" + nat.round(3).to_string())

    say("\n== CL61: native 4.8 m vs working 9.6 m (same profiles) ==")
    cl = v[v.itype == "CL61"]
    pw = cl[cl.grid == "working"].set_index(["key", "time_ns"])
    pn = cl[cl.grid == "native"].set_index(["key", "time_ns"])
    j = pw.join(pn, lsuffix="_w", rsuffix="_n", how="inner")
    nrow = []
    for p in PULSES:
        cw, cn = f"alpha_p{int(p)}_ag3.3_w", f"alpha_p{int(p)}_ag3.3_n"
        ok = (j[f"edge_p{int(p)}_ag3.3_w"] == 0) & (j[f"edge_p{int(p)}_ag3.3_n"] == 0)
        r = j.loc[ok, cn] / j.loc[ok, cw]
        nrow.append(dict(pulse_ns=p, n_pairs=int(ok.sum()), alpha_working_med=q(j.loc[ok, cw], 50),
                         alpha_native_med=q(j.loc[ok, cn], 50), ratio_native_working_med=q(r, 50),
                         ratio_p10=q(r, 10), ratio_p90=q(r, 90)))
    nw = pd.DataFrame(nrow)
    nw.to_csv(OUT / "cl61_native_vs_working.csv", index=False, float_format="%.4g")
    say(nw.to_string(index=False))
    per_unit = []
    for k, g in j.groupby(level=0):
        ok = (g["edge_p160_ag3.3_w"] == 0) & (g["edge_p160_ag3.3_n"] == 0)
        per_unit.append(dict(key=k, n=int(ok.sum()), alpha_w=q(g.loc[ok, "alpha_p160_ag3.3_w"], 50),
                             alpha_n=q(g.loc[ok, "alpha_p160_ag3.3_n"], 50)))
    say(pd.DataFrame(per_unit).round(2).to_string(index=False))

    say("\n== alpha_eff vs CBH (0.25 km bins) and vs month (nominal, no-edge, working) ==")
    vw2 = vw[vw.edge_nom == 0].copy()
    vw2["cbh_bin"] = (np.floor(vw2.cbh_m / 250.0) * 0.25).round(2)
    vw2["month"] = (vw2.date // 100).astype(int)
    tabs = []
    for gname, sel in GROUPS:
        g = vw2[sel(vw2)]
        for by in ("cbh_bin", "month"):
            t = g.groupby(by).alpha_nom.agg(n="size", med="median",
                                            p10=lambda x: q(x, 10), p90=lambda x: q(x, 90)).reset_index()
            t.insert(0, "by", by)
            t.insert(0, "group", gname)
            tabs.append(t.rename(columns={by: "bin"}))
    cb = pd.concat(tabs, ignore_index=True)
    cb.to_csv(OUT / "alpha_vs_cbh_month.csv", index=False, float_format="%.4g")
    say(cb[(cb.by == "cbh_bin") & (cb.group != "Vaisala excl. Payerne CL31")].pivot_table(
        index="bin", columns="group", values="med").round(1).to_string())

    say("\n== Cloudnet closure (alpha_shape nominal vs Cloudnet alpha, lowest 100 m, status 1, +-15 min) ==")
    cl = vw[np.isfinite(vw.alpha_cn_km)]
    crow = []
    for (k, it, site), g in cl.groupby(["key", "itype", "site"]):
        for lab, gg in (("all", g), ("no-edge", g[g.edge_nom == 0])):
            if len(gg) < 10:
                continue
            rs = spearmanr(gg.alpha_nom, gg.alpha_cn_km).statistic
            dm = gg.groupby("date")[["alpha_nom", "alpha_cn_km"]].median()
            rsd = spearmanr(dm.alpha_nom, dm.alpha_cn_km).statistic if len(dm) >= 5 else np.nan
            rf = spearmanr(gg.F100, gg.alpha_cn_km).statistic
            rk = spearmanr(gg.k_app, gg.alpha_cn_km, nan_policy="omit").statistic
            dk = gg.groupby("date")[["alpha_k", "alpha_cn_km"]].median()
            rkd = spearmanr(dk.alpha_k, dk.alpha_cn_km, nan_policy="omit").statistic if len(dk) >= 5 else np.nan
            crow.append(dict(unit=f"{it} {site}", key=k, sel=lab, n=len(gg), n_days=len(dm),
                             spearman_profile=rs, spearman_daily=rsd, spearman_F100_vs_cn=rf,
                             spearman_kapp_vs_cn=rk, spearman_alpha_k_daily=rkd,
                             alpha_shape_med=q(gg.alpha_nom, 50), alpha_k_med=q(gg.alpha_k, 50),
                             alpha_cn_med=q(gg.alpha_cn_km, 50),
                             median_ratio_shape_cn=q(gg.alpha_nom / gg.alpha_cn_km, 50),
                             median_ratio_alpha_k_cn=q(gg.alpha_k / gg.alpha_cn_km, 50),
                             aG_cn_med=q(gg.aG_cn_um, 50)))
    cs = pd.DataFrame(crow)
    cs.to_csv(OUT / "closure_cloudnet.csv", index=False, float_format="%.4g")
    say(cs.round(3).to_string(index=False))

    say("\n== FOV contrast at Payerne: CL31 (0.83 mrad, DISTORTED unit) vs CL61 (0.56 mrad) ==")
    a31 = vw[vw.key == "0-20000-0-06610_B"].sort_values("time_ns")
    a61 = vw[vw.key == "0-20000-0-06610_C"].sort_values("time_ns")
    if len(a31) and len(a61):
        m = pd.merge_asof(a61, a31, on="time_ns", direction="nearest", tolerance=60_000_000_000,
                          suffixes=("_61", "_31"))
        m = m[np.isfinite(m.cbh_m_31) & (np.abs(m.cbh_m_61 - m.cbh_m_31) <= 30.0)]
        say(f"pairs (+-1 min, |dCBH| <= 30 m): {len(m)} on {m.date_61.nunique()} days")
        if len(m):
            md31, md61 = models("CL31", "working"), models("CL61", "working")
            node = np.clip(np.round((m.cbh_m_61.to_numpy() - CBH_NODES[0]) / 50.0).astype(int), 0, CBH_NODES.size - 1)
            ai = np.abs(np.log(m.alpha_nom_61.to_numpy())[:, None] - np.log(ALPHA_FINE)[None, :]).argmin(1)
            p31, p61 = list(PULSES).index(110.0), list(PULSES).index(160.0)
            ka31, hh31 = model_kapp(md31["Mh"], md31["zg"])
            ka61, hh61 = model_kapp(md61["Mh"], md61["zg"])
            pred = dict(F50_31=md31["F50"][node, p31, 1, ai], F50_61=md61["F50"][node, p61, 1, ai],
                        F100_31=md31["F100"][node, p31, 1, ai], F100_61=md61["F100"][node, p61, 1, ai],
                        k_31=ka31[node, p31, 1, ai], k_61=ka61[node, p61, 1, ai])
            frow = []
            for nm, o31, o61, p_31, p_61 in (("F50", m.F50_31, m.F50_61, pred["F50_31"], pred["F50_61"]),
                                             ("F100", m.F100_31, m.F100_61, pred["F100_31"], pred["F100_61"]),
                                             ("k_app", m.k_app_31, m.k_app_61, pred["k_31"], pred["k_61"])):
                frow.append(dict(quantity=nm, obs_CL31_med=q(o31, 50), obs_CL61_med=q(o61, 50),
                                 obs_diff_31_minus_61_med=q(o31.to_numpy() - o61.to_numpy(), 50),
                                 pvc_diff_31_minus_61_med=q(p_31 - p_61, 50),
                                 pvc_CL31_med=q(p_31, 50), pvc_CL61_med=q(p_61, 50)))
            frow.append(dict(quantity="alpha_fit", obs_CL31_med=q(m.alpha_nom_31, 50), obs_CL61_med=q(m.alpha_nom_61, 50),
                             obs_diff_31_minus_61_med=q(m.alpha_nom_31.to_numpy() / m.alpha_nom_61.to_numpy(), 50)))
            fv = pd.DataFrame(frow)
            fv.to_csv(OUT / "fov_contrast_payerne.csv", index=False, float_format="%.4g")
            say(fv.round(4).to_string(index=False))
            say("(alpha_fit row: obs_diff column = median per-pair ratio CL31/CL61)")

    chk = OUT / "check_before_vs_after.csv"
    if chk.exists():
        c = pd.read_csv(chk)
        say("\n== beta_before vs beta_after (random profiles per unit) ==")
        say(f"profiles {int(c.n.sum())}; max |shape diff| (unit-sum) {c.max_abs_shape_diff.max():.2e}; "
            f"median dF100 {c.dF100_med.median():+.2e}; alpha ratio after/before: median "
            f"{c.alpha_ratio_after_before_med.median():.4f}, range {c.alpha_ratio_after_before_med.min():.4f}"
            f"..{c.alpha_ratio_after_before_med.max():.4f}")


def figures(df, say):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    df = add_nominal(df)
    v = df[df.valid]
    cols = {"CL31": "C0", "CL51": "C2", "CL61": "C3"}

    # (a) mean normalised profile per type + LUT curves
    fig, axes = plt.subplots(1, 4, figsize=(21, 6.2), sharey=True)
    for ax, (it, grid) in zip(axes, (("CL31", "working"), ("CL51", "working"), ("CL61", "working"), ("CL61", "native"))):
        md = models(it, grid)
        g = v[(v.itype == it) & (v.grid == grid)]
        if g.empty:
            continue
        P = []
        for k in g.key.unique():
            z = np.load(SRC / f"{k}.npz", allow_pickle=True)
            pre = "native_" if grid == "native" else ""
            if abs(float(z[f"{pre}dr_m"]) - md["dg"]) > 0.05:
                g = g[g.key != k]                                  # e.g. CL31 units at 20 m gates
                continue
            B = z[f"{pre}beta_before"].astype("float64")
            depth = z[f"{pre}range_rel_m"].astype("float64")
            kd = np.round(depth / md["dg"]).astype(int)
            cidx = np.array([int(np.nonzero(kd == kk)[0][0]) for kk in md["ks"]])
            W = B[:, cidx]
            sel = g[g.key == k]
            idx = np.searchsorted(z["time_ns"], sel.time_ns.to_numpy())
            W = W[idx]
            P.append(W / W.sum(1, keepdims=True))
        P = np.vstack(P)
        ax.plot(np.nanmean(P, 0), md["zg"], "k-", lw=2.4, label=f"observed mean (n={P.shape[0]})")
        ax.fill_betweenx(md["zg"], np.nanpercentile(P, 25, 0), np.nanpercentile(P, 75, 0), color="k", alpha=0.12,
                         label="observed P25-P75")
        cbh_med = float(np.nanmedian(g.cbh_m))
        c = int(np.argmin(np.abs(CBH_NODES - cbh_med)))
        pi = list(PULSES).index(NOMINAL[it])
        for al, cc in zip((5.0, 10.0, 20.0, 40.0), ("C0", "C1", "C2", "C3")):
            ai = int(np.argmin(np.abs(ALPHA_FINE - al)))
            ax.plot(md["Mh"][c, pi, 1, ai], md["zg"], "--", color=cc, lw=1.6,
                    label=f"PVC alpha={al:.0f} /km (a_G 3.3, {NOMINAL[it]:.0f} ns, CBH {CBH_NODES[c]/1e3:.2f} km)")
        ax.axhline(0, color="grey", lw=0.8)
        ax.set_title(f"{it} - {grid} gates ({md['dg']:.1f} m)")
        ax.set_xlabel("beta / sum over the window")
        ax.grid(alpha=0.3)
        ax.legend(fontsize=7)
    axes[0].set_ylabel("depth relative to the peak gate [m]")
    fig.suptitle("(a) Mean normalised in-cloud profile (beta before corrections) vs PVC forward model",
                 fontweight="bold")
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    fig.savefig(OUT / "fig_a_mean_profile_vs_lut.png", dpi=130)
    plt.close(fig)

    # (b) alpha_eff histograms per type
    fig, axes = plt.subplots(1, 4, figsize=(23, 5.2))
    bins = np.logspace(np.log10(3.0), np.log10(60.0), 31)
    for ax, (lab, sel) in zip(axes, (("full-shape fit, nominal pulse", (v.grid == "working")),
                                     ("full-shape fit, no pulse (0 ns)", (v.grid == "working")),
                                     ("tail slope alpha_k (response-independent)", (v.grid == "working")),
                                     ("CL61: native 4.8 m vs working 9.6 m", (v.itype == "CL61")))):
        g = v[sel]
        if "tail slope" in lab:
            for it in ("CL31", "CL51", "CL61"):
                x = g.loc[g.itype == it, "alpha_k"].dropna()
                ax.hist(x, bins=bins, weights=np.full(len(x), 1 / max(len(x), 1)), histtype="step", lw=2,
                        color=cols[it], label=f"{it} (med {x.median():.1f}, n={len(x)})")
        elif "no pulse" in lab:
            for it in ("CL31", "CL51", "CL61"):
                x = g.loc[g.itype == it, "alpha_p0_ag3.3"].dropna()
                ax.hist(x, bins=bins, weights=np.full(len(x), 1 / max(len(x), 1)), histtype="step", lw=2,
                        color=cols[it], label=f"{it} (med {x.median():.1f})")
        elif "CL61" in lab:
            for grid, ls in (("working", "-"), ("native", "--")):
                x = g.loc[g.grid == grid, "alpha_nom"].dropna()
                ax.hist(x, bins=bins, weights=np.full(len(x), 1 / max(len(x), 1)), histtype="step", lw=2,
                        ls=ls, color="C3", label=f"CL61 {grid} (med {x.median():.1f})")
        else:
            for it in ("CL31", "CL51", "CL61"):
                x = g.loc[g.itype == it, "alpha_nom"].dropna()
                ax.hist(x, bins=bins, weights=np.full(len(x), 1 / max(len(x), 1)), histtype="step", lw=2,
                        color=cols[it], label=f"{it} (med {x.median():.1f}, n={len(x)})")
            x = g.loc[(g.itype == "CL31") & ~g.distorted.astype(bool), "alpha_nom"].dropna()
            ax.hist(x, bins=bins, weights=np.full(len(x), 1 / max(len(x), 1)), histtype="step", lw=1.5, ls=":",
                    color="C0", label=f"CL31 excl. Payerne (med {x.median():.1f})")
        ax.set_xscale("log")
        ax.set_xlabel("alpha_eff [km$^{-1}$] (a_G 3.3 um; grid 3-60, ends = edge)")
        ax.set_ylabel("fraction of profiles per bin")
        ax.set_title(lab)
        ax.grid(alpha=0.3)
        ax.legend(fontsize=8)
    fig.suptitle("(b) alpha_eff histograms per instrument type", fontweight="bold")
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    fig.savefig(OUT / "fig_b_alpha_hist.png", dpi=130)
    plt.close(fig)

    # (c) alpha_shape vs alpha_Cloudnet
    cl = v[(v.grid == "working") & np.isfinite(v.alpha_cn_km)]
    keys = sorted(cl.key.unique())
    if keys:
        fig, axes = plt.subplots(2, len(keys), figsize=(5.2 * len(keys), 9.6))
        axes = np.atleast_2d(axes)
        from scipy.stats import spearmanr
        for row, (col, ylab) in enumerate((("alpha_nom", "alpha_shape, full-shape fit"),
                                           ("alpha_k", "alpha_k, tail slope"))):
            for ax, k in zip(axes[row], keys):
                g = cl[(cl.key == k) & np.isfinite(cl[col])]
                ax.scatter(g.alpha_cn_km, g[col], s=2, alpha=0.15, color="C0", rasterized=True)
                dm = g.groupby("date")[[col, "alpha_cn_km"]].median()
                ax.scatter(dm.alpha_cn_km, dm[col], s=14, color="C3", label="daily medians")
                rs = spearmanr(g[col], g.alpha_cn_km).statistic
                rsd = spearmanr(dm[col], dm.alpha_cn_km).statistic
                rr = np.nanmedian(g[col] / g.alpha_cn_km)
                ax.plot([0.3, 300], [0.3, 300], "k--", lw=1)
                ax.set_xscale("log")
                ax.set_yscale("log")
                ax.set_xlim(0.3, 300)
                ax.set_ylim(2, 80)
                ax.set_xlabel("Cloudnet alpha, lowest 100 m (radar+MWR) [km$^{-1}$]")
                ax.set_ylabel(f"{ylab} [km$^{{-1}}$]")
                g0 = g.iloc[0]
                ax.set_title(f"{g0.itype} {g0.site}{' (distorted)' if g0.distorted else ''}: n={len(g)}\n"
                             f"Spearman profile {rs:.2f} / daily {rsd:.2f}, median ratio {rr:.2f}", fontsize=9)
                ax.grid(alpha=0.3)
                ax.legend(fontsize=8)
        fig.suptitle("(c) Closure: ceilometer alpha (top: full-shape fit, bottom: tail slope) vs Cloudnet "
                     "alpha (+-15 min)", fontweight="bold")
        fig.tight_layout(rect=(0, 0, 1, 0.95))
        fig.savefig(OUT / "fig_c_closure_cloudnet.png", dpi=130)
        plt.close(fig)

    # (d) alpha_eff vs CBH (Y axis = CBH)
    fig, axes = plt.subplots(1, 2, figsize=(16, 6))
    g = v[(v.grid == "working") & (v.edge_nom == 0)].copy()
    g["cbh_bin"] = np.floor(g.cbh_m / 250.0) * 0.25 + 0.125
    for it in ("CL31", "CL51", "CL61"):
        for ax, sel, ls in ((axes[0], g.itype == it, "-"),
                            (axes[1], (g.itype == it) & ~g.distorted.astype(bool), "-")):
            t = g[sel].groupby("cbh_bin").alpha_nom
            st = t.agg(["median", "size"])
            st = st[st["size"] >= 100]
            p10 = t.quantile(0.1).reindex(st.index)
            p90 = t.quantile(0.9).reindex(st.index)
            ax.plot(st["median"], st.index, ls, color=cols[it], lw=2, marker="o", label=f"{it} full-shape fit")
            ax.fill_betweenx(st.index, p10, p90, color=cols[it], alpha=0.12)
            tk = g[sel].groupby("cbh_bin").alpha_k.median().reindex(st.index)
            ax.plot(tk, st.index, "--", color=cols[it], lw=2, marker="s", label=f"{it} tail slope alpha_k")
    for ax, t in zip(axes, ("all units", "excl. Payerne CL31 (distorted)")):
        ax.set_xscale("log")
        ax.set_xlabel("alpha [km$^{-1}$] (median; band = P10-P90 of the full-shape fit; a_G 3.3 um, nominal pulse)")
        ax.set_ylabel("cloud base height [km]")
        ax.set_title(t)
        ax.grid(alpha=0.3)
        ax.legend()
    fig.suptitle("(d) alpha_eff vs cloud base height (0.25 km bins, >= 100 profiles)", fontweight="bold")
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    fig.savefig(OUT / "fig_d_alpha_vs_cbh.png", dpi=130)
    plt.close(fig)
    say("figures (a)-(d) written")


def eta_preview(df, say):
    """Step-3 preview: eta(CBH) per type at a_G 3.3 and the median / P10 / P90 alpha_eff."""
    import pandas as pd
    from validation.multiple_scattering_eta import OPTICS, WAVELENGTH, eta_for_cloud_base
    from calibration.cloud._filters import (_ETA_CL31, _ETA_CL51, _ETA_CL61,
                                            _ETA_LEGACY_CL31, _ETA_LEGACY_CL51)
    df = add_nominal(df)
    v = df[df.valid & (df.grid == "working") & (df.edge_nom == 0)]
    oper = {"CL31": _ETA_CL31, "CL51": _ETA_CL51, "CL61": _ETA_CL61}
    hop = {"CL31": _ETA_LEGACY_CL31, "CL51": _ETA_LEGACY_CL51, "CL61": _ETA_LEGACY_CL51}
    rows = []
    for it in ("CL31", "CL51", "CL61"):
        sel = v.itype == it
        if it == "CL31":
            sel &= ~v.distorted.astype(bool)
        for est, col in (("full-shape fit", "alpha_nom"), ("tail slope alpha_k", "alpha_k")):
            x = v.loc[sel, col]
            if x.dropna().empty:
                continue
            a10, a50, a90 = q(x, 10), q(x, 50), q(x, 90)
            th = WAVELENGTH[it] / (np.pi * 3.3e-6)
            for zkm in (0.5, 1.0, 2.0):
                e = {lab: eta_for_cloud_base(zkm * 1000, *OPTICS[it], th, alpha=al * 1e-3)
                     for lab, al in (("p10", a10), ("med", a50), ("p90", a90))}
                rows.append(dict(itype=it, estimator=est, cbh_km=zkm, alpha_p10=a10, alpha_med=a50,
                                 alpha_p90=a90, eta_shape_med=e["med"], eta_at_alpha_p10=e["p10"],
                                 eta_at_alpha_p90=e["p90"],
                                 eta_operational=float(np.interp(zkm, oper[it][:, 0], oper[it][:, 1])),
                                 eta_hopkin_legacy=float(np.interp(zkm, hop[it][:, 0], hop[it][:, 1]))))
    t = pd.DataFrame(rows)
    t["ratio_shape_over_operational"] = t.eta_shape_med / t.eta_operational
    t.to_csv(OUT / "eta_preview_step3.csv", index=False, float_format="%.4g")
    say("\n== Step-3 preview: eta at a_G 3.3 um and the shape-derived alpha (CL31 excl. Payerne) ==")
    say(t.round(4).to_string(index=False))


if __name__ == "__main__":
    main()
