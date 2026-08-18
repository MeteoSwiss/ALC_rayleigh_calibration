"""Derive the APPLIED (onboard/manufacturer) overlap O_applied(z) from CLEAR-SKY NIGHT rcs_0 noise,
WITHOUT any hood measurement.

Physics (P view, P = rcs_0 / z^2 -> removes the range-correction; the firmware already divided by
O_applied so P carries a 1/O_applied factor in its NOISE):
  detector/electronic noise  sigma_elec_P(z) = sigma_elec_raw / O_applied(z)      (signal-INDEPENDENT)
  photon shot noise          sigma_shot_P(z)^2 = k_s * P(z) / O_applied(z)         (∝ signal)
  turbulent atmos. fluct.    sigma_turb_P(z)   ∝ P(z)                              (∝ signal)
Only the electronic term is INDEPENDENT of the signal P. Hence, at a fixed gate z, regressing the
white (successive-difference) variance against the mean signal,
      sigma_white^2(z) = beta(z) + alpha(z)*P + gamma(z)*P^2
isolates the electronic variance as the INTERCEPT beta(z)=sigma_elec_P^2 (shot -> alpha, turbulence
-> gamma). Then, since sigma_elec_P = sigma_elec_raw/O:
      O_applied(z) = sigma_elec_raw / sqrt(beta(z))   normalised so O -> 1 in the full-overlap far band.

Three estimators (as requested: our electronic-noise method + 2 adapted):
  M1  clean-night successive-difference floor  (assumes shot/turb negligible on the cleanest bins)
  M2  variance-vs-signal regression intercept  (rigorous shot/turbulence separation)  <-- principled
  M3  minimum-statistics envelope              (electronic floor = low quantile of sigma_white)
The three agree only if the noise model holds -> mutual validation.
"""
import numpy as np, netCDF4, glob, os, json, warnings, re
warnings.filterwarnings("ignore")
from scipy.signal import medfilt

L1 = "D:/E-PROFILE_L1_2026"
NORM_BAND = (1500.0, 3500.0)   # full-overlap band where O=1 and electronic noise dominates a clear night
CHM_GRIDS = (1024, 1535, 1536) # CHM15k range-gate configs (14.985 m for 1024; 9.99 m for 1535/1536) -
#                                excludes CL61(3276)/CL31(770)/CL51/Mini-MPL that share a station's WMO

def _orient(rcs, nr):
    return rcs.T if rcs.shape[0] == nr else rcs

def load_night(fp, h0=0, h1=6, cloud_ceiling=6000.0, wt_min=85.0):
    """Return dict with night, clear P(t,z)=rcs_0/z^2, time(min), rng, and a clearness score."""
    nc = netCDF4.Dataset(fp)
    t = np.asarray(nc['time'][:], float); rng = np.asarray(nc['range'][:], float)
    if len(rng) not in CHM_GRIDS:                                        # CHM15k only (skip CL61/CL31 under same WMO)
        nc.close(); return None
    rcs = _orient(np.asarray(nc['rcs_0'][:], float), len(rng))            # (time,range)
    cbh = np.asarray(nc['cloud_base_height'][:], float) if 'cloud_base_height' in nc.variables else None
    wt = np.asarray(nc['window_transmission'][:], float) if 'window_transmission' in nc.variables else None
    tres = float(nc['time_resol'][:]) if 'time_resol' in nc.variables else np.nan
    nc.close()
    import pandas as pd
    dt = pd.to_datetime(t, unit='D', origin='unix'); hh = dt.hour + dt.minute/60.0
    nm = np.asarray((hh >= h0) & (hh < h1))
    if nm.sum() < 60:
        return None
    P = rcs[nm, :] / (rng[None, :]**2)                                   # P = rcs_0 / z^2
    tmin = (t[nm] - np.floor(t[nm][0])) * 24 * 60                        # minutes in day
    # cloud / rain screen per profile
    if cbh is not None:
        cb = np.nanmin(np.where((cbh > 0) & (cbh < 20000), cbh, np.nan), axis=1) if cbh.ndim == 2 else cbh
        cb = cb[nm]
    else:
        cb = np.full(nm.sum(), np.nan)
    good = ~((cb > 0) & (cb < cloud_ceiling))                           # drop profiles with cloud below ceiling
    if wt is not None:                                                  # drop profiles with window transmission
        wtn = wt[nm]; med = np.nanmedian(wtn)                          # anomalously low vs the night's own norm
        if np.isfinite(med) and med > 0:                              # (rain/dew on the window) — instrument-agnostic
            good &= (wtn > 0.90 * med)
    if good.sum() < 60:
        return None
    P = P[good, :]; tmin = tmin[good]
    order = np.argsort(tmin); P = P[order]; tmin = tmin[order]
    clear_frac = float(good.mean())
    # atmospheric HOMOGENEITY (for the homogeneity filter), measured above the overlap where a clean,
    # stationary night should look molecular (smoothly decreasing, unchanging).
    #  vgrad = LAYEREDNESS = strongest signal-INCREASE-with-height of the range-corrected profile (an
    #          elevated aerosol layer; impossible for pure molecular which only decreases) -> catches layers.
    #  tchange = temporal drift (first vs last third of the night) -> catches MOVING layers (the real M2 killer).
    Pmed = np.nanmedian(P, axis=0)
    with np.errstate(all='ignore'):
        lnrc = np.log(np.clip(Pmed * rng ** 2, 1.0, None))              # range-corrected signal ln(rcs_0)
        g = np.gradient(medfilt(lnrc, 5), rng)
    vb = (rng >= 500) & (rng <= 2500)
    vgrad = float(np.nanmax(g[vb]) * 1000)                              # layeredness [per km]
    npf = P.shape[0]; a = np.nanmedian(P[:max(npf // 3, 1)], axis=0); bb = np.nanmedian(P[-max(npf // 3, 1):], axis=0)
    tb = (rng >= 500) & (rng <= 2500)
    tchange = float(np.nanmedian(np.abs(a[tb] - bb[tb]) / (0.5 * np.abs(a[tb] + bb[tb]) + 1e-30)))  # temporal drift
    return dict(P=P, tmin=tmin, rng=rng, tres=tres, clear_frac=clear_frac, vgrad=vgrad, tchange=tchange,
                nprof=P.shape[0], day=(re.search(r'(\d{8})(?=\.nc)', os.path.basename(fp)) or re.search(r'(\d{8})', os.path.basename(fp))).group(1))

def night_samples(night, bin_min=20.0, max_gap_factor=2.5):
    """Per (time-bin, gate): white-noise sigma from lag-1 successive differences and the mean signal.
    Returns arrays P_mean (nbin,ngate), sig_white (nbin,ngate)."""
    P = night['P']; tmin = night['tmin']; ng = P.shape[1]
    # lag-1 successive differences between temporally-adjacent profiles (skip gaps)
    dtp = np.diff(tmin); med = np.nanmedian(dtp[dtp > 0])
    okpair = dtp < max_gap_factor * med
    dP = (P[1:, :] - P[:-1, :])                                          # (npair,ngate)
    pmid = 0.5 * (tmin[1:] + tmin[:-1])
    dP = dP[okpair]; pmid = pmid[okpair]; Ppair = 0.5 * (P[1:, :] + P[:-1, :])[okpair]
    if dP.shape[0] < 20:
        return None
    edges = np.arange(pmid.min(), pmid.max() + bin_min, bin_min)
    idx = np.digitize(pmid, edges)
    Pm, Sw = [], []
    for b in np.unique(idx):
        sel = idx == b
        if sel.sum() < 8:
            continue
        d = dP[sel, :]
        # robust white sigma from successive diffs: MAD/0.6745 / sqrt(2)
        madd = np.nanmedian(np.abs(d - np.nanmedian(d, axis=0)), axis=0) * 1.4826 / np.sqrt(2.0)
        Pm.append(np.nanmedian(Ppair[sel, :], axis=0)); Sw.append(madd)
    if not Pm:
        return None
    return np.array(Pm), np.array(Sw)                                   # (nbin,ngate) each

def _norm_far(rng, y, band=NORM_BAND):
    m = (rng >= band[0]) & (rng <= band[1])
    ref = np.nanmedian(y[m])
    return y / ref if ref and np.isfinite(ref) else y * np.nan

def derive(files, verbose=True, norm_band=NORM_BAND, min_use=8, vgrad_max=6.0, tchange_max=0.4):
    """Aggregate clear nights -> O_applied by M1/M2/M3 on the L1 range grid.
    HOMOGENEITY FILTER: prefer temporally/vertically homogeneous nights (low vgrad & tchange); drop the
    inhomogeneous ones, but always keep at least `min_use` (the most homogeneous available)."""
    rng = None; nights = []
    for fp in files:
        n = load_night(fp)
        if n is None:
            continue
        s = night_samples(n)
        if s is None:
            continue
        if rng is None:
            rng = n['rng']
        if len(n['rng']) != len(rng):
            continue
        nights.append((s[0], s[1], n['vgrad'], n['tchange']))
    if not nights:
        raise RuntimeError("no usable clear nights")
    vg = np.array([x[2] for x in nights]); tc = np.array([x[3] for x in nights])
    homog_ok = (vg <= vgrad_max) & (tc <= tchange_max)                    # homogeneity filter
    if homog_ok.sum() < min_use:                                          # fallback: keep the most homogeneous min_use
        score = vg / vgrad_max + tc / tchange_max
        homog_ok = np.zeros(len(nights), bool); homog_ok[np.argsort(score)[:min_use]] = True
    keep = [nights[i] for i in range(len(nights)) if homog_ok[i]]
    PM = np.vstack([k[0] for k in keep]); SW = np.vstack([k[1] for k in keep])   # (Nbin, ngate)
    ndays = len(keep); nbins = PM.shape[0]
    ng = SW.shape[1]
    # ---- M1: clean-night floor = median sigma_white over the cleanest tercile of bins per gate
    sig1 = np.full(ng, np.nan)
    for g in range(ng):
        p = PM[:, g]; s = SW[:, g]; ok = np.isfinite(p) & np.isfinite(s) & (s > 0)
        if ok.sum() < 8:
            continue
        thr = np.nanpercentile(p[ok], 33)
        sig1[g] = np.nanmedian(s[ok & (p <= thr)])
    O1 = _norm_far(rng, 1.0 / sig1, norm_band)
    # ---- M2: per-gate NON-NEGATIVE linear regression of white variance vs mean signal, on P-quantile
    # bins (robust): sigma_white^2 = beta + alpha*P.  Intercept beta = signal-INDEPENDENT electronic
    # variance (extrapolated to P->0, removing shot+turbulence). O ∝ 1/sqrt(beta).
    from scipy.optimize import nnls
    beta = np.full(ng, np.nan)
    for g in range(ng):
        p = PM[:, g]; s2 = SW[:, g]**2; ok = np.isfinite(p) & np.isfinite(s2) & (s2 > 0)
        if ok.sum() < 15:
            continue
        x = p[ok]; y = s2[ok]
        if not (np.nanmax(x) > 3 * np.nanmax([np.nanmin(x), 1e-30])):    # need signal dynamic range
            beta[g] = np.nanpercentile(y, 10); continue
        # bin by signal quantiles, median variance per bin (down-weights outliers)
        qe = np.nanpercentile(x, np.linspace(0, 100, 9)); qe[-1] += 1e-9
        bx, by = [], []
        for k in range(len(qe) - 1):
            m = (x >= qe[k]) & (x < qe[k+1])
            if m.sum() >= 3:
                bx.append(np.nanmedian(x[m])); by.append(np.nanmedian(y[m]))
        if len(bx) < 4:
            beta[g] = np.nanpercentile(y, 10); continue
        bx = np.array(bx); by = np.array(by); sc = np.nanmedian(bx) or 1.0
        xs = bx / sc
        A = np.vstack([np.ones_like(xs), xs, xs**2]).T                    # beta + alpha P + gamma P^2 (turb)
        try:
            coef, _ = nnls(A, by)                                        # all >= 0
            beta[g] = max(coef[0], 1e-30)
        except Exception:
            beta[g] = np.nanpercentile(y, 10)
    from scipy.signal import medfilt
    O2 = _norm_far(rng, 1.0 / np.sqrt(medfilt(np.nan_to_num(beta, nan=np.nanmedian(beta)), 5)), norm_band)
    # ---- M3: minimum-statistics envelope = 10th percentile of sigma_white per gate
    sig3 = np.full(ng, np.nan)
    for g in range(ng):
        s = SW[:, g]; ok = np.isfinite(s) & (s > 0)
        if ok.sum() < 8:
            continue
        sig3[g] = np.nanpercentile(s[ok], 10)
    O3 = _norm_far(rng, 1.0 / sig3, norm_band)
    if verbose:
        print(f"  derived from {ndays} clear nights, {nbins} time-bins")
    return dict(rng=rng, O1=O1, O2=O2, O3=O3, sig1=sig1, beta=beta, sig3=sig3,
                ndays=ndays, nbins=nbins, PM=PM, SW=SW)

def select_clear_nights(wmo, months, max_nights=25, min_clear=0.80):
    """Scan L1 files, rank by night clear-fraction, return the clearest file paths."""
    cand = []
    for ym in months:
        y, m = ym.split('-')
        for fp in sorted(glob.glob(f"{L1}/{wmo}/{y}/{m}/L1_*.nc")):
            try:
                nc = netCDF4.Dataset(fp)
                if len(nc['range']) not in CHM_GRIDS:                    # CHM15k only
                    nc.close(); continue
                t = np.asarray(nc['time'][:], float)
                cbh = np.asarray(nc['cloud_base_height'][:], float) if 'cloud_base_height' in nc.variables else None
                nc.close()
                import pandas as pd
                hh = pd.to_datetime(t, unit='D', origin='unix'); hh = hh.hour + hh.minute/60.0
                nm = (hh >= 0) & (hh < 6)
                if nm.sum() < 100 or cbh is None:
                    continue
                cb = np.nanmin(np.where((cbh > 0) & (cbh < 20000), cbh, np.nan), axis=1) if cbh.ndim == 2 else cbh
                clear = np.mean(~((cb[nm] > 0) & (cb[nm] < 6000)))
                if clear >= min_clear:
                    cand.append((clear, fp))
            except Exception:
                continue
            if len(cand) >= 4 * max_nights:                             # early stop once plenty of clear nights
                break
        else:
            continue
        break
    cand.sort(reverse=True)
    return [fp for _, fp in cand[:max_nights]]

if __name__ == "__main__":
    S = os.path.dirname(__file__)
    # Payerne TUB140016 feasibility: winter+spring clearest nights (low aerosol -> electronic floor clean)
    months = [f"2025-{m:02d}" for m in [4,5,6,7,8,9,10,11,12]] + [f"2026-{m:02d}" for m in [1,2,3,4,5,6]]
    print("selecting clear Payerne nights...")
    files = select_clear_nights("0-20000-0-06610", months, max_nights=25, min_clear=0.85)
    print(f"  {len(files)} clear nights selected")
    R = derive(files)
    np.savez(S + "/clearsky_payerne.npz", **{k: R[k] for k in ['rng','O1','O2','O3','sig1','beta','sig3','PM','SW']},
             ndays=R['ndays'], nbins=R['nbins'])
    # compare to hood
    hood = np.load(S + "/tub140016_overlap_noise.npz"); Oh = hood['O_applied']; rh = hood['rng']
    print(f"\n{'h[m]':>6} | {'O1 clean':>8} {'O2 regr':>8} {'O3 min':>8} | {'O hood':>8}")
    for h in [150,225,300,450,600,750,900,1050,1200,1500]:
        f = lambda O: float(np.interp(h, R['rng'], O))
        print(f"{h:>6} | {f(R['O1']):>8.3f} {f(R['O2']):>8.3f} {f(R['O3']):>8.3f} | {float(np.interp(h,rh,Oh)):>8.3f}")
