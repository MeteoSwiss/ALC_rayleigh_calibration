# -*- coding: utf-8 -*-
"""ref23 step 3: the per-night pair observable ln rho and its predicted uncertainty.

Observable (per night n and pair X / R)
  rho(n) = median over the night's admitted pixels (common clear hours x gates of the band) of
           beta_X(h, z) / beta_R(h, z),
  beta = rcs_0 x 1e6 / C_L(n)  (each instrument divided by ITS OWN constant of night n), divided by
  the two-way water-vapour transmission of the 910 nm units (their own laser line), and -- when a
  910 nm unit is compared with a 1064 nm one -- converted to 1064 nm with the v3 dashboard's
  component-separated rule beta_1064 = beta_mol,1064 + (beta_910 - beta_mol,910) x (910/1064)^alpha
  (CAMS T/p molecular attenuated backscatter; alpha = Angstrom exponent, nominal 1).
  Pixels: hourly medians on the v3 display grid, hours of the calibration night (SZA > 100 deg)
  clear for EVERY co-located unit, and both units detected at SNR >= 3 (v3 "common" rule, 5-min
  windows, hour admitted when >= 50 % of its data are detected).  Bands 500-3000 m (v3) and
  1500-3000 m (above full overlap).

Calibration states of a Rayleigh channel
  v22pub  the published constant (calibration_version 220) and its published uncertainty
          (2 x robust sigma of the 125-member ensemble)
  v22     the v2.2 rerun (identical code path; reproduces v22pub) with the ensemble uncertainty
  v22b    the v2.2 constant with the v2.3 budget evaluated on its window (isolates the effect of
          the uncertainty model from that of the repaired tier 2)
  v23     the release candidate: eprof_v2.3 constant + v2.3 budget
  The budget totals are recomputed here from the stored 1-sigma components.
Cloud channels: see cloud_reference() (same-night day value, else bridged).

Usage: python validation/v23_reference/ref23_pairs.py   (needs ref23_tables.py first)
"""
from __future__ import annotations

import pickle
import sys
import warnings
from datetime import datetime, timedelta
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import ref23_common as RC  # noqa: E402

TAB = RC.OUT / "tables"
RAY_STATES = ("v22pub", "v22", "v23", "v23pre")  # distinct night sets (v22b shares v22's)
COMPONENTS = ("sig_ens", "sig_vert", "sig_mol", "sig_aer", "sig_noise")
# Payerne CL31: optical block replaced 2026-07-07 ~13:00 (cloud constant x2.2 over the 30 days on
# either side, ref23_report.cl31_step); the day itself mixes both.
CL31_STEP_DAY = "20260707"


# ------------------------------------------------------------------------------- constants
def load_tables():
    """rays = the REVIEWED v2.3 budget (main); rays_prereview.pkl is attached as rays.attrs['pre']."""
    rays = pd.read_pickle(TAB / "rays.pkl")
    pre = TAB / "rays_prereview.pkl"
    if pre.exists():
        rays.attrs["pre"] = pd.read_pickle(pre)
    clouds = pd.read_pickle(TAB / "clouds.pkl")
    pub = pd.read_pickle(TAB / "published.pkl")
    with open(TAB / "ingredients.pkl", "rb") as f:
        ing = pickle.load(f)
    return rays, clouds, pub, ing


def budget_sigma(row, drop=()):
    """1-sigma relative total of the v2.3 budget from its stored components (drop = common terms)."""
    v = [row[c] for c in COMPONENTS if c not in drop]
    v = np.asarray(v, float)
    if not np.all(np.isfinite(v)):
        return np.nan
    return float(np.sqrt(np.sum(v ** 2)))


def ray_states(rays, pub):
    """{(site, ident, date): {state: dict(lnc, sig, sig_nomol, tier, ...)}} for every Rayleigh night.

    States: v22pub (published), v22 (rerun, ensemble sigma), v22b (v2.2 constant + reviewed v2.3
    budget), v23 (eprof_v2.3 + reviewed budget and its QC) and, from the first pass,
    v23pre (eprof_v2.3 + the PRE-REVIEW budget and its QC: same constants, other sigma and nights)."""
    out = {}
    pre = rays.attrs.get("pre")
    if pre is not None:
        for row in pre[pre.variant == "v23"].itertuples():
            if row.flag in RC.SUCCESS and np.isfinite(row.c) and row.c > 0:
                r = row._asdict()
                st = out.setdefault((row.site, row.ident, row.date), {})
                st["v23pre"] = dict(lnc=np.log(row.c), sig=budget_sigma(r),
                                    sig_nomol=budget_sigma(r, drop=("sig_mol",)), tier=row.tier,
                                    zb=row.zb, zt=row.zt)
    rr = rays.set_index(["site", "ident", "date", "variant"]).sort_index()
    for (site, ident, ds, v), row in rr.iterrows():
        key = (site, ident, ds)
        st = out.setdefault(key, {})
        ok = row["flag"] in RC.SUCCESS and np.isfinite(row["c"]) and row["c"] > 0
        comps = {c: row[c] for c in COMPONENTS}
        if v == "v22":
            # the v2.2 'u > C' gate (-6) that the recording configuration does not apply
            if ok and np.isfinite(row["u_ens"]) and row["u_ens"] > row["c"]:
                ok = False
            if ok:
                st["v22"] = dict(lnc=np.log(row["c"]), sig=row["sig_ens"], sig_nomol=row["sig_ens"],
                                 tier=row["tier"], zb=row["zb"], zt=row["zt"], **comps)
                st["v22b"] = dict(lnc=np.log(row["c"]), sig=budget_sigma(row),
                                  sig_nomol=budget_sigma(row, drop=("sig_mol",)),
                                  tier=row["tier"], zb=row["zb"], zt=row["zt"], **comps)
            st["v22_flag"] = row["flag"]
            st["v22_tier"] = row["tier"]
        elif v == "v23":
            if ok:
                st["v23"] = dict(lnc=np.log(row["c"]), sig=budget_sigma(row),
                                 sig_nomol=budget_sigma(row, drop=("sig_mol",)),
                                 tier=row["tier"], zb=row["zb"], zt=row["zt"], g=row["g"], **comps)
            st["v23_flag"] = row["flag"]
            st["v23_tier"] = row["tier"]
    pr = pub[pub.method == "ray"]
    for row in pr.itertuples():
        if not (np.isfinite(row.c) and row.c > 0):
            continue
        st = out.setdefault((row.site, row.ident, row.date), {})
        s = row.u / (2.0 * row.c) if np.isfinite(row.u) else np.nan
        st["v22pub"] = dict(lnc=np.log(row.c), sig=s, sig_nomol=s,
                            tier=st.get("v22", {}).get("tier", np.nan))
    return out


# ------------------------------------------------------------------------------- cloud reference
def cloud_days(clouds, site, ident):
    """Successful cloud days of one channel -> DataFrame indexed by date (YYYYMMDD) with lnc and
    the candidate 1-sigma uncertainties of the day value."""
    g = clouds[(clouds.site == site) & (clouds.ident == ident)].copy()
    g = g[g.flag.isin(RC.SUCCESS) & np.isfinite(g.c) & (g.c > 0)]
    if site == "payerne" and ident == "B":
        g = g[g.date != CL31_STEP_DAY]
    g["lnc"] = np.log(g.c)
    g["sig_rep"] = g["std"] / g["coef"]                     # what is published (1 sigma spread)
    nsc = np.maximum(g["n_scenes"].fillna(1).values, 1)
    # standard error of the day's median coefficient, counting independent scenes, not profiles
    g["sig_sem"] = 1.2533 * g["p_mad_rel"] / np.sqrt(nsc)
    g["era"] = 0
    if site == "payerne" and ident == "B":
        g.loc[g.date > CL31_STEP_DAY, "era"] = 1
    return g.set_index("date").sort_index()


def _dnum(ds):
    return (datetime.strptime(ds, "%Y%m%d") - RC.EPOCH).days


BRIDGE_DAYS = 7          # how far a cloud value may be carried to a night (each side)


def _day_str(dn):
    return (RC.EPOCH + timedelta(days=int(dn))).strftime("%Y%m%d")


def cloud_reference(cd, ds, sig_model, max_bridge_days=BRIDGE_DAYS, era_of=None, mode="nearest"):
    """Cloud constant of night ds (the night D-1 -> D) and its 1-sigma uncertainty.

    mode "nearest" (primary, as the brief asks: the cloud value of the night, bridged if missing)
      same  : the cloud calibration of UTC day D (the night's morning) and/or D-1 (its evening),
              averaged in log space when both exist (sigma / sqrt(k));
      bridge: otherwise the nearest successful day(s) within +-max_bridge_days (same era), averaged
              in log space at the minimal distance L (days from the night's days); sigma is the day
              sigma inflated by the drift of the channel's cloud constant over L days,
              sqrt(sigma_day^2 + q L), q from its structure function (ref23_cloud).
    mode "local" (sensitivity): the median of ln C over every successful day within
      +-max_bridge_days of the night (same era, >= 2 days); sigma = 1.2533 sigma_day / sqrt(n)
      combined with the drift over the mean |lag| of the days used.
    Returns (lnc, sigma, kind, lag_days) or None.
    """
    d = _dnum(ds)
    era = era_of(ds) if era_of is not None else None

    def usable(xs):
        xs = [x for x in xs if x in cd.index]
        if era is not None:
            xs = [x for x in xs if cd.loc[x, "era"] == era]
        return xs

    if mode == "local":
        cand = usable([_day_str(d + k) for k in range(-max_bridge_days - 1, max_bridge_days + 1)])
        if len(cand) < 2:
            return None
        v = cd.loc[cand]
        lags = np.array([max(0, (d - 1) - _dnum(x), _dnum(x) - d) for x in cand], float)
        s_day = sig_model["day"](v)
        s = float(np.hypot(1.2533 * s_day / np.sqrt(len(cand)), sig_model["drift"](lags.mean())))
        return float(np.median(v.lnc)), s, "local", float(lags.mean())
    same = usable([ds, _day_str(d - 1)])
    if same:
        v = cd.loc[same]
        s = sig_model["day"](v) / np.sqrt(len(same))
        return float(v.lnc.mean()), float(s), "same", 0
    for lag in range(1, max_bridge_days + 1):
        # lag L: day D-1-L (before the evening) or day D+L (after the morning)
        cand = usable([_day_str(d - 1 - lag), _day_str(d + lag)])
        if cand:
            v = cd.loc[cand]
            s_day = sig_model["day"](v) / np.sqrt(len(cand))
            return (float(v.lnc.mean()), float(np.hypot(s_day, sig_model["drift"](lag))),
                    "bridge", lag)
    return None


# ------------------------------------------------------------------------------- the ratio
# Measured covered-telescope dark b(z) of the Payerne units (rcs_0 units; the v3 dashboard's file).
DARK_NPZ = Path("C:/DATA/Projects/202606_E-PROFILE_calibration/rayleigh_availability/"
                "dark_profiles_payerne.npz")
_DARK = {}


def payerne_dark(ident):
    """b(z) of a Payerne unit on the display grid (rcs_0 units), as build_l1_l2_dashboard.dark_for."""
    if ident not in _DARK:
        prof = None
        if DARK_NPZ.exists():
            with np.load(DARK_NPZ) as z:
                if f"{ident}_b_rcs" in z:
                    rd = np.asarray(z[f"{ident}_range"], "f8")
                    bd = np.asarray(z[f"{ident}_b_rcs"], "f8")
                    ok = np.isfinite(rd) & np.isfinite(bd)
                    prof = np.nan_to_num(np.interp(RC.Z_AGL, rd[ok], bd[ok], left=np.nan,
                                                   right=np.nan), nan=0.0)
        _DARK[ident] = prof
    return _DARK[ident]


def _beta(ingn, ident, lnc, convert, alpha, dark=None):
    v = ingn["inst"].get(ident)
    if v is None:
        return None
    H = v["H"].astype("f8")
    if dark is not None:
        H = H - np.asarray(dark, "f8")[None, :]      # a constant offset commutes with the median
    B = H * 1e6 / np.exp(lnc)
    is910 = abs(v["wl"] - RC.TARGET_NM) >= 1.0
    if is910:
        c = ingn["corr"].get(ident)
        if c is None or c.get("t2") is None:
            return None                           # 910 nm without its WV correction: excluded
        B = B / np.asarray(c["t2"], "f8")[None, :]
        if convert:
            f = (v["wl"] / RC.TARGET_NM) ** alpha
            B = np.asarray(c["bmt"], "f8")[None, :] + (B - np.asarray(c["bml"], "f8")[None, :]) * f
    return B


def needs_conversion(ingn, x, r):
    wx, wr = ingn["inst"][x]["wl"], ingn["inst"][r]["wl"]
    return abs(wx - wr) >= 1.0


CLEAN_RMAX = 1.5      # "clean air": gates whose night-median scattering ratio at 1064 nm is < 1.5


def pair_values(ingn, x, r, lnc_x, lnc_r, alpha=RC.ALPHA_NOMINAL, W=300, band="b500_3000",
                profile=False, clean=False, dark_x=None, dark_r=None):
    """ln rho of one night (and optionally the ln rho(z) profile).

    clean=True keeps only the gates whose night-median scattering ratio at 1064 nm (mean of the
    two units, beta / beta_mol,1064 from CAMS) is below CLEAN_RMAX -- where the 910 -> 1064 nm
    conversion hardly depends on the Angstrom exponent (mixed-wavelength pairs only).
    dark_x / dark_r: a measured dark b(z) subtracted from that unit's profile (Payerne)."""
    if x not in ingn["inst"] or r not in ingn["inst"]:
        return None
    key = f"{x}|{r}"
    if key not in ingn["pairs"]:
        return None
    conv = needs_conversion(ingn, x, r)
    if clean and not conv:
        return None
    Bx = _beta(ingn, x, lnc_x, conv, alpha, dark_x)
    Br = _beta(ingn, r, lnc_r, conv, alpha, dark_r)
    if Bx is None or Br is None:
        return None
    adm = ingn["pairs"][key][f"adm{W}"]
    with np.errstate(all="ignore"):
        R = Bx / Br
    ok = adm & np.isfinite(R) & (Bx > 0) & (Br > 0)
    z0, z1 = RC.BANDS[band]
    zb = (RC.Z_AGL >= z0) & (RC.Z_AGL <= z1)
    if clean:
        c910 = ingn["corr"].get(x if ingn["inst"][x]["wl"] < RC.TARGET_NM - 1 else r)
        if c910 is None or c910.get("bmt") is None:
            return None
        bmt = np.asarray(c910["bmt"], "f8")[None, :]
        with np.errstate(all="ignore"), warnings.catch_warnings():
            warnings.simplefilter("ignore")
            sr = 0.5 * (np.nanmedian(np.where(ok, Bx / bmt, np.nan), axis=0)
                        + np.nanmedian(np.where(ok, Br / bmt, np.nan), axis=0))
        zb = zb & np.isfinite(sr) & (sr < CLEAN_RMAX)
    m = ok & zb[None, :]
    out = dict(n_pix=int(m.sum()), n_hours=int(np.any(m, axis=1).sum()),
               lnr=(float(np.log(np.median(R[m]))) if m.sum() >= 10 else np.nan))
    if profile:
        with np.errstate(all="ignore"), warnings.catch_warnings():
            warnings.simplefilter("ignore")
            P = np.where(ok, R, np.nan)
            out["prof"] = np.log(np.nanmedian(P, axis=0)).astype("f4")
            out["prof_n"] = np.sum(ok, axis=0).astype("i2")
    return out


# ------------------------------------------------------------------------------- pair table
def pair_kind(site, x, mx, r, mr):
    itx = RC.SITES[site]["streams"][x][0]
    itr = RC.SITES[site]["streams"][r][0]
    if x == r:
        return "self"                              # same instrument, Rayleigh vs cloud
    if mx == "ray" and mr == "ray":
        return "RR_twin" if itx == itr else "RR_mixed"
    if mx == "ray" and mr == "cloud":
        return "RC"
    if mx == "cloud" and mr == "cloud":
        return "CC"
    return "other"


def build_pair_table(rays, clouds, pub, ing, sig_model, with_profiles=True, alt_models=None):
    """alt_models: {name: sig_model} -- alternative cloud-uncertainty models; their sigma of the
    cloud reference is stored as cr_sig_<name> / cx_sig_<name> (the constants are identical)."""
    alt_models = alt_models or {}
    rs = ray_states(rays, pub)
    rows, profs = [], {}
    cds = {}
    for site, plist in RC.PAIRS.items():
        nights = sorted({ds for (s, ds) in ing if s == site} |
                        {k[2] for k in rs if k[0] == site})
        for (x, mx), (r, mr) in plist:
            kind = pair_kind(site, x, mx, r, mr)
            for ident in {x, r}:
                if (site, ident) not in cds:
                    cds[(site, ident)] = cloud_days(clouds, site, ident)
            era_of = None
            if site == "payerne" and "B" in (x, r) and "cloud" in (mx, mr):
                era_of = (lambda ds_: 0 if ds_ < CL31_STEP_DAY else 1)
            for ds in nights:
                if era_of is not None and ds == CL31_STEP_DAY:
                    continue                          # the CL31 swap day: mixed eras
                row = dict(site=site, pair=f"{x}:{mx}/{r}:{mr}", x=x, mx=mx, r=r, mr=mr, kind=kind,
                           date=ds, month=int(ds[4:6]), season=RC.season_of(ds[4:6]))
                # --- constants of the night -------------------------------------------------
                sx = rs.get((site, x, ds), {})
                sr = rs.get((site, r, ds), {})
                cref = {}                                  # (tag, mode) -> reference tuple
                for tag, ident, mm in (("cx", x, mx), ("cr", r, mr)):
                    if mm == "cloud" and (site, ident) in sig_model:
                        # the era split only concerns the Payerne CL31 itself
                        eo = era_of if (site == "payerne" and ident == "B") else None
                        for mode in ("nearest", "local"):
                            cref[(tag, mode)] = cloud_reference(
                                cds[(site, ident)], ds, sig_model[(site, ident)], era_of=eo,
                                mode=mode)
                cref_x, cref_r = cref.get(("cx", "nearest")), cref.get(("cr", "nearest"))
                loc_x, loc_r = cref.get(("cx", "local")), cref.get(("cr", "local"))
                # (state, ln C_X, ln C_R, suffix): suffix "" = nearest cloud reference (primary),
                # "_cl" = local (+-7 d median) cloud reference (sensitivity)
                states = []
                if kind in ("RC", "self"):
                    for s_ in RAY_STATES:
                        if s_ in sx:
                            if cref_r is not None:
                                states.append((s_, sx[s_]["lnc"], cref_r[0], ""))
                            if loc_r is not None:
                                states.append((s_, sx[s_]["lnc"], loc_r[0], "_cl"))
                elif kind in ("RR_twin", "RR_mixed"):
                    for s_ in RAY_STATES:
                        if s_ in sx and s_ in sr:
                            states.append((s_, sx[s_]["lnc"], sr[s_]["lnc"], ""))
                elif kind == "CC":
                    if cref_x is not None and cref_r is not None:
                        states.append(("cloud", cref_x[0], cref_r[0], ""))
                    if loc_x is not None and loc_r is not None:
                        states.append(("cloud", loc_x[0], loc_r[0], "_cl"))
                if not states:
                    continue
                ingn = ing.get((site, ds))
                for s_, lx, lr, ref in states:
                    if kind == "self":
                        # same photons: beta_X / beta_R = C_R / C_X
                        row[f"lnr_{s_}{ref}"] = lr - lx
                        continue
                    if ingn is None:
                        continue
                    for band in RC.BANDS:
                        bs = "" if band == "b500_3000" else "_hi"
                        alphas = RC.ALPHAS if ref == "" else (RC.ALPHA_NOMINAL,)
                        for a in alphas:
                            pv = pair_values(ingn, x, r, lx, lr, alpha=a, W=300, band=band,
                                             profile=(with_profiles and band == "b500_3000"
                                                      and a == RC.ALPHA_NOMINAL and ref == ""))
                            if pv is None:
                                continue
                            sfx = "" if a == RC.ALPHA_NOMINAL else f"_a{a:g}"
                            row[f"lnr_{s_}{bs}{sfx}{ref}"] = pv["lnr"]
                            if a == RC.ALPHA_NOMINAL and ref == "":
                                row[f"npix_{s_}{bs}"] = pv["n_pix"]
                                row[f"nh_{s_}{bs}"] = pv["n_hours"]
                                if "prof" in pv:
                                    profs[(site, row["pair"], ds, s_)] = (pv["prof"], pv["prof_n"])
                        if ref == "":
                            pv = pair_values(ingn, x, r, lx, lr, alpha=RC.ALPHA_NOMINAL, W=3600,
                                             band=band)
                            if pv is not None:
                                row[f"lnr_{s_}{bs}_w60"] = pv["lnr"]
                            if site == "payerne" and "cloud" in (mx, mr):
                                # sensitivity: the measured dark subtracted from the profiles of the
                                # cloud-calibrated units (their constants are dark-immune, so this is
                                # the coherent state of a cloud channel); Rayleigh units untouched
                                dkx = payerne_dark(x) if mx == "cloud" else None
                                dkr = payerne_dark(r) if mr == "cloud" else None
                                pv = pair_values(ingn, x, r, lx, lr, alpha=RC.ALPHA_NOMINAL, W=300,
                                                 band=band, dark_x=dkx, dark_r=dkr)
                                if pv is not None:
                                    row[f"lnr_{s_}{bs}_dk"] = pv["lnr"]
                            if kind in ("RC", "RR_mixed"):
                                for a in (0.5, RC.ALPHA_NOMINAL, 1.5):
                                    pv = pair_values(ingn, x, r, lx, lr, alpha=a, W=300,
                                                     band=band, clean=True)
                                    if pv is None:
                                        continue
                                    sfx = "" if a == RC.ALPHA_NOMINAL else f"_a{a:g}"
                                    row[f"lnr_{s_}{bs}_clean{sfx}"] = pv["lnr"]
                                    if a == RC.ALPHA_NOMINAL:
                                        row[f"npix_{s_}{bs}_clean"] = pv["n_pix"]
                if not any(k.startswith("lnr_") for k in row):
                    continue
                # --- uncertainties of the night (1 sigma, relative) ---------------------------
                for s_ in ("v22pub", "v22", "v22b", "v23", "v23pre"):
                    if s_ in sx:
                        row[f"sx_{s_}"] = sx[s_]["sig"]
                        row[f"sxnm_{s_}"] = sx[s_]["sig_nomol"]
                    if s_ in sr:
                        row[f"sr_{s_}"] = sr[s_]["sig"]
                        row[f"srnm_{s_}"] = sr[s_]["sig_nomol"]
                for c in COMPONENTS + ("g",):
                    if "v23" in sx and c in sx["v23"]:
                        row[f"x_{c}"] = sx["v23"][c]
                    if "v23" in sr and c in sr["v23"]:
                        row[f"r_{c}"] = sr["v23"][c]
                row["tier_x"] = sx.get("v23", sx.get("v22", {})).get("tier", np.nan)
                row["tier_r"] = sr.get("v23", sr.get("v22", {})).get("tier", np.nan)
                row["tier_x_v22"] = sx.get("v22", {}).get("tier", np.nan)
                row["tier_r_v22"] = sr.get("v22", {}).get("tier", np.nan)
                row["v23_flag_x"] = sx.get("v23_flag", np.nan)
                row["v23_flag_r"] = sr.get("v23_flag", np.nan)
                row["v22_flag_x"] = sx.get("v22_flag", np.nan)
                row["v22_flag_r"] = sr.get("v22_flag", np.nan)
                for s_ in ("v22pub", "v22", "v23"):
                    if s_ in sx:
                        row[f"lncx_{s_}"] = sx[s_]["lnc"]
                    if s_ in sr:
                        row[f"lncr_{s_}"] = sr[s_]["lnc"]
                for tag, cr_, lc_, ident in (("cx", cref_x, loc_x, x), ("cr", cref_r, loc_r, r)):
                    if cr_ is not None:
                        row[f"{tag}_sig"] = cr_[1]
                        row[f"{tag}_kind"] = cr_[2]
                        row[f"{tag}_lag"] = cr_[3]
                        row[f"{tag}_lnc"] = cr_[0]
                        eo = era_of if (site == "payerne" and ident == "B") else None
                        for name, sm in alt_models.items():
                            alt = (cloud_reference(cds[(site, ident)], ds, sm[(site, ident)],
                                                   era_of=eo)
                                   if (site, ident) in sm else None)
                            row[f"{tag}_sig_{name}"] = alt[1] if alt is not None else np.nan
                    if lc_ is not None:
                        row[f"{tag}_sig_cl"] = lc_[1]
                        row[f"{tag}_lag_cl"] = lc_[3]
                        row[f"{tag}_lnc_cl"] = lc_[0]
                rows.append(row)
    return pd.DataFrame(rows), profs
