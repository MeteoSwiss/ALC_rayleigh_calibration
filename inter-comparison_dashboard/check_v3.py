# -*- coding: utf-8 -*-
"""Non-regression + approximation check for the v3 payload.

Decodes v3_<site>.json exactly the way the browser does (base64 -> log-quantised uint16 -> the
affine transform chain) and recomputes the statistics table, then compares it to the v2 page's
precomputed numbers in data_<site>.json.  Any deviation here is a deviation the operator would see.

Also quantifies the approximation the v3 SPEC allowed but v3 does NOT use: dividing a window-median
profile by the window-median C_L instead of dividing per hour.

Run:  python inter-comparison_dashboard/check_v3.py [site] [combo]
"""
from __future__ import annotations
import base64
import json
import re
import sys
from pathlib import Path

sys.path.insert(0, r"C:/Users/hervo/OneDrive/Documents/ALC_rayleigh_calibration")
sys.path.insert(0, str(Path(__file__).resolve().parent))
import numpy as np

import warnings
warnings.filterwarnings("ignore")

from validation.paper import intercompare as IC
import build_l1_l2_dashboard as BD
import l1_l2_io as IO
import nf_v3 as NF
import sites
import variants_v3 as V3

DATA = sites.DATA_ROOT
TARGET = 1064.0

# Mechanical gate (2026-08-16 review: the check used to have no threshold, no exit code, and a
# headline metric normalised by a near-zero denominator).  Healthy worst on the gated sweep is a
# fraction of a point (quantisation jitter); a real perturbation of the wiring moves rows by
# several points, so these absolute tolerances separate the two cleanly.
TOL_PP = 0.5      # percentage points, on medrelbias_pct / relbias_pct
# r_log is fragile on the CL31 with both corrections off (near-zero in-band signal: positive-pair
# selection flips on tiny numerical differences; healthy worst 0.12 on that one corner while
# medrel matches to 0.001 pt) — the tolerance sits above that but far below a wiring swap.
TOL_R = 0.15

# L1 rows of the dark-run variants are INFORMATIONAL, not gated: the v2 snapshot predates the
# 2026-08-16 dark refit + window extension of diag_v22_dark, so their deviation measures data
# vintage, not wiring (their L2 rows use no constants and stay gated — they match exactly).
STALE_CALS = {"v2.2dark", "v2.2dark_est"}

# v2 swept ONE variant dimension across all channels; v3 deliberately dropped the redundant
# cloud x pipeline-version corner (the cloud method does not depend on the molecular method), so
# only these cloud mappings are expressible from the v3 payload — the rest is skipped with a note.
CAL_MAP_CLOUD = {"v2.2": "cloudWV", "v2.2sansWV": "cloudNoWV"}
# v2's estimated-dark variant kept its own id; v3 collapsed the naming (the page now discloses
# the estimated provenance via dark_kind instead)
CAL_MAP_RAYLEIGH = {"v2.2dark_est": "v2.2dark"}
# for a 1064 nm channel the WV correction is an exact no-op, so v2's 'sansWV' variant is the
# ordinary v2.2 run — v3 dropped that physically meaningless entry (operator request), map it back
CAL_MAP_1064 = {"v2.2sansWV": "v2.2"}


def dequant(enc):
    """base64 log-quantised uint16 -> float64 (n_h x n_z).  The JS `decodeStream` twin."""
    q = np.frombuffer(base64.b64decode(enc["b"]), "<u2").reshape(enc["nh"], enc["nz"]).astype("i4")
    lo = np.asarray(enc["lo"], "f8")[None, :]
    hi = np.asarray(enc["hi"], "f8")[None, :]
    sign = np.where(q >= 32768, -1.0, 1.0)
    qq = q & 32767
    with np.errstate(all="ignore"):
        mag = lo * np.power(hi / lo, (qq - 1) / 32766.0)
    out = sign * mag
    out[qq == 0] = 0.0
    out[q == 0] = np.nan
    return out


def defloat(blk):
    if blk is None:
        return None
    return np.frombuffer(base64.b64decode(blk["b"]), "<f4").reshape(blk["n"], blk["m"]).astype("f8")


def transform(raw, c, corr, wv, wl, hday, dark=None, calibrated=True):
    v = raw.copy()
    if dark is not None:
        v = v - np.asarray(dark, "f8")[None, :]
    if calibrated:
        v = v * 1e6 / np.asarray(c, "f8")[:, None]
    if wv and corr.get("wv") is not None:
        v = v / defloat(corr["wv"])[hday]
    need = corr.get("lam") is not None and abs(corr["lam"] - TARGET) >= 1.0
    if need and wl == "molecular" and corr.get("bmt") is not None:
        v = defloat(corr["bmt"])[hday] + (v - defloat(corr["bml"])[hday]) * corr["f"]
    elif need and wl == "angstrom":
        v = v * corr["f"]
    return v


def _combo_channels(site, cal, P):
    """v3 (ident, method, variant, dark_vec|None) list reproducing the v2 combo's channels, in the
    v2 channel order.  variant=None marks a channel the v3 payload cannot express for this cal
    (stats row skipped; the channel still feeds the paired mask via its canonical variant,
    which is exact because dividing by a positive constant never changes finiteness)."""
    base = {"payerne": [("A", "rayleigh"), ("B", "cloud"), ("C", "cloud"), ("C", "rayleigh")],
            "amsterdam": [(i, "rayleigh") for i in "ABCD"],
            "lindenberg": [("0", "rayleigh"), ("C", "cloud"), ("C", "rayleigh")]}[site]
    out = []
    for ident, method in base:
        if method == "rayleigh":
            variant = CAL_MAP_RAYLEIGH.get(cal, cal)
            corr = P["corr"][P["corr_of"][f"L1|{ident}"]]
            if corr.get("wv") is None:              # 1064 nm channel: WV variants collapse
                variant = CAL_MAP_1064.get(cal, variant)
        else:
            variant = CAL_MAP_CLOUD.get(cal)
        if variant is not None:
            rec = P["calib"].get(f"{ident}|{method}|{variant}")
            if not (rec and rec.get("ok")):
                variant = None
        # v2 subtracted the measured dark from every L1 channel under a dark variant
        db = P.get("dark", {}).get(ident) if cal in V3.DARK_RUNS else None
        out.append((ident, method, variant, db))
    return out


# --------------------------------------------------------------------------- noise-filter checks
def _debits(blk, dtype):
    return np.frombuffer(base64.b64decode(blk["b"]), dtype).reshape(blk["nh"], blk["nz"])


def check_nf(P, site):
    """Verify the shipped nf blocks against an INDEPENDENT recomputation from the native L1
    streams: the 60-min per-instrument masks against retime_hourly-style hourly bins +
    intercompare.snr_mask (the operational reference implementation), the intersection and
    scene bits by direct reconstruction, the p2 sums by direct sums, and a brute-force 5-min
    spot check on two sample days.  Returns the worst mismatch fraction (gated by the caller).
    Also PRINTS (informational) the size of the sub-hourly all-or-nothing approximation."""
    nf = P.get("nf")
    if not nf:
        print("\nnoise filter: no nf block in payload (rebuild needed) — SKIPPED, gate fails")
        return 1.0
    v3 = V3.SITE_V3[site]
    s = sites.get_site(site)
    BD.set_site(site)
    z_agl = np.asarray(P["z"], "f8")
    nz = z_agl.size
    salt = float(P["alt"])
    alt_asl = z_agl + salt
    hours = np.array([np.datetime64(h) for h in P["hours"]])
    hour_epoch = hours.astype("datetime64[h]").astype("i8")
    t0 = np.datetime64(P["start"])
    t1 = np.datetime64(P["end"]) + np.timedelta64(1, "D")
    d0, d1 = P["start"].replace("-", ""), P["end"].replace("-", "")
    idents = [i["ident"] for i in P["instruments"]]
    print(f"\nnoise filter: recomputing masks from native L1 ({len(idents)} streams) ...")

    worst = 0.0
    per_inst = {}          # ident -> dict(kfrac60, rn60, hpos) for the common/scene checks
    natives = {}
    for ident in idents:
        d = IO.read_l1_native(P["wmo"], ident, d0, d1)
        if d is None:
            print(f"  {ident}: NO native data — cannot verify")
            worst = 1.0
            continue
        tarr = np.asarray(d["time"])
        m = (tarr >= t0) & (tarr <= t1)
        scr = NF.nf_screen(np.asarray(d["beta"], "f4")[m], np.asarray(d["cbh"])[m],
                           np.asarray(d["vv"])[m], tarr[m])
        tsec = tarr[m].astype("datetime64[s]").astype("i8")
        o = np.argsort(tsec, kind="stable")
        tsec, scr = tsec[o], scr[o]
        alt_nat = np.asarray(d["alt"], "f8")
        dark_nat = BD.dark_for(ident, alt_nat - salt) if s.get("dark_npz") else None
        natives[ident] = (tsec, scr, alt_nat, dark_nat)

        # -- 60 min: hourly bins + the OPERATIONAL snr_mask, gate by gate ----------------------
        hid = tsec // 3600
        change = np.flatnonzero(np.diff(hid)) + 1
        st = np.concatenate([[0], change])
        en = np.concatenate([change, [hid.size]])
        keep = np.ones((st.size, alt_nat.size), bool)
        nfin = np.zeros((st.size, alt_nat.size))
        for j in range(st.size):
            X = np.asarray(scr[st[j]:en[j]], "f8")
            if dark_nat is not None:
                X = X - dark_nat[None, :]
            keep[j] = IC.snr_mask(X)                       # the reference implementation
            nfin[j] = np.isfinite(X).sum(axis=0)
        rk = IC.regrid(keep * nfin, alt_nat, alt_asl)
        rn = IC.regrid(nfin, alt_nat, alt_asl)
        with np.errstate(all="ignore"):
            kfrac = np.nan_to_num(rk) / np.where(np.nan_to_num(rn) > 0, np.nan_to_num(rn), 1.0)
        admit_ref = (np.nan_to_num(rn) > 0) & (kfrac >= NF.NF_FMIN)
        hpos = {int(h): j for j, h in enumerate(hid[st])}
        rows = np.array([hpos.get(int(h), -1) for h in hour_epoch])
        ok = rows >= 0
        mine = (_debits(nf["inst"][ident], "u1") >> 2 & 1).astype(bool)
        ref = np.zeros_like(mine)
        ref[ok] = admit_ref[rows[ok]]
        mism = float((mine[ok] != ref[ok]).mean()) if ok.any() else 1.0
        print(f"  {ident}: 60-min admit vs retime/snr_mask reference — mismatch "
              f"{mism * 100:.3f} % of {int(ok.sum())} h x {nz} gates")
        worst = max(worst, mism)
        per_inst[ident] = dict(kfrac=(np.nan_to_num(rn) > 0) & (kfrac >= NF.NF_FMIN),
                               pres=np.nan_to_num(rn) > 0, rn=np.nan_to_num(rn),
                               med=IC.regrid(np.vstack([np.nanmedian(
                                   np.asarray(scr[st[j]:en[j]], "f8"), axis=0)
                                   for j in range(st.size)]), alt_nat, alt_asl),
                               sig=IC.regrid(np.vstack([
                                   NF.robust_std(np.asarray(scr[st[j]:en[j]], "f8"), 0)
                                   for j in range(st.size)]), alt_nat, alt_asl),
                               rows=rows, okrows=ok)

    # -- common (intersection), 60 min ---------------------------------------------------------
    if all(i in per_inst for i in idents):
        both = np.ones((hour_epoch.size, nz), bool)
        seen = np.zeros(hour_epoch.size, bool) | True
        for ident in idents:
            pi = per_inst[ident]
            contrib = np.zeros((hour_epoch.size, nz), bool)
            okr = pi["okrows"]
            contrib[okr] = (pi["pres"] & pi["kfrac"])[pi["rows"][okr]]
            both &= contrib
        minec = (_debits(nf["common"], "u1") >> 2 & 1).astype(bool)
        mism = float((minec != both).mean())
        print(f"  common: 60-min intersection reconstruction — mismatch {mism * 100:.3f} %")
        worst = max(worst, mism)

    # -- scene, 60 min, all three thresholds ---------------------------------------------------
    if nf.get("scene") and nf.get("ref") in per_inst:
        rid = nf["ref"]
        pr = per_inst[rid]
        method, variant = v3["default"][rid]
        rec = P["calib"][f"{rid}|{method}|{variant}"]
        kd = np.array([np.datetime64(x) for x in rec["kal"]["d"]])
        cref = IC.interp_calib(kd, np.asarray(rec["kal"]["v"], "f8"),
                               hours + np.timedelta64(30, "m"))
        dref = np.asarray(P["dark"][rid], "f8") if P["dark"].get(rid) else np.zeros(nz)
        okr = pr["okrows"]
        med_h = np.full((hour_epoch.size, nz), np.nan)
        sig_h = np.zeros((hour_epoch.size, nz))
        n_h = np.zeros((hour_epoch.size, nz))
        med_h[okr] = pr["med"][pr["rows"][okr]]
        sig_h[okr] = np.nan_to_num(pr["sig"][pr["rows"][okr]])
        n_h[okr] = pr["rn"][pr["rows"][okr]]
        with np.errstate(all="ignore"):
            beta_cal = (np.nan_to_num(med_h) - dref[None, :]) * 1e6 / cref[:, None]
            sig_ref = sig_h * 1e6 / cref[:, None] / np.sqrt(np.maximum(n_h, 1))
        hday = np.asarray(P["hour_day"])
        corrs = {k: P["corr"][k] for k in P["corr"]}
        bmt = next((defloat(c["bmt"]) for c in corrs.values() if c.get("bmt")), None)
        thr = {"t_ref": NF.SNR_MIN * sig_ref}
        fbmol = IC._molecular_beta(z_agl, salt, 1064.0)
        if bmt is not None:
            tr = bmt[hday]
            thr["t_ray"] = np.where(np.isfinite(tr), tr, fbmol[None, :])
        else:
            thr["t_ray"] = fbmol[None, :] * np.ones((hour_epoch.size, 1))
        t_all = np.zeros((hour_epoch.size, nz))
        for ident in idents:
            method_i, variant_i = v3["default"][ident]
            rec_i = P["calib"].get(f"{ident}|{method_i}|{variant_i}")
            if not (rec_i and rec_i.get("ok")) or ident not in per_inst:
                continue
            kdi = np.array([np.datetime64(x) for x in rec_i["kal"]["d"]])
            ci = IC.interp_calib(kdi, np.asarray(rec_i["kal"]["v"], "f8"),
                                 hours + np.timedelta64(30, "m"))
            pi = per_inst[ident]
            oki = pi["okrows"]
            sig_i = np.zeros((hour_epoch.size, nz))
            n_i = np.zeros((hour_epoch.size, nz))
            sig_i[oki] = np.nan_to_num(pi["sig"][pi["rows"][oki]])
            n_i[oki] = pi["rn"][pi["rows"][oki]]
            c = corrs[P["corr_of"][f"L1|{ident}"]]
            with np.errstate(all="ignore"):
                scal = sig_i * 1e6 / ci[:, None] * float(c.get("f") or 1.0) / \
                    np.sqrt(np.maximum(n_i, 1))
                wv = defloat(c["wv"]) if c.get("wv") else None
                if wv is not None:
                    wvv = np.where(np.isfinite(wv[hday]), wv[hday], 1.0)
                    scal = scal / wvv
            scal[n_i == 0] = 0.0
            t_all = np.maximum(t_all, scal)
        thr["t_all"] = NF.SNR_MIN * t_all
        sc = _debits(nf["scene"], "<u2")
        for t, key in enumerate(NF.THR_KEYS):
            ref_pass = (n_h > 0) & np.isfinite(np.nan_to_num(beta_cal)) & \
                (np.nan_to_num(beta_cal) >= thr[key])
            mineb = (sc >> (t * 4 + 2) & 1).astype(bool)
            comp = okr[:, None] & np.ones((1, nz), bool)
            # tolerate the exact threshold boundary (interp of C at bin center vs bin mean)
            with np.errstate(all="ignore"):
                nearline = np.abs(np.nan_to_num(beta_cal) - thr[key]) <= \
                    0.02 * np.abs(thr[key])
            hard = comp & ~nearline
            mism = float((mineb[hard] != ref_pass[hard]).mean()) if hard.any() else 0.0
            print(f"  scene {key}: 60-min reconstruction — mismatch {mism * 100:.3f} % "
                  f"(hors bande ±2 % du seuil)")
            worst = max(worst, mism)

    # -- p2 monthly sums, W = 3600 -------------------------------------------------------------
    hm = np.asarray(P["hour_month"])
    for ident in idents:
        if ident not in per_inst:
            continue
        pi = per_inst[ident]
        okr = pi["okrows"]
        sig_h = np.zeros((hour_epoch.size, nz))
        n_h = np.zeros((hour_epoch.size, nz))
        sig_h[okr] = np.nan_to_num(pi["sig"][pi["rows"][okr]])
        n_h[okr] = pi["rn"][pi["rows"][okr]]
        sn = defloat(P["nf"]["p2"][ident]["sn"])
        ss = defloat(P["nf"]["p2"][ident]["ss"])
        sn_ref = np.zeros_like(sn)
        ss_ref = np.zeros_like(ss)
        for mth in range(sn.shape[0]):
            sel = hm == mth
            sn_ref[mth] = n_h[sel].sum(axis=0)
            ss_ref[mth] = (n_h[sel] * sig_h[sel] ** 2).sum(axis=0)
        ok = sn > 0
        with np.errstate(all="ignore"):
            dn = float(np.nanmax(np.abs(sn_ref - sn) / np.maximum(sn, 1))) if ok.any() else 0.0
            rs = np.abs(ss_ref - ss) / np.where(ss > 0, ss, 1)
            ds = float(np.nanmax(rs[ok])) if ok.any() else 0.0
        print(f"  p2 {ident}: sums vs direct — max rel dev n {dn * 100:.2f} %, "
              f"n*sigma^2 {ds * 100:.2f} %")
        worst = max(worst, 0.0 if (dn < 0.01 and ds < 0.05) else 1.0)

    # -- 5-min brute-force spot check + admission approximation (informational) ---------------
    ident = nf.get("noisiest") or idents[0]
    if ident in natives:
        tsec, scr, alt_nat, dark_nat = natives[ident]
        daysec = np.unique(hour_epoch // 24 * 24)
        sample = [daysec[len(daysec) // 3], daysec[2 * len(daysec) // 3]]
        mineb = (_debits(nf["inst"][ident], "u1") & 1).astype(bool)
        n_cmp = n_bad = 0
        approx_dev = []
        for d0h in sample:
            for h in range(24):
                he = d0h + h
                hrow = np.flatnonzero(hour_epoch == he)
                if not hrow.size:
                    continue
                sel = (tsec >= he * 3600) & (tsec < (he + 1) * 3600)
                if not sel.any():
                    continue
                Xh = np.asarray(scr[sel], "f8")
                if dark_nat is not None:
                    Xh = Xh - dark_nat[None, :]
                tw = tsec[sel] // 300
                kn = np.zeros(alt_nat.size)
                nn = np.zeros(alt_nat.size)
                exact_pool = np.full((np.unique(tw).size, alt_nat.size), np.nan)
                for jw, wdd in enumerate(np.unique(tw)):
                    Xw = Xh[tw == wdd]
                    kw = IC.snr_mask(Xw)
                    nw = np.isfinite(Xw).sum(axis=0)
                    kn += kw * nw
                    nn += nw
                    with np.errstate(all="ignore"):
                        mw = np.nanmedian(Xw, axis=0)
                    exact_pool[jw] = np.where(kw, mw, np.nan)
                rk = IC.regrid(kn[None, :], alt_nat, alt_asl)[0]
                rn = IC.regrid(nn[None, :], alt_nat, alt_asl)[0]
                admit = (np.nan_to_num(rn) > 0) & \
                    (np.nan_to_num(rk) / np.where(np.nan_to_num(rn) > 0,
                                                  np.nan_to_num(rn), 1.0) >= NF.NF_FMIN)
                n_cmp += admit.size
                n_bad += int((admit != mineb[hrow[0]]).sum())
                # survivors-median vs all-or-nothing (both regridded), where the hour is admitted
                with np.errstate(all="ignore"):
                    ex = IC.regrid(np.nanmedian(exact_pool, axis=0)[None, :], alt_nat,
                                   alt_asl)[0]
                    un = IC.regrid(np.nanmedian(Xh, axis=0)[None, :], alt_nat, alt_asl)[0]
                good = admit & np.isfinite(ex) & np.isfinite(un) & (np.abs(un) > 0)
                if good.any():
                    approx_dev.append(np.abs(ex[good] - un[good]) /
                                      np.abs(un[good]) * 100)
        if n_cmp:
            mism = n_bad / n_cmp
            print(f"  {ident}: 5-min brute-force on 2 days — mismatch {mism * 100:.3f} % "
                  f"of {n_cmp} pixels")
            worst = max(worst, mism)
        if approx_dev:
            a = np.concatenate(approx_dev)
            print(f"  [info] admission f>={NF.NF_FMIN}: survivors-median vs unconditional "
                  f"hourly median, admitted pixels: median {np.median(a):.2f} %, "
                  f"p95 {np.percentile(a, 95):.2f} % (not gated — the documented "
                  f"approximation of the sub-hourly modes)")
    return worst


def main():
    site = sys.argv[1] if len(sys.argv) > 1 else "payerne"
    only = sys.argv[2] if len(sys.argv) > 2 else None
    P = json.loads((DATA / f"v3_{site}.json").read_text(encoding="utf-8"))
    V2 = json.loads((DATA / f"data_{site}.json").read_text(encoding="utf-8"))
    hday = np.asarray(P["hour_day"])
    b0, b1 = P["band"]
    band = slice(b0, b1 + 1)
    nz = len(P["z"])
    rng = "r0_" + str(len(V2["months"]) - 1)
    print(f"{P['name']}: {len(P['hours'])} paired hours, {nz} display gates, "
          f"band gates {b0}..{b1} ({b1-b0+1})")

    iref = 0
    zmask = np.zeros(nz, bool)
    zmask[band] = True
    worst_pp, worst_r, skipped = 0.0, 0.0, []
    combos = [c for c in V2["combos"] if only is None or c == only]
    print(f"sweeping {len(combos)} v2 combos against the v3 wiring "
          f"(gate: |d| <= {TOL_PP} pt / {TOL_R} on r_log)")
    for ckey in combos:
        m = re.match(r"(.+)_wv([01])_(molecular|none|angstrom)$", ckey)
        if not m:
            skipped.append((ckey, "unparseable combo key"))
            continue
        cal, wv, wl = m.group(1), m.group(2) == "1", m.group(3)
        v2combo = V2["combos"][ckey].get(rng)
        if not v2combo:
            skipped.append((ckey, "no full-range entry"))
            continue
        chans = _combo_channels(site, cal, P)
        if chans[iref][2] is None:
            skipped.append((ckey, "reference channel not expressible"))
            continue
        vals = {}
        for src in P["sources"]:
            for k, (ident, method, variant, db) in enumerate(chans):
                raw = dequant(P["streams"][f"{src}|{ident}"])
                corr = P["corr"][P["corr_of"][f"{src}|{ident}"]]
                # canonical variant stands in for an inexpressible one, ONLY for the mask
                vmask = variant or ("v2.2" if method == "rayleigh" else "cloudWV")
                rec = P["calib"].get(f"{ident}|{method}|{vmask}")
                if not (rec and rec.get("ok")):
                    rec = None
                vals[(src, k)] = transform(raw, rec["c"] if rec else np.ones(len(P["hours"])),
                                           corr, wv, wl, hday,
                                           dark=(db if src == "L1" else None),
                                           calibrated=(src == "L1" and rec is not None))
        have = np.logical_and.reduce([np.any(np.isfinite(vals[(s, k)][:, band]), axis=1)
                                      for s in P["sources"] for k in range(len(chans))])
        print(f"\n== {ckey}: rows {int(have.sum())} (v2 n_hours = {v2combo['n_hours']})")
        print(f"{'src':4s} {'channel':26s} {'medrelbias %':>26s} {'r_log':>16s} "
              f"{'relbias %':>22s} {'N':>14s}")
        for src in P["sources"]:
            for k, (ident, method, variant, db) in enumerate(chans):
                if k == iref:
                    continue
                if variant is None:
                    skipped.append((ckey, f"{src} {ident} {method}: no v3 record for '{cal}'"))
                    continue
                gated = not (src == "L1" and cal in STALE_CALS)
                new = IC._stats(vals[(src, k)][have], vals[(src, iref)][have], zmask)
                old = v2combo[src]["stats"][k]
                row = []
                for key, fmt in (("medrelbias_pct", "%+.4f"), ("r_log", "%.4f"),
                                 ("relbias_pct", "%+.4f"), ("n", "%d")):
                    o, n = old[key], new[key]
                    row.append(f"{fmt % n} / {fmt % o}")
                    if o is None or n is None or not gated:
                        continue
                    if key in ("medrelbias_pct", "relbias_pct"):
                        worst_pp = max(worst_pp, abs(n - o))          # absolute, in points
                    elif key == "r_log":
                        worst_r = max(worst_r, abs(n - o))
                lbl = f"{ident} {method} {variant}" + ("" if gated else "  [info]")
                print(f"{src:4s} {lbl:26s} {row[0]:>26s} {row[1]:>16s} {row[2]:>22s} "
                      f"{row[3]:>14s}")
    if skipped:
        print(f"\nnot expressible from the v3 payload ({len(skipped)}; by design, "
              f"see CAL_MAP_CLOUD):")
        for ckey, why in skipped:
            print(f"  - {ckey}: {why}")
    # ---- the approximation the spec allowed for the profile panel, measured but NOT used --------
    # Well-posed form: run the SAME transform chain but with C_L frozen at its window median,
    # i.e. "divide the window's median profile by the median C_L".  (Dividing an *uncalibrated*
    # transformed profile is not even defined once the molecular wavelength conversion adds its
    # additive term, which lives in calibrated units — that is one reason v3 does it per hour.)
    print("\nprofile panel: exact per-hour C_L division vs a frozen window-median C_L")
    mx = 0.0
    chans = _combo_channels(site, "v2.2", P)
    vals = {}
    for k, (ident, method, variant, db) in enumerate(chans):
        vals[k] = transform(dequant(P["streams"][f"L1|{ident}"]),
                            P["calib"][f"{ident}|{method}|{variant}"]["c"],
                            P["corr"][P["corr_of"][f"L1|{ident}"]], True, "molecular", hday)
    have = np.logical_and.reduce([np.any(np.isfinite(vals[k][:, band]), axis=1)
                                  for k in range(len(chans))])
    for k, (ident, method, variant, db) in enumerate(chans):
        raw = dequant(P["streams"][f"L1|{ident}"])
        corr = P["corr"][P["corr_of"][f"L1|{ident}"]]
        c = np.asarray(P["calib"][f"{ident}|{method}|{variant}"]["c"], "f8")
        cfix = np.full_like(c, float(np.median(c[have])))
        exact = np.nanmedian(transform(raw, c, corr, True, "molecular", hday)[have], axis=0)
        approx = np.nanmedian(transform(raw, cfix, corr, True, "molecular", hday)[have], axis=0)
        m = np.isfinite(exact) & np.isfinite(approx) & (np.abs(exact) > 1e-3)
        d = np.abs(approx[m] - exact[m]) / np.abs(exact[m]) * 100
        if d.size:
            print(f"  {ident} {method:8s} {variant:11s} max {d.max():7.3f} %   "
                  f"median {np.median(d):6.3f} %   (C_L spread "
                  f"{100*(c[have].max()/c[have].min()-1):.1f} %)")
            mx = max(mx, d.max())
    print(f"  -> worst-case profile approximation error {mx:.2f} % — v3 does NOT use it")

    # ---- noise-filter blocks vs independent recomputation ---------------------------------------
    nf_worst = check_nf(P, site)

    # ---- the mechanical gate, LAST so every section above always prints -------------------------
    print(f"\nworst absolute deviation new-vs-v2: {worst_pp:.4f} pt "
          f"(medrelbias/relbias), {worst_r:.5f} (r_log); noise-filter mismatch "
          f"{nf_worst * 100:.3f} %")
    if worst_pp > TOL_PP or worst_r > TOL_R or nf_worst > 0.001:
        print(f"FAIL: deviation exceeds tolerance ({TOL_PP} pt / {TOL_R} / nf 0.1 %)")
        sys.exit(1)
    print("PASS")


if __name__ == "__main__":
    main()
