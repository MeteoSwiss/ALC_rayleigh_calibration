# -*- coding: utf-8 -*-
"""q19 - injection self-test of the vertical-consistency gate (mandatory before the gate is trusted).

Clean tier-1 nights (Nov-Mar, fully clear, window bottom <= 3 km AGL, constant within +-10 % of the
unit's local tier-1 level; selection does NOT use any gate statistic) are perturbed on the NATIVE L1
profiles (before cloud screen, binning, noise estimate, window search, Klett and ensemble -- the whole
operational chain sees the injected signal), then calibrated again:

  (i)   white noise      rcs + a * sigma_night(z) * N(0,1), a = 1, 2, 3, 4, 6 (total noise x sqrt(1+a^2));
                         sigma_night = native first-difference sigma of the night profiles.
                         -> the level a_trig where tier 1 fails and tier 2 rescues ("natural tier 2").
  (ii)  aerosol          beta_a(z) = beta_a0 exp(-z/H), H = 1 or 2 km, beta_a(4 km)/beta_mol(4 km) = s
                         (0.05, 0.10, 0.20, 0.30), with its extinction (LR 50 sr): rcs x (1 + beta_a/beta_m)
                         x exp(-2 int alpha_a).
  (iii) baseline         rcs + B r^2 (a constant offset of the range-normalised signal), B / signal(4 km)
                         = -0.05, +0.02, +0.05, +0.10, +0.20 (the +9..+12 %/km units need ~0.1).
  (iv)  uniform aerosol  beta_a = s_u beta_mol at every height (s_u = 0.10, 0.20) + extinction: the BLIND
                         case (R(z) stays flat, C is biased by ~+s_u).
(ii)-(iv) are run twice: at the night's natural tier-2 noise level a_trig (same noise realisation as
the ladder), and at native noise with tier 1 switched off ("forced tier 2": the gate's discrimination
without the noise limitation). Every case is calibrated by operational v2.2 ("op"), v2.2 with the
measured noise ("fix") and both + the gate. Truth = the clean night's own tier-1 constant.
Outputs only under C:/DATA/Projects/202610_CHM15k_stability/q19/inject/.
"""
from __future__ import annotations

import argparse
import csv
import dataclasses
import json
import os
import sys
import time
import traceback
import warnings
import zlib
from concurrent.futures import ProcessPoolExecutor, as_completed
from datetime import datetime, timedelta
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1]))
sys.path.insert(0, str(HERE))
import q19_gradient_gate as Q  # noqa: E402

INJ = Q.Q19 / "inject"
NOISE_LADDER = (1.0, 2.0, 3.0, 4.0, 6.0)
AEROSOL = [(s, H) for H in (1000.0, 2000.0) for s in (0.05, 0.10, 0.20, 0.30)]
BASELINE = (-0.05, 0.02, 0.05, 0.10, 0.20)
UNIFORM = (0.10, 0.20)
Z_REF = 4000.0
LR_INJ = 50.0

FIELDS = ["wmo", "ident", "date", "case", "kind", "amp", "H", "a_noise", "cond", "variant", "gate",
          "flag", "c", "u", "zb", "zt", "tier", "c_clean", "bias", "s_win", "verdict_first", "verdict",
          "moved", "verdict_main", "lcb", "ucb_lo", "n_up_blk", "best_s", "g_ext", "sg_ext", "g_up",
          "sg_up", "eps_up", "seps_up", "eps_lo", "seps_lo", "z_lo", "z_hi", "zs", "ze", "message"]
MAIN = None


def main_gate():
    """The main rule (first entry of q19_gate_variants.json), used to label every operational selection."""
    global MAIN
    if MAIN is None:
        MAIN = Q.load_gate_variants()[0][1]["gate"]
    return MAIN


def clean_nights(per_unit=6, seed=19):
    """Clean tier-1 nights from the record run (op variant).

    Season/height/level criteria (Nov-Mar, fully clear, tier-1 window bottom <= 3 km AGL, constant within
    +-10 % of the unit's +-30 d tier-1 median) are NOT enough: winter nights can carry elevated layers or
    cirrus (e.g. Sofia 2025-02-12: R 1.8 in the window, 4-4.8 at 4.25 and 5.5-7 km). The injection base
    must be a molecular column, so in addition the tier-1 window must pass the main vertical-consistency
    rule at native noise and no block above it may be > 30 % brighter than the window (no layer). This
    uses the gate on the UNPERTURBED night only: the test then measures what noise / injected aerosol do
    to a verified-clean column (false rejection, detection power), not whether the base is clean."""
    import pickle
    import pandas as pd
    sys.path.insert(0, str(HERE))
    import q19_lib as L
    W = pickle.load(open(Q.Q19 / "q19_blocks_windows.pkl", "rb"))
    Wd = {(w["key"], w["date"]): w for w in W if w["variant"] == "op" and w["which"] == "sel"}
    G = Q.load_gate_variants()[0][1]["gate"]

    def strict(key, date):
        b = Wd.get((key, date))
        if b is None:
            return False
        up = b["is_up"]
        v = L.vc_verdict(b, k=G["k"], eps_tol=G["eps_tol"], tol_lo=G.get("tol_lo", np.inf))
        return v == "accept" and bool(up.any()) and float(np.min(b["exc"][up])) >= -0.30
    rec = Q.load_record()
    op = rec[(rec.variant == "op") & (rec.type == "CHM15k")].copy()
    op = op[op.key != "_".join(Q.PAYERNE)]
    op["t"] = pd.to_datetime(op.date)
    op["ok"] = op.flag.isin(Q.SUCCESS) & (op.c > 0)
    out = []
    rng = np.random.default_rng(seed)
    for key, g in op[op.ok & (op.tier == 1)].groupby("key"):
        g = g.sort_values("t")
        lc = np.log(g.c.values)
        loc = np.array([np.median(lc[np.abs((g.t - t).dt.days.values) <= 30]) for t in g.t])
        m = ((g.t.dt.month.isin([11, 12, 1, 2, 3])) & (g.flag == 1.0) & (g.zb <= 3000)
             & (np.abs(lc - loc) <= 0.10))
        cand = g[m]
        cand = cand[[strict(key, d) for d in cand.date]]
        if len(cand) == 0:
            continue
        pick = cand.iloc[np.sort(rng.choice(len(cand), size=min(per_unit, len(cand)), replace=False))]
        for _, r in pick.iterrows():
            out.append((r.wmo, r.ident, r.date, float(r.c)))
    return out


def night_sigma_rcs(s, d, native):
    """Native per-gate noise (rcs units) of the night profiles: first differences / sqrt(2)."""
    from calibration.io.data_loader import filter_time_range
    from calibration.sensitivity.noise import first_difference_sigma
    o = Q.make_options()
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        nn = filter_time_range(native, d, o)
    rcs = np.asarray(nn.rcs, float)
    t = np.asarray(nn.time, float) * 86400.0
    dt = float(np.median(np.diff(t)))
    sig, _ = first_difference_sigma(rcs, t, dt)
    sig = np.where(np.isfinite(sig), sig, np.nanmedian(sig))
    return sig


def mol_native(s, native):
    """(beta_mol, two-way molecular transmission) on the native range grid (the pipeline's own model)."""
    from calibration.rayleigh.atmosphere import (DEFAULT_STANDARD_ATMOSPHERE, load_standard_atmosphere,
                                                 calculate_molecular_properties)
    info = Q.make_info(s)
    ag = np.asarray(native.altitude_grid, float)
    atm = load_standard_atmosphere(DEFAULT_STANDARD_ATMOSPHERE, ag)
    mp = calculate_molecular_properties(atm.temperature, atm.pressure, np.asarray(native.range_alc, float),
                                        info.instrument_type.wavelength_nm * 1e-9)
    return np.asarray(mp.beta_mol, float), np.asarray(mp.transmission, float)


def aerosol_additive(rcs, z, beta_m, t2_mol, c_clean, s_ref, H, z_ref=Z_REF, lr=LR_INJ, uniform=False):
    """ADDITIVE injection: rcs' = [rcs + C beta_a(z) T_mol^2(z)] exp(-2 tau_a(z)), tau_a = LR int beta_a.

    rcs = C beta T^2 (range-corrected), so the added backscatter is exactly s(z) = beta_a/beta_mol relative
    to molecular in the free troposphere and does not scale the night's own boundary-layer aerosol (the
    first, multiplicative version did, without the matching extinction, and the LR-52 Klett correction then
    inflated C -- replaced). Returns (rcs', s(z))."""
    if uniform:
        ba = s_ref * beta_m
    else:
        i_ref = int(np.argmin(np.abs(z - z_ref)))
        ba0 = s_ref * beta_m[i_ref] / np.exp(-z[i_ref] / H)
        ba = ba0 * np.exp(-np.clip(z, 0, None) / H)
    tau = np.concatenate([[0.0], np.cumsum(0.5 * (ba[1:] + ba[:-1]) * np.diff(z))]) * lr
    out = (np.asarray(rcs, float) + (c_clean * ba * t2_mol)[None, :]) * np.exp(-2.0 * tau)[None, :]
    return out, ba / beta_m


def beta_mol_native(s, native):
    """US-1976 Bucholtz molecular backscatter on the native range grid (the pipeline's own model)."""
    from calibration.rayleigh.atmosphere import (DEFAULT_STANDARD_ATMOSPHERE, load_standard_atmosphere,
                                                 calculate_molecular_properties)
    info = Q.make_info(s)
    ag = np.asarray(native.altitude_grid, float)
    atm = load_standard_atmosphere(DEFAULT_STANDARD_ATMOSPHERE, ag)
    mp = calculate_molecular_properties(atm.temperature, atm.pressure, np.asarray(native.range_alc, float),
                                        info.instrument_type.wavelength_nm * 1e-9)
    return np.asarray(mp.beta_mol, float)


def _with_rcs(native, rcs):
    return dataclasses.replace(native, rcs=rcs)


def add_noise(native, sig_rcs, a, seed):
    if a <= 0:
        return native
    rng = np.random.default_rng(seed)
    rcs = np.asarray(native.rcs, float)
    return _with_rcs(native, rcs + a * sig_rcs[None, :] * rng.standard_normal(rcs.shape))


def aerosol_factor(z, beta_m, s_ref, H, z_ref=Z_REF, lr=LR_INJ):
    """Multiplicative factor (1 + beta_a/beta_m) exp(-2 int_0^z alpha_a) and the excess s(z)."""
    i_ref = int(np.argmin(np.abs(z - z_ref)))
    ba0 = s_ref * beta_m[i_ref] / np.exp(-z[i_ref] / H)
    ba = ba0 * np.exp(-np.clip(z, 0, None) / H)
    tau = np.concatenate([[0.0], np.cumsum(0.5 * (ba[1:] + ba[:-1]) * np.diff(z))]) * lr
    return (1.0 + ba / beta_m) * np.exp(-2.0 * tau), ba / beta_m


def uniform_factor(z, beta_m, s_u, lr=LR_INJ):
    ba = s_u * beta_m
    tau = np.concatenate([[0.0], np.cumsum(0.5 * (ba[1:] + ba[:-1]) * np.diff(z))]) * lr
    return (1.0 + s_u) * np.exp(-2.0 * tau), np.full(z.size, s_u)


def signal_at(native, s, d, z0=Z_REF):
    """Night-median range-normalised signal around z0 (+-150 m), for the baseline amplitude."""
    from calibration.io.data_loader import filter_time_range
    o = Q.make_options()
    nn = filter_time_range(native, d, o)
    z = np.asarray(nn.range_alc, float)
    m = np.abs(z - z0) <= 150
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return float(np.nanmedian(np.asarray(nn.rcs, float)[:, m] / z[m] ** 2))


def _s_win(sz, z, zs, ze):
    if not (np.isfinite(zs) and np.isfinite(ze)):
        return np.nan
    m = (z >= zs) & (z <= ze)
    return float(np.mean(sz[m])) if m.any() else np.nan


RECORDS = {"on": True}


def run_case(s, d, nat, base, c_clean, gates, sz, conds):
    """All variants of one injected native profile set."""
    rows = []
    z = np.asarray(nat.range_alc, float)
    for cond in conds:
        force = cond == "forced"
        runs = [("op", "record", False, None, None), ("fix", "record", True, None, None)]             if RECORDS["on"] else []
        for gname, spec in gates:
            runs.append((spec.get("base", "op"), "gate", spec.get("base", "op") == "fix", gname, spec["gate"]))
        for variant, mode, fix, gname, gpar in runs:
            res, cap = Q.run_one(s, d, nat, mode, gate=gpar, cells=False, sigma_fix=fix, force_t2=force)
            sel = cap.get("sel") or {}
            gsel = cap.get("gsel") or {}
            use = gsel if (mode == "gate" and gsel) else sel
            row = dict(base, cond=cond, variant=variant, gate=gname or "", flag=res["flag"], c=res["c"],
                       u=res["u"], zb=res["zb"], zt=res["zt"], tier=cap.get("tier", np.nan),
                       c_clean=c_clean, message=res["message"])
            if mode == "record" and sel and cap.get("profile") and cap["profile"][3] is not None:
                sg_, p_, z_, s_ = cap["profile"]
                P = Q.prep_profile(sg_, p_, z_, s_)
                bt = Q.candidate_stats(P, sel["s"], sel["e"], main_gate())
                row["verdict_main"] = Q.gate_decision(bt, main_gate())[1]
                for kk in ("lcb", "ucb_lo", "n_up_blk", "best_s"):
                    row[kk] = bt.get(kk, np.nan)
            elif mode == "gate" and gsel:
                for kk in ("lcb", "ucb_lo", "n_up_blk", "best_s"):
                    row[kk] = gsel.get(kk, np.nan)
            row["bias"] = res["c"] / c_clean - 1.0 if np.isfinite(res["c"]) else np.nan
            row["s_win"] = _s_win(sz, z, res["zb"], res["zt"]) if sz is not None else np.nan
            gd = cap.get("gate") or {}
            row.update(verdict=gd.get("verdict", ""), verdict_first=gd.get("first_verdict", ""),
                       moved=gd.get("moved", ""))
            for k in ("g_ext", "sg_ext", "g_up", "sg_up", "eps_up", "seps_up", "eps_lo", "seps_lo",
                      "z_lo", "z_hi", "zs", "ze"):
                row[k] = use.get(k, np.nan)
            rows.append(row)
    return rows


def night_job(wmo, ident, d, c_pub, gates, out_csv, records=True):
    if out_csv.exists():
        return str(out_csv), 0, 0.0
    RECORDS["on"] = bool(records)
    warnings.filterwarnings("ignore")
    import logging
    logging.disable(logging.CRITICAL)
    t0 = time.time()
    s = Q.census_entry(wmo, ident)
    nat0 = Q.load_native(s, d)
    rows = []
    if nat0 is None:
        Q._write_csv(out_csv, rows, FIELDS)
        return str(out_csv), 0, time.time() - t0
    base0 = dict(wmo=wmo, ident=ident, date=d)
    # clean reference (operational v2.2 on the unmodified night = tier 1 by selection)
    res, cap = Q.run_one(s, d, nat0, "record", cells=False)
    c_clean = res["c"]
    if not np.isfinite(c_clean):
        Q._write_csv(out_csv, rows, FIELDS)
        return str(out_csv), 0, time.time() - t0
    sig = night_sigma_rcs(s, d, nat0)
    bm = beta_mol_native(s, nat0)
    z = np.asarray(nat0.range_alc, float)
    seed0 = zlib.crc32(f"{wmo}_{ident}_{d}".encode())
    # forced tier 2 on the CLEAN night at native noise (false-rejection without added noise)
    rows += run_case(s, d, nat0, dict(base0, case="clean", kind="clean", amp=0.0, H=np.nan, a_noise=0.0),
                     c_clean, gates, None, ("forced",))
    # (i) noise ladder -> natural tier 2
    a_trig = None
    for a in NOISE_LADDER:
        nat = add_noise(nat0, sig, a, seed0)
        rr = run_case(s, d, nat, dict(base0, case=f"noise{a:g}", kind="noise", amp=a, H=np.nan, a_noise=a),
                      c_clean, gates, None, ("natural",))
        rows += rr
        if a_trig is None:
            if records:
                if any(r["variant"] == "op" and r["gate"] == "" and r["tier"] == 2 for r in rr):
                    a_trig = a
            else:
                # no record runs in this pass: the natural tier-2 trigger is where operational v2.2 is tier 2
                res_, cap_ = Q.run_one(s, d, nat, "record", cells=False)
                if cap_.get("tier") == 2 and res_["flag"] in Q.SUCCESS:
                    a_trig = a
    # (ii)-(iv) at native noise (forced tier 2) and at a_trig (natural tier 2, same noise realisation)
    cases = [("aerosol", sw, H) for sw, H in AEROSOL] + [("baseline", b, np.nan) for b in BASELINE] \
        + [("uniform", su, np.nan) for su in UNIFORM]
    s4 = signal_at(nat0, s, d)
    for kind, amp, H in cases:
        rcs0 = np.asarray(nat0.rcs, float)
        if kind == "aerosol":
            f, sz = aerosol_factor(z, bm, amp, H)
            rcs1 = rcs0 * f[None, :]
        elif kind == "uniform":
            f, sz = uniform_factor(z, bm, amp)
            rcs1 = rcs0 * f[None, :]
        else:
            sz = None
            rcs1 = rcs0 + amp * s4 * (z ** 2)[None, :]
        nat1 = _with_rcs(nat0, rcs1)
        case = f"{kind}{amp:+g}" + (f"_H{H / 1e3:g}" if np.isfinite(H) else "")
        base = dict(base0, case=case, kind=kind, amp=amp, H=H)
        rows += run_case(s, d, nat1, dict(base, a_noise=0.0), c_clean, gates, sz, ("forced",))
        if a_trig is not None:
            natn = add_noise(nat1, sig, a_trig, seed0)
            rows += run_case(s, d, natn, dict(base, a_noise=a_trig), c_clean, gates, sz, ("natural",))
    Q._write_csv(out_csv, rows, FIELDS)
    return str(out_csv), len(rows), time.time() - t0


def additive_job(wmo, ident, d, a_trig, gates, out_csv):
    """Pass 3: aerosol (8) and uniform (2) cases with ADDITIVE injection, forced (native noise) and natural
    (pass-1 trigger level a_trig, same noise realisation). Records op/fix + the given gate variants."""
    if out_csv.exists():
        return str(out_csv), 0, 0.0
    warnings.filterwarnings("ignore")
    import logging
    logging.disable(logging.CRITICAL)
    RECORDS["on"] = True
    t0 = time.time()
    s = Q.census_entry(wmo, ident)
    nat0 = Q.load_native(s, d)
    rows = []
    if nat0 is None:
        Q._write_csv(out_csv, rows, FIELDS)
        return str(out_csv), 0, time.time() - t0
    res, cap = Q.run_one(s, d, nat0, "record", cells=False)
    c_clean = res["c"]
    if not np.isfinite(c_clean):
        Q._write_csv(out_csv, rows, FIELDS)
        return str(out_csv), 0, time.time() - t0
    sig = night_sigma_rcs(s, d, nat0)
    bm, t2m = mol_native(s, nat0)
    z = np.asarray(nat0.range_alc, float)
    seed0 = zlib.crc32(f"{wmo}_{ident}_{d}".encode())
    base0 = dict(wmo=wmo, ident=ident, date=d)
    cases = [("aerosol", sw, H) for sw, H in AEROSOL] + [("uniform", su, np.nan) for su in UNIFORM]
    for kind, amp, H in cases:
        rcs1, sz = aerosol_additive(nat0.rcs, z, bm, t2m, c_clean, amp, H if np.isfinite(H) else 1.0,
                                    uniform=(kind == "uniform"))
        nat1 = _with_rcs(nat0, rcs1)
        case = f"{kind}{amp:+g}" + (f"_H{H / 1e3:g}" if np.isfinite(H) else "")
        base = dict(base0, case=case, kind=kind, amp=amp, H=H)
        rows += run_case(s, d, nat1, dict(base, a_noise=0.0), c_clean, gates, sz, ("forced",))
        if a_trig is not None and np.isfinite(a_trig):
            natn = add_noise(nat1, sig, a_trig, seed0)
            rows += run_case(s, d, natn, dict(base, a_noise=a_trig), c_clean, gates, sz, ("natural",))
    Q._write_csv(out_csv, rows, FIELDS)
    return str(out_csv), len(rows), time.time() - t0


def main_additive(workers=12):
    """Pass 3 driver: same clean nights; a_trig from pass 1 (inject/: first ladder level with op tier 2)."""
    import glob
    import pandas as pd
    setup = json.load(open(Q.Q19 / "inject" / "inject_setup.json"))
    nights = setup["nights"]
    fr = [pd.read_csv(f, dtype={"date": str}) for f in glob.glob(str(Q.Q19 / "inject" / "*.csv"))]
    P1 = pd.concat([f for f in fr if len(f)], ignore_index=True)
    lad = P1[(P1.kind == "noise") & (P1.variant == "op") & (P1.gate.fillna("") == "") & (P1.tier == 2)
             & P1.flag.isin(Q.SUCCESS)]
    trig = lad.groupby(["wmo", "ident", "date"]).amp.min().to_dict()
    gates = [g for g in Q.load_gate_variants() if g[0] in ("fix_vc2", "fix_vc2_k3", "op_vc2", "fix_blk", "op_gext")]
    out = Q.Q19 / "inject3"
    out.mkdir(parents=True, exist_ok=True)
    tasks = [(w, i, d, trig.get((w, str(i), d), np.nan), gates, out / f"{w}_{i}_{d}.csv") for w, i, d, c in nights]
    print(f"{len(tasks)} nights, gates {[g[0] for g in gates]}", flush=True)
    t0 = time.time()
    done = 0
    with ProcessPoolExecutor(max_workers=workers) as ex:
        futs = [ex.submit(additive_job, *t) for t in tasks]
        for f in as_completed(futs):
            try:
                p_, n, dt = f.result()
            except Exception:  # noqa: BLE001
                traceback.print_exc()
                continue
            done += 1
            print(f"[{done}/{len(tasks)}] {Path(p_).name}: {n} rows in {dt:.0f} s (elapsed {time.time() - t0:.0f} s)",
                  flush=True)


def main():
    if "--additive" in sys.argv:
        main_additive()
        return
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=12)
    ap.add_argument("--per-unit", type=int, default=6)
    ap.add_argument("--limit", type=int, default=0)
    ap.add_argument("--out", default="inject")
    ap.add_argument("--variants", default="")
    ap.add_argument("--no-records", action="store_true")
    args = ap.parse_args()
    global INJ
    INJ = Q.Q19 / args.out
    gates = Q.load_gate_variants()
    if args.variants:
        keep = args.variants.split(",")
        gates = [g for g in gates if g[0] in keep]
    nights = clean_nights(args.per_unit)
    if args.limit:
        nights = nights[:args.limit]
    INJ.mkdir(parents=True, exist_ok=True)
    json.dump(dict(gates=gates, nights=nights), open(INJ / "inject_setup.json", "w"), indent=1)
    print(f"{len(nights)} clean nights, gates {[g[0] for g in gates]}", flush=True)
    tasks = [(w, i, d, c, gates, INJ / f"{w}_{i}_{d}.csv", not args.no_records) for w, i, d, c in nights]
    t0 = time.time()
    done = 0
    with ProcessPoolExecutor(max_workers=args.workers) as ex:
        futs = [ex.submit(night_job, *t) for t in tasks]
        for f in as_completed(futs):
            try:
                p, n, dt = f.result()
            except Exception:  # noqa: BLE001
                traceback.print_exc()
                continue
            done += 1
            print(f"[{done}/{len(tasks)}] {Path(p).name}: {n} rows in {dt:.0f} s "
                  f"(elapsed {time.time() - t0:.0f} s)", flush=True)


if __name__ == "__main__":
    main()
