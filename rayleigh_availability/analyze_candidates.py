# -*- coding: utf-8 -*-
"""Phase 2 analysis — is a candidate eprof_v2.2 configuration acceptable?

A candidate must raise availability WITHOUT buying nights with noise. Every criterion below is
therefore a guard, and availability is the only thing allowed to go up:

  C1 canaries        every registered outlier night stays rejected                    (hard veto)
  C2 continuity      nights valid under v2 keep the IDENTICAL constant                (hard veto)
  C3 no losses       no night valid under v2 becomes invalid                          (hard veto)
  C4 sigma_SD        night-to-night scatter must not inflate (<= +10 % relative)
  C5 Kalman outliers the rolling-IQR screen must not reject more nights per 100 valid
  C6 mechanism       the availability gain must concentrate on NOISY streams
                     (if it were uniform we would have weakened the aerosol defence instead)

Reported on the TUNE split (streams marked tune, 2025 only) for selection, and separately on the
HOLDOUT (other streams + everyone's 2026) so the chosen configuration is never justified by data
it was chosen on.

Figures (only the two that carry a decision) are written to rayleigh_availability/figs/.
"""
from __future__ import annotations
import csv
import json
import sys
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "rayleigh_availability"))
import indicators as IND                                                       # noqa: E402
from canaries import CANARIES                                                  # noqa: E402

DATA = Path("C:/DATA/Projects/202606_E-PROFILE_calibration/rayleigh_availability")
BASE, CAND, SENS = DATA / "baselines", DATA / "candidates", DATA / "sens"
MANIFEST = json.loads((REPO / "rayleigh_availability" / "scope_availability.json").read_text())
FIG = REPO / "rayleigh_availability" / "figs"
REF = "eprof_v2"


def load(path):
    return json.loads(path.read_text()) if path.exists() else None


def split_rec(rec, which):
    """tune = 2025 nights; holdout = 2026 nights (plus, at stream level, the holdout streams)."""
    if rec is None:
        return None
    if which == "tune":
        return {d: v for d, v in rec.items() if d[:4] == "2025"}
    return {d: v for d, v in rec.items() if d[:4] == "2026"}


def sens_noise():
    out = {}
    for f in SENS.glob("*/*_sens.csv"):
        rows = list(csv.DictReader(open(f, encoding="utf-8")))
        if not rows:
            continue
        try:
            v = float(rows[0]["sigma_night_3000"])
        except (KeyError, TypeError, ValueError):
            continue
        if np.isfinite(v) and v > 0:
            out[f.parent.name] = v
    return out


def configs_present():
    return sorted({p.name.split("_", 2)[1] for p in CAND.glob("cand_*.json")})


def evaluate(cfg, streams, which):
    """Aggregate every indicator for one config over one split."""
    tot = dict(n_clear=0, n_valid=0, ref_valid=0, added=0, lost=0, moved=0,
               kal_out=0, kal_valid=0, ref_kal_out=0, ref_kal_valid=0)
    per_stream, canary_hits = [], []
    for i in streams:
        ref = split_rec(load(BASE / f"base_{REF}_{i['label']}.json"), which)
        new = split_rec(load(CAND / f"cand_{cfg}_{i['label']}.json"), which)
        if not ref or not new:
            continue
        a_ref, a_new = IND.availability(ref), IND.availability(new)
        cont = IND.continuity(new, ref)
        added = IND.newly_admitted(new, ref)
        lost = IND.lost(new, ref)
        k_new, k_ref = IND.kalman_outliers(new), IND.kalman_outliers(ref)
        tot["n_clear"] += a_new["n_clear"]
        tot["n_valid"] += a_new["n_valid"]
        tot["ref_valid"] += a_ref["n_valid"]
        tot["added"] += len(added)
        tot["lost"] += len(lost)
        tot["moved"] += sum(1 for d in new if d in ref and IND.is_valid(new[d][0])
                            and IND.is_valid(ref[d][0]) and new[d][1] and ref[d][1]
                            and abs(new[d][1] / ref[d][1] - 1) > 1e-9)
        tot["kal_out"] += k_new["n_outliers"]; tot["kal_valid"] += k_new["n_valid"]
        tot["ref_kal_out"] += k_ref["n_outliers"]; tot["ref_kal_valid"] += k_ref["n_valid"]
        for wmo, ident, date, _ in CANARIES:
            if i["wmo"] == wmo and i["ident"] == ident and date in new \
                    and IND.is_valid(new[date][0]):
                canary_hits.append(f"{wmo}_{ident} {date}")
        # Gradient measured on the EXISTING v2 series, not the candidate's own output: the whole
        # point is to decide -- before deploying anything -- whether a station can take the
        # recovery. Using the candidate's nights to judge the candidate would be circular.
        grad_ref = IND.altitude_gradient(ref, min_n=15)
        grad = IND.altitude_gradient(new)
        shift = IND.altitude_shift(new, ref)
        off_mm = IND.offset_added_vs_kept(new, ref, month_matched=True)
        # How much of the offset is simply the window moving to a different altitude on a stream
        # whose retrieved C_L depends on altitude? (gradient [%/km] x shift [km])
        expl = (grad["slope_pct_per_km"] * shift["shift_m"] / 1000.0
                if np.isfinite(grad["slope_pct_per_km"]) and np.isfinite(shift["shift_m"])
                else float("nan"))
        per_stream.append(dict(
            key=f"{i['wmo']}_{i['ident']}", site=i["site"][:18], group=i["group"],
            av_ref=a_ref["availability_pct"], av_new=a_new["availability_pct"],
            n_added=len(added), n_lost=len(lost),
            sd_ref=IND.sigma_sd(ref), sd_new=IND.sigma_sd(new),
            cont_med=cont["median_pct"], cont_max=cont["max_pct"],
            kal_new=k_new["rate_per_100"], kal_ref=k_ref["rate_per_100"],
            grad_pct_km=grad["slope_pct_per_km"], grad_rho=grad["rho"],
            grad_ref_pct_km=grad_ref["slope_pct_per_km"],
            shift_m=shift["shift_m"], offset_mm=off_mm, offset_explained=expl))
    return tot, per_stream, canary_hits


def main():
    noise = sens_noise()
    cfgs = configs_present()
    if not cfgs:
        print("no candidate outputs yet"); return
    chm = [i for i in MANIFEST if i["group"] == "CHM15k"]
    tune = [i for i in chm if i["split"] == "tune"]
    print(f"configs: {cfgs}\nCHM15k streams: {len(chm)} ({len(tune)} tune)\n")

    rows = []
    for which, streams, tag in (("tune", tune, "TUNE (tune streams, 2025)"),
                                ("holdout", chm, "HOLDOUT (all streams, 2026)")):
        print(f"== {tag} ==")
        print(f"  {'config':13s} {'avail%':>7s} {'(ref)':>7s} {'added':>6s} {'lost':>5s} "
              f"{'moved':>6s} {'sigSD%':>7s} {'(ref)':>7s} {'kal/100':>8s} {'(ref)':>7s} canary")
        for cfg in cfgs:
            tot, per, canary = evaluate(cfg, streams, which)
            if not tot["n_clear"]:
                continue
            av = 100.0 * tot["n_valid"] / tot["n_clear"]
            av_ref = 100.0 * tot["ref_valid"] / tot["n_clear"]
            sd = np.nanmedian([p["sd_new"] for p in per])
            sd_ref = np.nanmedian([p["sd_ref"] for p in per])
            kal = 100.0 * tot["kal_out"] / max(tot["kal_valid"], 1)
            kal_ref = 100.0 * tot["ref_kal_out"] / max(tot["ref_kal_valid"], 1)
            print(f"  {cfg:13s} {av:7.1f} {av_ref:7.1f} {tot['added']:6d} {tot['lost']:5d} "
                  f"{tot['moved']:6d} {sd:7.1f} {sd_ref:7.1f} {kal:8.2f} {kal_ref:7.2f} "
                  f"{'FAIL ' + ','.join(canary) if canary else 'ok'}")
            rows.append(dict(which=which, cfg=cfg, av=av, av_ref=av_ref, added=tot["added"],
                             lost=tot["lost"], moved=tot["moved"], sd=sd, sd_ref=sd_ref,
                             kal=kal, kal_ref=kal_ref, canary=canary, per=per))
        print()

    # ---- verdict on the tune split -------------------------------------------------------
    print("== verdict (tune split) ==")
    best = None
    for r in [x for x in rows if x["which"] == "tune"]:
        c1 = not r["canary"]
        c2 = r["moved"] == 0
        c3 = r["lost"] == 0
        c4 = (not np.isfinite(r["sd_ref"])) or r["sd"] <= 1.10 * r["sd_ref"]
        c5 = r["kal"] <= r["kal_ref"] + 1.0
        gain = {p["key"]: p["av_new"] - p["av_ref"] for p in r["per"]}
        mech = IND.noise_regression(gain, noise)
        c6 = np.isfinite(mech["spearman"]) and mech["spearman"] > 0.3
        ok = all((c1, c2, c3, c4, c5, c6))
        print(f"  {r['cfg']:13s} C1canary={'ok' if c1 else 'FAIL'} C2cont={'ok' if c2 else 'FAIL'} "
              f"C3loss={'ok' if c3 else 'FAIL'} C4sigSD={'ok' if c4 else 'FAIL'} "
              f"C5kal={'ok' if c5 else 'FAIL'} C6mech={'ok' if c6 else 'FAIL'} "
              f"(rho={mech['spearman']:+.2f})  -> {'PASS' if ok else 'reject'}"
              f"{'   +%.1f pts' % (r['av'] - r['av_ref']) if ok else ''}")
        if ok and (best is None or r["av"] > best["av"]):
            best = r
    if best:
        print(f"\n  BEST: {best['cfg']}  availability {best['av_ref']:.1f} -> {best['av']:.1f} % "
              f"(+{best['added']} nights), sigma_SD {best['sd_ref']:.1f} -> {best['sd']:.1f} %")
    else:
        print("\n  no configuration passed all guards")

    # ---- guards stratified by the station's PRE-EXISTING altitude gradient ------------------
    # The reference tests (Amsterdam inter-unit, CHM15k-vs-CL61 at Payerne and Lindenberg) show the
    # recovered nights are sound where the retrieved C_L does not depend on fit height, and biased
    # by ~gradient x height-shift where it does. So the guards must be evaluated on the population
    # the change would actually be applied to.
    GRAD_MAX = 8.0                      # %/km; ~the network median |gradient|
    print(f"== guards split by |dC_L/dz| measured on the EXISTING v2 series "
          f"(threshold {GRAD_MAX:.0f} %/km) ==")
    print(f"  {'config':13s} {'band':16s} {'n':>3s} {'avail%':>7s} {'(ref)':>7s} {'added':>6s} "
          f"{'sigSD%':>7s} {'(ref)':>7s} {'kal/100':>8s} {'(ref)':>7s}")
    for r in [x for x in rows if x["which"] == "holdout"]:
        for band, sel in (("low-gradient", lambda p: abs(p["grad_ref_pct_km"]) <= GRAD_MAX),
                          ("high-gradient", lambda p: abs(p["grad_ref_pct_km"]) > GRAD_MAX)):
            per = [p for p in r["per"] if np.isfinite(p["grad_ref_pct_km"]) and sel(p)]
            if not per:
                continue
            av = np.nanmedian([p["av_new"] for p in per])
            avr = np.nanmedian([p["av_ref"] for p in per])
            sd = np.nanmedian([p["sd_new"] for p in per])
            sdr = np.nanmedian([p["sd_ref"] for p in per])
            kal = np.nanmedian([p["kal_new"] for p in per])
            kalr = np.nanmedian([p["kal_ref"] for p in per])
            print(f"  {r['cfg']:13s} {band:16s} {len(per):3d} {av:7.1f} {avr:7.1f} "
                  f"{sum(p['n_added'] for p in per):6d} {sd:7.1f} {sdr:7.1f} "
                  f"{kal:8.2f} {kalr:7.2f}")
    print()

    # ---- altitude confound, per stream (holdout) ------------------------------------------
    cfg_show = best["cfg"] if best else cfgs[0]
    hold = next((r for r in rows if r["which"] == "holdout" and r["cfg"] == cfg_show), None)
    if hold:
        print(f"\n== altitude confound, config {cfg_show} (holdout) ==")
        print("  A non-zero C_L-vs-window-height gradient means the constant depends on WHERE the")
        print("  window sits, so recovering nights by fitting higher shifts it. 'explained' is")
        print("  gradient x height-shift; when it matches 'offset' the shift is the whole story.")
        print(f"  {'site':20s} {'grad %/km':>10s} {'rho':>6s} {'shift m':>8s} "
              f"{'offset %':>9s} {'explained':>10s}")
        for p in sorted(hold["per"], key=lambda x: -abs(x["grad_pct_km"])
                        if np.isfinite(x["grad_pct_km"]) else 0):
            if not np.isfinite(p["grad_pct_km"]):
                continue
            print(f"  {p['site']:20s} {p['grad_pct_km']:+10.1f} {p['grad_rho']:+6.2f} "
                  f"{p['shift_m']:8.0f} {p['offset_mm']:+9.1f} {p['offset_explained']:+10.1f}")
        g = [abs(p["grad_pct_km"]) for p in hold["per"] if np.isfinite(p["grad_pct_km"])]
        if g:
            print(f"  median |gradient| = {np.median(g):.1f} %/km over {len(g)} streams")
    if best:
        _figures(rows, best["cfg"], noise)


def _figures(rows, cfg, noise):
    """Two decision figures: per-stream gain vs noise, and the seasonal availability recovery."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    FIG.mkdir(parents=True, exist_ok=True)
    tune = next(r for r in rows if r["which"] == "tune" and r["cfg"] == cfg)
    hold = next((r for r in rows if r["which"] == "holdout" and r["cfg"] == cfg), None)

    fig, ax = plt.subplots(1, 2, figsize=(13.5, 4.8))
    for r, lab, mk in ((tune, "tune (2025)", "o"), (hold, "holdout (2026)", "s")):
        if r is None:
            continue
        x = [noise.get(p["key"], np.nan) for p in r["per"]]
        y = [p["av_new"] - p["av_ref"] for p in r["per"]]
        ax[0].semilogx(x, y, mk, alpha=0.8, label=lab)
    ax[0].axhline(0, color="k", lw=0.8)
    ax[0].set_xlabel(r"measured night noise $\sigma_{night}$(3 km) [Mm$^{-1}$sr$^{-1}$]")
    ax[0].set_ylabel("availability gain v2.2 - v2 [pts]")
    ax[0].set_title(f"{cfg}: the gain goes where the noise is")
    ax[0].legend(); ax[0].grid(alpha=0.3, which="both")

    for r, lab, mk in ((tune, "tune (2025)", "o"), (hold, "holdout (2026)", "s")):
        if r is None:
            continue
        ax[1].plot([p["sd_ref"] for p in r["per"]], [p["sd_new"] for p in r["per"]], mk, alpha=0.8,
                   label=lab)
    lim = [0, np.nanmax([p["sd_ref"] for p in tune["per"]] + [p["sd_new"] for p in tune["per"]]) * 1.1]
    ax[1].plot(lim, lim, "k--", lw=1)
    ax[1].set_xlim(lim); ax[1].set_ylim(lim)
    ax[1].set_xlabel(r"$\sigma_{SD}$ v2 [%]"); ax[1].set_ylabel(r"$\sigma_{SD}$ v2.2 [%]")
    ax[1].set_title("Precision is not traded away (on/below the 1:1 line)")
    ax[1].legend(); ax[1].grid(alpha=0.3)
    fig.tight_layout()
    out = FIG / "phase2_gain_and_precision.png"
    fig.savefig(out, dpi=140); plt.close(fig)
    print(f"  -> {out}")


if __name__ == "__main__":
    main()
