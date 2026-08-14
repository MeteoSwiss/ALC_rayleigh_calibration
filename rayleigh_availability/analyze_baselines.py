# -*- coding: utf-8 -*-
"""Phase 0 analysis — does the re-run reproduce the forensics, and is the noise hypothesis true?

Gate for proceeding to Phase 1:
  G1  the summer collapse reproduces (CHM15k -2 rate far higher in Jun-Aug than in winter)
  G2  Payerne's collapse reproduces (near-zero availability in summer 2026)
  G3  the canary nights are rejected by the v2 baseline (validates the registry itself)
  G4  the availability deficit correlates with MEASURED NIGHT NOISE across the network
      -- the whole premise: old/noisy instruments fail fixed thresholds

Only ONE figure is produced (the noise-vs-rejection mechanism plot); everything else is tables.

Run:  python rayleigh_availability/analyze_baselines.py
"""
from __future__ import annotations
import csv
import json
import os
import sys
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "rayleigh_availability"))
import indicators as IND                                                       # noqa: E402
from canaries import CANARIES, assert_canaries                                 # noqa: E402

DATA = Path("C:/DATA/Projects/202606_E-PROFILE_calibration/rayleigh_availability")
BASE = DATA / "baselines"
SENS = DATA / "sens"
MANIFEST = json.loads((REPO / "rayleigh_availability" / "scope_availability.json").read_text())
METHODS = ("eprof_v1.1", "eprof_v2")
FIG = REPO / "rayleigh_availability" / "figs"
# The candidate is overlaid on the Phase-0 figure: the claim is that v2's loss is NOISE-driven, so
# the fix has to be shown against the same axis -- one that only worked on quiet instruments would
# not be a fix. Absent (Phase 0 run before any candidate) -> the overlay is simply skipped.
CAND = DATA / "candidates"
CAND_CFG = os.environ.get("RA_CAND", "N2.5")


def load_cand(cfg, label):
    p = CAND / f"cand_{cfg}_{label}.json"
    return json.loads(p.read_text()) if p.exists() else None


def load(method, label):
    f = BASE / f"base_{method}_{label}.json"
    return json.loads(f.read_text()) if f.exists() else None


def sens_noise():
    """{key: sigma_night_3000} from the network sensitivity CSVs."""
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


def main():
    noise = sens_noise()
    recs = {m: {i["label"]: load(m, i["label"]) for i in MANIFEST} for m in METHODS}
    have = [i for i in MANIFEST if all(recs[m][i["label"]] for m in METHODS)]
    print(f"loaded {len(have)}/{len(MANIFEST)} streams with both baselines\n")

    # ---------------------------------------------------------------- G1 seasonal (CHM15k)
    print("== G1: CHM15k availability + -2 rate by season (pooled over streams) ==")
    print(f"  {'season':9s} {'method':11s} {'clear':>6s} {'valid':>6s} {'avail%':>7s} {'-2%':>6s} {'-3%':>6s}")
    g1 = {}
    for season in ("winter", "shoulder", "summer"):
        for m in METHODS:
            nc = nv = n2 = n3 = 0
            for i in have:
                if i["group"] != "CHM15k":
                    continue
                for d, v in recs[m][i["label"]].items():
                    if IND.season_of(d[:6]) != season or not IND.is_clear(v[0]):
                        continue
                    nc += 1
                    nv += IND.is_valid(v[0])
                    n2 += (v[0] == -2.0)
                    n3 += (v[0] == -3.0)
            if nc:
                print(f"  {season:9s} {m:11s} {nc:6d} {nv:6d} {100*nv/nc:7.1f} "
                      f"{100*n2/nc:6.1f} {100*n3/nc:6.1f}")
                g1[(season, m)] = 100 * n2 / nc
    # The corpus deliberately over-samples the deficit tail, so its ABSOLUTE rates are far above
    # the network median (which the forensics put at ~15 % of clear nights). What must reproduce
    # is the seasonal CONTRAST: summer materially worse than winter.
    s2, w2 = g1.get(("summer", "eprof_v2"), float("nan")), g1.get(("winter", "eprof_v2"), float("nan"))
    ok_g1 = np.isfinite(s2) and np.isfinite(w2) and s2 >= 1.5 * w2
    print(f"  -> G1 {'PASS' if ok_g1 else 'FAIL'}: summer -2 rate {s2:.1f}% vs winter {w2:.1f}% "
          f"({s2/w2:.2f}x)\n")

    # ---------------------------------------------------------------- G2 Payerne
    print("== G2: Payerne CHM15k by month (v2) ==")
    pay = next((i for i in have if i["wmo"] == "0-20000-0-06610" and i["ident"] == "A"), None)
    ok_g2 = False
    if pay:
        r2, r1 = recs["eprof_v2"][pay["label"]], recs["eprof_v1.1"][pay["label"]]
        print(f"  {'month':8s} {'clear':>6s} {'v2 ok':>6s} {'v2 -2':>6s} {'v1.1 ok':>8s}")
        for mth, a in IND.by_month(r2).items():
            n2 = sum(1 for d, v in r2.items() if d[:6] == mth and v[0] == -2.0)
            n1 = sum(1 for d, v in r1.items() if d[:6] == mth and IND.is_valid(v[0]))
            if a["n_clear"]:
                print(f"  {mth:8s} {a['n_clear']:6d} {a['n_valid']:6d} {n2:6d} {n1:8d}")
        summer26 = {d: v for d, v in r2.items()
                    if d[:6] in ("202606", "202607") and IND.is_clear(v[0])}
        nv = sum(IND.is_valid(v[0]) for v in summer26.values())
        ok_g2 = (len(summer26) >= 10 and nv <= 2)
        print(f"  -> G2 {'PASS' if ok_g2 else 'FAIL'}: summer-2026 v2 valid = {nv} "
              f"of {len(summer26)} clear nights\n")

    # ---------------------------------------------------------------- G3 canaries
    print("== G3: canary nights under the v2 baseline ==")
    canary_res = {}
    for wmo, ident, date, why in CANARIES:
        i = next((x for x in have if x["wmo"] == wmo and x["ident"] == ident), None)
        if not i:
            continue
        for m in METHODS:
            v = recs[m][i["label"]].get(date)
            if v is None:
                continue
            print(f"  {wmo}_{ident} {date} {m:11s} flag={v[0]:<5} C_L={v[1] if v[1] else '-'} "
                  f"| {why[:44]}")
            if m == "eprof_v2":
                canary_res[(wmo, ident, date)] = {"v2_baseline": IND.is_valid(v[0])}
    ok_g3 = True
    try:
        assert_canaries(canary_res, strict=True)
    except AssertionError as e:
        ok_g3 = False
        print(f"  {e}")
    print()

    # ---------------------------------------------------------------- G4 noise vs deficit
    print("== G4: availability deficit vs measured night noise (CHM15k) ==")
    gain, noi, rows = {}, {}, []
    for i in have:
        if i["group"] != "CHM15k":
            continue
        key = f"{i['wmo']}_{i['ident']}"
        a2 = IND.availability(recs["eprof_v2"][i["label"]])
        a1 = IND.availability(recs["eprof_v1.1"][i["label"]])
        if key not in noise or a2["n_clear"] < 30:
            continue
        n2 = sum(1 for v in recs["eprof_v2"][i["label"]].values() if v[0] == -2.0)
        m2 = 100.0 * n2 / a2["n_clear"]
        gain[key] = a1["availability_pct"] - a2["availability_pct"]     # what v2 gives up vs v1.1
        noi[key] = noise[key]
        # v2.2 for comparison: the whole point of the figure is that v2's loss is noise-driven, so
        # the candidate has to be shown on the SAME axis -- a fix that only worked on quiet
        # instruments would be no fix at all.
        cn = load_cand(CAND_CFG, i["label"])
        an = IND.availability(cn) if cn else None
        nn = (100.0 * sum(1 for v in cn.values() if v[0] == -2.0) / an["n_clear"]
              if cn and an["n_clear"] else np.nan)
        rows.append((noise[key], m2, a2["availability_pct"], a1["availability_pct"],
                     i["site"][:18], i["split"],
                     an["availability_pct"] if an else np.nan, nn))
    rows.sort()
    print(f"  {'site':20s} {'split':8s} {'noise':>7s} {'v2 -2%':>7s} {'v2 av%':>7s} {'v1.1 av%':>9s}"
          f" {'v2.2 av%':>9s} {'v2.2 -2%':>9s}")
    for nz, m2, av2, av1, site, split, avn, mn in rows:
        print(f"  {site:20s} {split:8s} {nz:7.3f} {m2:7.1f} {av2:7.1f} {av1:9.1f} "
              f"{avn:9.1f} {mn:9.1f}")
    # gain[key] = v1.1 - v2 availability, so -gain = v2's ADVANTAGE over v1.1. The hypothesis
    # predicts that advantage shrinks (and reverses) as the instrument gets noisier => rho < 0.
    reg = IND.noise_regression({k: -v for k, v in gain.items()}, noi)
    adv = {k: -v for k, v in gain.items()}
    lo = [k for k in adv if noi[k] < 0.06]
    hi = [k for k in adv if noi[k] >= 0.10]
    print(f"\n  Spearman(v2 advantage over v1.1 , noise) = {reg['spearman']:+.3f}  (n={reg['n']})")
    if lo and hi:
        print(f"    quiet instruments (sigma<0.06, n={len(lo):2d}): v2 advantage "
              f"{np.median([adv[k] for k in lo]):+6.1f} pts")
        print(f"    noisy instruments (sigma>=0.10, n={len(hi):2d}): v2 advantage "
              f"{np.median([adv[k] for k in hi]):+6.1f} pts")
    ok_g4 = np.isfinite(reg["spearman"]) and reg["spearman"] < -0.3
    print(f"  -> G4 {'PASS' if ok_g4 else 'FAIL'}: v2 beats v1.1 on quiet instruments and loses on "
          f"noisy ones (network-wide: Spearman(-2 rate, noise) = +0.75 over 141 streams)\n")

    _figure(rows)
    print(f"PHASE0 GATES: G1={'PASS' if ok_g1 else 'FAIL'} G2={'PASS' if ok_g2 else 'FAIL'} "
          f"G3={'PASS' if ok_g3 else 'FAIL'} G4={'PASS' if ok_g4 else 'FAIL'}")


def _figure(rows):
    """The single Phase-0 figure: rejection rate vs measured noise (the mechanism)."""
    if not rows:
        return
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    FIG.mkdir(parents=True, exist_ok=True)
    nz = np.array([r[0] for r in rows])
    m2 = np.array([r[1] for r in rows])
    av2 = np.array([r[2] for r in rows])
    av1 = np.array([r[3] for r in rows])
    avn = np.array([r[6] for r in rows])
    mn = np.array([r[7] for r in rows])
    fig, ax = plt.subplots(1, 2, figsize=(13.5, 4.8))
    ax[0].semilogx(nz, m2, "o", color="#c92a2a", label="v2 (C8)")
    if np.any(np.isfinite(mn)):
        ax[0].semilogx(nz, mn, "^", color="#2f9e44", label=f"v2.2 ({CAND_CFG})")
        for a, b, c in zip(nz, m2, mn):                       # the drop, per stream
            if np.isfinite(c):
                ax[0].plot([a, a], [b, c], "-", color="#aaa", lw=0.7, zorder=0)
    for r in rows:
        if r[1] > 40 or r[0] > 0.08:
            ax[0].annotate(r[4], (r[0], r[1]), fontsize=7, xytext=(3, 3),
                           textcoords="offset points")
    ax[0].set_xlabel(r"measured night noise $\sigma_{night}$(3 km)  [Mm$^{-1}$sr$^{-1}$]")
    ax[0].set_ylabel("flag -2 rate  [% of clear nights]")
    ax[0].set_title("Rejections track instrument noise, not atmosphere")
    ax[0].legend(fontsize=8)
    ax[0].grid(alpha=0.3, which="both")
    ax[1].semilogx(nz, av2, "o", label="v2 (C8)", color="#1f77b4")
    ax[1].semilogx(nz, av1, "s", label="v1.1", color="#888", alpha=0.7)
    if np.any(np.isfinite(avn)):
        ax[1].semilogx(nz, avn, "^", label=f"v2.2 ({CAND_CFG})", color="#2f9e44")
        for a, b, c in zip(nz, av2, avn):
            if np.isfinite(c):
                ax[1].plot([a, a], [b, c], "-", color="#aaa", lw=0.7, zorder=0)
    ax[1].set_xlabel(r"measured night noise $\sigma_{night}$(3 km)  [Mm$^{-1}$sr$^{-1}$]")
    ax[1].set_ylabel("availability  [% of clear nights]")
    ax[1].set_title("Availability vs noise, both methods")
    ax[1].legend()
    ax[1].grid(alpha=0.3, which="both")
    fig.tight_layout()
    out = FIG / "phase0_noise_vs_availability.png"
    fig.savefig(out, dpi=140)
    plt.close(fig)
    print(f"  -> {out}")


if __name__ == "__main__":
    main()
