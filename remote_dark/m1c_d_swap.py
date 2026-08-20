# -*- coding: utf-8 -*-
"""M1-c — full reproduction of the Schiphol-D optical-module swap from sky data alone.

L1 metadata says unit D changed optical module TUB150037 -> TUB160055 on 2026-07-11 (A/B/C
unchanged). The remote-dark question: is that hardware change VISIBLE, and quantifiable, from
clear-sky nights only — the exact capability the campaign is after (detect instrument changes
without a site visit)?

Three views, all on the merged cache (local 2025-01..2026-07-12 + balfrin extension
2026-06-01..2026-08-13; for duplicated dates the LOCAL night is kept):

  P1  per-night noise floor vs date, all four units — the detection view. The floor is the
      top-gate P-view noise (electronics-dominated); a module swap moves it as a STEP on the
      swap date, and the co-located controls say what the sky did.
  P2  D's noise profile, pre vs post, season-matched (pre = 2026-06-01..07-09 only) — where in
      altitude the new module differs.
  P3  era change of the v4 evidence intercept d(z) (per-gate OLS of S on M over nights,
      season-matched eras) for ALL four units. The sky's seasonal drift is common to the four;
      D's excess over the controls is the module's additive-baseline change.

Era buffer: the swap day and its neighbours (2026-07-10..07-11) belong to neither era.
"""
from __future__ import annotations

import numpy as np

from remote_dark.common import V4_DIR, ensure_out

WMO = "0-20000-0-06240"
IDENTS = "ABCD"
SWAP = 20260711
PRE_ERA = (20260601, 20260709)      #: season-matched pre era (same summer, adjacent weeks)
POST_ERA = (20260712, 20260813)


def load_merged(ident):
    """Local cache + balfrin extension, duplicated dates resolved in favour of the local night."""
    loc = np.load(V4_DIR / "cache" / f"{WMO}_{ident}_20250101_20260813_220.npz",
                  allow_pickle=True)
    ext = np.load(V4_DIR / "cache_ext" / f"{WMO}_{ident}_ext_20260601_20260813.npz",
                  allow_pickle=True)
    if loc["rng"].size != ext["rng"].size or not np.allclose(loc["rng"], ext["rng"]):
        raise RuntimeError(f"{ident}: range grids differ between cache and extension")
    new = ~np.isin(ext["dates"], loc["dates"])
    out = {"rng": loc["rng"]}
    for k in ("S", "sigma", "M", "n_prof", "dates"):
        out[k] = np.concatenate([loc[k], ext[k][new]], axis=0)
    order = np.argsort(out["dates"])
    for k in ("S", "sigma", "M", "n_prof", "dates"):
        out[k] = out[k][order]
    return out


def pview_noise(c):
    """Per-night per-gate noise per profile in P-view (rcs/z^2) — the electronics-plain view."""
    z2 = np.where(c["rng"] > 0, c["rng"], np.nan) ** 2
    n_prof = np.maximum(c["n_prof"].astype(float), 1.0)
    return c["sigma"] * np.sqrt(n_prof)[:, None] / z2[None, :]


def night_floor(c):
    noise = pview_noise(c)
    top = c["rng"] >= np.nanpercentile(c["rng"], 92)
    return np.nanmedian(noise[:, top], axis=1)


def evidence_intercept(c, era):
    """Per-gate OLS of S on M over the era's nights: S_n(z) = c(z) M_n(z) + d(z).
    Returns d(z) — the v4 evidence intercept (dark + common contamination)."""
    m = (c["dates"] >= era[0]) & (c["dates"] <= era[1])
    S, M = c["S"][m], c["M"][m]
    n = np.isfinite(S) & np.isfinite(M)
    cnt = n.sum(axis=0).astype(float)
    Sm = np.where(n, S, 0.0)
    Mm = np.where(n, M, 0.0)
    sM = Mm.sum(axis=0)
    sS = Sm.sum(axis=0)
    sMM = (Mm * Mm).sum(axis=0)
    sMS = (Mm * Sm).sum(axis=0)
    with np.errstate(invalid="ignore", divide="ignore"):
        den = cnt * sMM - sM * sM
        d = (sMM * sS - sM * sMS) / den
    d[cnt < 5] = np.nan
    return d, int(m.sum())


def dates_to_num(dates):
    import matplotlib.dates as mdates
    from datetime import datetime
    return mdates.date2num([datetime.strptime(str(int(x)), "%Y%m%d") for x in dates])


def main():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import matplotlib.dates as mdates
    from datetime import datetime
    FIG = ensure_out("figs")

    units = {i: load_merged(i) for i in IDENTS}
    fig, (ax1, ax2, ax3) = plt.subplots(
        1, 3, figsize=(15.5, 5.4), gridspec_kw={"width_ratios": [1.9, 1.0, 1.0]})

    # ---- P1: floor vs date, the detection view -----------------------------------------------
    swap_num = mdates.date2num(datetime(2026, 7, 11))
    for ident in IDENTS:
        c = units[ident]
        fl = night_floor(c)
        m = c["dates"] >= 20260401
        style = (dict(color="crimson", lw=1.6, marker="o", ms=3.5, zorder=5) if ident == "D"
                 else dict(color="0.65", lw=1.0, marker=".", ms=2.5))
        ax1.plot(dates_to_num(c["dates"][m]), fl[m], label=f"unit {ident}", **style)
    ax1.axvline(swap_num, color="k", lw=1.2, ls="--")
    ax1.text(swap_num, ax1.get_ylim()[1], "  module swap 2026-07-11\n  TUB150037 → TUB160055",
             fontsize=8, va="top")
    ax1.set_yscale("log")
    ax1.xaxis.set_major_formatter(mdates.DateFormatter("%b %d"))
    ax1.set_xlabel("night (2026)")
    ax1.set_ylabel("noise floor per night, P-view (top gates)")
    ax1.grid(alpha=0.25, which="both")
    ax1.legend(fontsize=8, loc="lower left", ncol=2)

    # ---- P2: D pre/post profile, season-matched ----------------------------------------------
    cD = units["D"]
    noise = pview_noise(cD)
    rng = cD["rng"]
    pre = (cD["dates"] >= PRE_ERA[0]) & (cD["dates"] <= PRE_ERA[1])
    post = (cD["dates"] >= POST_ERA[0]) & (cD["dates"] <= POST_ERA[1])
    m = (rng >= 100) & (rng <= 15000)
    for sel, lab, col in ((pre, f"TUB150037 (pre, {int(pre.sum())} n)", "k"),
                          (post, f"TUB160055 (post, {int(post.sum())} n)", "crimson")):
        med = np.nanmedian(noise[sel], axis=0)
        lo, hi = np.nanpercentile(noise[sel], [25, 75], axis=0)
        ax2.fill_betweenx(rng[m] / 1e3, lo[m], hi[m], alpha=0.25,
                          color=("0.5" if col == "k" else col))
        ax2.plot(med[m], rng[m] / 1e3, "-", color=col, lw=1.7, label=lab)
    ax2.set_xscale("log")
    ax2.set_xlabel("noise per profile, P-view")
    ax2.set_ylabel("altitude (km)")
    ax2.set_ylim(0, 15)
    ax2.grid(alpha=0.25, which="both")
    ax2.legend(fontsize=8, loc="upper right")
    top = rng >= np.nanpercentile(rng, 92)
    f_pre = float(np.nanmedian(np.nanmedian(noise[pre], axis=0)[top]))
    f_post = float(np.nanmedian(np.nanmedian(noise[post], axis=0)[top]))
    ax2.set_title(f"unit D noise profile\nfloor ×{f_post / f_pre:.2f} after swap", fontsize=10)

    # ---- P3: era change of the evidence intercept, all units ---------------------------------
    for ident in IDENTS:
        c = units[ident]
        d_pre, n_pre = evidence_intercept(c, PRE_ERA)
        d_post, n_post = evidence_intercept(c, POST_ERA)
        z2 = np.where(c["rng"] > 0, c["rng"], np.nan) ** 2
        dd = (d_post - d_pre) / z2
        # light median smoothing for display only
        k = 5
        pad = np.pad(dd, k // 2, mode="edge")
        dd_s = np.array([np.nanmedian(pad[i:i + k]) for i in range(dd.size)])
        mm = (c["rng"] >= 1000) & (c["rng"] <= 15000)
        style = (dict(color="crimson", lw=1.8, zorder=5) if ident == "D"
                 else dict(color="0.65", lw=1.1))
        ax3.plot(dd_s[mm], c["rng"][mm] / 1e3,
                 label=f"{ident} ({n_pre}/{n_post} n)", **style)
    ax3.axvline(0, color="k", lw=0.6)
    ax3.set_xlabel("Δ evidence intercept post − pre, P-view")
    ax3.set_ylabel("altitude (km)")
    ax3.set_ylim(0, 15)
    ax3.grid(alpha=0.25)
    ax3.legend(fontsize=8, loc="upper right")
    ax3.set_title("additive-baseline change\n(controls = the sky's seasonal part)", fontsize=10)

    fig.suptitle("Schiphol D optical-module swap (2026-07-11), reproduced from clear-sky nights "
                 "alone — detection (left), noise change (centre), dark change vs controls "
                 "(right)", fontsize=11)
    fig.tight_layout(rect=[0, 0, 1, 0.93])
    f = FIG / "fig12_d_swap_full.png"
    fig.savefig(f, dpi=140)
    plt.close(fig)
    print(f)
    print(f"D floor pre {f_pre:.3e} -> post {f_post:.3e}  = x{f_post / f_pre:.2f} "
          f"({int(pre.sum())} vs {int(post.sum())} nights, season-matched)")
    for ident in "ABC":
        c = units[ident]
        fl = night_floor(c)
        p1 = (c["dates"] >= PRE_ERA[0]) & (c["dates"] <= PRE_ERA[1])
        p2 = (c["dates"] >= POST_ERA[0]) & (c["dates"] <= POST_ERA[1])
        print(f"control {ident}: floor x{np.nanmedian(fl[p2]) / np.nanmedian(fl[p1]):.2f} "
              f"({int(p1.sum())}/{int(p2.sum())} nights)")


if __name__ == "__main__":
    main()
