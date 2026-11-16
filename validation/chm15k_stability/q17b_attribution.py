# -*- coding: utf-8 -*-
"""q17b (attribution): fraction of the CHM15k annual cycle that is aerosol vs instrument gain.

f_aer = s / r, where s = CL61 in-phase annual response per unit CHM15k annual response (q17b_bands.csv,
s_harm with its month-bootstrap CI; s_month as a check) and r = CL61/CHM15k aerosol-excess ratio:
uniform over [0.54 (AE 0), 0.78 (empirical: synoptic slope at 2-3 km, where aerosol dominates)].
Pooled over the free-tropospheric bands 4-5, 5-6, 6-8 km (per site: band mean, SE = mean band SE, not
reduced - bands are correlated) then inverse-variance across Aosta, Camborne, Lindenberg. Payerne is
reported separately (abnormal CHM15k, CL61 only Feb-Oct 2026).
Two CL61 versions: 'cl_wv' (operational WV correction) and 'cl_wv+cloudgain' (CL61 additionally divided
by its own liquid-cloud constant, removing CL61 gain drifts such as the Lindenberg spring artefact).
Two CHM15k targets: raw ln X (whole cycle) and kappa-normalised (the residual after the T_int term).

Also the height decomposition G(z) = A_CHM(z) - A_CL61,in(z)/r (gain part, should be flat) and the
window-level cross-check with the operational CL61 Rayleigh/cloud constants.
Outputs: q17/q17b_attribution.csv, q17b_decomposition.csv
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1]))
sys.path.insert(0, str(HERE))
from q17b_extract import Q17  # noqa: E402

R_RANGE = (0.536, 0.78)
FT_BANDS = ("4-5", "5-6", "6-8")
rng = np.random.default_rng(171)


def pooled(b, ycol, prefix, sites, bands):
    per_site = []
    for s in sites:
        d = b[(b.site == s) & (b.ycol == ycol) & (b.band.isin(bands))]
        if d.empty:
            continue
        sv = d[f"{prefix}s_harm"].to_numpy(float)
        se = (d[f"{prefix}s_harm_hi"] - d[f"{prefix}s_harm_lo"]).to_numpy(float) / 3.92
        ok = np.isfinite(sv) & np.isfinite(se)
        if not ok.any():
            continue
        per_site.append((s, float(np.mean(sv[ok])), float(np.mean(se[ok])),
                         float(np.nanmean(d[f"{prefix}s_month"])), float(np.nanmean(d[f"{prefix if prefix else ''}amp_x"]))))
    if not per_site:
        return None, []
    w = np.array([1 / p[2] ** 2 for p in per_site])
    s_bar = float(np.sum(w * np.array([p[1] for p in per_site])) / w.sum())
    se_bar = float(1 / np.sqrt(w.sum()))
    # inflate by the between-site scatter when it exceeds the formal SE (random-effects style)
    sv = np.array([p[1] for p in per_site])
    if len(sv) > 1:
        se_bar = float(max(se_bar, np.std(sv, ddof=1) / np.sqrt(len(sv))))
    sm = rng.normal(s_bar, se_bar, 200000)
    rr = rng.uniform(*R_RANGE, 200000)
    f = sm / rr
    return dict(s=s_bar, s_se=se_bar, f_med=float(np.median(f)), f_lo=float(np.percentile(f, 2.5)),
                f_hi=float(np.percentile(f, 97.5)), f_central=s_bar / 0.65), per_site


def main():
    b = pd.read_csv(Q17 / "q17b_bands.csv")
    rows = []
    for ycol in ("cl_wv", "cl_wv+cloudgain", "cl_wvdark", "cl_wv+aeronet", "cl_nowv"):
        for prefix, target in (("", "raw"), ("k_", "kappa")):
            for label, sites in (("pooled_3sites", ("Aosta", "Camborne", "Lindenberg")),
                                 ("pooled_AoCa", ("Aosta", "Camborne")),
                                 ("Aosta", ("Aosta",)), ("Camborne", ("Camborne",)),
                                 ("Lindenberg", ("Lindenberg",)), ("Payerne", ("Payerne",))):
                for bl, bands in (("FT_4-8km", FT_BANDS), ("3-4km", ("3-4",)), ("6-8km", ("6-8",))):
                    res, per = pooled(b, ycol, prefix, sites, bands)
                    if res is None:
                        continue
                    rows.append(dict(cl61=ycol, chm=target, sites=label, bands=bl, n_sites=len(per), **res))
    out = pd.DataFrame(rows)
    out.to_csv(Q17 / "q17b_attribution.csv", index=False)
    pd.set_option("display.width", 250)
    print(out[(out.bands == "FT_4-8km")].round(3).to_string())

    # height decomposition (cl_wv; Lindenberg also cloudgain)
    dec = []
    for (s, ycol), d in b.groupby(["site", "ycol"]):
        if ycol not in ("cl_wv", "cl_wv+cloudgain"):
            continue
        for _, r in d.iterrows():
            for tgt, ax in (("raw", r["amp_x"]), ("kappa", r.get("k_amp_x", np.nan))):
                yin = r["amp_y_inphase"] if tgt == "raw" else r.get("k_s_harm", np.nan) * r.get("k_amp_x", np.nan)
                dec.append(dict(site=s, cl61=ycol, chm=tgt, band=r["band"], A_chm=ax, A_cl61_inphase=yin,
                                aerosol_1064=yin / 0.65, gain=ax - yin / 0.65,
                                gain_lo=ax - yin / R_RANGE[0], gain_hi=ax - yin / R_RANGE[1]))
    pd.DataFrame(dec).to_csv(Q17 / "q17b_decomposition.csv", index=False)


if __name__ == "__main__" and len(sys.argv) == 1:
    main()


# ---------------------------------------------------------------------------------------------
# Additive in-phase budget of the CHM15k cycle: total = kappa*T_int + aerosol (CL61/r) + unexplained gain
# ---------------------------------------------------------------------------------------------
def budget():
    from q17b_analyse import (SITES, cl61_valid, cloud_gain, flag_iqr, harmonic, load_kappa_in, month_key)
    df = pd.read_csv(Q17 / "q17b_nights.csv")
    k_in = load_kappa_in()
    rows = []
    for site, (wmo, ich, _, _, _) in SITES.items():
        kap = k_in.get(f"{wmo}_{ich}", np.nan)
        d = df[df.site == site]
        for variant in ("cl_wv", "cl_wv+aeronet", "cl_wv+cloudgain", "cl_wv+aeronet+cloudgain"):
            acc = []
            for bn in FT_BANDS + ("3-4",):
                t = d["t"].to_numpy(float)
                x = d[f"chm_{bn}"].to_numpy(float)
                y = d[f"cl_wv_{bn}"].to_numpy(float)
                T = d["T_int"].to_numpy(float)
                if "aeronet" in variant:
                    x = x + 2 * d["tau1064"].to_numpy(float)
                    y = y + 2 * d["tau910"].to_numpy(float)
                if "cloudgain" in variant:
                    y = y - cloud_gain(site, t)
                ok = np.isfinite(x) & np.isfinite(y) & np.isfinite(T) & cl61_valid(site, t)
                t, x, y, T = t[ok], x[ok], y[ok], T[ok]
                bad = flag_iqr(x) | flag_iqr(y)
                t, x, y, T = t[~bad], x[~bad], y[~bad], T[~bad]
                if len(t) < 25:
                    continue
                kt = kap / 100 * (T - np.median(T))
                trend = np.ptp(t) > 540

                def comp(idx):
                    px, py, pk = harmonic(t[idx], [x[idx], y[idx], kt[idx]], trend and np.ptp(t[idx]) > 540)
                    u = px / abs(px)
                    a_tot = 100 * abs(px)
                    a_k = 100 * float(np.real(pk * np.conj(u)))
                    a_y = 100 * float(np.real(py * np.conj(u)))
                    return a_tot, a_k, a_y

                a_tot, a_k, a_y = comp(np.arange(len(t)))
                mk = month_key(t)
                um = np.unique(mk)
                bs = []
                for _ in range(400):
                    pick = rng.choice(um, len(um))
                    idx = np.concatenate([np.where(mk == m)[0] for m in pick])
                    if np.ptp(t[idx]) < 150:
                        continue
                    bs.append(comp(idx))
                bs = np.array(bs) if len(bs) >= 50 else np.full((1, 3), np.nan)
                rr = rng.uniform(*R_RANGE, len(bs))
                aer_bs = bs[:, 2] / rr
                unexpl_bs = bs[:, 0] - bs[:, 1] - aer_bs
                row = dict(site=site, variant=variant, band=bn, n=len(t), kappa=kap, A_total=a_tot, A_kappaT=a_k,
                           A_cl61_inphase=a_y, A_aerosol=a_y / 0.65,
                           A_aerosol_lo=float(np.percentile(aer_bs, 2.5)), A_aerosol_hi=float(np.percentile(aer_bs, 97.5)),
                           A_gain_unexplained=a_tot - a_k - a_y / 0.65,
                           A_unexpl_lo=float(np.percentile(unexpl_bs, 2.5)), A_unexpl_hi=float(np.percentile(unexpl_bs, 97.5)),
                           f_aer_total=(a_y / 0.65) / a_tot,
                           f_aer_total_lo=float(np.percentile(aer_bs / bs[:, 0], 2.5)),
                           f_aer_total_hi=float(np.percentile(aer_bs / bs[:, 0], 97.5)),
                           f_aer_residual=(a_y / 0.65) / (a_tot - a_k),
                           f_aer_residual_lo=float(np.percentile(aer_bs / (bs[:, 0] - bs[:, 1]), 2.5)),
                           f_aer_residual_hi=float(np.percentile(aer_bs / (bs[:, 0] - bs[:, 1]), 97.5)))
                rows.append(row)
    out = pd.DataFrame(rows)
    out.to_csv(Q17 / "q17b_budget.csv", index=False)
    pd.set_option("display.width", 250)
    print(out.round(2).to_string())


if __name__ == "__main__" and len(sys.argv) > 1 and sys.argv[1] == "budget":
    budget()


def budget_summary():
    """Free-tropospheric (4-8 km) budget per site and the 3-site mean (Payerne excluded)."""
    b = pd.read_csv(Q17 / "q17b_budget.csv")
    b = b[b.band.isin(FT_BANDS)]
    rows = []
    for v, dv in b.groupby("variant"):
        per = []
        for s, d in dv.groupby("site"):
            m = d[["A_total", "A_kappaT", "A_aerosol", "A_aerosol_lo", "A_aerosol_hi", "A_gain_unexplained",
                   "A_unexpl_lo", "A_unexpl_hi"]].mean()
            r = dict(variant=v, site=s, **m.to_dict())
            r["f_aer_total"] = m.A_aerosol / m.A_total
            r["f_aer_residual"] = m.A_aerosol / (m.A_total - m.A_kappaT)
            rows.append(r)
            if s != "Payerne":
                per.append(r)
        mm = pd.DataFrame(per)
        r = dict(variant=v, site="mean_AoLiCa", **mm.drop(columns=["variant", "site"]).mean().to_dict())
        # half-widths shrink by sqrt(3) for the mean of 3 independent sites
        for a, lo, hi in (("A_aerosol", "A_aerosol_lo", "A_aerosol_hi"), ("A_gain_unexplained", "A_unexpl_lo", "A_unexpl_hi")):
            hw_lo = (mm[a] - mm[lo]).mean() / np.sqrt(len(mm))
            hw_hi = (mm[hi] - mm[a]).mean() / np.sqrt(len(mm))
            r[lo], r[hi] = r[a] - hw_lo, r[a] + hw_hi
        r["f_aer_total"] = r["A_aerosol"] / r["A_total"]
        r["f_aer_total_lo"] = r["A_aerosol_lo"] / r["A_total"]
        r["f_aer_total_hi"] = r["A_aerosol_hi"] / r["A_total"]
        res = r["A_total"] - r["A_kappaT"]
        r["f_aer_residual"] = r["A_aerosol"] / res
        r["f_aer_residual_lo"] = r["A_aerosol_lo"] / res
        r["f_aer_residual_hi"] = r["A_aerosol_hi"] / res
        rows.append(r)
    out = pd.DataFrame(rows)
    out.to_csv(Q17 / "q17b_budget_summary.csv", index=False)
    pd.set_option("display.width", 250)
    print(out.round(2).to_string())


if __name__ == "__main__" and len(sys.argv) > 1 and sys.argv[1] == "summary":
    budget_summary()
