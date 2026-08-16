"""Regression du residu d ln C_L / dCBH sur le run natif PVC + CAMS 0,4 deg.

Reproduit l'estimateur legacy (effets fixes flux + bootstrap groupe par flux,
convention C_L) sur les scenes extraites par extract_cbh.py.
"""
import json
import os

import numpy as np
import pandas as pd

SCR = os.path.dirname(os.path.abspath(__file__))
FIGDIR = r"C:/Users/hervo/OneDrive/Documents/ALC_rayleigh_calibration/doc/reports/figs_altitude_audit"
os.makedirs(FIGDIR, exist_ok=True)

MIN_NWIN = 5       # profils nuageux dans la fenetre CBH pour retenir le jour
MIN_DAYS_FLUX = 30  # jours pour la regression par flux

rng = np.random.default_rng(42)


def load():
    df = pd.read_csv(os.path.join(SCR, "cbh_scenes.csv"),
                     dtype={"date": str}, low_memory=False)
    df["cal_value"] = pd.to_numeric(df["cal_value"], errors="coerce")
    df["cbh_med"] = pd.to_numeric(df["cbh_med"], errors="coerce")
    df["n_win"] = pd.to_numeric(df["n_win"], errors="coerce")
    df["frac_below800"] = pd.to_numeric(df["frac_below800"], errors="coerce")
    n_all = len(df)
    n_l1 = (df["l1"] == "1").sum() if df["l1"].dtype == object else (df["l1"] == 1).sum()
    df = df[(df["cal_value"] > 0) & df["cbh_med"].notna() & (df["n_win"] >= MIN_NWIN)]
    df = df.copy()
    df["lnC"] = np.log(df["cal_value"])
    df["cbh_km"] = df["cbh_med"] / 1000.0
    df["month"] = df["date"].str[:6]
    # ruptures instrumentales connues a Payerne: scinder les flux pour que les
    # effets fixes n'absorbent pas un saut de niveau correle a la saison
    m = (df["key"] == "0-20000-0-06610_B") & (df["date"] >= "20260707")
    df.loc[m, "key"] = "0-20000-0-06610_B_post"   # swap bloc optique CL31, C_L x2.8
    m = (df["key"] == "0-20000-0-06610_C") & (df["date"] >= "20260611")
    df.loc[m, "key"] = "0-20000-0-06610_C_post"   # perte regulation thermique diode CL61
    print(f"rows total {n_all}, with L1 {n_l1}, usable (cbh ok, n_win>={MIN_NWIN}) {len(df)}")
    return df


def ols(x, y):
    """slope, se, intercept (OLS classique)."""
    n = len(x)
    if n < 3:
        return np.nan, np.nan, np.nan
    X = np.column_stack([np.ones(n), x])
    beta, res, *_ = np.linalg.lstsq(X, y, rcond=None)
    yhat = X @ beta
    s2 = np.sum((y - yhat) ** 2) / (n - 2)
    cov = s2 * np.linalg.inv(X.T @ X)
    return beta[1], np.sqrt(cov[1, 1]), beta[0]


def fe_slope(df, extra_group=None):
    """Pente pooled apres demeaning par flux (ou flux x extra), bootstrap par flux."""
    g = df["key"] if extra_group is None else df["key"] + "|" + df[extra_group]
    d = df.copy()
    d["g"] = g
    # cellules avec >=3 jours seulement (sinon le demeaning n'apporte rien)
    cnt = d.groupby("g")["lnC"].transform("size")
    d = d[cnt >= 3]
    d["y"] = d["lnC"] - d.groupby("g")["lnC"].transform("mean")
    d["x"] = d["cbh_km"] - d.groupby("g")["cbh_km"].transform("mean")
    if len(d) < 10 or d["x"].std() == 0:
        return np.nan, np.nan, 0, 0
    s = np.sum(d["x"] * d["y"]) / np.sum(d["x"] ** 2)
    # bootstrap groupe par flux (sommes pre-agregees, pente = somme sxy / somme sxx)
    d["_xy"] = d["x"] * d["y"]
    d["_xx"] = d["x"] ** 2
    agg = d.groupby("key")[["_xy", "_xx"]].sum()
    sxy = agg["_xy"].values
    sxx = agg["_xx"].values
    nk = len(agg)
    bs = []
    for _ in range(1000):
        idx = rng.integers(0, nk, size=nk)
        den = sxx[idx].sum()
        if den > 0:
            bs.append(sxy[idx].sum() / den)
    return s, np.std(bs), len(d), nk


def main():
    df = load()
    types = ["CL31", "CL51", "CL61"]
    df = df[df["type"].isin(types + ["Mini-MPL"])]

    out = {}

    # ---- per-type FE estimates -------------------------------------------
    print("\n== Pentes par type (effets fixes flux, bootstrap par flux) ==")
    for t in types:
        d = df[df["type"] == t]
        s, se, n, nk = fe_slope(d)
        sm, sem, nm, nkm = fe_slope(d, "month")
        # bande basse / haute
        sl, sel, nl, _ = fe_slope(d[d["cbh_km"] < 0.8])
        sh, seh, nh, _ = fe_slope(d[d["cbh_km"] >= 0.8])
        out[t] = dict(fe=(s, se, n, nk), fe_month=(sm, sem, nm, nkm),
                      below08=(sl, sel, nl), above08=(sh, seh, nh))
        print(f"{t}: FE {100*s:+.2f} +/- {100*se:.2f} %/km (n={n} j, {nk} flux) | "
              f"FE+mois {100*sm:+.2f} +/- {100*sem:.2f} (n={nm}) | "
              f"CBH<0.8 {100*sl:+.2f} +/- {100*sel:.2f} (n={nl}) | "
              f">=0.8 {100*sh:+.2f} +/- {100*seh:.2f} (n={nh})")

    # ---- robustesse: flag 1.0 seulement ----------------------------------
    print("\n== Robustesse: flag 1.0 seulement ==")
    for t in types:
        d = df[(df["type"] == t) & (df["flag"].astype(str).isin(("1.0", "1")))]
        s, se, n, nk = fe_slope(d)
        print(f"{t}: FE {100*s:+.2f} +/- {100*se:.2f} %/km (n={n} j, {nk} flux)")

    # ---- first-difference (paires de jours proches, <=3 j) ---------------
    print("\n== Estimateur en premieres differences (paires <=3 jours) ==")
    for t in types:
        d = df[df["type"] == t].sort_values(["key", "date"])
        dx, dy, keys = [], [], []
        for k, g in d.groupby("key"):
            dates = pd.to_datetime(g["date"], format="%Y%m%d").values
            x = g["cbh_km"].values
            y = g["lnC"].values
            dt = np.diff(dates).astype("timedelta64[D]").astype(int)
            m = dt <= 3
            dx.append(np.diff(x)[m])
            dy.append(np.diff(y)[m])
            keys.append(np.repeat(k, m.sum()))
        if not dx or all(a.size == 0 for a in dx):
            out.setdefault(t, {})["firstdiff"] = (np.nan, np.nan, 0)
            continue
        dx = np.concatenate(dx); dy = np.concatenate(dy)
        keys = np.concatenate(keys)
        den = np.sum(dx ** 2)
        s = np.sum(dx * dy) / den
        # bootstrap par flux (sommes pre-agregees par flux -> rapide)
        uk = np.unique(keys)
        sxx = np.array([np.sum(dx[keys == k] ** 2) for k in uk])
        sxy = np.array([np.sum(dx[keys == k] * dy[keys == k]) for k in uk])
        bs = []
        for _ in range(1000):
            idx = rng.integers(0, len(uk), size=len(uk))
            d2 = sxx[idx].sum()
            if d2 > 0:
                bs.append(sxy[idx].sum() / d2)
        out[t]["firstdiff"] = (s, np.std(bs), len(dx))
        print(f"{t}: FD {100*s:+.2f} +/- {100*np.std(bs):.2f} %/km (n={len(dx)} paires)")

    # ---- per-flux slopes -------------------------------------------------
    rows = []
    for (k, t), g in df.groupby(["key", "type"]):
        if len(g) < MIN_DAYS_FLUX:
            continue
        s, se, _ = ols(g["cbh_km"].values, g["lnC"].values)
        rows.append(dict(key=k, type=t, slope=100 * s, se=100 * se, n=len(g),
                         cbh_med=g["cbh_km"].median(), cbh_p10=g["cbh_km"].quantile(.1),
                         cbh_p90=g["cbh_km"].quantile(.9)))
    per = pd.DataFrame(rows)
    per.to_csv(os.path.join(SCR, "per_flux_slopes.csv"), index=False)
    print("\n== Par flux (n>=%d) ==" % MIN_DAYS_FLUX)
    for t in types:
        p = per[per["type"] == t]
        pos = (p["slope"] > 0).sum()
        sig = (p["slope"] > 2 * p["se"]).sum()
        print(f"{t}: {len(p)} flux, mediane {p['slope'].median():+.2f} %/km "
              f"(IQR {p['slope'].quantile(.25):+.2f}..{p['slope'].quantile(.75):+.2f}), "
              f"{pos} positifs, {sig} > +2sigma")

    # ---- binned profile per type (normalise par flux) --------------------
    bins = np.array([500, 750, 1000, 1250, 1500, 1750, 2000, 2400], float)
    df["bin"] = pd.cut(df["cbh_med"], bins)
    prof = {}
    for t in types:
        d = df[df["type"] == t].copy()
        d["lnC_n"] = d["lnC"] - d.groupby("key")["lnC"].transform("median")
        gb = d.groupby("bin", observed=True)["lnC_n"]
        prof[t] = pd.DataFrame({"med": gb.median(), "n": gb.size(),
                                "lo": gb.quantile(.25), "hi": gb.quantile(.75)})
        print(f"\n{t} profil par bin CBH (mediane lnC normalisee par flux):")
        print(prof[t].to_string())

    # ---- level bias per flux --------------------------------------------
    # biais d'ancrage: pente type (FE) x (CBH mediane des scenes - reference).
    # On borne aussi par la dispersion P10-P90 de la CBH des scenes.
    print("\n== Biais de niveau implicite (pente FE type x distribution CBH) ==")
    lev = []
    for _, r in per.iterrows():
        s_t = out[r["type"]]["fe"][0]
        lev.append(dict(key=r["key"], type=r["type"],
                        bias_vs_1km=100 * s_t * (r["cbh_med"] - 1.0),
                        span_p10p90=100 * s_t * (r["cbh_p90"] - r["cbh_p10"]),
                        cbh_med_km=r["cbh_med"]))
    lev = pd.DataFrame(lev)
    lev.to_csv(os.path.join(SCR, "level_bias.csv"), index=False)
    for t in types:
        l = lev[lev["type"] == t]
        print(f"{t}: CBH mediane des scenes {l['cbh_med_km'].median():.2f} km "
              f"(range flux {l['cbh_med_km'].min():.2f}..{l['cbh_med_km'].max():.2f}); "
              f"biais vs ref 1 km: mediane {l['bias_vs_1km'].median():+.2f} % "
              f"(range {l['bias_vs_1km'].min():+.2f}..{l['bias_vs_1km'].max():+.2f}); "
              f"span P10-P90: mediane {l['span_p10p90'].median():.2f} %")

    with open(os.path.join(SCR, "fe_results.json"), "w") as fh:
        json.dump({t: {k: (list(v) if isinstance(v, tuple) else v)
                       for k, v in out[t].items()} for t in out}, fh, indent=1)

    # ---- CL61 fleet detail ----------------------------------------------
    sites = {"0-20000-0-03808_C": "Camborne", "0-20000-0-06418_B": "Zeebrugge",
             "0-20000-0-06447_B": "Uccle", "0-20000-0-06610_C": "Payerne",
             "0-20000-0-10393_C": "Lindenberg", "0-20000-0-11538_B": "Temelin",
             "0-20000-0-14015_B": "Ljubljana", "0-20008-0-BIR_A": "Birkenes",
             "0-20008-0-EDT_B": "Edmonton", "0-20008-0-LAU_A": "Lauder",
             "0-203-10-LNG_A": "Lanzhot", "0-380-5-1_B": "Aosta",
             "0-756-4-EERLCL61_A": "Sion"}
    print("\n== Flotte CL61 par flux (n>=15) ==")
    for k, name in sites.items():
        g = df[df["key"] == k]
        if len(g) >= 15:
            s, se, _ = ols(g["cbh_km"].values, g["lnC"].values)
            print(f"{name:12s} {k}: {100*s:+.2f} +/- {100*se:.2f} %/km, n={len(g)}, "
                  f"CBH med {g['cbh_km'].median():.2f} km")

    return df, per, prof, out, lev


if __name__ == "__main__":
    main()
