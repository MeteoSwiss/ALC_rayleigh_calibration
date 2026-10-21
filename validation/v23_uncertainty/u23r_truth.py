# -*- coding: utf-8 -*-
"""Review attack: is the truth leaking, and do the headline pulls reproduce?

1. Reproduce the doc's headline pulls from the review runs (exact shipped budget, eprof_v2.2).
2. Leakage controls on the time-series truth: exclusion of +-3 / +-7 days around the night (the u23
   truth excludes only the night itself, while a 7-night block shares weather regimes), a causal
   (past-only) window, no kappa removal, the unit's own sampled series instead of the published one.
3. The truth's own error: the 31-day running median has a standard error (1.2533 sigma / sqrt(n)),
   which inflates every anomaly; pulls corrected for it.
Same controls on the twin truth.

Writes u23r_prepared_<set>_<method>.pkl (runs + default truth) and u23r_truth.csv.
"""
from __future__ import annotations

import pickle
import sys
import warnings
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import u23r_common as C  # noqa: E402

warnings.filterwarnings("ignore")


def prepare(set_name, method, T=None, force=False):
    """Runs + default truths (cached)."""
    f = C.OUT / f"u23r_prepared_{set_name}_{method}.pkl"
    if f.exists() and not force:
        with open(f, "rb") as fh:
            return pickle.load(fh)
    T = T or C.Truth()
    df = C.load_runs(set_name, method)
    df["site"] = [T.cen[(w, i)].get("site", w) for w, i in zip(df.wmo, df.ident)]
    df["lat"] = [T.cen[(w, i)]["lat"] for w, i in zip(df.wmo, df.ident)]
    df["alt"] = [T.cen[(w, i)].get("alt", 0.0) for w, i in zip(df.wmo, df.ident)]
    df = C.add_truth(df, T)
    tw = C.twin_frame(df, T) if set_name in ("orig", "origx") else None
    co = None
    if set_name == "fresh":
        co = C.twin_frame(df, T, pairs=C.FRESH_COLOC)
    d = dict(df=df, tw=tw, coloc=co)
    with open(f, "wb") as fh:
        pickle.dump(d, fh)
    return d


def populations(df):
    base = df.ok & np.isfinite(df.anom)
    return {
        "CHM15k (Payerne excl.)": base & (df.type == "CHM15k") & (df.key != C.PAYERNE),
        "Payerne A": base & (df.key == C.PAYERNE),
        "CL61": base & (df.type == "CL61"),
        "Mini-MPL": base & (df.type == "Mini-MPL"),
    }


def twin_pulls(tw, df):
    s22 = C.twin_sigma(tw, df, "sig_ens")
    s23 = C.twin_sigma(tw, df, "sig23", drop_cols=("def_sig_mol", "def_sig_aer"))
    return tw.anom.values / s22, tw.anom.values / s23, s23


def main():
    import pandas as pd
    T = C.Truth()
    d = prepare("orig", "eprof_v2.2", T, force=True)
    df, tw = d["df"], d["tw"]
    rows = []

    def add(variant, pop, m, anom, se=None):
        a = anom[m]
        s22 = C.pull_stats(a / df.sig_ens[m])
        s23 = C.pull_stats(a / df.sig23[m])
        r = dict(variant=variant, population=pop, n=s23["n"], obs_sigma=100 * C.rstd(a),
                 pull22=s22["rstd"], out22=s22["out2"], pull23=s23["rstd"], out23=s23["out2"],
                 med_sig23=100 * np.nanmedian(df.sig23[m]))
        if se is not None:
            z = a / np.sqrt(df.sig23[m] ** 2 + se[m] ** 2)
            z22 = a / np.sqrt(df.sig_ens[m] ** 2 + se[m] ** 2)
            r.update(pull23_secorr=C.pull_stats(z)["rstd"], out23_secorr=C.pull_stats(z)["out2"],
                     pull22_secorr=C.pull_stats(z22)["rstd"], med_se=100 * np.nanmedian(se[m]))
        rows.append(r)

    # --- 1. reproduction (published nights only, as the doc) -------------------------------
    inpub = df.anom_inpub.astype(bool)
    for pop, m in populations(df).items():
        add("u23 truth (published nights)", pop, m & inpub, df.anom, df.anom_se)
    z22, z23, _ = twin_pulls(tw, df)
    rows.append(dict(variant="u23 truth", population="twins", n=int(np.isfinite(z23).sum()),
                     obs_sigma=100 * C.rstd(tw.anom), pull22=C.pull_stats(z22)["rstd"],
                     out22=C.pull_stats(z22)["out2"], pull23=C.pull_stats(z23)["rstd"],
                     out23=C.pull_stats(z23)["out2"]))

    # --- 2. leakage variants --------------------------------------------------------------
    variants = [("exclude +-3 d", dict(exclude_days=3.0)),
                ("exclude +-7 d", dict(exclude_days=7.0)),
                ("causal (past 31 d)", dict(causal=True)),
                ("no kappa removal", dict(use_kappa=False)),
                ("own sampled series", dict(own_series=True)),
                ("own sampled series, +-7 d", dict(own_series=True, exclude_days=7.0))]
    for lab, kw in variants:
        df = C.add_truth(df, T, col="a2", **kw)
        for pop, m in populations(df).items():
            mm = m & np.isfinite(df.a2)
            add(lab, pop, mm, df.a2, df.a2_se)
    for lab, kw in (("exclude +-3 d", dict(exclude_days=3.0)),
                    ("exclude +-7 d", dict(exclude_days=7.0)),
                    ("no kappa removal", dict(use_kappa=False))):
        t2 = C.twin_frame(df, T, **kw)
        s22 = C.twin_sigma(t2, df, "sig_ens")
        s23 = C.twin_sigma(t2, df, "sig23", drop_cols=("def_sig_mol", "def_sig_aer"))
        a, b = C.pull_stats(t2.anom / s22), C.pull_stats(t2.anom / s23)
        rows.append(dict(variant=lab, population="twins", n=b["n"], obs_sigma=100 * C.rstd(t2.anom),
                         pull22=a["rstd"], out22=a["out2"], pull23=b["rstd"], out23=b["out2"]))
    # twin truth error
    s23 = C.twin_sigma(tw, df, "sig23", drop_cols=("def_sig_mol", "def_sig_aer"))
    zc = tw.anom / np.sqrt(s23 ** 2 + tw.anom_se ** 2)
    rows.append(dict(variant="u23 truth, median error removed", population="twins",
                     n=int(np.isfinite(zc).sum()), pull23=C.pull_stats(zc)["rstd"],
                     out23=C.pull_stats(zc)["out2"], med_se=100 * np.nanmedian(tw.anom_se)))
    out = pd.DataFrame(rows)
    out.to_csv(C.OUT / "u23r_truth.csv", index=False)
    pd.set_option("display.width", 250)
    pd.set_option("display.max_columns", 30)
    print(out.round(3).to_string(index=False))


if __name__ == "__main__":
    main()
