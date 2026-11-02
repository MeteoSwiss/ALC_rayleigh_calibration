# -*- coding: utf-8 -*-
"""Histograms of the in-cloud extinction of the calibration scenes, two independent ways.

1. From the calibrated ceilometer profiles themselves (the stored in-cloud windows, `incloud/<key>.npz`,
   peak-60..+300 m, one window per accepted calibration profile). With the lidar ratio known
   (S = 18.8 sr) and the profile calibrated by the cloud method (beta_att = beta_before x c_profile /
   T2_wv, in 1/(m sr); its window integral is 0.9 x 1/(2 eta S), the missing 10 % being the cloud above
   the window), the extinction follows from the forward solution
       alpha_i = S x beta_att_i x exp(2 eta_i tau_i),   tau_(i+1) = tau_i + alpha_i dz,
   started at the window bottom (below the cloud onset, tau = 0), with the per-gate multiple-scattering
   factor eta_i the calibration applied. Reported per profile: the mean extinction over the first 50 m and
   100 m above the cloud onset (first gate with beta_att > 20 Mm-1 sr-1), the extinction at the signal
   peak and the optical depth at the peak. Profiles are kept up to tau = 3 (the forward solution amplifies
   noise beyond that).
2. From the Cloudnet microphysics matched to the same scenes (step 1d, `step1d/scenes_msfov.csv`):
   alpha = 2 pi N_scaled a_G^2 (Q_ext = 2) in the lowest 50 / 100 / 150 m above the liquid base
   (`alpha_km_<D>`), scenes with base_status 1.

Units: Payerne CL61 (C) and CL31 (B), Palaiseau CL31, Lindenberg CL61 (the step-1d sites). Also listed:
the per-unit tail-decay alpha_k of step 2 (`summary/unit_summary.csv`, a lower bound).

Outputs: C:/DATA/Projects/202611_eta_study/extinction/ (summary.txt, extinction_per_profile.csv.gz) and
doc/reports/figs_04v2/fig_extinction_histograms.png.

Run:  python validation/eta_study/extinction_histograms.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO))
STUDY = Path("C:/DATA/Projects/202611_eta_study")
OUT = STUDY / "extinction"
S_THEO = 18.8
ONSET = 20e-6            # 1/(m sr): cloud onset threshold on the calibrated attenuated backscatter
TAU_MAX = 3.0
UNITS = [("0-20000-0-06610_C", "Payerne CL61"), ("0-20000-0-06610_B", "Payerne CL31 (distorted unit)"),
         ("0-250-1001-07151_A", "Palaiseau CL31"), ("0-20000-0-10393_C", "Lindenberg CL61")]


def forward_extinction(key):
    z = np.load(STUDY / "incloud" / f"{key}.npz", allow_pickle=True)
    B = z["beta_before"].astype("f8")
    T2 = z["t2_wv"].astype("f8")
    eta = z["eta"].astype("f8")
    c = z["c_profile"].astype("f8")
    dz = float(z["dr_m"])
    rr = z["range_rel_m"].astype("f8")
    beta = B * c[:, None] / np.where(T2 > 0, T2, np.nan)              # calibrated beta_att, 1/(m sr)
    n, ng = beta.shape
    alpha = np.full_like(beta, np.nan)
    tau = np.zeros(n)
    tau_prof = np.zeros((n, ng))
    for i in range(ng):
        b = beta[:, i]
        a = S_THEO * np.where(np.isfinite(b) & (b > 0), b, 0.0) * np.exp(2.0 * eta[:, i] * tau)
        alpha[:, i] = a
        tau_prof[:, i] = tau
        tau = tau + a * dz
    # onset: first gate above the threshold; the peak gate is where rr == 0 (window is peak-relative)
    above = beta > ONSET
    has = above.any(axis=1)
    onset = np.where(has, above.argmax(axis=1), -1)
    ipeak = int(np.argmin(np.abs(rr)))
    rows = []
    for k in range(n):
        if onset[k] < 0 or not np.isfinite(c[k]) or c[k] <= 0:
            continue
        i0 = onset[k]
        i50 = i0 + int(round(50.0 / dz))
        i100 = i0 + int(round(100.0 / dz))
        if i100 >= ng:
            continue
        t = tau_prof[k]
        ok50 = t[i50] <= TAU_MAX
        ok100 = t[i100] <= TAU_MAX
        rows.append(dict(key=key, date=int(z["date"][k]), cbh_m=float(z["cbh_m"][k]),
                         onset_rel_m=float(rr[i0]),
                         alpha50=1e3 * (t[i50] - t[i0]) / (dz * (i50 - i0)) if ok50 else np.nan,
                         alpha100=1e3 * (t[i100] - t[i0]) / (dz * (i100 - i0)) if ok100 else np.nan,
                         alpha_peak=1e3 * alpha[k, ipeak] if t[ipeak] <= TAU_MAX else np.nan,
                         tau_peak=float(t[ipeak]), beta_peak=float(beta[k, ipeak]) * 1e6,
                         eta_onset=float(eta[k, i0])))
    return pd.DataFrame(rows)


def main():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    OUT.mkdir(parents=True, exist_ok=True)
    sc = pd.read_csv(STUDY / "step1d" / "scenes_msfov.csv", low_memory=False)
    sc = sc[(sc.status == "matched") & (sc.base_status == 1)]
    us = pd.read_csv(STUDY / "summary" / "unit_summary.csv").set_index("key")
    frames = []
    lines = ["unit | lidar forward alpha: 0-50 m, 0-100 m above onset, at peak [/km, median (P10-P90)], tau at peak, n | "
             "Cloudnet alpha 50 / 100 / 150 m [/km, median (P10-P90)], n | tail-decay alpha_k (step 2)"]
    fig, axes = plt.subplots(2, 4, figsize=(20, 9), sharex=True)
    bins = np.logspace(0, 2.7, 46)
    for j, (key, label) in enumerate(UNITS):
        df = forward_extinction(key)
        df["label"] = label
        frames.append(df)
        q = lambda s: f"{np.nanmedian(s):.1f} ({np.nanpercentile(s.dropna(), 10):.1f}-{np.nanpercentile(s.dropna(), 90):.1f})"
        cn = sc[sc.key == key]
        ak = us.alpha_k_med.get(key, np.nan)
        lines.append(f"{label}: {q(df.alpha50)}, {q(df.alpha100)}, {q(df.alpha_peak)}, tau_peak {np.nanmedian(df.tau_peak):.2f}, n {len(df)} | "
                     f"{q(cn.alpha_km_50)} / {q(cn.alpha_km_100)} / {q(cn.alpha_km_150)}, n {len(cn)} | {ak:.1f}")
        ax = axes[0, j]
        for col, c_, lab in (("alpha50", "#1f77b4", "0-50 m above onset"), ("alpha100", "#2ca02c", "0-100 m above onset"),
                             ("alpha_peak", "#7f7f7f", "at the signal peak")):
            v = df[col].dropna()
            v = v[v > 0]
            ax.hist(v, bins, histtype="step", lw=2, color=c_, density=True, label=f"{lab}: med {np.median(v):.0f} /km (n={len(v)})")
        ax.axvspan(15, 20, color="#f0c000", alpha=0.25, lw=0, label="O'Connor 2004 assumed 15-20")
        for a_ in (10, 20, 30):
            ax.axvline(a_, color="k", ls=":", lw=1)
        if np.isfinite(ak):
            ax.axvline(ak, color="#d62728", ls="--", lw=1.5, label=f"tail decay alpha_k {ak:.0f} (lower bound)")
        ax.set_xscale("log")
        ax.set_title(f"{label}: forward inversion of the calibrated profile", fontsize=10)
        ax.grid(alpha=0.3)
        ax.legend(fontsize=7, loc="upper left")
        ax = axes[1, j]
        for D, c_ in (("50", "#9467bd"), ("100", "#2ca02c"), ("150", "#ff7f0e")):
            v = cn[f"alpha_km_{D}"].dropna()
            v = v[v > 0]
            if len(v):
                ax.hist(v, bins, histtype="step", lw=2, color=c_, density=True, label=f"lowest {D} m: med {np.median(v):.0f} /km (n={len(v)})")
        ax.axvspan(15, 20, color="#f0c000", alpha=0.25, lw=0)
        for a_ in (10, 20, 30):
            ax.axvline(a_, color="k", ls=":", lw=1)
        ax.set_xscale("log")
        ax.set_title(f"{label}: Cloudnet 2 pi N a_G^2, same scenes", fontsize=10)
        ax.set_xlabel("extinction alpha [1/km]")
        ax.grid(alpha=0.3)
        ax.legend(fontsize=7, loc="upper left")
    axes[0, 0].set_ylabel("probability density (log bins)")
    axes[1, 0].set_ylabel("probability density (log bins)")
    fig.suptitle("In-cloud extinction of the liquid-cloud calibration scenes: ceilometer forward inversion (top) and Cloudnet microphysics (bottom); "
                 "dotted = the 10 / 20 / 30 /km eta-table ladder", fontsize=11)
    fig.tight_layout()
    fig.savefig(REPO / "doc" / "reports" / "figs_04v2" / "fig_extinction_histograms.png", dpi=130)
    plt.close(fig)
    allp = pd.concat(frames, ignore_index=True)
    allp.to_csv(OUT / "extinction_per_profile.csv.gz", index=False, float_format="%.5g")
    txt = "\n".join(lines)
    print(txt)
    (OUT / "summary.txt").write_text(txt, encoding="utf-8")
    print("->", REPO / "doc" / "reports" / "figs_04v2" / "fig_extinction_histograms.png")


if __name__ == "__main__":
    main()
