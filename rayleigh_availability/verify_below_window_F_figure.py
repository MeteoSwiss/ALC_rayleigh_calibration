# -*- coding: utf-8 -*-
"""Figure for the adversarial verification of mechanism (b) -- below-window transmission.

Three landscape panels, altitude on Y wherever altitude is plotted:

  1. The aerosol the mechanism REQUIRES versus the aerosol CAMS actually carries at Payerne on the
     v2.0-kept and v2.2-recovered nights (extinction profile, m^-1, altitude AGL on Y).
  2. The C_L bias ladder versus molecular-window centre (altitude AGL on Y) for the required scene
     and for the real median recovered night -- the signature test (level offset, zero slope) and
     the amplitude test in one panel.
  3. The (S_true, AOD_1064) a site would need to hit its published offset, with the CAMS
     per-night cloud of the same two quantities overlaid.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

REPO = Path("C:/Users/hervo/OneDrive/Documents/ALC_rayleigh_calibration")
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / "rayleigh_availability"))

from forward_model import (  # noqa: E402
    PAYERNE_CHM15K, default_grid, default_options, forward_signal, make_atmosphere, retrieve,
)
from forward_scan_transmission import scene  # noqa: E402

HERE = Path(__file__).parent   # verif_*.json live beside these scripts
FIGDIR = REPO / "doc" / "reports" / "figs_altitude_audit"
FIG = FIGDIR / "fwd_below_window_transmission_reality_check.png"
HALF = 490.0


def ladder(z, atm, C_true, options, aod, s_true, centres):
    ba, ea = scene(z, "below", aod=aod, lidar_ratio_true=s_true, top_m=1500.0)
    rcs = forward_signal(z, C_true, atm["beta_mol"], ba, ea)
    out = []
    for c in centres:
        r = retrieve(z, rcs, atm, options=options, window=(c - HALF, c + HALF))
        out.append((r["C_L"] / C_true - 1.0) * 100.0 if r["ok"] else np.nan)
    return np.array(out), ba, ea


def main():
    B = json.loads((HERE / "verif_B.json").read_text(encoding="utf-8"))
    C = json.loads((HERE / "verif_C.json").read_text(encoding="utf-8"))
    E = json.loads((HERE / "verif_E.json").read_text(encoding="utf-8"))

    z = default_grid()
    C_true = PAYERNE_CHM15K["C_true"]
    atm = make_atmosphere(z, station_alt_m=PAYERNE_CHM15K["station_alt_m"],
                          wavelength_nm=PAYERNE_CHM15K["wavelength_nm"])
    options = default_options()
    centres = np.arange(2500.0, 6501.0, 250.0)

    pay = B["PAYERNE_CHM15k_A"]
    med_rec = pay["summary"]["recovered"]["aod_below1500"]["median"]
    max_rec = pay["summary"]["recovered"]["aod_below1500"]["max"]
    med_kept = pay["summary"]["kept"]["aod_below1500"]["median"]
    s_rec = pay["summary"]["recovered"]["S_local_0_1km"]["median"]

    lad_req, ba_req, ea_req = ladder(z, atm, C_true, options, 0.38, 90.0, centres)
    lad_real, ba_real, ea_real = ladder(z, atm, C_true, options, med_rec, s_rec, centres)
    _, _, ea_max = ladder(z, atm, C_true, options, max_rec, 70.8, centres)
    _, _, ea_kept = ladder(z, atm, C_true, options, med_kept, 36.8, centres)

    fig, axes = plt.subplots(1, 3, figsize=(17.5, 6.6))

    # --- panel 1: extinction profiles -----------------------------------------------------------
    ax = axes[0]
    ax.plot(ea_req, z / 1000.0, color="C3", lw=2.4,
            label=f"REQUIRED by the claim\nAOD$_{{1064}}$ = 0.38, S = 90 sr")
    ax.plot(ea_max, z / 1000.0, color="C1", lw=1.6, ls="--",
            label=f"CAMS worst recovered night\nAOD = {max_rec:.4f}, S = 70.8 sr")
    ax.plot(ea_real, z / 1000.0, color="C0", lw=2.0,
            label=f"CAMS median recovered night\nAOD = {med_rec:.4f}, S = {s_rec:.0f} sr")
    ax.plot(ea_kept, z / 1000.0, color="0.45", lw=1.6, ls=":",
            label=f"CAMS median kept night\nAOD = {med_kept:.4f}, S = 36.8 sr")
    ax.set_xscale("log")
    ax.set_xlim(1e-7, 1e-3)
    ax.set_ylim(0, 3.0)
    ax.axhspan(2.0, 3.0, color="C2", alpha=0.08)
    ax.text(1.3e-7, 2.5, "bottom of the\nfit band", fontsize=8, color="C2")
    ax.set_xlabel(r"aerosol extinction $\alpha_{aer}$ at 1064 nm   [m$^{-1}$]")
    ax.set_ylabel("range above the instrument   [km AGL]")
    ax.set_title("1. The below-window layer:\nrequired vs what CAMS carries (Payerne)",
                 fontsize=11)
    ax.grid(alpha=0.3)
    ax.legend(loc="upper right", fontsize=7.6)

    # --- panel 2: the C_L ladder ---------------------------------------------------------------
    ax = axes[1]
    ax.axvline(0.0, color="0.4", lw=1.2, ls="-")
    ax.axvline(-22.5, color="C3", lw=1.8, ls="--",
               label="observed Payerne CHM15k\nrecovered $-$ kept = $-$22.5 %")
    ax.plot(lad_req, centres / 1000.0, "o-", color="C3", ms=5, lw=1.8,
            label=f"required scene (AOD 0.38, 90 sr)\nmean {np.nanmean(lad_req):+.1f} %, "
                  f"slope {np.polyfit(centres/1000, lad_req, 1)[0]:+.3f} %/km")
    ax.plot(lad_real, centres / 1000.0, "s-", color="C0", ms=5, lw=1.8,
            label=f"CAMS median recovered night\nmean {np.nanmean(lad_real):+.2f} %, "
                  f"slope {np.polyfit(centres/1000, lad_real, 1)[0]:+.3f} %/km")
    ax.set_xlim(-30, 10)
    ax.set_xlabel(r"$C_L / C_{true} - 1$   [%]")
    ax.set_ylabel("molecular-window centre   [km AGL]")
    ax.set_title("2. The signature: a pure LEVEL offset\n(zero slope -- it cannot make a gradient)",
                 fontsize=11)
    ax.grid(alpha=0.3)
    ax.legend(loc="lower left", fontsize=7.6)
    ax.text(0.97, 0.97,
            "closure of the shipped chain\non pure molecular air: 0.0046 %\n"
            "independent Fernald agrees to 0.23 pp",
            transform=ax.transAxes, fontsize=7.6, va="top", ha="right",
            bbox=dict(boxstyle="round", fc="white", ec="0.7", alpha=0.9))

    # --- panel 3: requirement vs the observed (S, AOD) cloud ------------------------------------
    ax = axes[2]
    colors = {"PAYERNE_CHM15k_A": "C3", "LINDENBERG_CHM15k_0": "C0", "PALAISEAU_CHM15k_B": "C2"}
    names = {"PAYERNE_CHM15k_A": "Payerne", "LINDENBERG_CHM15k_0": "Lindenberg",
             "PALAISEAU_CHM15k_B": "SIRTA/Palaiseau"}
    for label, col in colors.items():
        rows = B[label]["rows"]["recovered"]
        s = np.array([r["S_local_0_1km"] for r in rows], float)
        a = np.array([r["aod_below1500"] for r in rows], float)
        ax.scatter(s, a, s=13, color=col, alpha=0.55, edgecolor="none",
                   label=f"{names[label]} recovered nights (CAMS, n={np.isfinite(s).sum()})")
        req = C["required"][{"PAYERNE_CHM15k_A": "Payerne CHM15k",
                             "LINDENBERG_CHM15k_0": "Lindenberg CHM15k",
                             "PALAISEAU_CHM15k_B": "SIRTA CHM15k"}[label]]["required_aod_by_S"]
        ss = np.array(sorted(float(k) for k in req), float)
        aa = np.array([req[f"{v:.1f}"] if f"{v:.1f}" in req else req[str(v)] for v in ss], float)
        ax.plot(ss, aa, "-", color=col, lw=2.2,
                label=f"{names[label]}: AOD needed for {B[label]['observed_offset_pct']:+.1f} %")
    ax.set_yscale("log")
    ax.set_ylim(2e-3, 2.0)
    ax.set_xlim(10, 115)
    ax.axvline(52.0, color="0.3", lw=1.4, ls="--")
    ax.text(52.8, 1.2, "S assumed by\nthe retrieval\n= 52 sr", fontsize=8, color="0.3")
    ax.set_xlabel(r"true aerosol lidar ratio $S_{true}$ at 1064 nm   [sr]")
    ax.set_ylabel(r"AOD$_{1064}$ below 1500 m AGL   [-]")
    ax.set_title("3. The requirement sits far outside\nthe atmosphere CAMS reports", fontsize=11)
    ax.grid(alpha=0.3, which="both")
    ax.legend(loc="lower right", fontsize=6.8, ncol=1, framealpha=0.92)

    pred = E["PAYERNE_CHM15k_A"]["predicted_recovered_minus_kept_pct"]
    fig.suptitle(
        "Mechanism (b), below-window transmission with a lidar-ratio error: the simulation is "
        "right, the plausibility verdict is not.\n"
        f"Fed the real CAMS aerosol, it predicts recovered $-$ kept = {pred:+.3f} % at Payerne "
        f"against the observed $-$22.5 % (shortfall x"
        f"{E['PAYERNE_CHM15k_A']['shortfall_factor']:.0f}), "
        "and CAMS puts Payerne's S at 37 sr < 52 sr, i.e. the WRONG SIGN.",
        fontsize=11, y=0.985)
    fig.tight_layout(rect=(0, 0, 1, 0.925))
    FIGDIR.mkdir(parents=True, exist_ok=True)
    fig.savefig(FIG, dpi=140)
    plt.close(fig)
    print(f"-> {FIG}")


if __name__ == "__main__":
    main()
