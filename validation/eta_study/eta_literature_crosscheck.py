"""Cross-check of the published liquid-cloud multiple-scattering factors against the repo's PVC port.

Companion script of doc/reports/eta_alpha_literature_2026-11.md (section 5). It runs the literature
cloud inputs through validation.multiple_scattering_eta (the line-for-line port of Hogan's multiscatter
1.2.11 'original' PVC algorithm) and prints four tables:

A. O'Connor et al. (2004): CT75K optics (half-angle FOV 0.66 mrad, divergence 0.75 mrad, 905 nm, their
   Sect. 2), monodisperse droplet diameters 8 / 20 um (a_G = 4 / 10 um) and homogeneous extinction
   15 / 20 /km, with the optics scaled by 1, 1/sqrt(2) and 1/2, against their published eta
   (0.83 +- 0.09 at 1 km, 0.79 at 1.5 km, 0.74 +- 0.04 at 2.7 km, 0.73 +- 0.06 at 4 km; pp. 780-782).
   The swapped optics of their Sect. 3 (FOV 0.75, divergence 0.66 mrad) are run as well.
B. Le et al. (2026, AMT 19, 4923): the two eta*S curves of their supplement Fig. S8 (CL61, droplet
   effective diameters 8 and 20 um; digitised, see DIGITISED below), against PVC with the CL61 optics
   for alpha 5..30 /km and the optics scaled by 1, 1/sqrt(2), 1/2 (joint rms over both curves).
C. The operational CL61 tables of calibration/cloud/_filters.py (v2, pvc55, legacy ladder) against the
   same digitised envelope.
D. Our inputs (a_G 3.3 um, type optics) at alpha 10 / 15 / 20 /km against the adopted 30 /km: the
   C_L change (C_L scales with eta at fixed S).
E. The extinction just above cloud base of an adiabatic cloud (note section 2.9): alpha(z) for droplet
   numbers N 50..500 /cm3 and cloud-base temperatures 0 / 10 / 20 degC (900 hPa), the heights where the
   optical depth reaches 1 and 3, and the mean alpha over the tau <= 3 layer seen by a ground-based lidar.

None of this chooses an input: it only shows what the literature numbers correspond to in our model.

Digitising (only needed to re-derive DIGITISED): the supplement of Le et al. (2026),
https://amt.copernicus.org/articles/19/4923/2026/amt-19-4923-2026-supplement.pdf, page 6, holds Fig. S8
as one 9667 x 5469 raster. Save it (e.g. with PyMuPDF: Pixmap(doc, xref).save) and pass it with
``--digitise s8.png``: the 2023-12 panel (row 4, column 2) has few data points, so the red curves are
clean; x runs 0..20 sr between the inner spine edges (the midline gridline at x = 10 checks it), y is
set by the gridlines at 2000 and 4000 m. eta = (eta*S) / 18.8 sr.

Run: python validation/eta_study/eta_literature_crosscheck.py [--digitise s8.png]
         [--fig doc/reports/figs_eta_literature/fig_eta_literature.png]      (a few seconds)
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from calibration.cloud._filters import (  # noqa: E402
    _ETA_CL61, _ETA_LEGACY_CL51, _ETA_PVC55_CL61)
from validation.multiple_scattering_eta import (  # noqa: E402
    OPTICS, WAVELENGTH, eta_for_cloud_base)

S_LIQUID = 18.8          # sr, the lidar ratio behind every eta*S curve compared here
OPTICS_SCALES = (1.0, 1.0 / np.sqrt(2.0), 0.5)

# Le et al. (2026) supplement Fig. S8, panel 2023-12, digitised 2026-11-02 (eta*S in sr; the lower curve
# is the larger droplets, 20 um, i.e. the stronger multiple scattering). Reading error ~0.05 sr.
DIGITISED_Z_M = np.array([600.0, 1000.0, 1500.0, 2000.0, 3000.0, 4000.0])
DIGITISED_LOWER = np.array([16.55, 15.57, 14.70, 14.04, 13.09, 12.45])
DIGITISED_UPPER = np.array([18.00, 17.51, 16.98, 16.53, 15.79, 15.19])

# O'Connor et al. (2004): published CT75K eta (Sect. 2 p. 780, Sect. 3 p. 780, Sect. 4 p. 782).
OCONNOR_PUBLISHED = {1.0: "0.83 +- 0.09", 1.5: "0.79", 2.7: "0.74 +- 0.04 (printed +-0.4)",
                     4.0: "0.73 +- 0.06"}


def lobe_width(wavelength_m: float, a_g_um: float) -> float:
    """Forward-lobe 1/e half-width Theta = lambda / (pi a_G) (Hogan 2006, Eq. 5)."""
    return wavelength_m / (np.pi * a_g_um * 1e-6)


def eta_curve(z_m, rho_t, rho_l, wavelength_m, a_g_um, alpha_per_km):
    th = lobe_width(wavelength_m, a_g_um)
    return np.array([eta_for_cloud_base(z, rho_t, rho_l, th, alpha=alpha_per_km * 1e-3) for z in z_m])


def table_a_oconnor() -> None:
    print("A. O'Connor et al. (2004) cloud set (D 8/20 um x alpha 15/20 /km), CT75K, 905 nm")
    z = np.array(sorted(OCONNOR_PUBLISHED)) * 1e3
    runs = [(f"Sect. 2 optics x{s:.2f} (FOV {0.66 * s:.2f}, div {0.75 * s:.2f} mrad)", 0.66e-3 * s, 0.75e-3 * s)
            for s in OPTICS_SCALES]
    runs.append(("Sect. 3 optics (FOV 0.75, div 0.66 mrad)", 0.75e-3, 0.66e-3))
    for label, rho_t, rho_l in runs:
        cases = np.array([eta_curve(z, rho_t, rho_l, 905e-9, a, al) for a in (4.0, 10.0) for al in (15, 20)])
        print(f"  {label}")
        for k, cbh in enumerate(sorted(OCONNOR_PUBLISHED)):
            lo, hi = cases[:, k].min(), cases[:, k].max()
            print(f"    CBH {cbh:3.1f} km: PVC {lo:.3f}-{hi:.3f} (mid {0.5 * (lo + hi):.3f}; D8/a15 {cases[0, k]:.3f})"
                  f"   published {OCONNOR_PUBLISHED[cbh]}")


def table_b_le_fit(lower, upper) -> None:
    print("\nB. Le et al. (2026) Fig. S8 envelope vs PVC, CL61 optics (joint rms of eta over both curves)")
    rho_t, rho_l = OPTICS["CL61"]
    lam = WAVELENGTH["CL61"]
    obs = {20.0: lower / S_LIQUID, 8.0: upper / S_LIQUID}
    for s in OPTICS_SCALES:
        for alpha in (5, 10, 15, 20, 30):
            rms = {d: np.sqrt(np.mean((eta_curve(DIGITISED_Z_M, rho_t * s, rho_l * s, lam, d / 2.0, alpha)
                                       - o) ** 2)) for d, o in obs.items()}
            print(f"  optics x{s:.2f} (FOV {rho_t * s * 1e3:.2f}, div {rho_l * s * 1e3:.2f} mrad) alpha {alpha:2d}/km: "
                  f"rms D20 {rms[20.0]:.3f}  D8 {rms[8.0]:.3f}  joint {np.hypot(rms[20.0], rms[8.0]):.3f}")


def interp_table(tab: np.ndarray, z_km: np.ndarray) -> np.ndarray:
    """Linear interpolation with linear extrapolation above the last node (as the operational code)."""
    out = np.interp(z_km, tab[:, 0], tab[:, 1])
    above = z_km > tab[-1, 0]
    slope = (tab[-1, 1] - tab[-2, 1]) / (tab[-1, 0] - tab[-2, 0])
    out[above] = tab[-1, 1] + (z_km[above] - tab[-1, 0]) * slope
    return out


def table_c_operational(lower, upper) -> None:
    print("\nC. Operational CL61 tables vs the Le et al. (2026) envelope (eta; d = table - lower edge)")
    z_km = DIGITISED_Z_M / 1e3
    for name, tab in (("v2 (a_G 3.3 um, 30 /km)", _ETA_CL61), ("pvc55 (a_G 5.5 um, 10 /km)", _ETA_PVC55_CL61),
                      ("legacy ladder (CL51 rows)", _ETA_LEGACY_CL51)):
        eta = interp_table(tab, z_km)
        cells = "  ".join(f"{z:.1f} km {e:.3f} [{lo / S_LIQUID:.3f}-{hi / S_LIQUID:.3f}] d {e - lo / S_LIQUID:+.3f}"
                          for z, e, lo, hi in zip(z_km, eta, lower, upper))
        print(f"  {name:27s} {cells}")


def table_d_alpha_sensitivity() -> None:
    print("\nD. a_G 3.3 um, type optics: eta at alpha 10/15/20 /km vs 30 /km (C_L change at fixed S)")
    for inst in ("CL31", "CL51", "CL61"):
        rho_t, rho_l = OPTICS[inst]
        for cbh in (0.5, 1.0, 2.0):
            eta = {al: eta_curve([cbh * 1e3], rho_t, rho_l, WAVELENGTH[inst], 3.3, al)[0] for al in (10, 15, 20, 30)}
            cells = "  ".join(f"alpha {al}: {100 * (eta[al] / eta[30] - 1):+5.1f} %" for al in (20, 15, 10))
            print(f"  {inst} CBH {cbh:.1f} km: eta(30) {eta[30]:.4f}  {cells}")


def digitise(png: Path) -> tuple[np.ndarray, np.ndarray]:
    """Re-derive the two eta*S curves from the Fig. S8 raster (see module docstring)."""
    from PIL import Image

    im = np.asarray(Image.open(png).convert("RGB")).astype(int)
    h, w, _ = im.shape
    p = im[int(3 * h / 6):int(4 * h / 6), int(w / 6):int(2 * w / 6)]
    r, g, b = p[..., 0], p[..., 1], p[..., 2]
    black = (r < 60) & (g < 60) & (b < 60)
    gray = (abs(r - g) < 8) & (abs(g - b) < 8) & (r > 150) & (r < 205)
    red = (r > 200) & (g < 80) & (b < 80)
    cols = np.where(black.sum(axis=0) > 0.5 * black.sum(axis=0).max())[0]
    rows = np.where(black.sum(axis=1) > 0.5 * black.sum(axis=1).max())[0]
    mid_c, mid_r = 0.5 * (cols.min() + cols.max()), 0.5 * (rows.min() + rows.max())
    # inner spine edges; the crop also holds the colour-bar frame, hence the split at the middle
    x0, x1 = cols[cols < mid_c].max(), cols[cols > mid_c].min()
    y0, y1 = rows[rows < mid_r].max(), rows[rows > mid_r].min()
    inner = gray[y0 + 3:y1 - 2, x0 + 3:x1 - 2]
    grid = np.where(inner.sum(axis=1) > 0.4 * inner.shape[1])[0] + y0 + 3
    lines = [grp.mean() for grp in np.split(grid, np.where(np.diff(grid) > 2)[0] + 1) if grp.size]
    y4000, y2000 = sorted(lines)[:2]
    m_per_px = 2000.0 / (y2000 - y4000)
    lower, upper = [], []
    for z in DIGITISED_Z_M:
        row = int(round(y2000 + (2000.0 - z) / m_per_px))
        idx = np.where(red[row - 1:row + 2].any(axis=0))[0]
        segs = [s for s in np.split(idx, np.where(np.diff(idx) > 3)[0] + 1) if s.size >= 2]
        xs = sorted((s.mean() - x0) / (x1 - x0) * 20.0 for s in segs)
        lower.append(xs[-2])   # the two right-most red segments are the curves (the label is left)
        upper.append(xs[-1])
    return np.array(lower), np.array(upper)


def adiabatic_lwc_rate(t_base_k: float, p_hpa: float = 900.0) -> float:
    """Adiabatic increase of liquid water content with height at cloud base [g m-3 per m].

    Condensation along the moist adiabat: dq_l/dz = -dq_s/dz with dT/dz = -Gamma_m and the hydrostatic
    pressure drop; Bolton (1980) saturation vapour pressure. About 1.6 / 2.1 / 2.6 g m-3 km-1 at 0 / 10 /
    20 degC and 900 hPa.
    """
    g, rd, rv, cp, lv, eps = 9.81, 287.0, 461.5, 1004.0, 2.5e6, 0.622
    tc = t_base_k - 273.15
    es = 6.112 * np.exp(17.67 * tc / (tc + 243.5))
    qs = eps * es / (p_hpa - (1.0 - eps) * es)
    gamma_m = g * (1.0 + lv * qs / (rd * t_base_k)) / (cp + lv ** 2 * qs / (rv * t_base_k ** 2))
    dqs_dz = qs * (-(lv / (rv * t_base_k ** 2)) * gamma_m + g / (rd * t_base_k))
    rho_air = p_hpa * 100.0 / (rd * t_base_k)
    return -dqs_dz * rho_air * 1e3


def adiabatic_alpha_coeff(n_per_cm3: float, t_base_k: float, k: float = 0.8, f_ad: float = 1.0) -> float:
    """A in alpha(z) = A z^(2/3) [m^-5/3] for an adiabatic cloud with constant droplet number N.

    alpha = (3/2) LWC / (rho_w r_e), LWC = f_ad Gamma_l z, r_e = r_v k^(-1/3), r_v from LWC and N
    (Q_ext = 2; k = (r_v / r_e)^3, 0.8 for continental spectra, exp(-3 sigma^2) = 0.69 for a lognormal 0.35).
    """
    rho_w = 1e6                                   # g m-3
    n = n_per_cm3 * 1e6                           # m-3
    lwc_rate = f_ad * adiabatic_lwc_rate(t_base_k)   # g m-3 per m
    return 1.5 * (4.0 * np.pi * k / 3.0) ** (1.0 / 3.0) * rho_w ** (-2.0 / 3.0) * n ** (1.0 / 3.0) \
        * lwc_rate ** (2.0 / 3.0)


def table_e_adiabatic() -> None:
    print("\nE. Adiabatic cloud (k 0.8, f_ad 1, 900 hPa): alpha at 25/50/100 m above base, heights of tau = 1 and 3,"
          " mean alpha over tau <= 3")
    for t_c in (0.0, 10.0, 20.0):
        t_k = t_c + 273.15
        print(f"  cloud base {t_c:4.0f} degC: dLWC/dz = {adiabatic_lwc_rate(t_k) * 1e3:.2f} g m-3 km-1")
        for n in (50, 100, 200, 300, 500):
            a = adiabatic_alpha_coeff(n, t_k)
            alpha_at = {z: a * z ** (2.0 / 3.0) * 1e3 for z in (25.0, 50.0, 100.0)}
            z_tau = {tau: (5.0 * tau / (3.0 * a)) ** 0.6 for tau in (1.0, 3.0)}
            print(f"    N {n:3d}/cm3: alpha {alpha_at[25.0]:5.1f} / {alpha_at[50.0]:5.1f} / {alpha_at[100.0]:5.1f} /km;"
                  f" tau 1 at {z_tau[1.0]:4.0f} m, tau 3 at {z_tau[3.0]:4.0f} m;"
                  f" mean alpha over tau <= 3: {3.0 / z_tau[3.0] * 1e3:5.1f} /km")
    print("  (alpha scales as k^(1/3) f_ad^(2/3); the tau <= 3 mean as N^(1/5) (f_ad dLWC/dz)^(2/5) k^(1/5))")


def pvc_band(z_m, rho_t, rho_l, wavelength_m, sizes_um, alphas):
    """Min / max eta over the (droplet diameter x alpha) cases, as the literature envelopes are drawn."""
    cases = np.array([eta_curve(z_m, rho_t, rho_l, wavelength_m, d / 2.0, al) for d in sizes_um for al in alphas])
    return cases.min(axis=0), cases.max(axis=0)


def make_figure(path: Path, lower, upper) -> None:
    """Three landscape panels, cloud-base height on the y axis (section 5 of the note)."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    blue, orange, aqua, yellow = "#2a78d6", "#eb6834", "#1baf7a", "#eda100"   # fixed categorical order
    ink, muted = "#0b0b0b", "#52514e"
    # constrained layout keeps the legends, placed under each panel so they never cover data, inside the figure
    fig, axes = plt.subplots(1, 3, figsize=(17, 6.6), sharey=True, layout="constrained")
    z = np.arange(0.25, 4.01, 0.125) * 1e3
    pvc_runs = ((1.0, blue, None, "PVC, documented optics"), (0.5, orange, "//", "PVC, FOV and divergence × 0.5"))

    ax = axes[0]
    for scale, colour, hatch, label in pvc_runs:
        lo, hi = pvc_band(z, 0.66e-3 * scale, 0.75e-3 * scale, 905e-9, (8.0, 20.0), (15, 20))
        ax.fill_betweenx(z / 1e3, lo, hi, color=colour, alpha=0.22, hatch=hatch, edgecolor=colour, lw=0)
        ax.plot(lo, z / 1e3, color=colour, lw=1.5)
        ax.plot(hi, z / 1e3, color=colour, lw=1.5, label=label)
    pub_z = np.array([1.0, 1.5, 2.7, 4.0])
    pub_eta = np.array([0.83, 0.79, 0.74, 0.73])
    pub_err = np.array([0.09, 0.0, 0.04, 0.06])
    ax.errorbar(pub_eta, pub_z, xerr=pub_err, fmt="o", ms=7, color=ink, capsize=4, lw=1.2,
                label="published η (Eloranta model), O'Connor et al. 2004")
    ax.set_title("(a) CT75K, O'Connor et al. (2004) inputs\nD 8–20 µm, α 15–20 km⁻¹; FOV 0.66, div. 0.75 mrad; 905 nm",
                 fontsize=11)

    ax = axes[1]
    rho_t, rho_l = OPTICS["CL61"]
    zb = np.arange(0.5, 4.01, 0.125) * 1e3
    for scale, colour, hatch, label in pvc_runs:
        lo, hi = pvc_band(zb, rho_t * scale, rho_l * scale, WAVELENGTH["CL61"], (8.0, 20.0), (15, 20))
        ax.fill_betweenx(zb / 1e3, lo, hi, color=colour, alpha=0.22, hatch=hatch, edgecolor=colour, lw=0)
        ax.plot(lo, zb / 1e3, color=colour, lw=1.5)
        ax.plot(hi, zb / 1e3, color=colour, lw=1.5, label=label)
    for curve, marker, label in ((lower, "s", "Le et al. 2026 Fig. S8, digitised: D 20 µm"),
                                 (upper, "D", "Le et al. 2026 Fig. S8, digitised: D 8 µm")):
        ax.plot(curve / S_LIQUID, DIGITISED_Z_M / 1e3, color=ink, lw=1.2, ls="--", marker=marker, ms=6, label=label)
    ax.set_title("(b) CL61: Le et al. (2026) droplet sizes, O'Connor's α\n"
                 "D 8–20 µm, α 15–20 km⁻¹; FOV 0.56, div. 0.28 mrad; 910.55 nm", fontsize=11)

    ax = axes[2]
    ax.fill_betweenx(DIGITISED_Z_M / 1e3, lower / S_LIQUID, upper / S_LIQUID, color=muted, alpha=0.15, lw=0,
                     label="Le et al. 2026 envelope (digitised, D 8–20 µm)")
    for curve in (lower, upper):
        ax.plot(curve / S_LIQUID, DIGITISED_Z_M / 1e3, color=muted, lw=1.0, ls="--")
    for tab, colour, ls, marker, label in (
            (_ETA_CL61, blue, "-", "o", "operational v2 (2026-11): a_G 3.3 µm, α 30 km⁻¹"),
            (_ETA_PVC55_CL61, aqua, "-.", "^", "2026-07 table: a_G 5.5 µm, α 10 km⁻¹"),
            (_ETA_LEGACY_CL51, yellow, ":", "v", "legacy ladder (pre-2026-07; the CL61 used the CL51 rows)")):
        ax.plot(tab[:, 1], tab[:, 0], color=colour, lw=1.8, ls=ls, marker=marker, ms=6, label=label)
    ax.set_title("(c) CL61: tables of calibration/cloud/_filters.py\nagainst the Le et al. (2026) envelope", fontsize=11)

    for ax in axes:
        ax.set_xlim(0.5, 1.0)
        ax.set_ylim(0.2, 4.1)
        ax.set_xlabel("multiple-scattering factor η", color=ink)
        ax.grid(color="#d9d8d4", lw=0.6)
        ax.tick_params(colors=muted)
        for side in ("top", "right"):
            ax.spines[side].set_visible(False)
        ax.legend(fontsize=9, loc="upper left", bbox_to_anchor=(0.0, -0.12), frameon=False)
    axes[0].set_ylabel("cloud-base height (km)", color=ink)
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=150)
    plt.close(fig)
    print(f"[fig] {path}")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--digitise", type=Path, help="Fig. S8 raster (PNG) to re-derive the envelope from")
    ap.add_argument("--fig", type=Path, help="write the three-panel summary figure (PNG) here")
    args = ap.parse_args()
    lower, upper = DIGITISED_LOWER, DIGITISED_UPPER
    if args.digitise:
        lower, upper = digitise(args.digitise)
        print("digitised lower:", np.round(lower, 2), "\ndigitised upper:", np.round(upper, 2))
    table_a_oconnor()
    table_b_le_fit(lower, upper)
    table_c_operational(lower, upper)
    table_d_alpha_sensitivity()
    table_e_adiabatic()
    if args.fig:
        make_figure(args.fig, lower, upper)


if __name__ == "__main__":
    main()
