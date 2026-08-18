"""Reproduce the MATLAB dark-measurement experiment on Payerne L1 data.

Python port / validation of ``dark_measurement_cl61_chm_cl31.m``. It exercises
the reusable kernels in :mod:`calibration.sensitivity` so that the exact code
that will later feed the network-wide sensitivity product on the monitoring
dashboard is validated first on a single, well-understood site.

Data
----
Unlike the MATLAB script (which read hood-on raw files from the Payerne NAS),
this version reads the **operational L1 network files we already have on disk**
for the same dates extracted from the MATLAB code:

    primary window : 2026-05-12 09:35 .. 14:50 UTC
    Payerne streams: 0-20000-0-06610  A=CHM15k(1064nm) B=CL31(910nm) C=CL61(910nm)

Because these are operational (telescopes uncovered) profiles, the low-altitude
field contains real atmosphere, so the per-gate noise there is an *upper bound*
(atmosphere leaks in); the high-altitude floor and the whole detection-threshold
methodology reproduce the MATLAB experiment faithfully. Estimator (b) - the
temporal first difference - is used precisely because it cancels the static
atmosphere and slow variability and isolates the white detector noise.

Run:  python validation/dark_measurement_payerne.py
Outputs: doc/reports/dark_measurement/figs/*.png  and
         doc/reports/dark_measurement_payerne.md
"""

from __future__ import annotations

import os
import sys
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import matplotlib

matplotlib.use("Agg")
import matplotlib.dates as mdates  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
import netCDF4  # noqa: E402

from calibration.cloud.calibration import (  # noqa: E402
    _beta_conversion_factor,
    _is_raw_signal,
    INSTRUMENT_CAL_DEFAULT,
)
from calibration.sensitivity import detection as det  # noqa: E402
from calibration.sensitivity import noise as nz  # noqa: E402

# ----------------------------------------------------------------------------
# Configuration
# ----------------------------------------------------------------------------
L1_DIR = Path("D:/E-PROFILE_L1_2026/0-20000-0-06610/2026/05")
STATION = "0-20000-0-06610"
STATION_NAME = "Payerne"

# Window extracted from dark_measurement_cl61_chm_cl31.m (primary, hood-on in
# the original; here the same dates from the operational L1 archive).
T1 = np.datetime64("2026-05-12T09:35")
T2 = np.datetime64("2026-05-12T14:50")

UNIT_SCALE = 1e6  # m^-1 sr^-1 -> Mm^-1 sr^-1 for display (matches MATLAB)
RANGE_TOP = 15000.0  # m

# Detection-threshold parameters. Averaging time is 30 min (per the project
# decision to characterise sensitivity at SNR=3 with 30-min averaging).
SNR_MAIN = 3
TAU_MAIN = 1800  # [s] headline averaging time = 30 min
TAU_LIST = (TAU_MAIN,)
VOLC = det.VOLCANIC_SCENARIOS[0]  # LR=60 sr, MEC=0.60 m^2/g (headline)

# Altitude grid for binned noise (30 m bins as in the ambient script).
Z_EDGES = np.arange(0, 15390 + 1, 30, dtype=float)
Z_CTR = Z_EDGES[:-1] + np.diff(Z_EDGES) / 2
Z_PROBE = (500, 1000, 2000, 3000, 5000)
NOISE_BAND = (4000.0, 7000.0)  # clean high-altitude band for headline noise

OUT_DIR = Path("doc/reports/dark_measurement")
FIG_DIR = OUT_DIR / "figs"
REPORT = Path("doc/reports/dark_measurement_payerne.md")

# Per-instrument display colours (match the MATLAB script).
COL = {
    "CL31": (0.00, 0.45, 0.74),
    "CL61": (0.85, 0.33, 0.10),
    "CHM15k": (0.47, 0.67, 0.19),
}


@dataclass
class Inst:
    """Loaded, unit-harmonised data + derived noise for one instrument."""

    name: str
    ident: str
    wavelength: float
    dtime: np.ndarray = field(default=None)         # datetime64[ns], n_time
    dtime_s: np.ndarray = field(default=None)       # seconds from t0
    r: np.ndarray = field(default=None)             # range AGL [m], n_range
    beta: np.ndarray = field(default=None)          # Mm^-1 sr^-1, time x range
    dt: float = np.nan                               # nominal sampling step [s]
    # derived
    med: np.ndarray = field(default=None)
    p25: np.ndarray = field(default=None)
    p75: np.ndarray = field(default=None)
    std: np.ndarray = field(default=None)           # plain temporal std per gate
    rms: np.ndarray = field(default=None)
    sig_b: np.ndarray = field(default=None)         # estimator (b) per gate
    tmol2: np.ndarray = field(default=None)


# ----------------------------------------------------------------------------
# Loading
# ----------------------------------------------------------------------------
def load_instrument(ident: str, name: str) -> Inst:
    """Read one Payerne L1 file, convert rcs_0 to attenuated backscatter in
    Mm^-1 sr^-1 using the production unit logic, and subset to the window."""
    fp = L1_DIR / f"L1_{STATION}_{ident}20260512.nc"
    ds = netCDF4.Dataset(fp)
    try:
        units = ds.variables["rcs_0"].units
        wl = float(ds.variables["l0_wavelength"][...])
        days = np.asarray(ds.variables["time"][:], dtype=float)
        dtime = np.datetime64("1970-01-01") + (days * 86400.0 * 1e9).astype(
            "timedelta64[ns]"
        )
        r = np.asarray(ds.variables["range"][:], dtype=float)
        rcs = np.asarray(ds.variables["rcs_0"][:], dtype=float)  # time x range
    finally:
        ds.close()

    # Reuse the production conversion: physical beta scales by a unit factor,
    # raw range-corrected signals divide by the per-instrument constant.
    if _is_raw_signal(units):
        beta = rcs / INSTRUMENT_CAL_DEFAULT.get(name, 1.0)
    else:
        beta = rcs * _beta_conversion_factor(name, units)
    beta = beta * UNIT_SCALE  # -> Mm^-1 sr^-1

    sel = (dtime >= T1) & (dtime <= T2)
    dtime = dtime[sel]
    beta = beta[sel, :]

    t0 = dtime[0]
    dtime_s = (dtime - t0) / np.timedelta64(1, "s")
    dt = float(np.median(np.diff(dtime_s)))

    inst = Inst(name=name, ident=ident, wavelength=wl, dtime=dtime,
                dtime_s=dtime_s.astype(float), r=r, beta=beta, dt=dt)
    return inst


def derive(inst: Inst) -> None:
    """Compute median/IQR/std/RMS profiles, estimator (b) noise and T_mol^2."""
    b = inst.beta
    inst.med = np.nanmedian(b, axis=0)
    inst.p25 = np.nanpercentile(b, 25, axis=0)
    inst.p75 = np.nanpercentile(b, 75, axis=0)
    inst.std = nz.plain_std(b)
    with np.errstate(invalid="ignore"):
        inst.rms = np.sqrt(np.nanmean(b**2, axis=0))
    inst.sig_b, _ = nz.first_difference_sigma(b, inst.dtime_s, inst.dt)
    a0 = det.alpha_mol0_for_wavelength(inst.wavelength)
    inst.tmol2 = det.two_way_molecular_transmission(inst.r, a0)


def band_mean(values: np.ndarray, r: np.ndarray, band=NOISE_BAND) -> float:
    sel = (r >= band[0]) & (r <= band[1])
    return float(np.nanmean(values[sel])) if np.any(sel) else np.nan


def at_alt(values: np.ndarray, r: np.ndarray, z: float) -> float:
    iz = int(np.argmin(np.abs(r - z)))
    return float(values[iz])


# ----------------------------------------------------------------------------
# Figures
# ----------------------------------------------------------------------------
def _save(fig, name: str) -> str:
    fig.patch.set_facecolor("w")
    path = FIG_DIR / name
    fig.savefig(path, dpi=150, bbox_inches="tight", facecolor="w")
    plt.close(fig)
    return f"dark_measurement/figs/{name}"


def fig_timeheight(insts) -> str:
    fig, axes = plt.subplots(1, 3, figsize=(16, 5), constrained_layout=True)
    for ax, inst in zip(axes, insts):
        bb = inst.beta.copy()
        bb[bb <= 1e-2] = 1e-2  # clamp for log display (MATLAB c_min = -2)
        pcm = ax.pcolormesh(inst.dtime, inst.r, np.log10(bb).T,
                            cmap="jet", vmin=-2, vmax=2, shading="auto")
        ax.set_ylim(0, RANGE_TOP)
        ax.set_title(f"{inst.name}  ({inst.wavelength:.0f} nm)")
        ax.set_xlabel("Time (UTC)")
        ax.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
        cb = fig.colorbar(pcm, ax=ax)
        cb.set_label(r"log$_{10}\,\beta_{att}$ [Mm$^{-1}$ sr$^{-1}$]")
    axes[0].set_ylabel("Range [m]")
    fig.suptitle(f"Dark-window attenuated backscatter - {STATION_NAME} "
                 f"{np.datetime_as_string(T1, unit='m')} .. "
                 f"{np.datetime_as_string(T2, unit='m')}",
                 fontsize=13, fontweight="bold")
    return _save(fig, "01_timeheight.png")


def fig_profiles(insts) -> str:
    fig, axes = plt.subplots(1, 3, figsize=(16, 6), constrained_layout=True)
    for inst in insts:
        c = COL[inst.name]
        axes[0].fill_betweenx(inst.r, inst.p25, inst.p75, color=c, alpha=0.15)
        axes[0].plot(inst.med, inst.r, color=c, lw=1.5, label=inst.name)
        axes[1].plot(inst.std, inst.r, color=c, lw=1.5, label=inst.name)
        axes[2].plot(inst.rms, inst.r, color=c, lw=1.5, label=inst.name)
    axes[0].axvline(0, color="k", ls=":")
    axes[0].set_xlim(-0.5, 0.5)
    axes[0].set_title("Median profile (shaded: IQR 25/75 %)")
    axes[0].set_xlabel(r"$\beta_{att}$ [Mm$^{-1}$ sr$^{-1}$]")
    for ax, ttl in ((axes[1], "Temporal noise (plain std)"),
                    (axes[2], r"RMS = $\sqrt{bias^2+\sigma^2}$")):
        ax.set_xscale("log")
        ax.set_xlim(1e-4, 1e2)
        ax.set_title(ttl)
        ax.set_xlabel(r"[Mm$^{-1}$ sr$^{-1}$]")
    for ax in axes:
        ax.set_ylim(0, RANGE_TOP)
        ax.set_ylabel("Range [m]")
        ax.grid(True, which="both", alpha=0.3)
        ax.legend(loc="best")
    fig.suptitle(f"Noise profiles - {STATION_NAME} dark window",
                 fontsize=13, fontweight="bold")
    return _save(fig, "02_profiles.png")


def fig_estimatorb(insts) -> str:
    """Estimator (b) per-gate noise vs altitude - the headline sensitivity
    quantity. Binned to 30 m for readability."""
    fig, ax = plt.subplots(1, 1, figsize=(11, 6), constrained_layout=True)
    for inst in insts:
        sig_bin, _ = nz.bin_rms(inst.sig_b, np.isfinite(inst.beta).sum(0),
                                inst.r, Z_EDGES, min_n=10)
        ax.plot(sig_bin, Z_CTR, color=COL[inst.name], lw=1.8, label=inst.name)
    ax.axvspan(1e-4, 1e2, ymin=0, ymax=300 / RANGE_TOP, color="0.5", alpha=0.12)
    ax.set_xscale("log")
    ax.set_xlim(1e-4, 1e2)
    ax.set_ylim(0, RANGE_TOP)
    ax.set_xlabel(r"$\sigma(\beta_{att})$ estimator (b) [Mm$^{-1}$ sr$^{-1}$]")
    ax.set_ylabel("Altitude [m]")
    ax.set_title("Temporal first-difference noise (estimator b) - "
                 "isolates white detector noise")
    ax.grid(True, which="both", alpha=0.3)
    ax.legend(loc="best")
    return _save(fig, "03_estimator_b_noise.png")


def fig_allan(insts) -> str:
    """Overlapping Allan deviation of the range-collapsed (4-7 km) series, with
    a tau^-1/2 white-noise guide - validates the tau scaling used downstream."""
    fig, ax = plt.subplots(1, 1, figsize=(11, 6), constrained_layout=True)
    xt = [10, 30, 60, 300, 600, 1800, 3600]
    xtl = ["10s", "30s", "1min", "5min", "10min", "30min", "1h"]
    markers = {"CL31": "o", "CL61": "s", "CHM15k": "^"}
    anchor = None
    for inst in insts:
        ts = nz.band_timeseries(inst.beta, inst.r, np.mean(NOISE_BAND),
                                (NOISE_BAND[1] - NOISE_BAND[0]) / 2)
        span = inst.dtime_s[-1] - inst.dtime_s[0]
        taus = np.unique(np.round(np.logspace(
            np.log10(inst.dt), np.log10(span / 8), 22)))
        adev, lo, hi = [], [], []
        for tau in taus:
            a, l, h = nz.overlapping_adev_ci(ts, inst.dt, tau)
            adev.append(a); lo.append(l); hi.append(h)
        adev = np.array(adev); lo = np.array(lo); hi = np.array(hi)
        yerr = np.vstack([adev - lo, hi - adev])
        ax.errorbar(taus, adev, yerr=yerr, fmt="-" + markers[inst.name],
                    color=COL[inst.name], ms=4, lw=1.3, capsize=3,
                    label=inst.name)
        if inst.name == "CL61":
            ia = np.argmax(np.isfinite(adev))
            anchor = (taus[ia:], adev[ia] * np.sqrt(taus[ia] / taus[ia:]))
    if anchor is not None:
        ax.plot(anchor[0], anchor[1], "--k", lw=1, label=r"$\tau^{-1/2}$")
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xticks(xt); ax.set_xticklabels(xtl)
    ax.set_xlabel(r"Averaging duration $\tau$")
    ax.set_ylabel(r"Overlapping ADEV of $\beta_{att}$ (4-7 km) "
                  r"[Mm$^{-1}$ sr$^{-1}$]")
    ax.set_title(r"Allan deviation - white noise follows $\tau^{-1/2}$")
    ax.grid(True, which="both", alpha=0.3)
    ax.legend(loc="best")
    return _save(fig, "04_allan.png")


def fig_rawsignal(insts) -> str:
    """Non-range-corrected signal P = beta / r^2 (median, std, RMS)."""
    fig, axes = plt.subplots(1, 3, figsize=(16, 6), constrained_layout=True)
    for inst in insts:
        c = COL[inst.name]
        r2 = inst.r**2
        r2[r2 == 0] = np.nan
        praw = inst.beta / r2
        med = np.nanmedian(praw, axis=0)
        std = nz.plain_std(praw)
        with np.errstate(invalid="ignore"):
            rms = np.sqrt(np.nanmean(praw**2, axis=0))
        axes[0].plot(med, inst.r, color=c, lw=1.5, label=inst.name)
        axes[1].plot(np.abs(std), inst.r, color=c, lw=1.5, label=inst.name)
        axes[2].plot(rms, inst.r, color=c, lw=1.5, label=inst.name)
    axes[0].axvline(0, color="k", ls=":")
    axes[0].set_title("Median raw signal")
    axes[1].set_title("Std of raw signal"); axes[1].set_xscale("log")
    axes[2].set_title("RMS of raw signal"); axes[2].set_xscale("log")
    for ax in axes:
        ax.set_ylim(0, RANGE_TOP)
        ax.set_ylabel("Range [m]")
        ax.set_xlabel(r"[Mm$^{-1}$ sr$^{-1}$ m$^{-2}$]")
        ax.grid(True, which="both", alpha=0.3)
        ax.legend(loc="best")
    fig.suptitle(f"Non-range-corrected signal - {STATION_NAME} dark window",
                 fontsize=13, fontweight="bold")
    return _save(fig, "05_rawsignal.png")


def fig_detect_backscatter(insts) -> str:
    fig, axes = plt.subplots(1, len(TAU_LIST), figsize=(14, 6),
                             constrained_layout=True)
    for ax, tau in zip(np.atleast_1d(axes), TAU_LIST):
        for inst in insts:
            a0 = det.alpha_mol0_for_wavelength(inst.wavelength)
            bmin = det.min_detectable_backscatter(
                inst.sig_b, inst.r, inst.dt, tau, SNR_MAIN, a0)
            ax.plot(bmin, inst.r, color=COL[inst.name], lw=2, label=inst.name)
        ax.set_xscale("log")
        ax.set_xlim(1e-3, 1e1)
        ax.set_ylim(0, RANGE_TOP)
        ax.set_xlabel(r"$\beta_{min}$ aerosol [Mm$^{-1}$ sr$^{-1}$]")
        ax.set_ylabel("Altitude [m]")
        ax.set_title(f"SNR={SNR_MAIN}, "
                     fr"$\tau$ = {det_human(tau)}")
        ax.grid(True, which="both", alpha=0.3)
        ax.legend(loc="best")
    fig.suptitle(f"Aerosol detection threshold - backscatter - {STATION_NAME}",
                 fontsize=13, fontweight="bold")
    return _save(fig, "06_detect_backscatter.png")


def fig_detect_ext_mass(insts) -> str:
    tau = TAU_LIST[-1]  # longest tau headline (30 min)
    fig, axes = plt.subplots(1, 2, figsize=(14, 6), constrained_layout=True)
    for inst in insts:
        a0 = det.alpha_mol0_for_wavelength(inst.wavelength)
        bmin = det.min_detectable_backscatter(
            inst.sig_b, inst.r, inst.dt, tau, SNR_MAIN, a0)
        amin = det.min_detectable_extinction(bmin, VOLC.LR)
        mmin = det.min_detectable_mass(bmin, VOLC.LR, VOLC.MEC)
        axes[0].plot(amin, inst.r, color=COL[inst.name], lw=2, label=inst.name)
        axes[1].plot(mmin, inst.r, color=COL[inst.name], lw=2, label=inst.name)
    axes[0].set_title(f"Min detectable extinction (LR={VOLC.LR:g} sr)")
    axes[0].set_xlabel(r"$\alpha_{min}$ [Mm$^{-1}$]")
    axes[0].set_xlim(1e-2, 1e3)
    axes[1].set_title(f"Min detectable mass (MEC={VOLC.MEC:g} m$^2$/g)")
    axes[1].set_xlabel(r"$M_{min}$ [$\mu$g/m$^3$]")
    axes[1].set_xlim(1e-2, 1e3)
    styles = ["--", "-.", ":"]
    for lvl, lab, s in zip(det.ICAO_LEVELS_UG_M3, det.ICAO_LABELS, styles):
        axes[1].axvline(lvl, color="k", ls=s, lw=1, label=lab)
    for ax in axes:
        ax.set_xscale("log")
        ax.set_ylim(0, RANGE_TOP)
        ax.set_ylabel("Altitude [m]")
        ax.grid(True, which="both", alpha=0.3)
        ax.legend(loc="best", fontsize=8)
    fig.suptitle(f"Aerosol detection threshold - extinction & mass "
                 f"(SNR={SNR_MAIN}, tau={det_human(tau)}) - {STATION_NAME}\n"
                 "ICAO ash zones: EUR Doc 019 / NAT Doc 006 (2010)",
                 fontsize=12, fontweight="bold")
    return _save(fig, "07_detect_extinction_mass.png")


def det_human(s: float) -> str:
    if s < 60:
        return f"{s:.0f}s"
    if s < 3600:
        return f"{s / 60:.0f}min"
    return f"{s / 3600:.1f}h"


# ----------------------------------------------------------------------------
# Report
# ----------------------------------------------------------------------------
def md_table(header, rows) -> str:
    out = ["| " + " | ".join(header) + " |",
           "| " + " | ".join("---" for _ in header) + " |"]
    for row in rows:
        out.append("| " + " | ".join(str(c) for c in row) + " |")
    return "\n".join(out)


def build_report(insts, figs) -> None:
    lines = []
    lines.append("# Dark-measurement reproduction - Payerne L1\n")
    lines.append(
        f"Python port of `dark_measurement_cl61_chm_cl31.m` run on the "
        f"operational L1 files for **{STATION_NAME}** ({STATION}), window "
        f"**{np.datetime_as_string(T1, unit='m')} .. "
        f"{np.datetime_as_string(T2, unit='m')} UTC** (the primary date from "
        f"the MATLAB code).\n")
    lines.append(
        "> **Data caveat.** The MATLAB experiment used hood-on (covered) raw "
        "files; here we reuse the operational L1 network archive for the same "
        "dates, so the telescopes are uncovered and the low-altitude field "
        "contains real atmosphere. The per-gate noise below ~1 km is therefore "
        "an *upper bound*. Estimator (b) (temporal first difference) is used "
        "because it cancels the static atmosphere and slow variability and "
        "isolates the white detector noise; the 4-7 km floor and the whole "
        "detection-threshold methodology reproduce the MATLAB experiment.\n")
    lines.append(
        "> **Absolute-scale caveat.** CL61 L1 is already physical attenuated "
        "backscatter (m^-1 sr^-1), so its noise floor and thresholds are "
        "absolute. CL31 and CHM15k are operationally **uncalibrated**: their raw "
        "rcs_0 is divided by the default constants (CL31 1e8, CHM15k 3e11), so "
        "their beta, beta_min, M_min and ICAO altitudes scale inversely with "
        "those defaults - treat the CL31/CHM15k absolute numbers as "
        "default-constant-dependent (the MATLAB used CHM_CAL=5e11, a 5/3 shift). "
        "The relative behaviour and the methodology are unaffected.\n")

    lines.append("## Instruments\n")
    rows = [(i.name, i.ident, f"{i.wavelength:.1f}", i.beta.shape[0],
             i.r.size, f"{i.dt:.0f}", f"{i.r[-1]:.0f}") for i in insts]
    lines.append(md_table(
        ["Instrument", "Stream", "lambda [nm]", "Profiles", "Gates",
         "dt [s]", "Top [m]"], rows))
    lines.append("")

    lines.append("## 1. Dark-window attenuated backscatter\n")
    lines.append(f"![time-height]({figs['timeheight']})\n")

    lines.append("## 2. Noise profiles\n")
    lines.append(f"![profiles]({figs['profiles']})\n")
    lines.append(f"![estimator-b]({figs['estimatorb']})\n")

    lines.append("### Estimator (b) noise sigma(beta_att) at probe altitudes "
                 "[Mm^-1 sr^-1]\n")
    lines.append("*In this uncovered-L1 reproduction the 500-2000 m rows remain "
                 "**upper bounds** - the first difference does not fully cancel "
                 "fast atmospheric variability between consecutive profiles. The "
                 "**4-7 km mean** is the representative detector-noise floor.*\n")
    header = ["Altitude [m]"] + [i.name for i in insts]
    rows = []
    for z in Z_PROBE:
        rows.append([z] + [f"{at_alt(i.sig_b, i.r, z):.4g}" for i in insts])
    rows.append(["mean 4-7 km"] +
                [f"{band_mean(i.sig_b, i.r):.4g}" for i in insts])
    lines.append(md_table(header, rows))
    lines.append("")

    lines.append("## 3. Allan deviation (tau scaling)\n")
    lines.append(f"![allan]({figs['allan']})\n")

    lines.append("## 4. Non-range-corrected signal\n")
    lines.append(f"![rawsignal]({figs['rawsignal']})\n")

    lines.append("## 5. Aerosol detection thresholds\n")
    lines.append(f"![detect-backscatter]({figs['detectb']})\n")
    lines.append(f"![detect-extinction-mass]({figs['detectem']})\n")

    lines.append(f"### Detection thresholds at SNR={SNR_MAIN}, {VOLC.name}\n")
    for tau in TAU_LIST:
        lines.append(f"**tau = {det_human(tau)}** "
                     "(beta in Mm^-1 sr^-1, alpha in Mm^-1, M in ug/m^3)\n")
        header = ["Instrument", "z [m]", "beta_min", "alpha_min", "M_min"]
        rows = []
        for inst in insts:
            a0 = det.alpha_mol0_for_wavelength(inst.wavelength)
            for z in Z_PROBE:
                iz = int(np.argmin(np.abs(inst.r - z)))
                bmin = det.min_detectable_backscatter(
                    inst.sig_b, inst.r, inst.dt, tau, SNR_MAIN, a0)[iz]
                amin = det.min_detectable_extinction(bmin, VOLC.LR)
                mmin = det.min_detectable_mass(bmin, VOLC.LR, VOLC.MEC)
                rows.append([inst.name, z, f"{bmin:.3e}", f"{amin:.3e}",
                             f"{mmin:.1f}"])
        lines.append(md_table(header, rows))
        lines.append("")

    lines.append(f"### ICAO-threshold detection altitude (night, "
                 f"tau={det_human(TAU_MAIN)}, SNR={SNR_MAIN}, {VOLC.name})\n")
    header = ["Instrument"] + [f"{lvl:.0f} ug/m3" for lvl in
                               det.ICAO_LEVELS_UG_M3]
    rows = []
    for inst in insts:
        a0 = det.alpha_mol0_for_wavelength(inst.wavelength)
        bmin = det.min_detectable_backscatter(
            inst.sig_b, inst.r, inst.dt, TAU_MAIN, SNR_MAIN, a0)
        mmin = det.min_detectable_mass(bmin, VOLC.LR, VOLC.MEC)
        row = [inst.name]
        for lvl in det.ICAO_LEVELS_UG_M3:
            alt = det.detection_altitude(mmin, inst.r, lvl)
            row.append("n/a" if np.isnan(alt) else f"{alt:.0f} m")
        rows.append(row)
    lines.append(md_table(header, rows))
    lines.append("")
    lines.append("*This per-station ICAO detection altitude is the headline "
                 "scalar that will populate the network map on the dashboard.*\n")

    REPORT.write_text("\n".join(lines), encoding="utf-8")
    print(f"Wrote report: {REPORT}")


# ----------------------------------------------------------------------------
# Main
# ----------------------------------------------------------------------------
def main() -> None:
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    insts = []
    for ident, name in (("A", "CHM15k"), ("B", "CL31"), ("C", "CL61")):
        inst = load_instrument(ident, name)
        derive(inst)
        print(f"{name:8s} {inst.beta.shape[0]:5d} profiles  dt={inst.dt:.0f}s  "
              f"noise(b,4-7km)={band_mean(inst.sig_b, inst.r):.3g} Mm-1sr-1")
        insts.append(inst)

    # Order for plotting consistency: CL31, CL61, CHM15k (as in MATLAB legends).
    order = {"CL31": 0, "CL61": 1, "CHM15k": 2}
    insts.sort(key=lambda i: order[i.name])

    figs = {
        "timeheight": fig_timeheight(insts),
        "profiles": fig_profiles(insts),
        "estimatorb": fig_estimatorb(insts),
        "allan": fig_allan(insts),
        "rawsignal": fig_rawsignal(insts),
        "detectb": fig_detect_backscatter(insts),
        "detectem": fig_detect_ext_mass(insts),
    }
    build_report(insts, figs)
    print("Done.")


if __name__ == "__main__":
    main()
