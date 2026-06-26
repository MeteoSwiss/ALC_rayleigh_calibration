"""Network sensitivity over a period: per-day noise -> detection thresholds.

For each day in the window, the temporal first-difference noise (estimator b) is
estimated separately on clear-sky NIGHT and clear-sky DAY profiles, converted to
a minimum-detectable aerosol backscatter beta_min(z) at SNR and averaging time
tau, and assembled into an altitude x time field (the "evolution over time" the
user asked for). The period-mean profiles give the average minimum detectable
extinction and mass concentration (day & night), with the ICAO ash levels and
the headline detection altitude (highest z where the night threshold is below a
given ICAO concentration).

Input ``beta`` must already be on a physical, CALIBRATED scale **in Mm^-1 sr^-1**
(rcs / C_L * 1e6) to match the detection-threshold unit convention in
:mod:`calibration.sensitivity.detection` (the dark-measurement convention), so
beta_min comes out in Mm^-1 sr^-1, extinction in Mm^-1 and mass in ug/m^3.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, Optional

import numpy as np
from numpy.typing import NDArray

from . import detection as det
from . import noise as nz

__all__ = ["SensResult", "sensitivity_over_period", "combine_sens_results",
           "plot_sensitivity_station"]


@dataclass
class SensResult:
    wavelength: float
    z_ctr: NDArray                       # (n_z,) altitude bin centres AGL [m]
    dates: NDArray                       # (n_days,) datetime64[D]
    bmin_night: NDArray                  # (n_z, n_days) min detectable beta [m^-1 sr^-1]
    bmin_day: NDArray                    # (n_z, n_days)
    ext_night: NDArray = field(default=None)   # (n_z,) period-mean min ext [Mm^-1]
    ext_day: NDArray = field(default=None)
    mass_night: NDArray = field(default=None)  # (n_z,) period-mean min mass [ug/m^3]
    mass_day: NDArray = field(default=None)
    icao_alt: Dict[float, float] = field(default_factory=dict)  # level -> altitude [m]
    sigma_probe: Dict[int, float] = field(default_factory=dict)  # z[m] -> night sigma
    n_days_night: int = 0
    n_days_day: int = 0


def sensitivity_over_period(
    time: NDArray,
    beta: NDArray,
    range_agl: NDArray,
    cbh: Optional[NDArray],
    lat: float,
    lon: float,
    wavelength: float,
    *,
    dt: Optional[float] = None,
    snr: float = 3.0,
    tau: float = 1800.0,
    scenario: det.VolcanicScenario = det.VOLCANIC_SCENARIOS[0],
    z_edges: Optional[NDArray] = None,
    z_top: float = 15000.0,
    min_profiles: int = 10,
    icao_level: float = det.ICAO_LEVELS_UG_M3[0],
) -> SensResult:
    """Estimator-(b) noise -> detection thresholds per day (night & day classes)."""
    time = np.asarray(time)
    # Preserve the caller's dtype (float32) so a full-period multi-month window stays
    # within memory; the per-day first-difference subset is upcast to float64 internally.
    beta = np.asarray(beta)
    range_agl = np.asarray(range_agl, dtype=float)
    if z_edges is None:
        z_edges = np.arange(0.0, 15390.0 + 1, 30.0)
    z_ctr = z_edges[:-1] + np.diff(z_edges) / 2.0
    n_z = z_ctr.size

    dtime_s = (time - time[0]) / np.timedelta64(1, "s")
    if dt is None:
        dt = float(np.median(np.diff(dtime_s)))

    cls = nz.classify_day_night(time, lat, lon)          # 1 day, 2 night
    clear = ~np.isfinite(cbh) if cbh is not None else np.ones(time.size, bool)

    a0 = det.alpha_mol0_for_wavelength(wavelength)
    tmol2 = det.two_way_molecular_transmission(z_ctr, a0)

    day_keys = time.astype("datetime64[D]")
    uniq = np.unique(day_keys)
    bmin_night = np.full((n_z, uniq.size), np.nan)
    bmin_day = np.full((n_z, uniq.size), np.nan)

    def _bmin_for(sel: NDArray) -> NDArray:
        if np.count_nonzero(sel) < min_profiles:
            return np.full(n_z, np.nan)
        sig_gate, n_gate = nz.first_difference_sigma(beta, dtime_s, dt, sel)
        sig_bin, _ = nz.bin_rms(sig_gate, n_gate, range_agl, z_edges, min_n=10)
        return det.min_detectable_backscatter(sig_bin, z_ctr, dt, tau, snr, a0)

    for j, day in enumerate(uniq):
        on = day_keys == day
        bmin_night[:, j] = _bmin_for(on & clear & (cls == 2))
        bmin_day[:, j] = _bmin_for(on & clear & (cls == 1))

    res = SensResult(wavelength=wavelength, z_ctr=z_ctr, dates=uniq,
                     bmin_night=bmin_night, bmin_day=bmin_day)
    res.n_days_night = int(np.sum(np.any(np.isfinite(bmin_night), axis=0)))
    res.n_days_day = int(np.sum(np.any(np.isfinite(bmin_day), axis=0)))

    # Median daily detection limit (robust to bad days), not the mean.
    with np.errstate(invalid="ignore"):
        med_night = np.nanmedian(bmin_night, axis=1)
        med_day = np.nanmedian(bmin_day, axis=1)
    res.ext_night = det.min_detectable_extinction(med_night, scenario.LR)
    res.ext_day = det.min_detectable_extinction(med_day, scenario.LR)
    res.mass_night = det.min_detectable_mass(med_night, scenario.LR, scenario.MEC)
    res.mass_day = det.min_detectable_mass(med_day, scenario.LR, scenario.MEC)

    for lvl in det.ICAO_LEVELS_UG_M3:
        res.icao_alt[lvl] = det.detection_altitude(res.mass_night, z_ctr, lvl)
    for z in (500, 1000, 2000, 3000, 5000):
        iz = int(np.argmin(np.abs(z_ctr - z)))
        res.sigma_probe[z] = float(med_night[iz]) if np.isfinite(med_night[iz]) else float("nan")
    return res


def combine_sens_results(results, scenario: det.VolcanicScenario = det.VOLCANIC_SCENARIOS[0]):
    """Stitch per-chunk SensResults (same z_ctr grid) into one period-spanning result.

    Lets a long window be processed month-by-month within memory: each chunk's
    daily beta_min columns are concatenated, then the period-mean extinction/mass
    profiles, ICAO detection altitudes and probe values are recomputed from the
    combined field. Returns None if no chunk has data.
    """
    results = [r for r in results if r is not None and r.dates.size]
    if not results:
        return None
    z = results[0].z_ctr
    dates = np.concatenate([r.dates for r in results])
    bmin_night = np.concatenate([r.bmin_night for r in results], axis=1)
    bmin_day = np.concatenate([r.bmin_day for r in results], axis=1)
    order = np.argsort(dates)
    dates, bmin_night, bmin_day = dates[order], bmin_night[:, order], bmin_day[:, order]

    res = SensResult(wavelength=results[0].wavelength, z_ctr=z, dates=dates,
                     bmin_night=bmin_night, bmin_day=bmin_day)
    # Use the MEDIAN daily detection limit (robust to occasional bad days), not the mean.
    with np.errstate(invalid="ignore"):
        med_night = np.nanmedian(bmin_night, axis=1)
        med_day = np.nanmedian(bmin_day, axis=1)
    res.ext_night = det.min_detectable_extinction(med_night, scenario.LR)
    res.ext_day = det.min_detectable_extinction(med_day, scenario.LR)
    res.mass_night = det.min_detectable_mass(med_night, scenario.LR, scenario.MEC)
    res.mass_day = det.min_detectable_mass(med_day, scenario.LR, scenario.MEC)
    res.n_days_night = int(np.sum(np.any(np.isfinite(bmin_night), axis=0)))
    res.n_days_day = int(np.sum(np.any(np.isfinite(bmin_day), axis=0)))
    for lvl in det.ICAO_LEVELS_UG_M3:
        res.icao_alt[lvl] = det.detection_altitude(res.mass_night, z, lvl)
    for zp in (500, 1000, 2000, 3000, 5000):
        iz = int(np.argmin(np.abs(z - zp)))
        res.sigma_probe[zp] = float(med_night[iz]) if np.isfinite(med_night[iz]) else float("nan")
    return res


def plot_sensitivity_station(res: SensResult, instrument: str, save_path, *,
                             title: str = "",
                             snr: float = 3.0, tau: float = 1800.0,
                             scenario: det.VolcanicScenario = det.VOLCANIC_SCENARIOS[0],
                             z_top: float = 15000.0) -> str:
    """2x(2+1) figure: night/day beta_min(z,time) pcolors (left, 2 cols) + the
    period-mean extinction and mass profiles (right, 1 col) with ICAO lines."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.dates as mdates
    import matplotlib.pyplot as plt

    z = res.z_ctr
    t = res.dates.astype("datetime64[D]")
    col_n, col_d = (0.15, 0.20, 0.55), (0.93, 0.65, 0.10)
    fig = plt.figure(figsize=(16, 9))
    gs = fig.add_gridspec(2, 3)

    def _pcolor(ax, field, label):
        disp = np.log10(np.clip(field, 1e-4, None))
        pcm = ax.pcolormesh(t, z, disp, cmap="viridis", vmin=-3, vmax=1, shading="auto")
        ax.set_ylim(0, z_top); ax.set_ylabel("Altitude AGL [m]")
        ax.set_title(label)
        ax.xaxis.set_major_formatter(mdates.DateFormatter("%d-%b"))
        fig.colorbar(pcm, ax=ax,
                     label=r"log$_{10}\beta_{min}$ [Mm$^{-1}$sr$^{-1}$]")

    ax = fig.add_subplot(gs[0, 0:2])
    _pcolor(ax, res.bmin_night, f"Min detectable backscatter - NIGHT "
            f"(SNR={snr:.0f}, tau={tau/60:.0f} min)")
    ax2 = fig.add_subplot(gs[1, 0:2])
    _pcolor(ax2, res.bmin_day, "Min detectable backscatter - DAY")
    ax2.set_xlabel("Date")

    # right column: average extinction (top) and mass (bottom) profiles
    axe = fig.add_subplot(gs[0, 2])
    axe.plot(res.ext_night, z, "-", color=col_n, label="night")
    axe.plot(res.ext_day, z, "-", color=col_d, label="day")
    axe.set_xscale("log"); axe.set_xlim(1e-2, 1e3); axe.set_ylim(0, z_top)
    axe.set_xlabel(r"$\alpha_{min}$ [Mm$^{-1}$]"); axe.set_ylabel("Altitude AGL [m]")
    axe.set_title(f"Median min extinction (LR={scenario.LR:g} sr)")
    axe.grid(alpha=0.3, which="both"); axe.legend(fontsize=8)

    axm = fig.add_subplot(gs[1, 2])
    axm.plot(res.mass_night, z, "-", color=col_n, label="night")
    axm.plot(res.mass_day, z, "-", color=col_d, label="day")
    styles = ["--", "-.", ":"]
    for lvl, st in zip(det.ICAO_LEVELS_UG_M3, styles):
        axm.axvline(lvl, color="k", ls=st, lw=1, label=f"{lvl:.0f} ug/m3")
    axm.set_xscale("log"); axm.set_xlim(1e-2, 1e3); axm.set_ylim(0, z_top)
    axm.set_xlabel(r"$M_{min}$ [$\mu$g/m$^3$]"); axm.set_ylabel("Altitude AGL [m]")
    axm.set_title(f"Median min mass (MEC={scenario.MEC:g} m$^2$/g)")
    axm.grid(alpha=0.3, which="both"); axm.legend(fontsize=7)

    fig.suptitle(title or f"Sensitivity - {instrument}", fontsize=14, fontweight="bold")
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    Path(save_path).parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(save_path, dpi=150, facecolor="w")
    plt.close(fig)
    return str(save_path)
