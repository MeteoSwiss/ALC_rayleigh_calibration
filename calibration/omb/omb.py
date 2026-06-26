"""Core Observation-minus-Background computation against CAMS.

Python port of the per-station O-B logic in ``E_PROFILE_ALC_Monthly_OB.m``:

1. average the (already calibrated) attenuated backscatter to the CAMS temporal
   resolution (default 3 h) and a coarse vertical grid (default 150 m),
2. screen clouds (mask the gates around/above the cloud base, drop low-cloud
   profiles),
3. at 910 nm, divide the observation by the two-way water-vapour transmission
   (reusing :func:`calibration.cloud.calibration.compute_wv_transmission`),
4. read the CAMS aerosol backscatter at the instrument wavelength,
5. interpolate the observation onto the CAMS height grid (ASL) at each matching
   time step and form the bias = observation - background.

One or more observation *sources* can be passed at once (e.g. the operational
L2 constant and our Kalman best-estimate C_L), all sharing a single CAMS read
and water-vapour correction.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, Optional

import numpy as np
from numpy.typing import NDArray

from ..cloud.calibration import (
    CeiloData,
    CloudCalConfig,
    compute_wv_transmission,
    set_defaults,
)
from ..water_vapor_correction.water_vapor import in_water_vapor_band
from .cams_aerosol import cams_aerosol_backscatter

__all__ = ["OmBResult", "compute_omb"]

_DATENUM_UNIX = 719529.0  # datenum(1970,1,1)


@dataclass
class OmBResult:
    """OmB output. All backscatter/bias fields are in native **m^-1 sr^-1**
    (multiply by 1e6 for Mm^-1 sr^-1 display, as the report/figure layer does)."""

    wavelength: float
    range_mean: NDArray                       # (n_r,) AGL [m]
    z_cams: NDArray                           # (n_lev,) mean CAMS altitude ASL [m]
    time_cams: NDArray                        # (n_cams,) datetime64[ns]
    cams_beta: NDArray                        # (n_lev, n_cams) m^-1 sr^-1
    cams_med: NDArray                         # (n_lev,) median CAMS profile
    obs_mean: Dict[str, NDArray] = field(default_factory=dict)    # src -> (n_cams, n_r)
    obs_interp: Dict[str, NDArray] = field(default_factory=dict)  # src -> (n_lev, n_cams)
    bias: Dict[str, NDArray] = field(default_factory=dict)        # src -> (n_lev, n_cams)
    prof: Dict[str, dict] = field(default_factory=dict)           # src -> profile stats
    scalar: Dict[str, dict] = field(default_factory=dict)         # src -> {mean,median,rms,n}


def _datenum_to_dt64(datenum: NDArray) -> NDArray:
    days = np.asarray(datenum, dtype="float64") - _DATENUM_UNIX
    return np.datetime64("1970-01-01T00:00:00", "ns") + (days * 86400e9).astype(
        "timedelta64[ns]"
    )


def _interp_nan(x: NDArray, y: NDArray, xnew: NDArray) -> NDArray:
    """Linear interp with NaN handling; NaN outside the data range. ``x`` need
    not be sorted (sorted internally)."""
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    ok = np.isfinite(x) & np.isfinite(y)
    if ok.sum() < 2:
        return np.full(np.shape(xnew), np.nan)
    xs, ys = x[ok], y[ok]
    order = np.argsort(xs)
    return np.interp(xnew, xs[order], ys[order], left=np.nan, right=np.nan)


def _bin_to_grid(values: NDArray, gate_bin: NDArray, n_bins: int) -> NDArray:
    """Mean of ``values`` (per native gate) within each range bin (NaN-aware)."""
    out = np.full(n_bins, np.nan)
    ok = np.isfinite(values) & (gate_bin >= 0) & (gate_bin < n_bins)
    if not np.any(ok):
        return out
    sums = np.bincount(gate_bin[ok], weights=values[ok], minlength=n_bins)
    cnts = np.bincount(gate_bin[ok], minlength=n_bins)
    nz = cnts > 0
    out[nz] = sums[nz] / cnts[nz]
    return out


def _noise_keep_mask(prof: NDArray, r: NDArray, band: tuple, snr: float) -> NDArray:
    """Keep gates whose averaged signal exceeds ``snr`` * sigma, with sigma the
    high-altitude noise floor estimated over ``band`` [m AGL].

    Equivalent of the MATLAB ``remove_noise`` SNR gate (drop averaged gates below
    n_sigma * sigma). The 3 sigma floor from a clean high band suppresses the
    aloft detector noise that would otherwise dominate the O-B bias/RMS where the
    CAMS aerosol backscatter is ~0. Returns all-True if sigma cannot be estimated.
    """
    sel = (r >= band[0]) & (r <= band[1]) & np.isfinite(prof)
    if np.count_nonzero(sel) < 3:
        return np.ones(prof.shape, dtype=bool)
    sigma = np.nanstd(prof[sel])
    if not (sigma > 0):
        return np.ones(prof.shape, dtype=bool)
    return prof >= snr * sigma


def compute_omb(
    time: NDArray,
    range_agl: NDArray,
    beta_sources: Dict[str, NDArray],
    station_lat: float,
    station_lon: float,
    station_alt: float,
    wavelength: float,
    cams_file: str,
    instrument: str,
    *,
    cloud_base_height: Optional[NDArray] = None,
    hourly_resolution: float = 3.0,
    resol_z: float = 150.0,
    y_max: float = 15000.0,
    remove_lowclouds: bool = True,
    lowcloud_height: float = 1800.0,
    cbh_guard: float = 500.0,
    do_remove_noise: bool = False,
    noise_band: tuple = (6000.0, 8000.0),
    snr_noise: float = 3.0,
    apply_wv: Optional[bool] = None,
    abs_cs_lookup_table: str = "",
) -> OmBResult:
    """Compute OmB for one station over the period covered by ``cams_file``.

    Parameters
    ----------
    time : (n_time,) datetime64[ns]
    range_agl : (n_range,) range above ground [m]
    beta_sources : dict name -> (n_time, n_range) attenuated backscatter [m^-1 sr^-1]
        e.g. {"op": ..., "ours": ...}. 910 nm sources additionally get a
        water-vapour-corrected variant "<name>_wv".
    cloud_base_height : (n_time,) lowest cloud base AGL [m] (NaN = clear), or None.
    """
    time = np.asarray(time)
    range_agl = np.asarray(range_agl, dtype=float)
    if apply_wv is None:
        apply_wv = in_water_vapor_band(wavelength)

    # --- CAMS aerosol backscatter at the instrument wavelength --------------
    cams_dn_all, z_cams_all, cams_beta_all = cams_aerosol_backscatter(
        cams_file, station_lat, station_lon, wavelength)
    cams_t_all = _datenum_to_dt64(cams_dn_all)

    # Keep only CAMS steps that fall within the observation span (+ the window).
    half = np.timedelta64(int(hourly_resolution * 1800), "s")  # hourly_resolution/2 h
    in_span = (cams_t_all >= time.min() - half) & (cams_t_all <= time.max() + half)
    cams_t = cams_t_all[in_span]
    cams_dn = cams_dn_all[in_span]
    z_cams = z_cams_all[:, in_span]
    cams_beta = cams_beta_all[:, in_span]
    n_cams = cams_t.size
    if n_cams == 0:
        raise ValueError("No CAMS steps overlap the observation period.")

    # --- coarse vertical grid (lower-edge labels, anchored at first gate) ---
    # Matches average_ceilo_v5: range_mean = min(range):resol_z:max(range), each
    # bin averaging [range_mean(i), range_mean(i)+resol_z), and the lower edge is
    # used as the interpolation altitude.
    r0 = float(np.nanmin(range_agl))
    edges = np.arange(r0, y_max + resol_z, resol_z)
    range_mean = edges[:-1]
    n_r = range_mean.size
    gate_bin = np.digitize(range_agl, edges) - 1  # native gate -> bin index

    cbh = None if cloud_base_height is None else np.asarray(cloud_base_height, dtype=float)
    ref_key = next(iter(beta_sources))
    win = np.timedelta64(int(round(hourly_resolution * 3600)), "s")

    # --- cloud screening PER NATIVE PROFILE --------------------------------
    # The MATLAB script screens the 3 h-averaged bin on the bin-minimum CBH
    # (one low cloud nukes the whole bin). That was tuned for lower-density L3
    # data; at native L1 density (~hundreds of profiles per 3 h bin) almost every
    # bin contains at least one low cloud, so the whole-bin rule discards nearly
    # all data. We therefore screen each native profile on its own CBH and
    # average the survivors, which keeps the clear sub-windows of partly-cloudy
    # bins. The same screen applies to all sources (geometric -> identical gates).
    beta_scr = {k: np.array(v, dtype=float) for k, v in beta_sources.items()}
    if cbh is not None:
        for i in range(time.size):
            if not np.isfinite(cbh[i]):
                continue
            if remove_lowclouds and cbh[i] < lowcloud_height:
                for v in beta_scr.values():
                    v[i, :] = np.nan
            else:
                mask = range_agl >= (cbh[i] - cbh_guard)
                for v in beta_scr.values():
                    v[i, mask] = np.nan

    # --- average over the forward 3 h CAMS bin -----------------------------
    # remove_noise (drop averaged gates below snr*sigma) is available but OFF by
    # default: as a signal threshold it KEEPS only high-obs gates, which biases
    # the O-B high (it conditions the comparison on the observation detecting
    # aerosol). The MATLAB applies it on L3 data, but for an unbiased O-B we rely
    # on the median + the min-obs altitude filter to handle aloft noise instead.
    # When enabled, the SAME mask (reference source, scale-invariant) is applied
    # to every source so op-vs-ours stay on identical gates.
    obs_mean = {k: np.full((n_cams, n_r), np.nan) for k in beta_scr}
    for i in range(n_cams):
        sel = (time >= cams_t[i]) & (time < cams_t[i] + win)
        if not np.any(sel):
            continue
        prof_now = {
            k: _bin_to_grid(np.nanmean(v[sel, :], axis=0), gate_bin, n_r)
            for k, v in beta_scr.items()
        }
        drop = np.zeros(n_r, dtype=bool)
        if do_remove_noise:
            drop |= ~_noise_keep_mask(prof_now[ref_key], range_mean, noise_band, snr_noise)
        for k, p in prof_now.items():
            p[drop] = np.nan
            obs_mean[k][i, :] = p

    # --- water-vapour correction (910 nm): divide observation by trans2 ----
    sources = dict(obs_mean)
    if apply_wv:
        cfg = set_defaults(CloudCalConfig(
            instrument=instrument,
            cams_folder=str(Path(cams_file).parent),
            date_str=str(cams_t[0].astype("datetime64[D]")).replace("-", "")[:6],
            abs_cs_lookup_table=abs_cs_lookup_table,
        ))
        wv_data = CeiloData(
            time=cams_t,
            time_num=cams_dn,
            station_altitude=station_alt,
            station_latitude=station_lat,
            station_longitude=station_lon,
            range=range_mean,
            range_resol=resol_z,
            beta=np.zeros((n_r, n_cams)),
            cbh=np.full(n_cams, np.nan),
            quality_flag=None,
            window_transmission=None,
            laser_energy=None,
        )
        trans2 = compute_wv_transmission(wv_data, cfg)  # (n_r, n_cams)
        for k in list(obs_mean):
            sources[f"{k}_wv"] = obs_mean[k] / trans2.T

    # --- interpolate obs onto CAMS height (ASL) and form the bias ----------
    z_obs_asl = range_mean + station_alt
    obs_interp: Dict[str, NDArray] = {}
    bias: Dict[str, NDArray] = {}
    prof: Dict[str, dict] = {}
    scalar: Dict[str, dict] = {}
    # MATLAB rule: drop altitudes with < 3 days of data (nb_obs < 24/hourly_res*3),
    # capped at half the available steps so short windows still yield statistics.
    min_obs = max(3, min(int(round(3 * 24 / hourly_resolution)), n_cams // 2))
    for k, om in sources.items():
        oi = np.full((z_cams.shape[0], n_cams), np.nan)
        for i in range(n_cams):
            oi[:, i] = _interp_nan(z_obs_asl, om[i, :], z_cams[:, i])
        b = oi - cams_beta
        # keep only altitudes with enough observations across time
        keep = np.sum(np.isfinite(b), axis=1) >= min_obs
        b_filt = np.where(keep[:, None], b, np.nan)
        with np.errstate(invalid="ignore"):
            rms_profile = np.sqrt(np.nanmean(b_filt**2, axis=1))
        obs_interp[k] = oi
        bias[k] = b
        prof[k] = {
            "obs_med": np.nanmedian(oi, axis=1),
            "obs_p25": np.nanpercentile(oi, 25, axis=1),
            "obs_p75": np.nanpercentile(oi, 75, axis=1),
            "bias_med": np.nanmedian(b_filt, axis=1),
            "bias_p25": np.nanpercentile(b_filt, 25, axis=1),
            "bias_p75": np.nanpercentile(b_filt, 75, axis=1),
        }
        scalar[k] = {
            "mean_bias": float(np.nanmean(b_filt)),
            "median_bias": float(np.nanmedian(b_filt)),
            "rms": float(np.nanmean(rms_profile)),
            "n_obs": int(np.sum(np.isfinite(b_filt))),
        }

    return OmBResult(
        wavelength=wavelength,
        range_mean=range_mean,
        z_cams=np.nanmean(z_cams, axis=1),
        time_cams=cams_t,
        cams_beta=cams_beta,
        cams_med=np.nanmedian(cams_beta, axis=1),
        obs_mean=obs_mean,
        obs_interp=obs_interp,
        bias=bias,
        prof=prof,
        scalar=scalar,
    )
