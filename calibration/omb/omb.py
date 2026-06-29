"""Core Observation-minus-Background computation against CAMS.

Python port of the per-station O-B logic in ``E_PROFILE_ALC_Monthly_OB.m``,
updated to the E-PROFILE aerosol-processing recommendation (2026-06):

1. cloud-screen each native profile (drop low-cloud profiles, mask gates at/above
   the cloud base) -- O-B is a clear-sky comparison,
2. pre-average the (already calibrated) attenuated backscatter to a 5 min grid
   (reproducing the L2 cadence),
3. per 5 min profile, estimate the noise floor from the top ``noise_top_m`` of the
   de-range-corrected signal and drop gates below ``snr_min`` * sigma (SNR screen),
4. for each CAMS step, average the valid 5 min profiles within +/- ``agg_halfwin_min``
   onto a coarse vertical grid (default 150 m),
5. at 910 nm, divide the observation by the two-way water-vapour transmission
   (reusing :func:`calibration.cloud.calibration.compute_wv_transmission`),
6. read the CAMS aerosol backscatter at the instrument wavelength, interpolate the
   observation onto the CAMS height grid (ASL) and form bias = observation - background,
7. for the period statistics, drop altitudes valid in < ``valid_frac_min`` of steps.

One or more observation *sources* can be passed at once (e.g. the operational
L2 constant and our Kalman best-estimate C_L), all sharing a single CAMS read,
SNR screen (scale-invariant) and water-vapour correction.
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
    obs_full: Dict[str, NDArray] = field(default_factory=dict)    # src -> (n_cams, n_r) UNSCREENED obs (pcolor)
    cloud_base: NDArray = field(default_factory=lambda: np.empty(0))  # (n_cams,) cloud base AGL [m], NaN=clear


def _datenum_to_dt64(datenum: NDArray) -> NDArray:
    days = np.asarray(datenum, dtype="float64") - _DATENUM_UNIX
    return np.datetime64("1970-01-01T00:00:00", "ns") + (days * 86400e9).astype(
        "timedelta64[ns]"
    )


def _interp_nan(x: NDArray, y: NDArray, xnew: NDArray,
                max_gap: float | None = None) -> NDArray:
    """Linear interp with NaN handling; NaN outside the data range. ``x`` need
    not be sorted (sorted internally). If ``max_gap`` is given, any interpolated
    point farther than ``max_gap`` from the nearest finite source sample is set to
    NaN, so large flagged gaps are NOT bridged (keeps availability + bias honest)."""
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    xnew = np.asarray(xnew, dtype=float)
    ok = np.isfinite(x) & np.isfinite(y)
    if ok.sum() < 2:
        return np.full(np.shape(xnew), np.nan)
    xs, ys = x[ok], y[ok]
    order = np.argsort(xs)
    xs, ys = xs[order], ys[order]
    out = np.interp(xnew, xs, ys, left=np.nan, right=np.nan)
    if max_gap is not None:
        j = np.clip(np.searchsorted(xs, xnew), 1, len(xs) - 1)
        dist = np.minimum(np.abs(xnew - xs[j]), np.abs(xnew - xs[j - 1]))
        out = np.where(dist <= max_gap, out, np.nan)
    return out


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


def _preaverage_5min(time: NDArray, beta_scr: Dict[str, NDArray],
                     preavg_min: float):
    """Average native profiles into ``preavg_min``-minute bins (reproducing the L2
    cadence). Returns ``(t5, {src: (n5, n_range)})`` with ``t5`` the bin centres."""
    t = np.asarray(time)
    t0 = t.min()
    binw = np.timedelta64(int(round(preavg_min * 60)), "s")
    idx = ((t - t0) // binw).astype("int64")
    uniq = np.unique(idx)
    t5 = t0 + uniq * binw + binw // 2
    out: Dict[str, NDArray] = {}
    with np.errstate(invalid="ignore"):
        for k, v in beta_scr.items():
            arr = np.full((uniq.size, v.shape[1]), np.nan)
            for j, u in enumerate(uniq):
                sel = idx == u
                if np.any(sel):
                    arr[j, :] = np.nanmean(v[sel, :], axis=0)
            out[k] = arr
    return t5, out


def _snr_screen(beta5: Dict[str, NDArray], range_agl: NDArray, ref_key: str,
                noise_top_m: float, snr_min: float, morph_open: bool = True) -> None:
    """In-place SNR screen on the 5 min profiles. For each profile sigma is a robust
    (1.4826 * MAD) noise floor of the de-range-corrected signal (``beta/r**2`` ~ raw
    signal, C_L-invariant) over the top ``noise_top_m`` of the finite gates; gates with
    signal < ``snr_min`` * sigma are dropped. That keep-rule is one-sided, so on its own
    it admits only the *positive* noise tail at sub-noise altitudes (biasing O-B high);
    therefore when ``morph_open`` the 2-D (5min-time x range) keep mask is morphologically
    OPENED (binary opening, cross element) to delete isolated/transient survivors -- only
    signal coherent in BOTH time and range is kept. The mask is taken from ``ref_key`` and
    applied to every source (scale-invariant -> identical gates)."""
    r = np.asarray(range_agl, dtype=float)
    r2 = np.where(r > 0, r ** 2, np.nan)
    ref = beta5[ref_key]
    keep = np.zeros(ref.shape, dtype=bool)
    for i in range(ref.shape[0]):
        s = ref[i, :] / r2                          # de-range-corrected (raw-equivalent)
        fin = np.isfinite(s)
        if fin.sum() < 5:
            keep[i, :] = fin                        # too little to estimate noise -> don't filter
            continue
        top = fin & (r >= r[fin].max() - noise_top_m)
        if top.sum() < 5:
            top = fin
        v = s[top]
        sigma = 1.4826 * np.nanmedian(np.abs(v - np.nanmedian(v)))   # robust sigma (MAD)
        if not (sigma > 0):
            keep[i, :] = fin
            continue
        keep[i, :] = fin & (s >= snr_min * sigma)   # keep only signal >= snr_min * sigma
    if morph_open and keep.any():
        from scipy.ndimage import binary_opening
        keep = binary_opening(keep)                 # drop points isolated in time/range (speckle)
    for vv in beta5.values():
        vv[~keep] = np.nan


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
    preavg_min: float = 5.0,
    agg_halfwin_min: float = 15.0,
    snr_min: float = 3.0,
    noise_top_m: float = 2000.0,
    valid_frac_min: float = 0.25,
    morph_open: bool = True,
    interp_max_gap_m: float = 300.0,
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

    # --- 5 min pre-average + per-profile SNR screen (top-2km noise floor) ----
    # Reproduce the L2 5 min cadence, then keep only gates with genuine signal
    # (>= snr_min * sigma): the noise floor is a robust std of the de-range-corrected
    # signal over the top noise_top_m, so the long-period average cannot accumulate
    # sub-noise electronic distortion. The screen is scale-invariant (identical gates
    # for every source).
    t5, beta5 = _preaverage_5min(time, beta_scr, preavg_min)
    # unscreened 5 min field of the DISPLAYED source ('ours') -- cloud + SNR screens NOT
    # applied -- kept for the "all data" observation pcolor + its flagged-data overlay.
    disp_key = "ours" if "ours" in beta_sources else ref_key
    _, beta5_full = _preaverage_5min(
        time, {disp_key: np.asarray(beta_sources[disp_key], dtype=float)}, preavg_min)
    _snr_screen(beta5, range_agl, ref_key, noise_top_m, snr_min, morph_open=morph_open)

    # --- aggregate the valid 5 min profiles within +/- agg_halfwin of each CAMS step
    agg = np.timedelta64(int(round(agg_halfwin_min * 60)), "s")
    obs_mean = {k: np.full((n_cams, n_r), np.nan) for k in beta5}
    obs_full = {disp_key: np.full((n_cams, n_r), np.nan)}  # all data (unscreened) for the pcolor
    cloud_base = np.full(n_cams, np.nan)                    # representative cloud base per step
    with np.errstate(invalid="ignore"):
        for i in range(n_cams):
            sel = (t5 >= cams_t[i] - agg) & (t5 <= cams_t[i] + agg)
            if np.any(sel):
                for k, v in beta5.items():
                    obs_mean[k][i, :] = _bin_to_grid(
                        np.nanmean(v[sel, :], axis=0), gate_bin, n_r)
                obs_full[disp_key][i, :] = _bin_to_grid(
                    np.nanmean(beta5_full[disp_key][sel, :], axis=0), gate_bin, n_r)
            if cbh is not None:                            # cloud detections in the window (incl. Vaisala vis)
                cw = cbh[(time >= cams_t[i] - agg) & (time <= cams_t[i] + agg)]
                cw = cw[np.isfinite(cw)]
                if cw.size:
                    cloud_base[i] = float(np.median(cw))

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
    # period filter: keep altitudes valid in >= valid_frac_min of the CAMS steps
    # (replaces the old absolute min-obs rule). valid_frac is returned so the figure
    # can show data availability vs altitude on a second x-axis.
    for k, om in sources.items():
        oi = np.full((z_cams.shape[0], n_cams), np.nan)
        for i in range(n_cams):
            oi[:, i] = _interp_nan(z_obs_asl, om[i, :], z_cams[:, i],
                                   max_gap=interp_max_gap_m)
        b = oi - cams_beta
        valid_frac = (np.sum(np.isfinite(b), axis=1) / n_cams
                      if n_cams else np.zeros(b.shape[0]))
        keep = valid_frac >= valid_frac_min
        b_filt = np.where(keep[:, None], b, np.nan)
        oi_filt = np.where(keep[:, None], oi, np.nan)
        with np.errstate(invalid="ignore"):
            rms_profile = np.sqrt(np.nanmean(b_filt**2, axis=1))
        obs_interp[k] = oi
        bias[k] = b
        prof[k] = {
            "obs_med": np.nanmedian(oi_filt, axis=1),
            "obs_p25": np.nanpercentile(oi_filt, 25, axis=1),
            "obs_p75": np.nanpercentile(oi_filt, 75, axis=1),
            "bias_med": np.nanmedian(b_filt, axis=1),
            "bias_p25": np.nanpercentile(b_filt, 25, axis=1),
            "bias_p75": np.nanpercentile(b_filt, 75, axis=1),
            "valid_frac": valid_frac,
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
        obs_full=obs_full,
        cloud_base=cloud_base,
    )
