"""Correction & QC filters for the liquid-cloud calibration.

Split out of ``calibration.py`` (2026-07): multiple-scattering eta correction (+ per-type eta tables),
the instrument-health filters (window / laser energy / quality flag), the apparent lidar ratio, the
cloud filters (peak window / aerosol ratio / CBH range), the temporal-consistency filter and the
aerosol transmission correction, plus the small numeric helpers they use. Self-contained (duck-typed
``data``/``config``); re-exported from ``calibration`` for back-compat.
"""
from __future__ import annotations

import warnings
from typing import Optional, Tuple

import numpy as np
from numpy.typing import NDArray

S_THEORETICAL = 18.8                     # sr, liquid-water lidar ratio (O'Connor 2004)


_ETA_CL31 = np.array([
    [0.250, 0.93354], [0.375, 0.91224], [0.625, 0.87665], [0.875, 0.84739],
    [1.125, 0.82263], [1.375, 0.80128], [1.625, 0.78259], [1.875, 0.76604],
    [2.125, 0.75126], [2.375, 0.73795]])


_ETA_CL51 = np.array([
    [0.250, 0.95276], [0.375, 0.93630], [0.625, 0.90781], [0.875, 0.88331],
    [1.125, 0.86186], [1.375, 0.84282], [1.625, 0.82573], [1.875, 0.81027],
    [2.125, 0.79618], [2.375, 0.78327]])


_ETA_CL61 = np.array([
    [0.250, 0.95368], [0.375, 0.93759], [0.625, 0.90966], [0.875, 0.88566],
    [1.125, 0.86462], [1.375, 0.84593], [1.625, 0.82916], [1.875, 0.81397],
    [2.125, 0.80012], [2.375, 0.78742]])


_ETA_CHM15K = np.array([
    [0.250, 0.98139], [0.375, 0.97512], [0.625, 0.96241], [0.875, 0.95051],
    [1.125, 0.93943], [1.375, 0.92904], [1.625, 0.91925], [1.875, 0.90999],
    [2.125, 0.90122], [2.375, 0.89287], [2.750, 0.88110], [3.250, 0.86655],
    [3.750, 0.85318]])


_ETA_MINIMPL = np.array([
    [0.250, 0.98243], [0.375, 0.97661], [0.625, 0.96440], [0.875, 0.95288],
    [1.125, 0.94218], [1.375, 0.93214], [1.625, 0.92266], [1.875, 0.91367],
    [2.125, 0.90512], [2.375, 0.89697], [2.750, 0.88546], [3.250, 0.87119],
    [3.750, 0.85803]])


_ETA_MPL = np.array([
    [0.250, 0.98831], [0.375, 0.98667], [0.625, 0.98228], [0.875, 0.97700],
    [1.125, 0.97144], [1.375, 0.96590], [1.625, 0.96052], [1.875, 0.95533],
    [2.125, 0.95032], [2.375, 0.94547], [2.750, 0.93847], [3.250, 0.92954],
    [3.750, 0.92105]])


def apply_multiple_scattering_correction(
    beta: NDArray, range_data: NDArray, config: CloudCalConfig) -> NDArray:
    """beta *= eta(range): multiple-scattering correction of the O'Connor method.

    CL31/CL51/CL61 each use their OWN PVC table (Hogan 2006) at the measured calibration-scene
    droplet radius a_G = 5.5 um (see the tables above): the CL31's wider 0.83 mrad FOV gives a
    stronger correction than the 0.56 mrad CL51/CL61, and the CL61's slightly larger divergence
    (0.28 vs 0.21 mrad) gives it a marginally weaker correction than the CL51. This replaces the
    pre-2026-07 legacy Hopkin ladder (a_G = 8 um, CL61 borrowing the CL51 table). CHM15k, Mini-MPL
    and MPL use their own PVC tables -- the MATLAB applied NO correction for them (ones()). Unknown
    instrument types get no correction, which restores the MATLAB 'otherwise' behaviour."""
    range_km = range_data / 1000.0
    inst = config.instrument.upper()
    tables = {
        "CL31": _ETA_CL31,
        "CL51": _ETA_CL51,
        "CL61": _ETA_CL61,
        "CHM15K": _ETA_CHM15K,
        "CHM8K": _ETA_CHM15K,
        "MINI-MPL": _ETA_MINIMPL,
        "MINIMPL": _ETA_MINIMPL,
        "MPL": _ETA_MPL,
    }
    eta = tables.get(inst)
    if eta is None:
        return beta.copy()
    factor_profile = _interp1_linear_extrap(eta[:, 0], eta[:, 1], range_km)
    return beta * factor_profile[:, None]


def _interp1_linear_extrap(x: NDArray, y: NDArray, xi: NDArray) -> NDArray:
    """MATLAB interp1(x, y, xi, 'linear', 'extrap') with x ascending."""
    x = np.asarray(x, dtype="float64")
    y = np.asarray(y, dtype="float64")
    xi = np.asarray(xi, dtype="float64")
    out = np.empty_like(xi)
    # interior + edges via slopes
    idx = np.clip(np.searchsorted(x, xi) - 1, 0, len(x) - 2)
    x0 = x[idx]
    x1 = x[idx + 1]
    y0 = y[idx]
    y1 = y[idx + 1]
    slope = (y1 - y0) / (x1 - x0)
    out = y0 + slope * (xi - x0)
    return out


def apply_instrument_filters(
    beta: NDArray, data: CeiloData, config: CloudCalConfig
) -> Tuple[NDArray, Dict[str, int]]:
    """Quality-flag / window / energy filtering (vectorised over profiles). Window transmission is
    a reject-only gate (T < window_correction_threshold): the reported value is an arbitrary
    manufacturer scale and the calibration coefficient absorbs the real window attenuation, so no
    beta correction is applied (the legacy beta / (T/100)^2 correction double-counted it)."""
    beta_filtered = beta.copy()
    stats = {"window_rejected": 0, "energy_rejected": 0, "quality_flag_rejected": 0}

    def _alive():
        return ~np.all(np.isnan(beta_filtered), axis=0)

    # 1. Quality flag (only present if oriented [range x time]).
    if data.quality_flag is not None and data.quality_flag.size:
        lower_gate = _find_first(data.range >= config.cal_minheight)
        upper_gate = _find_last(data.range <= config.cal_maxheight)
        if lower_gate is not None and upper_gate is not None:
            qf = np.asarray(data.quality_flag)
            band = (qf[lower_gate:upper_gate + 1, :] if qf.shape[0] == beta.shape[0]
                    else qf[:, lower_gate:upper_gate + 1].T)            # (gate, time)
            bad = np.any(band > 0, axis=0)
            stats["quality_flag_rejected"] = int(np.sum(bad & _alive()))
            beta_filtered[:, bad] = np.nan

    # 2. Window transmission: reject-only gate (no beta correction, see docstring).
    if data.window_transmission is not None and data.window_transmission.size:
        wt = np.asarray(data.window_transmission, dtype=float)
        thr = float(config.window_correction_threshold)
        with np.errstate(invalid="ignore"):
            reject = wt < thr                                          # NaN -> False (kept), as before
        stats["window_rejected"] = int(np.sum(reject & _alive()))
        beta_filtered[:, reject] = np.nan

    # 3. Laser energy.
    if data.laser_energy is not None and data.laser_energy.size:
        le = np.asarray(data.laser_energy, dtype=float)
        with np.errstate(invalid="ignore"):
            rej = le < config.energy_threshold
        stats["energy_rejected"] = int(np.sum(rej & _alive()))
        beta_filtered[:, rej] = np.nan

    return beta_filtered, stats


def calculate_lidar_ratio(
    beta: NDArray, range_data: NDArray, config: CloudCalConfig
) -> Tuple[NDArray, NDArray]:
    """S = 1/(2*trapz(beta dz)) per profile (vectorised). Same skip rules (all-NaN / >10 % NaN / <3
    valid) and the same trapezoid as before; the common all-valid columns are integrated in one numpy
    op, only the rare partial-NaN columns fall back to the per-column loop -> bit-identical result."""
    n_profiles = beta.shape[1]
    S = np.full(n_profiles, np.nan)
    integrated_beta = np.full(n_profiles, np.nan)

    lower_gate = _find_first(range_data >= config.cal_minheight)
    upper_gate = _find_last(range_data <= config.cal_maxheight)
    if lower_gate is None or upper_gate is None:
        return S, integrated_beta

    bw = beta[lower_gate:upper_gate + 1, :]                       # (gate, time)
    x = np.asarray(range_data[lower_gate:upper_gate + 1], dtype="float64")
    dx = np.diff(x)[:, np.newaxis]                               # (gate-1, 1)
    nan_mask = np.isnan(bw)
    n_total = bw.shape[0]
    n_nan = nan_mask.sum(axis=0)
    eligible = (n_total - n_nan >= 3) & (n_nan <= n_total * 0.1)
    # all-valid columns -> one vectorised trapezoid (identical formula to _trapz)
    no_nan = eligible & (n_nan == 0)
    if np.any(no_nan):
        sub = bw[:, no_nan]
        integrated_beta[no_nan] = np.sum(0.5 * (sub[1:] + sub[:-1]) * dx, axis=0)
    # the few partial-NaN columns -> exact per-column trapz over the valid points
    for i in np.where(eligible & (n_nan > 0))[0]:
        m = ~nan_mask[:, i]
        integrated_beta[i] = _trapz(bw[m, i], x[m])
    pos = integrated_beta > 0
    S[pos] = 1.0 / (2.0 * integrated_beta[pos])
    return S, integrated_beta


def apply_cloud_filters(
    S: NDArray, beta: NDArray, data: CeiloData, config: CloudCalConfig
) -> Tuple[NDArray, Dict[str, int]]:
    """Port of ``apply_cloud_filters``: peak-based +/-300 m, aerosol ratio, CBH range."""
    n_profiles = S.size
    S_filtered = S.copy()
    stats = {"above_rejected": 0, "below_rejected": 0,
             "ratio_rejected": 0, "cbh_rejected": 0, "no_cloud_rejected": 0}

    range_resol = float(data.range[1] - data.range[0])
    gate_300m = int(round(300.0 / range_resol))
    if gate_300m < 1:
        gate_300m = 1

    lower_gate = _find_first(data.range >= config.cal_minheight)
    upper_gate = _find_last(data.range <= config.cal_maxheight)

    for i in range(n_profiles):
        if np.isnan(S_filtered[i]):
            continue
        beta_profile = beta[:, i]

        # ATTRIBUTION ONLY -- selection below is untouched (faithful MATLAB port). A profile with
        # no detected cloud base still runs the whole funnel: filters 1-3 test the pseudo-peak (the
        # strongest aerosol/noise gate) and filter 4 falls back to that pseudo-peak's HEIGHT. Its
        # rejections say nothing about cloud shape, so counting them under above/below/ratio/cbh
        # made cloudless profiles dominate the counters and the day was blamed on -22/-23/-24
        # "peak/aerosol" reasons when the truth was "no cloud in most profiles". Those rejections
        # are tallied under no_cloud_rejected, which the flag attribution maps to -1.
        has_cloud = (data.cbh is not None and i < np.size(data.cbh)
                     and np.isfinite(data.cbh[i]) and data.cbh[i] > 0)

        # peak within calibration range: zero out 1..lower_gate (MATLAB beta_roi(1:lower_gate)=0)
        beta_roi = beta_profile.copy()
        if lower_gate is not None:
            beta_roi[:lower_gate + 1] = 0.0  # MATLAB 1:lower_gate -> python [0:lower_gate+1)
        # max ignoring the fact that NaNs propagate: MATLAB max() ignores NaN
        if np.all(np.isnan(beta_roi)):
            S_filtered[i] = np.nan
            continue
        max_idx = int(np.nanargmax(beta_roi))
        max_beta = beta_roi[max_idx]

        if np.isnan(max_beta) or max_beta <= 0:
            S_filtered[i] = np.nan
            continue

        # Filter 1: 300 m above
        idx_above = max_idx + gate_300m
        if idx_above <= beta_profile.size - 1:
            beta_above = beta_profile[idx_above]
            if not np.isnan(beta_above) and beta_above * config.attenuation_factor > max_beta:
                S_filtered[i] = np.nan
                stats["above_rejected" if has_cloud else "no_cloud_rejected"] += 1
                continue

        # Filter 2: 300 m below
        idx_below = max_idx - gate_300m
        if idx_below >= 0:
            beta_below = beta_profile[idx_below]
            if not np.isnan(beta_below) and beta_below * config.attenuation_factor > max_beta:
                S_filtered[i] = np.nan
                stats["below_rejected" if has_cloud else "no_cloud_rejected"] += 1
                continue

        # Filter 3: aerosol contribution ratio
        if lower_gate is not None and max_idx > lower_gate + 5:
            idx_start_aero = max(0, lower_gate)  # MATLAB max(1,lower_gate) -> 0-based
            idx_end_aero = max_idx - 5
            if idx_end_aero > idx_start_aero:
                # MATLAB nansum(beta(idx_start:idx_end))*range_resol, inclusive end
                beta_below_cloud = np.nansum(
                    beta_profile[idx_start_aero:idx_end_aero + 1]) * range_resol
                beta_total = np.nansum(
                    beta_profile[lower_gate:upper_gate + 1]) * range_resol
                if beta_total > 0:
                    ratio = beta_below_cloud / beta_total
                    if ratio > config.ratio_filter:
                        S_filtered[i] = np.nan
                        stats["ratio_rejected" if has_cloud else "no_cloud_rejected"] += 1
                        continue

        # Filter 4: CBH range
        if data.cbh is not None and not np.isnan(data.cbh[i]):
            cbh = data.cbh[i]
        else:
            cbh = data.range[max_idx]
        if cbh < config.cbh_minheight or cbh > config.cbh_maxheight:
            S_filtered[i] = np.nan
            stats["cbh_rejected" if has_cloud else "no_cloud_rejected"] += 1
            continue

    return S_filtered, stats


def apply_temporal_consistency_filter(
    S: NDArray, config: CloudCalConfig
) -> Tuple[NDArray, Dict[str, int]]:
    """Port of ``apply_temporal_consistency_filter``: N consecutive within +/-X%."""
    n_profiles = S.size
    S_consistent = np.full(n_profiles, np.nan)
    n_consec = int(config.n_consecutive)
    range_percent = config.consistency_range
    plus_limit = 1.0 + range_percent / 100.0
    minus_limit = 1.0 - range_percent / 100.0

    for i in range(n_profiles - n_consec + 1):
        group = S[i:i + n_consec]
        if np.any(np.isnan(group)):
            continue
        mean_group = np.mean(group)
        if mean_group <= 0:
            continue
        is_consistent = np.all(
            (group >= mean_group * minus_limit) & (group <= mean_group * plus_limit))
        if is_consistent:
            S_consistent[i:i + n_consec] = S[i:i + n_consec]

    stats = {"n_rejected": int(np.sum(~np.isnan(S) & np.isnan(S_consistent)))}
    return S_consistent, stats


def apply_transmission_correction(
    beta: NDArray, data: CeiloData, S: NDArray, config: CloudCalConfig
) -> Tuple[NDArray, NDArray, NDArray]:
    """Port of ``apply_transmission_correction``: aerosol two-way transmission below cloud."""
    n_profiles = S.size
    C_corrected = np.full(n_profiles, np.nan)
    C_low = np.full(n_profiles, np.nan)
    C_high = np.full(n_profiles, np.nan)
    range_resol = float(data.range[1] - data.range[0])

    if config.aerosol_lidar_ratio is None:
        raise ValueError(
            "config.aerosol_lidar_ratio is required when apply_transmission_correction "
            "is on (MATLAB set_defaults does not define it; the runner sets it to 50).")

    for i in range(n_profiles):
        if np.isnan(S[i]):
            continue
        beta_profile = beta[:, i]
        # [~, max_idx] = max(beta_profile)  (over full profile; NaN ignored)
        if np.all(np.isnan(beta_profile)):
            continue
        max_idx = int(np.nanargmax(beta_profile))

        # MATLAB max_idx < 10 (1-based) -> first 9 gates. 0-based: max_idx < 9.
        if max_idx < 9:
            C_corrected[i] = S[i] / S_THEORETICAL
            C_low[i] = S[i] / S_THEORETICAL
            C_high[i] = S[i] / S_THEORETICAL
            continue

        # Integrated below-cloud backscatter, in the INSTRUMENT (uncalibrated) beta units:
        #   nansum(beta(5:max_idx-5))*range_resol     (MATLAB 1-based 5:(max_idx-5) inclusive
        #   -> 0-based [4 : max_idx-5] inclusive).
        seg = beta_profile[4:(max_idx - 5) + 1]
        B_aerosol_raw = np.nansum(seg) * range_resol
        # Convert to a PHYSICAL aerosol optical depth before the Beer-Lambert transmission.
        # ``beta`` here is the uncalibrated attenuated backscatter (e.g. L2 stored in
        # 1E-6*1/(m*sr) units, ~O(1)); the per-profile calibration coefficient
        # C = S/S_THEORETICAL is exactly what scales it to physical 1/(m*sr), so
        #   AOD = LR * integral(C * beta) dr = LR * C * B_aerosol_raw.
        # Without the C factor the "AOD" is ~1e6x too large for L2 input and T2 underflows
        # to 0 (the MATLAB reference shares this latent bug — its parity test is skipped).
        C_base = S[i] / S_THEORETICAL
        B_aerosol = C_base * B_aerosol_raw
        if B_aerosol <= 0:
            C_corrected[i] = C_base
            C_low[i] = C_base
            C_high[i] = C_base
            continue

        T2 = np.exp(-2.0 * config.aerosol_lidar_ratio * B_aerosol)
        T2_low = np.exp(-2.0 * config.aerosol_lidar_ratio_low * B_aerosol)
        T2_high = np.exp(-2.0 * config.aerosol_lidar_ratio_high * B_aerosol)
        # Same form as the MATLAB reference (C_corrected = C_base * T2); only B_aerosol's
        # units were corrected above (the scale bug). T2 <= 1.
        C_corrected[i] = C_base * T2
        C_low[i] = C_base * T2_low
        C_high[i] = C_base * T2_high

    return C_corrected, C_low, C_high


def _trapz(y: NDArray, x: NDArray) -> float:
    """Trapezoidal integral matching MATLAB ``trapz(x, y)`` (numpy.trapezoid, no deprecation)."""
    y = np.asarray(y, dtype="float64")
    x = np.asarray(x, dtype="float64")
    return float(np.sum(0.5 * (y[1:] + y[:-1]) * np.diff(x)))


def _find_first(mask: NDArray) -> Optional[int]:
    idx = np.where(mask)[0]
    return int(idx[0]) if idx.size else None


def _find_last(mask: NDArray) -> Optional[int]:
    idx = np.where(mask)[0]
    return int(idx[-1]) if idx.size else None
