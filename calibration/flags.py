"""Homogenized calibration-quality flags shared by the Rayleigh and cloud methods.

Single source of truth for flag values + human-readable meanings, imported by
``calibration.config.CalibrationResult`` (Rayleigh), the cloud runner, and the
monitoring dashboard. Keeping one table means a CL61 file that carries BOTH a Rayleigh
(method 0) and a cloud (method 1) calibration uses the same flag vocabulary for both.

Convention
----------
A calibration "succeeded" for a period when its flag is ``1`` (full) or ``0.5`` (partial);
``0`` means no data and negatives are failures. The labels are deliberately method-neutral;
the specific reason for a given row is carried in its free-text ``message``.

Code classification
-------------------
GENERAL (either method can emit): 1, 0.5, 0, -1, -4, -5, -6, -10, -99.
RAYLEIGH-SPECIFIC molecular-fit diagnostics: -2, -3, -7, -8, -9 (the cloud method never emits these).
"""
from __future__ import annotations

import math

#: Canonical flag -> short, method-neutral label. The exact numeric values match the
#: long-standing Rayleigh flags so existing NetCDF/CSV stay valid; -99 is added for an
#: explicit "exception/driver error" bucket (previously such cases collapsed to 0).
FLAG_MEANINGS = {
    1: "Successful",
    0.5: "Partial success",
    0: "No data",
    -1: "Unsuitable conditions",                    # Rayleigh: not a clear night; Cloud: no valid liquid-cloud profiles
    -2: "Signal not proportional to molecular",     # Rayleigh-specific
    -3: "Method disagreement",                       # Rayleigh-specific
    -4: "Missing model data (CAMS)",                 # both (WV / molecular auxiliary data)
    -5: "Signal all-NaN",                            # both
    -6: "Uncertainty exceeds value",                 # both
    -7: "Negative fit slope",                        # Rayleigh-specific
    -8: "Fit issue: |b| > a",                        # Rayleigh-specific
    -9: "Another layer with lower signal",           # Rayleigh-specific
    -10: "Closest CAMS data too far",                # both (910 nm WV: station outside CAMS domain)
    # Cloud-specific rejection reasons: when no profile survives the cloud filters, the dominant
    # rejection (the filter that removed the most profiles) is reported instead of the generic -1.
    -20: "Cloud: window transmission too low",       # Cloud-specific
    -21: "Cloud: laser energy too low",              # Cloud-specific
    -22: "Cloud: peak not sharp above (300 m)",      # Cloud-specific
    -23: "Cloud: peak not sharp below (300 m)",      # Cloud-specific
    -24: "Cloud: aerosol below cloud > limit",       # Cloud-specific
    -25: "Cloud: cloud base out of range",           # Cloud-specific
    -26: "Cloud: inconsistent neighbours",           # Cloud-specific
    -99: "Exception during calibration",             # both
}

#: Cloud filter-stat key -> rejection flag (the most-rejected filter is reported for a failed cloud
#: calibration). window/quality both map to -20 (instrument-window family).
CLOUD_REJECT_FLAG = {
    "window_rejected": -20,
    "quality_flag_rejected": -20,
    "energy_rejected": -21,
    "above_rejected": -22,
    "below_rejected": -23,
    "ratio_rejected": -24,
    "cbh_rejected": -25,
    "n_rejected": -26,
}


def dominant_cloud_reject_flag(*stat_dicts):
    """From the cloud filter-stat dicts (filter_stats, cloud_stats, consistency_stats), return
    (flag, reason, counts) for a failed cloud calibration. If nothing was rejected (genuinely no
    liquid cloud / clear sky) returns (-1.0, 'no liquid cloud', counts)."""
    counts: dict = {}
    for d in stat_dicts:
        if d:
            for k, v in d.items():
                try:
                    counts[k] = counts.get(k, 0) + int(v)
                except (TypeError, ValueError):
                    continue
    if sum(counts.values()) <= 0:
        return -1.0, "no liquid cloud", counts
    reason = max(counts, key=counts.get)
    return float(CLOUD_REJECT_FLAG.get(reason, -1)), reason, counts

#: Flags that count as a usable calibration.
SUCCESS_FLAGS = (1.0, 0.5)


def is_success(flag) -> bool:
    """True if ``flag`` is a usable calibration (1 or 0.5)."""
    try:
        return float(flag) in SUCCESS_FLAGS
    except (TypeError, ValueError):
        return False


def flag_label(flag) -> str:
    """Method-neutral label for a flag value, tolerant of float/NaN/unknown."""
    try:
        f = float(flag)
    except (TypeError, ValueError):
        return "Unknown"
    if f != f:  # NaN
        return "No calibration"
    key = f if f == 0.5 else int(f)  # 0.5 is the only non-integer flag
    return FLAG_MEANINGS.get(key, f"Unknown ({flag})")


def cloud_flag(n_profiles, cal_median, cal_std, *, min_profiles_success: int = 3) -> float:
    """Homogenized flag for one cloud (O'Connor/Hopkin) calibration from its summary stats.

    The cloud method has no native flag, so success is derived from the result: enough
    valid in-cloud profiles and a finite, positive, not-too-noisy calibration coefficient.
    CAMS-missing / exception cases are handled by the caller (they raise) and map to -4 / -99.
    """
    if (n_profiles is None or n_profiles <= 0 or cal_median is None
            or not math.isfinite(cal_median) or cal_median <= 0):
        return -1.0  # no usable liquid-cloud calibration scene
    if cal_std is not None and math.isfinite(cal_std) and cal_std > cal_median:
        return -6.0  # spread exceeds the value
    if n_profiles < min_profiles_success:
        return 0.5   # a calibration, but from few profiles
    return 1.0
