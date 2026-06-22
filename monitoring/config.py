"""Shared constants and default paths for the monitoring dashboard.

Defaults point at the live full-network Rayleigh outputs; every path can be overridden
on the build_dashboard.py command line when running against a different archive.
"""
from __future__ import annotations

from pathlib import Path

# --- Default I/O locations (override via the CLI) ---------------------------
_PROJECT = Path("C:/DATA/Projects/202606_E-PROFILE_calibration")
DEFAULT_FULLCAL_DIR = _PROJECT / "fullcal_all"          # <key>/<key>_cl.csv per instrument
DEFAULT_MANIFEST = _PROJECT / "stations_l2_manifest.json"  # lat/lon/type/n_months per station
DEFAULT_L2_DIR = Path("A:/E-PROFILE_L2_monthly")        # source of station name/country/institution
DEFAULT_OUT_DIR = _PROJECT / "dashboard"                # generated static site
DB_NAME = "calib_index.sqlite"

# --- Flag meanings ----------------------------------------------------------
# MIRROR of calibration/flags.py :: FLAG_MEANINGS (the homogenized cloud/Rayleigh table).
# Kept local on purpose so the dashboard runs without importing (or installing) the heavy
# calibration package. Keep in sync with calibration/flags.py.
FLAG_MEANINGS = {
    1: "Successful",
    0.5: "Partial success",
    0: "No data",
    -1: "Unsuitable conditions",
    -2: "Signal not proportional to molecular",
    -3: "Method disagreement",
    -4: "Missing model data (CAMS)",
    -5: "Signal all-NaN",
    -6: "Uncertainty exceeds value",
    -7: "Negative fit slope",
    -8: "Fit issue: |b| > a",
    -99: "Exception during calibration",
}

# A calibration "succeeded" on a night if its flag is one of these.
SUCCESS_FLAGS = (1.0, 0.5)

# Discrete color per flag: greens = usable, warm = atmospheric/fit reject, greys = no data.
FLAG_COLORS = {
    1: "#1a9850",
    0.5: "#66bd63",
    0: "#bdbdbd",
    -1: "#fee08b",
    -2: "#fdae61",
    -3: "#f46d43",
    -4: "#762a83",
    -5: "#9e9e9e",
    -6: "#d6604d",
    -7: "#d73027",
    -8: "#a50026",
    -99: "#000000",
}

# Instrument-type palette (E-PROFILE ALC families).
TYPE_COLORS = {
    "CHM15k": "#1f77b4",
    "CL31": "#ff7f0e",
    "CL51": "#2ca02c",
    "CL61": "#d62728",
    "Mini-MPL": "#9467bd",
}
TYPE_ORDER = ["CHM15k", "CL31", "CL51", "CL61", "Mini-MPL"]

# Theoretical (reference) lidar constant per instrument type, on the C_L scale. Used to express a
# station's median C_L as a percent of the nominal value. Mirrors INSTRUMENT_CAL_DEFAULT in the
# cloud calibration core.
THEORETICAL_CL = {"CL31": 1e8, "CL51": 1e8, "CL61": 1.0, "CHM15k": 3e11, "Mini-MPL": 5e5}


def theoretical_cl(itype):
    """Nominal C_L for an instrument type, or None if unknown."""
    return THEORETICAL_CL.get(str(itype))

# Calibration methods (a station can carry one or both; CL61 carries both).
METHOD_LABELS = {"rayleigh": "Rayleigh", "cloud": "Liquid-cloud"}
METHOD_COLORS = {"rayleigh": "#1f77b4", "cloud": "#2ca02c"}
METHOD_ORDER = ["rayleigh", "cloud"]


def method_label(method) -> str:
    return METHOD_LABELS.get(str(method), str(method))

# --- Watchlist thresholds (tunable) -----------------------------------------
RECENT_WINDOW_DAYS = 60       # span of the "recent" success-rate metric
DRIFT_SIGMA = 3.0             # robust-sigma multiple for the C_L drift alert
FAILURE_STREAK_DAYS = 21      # consecutive most-recent days with no successful cal
LOW_SUCCESS_FRAC = 0.15       # recent success fraction below this is flagged

# Approximate map extent for the network overview (Europe-centred).
MAP_LAT_RANGE = (33.0, 72.0)
MAP_LON_RANGE = (-25.0, 45.0)


# Method-specific wording for flags whose generic meaning is ambiguous, so the dashboard
# distinguishes cloud "No data" / "No liquid cloud" from Rayleigh "Not a clear night".
METHOD_FLAG_OVERRIDES = {
    "cloud":    {0: "No data", -1: "No liquid cloud"},
    "rayleigh": {-1: "Not a clear night"},
}


def flag_label(flag, method=None) -> str:
    """Human-readable label for a flag, optionally method-aware (cloud vs Rayleigh wording)."""
    try:
        f = float(flag)
    except (TypeError, ValueError):
        return "Unknown"
    if f != f:  # NaN (e.g. a station with no parseable calibration nights)
        return "No calibration"
    # 0.5 is the only non-integer flag; everything else keys as an int.
    key = f if f == 0.5 else int(f)
    ov = METHOD_FLAG_OVERRIDES.get(str(method), {})
    if key in ov:
        return ov[key]
    return FLAG_MEANINGS.get(key, f"Unknown ({flag})")


def flag_color(flag) -> str:
    """Chart color for a flag value (falls back to grey)."""
    try:
        f = float(flag)
    except (TypeError, ValueError):
        return "#cccccc"
    if f != f:  # NaN
        return "#cccccc"
    key = f if f == 0.5 else int(f)
    return FLAG_COLORS.get(key, "#cccccc")


def flag_anchor(flag) -> str:
    """Stable id/anchor fragment for a flag value: 1 -> '1', -1 -> 'm1', 0.5 -> '0_5'."""
    try:
        f = float(flag)
    except (TypeError, ValueError):
        return "unknown"
    if f != f:  # NaN
        return "nan"
    s = str(int(f)) if f == int(f) else ("0_5" if f == 0.5 else str(f))
    return s.replace("-", "m").replace(".", "_")


#: Long-form flag reference for the explanation page (flags.html). `methods` is which
#: calibration(s) can emit the flag; `example` (optional) is a site-relative diagnostic image.
FLAG_DOCS = [
    {"value": 1, "methods": "Both",
     "summary": "Calibration succeeded — a lidar constant C_L was produced and passed every check.",
     "detail": "Rayleigh: the range-corrected signal followed the molecular (Rayleigh) profile in a clean reference region and the sensitivity-perturbation ensemble agreed. Cloud: at least three fully-attenuating liquid-cloud profiles gave a consistent O'Connor coefficient.",
     "recognize": "Green point on the time series; a diagnostic image exists for the date."},
    {"value": 0.5, "methods": "Both",
     "summary": "Partial success — usable, but from limited evidence.",
     "detail": "A calibration was produced from fewer profiles / a shorter clean period than the full-success threshold. Cloud: 1–2 valid liquid-cloud profiles (below the 3-profile bar). Rayleigh: only part of the night met the clear-sky criterion. Treat the value as lower-confidence.",
     "recognize": "Light-green point on the time series."},
    {"value": 0, "methods": "Both",
     "summary": "No data — nothing could be attempted.",
     "detail": "No usable signal for that period: the file was missing or empty, or it carried no finite backscatter in the working range. Distinct from −5 (a file existed but every sample was NaN).",
     "recognize": "No point; the Message column reads 'No data'."},
    {"value": -1, "methods": "Both",
     "summary": "Unsuitable conditions — the scene is the wrong type for this method.",
     "detail": "Rayleigh needs a CLEAR night so the signal can be matched to the molecular backscatter aloft; a cloudy / aerosol-laden night fails this. Cloud needs a fully-attenuating LIQUID cloud in the search window; a clear night legitimately has none. So −1 is expected and benign for the 'off' method, and the dashboard excludes it from success-rate denominators.",
     "recognize": "The most common non-success flag. Cloud → 'No liquid cloud'; Rayleigh → 'Not a clear night'."},
    {"value": -2, "methods": "Rayleigh",
     "summary": "Signal not proportional to molecular backscatter.",
     "detail": "In the reference region aloft the range-corrected signal should decay like the known molecular (Rayleigh) profile. When it does not track the molecular curve — residual aerosol, thin cloud, or instrument issues — there is no clean region to anchor the constant, and the fit is rejected.",
     "recognize": "Rayleigh diagnostic: signal and molecular curves diverge inside the fit window."},
    {"value": -3, "methods": "Rayleigh",
     "summary": "Method disagreement.",
     "detail": "The lidar constant is estimated by more than one independent route (different reference windows / sub-methods). When those estimates disagree by more than the allowed tolerance the result is treated as unreliable and rejected.",
     "recognize": "Rayleigh diagnostic: candidate constants from different windows spread far apart."},
    {"value": -4, "methods": "Both",
     "summary": "Missing model data (CAMS).",
     "detail": "The calibration needs auxiliary model fields — molecular density (temperature/pressure) for Rayleigh, and water vapour for the transmission correction in the cloud method. When the required CAMS file or timestep is absent the calibration cannot proceed. Usually fixable by downloading the missing CAMS day.",
     "recognize": "Message mentions CAMS / model data."},
    {"value": -5, "methods": "Both",
     "summary": "Signal is all-NaN.",
     "detail": "Every backscatter sample in the working range was NaN / fill value — a data-quality failure distinct from 'No data': a file existed, but carried no finite values.",
     "recognize": "Message reads 'all-NaN'."},
    {"value": -6, "methods": "Both",
     "summary": "Uncertainty exceeds the value.",
     "detail": "The estimated uncertainty on C_L came out larger than C_L itself — the result is not significantly different from noise. Rayleigh: the perturbation-ensemble spread was huge. Cloud: the per-profile coefficients scattered more than their median. The value is withheld.",
     "recognize": "A value was computed but rejected; its relative uncertainty would exceed 100 %."},
    {"value": -7, "methods": "Rayleigh",
     "summary": "Negative fit slope.",
     "detail": "The linear fit of signal vs molecular profile returned a negative slope, which is unphysical (the constant must be positive). Indicates the fit window did not contain a valid molecular region.",
     "recognize": "Rayleigh-only molecular-fit diagnostic."},
    {"value": -8, "methods": "Rayleigh",
     "summary": "Ill-conditioned molecular fit (|b| > a).",
     "detail": "A consistency check on the molecular-fit coefficients failed — the offset term dominates the slope term, flagging an ill-conditioned fit rather than a clean molecular signal.",
     "recognize": "Rayleigh-only molecular-fit diagnostic."},
    {"value": -99, "methods": "Both",
     "summary": "Exception during calibration.",
     "detail": "The calibration code raised an unexpected error for that day — a driver / IO / edge-case bug rather than a physical rejection. The Message column carries the exception type; these are worth investigating as code issues.",
     "recognize": "Message shows an exception name (e.g. 'ValueError: ...')."},
]
