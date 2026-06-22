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
