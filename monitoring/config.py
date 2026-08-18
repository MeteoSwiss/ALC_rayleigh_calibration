"""Shared constants and default paths for the monitoring dashboard.

Defaults point at the live full-network Rayleigh outputs; every path can be overridden
on the build_dashboard.py command line when running against a different archive.
"""
from __future__ import annotations

import os
from pathlib import Path

# --- Default I/O locations (env-overridable; also accept a --flag on build_dashboard.py) ------------
# Set the ALC_* vars in ops/config.sh on the server; the fallbacks are the local Windows dev paths and
# now match what the runner writes (fullcal_l1_2026 / dashboard_l1_2026) and the in-repo census, so a
# bare `python scripts/build_dashboard.py` resolves to the right places.
_REPO = Path(__file__).resolve().parent.parent
_PROJECT = Path(os.environ.get("ALC_PROJECT_DIR", "C:/DATA/Projects/202606_E-PROFILE_calibration"))
DEFAULT_FULLCAL_DIR = Path(os.environ.get("ALC_FULLCAL_DIR", str(_PROJECT / "fullcal_l1_2026")))   # <key>/<key>_cal.csv
DEFAULT_MANIFEST = Path(os.environ.get("ALC_MANIFEST", str(_REPO / "validation" / "scope_l1_2026_census.json")))  # lat/lon/type
DEFAULT_L2_DIR = Path(os.environ.get("ALC_L2_DIR", "A:/E-PROFILE_L2_monthly"))   # station name/country/institution
DEFAULT_OUT_DIR = Path(os.environ.get("ALC_DASHBOARD_DIR", str(_PROJECT / "dashboard_l1_2026")))   # generated static site
DB_NAME = "calib_index.sqlite"


def _norm_base_url(u: str) -> str:
    """Normalize an image base URL: strip whitespace and ensure a single trailing slash. '' stays ''."""
    u = (u or "").strip()
    if not u:
        return ""
    return u if u.endswith("/") else u + "/"


# Public base URL for the bulky per-image assets (per-night diagnostic PNGs, OmB / sensitivity figures,
# flag examples). When set (ALC_IMG_BASE_URL, or --img-base-url), the site references those images as
# ABSOLUTE URLs under this base instead of site-relative paths -- so they can live in object storage
# (e.g. an EWC S3 bucket https://object-store.os-api.cci2.ecmwf.int/<bucket>/) while the HTML is served
# from a small web server. Empty (default) -> unchanged site-relative behaviour. The bucket layout
# mirrors the site: <base>diag/<key>/<method>_<date>.png, <base>ombsens/<key>/..., <base>flagex/...
IMG_BASE_URL = _norm_base_url(os.environ.get("ALC_IMG_BASE_URL", ""))

# "Images in bucket" mode: the per-image PNG assets (per-night diagnostics, OmB/sensitivity) are already
# mirrored in the object-storage bucket pointed at by IMG_BASE_URL, so the local diagnostic PNGs may be
# deleted to reclaim disk. When ON, the dashboard enumerates the diagnostic/OmB/sens panels from the
# persistent CSVs instead of the on-disk PNGs (whose existence is no longer required), and emits the
# bucket URLs -- so a rebuild keeps the full diag history even after the local PNGs are gone. Fresh local
# PNGs that still exist are staged for publishing as usual; deleted ones are simply not staged. Default
# OFF -> unchanged behaviour (panels enumerated by scanning the local PNGs).
IMAGES_IN_BUCKET = os.environ.get("ALC_IMAGES_IN_BUCKET", "").strip() not in ("", "0", "false", "False")


def nc_url(key: str, fname: str) -> str:
    """Public URL for a per-station calibration NetCDF. Absolute under the bucket base when set
    (mirrors diag/ombsens), else site-relative from a station page (stations/<key>.html)."""
    return (f"{IMG_BASE_URL}nc/{key}/{fname}" if IMG_BASE_URL else f"../nc/{key}/{fname}")

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
    -9: "Another layer with lower signal",
    -10: "Closest CAMS data too far",
    -11: "Rayleigh window contaminated (classification)",
    -20: "Cloud: window transmission too low",
    -21: "Cloud: laser energy too low",
    -22: "Cloud: peak not sharp above",
    -23: "Cloud: peak not sharp below",
    -24: "Cloud: aerosol below cloud",
    -25: "Cloud: cloud base out of range",
    -26: "Cloud: inconsistent neighbours",
    -99: "Exception during calibration",
}

#: Cloud rejection flags (a cloud was present but a filter rejected it). These COUNT as failures in
#: the success-rate denominator (EVERY non-success counts, including 0 and -1 -- the rate is the
#: true daily yield); kept here for the
#: per-reason outcome breakdown.
CLOUD_REJECT_FLAGS = (-20, -21, -22, -23, -24, -25, -26)

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
    -9: "#e07b39",
    -10: "#5e3c99",
    -11: "#e8a838",
    -20: "#fff2cc", -21: "#ffe699", -22: "#ffd966", -23: "#f1c232",
    -24: "#e69138", -25: "#d79b00", -26: "#bf9000",
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

# Plotly marker symbol per instrument type — used on EVERY network map so the symbol encodes the
# instrument family (colour is then free to encode the data value). scattergeo honours a per-POINT
# marker.symbol array (verified), so the maps stay single-trace and filter.js needs no restructure.
TYPE_SYMBOLS = {
    "CHM15k": "circle",
    "CL31": "square",
    "CL51": "diamond",
    "CL61": "triangle-up",
    "Mini-MPL": "star",
}


def type_symbol(itype) -> str:
    """Plotly marker symbol for an instrument type (falls back to circle)."""
    return TYPE_SYMBOLS.get(str(itype), "circle")

# --- Instrument health / availability classes (Cloudnet-style daily bar) ---------------------------
# Per-day instrument health decoded from the status strings (calibration/status/decode.py) + data
# availability + housekeeping thresholds. Cloudnet's legend is Pass/Info/Warning/Error/No data; we
# drop "info" (informational Vaisala (S) / CHM15k Note bits are ignored per operator decision) -> the
# four classes below. Written per day to <key>_status.csv (quality column) by the runner.
QUALITY_ORDER = ["error", "warning", "pass", "nodata"]
QUALITY_COLORS = {"pass": "#3cb371", "warning": "#f4c430", "error": "#e03131", "nodata": "#e9ecef"}
QUALITY_LABELS = {"pass": "Pass", "warning": "Warning", "error": "Error", "nodata": "No data"}


def quality_color(q) -> str:
    """Colour for a daily quality class (falls back to the 'no data' grey)."""
    return QUALITY_COLORS.get(str(q), QUALITY_COLORS["nodata"])

# --- Daily calibration outcome classes (availability card, row 3/4) -------------------------------
# The 26-value FLAG_MEANINGS table is far too fine for a ~3 px daily cell, so the flags are folded
# into six classes chosen to answer the operator's actual question -- "did we calibrate, and if not,
# whose fault?". The boundaries follow the measured network distribution (433 streams, 112,934
# (date, method) rows): calibrated 48.4 %, no usable scene 16.5 %, atmosphere 28.6 %, instrument
# 0.94 %, retrieval 0.57 %, no data 1.5 %. "Instrument" gets a colour that stands out precisely
# because it is rare AND actionable (dirty window / weak laser).
CAL_CLASS_ORDER = ["ok", "noscene", "nonight", "atmos", "retrieval", "instrument", "nodata"]
CAL_CLASS_LABELS = {
    "ok": "Calibrated",
    "noscene": "No usable scene",
    "nonight": "No night at this latitude",
    "atmos": "Rejected — atmosphere",
    "retrieval": "Rejected — retrieval",
    "instrument": "Rejected — instrument",
    "nodata": "No data / missing input",
}
#: Legend forms. The full labels above are what the hover shows; a legend of six of them wraps to
#: three lines on a narrow window, so the chips use these instead.
CAL_CLASS_SHORT = {
    "ok": "calibrated",
    "noscene": "no scene",
    "nonight": "no night",
    "atmos": "atmosphere",
    "retrieval": "retrieval",
    "instrument": "instrument",
    "nodata": "no data",
}
CAL_CLASS_COLORS = {
    "ok": "#1a9850",          # same green as FLAG_COLORS[1]
    "noscene": "#f6e3a1",     # pale sand: expected, the sky's fault, not a defect
    #: Twilight blue. A polar-summer day is NOT a failure and not missing data: the instrument
    #: measured all day, there was simply no astronomical night to fit. At Hopen (76.5 N) this is
    #: 300 days of 590, at Oslo 127 -- painting them the same grey as "instrument off" made a
    #: latitude read as an outage.
    "nonight": "#8fa4c8",
    "atmos": "#f0932b",
    "retrieval": "#d64545",
    "instrument": "#8e44ad",  # rare + actionable -> deliberately unlike the others
    "nodata": "#dfe3e8",
}
#: The pipeline writes flag 0 for BOTH "no measurement reached the retrieval" and "there was no
#: night to work with"; only the message separates them, so cal_class() looks at it.
NO_NIGHT_MARKER = "nighttime window"

_CAL_CLASS_OF_FLAG = {
    1: "ok", 0.5: "ok",
    -1: "noscene",
    -2: "atmos", -9: "atmos", -11: "atmos",
    -22: "atmos", -23: "atmos", -24: "atmos", -25: "atmos", -26: "atmos",
    -3: "retrieval", -6: "retrieval", -7: "retrieval", -8: "retrieval", -99: "retrieval",
    -20: "instrument", -21: "instrument",
    0: "nodata", -4: "nodata", -5: "nodata", -10: "nodata",
}


def cal_class(flag, message=None) -> str:
    """Display class for a daily calibration flag. Unknown flags read as a retrieval failure rather
    than silently as 'calibrated' -- a new flag must never be able to look like a success.

    *message* separates the two things the pipeline both files under flag 0: a day where nothing
    usable reached the retrieval, and a day where the sun never set far enough to give a night. The
    second is geography, not a fault, and it dominates whole seasons at high latitude -- so it gets
    its own class rather than being folded into "no data" beside genuinely dead instruments.
    """
    try:
        f = float(flag)
    except (TypeError, ValueError):
        return "nodata"
    key = int(f) if float(f).is_integer() else f
    if key == 0 and message and NO_NIGHT_MARKER in str(message):
        return "nonight"
    return _CAL_CLASS_OF_FLAG.get(key, "retrieval")

# --- Mean cloud cover (availability card, row 2) ---------------------------------------------------
# Octas 0..8 from the L1 per-profile cloud_amount (or the cbh-presence fallback). Clear sky reads as
# a saturated blue and overcast as a dark grey, so the row reads like a sky, not like a quality
# scale. Days the instrument does not report cloud cover at all must be z=None (paper background) --
# visually distinct from 0 octa, which is a real measurement of a clear sky.
CLOUD_COVER_SCALE = [
    [0.00, "#2b6cb0"],
    [0.25, "#7fa9d4"],
    [0.50, "#b9c4cf"],
    [0.75, "#8b929a"],
    [1.00, "#5a6068"],
]

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

# --- Instrument monitoring panel (station pages) ---------------------------------------------------
# Canonical housekeeping field (column in <key>_hk.csv) -> (label, unit, axis group). Two y-axes:
# laser power/energy and window transmission share the "pct" (%) axis; temperatures use the "temp"
# (degC) axis. The panel draws only the fields a stream actually reports (others are blank -> skipped).
HK_PANEL = [
    ("laser",         "Laser power/energy",  "%",  "pct"),
    ("window",        "Window transmission", "%",  "pct"),
    ("temp_optics",   "Optics/laser temp",   "degC", "temp"),
    ("temp_internal", "Internal temp",       "degC", "temp"),
    ("temp_detector", "Detector temp",       "degC", "temp"),
]
HK_COLORS = {"laser": "#d62728", "window": "#1f77b4", "temp_optics": "#ff7f0e",
             "temp_internal": "#2ca02c", "temp_detector": "#9467bd"}


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
     "detail": "Rayleigh needs a CLEAR night so the signal can be matched to the molecular backscatter aloft; a cloudy / aerosol-laden night fails this. Cloud needs a fully-attenuating LIQUID cloud in the search window; a clear night legitimately has none. So −1 is expected and benign for the 'off' method. Success rates on this dashboard are computed over ALL days with data — −1 days count in the denominator, so the rate reads as the true daily yield.",
     "recognize": "The most common non-success flag. Cloud → 'No liquid cloud'; Rayleigh → 'Not a clear night'."},
    {"value": -2, "methods": "Rayleigh",
     "summary": "Signal not proportional to molecular backscatter.",
     "detail": "In the reference region aloft the range-corrected signal should decay like the known molecular (Rayleigh) profile. When it does not track the molecular curve — residual aerosol, thin cloud, or instrument issues — there is no clean region to anchor the constant, and the fit is rejected.",
     "recognize": "Rayleigh diagnostic: signal and molecular curves diverge inside the fit window."},
    {"value": -3, "methods": "Rayleigh",
     "summary": "Method disagreement — two independent estimates of the lidar constant don't agree.",
     "detail": "The Rayleigh constant is computed two ways and cross-checked. "
               "(1) C_L^window = the median of signal / molecular-backscatter taken directly inside "
               "the aerosol-free molecular reference window. "
               "(2) C_L^slope = the slope of the Rayleigh fit (signal vs molecular backscatter), which "
               "equals C_L attenuated by the aerosol two-way transmittance up to the window "
               "(C_L·T_a²); the transmittance T_a² is then estimated independently from a Klett "
               "aerosol-extinction inversion integrated from the ground to the window and divided out. "
               "The relative gap |C_L^slope − C_L^window| / C_L^window is compared to the configured "
               "tolerance (threshold_quality, ≈15 %); above it, the night is rejected as unreliable. "
               "The two routes agree in clean air but diverge when aerosol sits between the ground and "
               "the window — the Klett transmittance correction and the in-window value respond to it "
               "differently — so −3 is in practice the calibration's aerosol-contamination guard. "
               "Real example: STORNOWAY (0-20000-0-03018) on 2026-03-06 — the slope and window "
               "estimates disagreed by 17.9 % (> tolerance) and the night was rejected. Across "
               "Mar–May 2026 there are 183 such rejections, typically 15–18 %.",
     "recognize": "Message reads 'Method disagreement: X%'. In the Rayleigh diagnostic, the 'Lidar "
                  "constant C_L spread' panel shows the slope-method estimate sitting away from the "
                  "window/perturbation ensemble (see example below)."},
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
    {"value": -9, "methods": "Rayleigh",
     "summary": "Another candidate window has much lower signal than the one chosen.",
     "detail": "A molecular-window selection QC. Each candidate window's fit slope is a lidar-constant "
               "proxy; in clean air these are nearly equal across the 2–6 km search range (molecular "
               "two-way transmittance varies < 1 %). When the chosen window carries excess backscatter "
               "(typically an aerosol layer in or near it), its slope is inflated relative to the "
               "cleanest (lowest-signal) candidate window — meaning a cleaner molecular layer was "
               "available and the selection is suspect. The night is rejected when the chosen window's "
               "slope exceeds a robust cleanest reference (10th percentile of the clean-window "
               "slopes) by more than the threshold (default 2.0). Validated on L1 Mar–May 2026: "
               "clean nights cluster at ~1.1–1.4 (p95 = 1.4) with a clear gap to the aerosol tail "
               "above ~2.8 (e.g. STORNAWAY-type residual free-tropospheric aerosol; 0-20000-0-10838 "
               "on 2026-03-03 scores ~2.8). It catches contamination that the slope-vs-window "
               "cross-check (−3) lets through.",
     "recognize": "Message reads 'Another layer with lower signal found (signal ratio X)'. In the "
                  "Rayleigh diagnostic the range-corrected-signal panel shows a cleaner (lower-signal) "
                  "layer than the chosen molecular window."},
    {"value": -20, "methods": "Cloud",
     "summary": "Cloud rejected — window transmission too low.",
     "detail": "Cloud rejection reasons (−20…−26) replace the generic 'no liquid cloud' when a cloud "
               "WAS present but every candidate profile failed a filter; the dominant filter (the one "
               "that rejected the most profiles) is reported. −20: the ceilometer window/blower "
               "transmission was below threshold for the cloud profiles, so they are untrustworthy.",
     "recognize": "Message 'Cloud: window transmission too low'. Genuine clear sky stays −1."},
    {"value": -21, "methods": "Cloud",
     "summary": "Cloud rejected — laser pulse energy too low.",
     "detail": "The laser pulse energy was below threshold for the cloud profiles (degraded emission), "
               "so the integrated backscatter cannot be trusted for calibration.",
     "recognize": "Message 'Cloud: laser energy too low'."},
    {"value": -22, "methods": "Cloud",
     "summary": "Cloud rejected — peak not sharp enough above.",
     "detail": "The peak-sharpness test failed on the upper side: the backscatter 300 m ABOVE the "
               "peak was not a factor ≥ 20 smaller. The cloud did not fully attenuate the beam (not a "
               "suitable opaque stratocumulus), so the O'Connor integral constraint does not hold.",
     "recognize": "Message 'Cloud: peak not sharp above'."},
    {"value": -23, "methods": "Cloud",
     "summary": "Cloud rejected — peak not sharp enough below.",
     "detail": "The peak-sharpness test failed on the lower side: the backscatter 300 m BELOW the peak "
               "was not a factor ≥ 20 smaller — drizzle or a diffuse base, not a sharp liquid cloud "
               "base, so the profile is unsuitable.",
     "recognize": "Message 'Cloud: peak not sharp below'."},
    {"value": -24, "methods": "Cloud",
     "summary": "Cloud rejected — too much aerosol below the cloud.",
     "detail": "The aerosol below the cloud contributed more than the allowed fraction (~5 %) of the "
               "total integrated backscatter, so the integral constraint B = 1/(2ηS) would be biased "
               "by the sub-cloud aerosol.",
     "recognize": "Message 'Cloud: aerosol below cloud'."},
    {"value": -25, "methods": "Cloud",
     "summary": "Cloud rejected — cloud base out of range.",
     "detail": "The detected cloud base fell outside the allowed [cbh_minheight, cbh_maxheight] band "
               "(too low — overlap/near-range artefacts — or too high), so the profile is excluded.",
     "recognize": "Message 'Cloud: cloud base out of range'."},
    {"value": -26, "methods": "Cloud",
     "summary": "Cloud rejected — inconsistent neighbouring profiles.",
     "detail": "The temporal-consistency filter removed the profiles: the lidar ratio of neighbouring "
               "in-cloud profiles did not agree within tolerance, so no stable run of consistent "
               "profiles remained — the cloud field was too variable to calibrate.",
     "recognize": "Message 'Cloud: inconsistent neighbours'."},
    {"value": -10, "methods": "Both",
     "summary": "Station outside the CAMS domain",
     "detail": "The nearest CAMS grid point is too far from the station for the water-vapour / "
               "molecular reference (distinct from -4, a missing CAMS FILE). Out-of-domain "
               "stations are served by their own regional CAMS boxes once configured.",
     "recognize": "Persistent for a station until its regional CAMS box exists; never sporadic."},
    {"value": -11, "methods": "Rayleigh",
     "summary": "Rayleigh window contaminated (classification)",
     "detail": "The eprof_v2.2 pre-fit target-classification mask flagged aerosol/cloud INSIDE "
               "the chosen molecular window, so the night is rejected rather than calibrated "
               "against a contaminated reference. An atmospheric rejection, not an instrument "
               "fault.",
     "recognize": "Appears with v2.2; the classification curtain for the night shows the layer "
                  "inside the fitted window."},
    {"value": -99, "methods": "Both",
     "summary": "Exception during calibration.",
     "detail": "The calibration code raised an unexpected error for that day — a driver / IO / edge-case bug rather than a physical rejection. The Message column carries the exception type; these are worth investigating as code issues.",
     "recognize": "Message shows an exception name (e.g. 'ValueError: ...')."},
]
