"""Decode ceilometer status strings and classify a day's instrument health.

Bit tables are transcribed from the manufacturer handbooks and the MeteoSwiss
MATLAB reference decoders, and cross-checked against the online manuals:

* **CHM15k** ``error_ext`` — 32-bit service code. Bit *b* (0 = LSB) set → a fixed
  "Error:/Warning:/Note:" message. Source: ``plot_error_chm15k.m`` +
  ``display_error_text.m`` (Y. Poltera / M. Hervo); the Lufft CHM15k user manual.
* **Vaisala CL31 / CL51** ``error_string`` — the instrument status word, either the
  native **48-bit** (12 hex chars) or the CT25K-legacy **32-bit** (8 hex chars)
  format. Manual-verified against the Vaisala CL31 User's Guide M210482EN and
  ``FVaisala_CL31_Line2.m`` / ``FVaisala_CL31_ErrorCode.m``.
* **Vaisala CL61** ``/status`` group — ~38 per-subsystem integer fields. Value→
  severity semantics follow the Vaisala CL61 traffic-light convention
  (0 = OK, 1 = warning, 2 = error); *provisional*, to confirm vs CL61 User Guide
  M212475EN. Not yet present in E-PROFILE L1 — decoder is ready for when it lands.

Severity is one of ``error`` / ``warning`` / ``info``. Per user decision, **info**
(Vaisala internal status ``(S)`` and CHM15k ``Note:``) is neither colored nor
listed — it is decoded but dropped from the day classification and summary text.

The day quality class is one of ``pass`` / ``warning`` / ``error`` / ``nodata``
(Cloudnet-style, minus Cloudnet's "info"). It is the worst of the flag class,
the data-availability class, and the housekeeping-threshold class.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, List, Optional, Sequence, Tuple

import math


# --------------------------------------------------------------------------- #
# Severity + quality class
# --------------------------------------------------------------------------- #
class Severity:
    ERROR = "error"
    WARNING = "warning"
    INFO = "info"  # decoded but ignored for color + summary (per user)


# Quality classes, ordered worst-last so ``max(..., key=QUALITY_ORDER.index)`` works.
QUALITY_ORDER = ["nodata", "pass", "warning", "error"]

# Severity -> quality class it contributes (info contributes nothing).
_SEV_TO_QUALITY = {Severity.ERROR: "error", Severity.WARNING: "warning"}

# Hourly prevalence gate: a decoded warning/alarm is credited to a clock-hour only if it appears in
# MORE than this fraction of that hour's profiles (default 5%). This filters transient single-profile
# blips (e.g. the CHM15k bit-12 laser-controller-temperature firmware quirk, present in ~0.1% of
# profiles) from painting a whole day 'warning'; a genuine sustained fault still exceeds it.
PREVALENCE_THRESHOLD = 0.05


def _worst(*classes: str) -> str:
    """Return the most severe of several quality classes ('error' beats 'warning' …)."""
    present = [c for c in classes if c]
    if not present:
        return "pass"
    return max(present, key=lambda c: QUALITY_ORDER.index(c) if c in QUALITY_ORDER else 0)


# --------------------------------------------------------------------------- #
# CHM15k — error_ext (32-bit service code); bit b (0=LSB) -> (name, severity)
# --------------------------------------------------------------------------- #
CHM15K_ERROR_EXT: Dict[int, Tuple[str, str]] = {
    0: ("Error: Signal quality", Severity.ERROR),
    1: ("Error: Signal recording", Severity.ERROR),
    2: ("Error: Signal values null or void", Severity.ERROR),
    3: ("Error: Signal recording error channel 2", Severity.ERROR),
    4: ("Error: Create new NetCDF file", Severity.ERROR),
    5: ("Error: Write / add to NetCDF", Severity.ERROR),
    6: ("Error: RS485 telegram cannot be generated/transmitted", Severity.ERROR),
    7: ("Error: Mount SD card failed", Severity.ERROR),
    8: ("Error: Detector high voltage control failed / cable defect", Severity.ERROR),
    9: ("Error: Inner housing temperature out of range", Severity.ERROR),
    10: ("Error: Laser optical unit temperature error", Severity.ERROR),
    11: ("Error: Laser trigger not detected", Severity.ERROR),
    12: ("Warning: Laser driver board temperature", Severity.WARNING),
    13: ("Error: Laser interlock", Severity.ERROR),
    14: ("Error: Laser head temperature", Severity.ERROR),
    15: ("Warning: Replace laser - ageing", Severity.WARNING),
    16: ("Warning: Signal quality - low signal/noise level", Severity.WARNING),
    17: ("Warning: Windows contaminated", Severity.WARNING),
    18: ("Warning: Signal processing", Severity.WARNING),
    19: ("Warning: Max. detection range cannot be determined", Severity.WARNING),
    20: ("Warning: File system, fsck repaired bad sectors", Severity.WARNING),
    21: ("Warning: RS485 baud rate/transfer mode reset", Severity.WARNING),
    22: ("Warning: AFD", Severity.WARNING),
    23: ("Warning: Configuration problem", Severity.WARNING),
    24: ("Warning: Laser optical unit temperature", Severity.WARNING),
    25: ("Warning: External temperature", Severity.WARNING),
    26: ("Warning: Detector temperature out of range", Severity.WARNING),
    27: ("Warning: General laser issue", Severity.WARNING),
    28: ("Note: NOL > 3 and standard telegram selected", Severity.INFO),
    29: ("Note: Power save mode on", Severity.INFO),
    30: ("Note: Standby mode on", Severity.INFO),
    # bit 31 = "not used"
}


# --------------------------------------------------------------------------- #
# Vaisala CL31/CL51 — error_string (48-bit, native); bit b (0=LSB) -> (name, sev)
# Manual bit numbering (b47 = MSB) == integer bit position, so a packed integer
# decodes by plain (value >> b) & 1. Suffix (A)=alarm/error, (W)=warning, (S)=status.
# --------------------------------------------------------------------------- #
CL31_STATUS_48: Dict[int, Tuple[str, str]] = {
    47: ("Transmitter shut-off (A)", Severity.ERROR),
    46: ("Transmitter failure (A)", Severity.ERROR),
    45: ("Receiver failure (A)", Severity.ERROR),
    44: ("Voltage failure (A)", Severity.ERROR),
    43: ("Alignment failure (A)", Severity.ERROR),
    42: ("Memory error (A)", Severity.ERROR),
    41: ("Light path obstruction (A)", Severity.ERROR),
    40: ("Receiver saturation (A)", Severity.ERROR),
    33: ("Coaxial cable failure (A)", Severity.ERROR),
    32: ("Ceilometer engine board failure (A)", Severity.ERROR),
    31: ("Window contamination (W)", Severity.WARNING),
    30: ("Battery voltage low (W)", Severity.WARNING),
    29: ("Transmitter expires (W)", Severity.WARNING),
    28: ("High humidity (W)", Severity.WARNING),
    26: ("Blower failure (W)", Severity.WARNING),
    24: ("Humidity sensor failure (W)", Severity.WARNING),
    23: ("Heater fault (W)", Severity.WARNING),
    22: ("High background radiance (W)", Severity.WARNING),
    21: ("Ceilometer engine board failure (W)", Severity.WARNING),
    20: ("Battery failure (W)", Severity.WARNING),
    19: ("Laser monitor failure (W)", Severity.WARNING),
    18: ("Receiver warning (W)", Severity.WARNING),
    17: ("Tilt angle > 45 degrees (W)", Severity.WARNING),
    15: ("Blower is on (S)", Severity.INFO),
    14: ("Blower heater is on (S)", Severity.INFO),
    13: ("Internal heater is on (S)", Severity.INFO),
    12: ("Working from battery (S)", Severity.INFO),
    11: ("Standby mode is on (S)", Severity.INFO),
    10: ("Self test in progress (S)", Severity.INFO),
    9: ("Manual data acquisition settings effective (S)", Severity.INFO),
    7: ("Units are meters (S)", Severity.INFO),
    6: ("Manual blower control (S)", Severity.INFO),
    5: ("Polling mode is on (S)", Severity.INFO),
}

# CT25K-legacy 32-bit (8 hex char) format — a DIFFERENT map (Vaisala CL31 UG, sec.
# "Data Message ... 4-byte hex"). Used when a station outputs the condensed status.
CL31_STATUS_32: Dict[int, Tuple[str, str]] = {
    31: ("Transmitter shut-off (A)", Severity.ERROR),
    30: ("Transmitter failure (A)", Severity.ERROR),
    29: ("Receiver or coaxial cable failure (A)", Severity.ERROR),
    28: ("Engine, voltage or memory failure (A)", Severity.ERROR),
    23: ("Window contaminated (W)", Severity.WARNING),
    22: ("Battery low (W)", Severity.WARNING),
    21: ("Transmitter expires (W)", Severity.WARNING),
    20: ("Heater or humidity sensor failure (W)", Severity.WARNING),
    19: ("High background radiance (W)", Severity.WARNING),
    18: ("Engine/receiver/laser monitor failure (W)", Severity.WARNING),
    17: ("High relative humidity > 85% (W)", Severity.WARNING),
    16: ("Light path obstruction or receiver saturation (A)", Severity.ERROR),
    15: ("Blower failure (W)", Severity.WARNING),
}


# --------------------------------------------------------------------------- #
# CHM15k sky-condition index (environmental; not an instrument fault)
# --------------------------------------------------------------------------- #
SCI_MEANINGS = {
    0: "nothing",
    1: "rain",
    2: "fog",
    3: "snow",
    4: "precipitation or particles on window",
}


def _as_int(value) -> Optional[int]:
    """Coerce a decimal status value (often a float holding an int) to int, else None."""
    try:
        if value is None or value == "":
            return None
        f = float(value)
    except (TypeError, ValueError):
        return None
    if not math.isfinite(f) or f < 0:  # NaN / fill (-999) -> no status
        return None
    return int(round(f))


def _hex_status_to_int_fmt(value, fmt: str) -> Tuple[Optional[int], str]:
    """Parse a Vaisala status value to (int, fmt).

    The E-PROFILE L1 ``error_string`` is the raw hex status string (e.g.
    ``"00000000C080"`` = 12 hex chars = 48-bit; ``"C0800000"`` = 8 hex = 32-bit).
    The string length selects the format. A plain number is also accepted and
    decoded with ``fmt`` (default 48-bit); ``fmt='auto'`` on a number falls to 48.
    """
    if isinstance(value, bytes):
        value = value.decode("ascii", "ignore")
    if isinstance(value, str):
        s = value.strip().strip("\x00").replace(" ", "")
        if not s:
            return None, fmt
        try:
            iv = int(s, 16)
        except ValueError:
            iv = _as_int(s)  # tolerate a decimal string
            return iv, ("48" if fmt == "auto" else fmt)
        length = len(s)
        detected = "32" if length <= 8 else "48"
        return (iv if iv >= 0 else None), detected
    # numeric input
    return _as_int(value), ("48" if fmt == "auto" else fmt)


def _decode_int(iv: Optional[int], table: Dict[int, Tuple[str, str]]) -> List[Tuple[str, str]]:
    if iv is None or iv == 0:
        return []
    return [(name, sev) for bit, (name, sev) in table.items() if (iv >> bit) & 1]


def decode_chm15k_error_ext(value) -> List[Tuple[str, str]]:
    """Decode a CHM15k ``error_ext`` service code into [(name, severity), ...]."""
    return _decode_int(_as_int(value), CHM15K_ERROR_EXT)


def decode_vaisala_error_string(value, fmt: str = "auto") -> List[Tuple[str, str]]:
    """Decode a Vaisala CL31/CL51 ``error_string`` into [(name, severity), ...].

    Accepts the raw hex string (E-PROFILE L1 form) or a number. ``fmt`` selects
    the bit map when it can't be inferred: ``"auto"`` (default; hex-string length
    decides, numbers → 48-bit), ``"48"`` (native 12-hex), or ``"32"`` (CT25K 8-hex).
    """
    iv, detected = _hex_status_to_int_fmt(value, str(fmt))
    table = CL31_STATUS_32 if detected == "32" else CL31_STATUS_48
    return _decode_int(iv, table)


def decode_sci(value) -> Optional[str]:
    """Decode a CHM15k sky-condition index to text (environmental, not a fault)."""
    iv = _as_int(value)
    if iv is None:
        return None
    return SCI_MEANINGS.get(iv)


# --------------------------------------------------------------------------- #
# CL61 /status group — subsystem int fields (provisional; confirm vs M212475EN)
# --------------------------------------------------------------------------- #
# Vaisala CL61 traffic-light convention (to confirm): 0=OK, 1=warning, 2=error.
CL61_VALUE_SEVERITY = {0: None, 1: Severity.WARNING, 2: Severity.ERROR}


def decode_cl61_status(
    fields: Dict[str, object],
    value_severity: Optional[Dict[int, Optional[str]]] = None,
) -> List[Tuple[str, str]]:
    """Decode a CL61 ``/status`` subsystem snapshot into [(name, severity), ...].

    ``fields`` maps subsystem name -> integer status value for one profile
    (e.g. ``{"Transmitter_overall": 0, "Window_condition": 1, ...}``). Value
    semantics follow ``value_severity`` (default: 0 OK / 1 warning / 2 error);
    unknown non-zero values fall back to warning. A trailing ``_overall``/
    ``_failure`` field escalates to error when non-zero even if the value is 1.
    """
    vmap = value_severity or CL61_VALUE_SEVERITY
    out: List[Tuple[str, str]] = []
    for name, raw in (fields or {}).items():
        iv = _as_int(raw)
        if iv is None or iv == 0:
            continue
        sev = vmap.get(iv)
        if sev is None:
            sev = Severity.WARNING  # unknown non-zero code -> at least a warning
        pretty = name.replace("_", " ")
        out.append((pretty, sev))
    return out


# --------------------------------------------------------------------------- #
# Housekeeping thresholds -> degraded-HK warning
# --------------------------------------------------------------------------- #
@dataclass
class HKThresholds:
    """Daily-mean housekeeping thresholds that raise a 'warning' quality class.

    Conservative defaults chosen to flag clearly-degraded hardware without firing
    on normal seasonal variation (window transmission commonly 70-100 %).
    """
    window_min_pct: float = 60.0      # window transmission below this = dirty/degraded
    laser_min_pct: float = 80.0       # laser power/energy below this = weak transmitter
    temp_internal_max_c: float = 45.0  # internal electronics running hot
    temp_internal_min_c: float = -25.0

    def evaluate(self, hk: Optional[Dict[str, float]]) -> List[str]:
        """Return a list of human-readable HK degradation notes (empty if healthy)."""
        if not hk:
            return []
        notes: List[str] = []
        w = _safe_float(hk.get("window"))
        if w is not None and w < self.window_min_pct:
            notes.append(f"Low window transmission ({w:.0f}%)")
        la = _safe_float(hk.get("laser"))
        if la is not None and la < self.laser_min_pct:
            notes.append(f"Low laser power ({la:.0f}%)")
        ti = _safe_float(hk.get("temp_internal"))
        if ti is not None and (ti > self.temp_internal_max_c or ti < self.temp_internal_min_c):
            notes.append(f"Internal temperature out of range ({ti:.0f}°C)")
        return notes


def _safe_float(v) -> Optional[float]:
    try:
        if v is None or v == "":
            return None
        f = float(v)
    except (TypeError, ValueError):
        return None
    return f if math.isfinite(f) else None


# --------------------------------------------------------------------------- #
# Per-day aggregation + classification
# --------------------------------------------------------------------------- #
@dataclass
class DayStatus:
    """One day of decoded instrument health for a stream."""
    quality: str = "nodata"            # pass | warning | error | nodata
    summary: str = "No data"           # ready-to-show reporting string
    n_profiles: int = 0
    coverage_pct: float = float("nan")  # profiles / expected * 100
    gap_hours: int = 0                  # clock-hours (0-23) with no profile
    n_alarm_h: int = 0                  # clock-hours with >=1 alarm/error flag
    n_warn_h: int = 0                   # clock-hours with >=1 warning flag
    flags: Dict[str, int] = field(default_factory=dict)  # flag name -> hours present

    def to_row(self) -> Dict[str, object]:
        """Flatten to a CSV row (flags as a compact JSON string)."""
        import json
        return {
            "quality": self.quality,
            "n_profiles": self.n_profiles,
            "coverage_pct": "" if math.isnan(self.coverage_pct) else f"{self.coverage_pct:.1f}",
            "gap_hours": self.gap_hours,
            "n_alarm_h": self.n_alarm_h,
            "n_warn_h": self.n_warn_h,
            "flags_json": json.dumps(self.flags, ensure_ascii=False),
            "summary": self.summary,
        }


def _infer_interval_s(times_epoch: Sequence[float]) -> Optional[float]:
    """Median sample interval (s) from monotonic profile times; None if <2 points."""
    ts = sorted(float(t) for t in times_epoch if t is not None and math.isfinite(float(t)))
    if len(ts) < 2:
        return None
    diffs = [b - a for a, b in zip(ts, ts[1:]) if b > a]
    if not diffs:
        return None
    diffs.sort()
    return diffs[len(diffs) // 2] or None


def _hour_of(t: float) -> int:
    return int((int(t) // 3600) % 24)


def summarize_day(
    times_epoch: Optional[Sequence[float]],
    itype: str,
    *,
    status_values: Optional[Sequence[float]] = None,
    cl61_status: Optional[Dict[str, Sequence[float]]] = None,
    hk: Optional[Dict[str, float]] = None,
    thresholds: Optional[HKThresholds] = None,
    error_string_fmt: str = "auto",
    prevalence_threshold: float = PREVALENCE_THRESHOLD,
    expected_interval_s: Optional[float] = None,
) -> DayStatus:
    """Aggregate one day of profiles into a :class:`DayStatus`.

    Parameters
    ----------
    times_epoch : profile timestamps (seconds since epoch). Empty/None -> nodata.
    itype       : instrument type string ("CHM15k" | "CL31" | "CL51" | "CL61" | ...).
    status_values : per-profile scalar status word (``error_ext`` for CHM15k,
                    ``error_string`` for CL31/CL51), aligned to ``times_epoch``.
    cl61_status : {subsystem_name: per-profile array} for CL61 (ready-for-later).
    hk          : daily-mean housekeeping ({window, laser, temp_internal, ...}).
    thresholds  : HK thresholds (defaults if None).
    error_string_fmt : "48" (default) or "32" for the Vaisala decoder.
    prevalence_threshold : a flag counts toward a clock-hour only if present in > this fraction of
                    that hour's profiles (default 5%); the reported count is the number of such hours.

    Occurrence counting is in **prevalence-gated clock-hours** (0-24): the number of hours in which a
    code exceeded ``prevalence_threshold`` — matching the operator's "N hourly occurrences" and
    suppressing transient single-profile blips.
    """
    thresholds = thresholds or HKThresholds()
    itype = str(itype or "")

    seq = [] if times_epoch is None else list(times_epoch)
    times = [float(t) for t in seq if t is not None and math.isfinite(float(t))]
    n_profiles = len(times)
    if n_profiles == 0:
        return DayStatus(quality="nodata", summary="No data", n_profiles=0)

    hours_with_data = {_hour_of(t) for t in times}
    gap_hours = 24 - len(hours_with_data)

    interval = expected_interval_s or _infer_interval_s(times)
    coverage = float("nan")
    if interval and interval > 0:
        expected = 86400.0 / interval
        coverage = min(100.0, 100.0 * n_profiles / expected) if expected > 0 else float("nan")

    # --- decode each profile, tallying per-(clock-hour, flag) counts + per-hour profile totals ---
    hour_total: Dict[int, int] = {}          # hour -> # profiles (prevalence denominator)
    hour_flag: Dict[int, Dict[str, int]] = {}  # hour -> {flag: # profiles with it}
    flag_sev: Dict[str, str] = {}

    def _accumulate(decoded: List[Tuple[str, str]], hour: int) -> None:
        hour_total[hour] = hour_total.get(hour, 0) + 1
        if not decoded:
            return
        hf = hour_flag.setdefault(hour, {})
        for name, sev in decoded:
            if sev == Severity.INFO:
                continue  # informational (S)/Note -> ignored per user
            hf[name] = hf.get(name, 0) + 1
            flag_sev[name] = sev

    if cl61_status:
        names = list(cl61_status.keys())
        for i, t in enumerate(times):
            snapshot = {nm: cl61_status[nm][i] for nm in names
                        if i < len(cl61_status[nm])}
            _accumulate(decode_cl61_status(snapshot), _hour_of(t))
    elif status_values is not None:
        is_chm = itype.upper().startswith("CHM")
        n = min(len(times), len(status_values))
        for i in range(n):
            val = status_values[i]
            decoded = (decode_chm15k_error_ext(val) if is_chm
                       else decode_vaisala_error_string(val, fmt=error_string_fmt))
            _accumulate(decoded, _hour_of(times[i]))

    # Prevalence gate: a flag "occurred" in an hour only if it exceeds the threshold there; its
    # reported count is the number of such hours. Alarm/warning hours are the gated hours.
    flag_hours: Dict[str, set] = {}
    for h, fl in hour_flag.items():
        tot = hour_total.get(h, 0)
        if tot <= 0:
            continue
        for name, cnt in fl.items():
            if cnt / tot > prevalence_threshold:
                flag_hours.setdefault(name, set()).add(h)
    alarm_hours = {h for nm, hrs in flag_hours.items()
                   if flag_sev.get(nm) == Severity.ERROR for h in hrs}
    warn_hours = {h for nm, hrs in flag_hours.items()
                  if flag_sev.get(nm) == Severity.WARNING for h in hrs}
    flags = {name: len(hrs) for name, hrs in flag_hours.items()}

    # --- quality class = worst of flags / availability / housekeeping ---
    flag_class = "pass"
    if alarm_hours:
        flag_class = "error"
    elif warn_hours:
        flag_class = "warning"

    avail_class = "pass"
    if not math.isnan(coverage) and coverage < 90.0:
        avail_class = "warning"
    elif gap_hours >= 1 and interval is None:
        # no cadence to judge coverage, but whole clock-hours are missing
        avail_class = "warning"

    hk_notes = thresholds.evaluate(hk)
    hk_class = "warning" if hk_notes else "pass"

    quality = _worst(flag_class, avail_class, hk_class)

    # --- reporting summary string (warnings + alarms + gaps + HK; no info) ---
    parts: List[str] = []
    if gap_hours >= 1:
        parts.append(f"Missing data for {gap_hours} h")
    # alarms first, then warnings; each sorted by hours-present desc
    for want in (Severity.ERROR, Severity.WARNING):
        items = sorted(((nm, flags[nm]) for nm in flags if flag_sev.get(nm) == want),
                       key=lambda kv: (-kv[1], kv[0]))
        for nm, hrs in items:
            parts.append(f"'{nm}' {hrs} h")
    parts.extend(hk_notes)
    summary = " · ".join(parts) if parts else "No warning or error recorded"

    return DayStatus(
        quality=quality,
        summary=summary,
        n_profiles=n_profiles,
        coverage_pct=coverage,
        gap_hours=gap_hours,
        n_alarm_h=len(alarm_hours),
        n_warn_h=len(warn_hours),
        flags=flags,
    )
