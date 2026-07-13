"""Unit tests for calibration.status.decode — CHM15k/CL31/CL51/CL61 status decoding
and the per-day quality classification.

Bit meanings are manual-verified (Vaisala CL31 UG M210482EN; MeteoSwiss MATLAB
plot_error_chm15k.m / FVaisala_CL31_Line2.m)."""
import math

from calibration.status.decode import (
    Severity,
    HKThresholds,
    decode_chm15k_error_ext,
    decode_vaisala_error_string,
    decode_cl61_status,
    decode_sci,
    summarize_day,
)


# --------------------------------------------------------------------------- #
# CHM15k error_ext
# --------------------------------------------------------------------------- #
def test_chm15k_clean_and_invalid():
    assert decode_chm15k_error_ext(0) == []
    assert decode_chm15k_error_ext(float("nan")) == []
    assert decode_chm15k_error_ext(-999) == []
    assert decode_chm15k_error_ext(None) == []


def test_chm15k_known_bits():
    # bit 17 = Windows contaminated (warning)
    got = decode_chm15k_error_ext(1 << 17)
    assert got == [("Warning: Windows contaminated", Severity.WARNING)]
    # bit 10 = Laser optical unit temperature error (error)
    got = decode_chm15k_error_ext(1 << 10)
    assert got == [("Error: Laser optical unit temperature error", Severity.ERROR)]
    # bit 30 = Standby mode on (info/Note)
    assert decode_chm15k_error_ext(1 << 30) == [("Note: Standby mode on", Severity.INFO)]
    # float holding an int decodes the same
    assert decode_chm15k_error_ext(float(1 << 17))[0][1] == Severity.WARNING


def test_chm15k_multiple_bits():
    got = dict(decode_chm15k_error_ext((1 << 10) | (1 << 17)))
    assert got["Error: Laser optical unit temperature error"] == Severity.ERROR
    assert got["Warning: Windows contaminated"] == Severity.WARNING


# --------------------------------------------------------------------------- #
# Vaisala CL31/CL51 error_string — 48-bit + 32-bit
# --------------------------------------------------------------------------- #
def test_vaisala_48bit_known_bits():
    assert decode_vaisala_error_string(0) == []
    # bit 46 = Transmitter failure (A) -> error
    assert decode_vaisala_error_string(1 << 46) == [("Transmitter failure (A)", Severity.ERROR)]
    # bit 29 = Transmitter expires (W) -> warning
    assert decode_vaisala_error_string(1 << 29) == [("Transmitter expires (W)", Severity.WARNING)]
    # bit 31 = Window contamination (W)
    assert decode_vaisala_error_string(1 << 31) == [("Window contamination (W)", Severity.WARNING)]


def test_vaisala_48bit_real_sample_is_status_only():
    # Real Payerne-style sample value 0xC020 = bits 5,14,15 -> all internal status (S)
    decoded = decode_vaisala_error_string(0xC020)
    assert decoded, "expected some (S) bits set"
    assert all(sev == Severity.INFO for _, sev in decoded)


def test_vaisala_hex_string_input():
    # E-PROFILE L1 form: raw hex string. Length selects the format.
    assert decode_vaisala_error_string("000000000000") == []
    assert decode_vaisala_error_string("00000000C080"), "status bits should decode"
    assert all(sev == Severity.INFO for _, sev in decode_vaisala_error_string("00000000C080"))
    # 12-hex with Transmitter failure (A) = bit 46 -> 0x400000000000
    assert decode_vaisala_error_string("400000000000") == [
        ("Transmitter failure (A)", Severity.ERROR)
    ]
    # 8-hex string auto-detects the 32-bit map: bit 30 = Transmitter failure (A)
    assert decode_vaisala_error_string("40000000") == [
        ("Transmitter failure (A)", Severity.ERROR)
    ]
    # bytes + whitespace tolerated
    assert decode_vaisala_error_string(b"  000000000000  ") == []


def test_summarize_accepts_hex_strings():
    times = _times([1, 2, 3])
    vals = ["400000000000"] * len(times)  # Transmitter failure (A)
    ds = summarize_day(times, "CL31", status_values=vals)
    assert ds.quality == "error"
    assert "Transmitter failure (A)" in ds.summary


def test_vaisala_32bit_format_differs():
    # In the 32-bit CT25K map, bit 30 = Transmitter failure (A)
    assert decode_vaisala_error_string(1 << 30, fmt="32") == [
        ("Transmitter failure (A)", Severity.ERROR)
    ]
    # bit 23 = Window contaminated (W)
    assert decode_vaisala_error_string(1 << 23, fmt="32") == [
        ("Window contaminated (W)", Severity.WARNING)
    ]


# --------------------------------------------------------------------------- #
# CL61 /status group
# --------------------------------------------------------------------------- #
def test_cl61_status_traffic_light():
    assert decode_cl61_status({"Transmitter_overall": 0}) == []
    assert decode_cl61_status({"Window_condition": 1}) == [("Window condition", Severity.WARNING)]
    assert decode_cl61_status({"Transmitter_overall": 2}) == [
        ("Transmitter overall", Severity.ERROR)
    ]
    # unknown non-zero code -> at least a warning
    assert decode_cl61_status({"Foo_bar": 7})[0][1] == Severity.WARNING


# --------------------------------------------------------------------------- #
# sky-condition index
# --------------------------------------------------------------------------- #
def test_decode_sci():
    assert decode_sci(0) == "nothing"
    assert decode_sci(2) == "fog"
    assert decode_sci(4) == "precipitation or particles on window"
    assert decode_sci(float("nan")) is None


# --------------------------------------------------------------------------- #
# summarize_day
# --------------------------------------------------------------------------- #
DAY0 = 1706745600  # 2024-02-01 00:00:00 UTC


def _times(hours, per_hour=1):
    out = []
    for h in hours:
        for k in range(per_hour):
            out.append(DAY0 + h * 3600 + k * 60)
    return out


def _uniform(n_hours, interval_s=900):
    """Uniform cadence over the first ``n_hours`` of the day (realistic sampling)."""
    return [DAY0 + s for s in range(0, int(n_hours * 3600), interval_s)]


def test_summarize_nodata():
    ds = summarize_day([], "CL31")
    assert ds.quality == "nodata"
    assert ds.summary == "No data"
    ds2 = summarize_day(None, "CHM15k", status_values=None)
    assert ds2.quality == "nodata"


def test_summarize_clean_full_day():
    times = _uniform(24)  # 15-min cadence, full coverage
    vals = [0] * len(times)
    ds = summarize_day(times, "CL31", status_values=vals)
    assert ds.quality == "pass"
    assert ds.summary == "No warning or error recorded"
    assert ds.gap_hours == 0
    assert ds.flags == {}
    assert ds.coverage_pct > 95


def test_summarize_warning_flag_counts_hours():
    # Transmitter expires (W) present in 3 distinct hours
    hours = [2, 3, 4]
    times = _times(hours)
    vals = [1 << 29] * len(times)
    ds = summarize_day(times, "CL51", status_values=vals)
    assert ds.flags == {"Transmitter expires (W)": 3}
    assert ds.n_warn_h == 3
    assert ds.quality == "warning"
    assert "'Transmitter expires (W)' 3 h" in ds.summary


def test_summarize_alarm_is_error():
    times = _times([10, 11])
    vals = [1 << 46] * len(times)  # Transmitter failure (A)
    ds = summarize_day(times, "CL31", status_values=vals)
    assert ds.quality == "error"
    assert ds.n_alarm_h == 2
    assert "'Transmitter failure (A)' 2 h" in ds.summary


def test_summarize_info_bits_ignored():
    # Blower on (S) bit 15 all day -> still pass, nothing reported
    times = _uniform(24)
    vals = [1 << 15] * len(times)
    ds = summarize_day(times, "CL31", status_values=vals)
    assert ds.quality == "pass"
    assert ds.flags == {}
    assert ds.summary == "No warning or error recorded"


def test_summarize_missing_data_reported():
    times = _uniform(18)  # 18 h present, 6 h missing
    vals = [0] * len(times)
    ds = summarize_day(times, "CHM15k", status_values=vals)
    assert ds.gap_hours == 6
    assert ds.quality == "warning"
    assert "Missing data for 6 h" in ds.summary


def test_summarize_hk_threshold_warns():
    times = _uniform(24)
    vals = [0] * len(times)
    ds = summarize_day(times, "CL61", status_values=vals, hk={"window": 30.0})
    assert ds.quality == "warning"
    assert "Low window transmission" in ds.summary


def test_summarize_chm_uses_error_ext_table():
    times = _times([5, 6])
    vals = [1 << 10] * len(times)  # CHM error bit
    ds = summarize_day(times, "CHM15k", status_values=vals)
    assert ds.quality == "error"
    assert any("Laser optical unit temperature error" in k for k in ds.flags)


def test_prevalence_gate_transient_does_not_warn_day():
    # Full day at 60 s cadence (1440 profiles, 60/hour). A single transient warning profile
    # (1/60 = 1.7% of its hour, < 5%) is gated out -> the day stays 'pass'.
    times = _uniform(24, interval_s=60)
    vals = [0] * len(times)
    vals[100] = 1 << 29                      # one profile, Transmitter expires (W)
    ds = summarize_day(times, "CL31", status_values=vals)
    assert ds.quality == "pass"
    assert ds.flags == {}
    assert ds.n_warn_h == 0
    assert ds.summary == "No warning or error recorded"


def test_prevalence_gate_sustained_warns_day():
    # 10 of 60 profiles in one hour (16.7% > 5%) is kept -> that hour counts.
    times = _uniform(24, interval_s=60)
    vals = [0] * len(times)
    for k in range(10):
        vals[5 * 60 + k] = 1 << 46           # Transmitter failure (A), sustained in hour 5
    ds = summarize_day(times, "CL31", status_values=vals)
    assert ds.quality == "error"
    assert ds.flags.get("Transmitter failure (A)") == 1   # exactly 1 gated hour
    assert ds.n_alarm_h == 1
    assert "'Transmitter failure (A)' 1 h" in ds.summary


def test_prevalence_threshold_is_tunable():
    times = _uniform(24, interval_s=60)
    vals = [0] * len(times)
    vals[100] = 1 << 29                      # 1.7% of its hour
    # with a 1% threshold the same blip now counts
    ds = summarize_day(times, "CL31", status_values=vals, prevalence_threshold=0.01)
    assert ds.flags.get("Transmitter expires (W)") == 1


def test_summarize_accepts_numpy_arrays():
    np = __import__("numpy")
    times = np.array(_uniform(24), dtype=float)
    vals = np.zeros(len(times))
    ds = summarize_day(times, "CL31", status_values=vals)
    assert ds.quality == "pass"
    assert ds.n_profiles == len(times)


def test_hk_thresholds_healthy():
    assert HKThresholds().evaluate({"window": 95, "laser": 100, "temp_internal": 20}) == []
    assert HKThresholds().evaluate(None) == []
