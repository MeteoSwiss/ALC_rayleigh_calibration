"""Instrument status-string decoding + per-day health classification.

Decodes the per-profile status/error/alarm words that ceilometers report
(CHM15k ``error_ext`` 32-bit service code; Vaisala CL31/CL51 ``error_string``
48-bit / 32-bit status; CL61 ``/status`` subsystem group) into named
warnings/errors, and aggregates a day into a Cloudnet-style quality class
(pass / warning / error / nodata) for the monitoring dashboard.

See ``calibration/status/decode.py`` for the bit tables (manual-verified against
the Vaisala CL31 User's Guide M210482EN and the MeteoSwiss MATLAB reference
``plot_error_chm15k.m`` / ``FVaisala_CL31_Line2.m``).
"""

from .decode import (  # noqa: F401
    Severity,
    DayStatus,
    HKThresholds,
    decode_chm15k_error_ext,
    decode_vaisala_error_string,
    decode_cl61_status,
    decode_sci,
    summarize_day,
    QUALITY_ORDER,
)
