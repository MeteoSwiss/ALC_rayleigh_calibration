"""Cloudnet-style target classification (``ceiloclass``) on the E-PROFILE ALC network.

Bridges the E-PROFILE L1 archive to the ``actris-cloudnet/ceiloclass`` classifier:

* :func:`read_eprofile_l1` / :func:`ceilo_from_shared` — L1 rcs_0 (file or the read-once shared
  working grid) -> ceilopyter ``Ceilo`` (screened attenuated backscatter + CL61 depolarization);
* :func:`cams_to_model` — the operational CAMS temperature -> ceiloclass ``Model``.

``ceiloclass`` + ``ceilopyter`` are optional dependencies (not needed for calibration); importing
this subpackage requires them.
"""

from .cams_model import cams_to_model
from .eprofile_l1 import (
    ceilo_from_shared,
    instrument_type,
    read_eprofile_l1,
    station_altitude,
)

__all__ = [
    "cams_to_model",
    "ceilo_from_shared",
    "instrument_type",
    "read_eprofile_l1",
    "station_altitude",
]
