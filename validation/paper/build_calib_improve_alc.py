"""build_calib_improve_alc.py - daily best-estimate lidar constants for the paper (manuscript v0.9 onwards).

The paper uses the nightly calibration constants that the operational dashboard publishes (Rayleigh
eprof_v2.2 for the CHM15k / Mini-MPL / CL61, O'Connor liquid-cloud for the Vaisala ceilometers) and
smooths them with the OPERATIONAL E-Profile Kalman filter improve_alc_calib
(validation/paper/kalman_improve_alc.py). The dashboard's own Kalman curve is for display only and is
never used for paper numbers (operator decision 2026-10-10).

Input  : <calib_dashboard>/<WMO>_<ident>_<rayleigh|cloud>.csv  (time, C_daily, C_daily_std, C_kalman,
         C_kalman_std), written by validation/paper/fetch_dashboard_calib.py and
         rayleigh_availability/background_bias_network.py (network CHM15k), fetched 2026-10-09.
Output : <calib_improve_alc>/<same file name>.csv with the same columns; C_daily/C_daily_std copied,
         C_kalman/C_kalman_std replaced by the improve_alc state and its 1-sigma uncertainty, plus
         SOURCE.txt (one line per stream: type, nights, outliers, model, filter segments).

Filter details (see kalman_improve_alc.py): one global 1.5 IQR outlier test, annual-cycle prediction
model when > 20 valid nights (constant otherwise), one scalar measurement variance, type-specific
process noise. The filter is restarted at known optical-module replacements (SWAPS below), because the
Vaisala L1 files carry no optical_module_id the operational driver could use.

Streams with fewer than 5 valid nights are skipped. The instrument type comes from the L1 census
(validation/scope_l1_2026_census.json); streams missing from it default to CHM15k (Rayleigh) or CL31
(cloud).

Usage (repository root, project venv):
    python -m validation.paper.build_calib_improve_alc [<calib_dashboard_dir> [<output_dir>]]

History: this is the code that built paper_v09/calib_improve_alc on 2026-10-10 (run then as an inline
snippet, transcribed here unchanged apart from layout and the command-line arguments).
"""
from __future__ import annotations

import json
import re
import sys
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd

from validation.paper.kalman_improve_alc import kalman_improve_alc

PAPER = Path("C:/DATA/Projects/202606_E-PROFILE_calibration/figs_paper_validation")
DEFAULT_SRC = PAPER / "paper_v08" / "calib_dashboard"
DEFAULT_DST = PAPER / "paper_v09" / "calib_improve_alc"
CENSUS = Path(__file__).resolve().parents[1] / "scope_l1_2026_census.json"

# known optical-module replacements: the filter is restarted at each date
SWAPS = {"0-20000-0-06610_B": [datetime(2026, 7, 7)]}       # Payerne CL31, 2026-07-07 ~13:00 UTC
MIN_NIGHTS = 5


def instrument_types(census_path=CENSUS):
    """{'<WMO>_<ident>': instrument type} from the L1 census."""
    census = json.loads(Path(census_path).read_text(encoding="utf-8"))
    return {f"{e['wmo']}_{e['ident']}": e.get("type") for e in census}


def smooth_stream(path, itype, swaps):
    """Run improve_alc on one dashboard CSV; return (output frame, n_nights, n_outliers, models)."""
    d = pd.read_csv(path, parse_dates=["time"])
    valid = d.C_daily.notna()
    out = pd.DataFrame({"time": d.time, "C_daily": d.C_daily, "C_daily_std": d.C_daily_std,
                        "C_kalman": np.nan, "C_kalman_std": np.nan})
    bounds = [datetime(1900, 1, 1)] + list(swaps) + [datetime(2100, 1, 1)]
    n_outliers = 0
    models = []
    for start, end in zip(bounds[:-1], bounds[1:]):
        seg = valid & (d.time >= start) & (d.time < end)
        if seg.sum() < MIN_NIGHTS:
            continue
        r = kalman_improve_alc(list(d.time[seg]), d.C_daily[seg].values, itype)
        state = pd.Series(r["state"], index=pd.to_datetime(r["grid"]))
        sigma = pd.Series(r["std"], index=state.index)
        sel = (out.time >= max(start, state.index[0])) & (out.time < end)
        out.loc[sel, "C_kalman"] = state.reindex(out.time[sel]).values
        out.loc[sel, "C_kalman_std"] = sigma.reindex(out.time[sel]).values
        n_outliers += len(r["outlier_dates"])
        models.append(r["model"])
    # prediction-only days at the end of the dashboard file keep the last state
    out["C_kalman"] = out["C_kalman"].ffill()
    out["C_kalman_std"] = out["C_kalman_std"].ffill()
    return out, int(valid.sum()), n_outliers, models


def main(argv):
    src = Path(argv[0]) if len(argv) > 0 else DEFAULT_SRC
    dst = Path(argv[1]) if len(argv) > 1 else DEFAULT_DST
    dst.mkdir(parents=True, exist_ok=True)
    types = instrument_types()
    n_ok = n_skip = 0
    log = []
    for f in sorted(src.glob("*_rayleigh.csv")) + sorted(src.glob("*_cloud.csv")):
        key = re.sub(r"_(rayleigh|cloud)\.csv$", "", f.name)
        method = f.name.rsplit("_", 1)[1][:-4]
        nights = pd.read_csv(f).C_daily.notna().sum()
        if nights < MIN_NIGHTS:
            n_skip += 1
            continue
        itype = types.get(key) or ("CHM15k" if method == "rayleigh" else "CL31")
        out, n, n_out, models = smooth_stream(f, itype, SWAPS.get(key, []))
        out.to_csv(dst / f.name, index=False, float_format="%.6e", date_format="%Y-%m-%d", na_rep="")
        n_ok += 1
        log.append(f"{f.name}: {itype}, {n} nights, {n_out} outliers, model {'/'.join(models)}, "
                   f"segments {len(models)}")
    (dst / "SOURCE.txt").write_text(
        "Nightly values (C_daily) = dashboard v2.2 / O'Connor constants fetched 2026-10-09 (paper_v08/calib_dashboard).\n"
        "C_kalman = improve_alc_calib operational Kalman (validation/paper/kalman_improve_alc.py), global IQR outliers, "
        "annual-cycle model\nwhen > 20 nights; Payerne CL31 (0-20000-0-06610_B) restarted at the optical-module "
        "replacement of 2026-07-07.\n" + "\n".join(log) + "\n", encoding="utf-8")
    print(f"written {n_ok} streams, skipped {n_skip} (< {MIN_NIGHTS} nights) -> {dst}")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
