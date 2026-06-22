"""Generate a few curated flag-example diagnostics for the dashboard flag page (flags.html).

Output files are named '<anchor>__<caption>.png' so render._copy_flag_examples picks them up:
  1__...   -> flag 1  (success)
  m1__...  -> flag -1 (unsuitable / no liquid cloud)
Run:  python scripts/make_flag_examples.py [OUTDIR]
"""
import json
import sys
from datetime import datetime, timedelta
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))          # for run_all_l1_2026
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))      # repo root

import numpy as np
import run_all_l1_2026 as R  # noqa: E402  (reuse its paths + per-stream helpers)
from calibration.cloud.calibration import (  # noqa: E402
    CloudCalConfig, liquid_cloud_calibration_from_data, read_ceilometer_data, set_defaults)
from calibration.plotting import plot_cloud_diagnostics_compact  # noqa: E402

OUTDIR = Path(sys.argv[1]) if len(sys.argv) > 1 else Path(
    "C:/DATA/Projects/202606_E-PROFILE_calibration/flag_examples")
OUTDIR.mkdir(parents=True, exist_ok=True)
census = json.loads(R.CENSUS.read_text(encoding="utf-8"))


def _cfg(s, fp):
    return set_defaults(CloudCalConfig(
        nc_file=str(fp), instrument=s["type"], apply_wv_correction=True,
        apply_transmission_correction=True, aerosol_lidar_ratio=50.0,
        cams_folder=str(R.CAMS), abs_cs_lookup_table=str(R.WV_LUT),
        station_latitude=s["lat"], station_longitude=s["lon"],
        average_time_s=300.0, average_range_m=10.0))


def _run(s, d):
    """Return (data, res) for one station-day, or None if the file/data is unusable."""
    fp = R._l1_file(s["wmo"], s["ident"], d)
    if not fp.exists():
        return None
    cfg = _cfg(s, fp)
    data, status = read_ceilometer_data(cfg.nc_file, cfg)
    beta = getattr(data, "beta", None) if data is not None else None
    if status != 0 or beta is None or not np.any(np.isfinite(np.asarray(beta, dtype=float))):
        return None
    return data, liquid_cloud_calibration_from_data(data, cfg)


# 1) flag 1 (success): KLEINE_BROGEL 0-20000-0-06479 on 2026-03-01 (the worked example).
s = next(x for x in census if R._key(x).startswith("0-20000-0-06479"))
out = _run(s, datetime(2026, 3, 1))
if out and int(getattr(out[1], "n_profiles", 0)) > 0:
    data, res = out
    site = s.get("site", R._key(s))
    plot_cloud_diagnostics_compact(
        data, res, title=f"{site} — 20260301 — cloud SUCCESS (flag 1)",
        save_path=OUTDIR / "1__cloud_calibration_success_(KLEINE-BROGEL_2026-03-01).png")
    print(f"[1] cloud success: {R._key(s)} n_profiles={res.n_profiles}", flush=True)
else:
    print("[1] could not build cloud-success example", flush=True)

# 2) flag -1 (no liquid cloud): first day with data but zero selected profiles.
done = False
for s in [x for x in census if x.get("type") in ("CL31", "CL51", "CL61")]:
    if done:
        break
    for k in range(0, 75):
        out = _run(s, datetime(2026, 3, 1) + timedelta(days=k))
        if not out:
            continue
        data, res = out
        if int(getattr(res, "n_profiles", 0)) == 0:
            d = datetime(2026, 3, 1) + timedelta(days=k)
            site = s.get("site", R._key(s))
            plot_cloud_diagnostics_compact(
                data, res, title=f"{site} — {d:%Y%m%d} — NO liquid cloud in view (flag -1)",
                save_path=OUTDIR / "m1__no_liquid_cloud_in_view_(data_present,_clear_or_ice_only).png")
            print(f"[-1] no-liquid-cloud: {R._key(s)} {d:%Y%m%d}", flush=True)
            done = True
            break
if not done:
    print("[-1] could not build a no-liquid-cloud example", flush=True)
print("DONE", flush=True)
