"""Slide 9 of the E-PROFILE Annual Meeting 2026 deck, 4th configuration: legacy Hopkin eta ladder
WITH the operational WV correction (CL31 909.7/6.0 nm = code default; CL61 910.55/0.188 nm, as in
ops/config.sh ALC_WV_SPECTRUM), on the same Payerne day list as slide9_eta_ladder.py.

Together with hopkin_nowv / nowv / operational this completes the 2x2 {Hopkin, PVC} x {no WV, WV},
so the slide can compare the eta tables in the operational (WV-corrected) setup.

Run from the repo root:  python presentation/slide9_eta_ladder_wv.py [--workers 12]
"""
import argparse
import os
import shutil
import sys
from pathlib import Path

# Set BEFORE importing cloud_profile_dump (it reads them at import) -- also inherited by workers.
os.environ["ALC_DUMP_UNITS"] = "0-20000-0-06610_B,0-20000-0-06610_C"
os.environ["ALC_DUMP_TAG"] = "hopkin_wv"
os.environ.pop("ALC_WV_DISABLE", None)
os.environ["ALC_ETA_TABLES"] = "hopkin_legacy"
os.environ["ALC_WV_SPECTRUM"] = '{"CL61": [910.55, 0.188]}'

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / "rayleigh_availability"))

# The dump's CAMS_FOLDERS points at A:/CAMS_Monthly_04, which no longer exists on this PC (the 0.4-deg
# monthlies are on D:). Overridden at MODULE level so the spawned workers (which re-run this module as
# __mp_main__) get it too. Never D:/CAMS or D:/CAMS_run_v20 (1-degree).
import cloud_profile_dump as _D  # noqa: E402
_D.CAMS_FOLDERS = "D:/CAMS_Monthly_04;D:/CAMS_daily"

DATA = Path("C:/DATA/Projects/202606_E-PROFILE_calibration")
OUT = REPO / "presentation" / "figs_annual2026" / "eta_ladder"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=12)
    a = ap.parse_args()
    import cloud_profile_dump as D
    assert str(D.DUMP_DIR).endswith("cloud_profile_dump_hopkin_wv"), D.DUMP_DIR
    from calibration.cloud._filters import _ETA_LEGACY_CL31  # noqa: F401  (switch present)
    from calibration.water_vapor_correction.water_vapor import laser_spectrum_for
    print("CL31 spectrum:", laser_spectrum_for("CL31", 910.0), " CL61 spectrum:", laser_spectrum_for("CL61", 910.0),
          flush=True)
    print("CAMS_FOLDERS:", D.CAMS_FOLDERS, flush=True)
    D.run_dump(a.workers)
    OUT.mkdir(parents=True, exist_ok=True)
    for i in ("B", "C"):
        f = DATA / "cloud_profile_dump_hopkin_wv" / f"0-20000-0-06610_{i}_profiles.npz"
        print(f"hopkin_wv {i}  exists={f.exists()}  {f}", flush=True)
        if f.exists():
            shutil.copy2(f, OUT / f"hopkin_wv_06610_{i}.npz")
    print("SLIDE9WV_DONE", flush=True)


if __name__ == "__main__":
    main()
