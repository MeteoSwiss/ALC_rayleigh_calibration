"""Slide 9 of the E-PROFILE Annual Meeting 2026 deck: Payerne cloud calibration C_L vs cloud-base
height under three configurations (per-profile dumps, same day list = cbh_scenes.csv):

  1. hopkin_nowv : legacy Hopkin eta ladder (a_G 8 um), no WV correction      -> computed here
  2. nowv        : operational PVC eta (a_G 5.5 um), no WV correction          -> existing dump
  3. operational : PVC eta + WV (CL31 909.7/6.0 = tag 'nominal';
                   CL61 910.55/0.188 = tag 'l55s008')                          -> existing dumps

Only the dump is run (NOT cloud_profile_dump.main, whose plots would overwrite the report figure
and the gallery). The npz of all three configurations are then copied to
presentation/figs_annual2026/eta_ladder/ for the figure.

Run from the repo root:  python presentation/slide9_eta_ladder.py [--workers 12]
"""
import argparse
import os
import shutil
import sys
from pathlib import Path

# Set BEFORE importing cloud_profile_dump (it reads them at import) -- also inherited by workers.
os.environ["ALC_DUMP_UNITS"] = "0-20000-0-06610_B,0-20000-0-06610_C"
os.environ["ALC_DUMP_TAG"] = "hopkin_nowv"
os.environ["ALC_WV_DISABLE"] = "1"
os.environ["ALC_ETA_TABLES"] = "hopkin_legacy"
os.environ.pop("ALC_WV_SPECTRUM", None)

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / "rayleigh_availability"))

DATA = Path("C:/DATA/Projects/202606_E-PROFILE_calibration")
OUT = REPO / "presentation" / "figs_annual2026" / "eta_ladder"
SRC = {"hopkin_nowv": ["B", "C"], "nowv": ["B", "C"], "nominal": ["B"], "l55s008": ["C"]}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=12)
    ap.add_argument("--copy-only", action="store_true")
    a = ap.parse_args()
    if not a.copy_only:
        import cloud_profile_dump as D
        assert str(D.DUMP_DIR).endswith("cloud_profile_dump_hopkin_nowv"), D.DUMP_DIR
        from calibration.cloud._filters import _ETA_LEGACY_CL31  # noqa: F401  (switch present)
        D.run_dump(a.workers)
    OUT.mkdir(parents=True, exist_ok=True)
    for tag, idents in SRC.items():
        for i in idents:
            f = DATA / f"cloud_profile_dump_{tag}" / f"0-20000-0-06610_{i}_profiles.npz"
            print(f"{tag:12s} {i}  exists={f.exists()}  {f}", flush=True)
            if f.exists():
                shutil.copy2(f, OUT / f"{tag}_06610_{i}.npz")
    print("SLIDE9_DONE", flush=True)


if __name__ == "__main__":
    main()
