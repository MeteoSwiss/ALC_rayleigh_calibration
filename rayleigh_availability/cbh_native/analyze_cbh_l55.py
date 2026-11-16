"""`_l55` re-run of analyze_cbh.py on the CL61 constants re-expressed at 910.55 nm / 0.188 nm.

Input:  C:/DATA/Projects/202606_E-PROFILE_calibration/cbh_native_l55/cbh_scenes.csv
        (from rescale_cl61_l55.py; CL31/CL51 rows unchanged).
Output: per_flux_slopes.csv, level_bias.csv, fe_results.json in that same folder (the originals in
        rayleigh_availability/cbh_native/ are untouched); printed tables = phase4 section 3 numbers.
Run:    python rayleigh_availability/cbh_native/analyze_cbh_l55.py
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import analyze_cbh as A  # noqa: E402

A.SCR = r"C:/DATA/Projects/202606_E-PROFILE_calibration/cbh_native_l55"

if __name__ == "__main__":
    import matplotlib
    matplotlib.use("Agg")
    from matplotlib.figure import Figure
    _orig = Figure.savefig

    def _savefig(self, fname, *a, **k):  # never overwrite an original figure
        from pathlib import Path
        p = Path(fname)
        if not p.stem.endswith("_l55"):
            p = p.with_name(f"{p.stem}_l55{p.suffix}")
        return _orig(self, p, *a, **k)
    Figure.savefig = _savefig
    A.main()
