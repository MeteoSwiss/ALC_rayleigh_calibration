# -*- coding: utf-8 -*-
"""Is the propagated night-mean noise (calibrate_rayleigh: sigma_signal) right?

sigma_signal = RMS(native per-gate sigma in the fit gate) / sqrt(native gates per fit gate)
               / sqrt(n_profiles), with n_profiles = the number of BINNED (300 s) profiles kept.
The native sigma (temporal first differences / sqrt 2) is the noise of ONE native profile (15 s
CHM15k), so dividing by the number of 300-s profiles misses the time binning (x sqrt(20) on CHM15k).

Empirical reference: the same first-difference estimator applied to the BINNED profile stack
(consecutive 300-s profiles; the atmosphere is steady over 5 min at 3-6 km), divided by
sqrt(n_binned) -> noise of the night mean. Ratio propagated / empirical per gate, median 3-6 km.
"""
from __future__ import annotations

import sys
import warnings
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from u23_validation import OUT, _info, _install, _options, census  # noqa: E402

NIGHTS = [("0-20000-0-06240", "B", "20250412"), ("0-20000-0-06240", "B", "20250316"),
          ("0-20000-0-10393", "0", "20250425"), ("0-20000-0-06610", "A", "20250409"),
          ("0-20000-0-06610", "A", "20260912"), ("0-20008-0-UGR", "A", "20250330"),
          ("0-20000-0-10393", "C", "20250425"), ("0-20000-0-07014", "A", "20250425")]


def main():
    import logging
    from calibration import calibrate_rayleigh
    from calibration.io.data_loader import build_file_paths, load_data
    from calibration.sensitivity.noise import first_difference_sigma
    warnings.filterwarnings("ignore")
    logging.disable(logging.CRITICAL)
    _install()
    cen = census()
    rows = []
    for wmo, ident, d in NIGHTS:
        s = cen[(wmo, ident)]
        info = _info(s)
        o = _options({}, "ensemble")
        files = [f for f in build_file_paths(d, info, o) if f.exists()]
        native = load_data(files, info.instrument_type, o.data_level)
        fi = {}
        r = calibrate_rayleigh(d, info, o, preloaded_data=native, fit_inputs_out=fi)
        if "signal_stack" not in fi or fi.get("sigma_signal") is None:
            print(wmo, ident, d, "no fit inputs (flag", r.flag, ")")
            continue
        stack = np.asarray(fi["signal_stack"], float)
        hours = np.asarray(fi["hours"], float)
        rng = np.asarray(fi["range_alc"], float)
        prop = np.asarray(fi["sigma_signal"], float)
        dt_b = float(np.median(np.diff(hours))) * 3600.0
        sig_b, npair = first_difference_sigma(stack, hours * 3600.0, dt_b)
        emp = sig_b / np.sqrt(stack.shape[0])
        band = (rng >= 3000) & (rng <= 6000) & np.isfinite(prop) & np.isfinite(emp) & (emp > 0)
        ratio = float(np.median(prop[band] / emp[band]))
        sn = fi.get("sigma_native") or {}
        dt_n = float(sn.get("dt_s", np.nan))
        k = dt_b / dt_n if np.isfinite(dt_n) and dt_n > 0 else np.nan
        rows.append((wmo, ident, d, s["type"], r.flag, stack.shape[0], dt_n, dt_b, k, ratio))
        print(f"{s['site'][:14]:14s} {ident} {d} {s['type']:8s} flag {r.flag:4} n_binned "
              f"{stack.shape[0]:3d} dt native {dt_n:5.1f} s binned {dt_b:5.0f} s (k={k:5.1f}, "
              f"sqrt k={np.sqrt(k):4.2f}) | propagated / empirical noise of the night mean, "
              f"3-6 km: {ratio:5.2f}", flush=True)
    np.save(OUT / "u23_noise_check.npy", np.array(rows, dtype=object), allow_pickle=True)


if __name__ == "__main__":
    main()
