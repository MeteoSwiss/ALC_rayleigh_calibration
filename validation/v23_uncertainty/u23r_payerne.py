# -*- coding: utf-8 -*-
"""Review (claim 4): Payerne A, 2026-09-07 and its neighbours, eprof_v2.2 AND eprof_v2.3.

Recalibrates the four nights of doc section 3.8 with the classification screen (curtains computed
by u23_analysis.cmd_payerne with the runner's own classification, read here read-only), for both
window methods, with the ensemble (v2.2) uncertainty, the shipped v2.3 budget and any extra
uncertainty_params variants given on the command line as JSON.

Usage: python u23r_payerne.py ['{"label": {...params...}}']
"""
from __future__ import annotations

import json
import sys
import warnings
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import u23_validation as V  # noqa: E402

NIGHTS = ("20260906", "20260907", "20260912", "20260923")
CLS_ROOT = Path("C:/DATA/Projects/202610_v23/uncertainty/payerne_cls")
UNIT_MEDIAN_2026 = 5.26e11          # doc section 3.8 (published 2026 unit median)


def main():
    import logging
    from calibration import calibrate_rayleigh
    from calibration.io.classification import (classification_files_for_night,
                                               read_classification_curtain)
    from calibration.io.data_loader import build_file_paths, load_data
    warnings.filterwarnings("ignore")
    logging.disable(logging.CRITICAL)
    V._install()
    extra = json.loads(sys.argv[1]) if len(sys.argv) > 1 else {}
    s = V.census()[("0-20000-0-06610", "A")]
    info = V._info(s)
    key = "0-20000-0-06610_A"
    pre_review = dict(noise_estimator="per_gate", molecular_height_reference="agl",
                      vertical_cap_rel=None)
    models = [("ensemble", "ensemble", {}), ("v2.3 pre-review", "v2.3", pre_review),
              ("v2.3 reviewed", "v2.3", {})]
    models += [(k, "v2.3", v) for k, v in extra.items()]
    lines = []
    for d in NIGHTS:
        o0 = V._options({}, "ensemble")
        files = [f for f in build_file_paths(d, info, o0) if f.exists()]
        native = load_data(files, info.instrument_type, o0.data_level)
        cls = read_classification_curtain(classification_files_for_night(CLS_ROOT, key, s["wmo"], d))
        for method in ("eprof_v2.2", "eprof_v2.3"):
            for label, model, params in models:
                o = V._options(params, model)
                o.molecular_method = method
                fi = {}
                r = calibrate_rayleigh(d, info, o, preloaded_data=native, fit_inputs_out=fi,
                                       classification=cls)
                b = fi.get("uncertainty_budget")
                c = b.lidar_constant if b is not None else r.lidar_constant
                txt = (f"{d} {method:10s} {label:22s} flag {r.flag:5} C/median "
                       f"{c / UNIT_MEDIAN_2026 if c > 0 else np.nan:5.2f}")
                if b is not None:
                    txt += (f" | 2s total {100 * b.rel_total:5.0f} % (ens {200 * b.sigma_ensemble_rel:4.0f}, "
                            f"vert {200 * b.sigma_vertical_rel:4.0f}, mol {200 * b.sigma_molecular_rel:4.1f}, "
                            f"noise {200 * b.sigma_noise_rel:4.1f}) g {100 * b.gradient_per_km:+5.0f} %/km "
                            f"n_sub {b.n_subwindows} {b.z_lo_m:.0f}-{b.z_hi_m:.0f} m")
                elif r.flag in (1, 0.5):
                    txt += f" | 2s {100 * r.uncertainty / r.lidar_constant:5.0f} %"
                txt += f" | {r.message[:60]}"
                print(txt, flush=True)
                lines.append(txt)
    out = Path("C:/DATA/Projects/202610_v23/uncertainty_review/u23r_payerne.txt")
    out.write_text("\n".join(lines), encoding="utf-8")


if __name__ == "__main__":
    main()
