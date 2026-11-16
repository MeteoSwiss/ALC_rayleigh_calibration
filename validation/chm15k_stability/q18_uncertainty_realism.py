# -*- coding: utf-8 -*-
"""Is the published Rayleigh uncertainty realistic? (report 19, q18)

The calibration writes `lidar_constant_uncertainty` = 2 x robust sigma of a 125-member ensemble
(5 lidar ratios x 5 window shifts of +-200 m x 5 time samples). For every CHM15k unit (main module
era) this compares the reported sigma (u/2) with the observed night-to-night scatter of ln C around a
31-day running median, and computes the pull z = (ln C - running median) / (u / 2C). A realistic
uncertainty gives a robust std(z) of ~1.
Output: C:/DATA/Projects/202610_CHM15k_stability/q18/q18_uncertainty_units.csv
"""
import csv
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import OUT, chm15k_streams, load_constants, stream_key  # noqa: E402
from q1_temperature import ANOM_WINDOW, running_median  # noqa: E402

Q18 = OUT / "q18"


def robust_std(x):
    x = x[np.isfinite(x)]
    return 1.4826 * np.median(np.abs(x - np.median(x))) if x.size > 5 else np.nan


def main():
    rows = []
    for wmo, ident, site, lat, lon in chm15k_streams():
        cons = load_constants(wmo, ident)
        if cons is None or len(cons["c"]) < 40:
            continue
        mod = max(set(cons["module"]), key=list(cons["module"]).count)
        m = (cons["module"] == mod) & (cons["c"] > 0) & (cons["u"] > 0)
        t, c, u = cons["t"][m], cons["c"][m], cons["u"][m]
        o = np.argsort(t)
        t, lc, rel_sigma = t[o], np.log(c[o]), (u[o] / c[o]) / 2.0     # u is a 2-sigma value
        if t.size < 40:
            continue
        anom = lc - running_median(t, lc, ANOM_WINDOW)
        z = anom / rel_sigma
        rows.append(dict(stream=stream_key(wmo, ident), site=site, n=int(t.size),
                         reported_sigma_pct=100 * float(np.median(rel_sigma)),
                         observed_sigma_pct=100 * float(robust_std(anom)),
                         pull_std=float(robust_std(z)),
                         frac_outside_2sigma=float(np.mean(np.abs(z) > 2))))
    Q18.mkdir(parents=True, exist_ok=True)
    with open(Q18 / "q18_uncertainty_units.csv", "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    rep = np.array([r["reported_sigma_pct"] for r in rows])
    obs = np.array([r["observed_sigma_pct"] for r in rows])
    pull = np.array([r["pull_std"] for r in rows])
    out = np.array([r["frac_outside_2sigma"] for r in rows])
    print(f"{len(rows)} CHM15k units (main module era)")
    print(f"reported sigma  median {np.median(rep):.2f} % (IQR {np.percentile(rep, 25):.2f}-{np.percentile(rep, 75):.2f})")
    print(f"observed sigma  median {np.median(obs):.2f} % (IQR {np.percentile(obs, 25):.2f}-{np.percentile(obs, 75):.2f})")
    print(f"pull std        median {np.median(pull):.2f} (IQR {np.percentile(pull, 25):.2f}-{np.percentile(pull, 75):.2f}); "
          f"nights outside 2 sigma: median {100 * np.median(out):.0f} % (Gaussian 4.6 %)")
    print(f"Spearman(reported, observed) across units: {np.corrcoef(np.argsort(np.argsort(rep)), np.argsort(np.argsort(obs)))[0, 1]:+.2f}")
    for r in sorted(rows, key=lambda r: -r["observed_sigma_pct"])[:8]:
        print(f"  {r['site'][:22]:22s} reported {r['reported_sigma_pct']:5.2f} %  observed {r['observed_sigma_pct']:5.2f} %  pull {r['pull_std']:4.1f}")
    pay = [r for r in rows if r["stream"].startswith("0-20000-0-06610")]
    for r in pay:
        print(f"  PAYERNE {r['stream']}: reported {r['reported_sigma_pct']:.2f} %  observed {r['observed_sigma_pct']:.2f} %  pull {r['pull_std']:.1f}")


if __name__ == "__main__":
    main()
