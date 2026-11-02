"""Slide 10 of the E-PROFILE Annual Meeting 2026 deck: availability, outliers and night-to-night
spread of Rayleigh v1 (eprof_v1.1 = the v1.0 window without the Klett sign error, which does not
affect C_L), v2 (C8) and v2.2 (candidate N2.5) on the SAME corpus streams and nights (the
rayleigh_availability corpus: 22 CHM15k + 10 CL61 + 4 Mini-MPL, 2025-01 .. 2026-07).

Per stream and method:
  availability = valid / clear nights (indicators.availability)
  sigma_SD     = 1.4826*median|diff C|/sqrt(2)/median C   (indicators.sigma_sd, the C8 definition)
  out_kalman   = operational Kalman rolling-IQR screen     (indicators.kalman_outliers)
  out_drift3s  = report-09 drift-aware rate: |C - rolling median(9 nights)| > 3 sigma_rob
  out_local50  = phase-4 rate: |C / median(valid nights within +/-30 d) - 1| > 50 %, >= 5 neighbours
Prints per-type pooled availability and the median over streams of the other indicators, as JSON.

Run from the repo root:  python presentation/slide10_outliers_spread.py [cand_cfg=N2.5]
"""
import json
import sys
from datetime import datetime
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / "rayleigh_availability"))
import indicators as IND  # noqa: E402

DATA = Path("C:/DATA/Projects/202606_E-PROFILE_calibration/rayleigh_availability")
BASE, CAND = DATA / "baselines", DATA / "candidates"
MANIFEST = json.loads((REPO / "rayleigh_availability" / "scope_availability.json").read_text())


def load(p):
    return json.loads(p.read_text()) if p.exists() else None


def drift3s(rec):
    d, c = IND.constants(rec)
    if c.size < 12:
        return float("nan")
    med = np.array([np.median(c[max(0, k - 4):k + 5]) for k in range(c.size)])
    r = c - med
    s = 1.4826 * np.median(np.abs(r - np.median(r)))
    return 100.0 * float(np.mean(np.abs(r - np.median(r)) > 3 * s)) if s > 0 else float("nan")


def local50(rec):
    d, c = IND.constants(rec)
    if c.size < 20:
        return float("nan")
    t = [datetime.strptime(x, "%Y%m%d") for x in d]
    n = bad = 0
    for k in range(c.size):
        near = [c[j] for j in range(c.size) if j != k and abs((t[j] - t[k]).days) <= 30]
        if len(near) < 5:
            continue
        n += 1
        bad += abs(c[k] / np.median(near) - 1) > 0.5
    return 100.0 * bad / n if n else float("nan")


def main():
    cfg = sys.argv[1] if len(sys.argv) > 1 else "N2.5"
    out = {}
    for grp in ("CHM15k", "CL61", "Mini-MPL"):
        ins = [i for i in MANIFEST if i["group"] == grp]
        recs = {}
        for i in ins:
            r = dict(v1=load(BASE / f"base_eprof_v1.1_{i['label']}.json"),
                     v2=load(BASE / f"base_eprof_v2_{i['label']}.json"),
                     v22=load(CAND / f"cand_{cfg}_{i['label']}.json"))
            if all(r.values()):
                recs[i["label"]] = r
            else:
                print(f"skip {i['label']}: missing {[k for k, v in r.items() if not v]}", flush=True)
        res = {"n_streams": len(recs)}
        for m in ("v1", "v2", "v22"):
            av = [IND.availability(r[m]) for r in recs.values()]
            nc, nv = sum(a["n_clear"] for a in av), sum(a["n_valid"] for a in av)
            ko = [IND.kalman_outliers(r[m]) for r in recs.values()]
            sd = [IND.sigma_sd(r[m]) for r in recs.values()]
            d3 = [drift3s(r[m]) for r in recs.values()]
            l5 = [local50(r[m]) for r in recs.values()]
            res[m] = dict(
                n_clear=nc, n_valid=nv, avail_pooled=round(100.0 * nv / nc, 1) if nc else None,
                sigma_sd_median=round(float(np.nanmedian(sd)), 1),
                sigma_sd_iqr=[round(float(x), 1) for x in np.nanpercentile(sd, [25, 75])],
                out_kalman_pooled=round(100.0 * sum(k["n_outliers"] for k in ko)
                                        / max(1, sum(k["n_valid"] for k in ko)), 2),
                out_kalman_median=round(float(np.nanmedian([k["rate_per_100"] for k in ko])), 2),
                out_drift3s_median=round(float(np.nanmedian(d3)), 2),
                out_local50_median=round(float(np.nanmedian(l5)), 2),
                n_out_local50=None)
        out[grp] = res
        print(grp, json.dumps(res), flush=True)
    (REPO / "presentation" / "figs_annual2026").mkdir(parents=True, exist_ok=True)
    (REPO / "presentation" / "figs_annual2026" / "slide10_outliers_spread.json").write_text(
        json.dumps(out, indent=1))
    print("STATS", json.dumps(out))


if __name__ == "__main__":
    main()
