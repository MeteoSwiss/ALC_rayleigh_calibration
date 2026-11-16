# -*- coding: utf-8 -*-
"""q16e check 5: is the published (calibration_version 220 = eprof_v2.2) series seasonally biased by
its noise-aware FALLBACK tier?

eprof_v2.2 = eprof_v2 strict gates, plus a relaxed noise-relative tier that is used ONLY when v2
rejects the night (molecular_methods._select_optimal, `fallback`). Re-run sampled published CHM15k
nights with eprof_v2: a night that v2 rejects but v2.2 calibrates is a tier-2 night. Measure, per
month, the tier-2 share and the offset of tier-2 constants from the unit's tier-1 constants, and the
per-unit annual amplitude of ln C with and without tier-2 nights. Package untouched; outputs only in
C:/DATA/Projects/202610_CHM15k_stability/q16/.
"""
from __future__ import annotations

import csv
import logging
import os
import sys
import warnings
from concurrent.futures import ProcessPoolExecutor
from datetime import datetime, timedelta
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import L1_ROOT, OUT, chm15k_streams, load_constants  # noqa: E402
from q16e_klett_rerun import annual_fit  # noqa: E402

Q16 = OUT / "q16"
EPOCH = datetime(1970, 1, 1)
N_UNITS = 16
NIGHTS_PER_MONTH = 4
PAYERNE = "0-20000-0-06610_A"


def run(args):
    warnings.filterwarnings("ignore")
    logging.disable(logging.CRITICAL)
    from calibration import CalibrationOptions, InstrumentInfo, calibrate_rayleigh
    from calibration.config import DataLevel, InstrumentType
    s, ds, c_pub, zb_pub = args
    info = InstrumentInfo(site_name=s["site"], wmo_id=s["wmo"], identifier=s["ident"],
                          instrument_type=InstrumentType.CHM15k, latitude=s["lat"], longitude=s["lon"],
                          altitude=s.get("alt", 0.0))
    out = dict(stream=f"{s['wmo']}_{s['ident']}", date=ds, c_pub=c_pub, zb_pub=zb_pub)
    for meth in ("eprof_v2", "eprof_v2.2"):
        o = CalibrationOptions.from_json(REPO / "options.json")
        o.folder_root = L1_ROOT
        o.data_level = DataLevel.L1
        o.folder_output = Q16 / "q16e_recal" / f"w{os.getpid()}"   # one writer per process
        o.plot_main = False
        o.plot_all = False
        o.molecular_method = meth
        try:
            r = calibrate_rayleigh(ds, info, o)
        except Exception as exc:  # noqa: BLE001
            out[f"err_{meth}"] = str(exc)[:60]
            break
        out[f"c_{meth}"] = r.lidar_constant
        out[f"flag_{meth}"] = r.flag
        out[f"zb_{meth}"] = r.calibration_bottom_height
        if r.lidar_constant > 0:                 # v2 succeeded -> v2.2 is identical by construction
            break
    return out


def jobs():
    import json
    census = {f"{s['wmo']}_{s['ident']}": s for s in json.load(open(REPO / "validation/scope_l1_2026_census.json",
                                                                  encoding="utf-8"))}
    cands = []
    for wmo, ident, site, lat, lon in chm15k_streams():
        key = f"{wmo}_{ident}"
        if key == PAYERNE:
            continue
        c = load_constants(wmo, ident)
        if c is None:
            continue
        ok = (c["c"] > 0) & (c["t"] >= (datetime(2025, 1, 1) - EPOCH).days) & (c["t"] < (datetime(2026, 1, 1) - EPOCH).days)
        mon = np.array([(EPOCH + timedelta(days=float(x))).month for x in c["t"][ok]])
        if ok.sum() >= 150 and len(np.unique(mon)) == 12:
            cands.append((lon, key, c, ok))
    cands.sort()
    idx = np.linspace(0, len(cands) - 1, N_UNITS).round().astype(int)
    rng = np.random.default_rng(1616)
    out = []
    for k in idx:
        lon, key, c, ok = cands[k]
        t, t0, cc, zb = c["t"][ok], c["t0"][ok], c["c"][ok], c["zb"][ok]
        mon = np.array([(EPOCH + timedelta(days=float(x))).month for x in t])
        for m in range(1, 13):
            sel = np.where(mon == m)[0]
            for j in rng.choice(sel, size=min(NIGHTS_PER_MONTH, sel.size), replace=False):
                d = EPOCH + timedelta(days=float(np.floor(t0[j] + 0.5)))
                out.append((census[key], d.strftime("%Y%m%d"), float(cc[j]), float(zb[j])))
    return out


def main():
    J = jobs()
    print(f"{len(J)} nights")
    with ProcessPoolExecutor(max_workers=4) as ex:
        rows = list(ex.map(run, J, chunksize=4))
    keys = []
    for r in rows:
        keys += [k for k in r if k not in keys]
    with open(Q16 / "q16e_tier2_selection.csv", "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=keys)
        w.writeheader()
        w.writerows(rows)
    analyse(rows)


def analyse(rows):
    def f(r, k):
        try:
            return float(r.get(k, "nan"))
        except (TypeError, ValueError):
            return np.nan
    R = []
    for r in rows:
        c1, c2 = f(r, "c_eprof_v2"), f(r, "c_eprof_v2.2")
        tier = 1 if c1 > 0 else (2 if c2 > 0 else 0)
        rep = (c1 if tier == 1 else c2) / f(r, "c_pub") if tier else np.nan
        R.append((r["stream"], datetime.strptime(str(r["date"]), "%Y%m%d"), f(r, "c_pub"), tier, rep))
    tier = np.array([x[3] for x in R])
    mon = np.array([x[1].month for x in R])
    rep = np.array([x[4] for x in R])
    print(f"tier1 {np.sum(tier == 1)}, tier2 {np.sum(tier == 2)}, not reproduced {np.sum(tier == 0)}; "
          f"re-run/published within 2 %: {np.mean(np.abs(rep[tier > 0] - 1) < 0.02) * 100:.0f} %")
    print("month: tier-2 share of the published nights [%]")
    print("  " + " ".join(f"{m}:{100 * np.mean(tier[mon == m] == 2):.0f}" for m in range(1, 13)))
    # offset of tier-2 constants vs the unit's tier-1 constants within +-45 days
    offs, offm = [], []
    for i, x in enumerate(R):
        if x[3] != 2:
            continue
        nb = [y[2] for y in R if y[0] == x[0] and y[3] == 1 and abs((y[1] - x[1]).days) <= 45]
        if len(nb) >= 3:
            offs.append(np.log(x[2] / np.median(nb)))
            offm.append(x[1].month)
    offs = np.array(offs)
    print(f"tier-2 offset vs neighbouring tier-1 nights: n {offs.size}, median {100 * np.median(offs):+.1f} %, "
          f"mean {100 * np.mean(offs):+.1f} %, IQR {100 * np.percentile(offs, 25):+.1f}..{100 * np.percentile(offs, 75):+.1f} %")
    units = sorted({x[0] for x in R})
    a_all, a_t1, p_all, p_t1 = [], [], [], []
    for u in units:
        s = [x for x in R if x[0] == u and x[3] > 0]
        doy = np.array([x[1].timetuple().tm_yday for x in s], float)
        y = np.log([x[2] for x in s])
        t1 = np.array([x[3] == 1 for x in s])
        a, p = annual_fit(doy, y)
        b, q = annual_fit(doy[t1], y[t1])
        a_all.append(a)
        a_t1.append(b)
        p_all.append(p)
        p_t1.append(q)
    a_all, a_t1 = np.array(a_all), np.array(a_t1)
    print(f"per-unit annual amplitude of ln C (published nights): all {100 * np.nanmedian(a_all):.1f} % "
          f"(peak {np.nanmedian(p_all):.0f}), tier-1 only {100 * np.nanmedian(a_t1):.1f} % (peak {np.nanmedian(p_t1):.0f}); "
          f"median per-unit difference {100 * np.nanmedian(a_all - a_t1):+.2f} %-pt")


if __name__ == "__main__":
    main()
