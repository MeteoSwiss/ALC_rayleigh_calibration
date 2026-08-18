"""Extract per-day calibration-scene CBH for every successful cloud calibration.

Pour chaque (flux, jour) avec calibration nuage flag 1.0/0.5 dans calout_v22_04,
lit le L1 local (D:/E-PROFILE_L1_2026) et calcule la CBH mediane des profils dont
la CBH (couche la plus basse) tombe dans la fenetre de calibration [500, 2400] m
(calibration/cloud/calibration.py: cbh_minheight=500, cbh_maxheight=2400).
"""
import csv
import json
import os
import sys
from concurrent.futures import ProcessPoolExecutor

import numpy as np

CAL_ROOT = r"C:/DATA/Projects/202606_E-PROFILE_calibration/calout_v22_04"
L1_ROOT = r"D:/E-PROFILE_L1_2026"
CENSUS = r"C:/Users/hervo/OneDrive/Documents/ALC_rayleigh_calibration/validation/scope_l1_2026_census.json"
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "cbh_scenes.csv")

CBH_MIN, CBH_MAX = 500.0, 2400.0


def list_jobs():
    """One job per flux: (key, wmo, ident, list of (date, cal, unc, nprof, flag))."""
    census = json.load(open(CENSUS))
    types = {(e["wmo"], e["ident"]): e["type"] for e in census}
    jobs = []
    for d in sorted(os.listdir(CAL_ROOT)):
        f = os.path.join(CAL_ROOT, d, d + "_cal.csv")
        if not os.path.exists(f):
            continue
        wmo, ident = d.rsplit("_", 1)
        itype = types.get((wmo, ident), "?")
        rows = []
        with open(f) as fh:
            for r in csv.DictReader(fh):
                if r["method"] == "cloud" and r["flag"] in ("1.0", "1", "0.5"):
                    rows.append((r["date"], r["cal_value"], r["uncertainty"],
                                 r["n_profiles"], r["flag"]))
        if rows:
            jobs.append((d, wmo, ident, itype, rows))
    return jobs


def process_flux(job):
    key, wmo, ident, itype, rows = job
    import netCDF4 as nc  # import in worker
    out = []
    for date, cal, unc, nprof, flag in rows:
        yyyy, mm = date[:4], date[4:6]
        path = os.path.join(L1_ROOT, wmo, yyyy, mm, f"L1_{wmo}_{ident}{date}.nc")
        rec = dict(key=key, type=itype, date=date, cal_value=cal, uncertainty=unc,
                   n_profiles=nprof, flag=flag, cbh_med="", cbh_p25="", cbh_p75="",
                   n_win="", frac_below800="", l1="0")
        if os.path.exists(path):
            try:
                with nc.Dataset(path) as ds:
                    v = ds.variables.get("cloud_base_height")
                    if v is not None:
                        # cast AVANT filled: Vaisala stocke cbh en int32
                        a = np.ma.asarray(v[:]).astype("float64").filled(np.nan)
                        if a.ndim == 2:
                            a = a[:, 0]  # lowest layer, same as calibration code
                        a[(a < 0) | (a > 20000)] = np.nan
                        w = a[(a >= CBH_MIN) & (a <= CBH_MAX)]
                        rec["l1"] = "1"
                        if w.size:
                            rec["cbh_med"] = f"{np.median(w):.1f}"
                            rec["cbh_p25"] = f"{np.percentile(w, 25):.1f}"
                            rec["cbh_p75"] = f"{np.percentile(w, 75):.1f}"
                            rec["n_win"] = str(w.size)
                            rec["frac_below800"] = f"{np.mean(w < 800.0):.3f}"
            except Exception as e:  # corrupt file etc.
                rec["l1"] = "err:" + type(e).__name__
        out.append(rec)
    return out


def main():
    jobs = list_jobs()
    print(f"{len(jobs)} flux, {sum(len(j[4]) for j in jobs)} cloud-success days")
    cols = ["key", "type", "date", "cal_value", "uncertainty", "n_profiles",
            "flag", "cbh_med", "cbh_p25", "cbh_p75", "n_win", "frac_below800", "l1"]
    n_done = 0
    with open(OUT, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols)
        w.writeheader()
        with ProcessPoolExecutor(max_workers=8) as ex:
            for recs in ex.map(process_flux, jobs, chunksize=1):
                for r in recs:
                    w.writerow(r)
                n_done += 1
                if n_done % 25 == 0:
                    print(f"  {n_done}/{len(jobs)} flux", flush=True)
    print("wrote", OUT)


if __name__ == "__main__":
    main()
