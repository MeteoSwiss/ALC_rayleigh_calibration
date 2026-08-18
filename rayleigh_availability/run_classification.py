# -*- coding: utf-8 -*-
"""Phase 3 — produce Cloudnet target classifications for the availability corpus.

The pre-fit cell mask needs a classification curtain per stream-day; this is the (expensive) step
that makes one. It is deliberately a SEPARATE run from the calibration sweep so the classification
is computed ONCE and every candidate configuration reads the same product.

Only a subset of the corpus is classified -- the streams that can answer the Phase-3 questions:
reference sites (a co-located CL61), the co-located Amsterdam quad, the healthy streams that show
whether the mask costs clean nights, and the worst-rejection streams that show whether it converts
steady-layer failures. Everything else would only add cost.

Run:  python rayleigh_availability/run_classification.py [--force] [label-substring ...]
Out:  <DATA>/rayleigh_availability/classification/<key>/classification/<wmo>/<year>/*.nc
"""
from __future__ import annotations
import json
import logging
import os
import sys
import warnings
from concurrent.futures import ProcessPoolExecutor, as_completed
from datetime import datetime, timedelta
from pathlib import Path

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[_v] = "1"
try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))

from calibration import CalibrationOptions                                      # noqa: E402
from calibration.io.cams import find_cams_file                                  # noqa: E402
from calibration.io.instrument_day import load_instrument_day                   # noqa: E402

MANIFEST = json.loads((REPO / "rayleigh_availability" / "scope_availability.json").read_text())
DATA = Path("C:/DATA/Projects/202606_E-PROFILE_calibration/rayleigh_availability")
OUT = DATA / "classification"
L1_ROOT = Path("D:/E-PROFILE_L1_2026")

# The classified subset, by role in the Phase-3 argument.
SUBSET = {
    "PAYERNE_CHM15k_A":            "reference (co-located CL61, dC_L/dz < 0)",
    "LINDENBERG_CHM15k_0":         "reference (co-located CL61, dC_L/dz > 0)",
    "AMSTERDAM_AP_SCHIP_CHM15k_A": "inter-unit quad",
    "AMSTERDAM_AP_SCHIP_CHM15k_B": "inter-unit quad",
    "AMSTERDAM_AP_SCHIP_CHM15k_C": "inter-unit quad",
    "AMSTERDAM_AP_SCHIP_CHM15k_D": "inter-unit quad",
    "HAMBURG_CHM15k_0":            "healthy (does the mask cost clean nights?)",
    "HOHENPEISSENBERG_CHM15k_0":   "healthy (does the mask cost clean nights?)",
    "GOTTFRIEDING_CHM15k_0":       "worst rejection (steady-layer conversion)",
    "GUADIANA_UGR_CHM15k_A":       "worst rejection (steady-layer conversion)",
    "MONTSEC_CHM15k_A":            "worst rejection (steady-layer conversion)",
    "KLIPPENECK_CHM15k_0":         "worst rejection (steady-layer conversion)",
}

# The co-located CL61s, which are NOT calibration streams here but the better AEROSOL DETECTOR: they
# carry depolarization, so their classification distinguishes layers a single-channel CHM15k cannot.
# Classifying them lets the CHM15k's recovered nights be screened by an INDEPENDENT, better-informed
# instrument -- the difference between "our mask found no aerosol" and "there is no aerosol".
EXTRA = [
    dict(label="PAYERNE_CL61_C", site="PAYERNE", wmo="0-20000-0-06610", ident="C", type="CL61",
         twin="PAYERNE_CHM15k_A"),
    dict(label="LINDENBERG_CL61_C", site="LINDENBERG", wmo="0-20000-0-10393", ident="C",
         type="CL61", twin="LINDENBERG_CHM15k_0"),
]


def extra_streams():
    """The EXTRA reference streams, taking position and period from their co-located twin."""
    out = []
    for e in EXTRA:
        twin = next((i for i in MANIFEST if i["label"] == e["twin"]), None)
        if twin is None:
            continue
        out.append({**{k: twin[k] for k in ("lat", "lon", "alt", "periods", "group", "split")},
                    **{k: e[k] for k in ("label", "site", "wmo", "ident", "type")}})
    return out


CHUNK_DAYS = 45          # work unit: one stream x ~45 days, so 22 streams do not idle 20 cores


def dates_of(inst, period=None):
    """Every date in the stream's periods, plus the day BEFORE each period start.

    A Rayleigh night for date d draws on the evening of d-1, so classifying only the corpus dates
    would leave the first night of each period half-covered.
    """
    seen = set()
    for first, last in ([period] if period else inst["periods"]):
        d = datetime.strptime(first, "%Y%m%d") - timedelta(days=1)
        d1 = datetime.strptime(last, "%Y%m%d")
        while d <= d1:
            ds = d.strftime("%Y%m%d")
            if ds not in seen:
                seen.add(ds)
                yield ds
            d += timedelta(days=1)


def chunks_of(inst):
    """The stream's periods split into CHUNK_DAYS work units."""
    out = []
    for first, last in inst["periods"]:
        d = datetime.strptime(first, "%Y%m%d")
        d1 = datetime.strptime(last, "%Y%m%d")
        while d <= d1:
            e = min(d + timedelta(days=CHUNK_DAYS - 1), d1)
            out.append((d.strftime("%Y%m%d"), e.strftime("%Y%m%d")))
            d = e + timedelta(days=1)
    return out


def l1_file(wmo, ident, ds):
    return L1_ROOT / wmo / ds[:4] / ds[4:6] / f"L1_{wmo}_{ident}{ds}.nc"


def run_one(job):
    inst, period = job
    warnings.filterwarnings("ignore")
    logging.getLogger().setLevel(logging.ERROR)
    try:
        from calibration.classify import cams_to_model, ceilo_from_shared
        from ceiloclass.classification import classify
        from ceiloclass.write import write_classification
    except ImportError as exc:
        return inst["label"], 0, 0, f"ceiloclass missing: {exc}"

    o = CalibrationOptions.from_json(REPO / "options.json")
    key = inst["label"]
    # Monthly first (what the calibration baselines used), then the daily cache, which is the only
    # source covering Jul-2026 -- the classification only needs the temperature profile.
    cams_dirs = [Path(str(o.cams_folder)), Path("D:/CAMS_daily")]
    force = "--force" in sys.argv
    n_ok = n_skip = 0
    err = ""
    for ds in dates_of(inst, period):
        fp = l1_file(inst["wmo"], inst["ident"], ds)
        if not fp.exists():
            continue
        cdir = OUT / key / "classification" / inst["wmo"] / ds[:4]
        out_nc = cdir / f"{key}_{ds}_classification.nc"
        if out_nc.exists() and not force:
            n_skip += 1
            continue
        cams = next((c for c in (find_cams_file(p, ds) for p in cams_dirs if p.is_dir())
                     if c is not None), None)
        if cams is None:
            continue                       # no temperature -> no melting layer -> no classification
        try:
            idd = load_instrument_day([str(fp)], inst["type"], str(cams_dirs[0]),
                                      read_cams=False, build_working=True)
            if idd is None:
                continue
            sub = idd.slice_to_date(datetime.strptime(ds, "%Y%m%d").date())
            ceilo = ceilo_from_shared(sub)
            alt = float(sub.altitude)
            model = cams_to_model(str(cams), inst["lat"], inst["lon"], ceilo.time, ceilo.range, alt)
            result = classify(ceilo, model, altitude=alt, use_wet_bulb=False)
            cdir.mkdir(parents=True, exist_ok=True)
            write_classification(result, out_nc, wavelength=float(ceilo.wavelength), altitude=alt,
                                 latitude=inst["lat"], longitude=inst["lon"],
                                 location=inst["site"], source_files=[str(fp)])
            n_ok += 1
        except Exception as exc:           # one bad day must not lose the stream
            err = f"{ds}: {type(exc).__name__}: {exc}"[:110]
    return key, n_ok, n_skip, err


def main():
    pats = [a.upper() for a in sys.argv[1:] if not a.startswith("--")]
    insts = [i for i in MANIFEST if i["label"] in SUBSET] + extra_streams()
    if pats:
        insts = [i for i in insts if any(p in i["label"].upper() for p in pats)]
    jobs = [(i, p) for i in insts for p in chunks_of(i)]
    workers = int(os.environ.get("RA_WORKERS", str(max(1, (os.cpu_count() or 8) - 2))))
    print(f"classifying {len(insts)} streams as {len(jobs)} chunks over {workers} workers -> {OUT}",
          flush=True)
    done = 0
    with ProcessPoolExecutor(max_workers=workers) as ex:
        futs = {ex.submit(run_one, j): j for j in jobs}
        for fut in as_completed(futs):
            done += 1
            try:
                key, n_ok, n_skip, err = fut.result()
                print(f"  [{done}/{len(jobs)}] {key:30s} {futs[fut][1][0]} {n_ok:3d} new, "
                      f"{n_skip:3d} cached{('  last error: ' + err) if err else ''}", flush=True)
            except Exception as e:
                print(f"  [{done}/{len(jobs)}] FAILED {futs[fut][0]['label']}: {e}", flush=True)
    print("CLASSIFICATION_DONE", flush=True)


if __name__ == "__main__":
    main()
