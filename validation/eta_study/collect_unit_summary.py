"""Per-unit summary of the eta-study array (CSCS) or of a local run: one row per unit with the
tail-slope alpha_k and the full-shape alpha (nominal pulse), plus the metadata needed to tell
"cloud" from "receiver" (type, gate spacing, firmware/serial from the L1 global attributes,
lat/lon/alt, country from the L1 `site_location`), and a long table of monthly alpha_k medians.

Reads  $ETA_OUT/units/*.json (or units.json), $ETA_OUT/<key>.npz (gate spacing),
       $ETA_BASE/step2/profiles/<key>_working.(npz|parquet), and one L1 file per unit-month.
Writes $ETA_BASE/summary/unit_summary.csv, unit_month_alpha_k.csv (small; the only files to copy back).

Run:  python validation/eta_study/collect_unit_summary.py
"""
from __future__ import annotations

import json
import os
import sys
from collections import Counter
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import fit_alpha_shape as F  # noqa: E402
import extract_incloud_profiles as E  # noqa: E402


def load_profiles(key):
    import pandas as pd
    p = F.OUT / "profiles" / f"{key}_working.parquet"
    if p.exists():
        return pd.read_parquet(p)
    p = F.OUT / "profiles" / f"{key}_working.npz"
    if p.exists():
        return pd.DataFrame(dict(np.load(p, allow_pickle=True)))
    return None


def l1_meta(key, months):
    """Firmware / serial versions (with the months they cover) and the country, from the L1."""
    from netCDF4 import Dataset
    wmo, ident = key.rsplit("_", 1)
    fw, sn, country = Counter(), Counter(), ""
    for ym in months:
        d = E.L1_ROOT / wmo / ym[:4] / ym[4:]
        files = sorted(d.glob(f"L1_{wmo}_{ident}{ym}*.nc")) if d.is_dir() else []
        if not files:
            continue
        try:
            with Dataset(files[0]) as nc:
                fw[str(getattr(nc, "instrument_firmware_version", "?"))] += 1
                sn[str(getattr(nc, "instrument_serial_number", "?"))] += 1
                loc = str(getattr(nc, "site_location", ""))
                country = loc.split(",")[-1].strip() if "," in loc else country
        except OSError:
            continue
    return ";".join(fw), ";".join(sn), country


def q(x, p):
    x = np.asarray(x, float)
    x = x[np.isfinite(x)]
    return float(np.percentile(x, p)) if x.size else np.nan


def main():
    import pandas as pd
    src = Path(os.environ.get("ETA_OUT", str(F.SRC)))
    if (src / "units.json").exists():
        units = json.load(open(src / "units.json", encoding="utf-8"))
    else:
        units = {p.stem: json.load(open(p, encoding="utf-8")) for p in sorted((src / "units").glob("*.json"))}
    out = F.BASE / "summary"
    out.mkdir(parents=True, exist_ok=True)
    rows, mrows = [], []
    for key, u in sorted(units.items()):
        df = load_profiles(key)
        npz = src / f"{key}.npz"
        if df is None or not npz.exists():
            rows.append(dict(key=key, itype=u["itype"], site=u["site"], status="no output"))
            continue
        z = np.load(npz, allow_pickle=True)
        v = df[df.valid.astype(bool)].copy()
        p = int(F.NOMINAL[u["itype"]])
        full = v.loc[v[f"edge_p{p}_ag3.3"] == 0, f"alpha_p{p}_ag3.3"]
        v["month"] = (v.date.astype(int) // 100).astype(str)
        fw, sn, country = l1_meta(key, sorted(v.month.unique()))
        rows.append(dict(key=key, itype=u["itype"], site=u["site"], country=country, lat=u["lat"], lon=u["lon"],
                         alt=u["alt"], distorted=u["distorted"], gate_m=float(z["dr_m"]),
                         coarsen=";".join(f"{a}x{b}" for a, b in np.unique(z["coarsen_factors"], axis=0)),
                         firmware=fw, serial=sn, n=len(v), n_days=v.date.nunique(),
                         alpha_k_med=q(v.alpha_k, 50), alpha_k_p10=q(v.alpha_k, 10), alpha_k_p90=q(v.alpha_k, 90),
                         frac_alpha_k_outside=float(v.alpha_k.isna().mean()),
                         alpha_full_med=q(full, 50), alpha_full_p10=q(full, 10), alpha_full_p90=q(full, 90),
                         kapp_med=q(v.k_app, 50), half_height_med=q(v.half_height_m, 50),
                         F50_med=q(v.F50, 50), cbh_med=q(v.cbh_m, 50), status="ok"))
        for m, g in v.groupby("month"):
            mrows.append(dict(key=key, itype=u["itype"], site=u["site"], month=m, n=len(g),
                              n_days=g.date.nunique(), alpha_k_med=q(g.alpha_k, 50),
                              alpha_full_med=q(g[f"alpha_p{p}_ag3.3"], 50), cbh_med=q(g.cbh_m, 50)))
    pd.DataFrame(rows).to_csv(out / "unit_summary.csv", index=False, float_format="%.4g")
    pd.DataFrame(mrows).to_csv(out / "unit_month_alpha_k.csv", index=False, float_format="%.4g")
    print(f"{len(rows)} units -> {out / 'unit_summary.csv'}", flush=True)


if __name__ == "__main__":
    main()
