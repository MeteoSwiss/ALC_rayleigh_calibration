"""fetch_dashboard_calib.py — take the calibration constants for the paper validation straight from the
LIVE operational dashboard (https://alc-calib.ch-meteoswiss-emermet.f.ewcloud.host/), so the validation
uses exactly the numbers the network publishes and nothing else (no local reruns, no dark variants).

Each station page embeds its time-series chart as Plotly JSON; the traces named
"Kalman estimate (v2.2)" (Rayleigh) and "Kalman estimate (O'Connor)" (liquid cloud) carry the daily
Kalman lidar constant C_L (Wiegner & Geiss 2012 convention: beta_att = rcs_0 / C_L) that the
dashboard displays. They are written in the paper-validation calibration format
    time, C_daily, C_daily_std, C_kalman, C_kalman_std
to <out>/<WMO>_<ident>_<rayleigh|cloud>.csv, where intercompare.load_calib_series reads them once
intercompare.CALIB points at <out>. The per-night values are taken from the yearly
ALC_calibration_<key><YYYY>.nc files the dashboard links (EWC bucket) and fill C_daily.

Usage:  python -m validation.paper.fetch_dashboard_calib <out_dir> <WMO_ident> [<WMO_ident> ...]
"""
from __future__ import annotations
import base64
import json
import re
import sys
import urllib.request
from datetime import datetime, timedelta, timezone
from pathlib import Path

import numpy as np

DASHBOARD = "https://alc-calib.ch-meteoswiss-emermet.f.ewcloud.host"
BUCKET = "https://object-store.os-api.cci2.ecmwf.int/eprofile-alc-dashboard/nc"
METHOD_CODE = {"rayleigh": 0, "cloud": 1}     # calibration_method flag values in the NetCDF


def _get(url, timeout=120):
    with urllib.request.urlopen(url, timeout=timeout) as r:
        return r.read()


def plotly_traces(html):
    """All trace dicts of every Plotly.newPlot(...) call on a page."""
    out = []
    for m in re.finditer(r'Plotly\.newPlot\(\s*"[^"]+",\s*(\[.*?\])\s*,\s*\{', html, re.S):
        try:
            out.extend(json.loads(m.group(1)))
        except json.JSONDecodeError:
            continue
    return out


def _values(y):
    """Plotly y array, plain list or binary-encoded {'dtype','bdata'}."""
    if isinstance(y, dict):
        return np.frombuffer(base64.b64decode(y["bdata"]), dtype=y["dtype"]).astype("f8")
    return np.array([np.nan if v is None else v for v in y], "f8")


def kalman_series(stream):
    """{method: (dates[str], C_L[float])} for one stream, from its dashboard station page."""
    html = _get(f"{DASHBOARD}/stations/{stream}.html").decode("utf-8")
    out = {}
    for t in plotly_traces(html):
        name = str(t.get("name", ""))
        if not name.startswith("Kalman estimate"):
            continue
        method = "cloud" if "O'Connor" in name else "rayleigh"
        out[method] = ([str(x)[:10] for x in t["x"]], _values(t["y"]))
    return out


def nightly_values(stream, years, cache_dir):
    """{method: {date: (C_L, unc)}} from the yearly calibration NetCDFs the dashboard links."""
    from netCDF4 import Dataset
    out = {"rayleigh": {}, "cloud": {}}
    for y in years:
        fn = cache_dir / f"ALC_calibration_{stream}{y}.nc"
        if not fn.exists():
            try:
                fn.write_bytes(_get(f"{BUCKET}/{stream}/ALC_calibration_{stream}{y}.nc"))
            except Exception:
                continue
        with Dataset(fn) as nc:
            t = np.asarray(nc["time"][:], "f8")
            c = np.ma.filled(nc["lidar_constant"][:].astype("f8"), np.nan)
            u = np.ma.filled(nc["lidar_constant_uncertainty"][:].astype("f8"), np.nan)
            m = np.ma.filled(nc["calibration_method"][:], -1)
        for ti, ci, ui, mi in zip(t, c, u, m):
            meth = {v: k for k, v in METHOD_CODE.items()}.get(int(mi))
            if meth is None or not np.isfinite(ci) or ci <= 0:
                continue
            d = (datetime(1970, 1, 1) + timedelta(days=float(ti))).strftime("%Y-%m-%d")
            out[meth][d] = (ci, ui)
    return out


def write_stream(stream, out_dir, years=(2025, 2026)):
    out_dir = Path(out_dir)
    nc_dir = out_dir / "nc"
    nc_dir.mkdir(parents=True, exist_ok=True)
    kal = kalman_series(stream)
    night = nightly_values(stream, years, nc_dir)
    written = []
    for method, (dates, vals) in kal.items():
        daily = night.get(method, {})
        alld = sorted(set(dates) | set(daily))
        kd = dict(zip(dates, vals))
        lines = ["time,C_daily,C_daily_std,C_kalman,C_kalman_std"]
        for d in alld:
            cd, cu = daily.get(d, (np.nan, np.nan))
            ck = kd.get(d, np.nan)
            f = lambda v: "" if not np.isfinite(v) else f"{v:.6e}"
            lines.append(f"{d},{f(cd)},{f(cu)},{f(ck)},")
        (out_dir / f"{stream}_{method}.csv").write_text("\n".join(lines) + "\n", encoding="utf-8")
        written.append((method, int(np.isfinite(vals).sum()), len(daily)))
    return written


def main(argv):
    if len(argv) < 2:
        print(__doc__)
        return 1
    out_dir = Path(argv[0])
    out_dir.mkdir(parents=True, exist_ok=True)
    for stream in argv[1:]:
        for method, nk, nd in write_stream(stream, out_dir):
            print(f"  {stream:24s} {method:8s} kalman={nk:4d} days  nightly={nd:4d}")
    (out_dir / "SOURCE.txt").write_text(
        f"Kalman C_L from {DASHBOARD} station pages, fetched {datetime.now(timezone.utc):%Y-%m-%d %H:%M} UTC\n",
        encoding="utf-8")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
