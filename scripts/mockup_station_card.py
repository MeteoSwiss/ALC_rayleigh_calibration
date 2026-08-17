#!/usr/bin/env python3
"""Standalone mockup of the redesigned station availability card, on REAL data.

Deliberately a thin shell over the production functions -- it contains argument parsing, an
HTML wrapper and the ad-hoc L1 cloud-cover probe, and nothing else. Every pixel it shows is
drawn by the code that will ship (monitoring.charts / monitoring.config / the real style.css),
so the mockup cannot flatter the design.

What it shows, per station:
  * the 3-row card (instrument status / mean cloud cover / calibration), 4 rows for a CL61
    which carries both a Rayleigh and a liquid-cloud calibration;
  * the re-encoded instrument-monitoring figure (x0/dx + float32);
  * the compact header line.

Row 2 needs `mean_cloud_cover`, which the operational producer does not write yet: pass
--l1-root and --l1-month to compute it live from the L1 `cloud_amount` octas, using the exact
0..8 screen the producer will use. Without it the row is simply absent (the card degrades the
way an old archive would), and the page says so.

Examples
--------
  python scripts/mockup_station_card.py ^
      --cal-dir  C:/DATA/Projects/202606_E-PROFILE_calibration/calout_v22_04 ^
      --status-dir C:/DATA/Projects/202606_E-PROFILE_calibration/diag_v22_04 ^
      --key 0-20000-0-06610_A --key 0-20000-0-10393_C ^
      --l1-root D:/E-PROFILE_L1_2026 --l1-month 202607 ^
      --out C:/tmp/mock_station_card.html
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import plotly.offline as pyo

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from monitoring import charts, config  # noqa: E402
from monitoring.render import _load_hk, _load_status  # noqa: E402

#: L1 cloud_amount / tcc are octas 0..8. The observed missing-value sentinels are -9, -99 and
#: -32767 (so the housekeeping screen `> -990` would let -9/-99 through as negative cloud cover),
#: and the value 9 appears on CHM15k with no documented meaning -- presumably "obscured / not
#: determinable". Both are treated as missing until confirmed; averaging 9 in would make every
#: foggy day read as more-than-overcast.
OCTA_MIN, OCTA_MAX = 0.0, 8.0


def _wmo_ident(key: str):
    """'0-20000-0-06610_A' -> ('0-20000-0-06610', 'A')."""
    wmo, _, ident = key.rpartition("_")
    return wmo, ident


def _month_with_data(l1_root: Path, key: str, preferred: str) -> str | None:
    """The preferred YYYYMM if that stream has L1 files there, else the most recent month that does.

    A stream's identifier can change over its life (the Lindenberg CL61 files are '_C' until
    mid-2026 and '_0' after), so a fixed month is not a safe assumption for every station."""
    wmo, ident = _wmo_ident(key)
    root = Path(l1_root) / wmo
    if any((root / preferred[:4] / preferred[4:6]).glob(f"L1_{wmo}_{ident}*.nc")):
        return preferred
    months = sorted({f"{p.parent.name}{p.name}"
                     for p in root.glob("*/*")
                     if p.is_dir() and any(p.glob(f"L1_{wmo}_{ident}*.nc"))})
    return months[-1] if months else None


def probe_cloud_cover(l1_root: Path, key: str, month: str) -> pd.DataFrame:
    """Daily mean cloud cover (octas) for one month, straight off the L1 archive.

    This is the ad-hoc stand-in for the producer change (a new `mean_cloud_cover` column in
    <key>_status.csv computed inside the existing single-open monitoring pass). Returns an empty
    frame when the month or the variable is absent, so the caller can say so honestly."""
    import netCDF4

    wmo, ident = _wmo_ident(key)
    mdir = Path(l1_root) / wmo / month[:4] / month[4:6]
    rows = []
    for f in sorted(mdir.glob(f"L1_{wmo}_{ident}*.nc")):
        date = f.stem[-8:]
        try:
            with netCDF4.Dataset(f) as ds:
                rec = _day_cloud_cover(ds)
        except Exception as exc:                       # a single unreadable file must not abort
            print(f"    ! {f.name}: {exc}", file=sys.stderr)
            continue
        if rec is not None:
            rows.append({"date": date, **rec})
    return pd.DataFrame(rows)


def _day_cloud_cover(ds) -> dict | None:
    """One open L1 file -> {mean_cloud_cover, cloud_cover_n, cloud_src}, or None.

    Prefers the reported octas and falls back to cloud-base presence. The fallback is NOT cosmetic:
    verified on the real archive, the Lindenberg CL61 declares `cloud_amount` and fills it entirely
    with the sentinel -9, so a producer that trusted the variable's presence (or reused the
    housekeeping `> -990` screen) would publish a permanent -9 octa row for that stream."""
    n_prof, amount = None, None
    if "cloud_amount" in ds.variables:
        a = np.ma.filled(np.asarray(ds.variables["cloud_amount"][:], "f8"), np.nan)
        if a.ndim > 1:              # (time, layer) on the CL61 -> the BASE layer is the octa value
            a = a[:, 0]
        n_prof = a.size
        amount = a[np.isfinite(a) & (a >= OCTA_MIN) & (a <= OCTA_MAX)]
    # Trust the octas only if the instrument actually reports them for a fair share of the day.
    if amount is not None and n_prof and amount.size >= 0.2 * n_prof:
        return {"mean_cloud_cover": float(amount.mean()), "cloud_cover_n": int(amount.size),
                "cloud_src": "cloud_amount"}
    if "cloud_base_height" in ds.variables:
        cbh = np.ma.filled(np.asarray(ds.variables["cloud_base_height"][:], "f8"), np.nan)
        if cbh.ndim > 1:
            cbh = cbh[:, 0]
        ok = np.isfinite(cbh) & (cbh > 0)
        if ok.size:
            # Fraction of profiles carrying a cloud base, scaled to octas. NOT the same physical
            # quantity as reported octas -- recorded in cloud_src so the hover stays honest.
            return {"mean_cloud_cover": float(8.0 * ok.mean()), "cloud_cover_n": int(ok.size),
                    "cloud_src": "cbh"}
    return None


def build(key: str, cal_dir: Path, status_dir: Path, l1_root: Path | None, l1_month: str | None):
    """One station -> (html fragment, note). Uses the production loaders and chart builders."""
    status_df = _load_status(status_dir, key)
    if status_df is None:
        return None, f"{key}: no {key}_status.csv under {status_dir}"
    hk_df = _load_hk(status_dir, key)

    cal_path = Path(cal_dir) / key / f"{key}_cal.csv"
    if not cal_path.exists():
        return None, f"{key}: no {cal_path.name} under {cal_dir}"
    cal = pd.read_csv(cal_path, dtype={"date": str})
    cal["key"] = key
    methods = [m for m in config.METHOD_ORDER if len(cal[cal["method"] == m])]

    note = ""
    if l1_root and l1_month:
        month = _month_with_data(Path(l1_root), key, l1_month) or l1_month
        cc = probe_cloud_cover(Path(l1_root), key, month)
        if len(cc):
            status_df = status_df.merge(cc, on="date", how="left")
            fallback = "" if month == l1_month else f" (no L1 for this stream in {l1_month})"
            note = (f"row 2 = REAL octas probed live from {len(cc)} L1 files of {month}{fallback}, "
                    f"source '{cc['cloud_src'].iloc[0]}'; the rest of the record is blank because "
                    f"the operational producer does not write mean_cloud_cover yet")
        else:
            note = f"no L1 cloud_amount found for {month} — row 2 omitted"
    else:
        note = "no --l1-root/--l1-month given — row 2 omitted (pass them to see real octas)"

    fig = charts.daily_availability_rows(status_df, cal, methods)
    if fig is None:
        return None, f"{key}: availability figure could not be built"

    meta_line = (f"{key} · {len(cal)} calibration rows · methods: "
                 f"{', '.join(config.method_label(m) for m in methods)} · "
                 f"{len(status_df)} status days")
    parts = [f"<h2>{key}</h2>", f'<p class="subtitle">{meta_line}</p>']
    parts.append('<div class="avail-legend"><b>Status</b>'
                 '<span class="ak ak-pass"></span> Pass'
                 '<span class="ak ak-warning"></span> Warning'
                 '<span class="ak ak-error"></span> Error'
                 '<span class="ak ak-nodata"></span> No data'
                 '<span class="legsep"></span><b>Cloud cover</b>'
                 '<span class="ak ak-cc0"></span> 0'
                 '<span class="ak ak-cc4"></span> 4'
                 '<span class="ak ak-cc8"></span> 8 octas'
                 '<span class="legsep"></span><b>Calibration</b>'
                 + "".join(f'<span class="ak" style="background:{config.CAL_CLASS_COLORS[c]}" '
                           f'title="{config.CAL_CLASS_LABELS[c]}"></span>'
                           f' {config.CAL_CLASS_SHORT[c]}' for c in config.CAL_CLASS_ORDER)
                 + "</div>")
    parts.append(f'<p class="muted">{note}</p>')
    parts.append(f'<div class="card">{charts.fig_to_div(fig, f"fig-avail-{key}")}</div>')
    if hk_df is not None:
        parts.append("<h3>Instrument monitoring "
                     '<span class="muted">— re-encoded x0/dx + float32</span></h3>')
        parts.append(f'<div class="card">'
                     f'{charts.fig_to_div(charts.monitoring_timeseries(hk_df), f"fig-hk-{key}")}'
                     f'</div>')
    return "\n".join(parts), note


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--cal-dir", type=Path, required=True,
                    help="tree holding <key>/<key>_cal.csv (e.g. the v2.2 network run)")
    ap.add_argument("--status-dir", type=Path, default=None,
                    help="tree holding <key>/<key>_status.csv and _hk.csv (default: --cal-dir)")
    ap.add_argument("--key", action="append", required=True, help="station key (repeatable)")
    ap.add_argument("--l1-root", type=Path, default=None, help="L1 archive root for the cloud probe")
    ap.add_argument("--l1-month", default=None, help="YYYYMM to probe for cloud cover")
    ap.add_argument("--out", type=Path, required=True, help="output HTML file")
    args = ap.parse_args()

    status_dir = args.status_dir or args.cal_dir
    css = (Path(__file__).resolve().parent.parent / "monitoring/static/style.css").read_text(
        encoding="utf-8")

    blocks, notes = [], []
    for key in args.key:
        print(f"  {key} ...", flush=True)
        html, note = build(key, args.cal_dir, status_dir, args.l1_root, args.l1_month)
        notes.append(f"{key}: {note}")
        blocks.append(html if html else f'<p class="muted">{note}</p>')

    page = f"""<!doctype html><html lang="en"><head><meta charset="utf-8">
<title>Mockup — station availability card</title>
<style>{css}
body {{ padding: 22px 28px; }} h3 {{ margin: 18px 0 6px; }}
.mockhdr {{ background:#fff3cd; border:1px solid #ffe08a; border-radius:8px;
  padding:10px 14px; margin-bottom:18px; }}</style>
<script>{pyo.get_plotlyjs()}</script></head><body>
<div class="mockhdr"><b>MOCKUP</b> — drawn by the production functions
(<code>monitoring.charts.daily_availability_rows</code>,
<code>monitoring.charts.monitoring_timeseries</code>, the real <code>style.css</code>).
Plotly and the CSS are inlined here so the file opens by double-click; in production they are
shared assets.</div>
{"<hr>".join(blocks)}
</body></html>"""
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(page, encoding="utf-8")
    size = args.out.stat().st_size / 1e6
    print(f"-> {args.out}  ({size:.2f} MB)")
    for n in notes:
        print(f"   {n}")


if __name__ == "__main__":
    main()
