# -*- coding: utf-8 -*-
"""Time-height PCOLOR quicklook for EVERY scanned instrument (all 31), one panel per station on a
per-type grid, so the actual clear-day signal (and any visible banding) is auditable network-wide.
Height on the Y axis (convention); time on X; log attenuated backscatter. Panels are framed green
(correctable) / grey (not) and titled with the offset amplitude and split-half reproducibility."""
import sys, glob, csv, warnings
warnings.filterwarnings("ignore")
sys.path.insert(0, r"C:/Users/hervo/OneDrive/Documents/ALC_rayleigh_calibration")
try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import matplotlib.dates as mdates
from validation.paper._offset_lib import read_file, list_files

D = Path("C:/DATA/Projects/202606_E-PROFILE_calibration/figs_paper_validation/paper_python")
NPZ = D / "network_offset"
REPORT_FIG = Path("C:/Users/hervo/OneDrive/Documents/ALC_rayleigh_calibration/doc/reports/figs_paper_report")
ZTOP = 6000.0

meta = {}
for r in csv.DictReader(open(D / "network_offset_coeffs_final.csv", encoding="utf-8")):
    meta[(r["wmo"], r["ident"], r["itype"])] = r


def clear_day(wmo, ident):
    best = None
    for f in list_files(wmo, ident, ("05",)):
        try:
            tt, r, x, cbh, _ = read_file(f)
        except Exception:
            continue
        if x.ndim != 2 or tt.size < 200:
            continue
        nc = np.nanmax(cbh, axis=1) < 0 if cbh is not None else np.ones(tt.size, bool)
        if best is None or int(nc.sum()) > best[0]:
            best = (int(nc.sum()), tt, r, x)
    return best[1:] if best else None


def stations(itype):
    out = []
    for f in sorted(glob.glob(str(NPZ / f"*_{itype}.npz"))):
        name = Path(f).stem
        wmo = "_".join(name.split("_")[:-2]); ident = name.split("_")[-2]
        z = np.load(f, allow_pickle=True); m = meta.get((wmo, ident, itype), {})
        out.append((str(z["site"]), wmo, ident, m))
    out.sort(key=lambda t: (-(t[3].get("correctable") == "1"), -float(t[3].get("rel_amp_pct", 0) or 0)))
    return out


TYPES = tuple(sys.argv[1:]) or ("CL51", "CL31", "CL61")     # optional: process only listed types
for itype in TYPES:
    rows = stations(itype)
    n = len(rows); ncol = 3; nrow = int(np.ceil(n / ncol))
    fig, ax = plt.subplots(nrow, ncol, figsize=(16, 2.9 * nrow), squeeze=False)
    for k, (site, wmo, ident, m) in enumerate(rows):
        a = ax[k // ncol][k % ncol]
        got = clear_day(wmo, ident)
        corr = m.get("correctable", "0") == "1"
        col = "#2ca02c" if corr else "0.6"
        for s in a.spines.values():
            s.set_color(col); s.set_linewidth(2.4 if corr else 1.0)
        if got is None:
            a.text(0.5, 0.5, "no data", ha="center"); continue
        tt, rng, X = got
        beta = np.abs(X) * (1.0 if itype == "CL61" else 1e-2)     # nominal β_att scale
        zm = rng <= ZTOP; tnum = mdates.date2num(tt)
        pos = beta[:, zm][beta[:, zm] > 0]
        # focus the colour scale on the bulk signal (25-99.5 pct) so weak, uniform profiles (CL61)
        # are not washed to the dark floor by a few bright surface pixels
        vmin = np.nanpercentile(pos, 25) if pos.size else 1e-3
        vmax = np.nanpercentile(pos, 99.5) if pos.size else 1.0
        if not np.isfinite(vmin) or vmin <= 0:
            vmin = 1e-4
        if not np.isfinite(vmax) or vmax <= vmin:
            vmax = vmin * 100                                      # guard degenerate LogNorm limits
        a.imshow(np.clip(beta[:, zm].T, vmin, vmax), aspect="auto", origin="lower",
                 extent=[tnum[0], tnum[-1], 0, rng[zm].max() / 1000], cmap="turbo",
                 norm=LogNorm(vmin=vmin, vmax=vmax))
        a.xaxis_date(); a.xaxis.set_major_formatter(mdates.DateFormatter("%H")); a.tick_params(labelsize=7)
        if k % ncol == 0:
            a.set_ylabel("height [km]", fontsize=8)
        amp = float(m.get("rel_amp_pct", 0) or 0); rep = float(m.get("split_half_repro", 0) or 0)
        a.set_title(f"{site[:20]}  amp={amp:.1f}% repro={rep:+.2f}" + ("  ✔" if corr else ""),
                    fontsize=8.2, color=("#1a7a1a" if corr else "0.25"))
        print(f"  {itype} {site}", flush=True)
    for k in range(n, nrow * ncol):
        ax[k // ncol][k % ncol].axis("off")
    ncorr = sum(1 for r in rows if r[3].get("correctable") == "1")
    fig.suptitle(f"{itype} clear-day time-height quicklook — all {n} stations ({ncorr} correctable, green). "
                 f"Time (x, hour UTC) vs height (y, 0–{ZTOP/1000:.0f} km); log β_att (nominal calibration).",
                 fontweight="bold", fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    out = D / f"fig_offset_pcolor_{itype}.png"
    fig.savefig(out, dpi=135); fig.savefig(REPORT_FIG / out.name, dpi=135)
    print(f"saved {out.name}", flush=True)
print("OFFSET_PCOLOR_DIAGNOSTICS_DONE")
