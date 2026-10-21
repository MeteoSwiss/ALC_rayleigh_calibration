# -*- coding: utf-8 -*-
"""Review (code): the molecular climatology where it is extrapolated -- the southern-hemisphere
season mirroring (Lauder, 45 S) and the low-latitude clip (Bonaire, 12 N).

For every calibrated night of these fresh units (windows from the review run, eprof_v2.3) the
realised ln-bias b of C_L from US-1976 instead of the real T/p is computed exactly as q15d does
(slope of the CAMS molecular signal against the US-1976 one inside the night's window), with T/p
from the regional CAMS boxes (calibration.rayleigh.atmosphere.load_cams_atmosphere, available from
2026-06-24), and compared with the climatology the budget uses.

Output: u23r_molecular_south.txt
"""
from __future__ import annotations

import sys
import warnings
from datetime import datetime, timedelta
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / "validation" / "chm15k_stability"))
sys.path.insert(0, str(HERE))
warnings.filterwarnings("ignore")

UNITS = [("0-20008-0-LAU", "A", 910e-9), ("0-20000-0-78990", "A", 1064e-9)]


def main():
    import u23r_common as C
    from calibration.io.cams import ensure_cams_file
    from calibration.rayleigh.atmosphere import (calculate_molecular_properties, load_cams_atmosphere,
                                                 load_standard_atmosphere)
    from calibration.rayleigh import uncertainty as U
    d = C.load_runs("fresh", "eprof_v2.3")
    cen = C.A.census()
    lines = []
    for wmo, ident, wl in UNITS:
        s = cen[(wmo, ident)]
        g = d[(d.key == f"{wmo}_{ident}") & d.ok & (d.date >= "20260624")]
        rows = []
        for r in g.itertuples():
            f = ensure_cams_file("D:/CAMS_Monthly_04;D:/CAMS_daily", r.date, auto_download=False,
                                 latitude=s["lat"], longitude=s["lon"])
            if f is None:
                continue
            night = np.datetime64(f"{r.date[:4]}-{r.date[4:6]}-{r.date[6:]}")
            z = s["alt"] + np.arange(15.0, 9000.0, 30.0)
            atm = load_cams_atmosphere(f, s["lat"], s["lon"], night - np.timedelta64(4, "h"),
                                       night + np.timedelta64(4, "h"), z)
            if atm is None:
                continue
            std = load_standard_atmosphere(None, z)
            rng = z - s["alt"]
            mc = calculate_molecular_properties(atm.temperature, atm.pressure, rng, wl)
            ms = calculate_molecular_properties(std.temperature, std.pressure, rng, wl)
            w = (rng >= r.zb_agl) & (rng <= r.zt_agl)
            if w.sum() < 5:
                continue
            x, y = ms.p_mol[w], mc.p_mol[w]
            xm = x - x.mean()
            b = np.log(float(np.sum(xm * (y - y.mean())) / np.sum(xm * xm)))
            doy = datetime.strptime(r.date, "%Y%m%d").timetuple().tm_yday
            zc_asl = (0.5 * (r.zb_agl + r.zt_agl) + s["alt"]) / 1000.0
            clim = U.molecular_bias_climatology(s["lat"], doy, zc_asl, s["alt"] / 1000.0)
            rows.append((b, clim))
        if not rows:
            lines.append(f"{s['site']}: no night with CAMS T/p")
            continue
        b, clim = (np.array(v) for v in zip(*rows))
        sig = np.sqrt(clim ** 2 + 0.0135 ** 2)
        lines.append(f"{s['site']} ({s['lat']:.1f}): {len(b)} nights 2026-06-24.., realised b mean "
                     f"{100 * b.mean():+.2f} % (sd {100 * b.std():.2f} %), climatology mean "
                     f"{100 * clim.mean():+.2f} %, error {100 * (b - clim).mean():+.2f} %; |b| <= 2 sigma_mol on "
                     f"{100 * np.mean(np.abs(b) <= 2 * sig):.0f} % of the nights")
    print("\n".join(lines))
    (C.OUT / "u23r_molecular_south.txt").write_text("\n".join(lines), encoding="utf-8")


if __name__ == "__main__":
    main()
