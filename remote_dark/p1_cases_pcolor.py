# -*- coding: utf-8 -*-
"""Fig 24 — time-height backscatter quicklooks of the P1 case families on their transition
days, with the three data discriminators annotated per panel:

  1. the instrument's own error_ext service code (operational decoder);
  2. window_transmission housekeeping;
  3. the FIRST-GATES brightness ratio — the physical separator the operator asked for:
     fog and snow-on-window SCATTER the transmit pulse into the receiver (bright first
     gates, nothing above), a cover/hood absorbs it (dark everywhere).

Cases: BERUS true cover; EXETER detector fault; PAYERNE hood (positive control);
Schiphol on the synchronised NL snow day; and Kleine Scheidegg twice — operator ground
truth: NO hood was ever placed there, every flagged day is snow (negative control).
"""
from __future__ import annotations

import sys
from collections import Counter
from datetime import datetime, timedelta
from pathlib import Path

import numpy as np
import netCDF4
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

from remote_dark.common import ensure_out, l1_file
from calibration.status.decode import decode_chm15k_error_ext

FIG = ensure_out("figs")
SCR = Path(r"C:/Users/hervo/AppData/Local/Temp/claude"
           r"/C--Users-hervo-OneDrive-Documents-ALC-rayleigh-calibration"
           r"/265f36ab-1a40-4592-bae6-6c3550f2fac1/scratchpad")

CASES = [
    (SCR / "L1_0-20000-0-10704_020260203.nc",
     "BERUS — TRUE COVER goes on (2026-02-03)"),
    (SCR / "L1_0-20000-0-03838_A20250724.nc",
     "EXETER — detector fault begins (2025-07-24)"),
    (l1_file("0-20000-0-06610", "A", datetime(2026, 5, 26)),
     "PAYERNE — hood session, positive control (2026-05-26)"),
    (l1_file("0-20000-0-06240", "A", datetime(2025, 12, 25)),
     "SCHIPHOL A — synchronised NL snow day (2025-12-25)"),
    (SCR / "L1_0-20000-0-06736_A20231203.nc",
     "KLEINE SCHEIDEGG — snow (op. ground truth), 2023-12-03"),
    (SCR / "L1_0-20000-0-06736_A20190408.nc",
     "KLEINE SCHEIDEGG — snow (op. ground truth), 2019-04-08"),
]


def main():
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    fig, axes = plt.subplots(2, 3, figsize=(16.5, 9.2), sharey=True)
    for ax, (f, lab) in zip(axes.ravel(), CASES):
        with netCDF4.Dataset(f) as ds:
            t = np.ma.filled(ds.variables["time"][:].astype("f8"), np.nan)
            rng = np.ma.filled(ds.variables["range"][:].astype("f8"), np.nan)
            rcs = np.ma.filled(ds.variables["rcs_0"][::2, :].astype("f4"), np.nan)
            err = (np.ma.filled(ds.variables["error_ext"][::2].astype("f8"), np.nan)
                   if "error_ext" in ds.variables else None)
            wt = (float(np.nanmedian(np.ma.filled(
                ds.variables["window_transmission"][:].astype("f8"), np.nan)))
                if "window_transmission" in ds.variables else np.nan)
        hours = np.array([(datetime(1970, 1, 1) + timedelta(days=float(x))).hour
                          + (datetime(1970, 1, 1) + timedelta(days=float(x))).minute / 60.0
                          for x in t[::2]])
        mz = (rng > 0) & (rng <= 8000)
        pos = rcs[:, mz]
        finite = np.isfinite(pos) & (pos < 1e30)
        vals = pos[finite & (pos > 0)]
        vmax = float(np.percentile(vals, 99.5)) if vals.size else 1.0
        if not np.isfinite(vmax) or vmax <= 0:
            vmax = 1.0
        cm = plt.cm.viridis.copy()
        cm.set_bad("0.25")
        pc = ax.pcolormesh(hours, rng[mz] / 1e3,
                           np.where(finite & (pos > 0), pos, np.nan).T,
                           norm=LogNorm(vmin=vmax / 1e4, vmax=vmax),
                           cmap=cm, shading="nearest", rasterized=True)
        fig.colorbar(pc, ax=ax, pad=0.01)
        ax.set_xlabel("hour (UTC)")
        ax.set_xlim(0, 24)

        # discriminator 3: first-gates vs mid-band brightness (fog/snow scatter vs cover dark)
        z = np.where(rng > 0, rng, np.nan)
        P = np.abs(np.where(np.isfinite(rcs) & (np.abs(rcs) < 1e30), rcs, np.nan)
                   / z[None, :] ** 2)
        first = np.nanmedian(P[:, (rng > 0) & (rng <= 200)])
        band = np.nanmedian(P[:, (rng >= 400) & (rng <= 3000)])
        ratio = first / band if band > 0 else np.nan

        estr = ""
        if err is not None and np.isfinite(err).any():
            cnt = Counter()
            for v in err[np.isfinite(err)]:
                names = decode_chm15k_error_ext(v)
                if not names:
                    cnt["OK"] += 1
                for nm, _sev in names:
                    cnt[nm.replace("Error: ", "")] += 1
            tot = int(np.isfinite(err).sum())
            estr = "; ".join(f"{nm} {100 * c / tot:.0f}%" for nm, c in cnt.most_common(2))
        ax.set_title(f"{lab}\nerr: {estr}\n"
                     f"window_trans={wt:.0f}   first-gates/band ×{ratio:.0f}", fontsize=8)
        print(f"{lab}\n  err: {estr}\n  wt={wt:.0f}  first/band x{ratio:.0f}")
    for ax in axes[:, 0]:
        ax.set_ylabel("altitude (km)")
    fig.suptitle("P1 case families — the three discriminators: cover = dark at ALL gates; "
                 "snow/fog = BRIGHT first gates (pulse scattered back), dark above; "
                 "fault = flags + dead receiver", fontsize=11)
    fig.tight_layout(rect=[0, 0, 1, 0.93])
    f = FIG / "fig24_cases_pcolor.png"
    fig.savefig(f, dpi=140)
    plt.close(fig)
    print(f)


if __name__ == "__main__":
    main()
