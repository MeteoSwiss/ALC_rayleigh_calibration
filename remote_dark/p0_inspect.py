# -*- coding: utf-8 -*-
"""P0 triage inspection — one flagged stream, four diagnostic views (fig21_<key>.png).

Usage:  python -m remote_dark.p0_inspect <path-to-cache-npz> [label]

Views: (1) the evidence intercept d(z) against Payerne A's evidence and hood (the scale and
shape references); (2) the nightly constants vs time (steps = hardware eras, drift = optics);
(3) the median measured signal vs the predicted molecular signal (where the atmosphere
departs from the model — persistent haze shows here); (4) the per-night noise floor
(hardware health / swap tripwire). Numbers: the P0 bound recomputed locally + the evidence
amplitude relative to Payerne's.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from datetime import datetime

from remote_dark.common import DATA, PAYERNE, ensure_out
from remote_dark.p0_net_cache_balfrin import nightly_c, evidence_and_bound
from remote_dark import hood

FIG = ensure_out("figs")
M1_DIR = DATA / "remote_dark" / "m1"


def smooth(x, k=7):
    pad = np.pad(x, k // 2, mode="edge")
    return np.array([np.nanmedian(pad[i:i + k]) for i in range(x.size)])


def main():
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    f = Path(sys.argv[1])
    label = sys.argv[2] if len(sys.argv) > 2 else f.stem
    z = np.load(f, allow_pickle=True)
    S, sig, M, rng, dates = z["S"], z["sigma"], z["M"], z["rng"], z["dates"]
    z2 = np.where(rng > 0, rng, np.nan) ** 2
    r = evidence_and_bound(z)
    c = nightly_c(S, sig, M, rng)
    ok = np.isfinite(c) & (c > 0)

    # Payerne references: sky evidence (same estimator class) + hood truth
    m1 = np.load(M1_DIR / f"{PAYERNE['wmo']}_A_m1.npz", allow_pickle=True)
    pay_rng = m1["rng"]
    pay_free = m1["b_free"] / np.where(pay_rng > 0, pay_rng, np.nan) ** 2
    pay_hood = m1["hood_b"] / np.where(pay_rng > 0, pay_rng, np.nan) ** 2

    fig, (a1, a2, a3, a4) = plt.subplots(1, 4, figsize=(16.5, 6.0))

    mm = (rng >= 500) & (rng <= 10000)
    pm = (pay_rng >= 500) & (pay_rng <= 10000)
    a1.plot(smooth(r["d"] / z2)[mm] * 1e3, rng[mm] / 1e3, "-", color="crimson", lw=1.9,
            label=f"{label} evidence d(z)")
    a1.plot(smooth(pay_free)[pm] * 1e3, pay_rng[pm] / 1e3, "-", color="steelblue", lw=1.2,
            label="Payerne A evidence (sky)")
    a1.plot(smooth(pay_hood)[pm] * 1e3, pay_rng[pm] / 1e3, "k--", lw=1.6,
            label="Payerne A hood truth")
    a1.axvline(0, color="k", lw=0.5)
    a1.set_xlabel("intercept, P-view (×1e-3)")
    a1.set_ylabel("altitude (km)")
    ref = np.concatenate([smooth(r["d"] / z2)[mm], smooth(pay_free)[pm]]) * 1e3
    span = float(np.nanpercentile(np.abs(ref), 97))
    a1.set_xlim(-2.5 * span, 2.5 * span)
    a1.legend(fontsize=8, loc="best")
    a1.grid(alpha=0.25)
    a1.set_title("evidence vs the Payerne references", fontsize=10)

    tnum = mdates.date2num([datetime.strptime(str(int(x)), "%Y%m%d") for x in dates])
    a2.plot(tnum[ok], c[ok], ".", color="k", ms=4)
    a2.plot(tnum[ok], smooth(c[ok], 15), "-", color="crimson", lw=1.5)
    a2.xaxis.set_major_formatter(mdates.DateFormatter("%y-%m"))
    a2.set_ylabel("nightly Rayleigh constant (4.5-6.5 km)")
    a2.set_yscale("log")
    a2.grid(alpha=0.25)
    a2.set_title(f"nightly constants, {int(ok.sum())} nights", fontsize=10)

    medS = np.nanmedian(S / z2[None, :], axis=0)
    medCM = np.nanmedian(c[ok, None] * M[ok] / z2[None, :], axis=0)
    mm3 = (rng >= 300) & (rng <= 15000)
    a3.plot(smooth(medS)[mm3], rng[mm3] / 1e3, "-", color="k", lw=1.6, label="measured (median)")
    a3.plot(smooth(medCM)[mm3], rng[mm3] / 1e3, "--", color="crimson", lw=1.6,
            label="predicted molecular C·M")
    a3.set_xscale("log")
    a3.set_xlabel("signal, P-view")
    a3.legend(fontsize=8, loc="best")
    a3.grid(alpha=0.25, which="both")
    a3.set_title("signal vs molecular model\n(departure = haze/overlap/dark)", fontsize=10)

    n_prof = np.maximum(z["n_prof"].astype(float), 1.0)
    top = rng >= np.nanpercentile(rng, 92)
    fl = np.nanmedian((sig * np.sqrt(n_prof)[:, None] / z2[None, :])[:, top], axis=1)
    a4.plot(tnum, fl, ".-", color="steelblue", ms=3, lw=0.8)
    a4.set_yscale("log")
    a4.xaxis.set_major_formatter(mdates.DateFormatter("%y-%m"))
    a4.set_ylabel("noise floor per night, P-view")
    a4.grid(alpha=0.25, which="both")
    a4.set_title("noise floor (hardware tripwire)", fontsize=10)

    # numbers
    amp = float(np.nanmedian(np.abs(smooth(r["d"] / z2)[mm]))
                / np.nanmedian(np.abs(smooth(pay_free)[pm])))
    print(f"{label}: nights={r['n_nights']}  bound={100 * r['bound']:+.1f} %  "
          f"evidence amplitude = x{amp:.1f} Payerne-A-evidence")
    fig.suptitle(f"P0 triage inspection — {label}  (bound {100 * r['bound']:+.1f} %, "
                 f"evidence ×{amp:.1f} the Payerne-A sky evidence)", fontsize=11)
    fig.tight_layout(rect=[0, 0, 1, 0.93])
    out = FIG / f"fig21_{label}.png"
    fig.savefig(out, dpi=140)
    plt.close(fig)
    print(out)


if __name__ == "__main__":
    main()
