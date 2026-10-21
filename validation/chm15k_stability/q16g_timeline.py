# -*- coding: utf-8 -*-
"""q16g (task 1b): timeline of the network-median CHM15k residual 2025-2026 against the processing /
firmware changes found in the L1 attributes (q16g_l1_provenance.py).

* Network-median monthly residual (ln C - kappa_in (T_int - median), unit mean removed) and the same for
  raw ln C and for the aerosol-free clear-sky signal X(7-8 km) (q15e series), with the per-unit annual
  harmonic removed too (= the non-seasonal "steps").
* Change dates (network-wide): raw2l1 3.1.2 -> 3.2.2 (2025-03), hermes-eprofile-alc-science 2.0.11 -> .12
  (2025-01), .13/.14 (2025-07), .15 (2025-10); DWD firmware 1.14 -> 1.15 on ~57 units (2026-04-29 .. 06-03).
* Difference-in-differences for the firmware switch: per unit, mean residual 30 nights after minus 30 nights
  before the switch date in 2026, minus the same calendar windows in 2025 (removes the season); switched units
  vs the others (pseudo-date 2026-05-08).

Outputs: q16g_timeline_months.csv, q16g_firmware_did.csv, q16g_timeline.png
"""
from __future__ import annotations

import csv
import sys
from datetime import datetime, timedelta
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1]))
sys.path.insert(0, str(HERE))
from q16b_insolation import PAYERNE, Q15, Q16, harmonic, load_kappas, residuals  # noqa: E402

EPOCH = datetime(1970, 1, 1)
CHANGES = [("2025-01-08", "hermes-science 2.0.12"), ("2025-03-12", "raw2l1 3.2.2 / hermes-raw2l1 2.0.4"),
           ("2025-07-09", "hermes 2.0.5 / science 2.0.13-14"), ("2025-10-08", "hermes-science 2.0.15"),
           ("2026-05-06", "DWD firmware 1.14->1.15 (57 units)")]


def day(s):
    return (datetime.strptime(s, "%Y-%m-%d") - EPOCH).days


def x78_series(ser, skeys, s, r, col="x78"):
    if s not in skeys:
        return np.full(len(r["t"]), np.nan)
    idx = {round(float(v), 9): i for i, v in enumerate(ser[f"{s}__lnc"])}
    j = np.array([idx.get(round(float(v), 9), -1) for v in r["lnC"]])
    src = ser[f"{s}__lnx"][:, -1] if col == "x78" else ser[f"{s}__lnx_win"]
    return np.where(j >= 0, src[np.clip(j, 0, None)], np.nan)


def firmware_dates():
    out = {}
    for c in csv.DictReader(open(Q16 / "q16g_l1_attr_changes.csv", encoding="utf-8")):
        if c["attribute"] == "instrument_firmware_version" and c["old"] == "1.14" and c["new"] == "1.15":
            out[c["stream"]] = (datetime.strptime(c["date"], "%Y%m%d") - EPOCH).days
    return out


def main():
    units = list(np.load(Q15 / "q15f_nightly.npz", allow_pickle=True)["units"])
    k_in, k_m2, _ = load_kappas()
    ser = np.load(Q15 / "q15e_series.npz")
    skeys = set(k.split("__")[0] for k in ser.files)
    fw = firmware_dates()
    months = [(y, m) for y in (2025, 2026) for m in range(1, 13) if (y, m) <= (2026, 9)]
    mstart = np.array([(datetime(y, m, 1) - EPOCH).days for y, m in months] + [(datetime(2026, 10, 1) - EPOCH).days])
    series_names = ("raw", "res_in", "res_in_noharm", "x78", "x78_noharm")
    acc = {k: [[] for _ in months] for k in series_names}
    did_rows = []
    for m, r in units:
        s = m["stream"]
        if s == PAYERNE:
            continue
        t = r["t"]
        kin = k_in.get(s, np.nan)
        res = residuals(r, kin) if np.isfinite(kin) else np.full(len(t), np.nan)
        x78 = x78_series(ser, skeys, s, r)
        ys = {"raw": r["lnC"], "res_in": res, "x78": x78}
        for k in ("res_in", "x78"):
            _, e, _ = harmonic(t, ys[k])
            ys[f"{k}_noharm"] = e
        for k in series_names:
            y = ys[k]
            if np.isfinite(y).sum() < 60:
                continue
            y = y - np.nanmean(y)
            for i in range(len(months)):
                sel = (t >= mstart[i]) & (t < mstart[i + 1]) & np.isfinite(y)
                if sel.sum() >= 4:
                    acc[k][i].append(100 * float(np.mean(y[sel])))
        # firmware difference-in-differences (2026 step minus the same calendar step in 2025)
        d0 = fw.get(s, day("2026-05-08"))
        for win in (30, 45):
            row = dict(stream=s, switched=int(s in fw), date=f"{EPOCH + timedelta(days=d0):%Y-%m-%d}", win=win)
            for k in ("raw", "res_in", "x78"):
                y = ys[k]
                vals = []
                for shift in (0, -365):
                    a = (t >= d0 + shift - win) & (t < d0 + shift) & np.isfinite(y)
                    b = (t >= d0 + shift) & (t < d0 + shift + win) & np.isfinite(y)
                    vals.append(100 * (y[b].mean() - y[a].mean()) if (a.sum() >= 5 and b.sum() >= 5) else np.nan)
                row[f"step2026_{k}"], row[f"step2025_{k}"] = vals
                row[f"did_{k}"] = vals[0] - vals[1]
            did_rows.append(row)
    rows = []
    rng = np.random.default_rng(3)
    for i, (y, mo) in enumerate(months):
        row = dict(year=y, month=mo)
        for k in series_names:
            v = np.array(acc[k][i])
            if v.size >= 10:
                med = float(np.median(v))
                bs = [np.median(rng.choice(v, v.size)) for _ in range(300)]
                row.update({f"{k}_med": med, f"{k}_se": float(np.std(bs)), f"{k}_n": int(v.size),
                            f"{k}_frac_same_sign": float(np.mean(np.sign(v) == np.sign(med)))})
        rows.append(row)
        print({kk: (round(vv, 2) if isinstance(vv, float) else vv) for kk, vv in row.items()
               if kk.endswith("_med") or kk in ("year", "month") or kk == "res_in_n"})
    write_csv(rows, Q16 / "q16g_timeline_months.csv")
    write_csv(did_rows, Q16 / "q16g_firmware_did.csv")
    from scipy.stats import mannwhitneyu
    for win in (30, 45):
        for k in ("raw", "res_in", "x78"):
            a = np.array([d[f"did_{k}"] for d in did_rows if d["win"] == win and d["switched"]], float)
            b = np.array([d[f"did_{k}"] for d in did_rows if d["win"] == win and not d["switched"]], float)
            a, b = a[np.isfinite(a)], b[np.isfinite(b)]
            if a.size > 5 and b.size > 5:
                p = mannwhitneyu(a, b).pvalue
                print(f"firmware DiD win {win} {k}: switched {np.median(a):+.2f} (n={a.size}) vs others "
                      f"{np.median(b):+.2f} (n={b.size}); diff {np.median(a) - np.median(b):+.2f} %, MWU p={p:.3f}")
    figure(rows, months)


def figure(rows, months):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    x = np.array([(datetime(y, m, 15) - EPOCH).days for y, m in months])
    fig, ax = plt.subplots(1, 2, figsize=(17, 5.2))
    for k, lab, c in (("raw", "raw ln C_L", "0.5"), ("res_in", "residual (kappa_in)", "C0"), ("x78", "clear-sky X(7-8 km)", "C2")):
        v = np.array([r.get(f"{k}_med", np.nan) for r in rows])
        e = np.array([r.get(f"{k}_se", np.nan) for r in rows])
        ax[0].errorbar([EPOCH + timedelta(days=float(d)) for d in x], v, yerr=e, fmt="-o", ms=4, color=c, label=lab)
    for k, lab, c in (("res_in_noharm", "residual - annual harmonic", "C0"), ("x78_noharm", "X(7-8 km) - annual harmonic", "C2")):
        v = np.array([r.get(f"{k}_med", np.nan) for r in rows])
        e = np.array([r.get(f"{k}_se", np.nan) for r in rows])
        ax[1].errorbar([EPOCH + timedelta(days=float(d)) for d in x], v, yerr=e, fmt="-o", ms=4, color=c, label=lab)
    for a in ax:
        for d, lab in CHANGES:
            dd = datetime.strptime(d, "%Y-%m-%d")
            a.axvline(dd, color="r", lw=0.8, ls="--")
            a.text(dd, a.get_ylim()[1] if False else 0, " " + lab, rotation=90, fontsize=7, color="r", va="bottom")
        a.axhline(0, color="k", lw=0.5)
        a.set_ylabel("network median, unit mean removed [%]")
        a.legend(fontsize=8, loc="lower left")
    ax[0].set_title("(a) monthly network median (135 CHM15k, Payerne A excluded)")
    ax[1].set_title("(b) after removing each unit's annual harmonic: non-seasonal steps")
    fig.suptitle("q16g - CHM15k monthly residual vs processing / firmware change dates (red)", fontsize=10)
    fig.autofmt_xdate()
    fig.tight_layout()
    fig.savefig(Q16 / "q16g_timeline.png", dpi=130)


def write_csv(rows, path):
    keys = []
    for r in rows:
        for k in r:
            if k not in keys:
                keys.append(k)
    with open(path, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=keys)
        w.writeheader()
        w.writerows(rows)


if __name__ == "__main__":
    main()
