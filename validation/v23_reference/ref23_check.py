# -*- coding: utf-8 -*-
"""ref23: completeness check of the harness output against the L1 archive as it is NOW.

For every site-month pickle and every night, each stream whose day-D L1 file exists must have been
read (t_read) and processed (Rayleigh rows for Rayleigh streams, a cloud row for cloud streams).
A night where a stream is missing although its file exists now was processed while the file was
unavailable (e.g. a temporarily moved folder) and its month must be re-run.  The task start times
are reconstructed from the run log (completion time = pickle mtime, duration from the log line).

Usage: python validation/v23_reference/ref23_check.py [--delete]
       --delete removes the pickles of incomplete months so that `ref23_run.py run` redoes them.
"""
from __future__ import annotations

import pickle
import re
import sys
from datetime import datetime, timedelta
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import ref23_common as RC  # noqa: E402


def task_times():
    """{(site, month): (start, end)} from run.log lines '[k/N] site month: n nights in T s'."""
    out = {}
    log = RC.OUT / "run.log"
    if not log.exists():
        return out
    for ln in log.read_text(encoding="utf-8", errors="ignore").splitlines():
        m = re.match(r"\[\d+/\d+\] (\w+) (\d{6}): \d+ nights in (\d+) s", ln)
        if not m:
            continue
        site, month, dur = m.group(1), m.group(2), int(m.group(3))
        f = RC.RUN_DIR / site / f"{month}.pkl"
        if f.exists():
            end = datetime.fromtimestamp(f.stat().st_mtime)
            out[(site, month)] = (end - timedelta(seconds=dur), end)
    return out


def main(delete=False):
    times = task_times()
    bad = {}
    n_nights = 0
    for f in sorted(RC.RUN_DIR.glob("*/*.pkl")):
        site, month = f.parent.name, f.stem
        s = RC.SITES[site]
        with open(f, "rb") as fh:
            recs = pickle.load(fh)
        for r in recs:
            n_nights += 1
            d = datetime.strptime(r["date"], "%Y%m%d").date()
            for ident, (itype, methods) in s["streams"].items():
                if not RC.l1_file(s["wmo"], ident, d).exists():
                    continue
                miss = []
                if "fatal" in r:
                    miss.append("fatal")
                elif ident not in r.get("t_read", {}):
                    miss.append("not read")
                else:
                    if "ray" in methods and ident not in r.get("ray", {}):
                        miss.append("no Rayleigh")
                    if "cloud" in methods and ident not in r.get("cloud", {}):
                        miss.append("no cloud")
                    if r.get("err", {}).get(ident):
                        miss.append(r["err"][ident][:60])
                if miss:
                    bad.setdefault((site, month), []).append((r["date"], ident, ",".join(miss)))
    print(f"{n_nights} site-nights checked")
    for (site, month), lst in sorted(bad.items()):
        t = times.get((site, month))
        print(f"INCOMPLETE {site} {month} ({len(lst)} stream-nights)"
              + (f" task {t[0]:%H:%M:%S}-{t[1]:%H:%M:%S}" if t else ""))
        for x in lst[:5]:
            print("    ", x)
        if delete:
            (RC.RUN_DIR / site / f"{month}.pkl").unlink()
    for (site, month), (t0, t1) in sorted(times.items(), key=lambda kv: kv[1][0]):
        if site == "schiphol":
            print(f"schiphol {month}: ran {t0:%H:%M:%S} - {t1:%H:%M:%S}")
    if not bad:
        print("all site-months complete")


if __name__ == "__main__":
    main(delete="--delete" in sys.argv)
