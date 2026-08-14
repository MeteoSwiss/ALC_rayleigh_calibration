# -*- coding: utf-8 -*-
"""Canary nights — calibrations that MUST stay rejected, whatever we do to the gates.

The whole risk of this work is that loosening/reformulating the QC lets outliers back in. These are
nights where the old operational v1 produced a constant that is demonstrably wrong and the current
v2 correctly rejects; any candidate configuration that accepts one of them is disqualified, no
matter how good its availability numbers look.

Registry (Payerne CHM15k 0-20000-0-06610_A, from the v1-vs-v2 forensics):
  20260531  v1 C_L = 3.163e12  ~6x the station norm (~5.3e11), v1 claimed only 2.7 % uncertainty
  20260608  v1 C_L = 9.781e11  ~1.9x the norm, claimed 4.4 %

Add to CANARIES whenever a new confirmed outlier is found. `assert_canaries` is called by every
sweep-analysis script BEFORE any configuration may be reported.
"""
from __future__ import annotations

# (wmo, ident, date, why) — the night must NOT yield a valid calibration
CANARIES = [
    ("0-20000-0-06610", "A", "20260531", "v1 C_L 3.163e12 = 6x station norm (v1 claimed 2.7% unc)"),
    ("0-20000-0-06610", "A", "20260608", "v1 C_L 9.781e11 = 1.9x station norm (v1 claimed 4.4% unc)"),
]

CANARY_KEYS = {(w, i, d) for w, i, d, _ in CANARIES}


def is_canary(wmo, ident, date_str):
    return (wmo, ident, str(date_str)) in CANARY_KEYS


def check(results):
    """results: {(wmo, ident, date): {config: accepted_bool}} -> {config: [failed canary, ...]}.

    A canary "fails" when a config accepts it. Returns only the configs with failures.
    """
    bad = {}
    for (wmo, ident, date), per_cfg in results.items():
        if not is_canary(wmo, ident, date):
            continue
        for cfg, accepted in per_cfg.items():
            if accepted:
                bad.setdefault(cfg, []).append(f"{wmo}_{ident} {date}")
    return bad


def assert_canaries(results, *, strict=True):
    """Print the canary verdict; raise (strict) if any config accepted a canary night.

    Also warns when a canary was never evaluated — an unevaluated canary is not a pass.
    """
    seen = {k for k in results if is_canary(*k)}
    missing = CANARY_KEYS - seen
    for wmo, ident, date in sorted(missing):
        print(f"  [canary] NOT EVALUATED: {wmo}_{ident} {date} (corpus/period does not cover it)")
    bad = check(results)
    if not bad:
        print(f"  [canary] OK - all {len(seen)}/{len(CANARY_KEYS)} evaluated canaries rejected "
              f"by every config")
        return True
    for cfg, nights in sorted(bad.items()):
        print(f"  [canary] FAIL: config {cfg!r} ACCEPTED {', '.join(nights)}")
    if strict:
        raise AssertionError(f"{len(bad)} config(s) accepted a canary night: {sorted(bad)}")
    return False


if __name__ == "__main__":
    print(f"{len(CANARIES)} canary nights registered:")
    for wmo, ident, date, why in CANARIES:
        print(f"  {wmo}_{ident} {date}  {why}")
