# Sample data

Small, trimmed **real** E-PROFILE L1/L2 NetCDF fixtures — one or two per instrument
type — committed so the package can be exercised without the full archive. They drive
`tests/test_sample_data.py` and `examples/03_sample_data_calibration.ipynb`.

The layout mirrors the operational archive, so point `folder_root` here and set
`data_level` accordingly:

```
examples/data/
├── L1/<WMO>/2026/<MM>/L1_<WMO>_A<YYYYMMDD>.nc      # data_level = L1
└── L2/<WMO>/2026/<MM>/L2_<WMO>_A<YYYYMMDD>.nc      # data_level = L2_daily
```

| Instrument (λ) | WMO id (site) | L1 day | L2 day(s) | condition | calibration on the L2 fixture |
|---|---|---|---|---|---|
| CHM15k (1064 nm) | `0-20000-0-06610` (Payerne) | 2026-02-25 | 2026-02-24 eve + 2026-02-25 | clear night | **Rayleigh → flag 1** (self-contained) |
| Mini-MPL (532 nm) | `0-20000-0-07014` (Lille) | 2026-04-23 | 2026-04-23 eve + **2026-04-24** | clear night | **Rayleigh → flag 1** (self-contained) |
| CL61 (910 nm) | `0-756-4-EERLCL61` (Sion) | 2026-03-04 | 2026-03-03 eve + 2026-03-04 | clear night | **Rayleigh → flag 1** (needs CAMS) |
| CL61 (910 nm) | `0-756-4-EERLCL61` (Sion) | — | 2026-04-14 | overcast | **Cloud →** runs, n≈184 in-cloud profiles (needs CAMS) |
| CL31 (910 nm) | `0-20000-0-06602` (Delémont) | 2026-02-20 | 2026-02-20 | overcast | **Cloud → C ≈ 1.6** (needs CAMS) |
| CL51 (910 nm) | `0-20000-0-02998` (Helsinki) | 2026-04-09 | 2026-01-16 | overcast | **Cloud → C ≈ 3.6** (needs CAMS) |

**Trimming** (to keep the repo small): L2 Rayleigh fixtures are the FULL day D plus the
≥16 UTC evening slice of day D−1 (the calibration night is D−1 evening → D morning, and
the darkness-adaptive SZA window needs both files), range ≤ 12 km; cloud days are kept
full-day. L1 fixtures are decimated (they back the reader test only). All are real data,
just cropped. Re-curated 2026-07: nights re-selected to pass the current production
settings (eprof_v2 + SZA night + post-0fb4cc4 cloud screening on coarse L2 data).

**Water vapour (910 nm).** CL61/CL31/CL51 need a monthly CAMS file for the mandatory WV
correction, so their tests **skip when CAMS is absent**. The HITRAN WV cross-section LUT
is bundled (`calibration/data/abs_cross_wv_910nm.nc`) — no external LUT required.

**Notes.** CL51 L1 only exists as short sub-daily chunks, so its L1 fixture is one such
chunk (a different day from its cloud L2). CL61 (Sion) reliably Rayleigh-calibrates but
rarely meets the warm-liquid-stratocumulus conditions the O'Connor cloud method needs, so
its cloud fixture asserts only that the cloud path runs and finds in-cloud profiles.
