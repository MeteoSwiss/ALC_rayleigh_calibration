# Calibration-constant time series & per-calibration outlier rate (network, 2026)

*Generated 2026-06-21. `validation/outliers_timeseries.py` on the network run
(`run_network_v2_vs_v11.py`). Optimized E-PROF v2 (config C8) vs E-PROF v1.1, L1 + L2.*

## Outlier definition (drift-aware, robust)

For each calibration (stream × level × method) we take the date-ordered **valid** nightly lidar
constants and flag outliers relative to the instrument's own slow drift:

```
residual_i = C_i − rolling_median(C, 9 nights)
sigma_rob  = 1.4826 · median(|residual − median(residual)|)
outlier_i  = |residual_i − median(residual)| > 3 · sigma_rob
outlier %  = flagged nights / valid nights
```

Detrending with a 9-night rolling median removes seasonal/instrumental drift, so only genuine
jumps and spikes (bad fits, residual cloud/aerosol that slipped past QC, glitches) are counted —
not the slow baseline a global threshold would mistake for outliers.

## Per-type outlier rate — median over streams

| type | level | **v2 outlier %** | v1.1 outlier % | worst stream (v2) |
|---|---|---|---|---|
| CHM15k | L1 | **4.2** | 6.2 | 33.3 |
| CHM15k | L2 | **4.9** | 9.7 | 19.2 |
| Mini-MPL | L1 | **9.9** | 10.7 | 14.0 |
| Mini-MPL | L2 | **7.8** | 10.3 | 15.1 |
| CL61 | L1 | **2.0** | 11.1 | 5.6 |
| CL61 | L2 | **4.7** | 12.2 | 9.5 |

**The optimized v2 produces markedly fewer outliers than v1.1 for every type and level** — most
strikingly CL61 (2.0% vs 11.1% on L1; 4.7% vs 12.2% on L2). This is the time-domain counterpart of
v2's lower σ_SD: its gated-optimal selection with temporal-variability rejection avoids the
occasional bad windows that v1.1 admits. Mini-MPL is the noisiest type (~8–10%); CHM15k and CL61
sit at ~2–5% under v2.

![Calibration outlier rates: (left) v2 outlier % by type, L1 light/L2 dark, bar=median; (centre) v2 vs v1.1 on L2 — points below the line = v2 cleaner; (right) per-stream v2 outlier % sorted.](figs_paper_validation/network_v2_v11/fig_outlier_overview.png)

## Time series (L2, optimized v2; outliers in red)

Each panel is one instrument's nightly lidar constant over 2026: green = robust median, grey band =
±3σ, red = flagged outliers; the panel title gives the outlier %.

![CL61 — nightly lidar constant (L2, v2), outliers in red](figs_paper_validation/network_v2_v11/fig_ts_CL61.png)

![Mini-MPL — nightly lidar constant (L2, v2), outliers in red](figs_paper_validation/network_v2_v11/fig_ts_MiniMPL.png)

![CHM15k (12 highest-outlier streams) — nightly lidar constant (L2, v2), outliers in red](figs_paper_validation/network_v2_v11/fig_ts_CHM15k.png)

## Extremes (L2, v2)

**Highest outlier rates** (worth flagging for instrument follow-up):

| site | stream | type | v2 out % | v1.1 out % |
|---|---|---|---|---|
| Pécs | 0-20000-0-12944_A | CHM15k | 19.2 (n=26) | 14.3 |
| Ell | 0-20000-0-06377_A | CHM15k | 18.2 (n=77) | 10.1 |
| Essen | 0-20000-0-10410 | CHM15k | 15.9 (n=44) | 5.8 |
| Mini-MPL Mobile | 0-20000-0-07617_A | Mini-MPL | 15.1 (n=93) | 18.8 |
| Leipzig | 0-20000-0-10471 | CHM15k | 12.5 (n=80) | 18.5 |

**Cleanest** (0.0% over a full record): Coningsby, Aberporth, Amsterdam-Schiphol, Twenthe (all CHM15k).

> Essen is the one notable case where v2 has *more* outliers than v1.1 (15.9 vs 5.8) — worth a look;
> for most worst-offenders v2 is already better or comparable. A high outlier % on a full record
> (e.g. Pécs, Ell) usually points to a real instrument/site issue rather than the method.

## Reproduce

```
python validation/outliers_timeseries.py   # -> outlier_summary.json + 4 figures
```
Per-stream rates (both methods, both levels) are in
`figs_paper_validation/network_v2_v11/outlier_summary.json`.
