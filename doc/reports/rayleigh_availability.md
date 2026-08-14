# CHM15k Rayleigh calibration — availability study

*Branch `rayleigh-availability`, generated 2026-08-14 17:07. Candidate: **`eprof_v2.2` / N2.5**. Every number below is recomputed from the run outputs by `rayleigh_availability/make_report.py`.*

## 1. The gates were measuring instrument age, not atmosphere

v2 rejects CHM15k nights with flag -2 ("no molecular window passed the validity gates") in proportion to how NOISY the instrument is. Four of v2's gates — `residual_pct`, `min_r2`, `ratio_std`, `temporal_cv` — are functions of SNR compared against FIXED thresholds (v1.1 has no such gates at all), so an ageing laser fails them on perfectly clean nights.

![Rejection rate and availability versus measured night noise](figs_rayleigh_availability/phase0_noise_vs_availability.png)

| evidence | statistic |
|---|---|
| Network, 141 CHM15k streams: flag -2 rate vs measured night noise | Spearman **+0.75** |
| Corpus, per-stream v2-minus-v1.1 availability vs noise | Spearman **-0.87** |

The figure also carries the candidate, because that is the test the fix has to pass: a gate change that only worked on quiet instruments would be no fix at all. **v2.2 is flat across the whole noise range** (right panel, green) where v2 falls from 84 % to 3 % of clear nights, and the flag -2 collapse is largest exactly on the worst streams (left panel: Vasarosnameny 92 -> 14 %, Gottfrieding 90 -> 10 %, Payerne 77 -> 5 %). The availability deficit is no longer a function of instrument condition.

v2's advantage over v1.1 reverses sign with instrument condition: **+13.2 pts** on the 9 quiet instruments (sigma < 0.06) versus **-16.3 pts** on the 7 noisy ones (sigma >= 0.10). v2 is genuinely better on healthy hardware and progressively worse as the laser ages.

The failure is also strongly seasonal — a regime the C8 gate tuning (Feb-May 2026 only) never contained. Pooled over the corpus CHM15k streams, per CLEAR night (the denominator the gates actually act on):

| season | clear nights | v1.1 avail | v2 avail | v2.2 avail | v2 flag -2 | v2.2 flag -2 |
|---|---|---|---|---|---|---|
| winter | 1177 | 67.2 % | 64.4 % | **82.9 %** | 29.9 % | 5.1 % |
| shoulder | 2278 | 56.7 % | 57.2 % | **84.3 %** | 39.6 % | 5.6 % |
| summer | 1240 | 32.5 % | 34.0 % | **69.1 %** | 63.5 % | 21.0 % |

*(Winter = Nov-Feb. An earlier version of this corpus stopped at 30 September, so winter was represented by Jan/Feb alone; the table above uses the full year.)*

## 2. What `eprof_v2.2` changes

The three SNR-driven gates are replaced by their **noise-relative** forms, using a per-night photon-noise profile measured from consecutive-profile differences on the NATIVE resolution (before the L2-grid binning averages it away). A window is kept when its departure from a Rayleigh shape is no larger than its OWN measured noise explains (reduced chi-square), rather than smaller than a fixed percentage.

The scattering-ratio reference is also de-biased: it was the MINIMUM over ~136 candidate windows, and a minimum-of-N is biased low by an amount that grows with noise — so a noisy night mechanically inflated every window's ratio. It is now a low percentile over windows that are individually noise-consistent. **The 1.15 threshold itself is unchanged.**

`eprof_v2.2` is a strict **superset** of `eprof_v2`: the noise-relative tier is consulted only when the strict v2 gates leave no eligible window. Every night v2 calibrates is calibrated identically — same window, same constant — so the change can only ADD nights.

**Ablation.** Re-running with the ORIGINAL min-based scattering reference isolates what each change buys (22 CHM15k streams, 4695 clear nights):

| variant | valid nights | availability |
|---|---|---|
| v2 | 2481 | 52.8 % |
| v2.2, full | 3738 | **79.6 %** |
| v2.2, raw (min-based) scattering reference | 3662 | 78.0 % |

The de-biased reference accounts for only **6 %** of the recovery (76 of 1257 nights). The effect is real but small: essentially all of the gain comes from replacing the three direct SNR gates with their noise-relative forms. The de-biasing can therefore be dropped if simplicity is preferred — it costs a second scattering array and two extra parameters.

## 3. Availability

Split by the station's C_L-vs-window-height gradient (section 4 explains why this matters). Holdout = unseen streams and 2026. **Caveat:** the band assignment itself is unreliable when taken from v2's own series -- see section 4.1 -- so this table indicates the shape of the effect, not a deployable rule.

| band (holdout) | n | availability | sigma_SD | Kalman outliers /100 |
|---|---|---|---|---|
| low-gradient (<= 8 %/km) | 12 | 49.2 -> **86.2 %** | 9.5 -> 11.2 % | 6.12 -> 6.46 |
| high-gradient (> 8 %/km) | 9 | 79.6 -> **87.8 %** | 9.2 -> 9.2 % | 10.53 -> 9.33 |

![Calibration time series, Gottfrieding](figs_rayleigh_availability/timeseries_GOTTFRIEDING.png)

Gottfrieding is the case that matters most: v2 produced **19 nights in 18 months**, so the operational Kalman was effectively interpolating a constant across month-long gaps. v2.2 adds ~115 nights, the recovered points sit ON the existing Kalman line, and sigma_SD *improves* by 40 %.

![Calibration time series, Amsterdam](figs_rayleigh_availability/timeseries_AMSTERDAM.png)

Amsterdam (already densely sampled, near-zero gradient): the recovered nights are indistinguishable from the retained ones.

## 4. Validation against independent references

Internal guards cannot see a *coherent* bias: sigma_SD, continuity and the Kalman-outlier count are all difference statistics, so a set of nights that is uniformly 25 % low but internally consistent passes all of them. Only a co-located independent instrument can detect it — which is why this section, not section 3, is the decisive one.

![Amsterdam inter-unit consistency](figs_rayleigh_availability/phase2_amsterdam_interunit.png)

**Amsterdam quad** (four co-located CHM15k — same night, same atmosphere, no wavelength or water-vapour correction involved): the pair-ratio scatter on recovered nights is **1.13x** that of the nights v2 already had (three pairs better, three worse), with no systematic offset. At a near-zero-gradient site the recovery is sound.

**CHM15k vs co-located CL61** (CL61 calibrated by the independent liquid-cloud route). The offset of the recovered nights tracks the station's altitude gradient in sign and approximate magnitude — a directional, falsifiable prediction that held at two sites with gradients of OPPOSITE sign:

| site | gradient | predicted offset | observed vs CL61 |
|---|---|---|---|
| Payerne | **-14.3 %/km** | -20.6 % | **-26.5 %** |
| Lindenberg | **+7.8 %/km** | +14.9 % | **+27.9 %** |

![Calibration time series, Payerne](figs_rayleigh_availability/timeseries_PAYERNE.png)

Payerne shows the failure directly: v2.2 fills an 8-month hole in which the Kalman had been holding a flat constant, but the recovered points sit visibly below the retained ones. **A true lidar constant cannot depend on the fit altitude**, so a non-zero gradient is an unmodelled profile defect (background subtraction, overlap residual, or aerosol the molecular model does not capture) — and on nights v2 rejects, the only eligible windows are higher up, so any recovery inherits gradient x height-shift.

### 4.1 The gradient cannot be used as a pre-flight gate

The obvious safeguard — measure each station's gradient from its EXISTING v2 series and only enable the recovery where it is small — does not work, for a structural reason. v2's accepted nights are exactly the clean LOW windows, so they span too little height to fit a slope against height. Bootstrapping the gradient (90 % CI) from each source:

| station | from v2 nights | from v2.2 nights |
|---|---|---|
| Payerne | 43 n, 1439 m span, **+1.7 [-5.6, +8.4]** | 150 n, 2398 m, **-14.7 [-18.8, -10.9]** |
| Gottfrieding | 20 n, 959 m, -1.0 [-26.5, +12.6] | 140 n, 2158 m, **-5.7 [-10.9, -1.2]** |
| Vasarosnamaony | 25 n, 719 m, +5.6 [-14.0, +24.5] | 192 n, 2158 m, **+11.4 [+5.8, +17.4]** |
| Bonaire | 7 n, not computable | 33 n, 1439 m, **-12.1 [-22.0, -3.5]** |
| Guadiana | 65 n, 1439 m, +43.2 [+30.7, +56.0] | 189 n, 2398 m, **+6.8 [+1.4, +12.9]** |

Only 6 of 22 streams get a stable band assignment from the v2 series, and Payerne's reads as consistent with ZERO (+1.7 +/- 7) when the truth from the wider span is -14.7 +/- 4. Guadiana is worse than imprecise: the narrow-span estimate (+43) is simply wrong (+6.8). Where v2 does happen to span a wide height range (Amsterdam, Magurele, Aosta) the two sources agree — confirming it is the SPAN, not the method, that governs.

So the gradient is only observable once the high-altitude windows are taken, which is the very thing it was meant to gate. It has to be a second pass, not a pre-flight check.

## 5. What the recovered nights actually are

Phase 2 left the recovered nights described only by what the gates did to them. This section characterises them from the data, using the quantity every method fits and **no window selection at all**:

> `R(z) = signal(z) / p_mol(z)`, which equals C_L wherever the atmosphere is purely molecular. Normalised per night at 3 km, then combined over nights.

![Signal over molecular, kept versus recovered nights](figs_rayleigh_availability/ratio_profile.png)

| stream | slope of R(z), 2-6 km — kept by v2 | recovered by v2.2 |
|---|---|---|
| PAYERNE_CHM15k_A | -3.7 %/km | **-19.5 %/km** |
| LINDENBERG_CHM15k_0 | -3.5 %/km | **-18.6 %/km** |
| GOTTFRIEDING_CHM15k_0 | +4.5 %/km | **-18.5 %/km** |

On the nights v2 already had, R(z) is **flat** through the fit band. On the nights v2.2 recovers it **falls by about 19 %/km, at every site tested** — and it is the same ~19 %/km whether the station's dC_L/dz is negative (Payerne), positive (Lindenberg) or small (Gottfrieding). That is a property of the NIGHTS, not of the station.

A decline of that size is what two-way aerosol extinction looks like: `R ∝ exp(-2∫α dz)`, so 19 %/km implies α ≈ 0.1 /km — a real, ordinary aerosol column. So these are **not clean nights that noise falsely rejected**. v2's verdict on them ("signal not proportional to molecular") is literally true; what v2 lacks is the ability to find the sub-window where the departure is no larger than the night's own noise explains. The pipeline already knows they are weaker evidence — their median reported uncertainty is **12.4 %** against **5.0 %** for the retained nights.

### 5.1 ...but the aerosol load does not PREDICT the constant

If the offset of the recovered nights were simply the two-way transmission of the aerosol below the window, it would scale with that column and be correctable per night. Testing it directly — each night's R(z) slope (the aerosol proxy) against its constant relative to the stream level — the correlation is inconsistent in SIGN across sites:

| stream | r (aerosol proxy vs C_L) | median slope, kept | recovered |
|---|---|---|---|
| PAYERNE_CHM15k_A | **+0.33** | -4.2 %/km | -20.9 %/km |
| LINDENBERG_CHM15k_0 | **-0.26** | -3.8 %/km | -29.9 %/km |
| GOTTFRIEDING_CHM15k_0 | **-0.01** | +2.0 %/km | -13.9 %/km |

![Aerosol load below the window versus the retrieved constant](figs_rayleigh_availability/aerosol_load_vs_cl.png)

Payerne goes the way transmission predicts (more aerosol -> lower C_L), Lindenberg goes the opposite way and Gottfrieding not at all. So the aerosol column reliably IDENTIFIES the recovered nights but does not PREDICT their constant, and a per-night transmission correction is not supported by this evidence. The two-pass, per-station altitude route of §8 remains the one with evidence behind it.

### 5.2 The classification mask does not address this

The plan's Phase 3 assumed steady elevated aerosol layers were forcing the fit high, and that masking them with the Cloudnet target classification would bring the windows back down. Measured on the classified corpus, it does not:

*(measured by `rayleigh_availability/cl61_crossmask.py`, which also repeats the test with the co-located CL61 as the screening instrument.)*

| nights | own classifier, 2-4 km | co-located CL61, 2-4 km |
|---|---|---|
| Payerne, kept by v2 | 0.0 % | 0.0 % |
| Payerne, recovered by v2.2 | **0.0 %** | **0.2 %** |
| Payerne, still rejected | 0.1 % | 9.6 % |
| Lindenberg, kept by v2 | 0.0 % | 0.0 % |
| Lindenberg, recovered by v2.2 | **0.0 %** | **2.8 %** |
| Lindenberg, still rejected | 0.0 % | 7.9 % |

The classifier is working — it flags the nights that stay rejected — but on the recovered nights there is nothing in the relevant band for it to mask. The obvious objection is that a single-channel CHM15k classifier could be BLIND to a layer that still biases a fit, so the test was repeated with the co-located CL61, which has depolarisation, as the screening instrument. It is indeed more sensitive — it finds 8-10 % on the nights that STAY rejected, where the CHM15k's own classifier finds ~0 % — but on the RECOVERED nights it finds 0.2 % and 2.8 %. The CHM15k classifier is not blind; there is no classifiable layer there. That is exactly what §5 predicts: distributed extinction, not a layer.

Run over the classified streams, the mask changes little and moves nothing it should not:

| gates | valid nights without -> with mask | lost | gained | kept nights moved >1 % | flag -11 |
|---|---|---|---|---|---|
| operational gates | 1378 -> 1262 | 117 | 1 | 1 / 1261 | 0 |
| noise-aware gates | 1850 -> 1756 | 98 | 4 | 0 / 1096 | 0 |

The mask is implemented and tested, and is a genuine improvement to the pre-fit cleaning in general (it is unioned with the temporal MAD screen, which by construction cannot see anything that persists all night). It is simply not the lever for this problem.

### 5.3 Flag -11 repaired

Independently of the mask, the classification veto was broken. Its contaminated fraction was computed over the **whole 48 h** spanned by the two classification files, while a Rayleigh night is 6-10 h — so a cloud filling the fit window for an entire night scored ~15 % against a 30 % threshold. The flag had never fired network-wide. It is now computed over the night the fit actually used, and the threshold is a real setting (`classification_veto_fraction`) rather than a literal.

## 6. Outliers among the recovered nights

Measured against the **local** level (median of valid nights within +/-30 days), so a genuine hardware step is not mistaken for an outlier — a station-wide median would flag a whole pre-step era, which is exactly what happens at Guadiana (its constant steps ~2.5x at the 2026 boundary; see figure).

| night origin | outside 0.6-1.67x local | rate |
|---|---|---|
| kept by v2 | 19 / 2090 | **0.9 %** |
| recovered by v2.2 | 60 / 1152 | **5.2 %** |

![Calibration time series, Guadiana](figs_rayleigh_availability/timeseries_GUADIANA.png)

The recovered nights carry a materially higher local-outlier rate. The tail is **one-sided (low)**, its windows sit high, and it concentrates in the noisiest stations — i.e. it is the extreme end of the same altitude mechanism, not a separate failure mode.

## 7. The best estimate: the added nights must not weigh the same

The operational Kalman (`improve_alc_calib/cal_best_estimate.py`, vendored in `monitoring/kalman.py`) uses **one scalar measurement variance for the whole series** — `(res**2).mean()` of the fit residuals, reused for every update. A night the pipeline reports at 25 % uncertainty therefore pulls the best estimate exactly as hard as one at 5 %. That is defensible while the gates admit only the cleanest nights; §5 shows it stops being defensible once availability is raised, because the added nights are systematically the less certain ones. The reference author foresaw this — `cal_best_estimate.py` carries the TODOs *"check whether lidar_constant_uncertainty can be used instead of residuals"* and *"launch kalman — configure with uncertainties"*.

`kalman_best_estimate(..., uncertainties=...)` implements it: the empirical variance scale is kept (it captures error the formal uncertainty does not) and only the RATIO between nights is taken from the reported values, clipped to a factor 4. Over 21 CHM15k streams (2474 -> 3252 nights):

| | median movement of the best estimate from v2, on days v2 already covered | p90 |
|---|---|---|
| flat (operational) | 1.86 % | 13.33 % |
| uncertainty-weighted | 2.79 % | **11.52 %** |

The medians favour the flat filter and the tail favours the weighted one, and that is the whole point: weighting halves the swing exactly where availability jumps most — Payerne 13.3 -> 9.5 %, Vasarosnameny 19.8 -> 13.2 %, Montsec 7.7 -> 4.2 % — at the cost of 1-2 % on streams that never needed help. It does NOT change the outlier count: the IQR screen runs before the filter.

**Outlier indicator, both definitions.** The operational screen (`flag_outliers_lom`) is a SINGLE IQR over the whole series; the dashboard's vendored copy used a rolling 30-day window. They differ by 2.5x, so the definition has to be stated:

| Kalman outlier screen | v2 | v2.2 | rate |
|---|---|---|---|
| operational (one global IQR) | 73 | 102 | 2.95 % -> **3.14 %** |
| rolling 30-day (dashboard) | 179 | 269 | 7.2 % -> 8.3 % |

Nights rise 31 %; the operational outlier RATE rises 0.2 pp.

## 8. Recommendation

In decreasing order of confidence:

1. **Adopt the noise-aware gates.** They are a strict superset of v2 (nothing moves or is lost), they lift availability in every season — most where it is worst (summer 34 -> 69 % of clear nights) — and on the operational outlier definition the rate rises only 2.95 -> 3.14 % for 31 % more nights.
2. **Weight the Kalman by the per-night uncertainty at the same time.** The two changes belong together: v2.2 admits nights the pipeline itself rates 2.5x less certain, and the operational filter currently ignores that. This is the reference implementation's own unimplemented TODO, and it damps the best-estimate swing precisely on the streams the availability gain is largest for.
3. **Deploy in TWO PASSES, not behind a pre-flight gate.** Run v2.2, then measure each station's C_L-vs-height gradient FROM ITS OWN OUTPUT (where it is tightly constrained, +/-2-5 %/km) and either correct the recovered constants back to the station's reference altitude (`C_corr = C * exp(-grad * dz)`) or flag them. Section 4.1 shows why the intuitive one-pass version is not implementable. This stays safe because v2's own nights are never touched: the worst case is that a station's recovered nights end up flagged rather than corrected.
4. **Apply a local-consistency check to recovered nights only** (reject beyond ~1.6x the local level). Verified not to clip Guadiana's real 2.5x step; restricting it to recovered nights means it can never remove a night from the existing series.
5. **Do not deploy the classification pre-fit mask for this purpose** (§5.1) — it does not touch the problem. **Do** take the flag -11 repair (§5.2), which is an outright bug fix.

**Still open:** the altitude correction in (3) is proposed, not validated. It should be tested exactly as the gates were — against the co-located CL61 pairs, checking that the corrected recovered nights land on the same level as the retained ones at Payerne (-26.5 % today) and Lindenberg (+27.9 % today). The obvious physical shortcut — treating the offset as the two-way transmission of the aerosol column below the window, and correcting it per night from R(z)'s own slope — was tested and is NOT supported (§5.3).

`options.json` still selects `eprof_v2`; nothing is deployed by this branch.

## 9. Limits and caveats

- **The recovered nights are aerosol-loaded, not clean** (§5). Framing v2's gates as "measuring noise, not atmosphere" is right about WHICH INSTRUMENTS reject most (Spearman +0.75 against measured noise) but wrong if read as "the rejected nights were clean". They were not. What v2.2 adds is the ability to fit them where the departure from Rayleigh is within the night's own noise — with a correspondingly larger uncertainty.
- Not every summer rejection is a noise victim: an aerosol-laden night can have *higher* SNR and still be correctly rejected. The recovery is partial by design.
- The chi-square tolerance is **inert** over the range tested (1.5-3.0 give near-identical results, differing on 4 of 36 streams); the active ingredients are the de-biased scattering reference and the noise-excess gates.
- The noise propagation from native to binned grid assumes white noise; correlated components (afterpulse, background drift) would break it, which is why `max_chi2red` is held above 1.
- The clear-night screen (flag -1) is the single largest v1-to-v2 divergence bucket and is deliberately **out of scope** here.
- Costs ~2.5x v2 per night (the native-resolution noise estimate runs on every night). Acceptable for the daily operational run; the estimate could be restricted to the fit band.
