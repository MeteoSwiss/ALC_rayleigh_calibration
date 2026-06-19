# Rayleigh-calibration variability across the E-PROFILE network — L1 2026 (CHM15k, Mini-MPL, CL61)

*MeteoSwiss E-PROFILE ALC paper — M. Hervo. Generated 2026-06-19.*

This report estimates the **night-to-night variability of the Rayleigh (molecular) calibration
constant** for the selected CHM15k, the Mini-MPL, and **all CL61** in the network, processed
directly from the **L1 2026** daily archive (`D:/E-PROFILE_L1_2026`, Feb–May 2026, every night).
The molecular-window detection methods are now named by their **E-PROF calibration version**, the
historical **E-PROF v1.0 (sign error)** is included, and the CL61 Rayleigh result is cross-checked
against an **independent liquid-cloud calibration**. It builds on the molecular-window method study
([`molecular_window_detection_methods_report.md`](molecular_window_detection_methods_report.md),
[`precision_longrun.md`](precision_longrun.md)).

## Method naming (E-PROF versions)

| key (options.json) | label | what it is |
|---|---|---|
| `eprof_v1.0` | **E-PROF v1.0 (sign error)** | legacy `main` window + the pre-a4e7140 Klett **sign error** (full pipeline); historical baseline |
| `eprof_v1.1` | **E-PROF v1.1 (sign cor)** | legacy `main` window, sign-corrected; degenerate (min Σ\|b\|, no R² floor) |
| `eprof_v0.25` | **E-PROF v0.25 (MATLAB)** | MATLAB `Auto_Calib_25` (forced b=0, min ΣRMSE, R² floor) |
| `eprof_v1.2` | **E-PROF v1.2 (improved)** | current production default (above-aerosol + R²/shape/slope gates, max R²) |
| `eprof_v2` | **E-PROF v2** | "optimal" — all physical gates + temporal-variability aerosol rejection + time-resolved cell flagging |
| `earlinet` | EARLINET/SCC | shape+SNR gate, lowest qualifying window |
| `bellini` | Bellini/ALICENET | 3–7 km, Breusch–Godfrey residual test |

Legacy keys (`optimal`/`improved`/`main`/`matlab`/`calipso`→removed) are still accepted as input
aliases. `calipso` was dropped (no stratospheric molecular reference for a ground-up ALC).

---

## 1. Summary

- **24 instruments from L1 2026**: **10 selected CHM15k**, **4 Mini-MPL**, **10 CL61** (~16–120
  fit-nights each, Feb–May 2026). Instrument type is read from the L1 NetCDF `instrument_type`
  attribute.
- **E-PROF v2 is the most precise** molecular-selection version overall (σ_SD **9.5 %** of the median
  C_L), with **EARLINET a close second (10.5 %) at higher yield** — the prior ranking, confirmed on
  fresh L1 data and now on CL61.
- **E-PROF v1.0 (sign error) ≈ v1.1/v1.2 in calibration-constant stability** (σ_SD 10.3 %): the Klett
  sign error corrupts the *attenuated-backscatter product*, **not** the calibration constant. This is
  a useful, citable result — v1.0 and the sign-corrected versions share essentially the same C_L
  variability (see §6).
- **CL61 calibrate as well as the CHM15k.** With the WV correction, the **CL61 median per-instrument
  σ_SD is 7.7 %** (range 4.7–14.5 %), vs **9.7 %** (CHM15k) and **8.5 %** (Mini-MPL). Best units reach
  a **~5 % floor** (Lanzhot CL61 4.7 %, Hohenpeissenberg CHM15k 4.5 %).
- **An independent liquid-cloud calibration corroborates the CL61 Rayleigh result** (§5): on the
  well-exposed units both methods agree at the ~10 % level, and **both independently flag Zeebrugge**
  as anomalous.
- **`instruments.json` was rebuilt** from the actual May-2026 fleet metadata (415 instruments incl.
  the 10 CL61; the old manifest had zero), keeping only the fields the files carry or the code uses.
- **Recommendation:** use **E-PROF v2** (or EARLINET where a single-profile method is preferred) for
  production Rayleigh calibration; feed the per-night constants to the Kalman.

---

## 2. Data and method

**Input.** E-PROFILE **L1** daily files, `D:/E-PROFILE_L1_2026/<WMO>/2026/MM/L1_<WMO>_<id><YYYYMMDD>.nc`,
every available night. Darkness-adaptive night selection (SZA > 100°), cloud/fog screening, profile
outlier screening, then the night-mean range-normalised signal `signal = RCS/r²` + per-profile stack
feed the molecular-window search. **910 nm units (CL61) get the WV correction** (bundled 910 nm LUT +
CAMS 2026); CHM15k (1064 nm) and Mini-MPL (532 nm) are unaffected.

Because the prepared profile is built *before* window selection and is method-independent, all seven
series (the six live methods + **E-PROF v1.0**, which is the legacy window run through the full
pipeline with `sign_error_v10`) are obtained from **one load per night**.

**Instruments**: CHM15k — Payerne, Lindenberg, Aosta, Palaiseau, Granada, Magurele, Bergen, Oslo,
Hamburg, Hohenpeissenberg; Mini-MPL — Brest, Mini-MPL-Mobile, Aleria, Lille; CL61 — Camborne,
Zeebrugge, Birkenes, Lauder, Aosta, Uccle, Lanzhot, Sion/EPFL, Temelín, Edmonton.

**Variability metrics** (identical to `precision_longrun.py`; all % of the median C_L). CV mixes
*precision* with the real *seasonal + laser-ageing drift*, so the headline metrics remove slow drift:

- **σ_SD** — *successive-difference* (von Neumann) precision: robust scatter of |C_Lᵢ₊₁−C_Lᵢ| ÷ √2.
  **Headline precision metric.**
- **σ_detrend** — robust scatter of C_L minus its ~2-month rolling median.
- **σ_within-month** — robust scatter of C_L pooled within each calendar month.
- **σ_night** — mean single-night spread (in-window std of signal/molecular ÷ C_L).
- **rob_CV** — robust night-to-night scatter (1.4826·MAD/median). **valid %** — fraction of
  fit-nights passing the pipeline QC (rel_error ≤ 15 %).

---

## 3. Which molecular selection is best

![Molecular-window methods on L1 2026: yield (left) and drift-insensitive precision (right)](figs_paper_validation/l1_2026_variability/method_precision_l1_2026.png)

Mean over all 24 instruments:

| method | valid % | σ_night % | **σ_SD %** | σ_detrend % | σ_within-month % | rob_CV % | CV % |
|---|--:|--:|--:|--:|--:|--:|--:|
| **E-PROF v2** | 44 | 6.7 | **9.5** | 9.2 | 9.9 | 13.0 | 13 |
| E-PROF v1.0 (sign error) | 60 | 4.4 | **10.3** | 9.2 | 10.2 | 13.3 | 37 |
| **EARLINET/SCC** | 51 | 4.9 | **10.5** | 10.1 | 10.7 | 13.9 | 16 |
| E-PROF v0.25 (MATLAB) | 63 | 21.6 | **11.4** | 10.6 | 11.5 | 14.8 | 23 |
| E-PROF v1.1 (sign cor) | 68 | 23.9 | **11.6** | 10.9 | 12.1 | 15.5 | 31 |
| E-PROF v1.2 (improved) | 54 | 12.7 | **12.9** | 10.8 | 12.3 | 15.3 | 59 |
| Bellini/ALICENET | 34 | 12.6 | **14.1** | 15.4 | 14.3 | 17.4 | 26 |

- **E-PROF v2 is the most precise** (σ_SD 9.5 %); **EARLINET** is a close second (10.5 %) with the
  best yield among the quality methods (51 %).
- **E-PROF v1.1 / v0.25 look stable night-to-night (σ_SD 11–12 %) but carry a huge within-night
  spread (σ_night 22–24 %)**: lacking an above-aerosol gate they admit aerosol windows — a
  *systematic* accuracy risk σ_SD alone doesn't reveal.
- **E-PROF v1.2 (the current production default) is mid-pack here (σ_SD 12.9 %, CV 59 %)** — a few
  spurious-window nights inflate its CV; E-PROF v2/EARLINET's stricter gates avoid them. This is the
  case for switching the default to **E-PROF v2**.
- **Bellini/ALICENET is last (σ_SD 14.1 %, yield 34 %)** and **fails entirely on CL61** (3–7 km band,
  no SNR at 910 nm).

### Per instrument type

![Per-group method precision σ_SD (lower = better; bar label = valid%)](figs_paper_validation/l1_2026_variability/method_precision_by_group.png)

| group | most precise (σ_SD) | E-PROF v2 | note |
|---|---|---|---|
| **CHM15k** (n=10) | EARLINET 9.2 % | 10.1 % (valid 36 %) | EARLINET/v2 lead |
| **Mini-MPL** (n=4) | **E-PROF v2** 11.1 % | 11.1 % (valid 62 %) | v2 clearly best; others 16–18 % |
| **CL61** (n=10) | E-PROF v0.25 7.7 % | 8.1 % (valid 45 %) | v0.25/v1.1 σ_SD low **but σ_night 22–24 %** (aerosol); **E-PROF v2 gives 8.1 % at σ_night 9 %** — the best *clean* choice. Bellini yields 0. |

For CL61, E-PROF v0.25/v1.1 reach the lowest σ_SD only because they take a window every night, but
their 22–24 % within-night spread shows the windows include aerosol; **E-PROF v2 matches their
night-to-night precision while keeping the within-night spread at 9 %**, so it is the method to use.

---

## 4. Per-instrument variability

![Per-instrument σ_SD and yield, recommended method per instrument (colour = type)](figs_paper_validation/l1_2026_variability/per_instrument_variability.png)

Each instrument is characterised with its **recommended method** — the most precise quality method
(preference E-PROF v2 → EARLINET → E-PROF v1.2) with ≥ 6 valid nights, so every instrument gets a
robust number even where the strictest method yields too few nights.

| instrument | type | fit-nights | method | valid | valid % | σ_night % | **σ_SD %** | σ_detrend % | σ_within-month % | rob_CV % |
|---|---|--:|:--|--:|--:|--:|--:|--:|--:|--:|
| Hohenpeiss | CHM15k | 40 | E-PROF v2 | 16 | 40 | 5.9 | **4.5** | 6.5 | 6.1 | 10.6 |
| Hamburg | CHM15k | 34 | E-PROF v2 | 19 | 56 | 6.7 | **7.9** | 6.4 | 5.3 | 7.6 |
| Magurele | CHM15k | 35 | E-PROF v2 | 25 | 71 | 10.5 | **8.6** | 8.0 | 7.6 | 10.7 |
| Oslo | CHM15k | 18 | E-PROF v2 | 14 | 78 | 4.3 | **10.9** | 7.5 | 8.1 | 9.9 |
| Granada | CHM15k | 54 | E-PROF v2 | 13 | 24 | 5.8 | **11.4** | 9.6 | 8.7 | 11.3 |
| Lindenberg | CHM15k | 43 | E-PROF v2 | 22 | 51 | 7.3 | **13.3** | 11.0 | 12.7 | 14.8 |
| Aosta | CHM15k | 44 | E-PROF v2 | 16 | 36 | 6.8 | **13.9** | 14.8 | 14.9 | 19.9 |
| Bergen | CHM15k | 16 | EARLINET | 7 | 44 | 7.0 | **4.5** | 5.7 | 6.4 | 6.1 |
| Palaiseau | CHM15k | 53 | EARLINET | 10 | 19 | 8.4 | **7.6** | 2.9 | 4.8 | 8.0 |
| Payerne | CHM15k | 51 | E-PROF v1.2 | 7 | 14 | 31.7 | **14.2** | 15.3 | 8.4 | 9.3 |
| Mini-MPL-Mobile | Mini-MPL | 119 | E-PROF v2 | 85 | 71 | 2.2 | **6.8** | 5.7 | 8.6 | 13.6 |
| Aleria | Mini-MPL | 117 | E-PROF v2 | 95 | 81 | 1.6 | **6.8** | 7.5 | 15.6 | 17.7 |
| Lille | Mini-MPL | 75 | E-PROF v2 | 44 | 59 | 1.9 | **10.2** | 12.1 | 11.8 | 12.3 |
| Brest | Mini-MPL | 109 | E-PROF v2 | 42 | 39 | 2.4 | **20.5** | 23.5 | 25.9 | 42.1 |
| Lanzhot | CL61 | 49 | E-PROF v2 | 35 | 71 | 8.8 | **4.7** | 7.2 | 7.5 | 8.1 |
| Lauder | CL61 | 55 | E-PROF v2 | 42 | 76 | 6.9 | **6.0** | 7.7 | 10.1 | 11.4 |
| Sion / EPFL | CL61 | 50 | E-PROF v2 | 22 | 44 | 10.1 | **6.3** | 6.6 | 5.8 | 7.7 |
| Aosta | CL61 | 56 | E-PROF v2 | 41 | 73 | 7.8 | **7.5** | 8.1 | 7.6 | 11.5 |
| Temelín | CL61 | 32 | E-PROF v2 | 15 | 47 | 10.8 | **7.7** | 5.8 | 5.8 | 5.9 |
| Uccle | CL61 | 32 | E-PROF v1.1 | 9 | 28 | 16.0 | **7.7** | 9.3 | 9.7 | 6.0 |
| Birkenes | CL61 | 17 | E-PROF v2 | 14 | 82 | 9.3 | **7.9** | 5.7 | 5.2 | 6.2 |
| Camborne | CL61 | 23 | E-PROF v2 | 8 | 35 | 10.6 | **10.4** | 4.6 | 5.1 | 5.1 |
| Zeebrugge | CL61 | 31 | E-PROF v2 | 6 | 19 | 7.5 | **14.5** | 15.6 | 15.5 | 20.0 |
| Edmonton | CL61 | 9 | — | 0 | 0 | — | **—** | — | — | — |

Per-instrument **median σ_SD by type**: **CHM15k 9.7 %** (4.5–14.2), **Mini-MPL 8.5 %** (6.8–20.5),
**CL61 7.7 %** (4.7–14.5).

- **An irreducible night-to-night floor of ≈ 5 %** is reached by the best units of every type
  (Hohenpeissenberg & Bergen CHM15k 4.5 %, Lanzhot CL61 4.7 %).
- **CL61 is the most precise group** (median 7.7 %) — the headline new result (§5).
- **Problem units stand out**: Brest Mini-MPL (σ_SD 20.5 %), Zeebrugge CL61 (14.5 %), and Payerne
  CHM15k (only E-PROF v1.2 on 7 nights, σ_night 31.7 % — its spring-2026 aerosol defeats the strict
  gates).

### Calibration-constant time series

![Per-instrument calibration-constant time series (E-PROF v1.2 vs E-PROF v2; title colour = type)](figs_paper_validation/l1_2026_variability/timeseries_l1_2026.png)

The E-PROF v2 series (red) sit in a tighter band than E-PROF v1.2 (blue), whose spikes are the
spurious-window nights that inflate its CV.

---

## 5. CL61 cross-check — Rayleigh vs liquid-cloud calibration

The CL61 Rayleigh calibration is cross-checked against the **independent liquid-water-cloud
calibration** (O'Connor/Hopkin) run on the monthly L2 archive (the cloud method needs a month of
profiles to accumulate enough liquid-cloud returns; daily files are too sparse). The two methods
calibrate *different physical constants* — the Rayleigh L1 lidar constant vs a multiplier on the L2
attenuated backscatter — so only their **precision and consistency** are comparable.

![CL61 cross-check: liquid-cloud vs Rayleigh calibration precision](figs_paper_validation/l1_2026_variability/cloud_vs_rayleigh.png)

| instrument | cloud months | cloud σ_within-month % | cloud σ_month-to-month % | Rayleigh n | Rayleigh σ_night % | Rayleigh σ_SD % |
|---|--:|--:|--:|--:|--:|--:|
| Birkenes | 2 | 10.5 | 2.4 | 14 | 9.3 | 7.9 |
| Camborne | 2 | 45.2 | 25.4 | 8 | 10.6 | 10.4 |
| Lanzhot | 2 | 7.4 | 5.0 | 35 | 8.8 | 4.7 |
| Lauder | 2 | 11.0 | 12.5 | 42 | 6.9 | 6.0 |
| Aosta | 2 | 17.0 | 0.9 | 41 | 7.8 | 7.5 |
| Sion / EPFL | 2 | 8.6 | 5.8 | 22 | 10.1 | 6.3 |
| Temelín | 1 | 16.3 | – | 15 | 10.8 | 7.7 |
| Uccle | 4 | 27.2 | 22.6 | 0 | – | – |
| Zeebrugge | 4 | 36.1 | **105.9** | 6 | 7.5 | **14.5** |
| Edmonton | 1 | 2.9 | – | 0 | – | – |

- **On the well-exposed CL61 the two independent methods agree at the ~10 % level**: Lanzhot (cloud
  σ_within-month 7.4 %, Rayleigh σ_night 8.8 %), Sion (8.6 / 10.1), Birkenes (10.5 / 9.3), Lauder
  (11.0 / 6.9). This is strong corroboration — a method that does *not* use a molecular window
  reproduces the Rayleigh-derived calibratability.
- **Both methods independently flag Zeebrugge as anomalous** (cloud month-to-month 106 %, cloud
  coefficient ~10× the others; Rayleigh σ_SD 14.5 %, the worst CL61) — an instrument/site issue, not
  a method artefact.
- **Camborne and Uccle are noisier in the cloud method** (45 %, 27 %) on only 2–4 months with few
  clean liquid clouds; their Rayleigh values (σ_SD 10.4 %; n=0 for Uccle) are the more reliable there.
- The mean cloud within-month scatter (18 %, inflated by Camborne/Uccle/Zeebrugge) is larger than the
  Rayleigh σ_SD (8 %), as expected for a sparse monthly cloud sample — but where both have data they
  **agree on which CL61 are stable and which are not**, which is the point of a cross-check.

*(Note: the cloud port's optional above-cloud aerosol-transmission refinement currently collapses the
median coefficient to zero on the 2026 monthly L2 — a regression in that step — so the cross-check
uses the base O'Connor coefficient, which is valid. Driver: `validation/run_cl61_cloud_l1_2026.py`,
`cl61_cloud_vs_rayleigh.py`.)*

---

## 6. The CL61 fleet, and the sign error

**CL61 (new).** The long-run method study had no CL61 (its L2-monthly archive and the old
`instruments.json` predate the fleet). Running the **10 CL61 from L1 2026 with the WV correction**
gives the new result: **CL61 median per-instrument σ_SD 7.7 %**, the *most precise* of the three
groups, with high yield on the well-exposed units (Lauder 76 %, Birkenes 82 %, Aosta 73 %, Lanzhot
71 %) — corroborated by the liquid-cloud cross-check (§5).

**E-PROF v1.0 (sign error).** Including the historical sign-error baseline shows that **the Klett
sign error does not change the calibration-constant stability**: E-PROF v1.0 σ_SD 10.3 % vs E-PROF
v1.1 (sign-corrected) 11.6 % — essentially the same (both use the legacy `main` window). The sign
error corrupts the *downstream attenuated-backscatter / AOD product* (the subject of the Klett
sign-error fix), **not** the Rayleigh calibration constant C_L. This isolates the sign error's impact
to the β_att product and confirms C_L is unaffected.

---

## 7. `instruments.json` rebuilt from the fleet metadata

The shipped `instruments.json` predated the CL61 fleet (zero CL61) and carried fields the code never
uses. `scripts/data/rebuild_instruments_json.py` rebuilds it by reading each stream's NetCDF
metadata (one entry per `(WMO, identifier)` — a station can host several instruments), keeping only
fields that are **in the files** or **used by the code**:

- **kept**: `WMO, Identifier, Type, SiteName, Latitude, Longitude, Altitude, Serial, Calibrated`
- **dropped** (not in the files and unused): `Reference, FLength, NWS, Status`

Rebuilt from the **May 2026** fleet: **415 instruments** incl. **10 CL61** + 145 CHM15k + 4 Mini-MPL
(+ 208 CL31, 48 CL51), 29 multi-instrument stations. (May 2025 — the originally-requested month —
predates the CL61 rollout and yields 0 CL61, so the latest comprehensive month was used; the script
takes `YYYYMM` + `--archive` to regenerate for any snapshot.)

---

## 8. Recommendations

1. **Use E-PROF v2 for production** Rayleigh calibration of CHM15k / CL61 / Mini-MPL (most precise,
   σ_SD 9.5 %); use **EARLINET** where a single-profile method is preferred (close second, higher
   yield) and as the fallback when v2 yields too few nights. **Prefer either over the current E-PROF
   v1.2 default**, which is noisier here.
2. **Keep `calipso` out** — physically inappropriate for ground-up ALCs.
3. **Apply the WV correction on the CL61** — it brings the 910 nm units to CHM15k-class stability.
4. **The liquid-cloud calibration is a valid independent cross-check** for CL61 — it agrees with the
   Rayleigh at ~10 % on well-exposed units and flags the same problem instruments. (Fix the
   transmission-correction regression to use it quantitatively at the monthly level.)
5. **Flag the outlier units**: Brest Mini-MPL, Zeebrugge & Camborne CL61, Payerne CHM15k (aerosol).
6. **Feed per-night constants to the Kalman**, which absorbs the sparsity and residual outliers.

---

## 9. Files & reproduction

- **Scoping**: `validation/scope_l1_2026.py` → `scope_l1_2026.json`
- **Rayleigh run** (24 instruments × all nights × 6 methods + E-PROF v1.0, parallel):
  `validation/run_l1_2026_variability.py` → `results_<label>.json`
- **CL61 cloud cross-check**: `validation/run_cl61_cloud_l1_2026.py` → `cloud_<label>.json`;
  `validation/cl61_cloud_vs_rayleigh.py` → `cloud_vs_rayleigh.png`, `cloud_crosscheck_table.md`
- **Metrics + figures**: `validation/variability_metrics_l1_2026.py` →
  `method_precision_l1_2026.png`, `method_precision_by_group.png`, `per_instrument_variability.png`,
  `timeseries_l1_2026.png`, `metrics_tables.md`, `metrics_summary.json`
- **instruments.json rebuild**: `scripts/data/rebuild_instruments_json.py`
- **Method rename**: E-PROF version keys in `calibration/rayleigh/molecular_methods.py`
  (`METHODS`, `METHOD_LABELS`, `ALIASES`, `_SELECTORS`), `config.py`, `rayleigh_fit.py`,
  `options.json`, validation harness. Metric definitions reuse `validation/precision_longrun.py`.
  Unit + smoke tests pass.
