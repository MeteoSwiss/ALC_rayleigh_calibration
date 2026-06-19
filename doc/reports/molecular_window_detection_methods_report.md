# Molecular-window detection for Rayleigh calibration — pluggable methods, comparison & recommendations

*MeteoSwiss E-PROFILE ALC paper — M. Hervo. Generated 2026-06-17.*

Molecular-window detection is the core of elastic-lidar Rayleigh calibration: an elastic
ceilometer (CL31/CL51/CL61/CHM15k) measures only *total* attenuated backscatter, so to get an
absolute constant it must find an altitude window it can assume is **aerosol-free** and force the
range-corrected signal to match the molecular (Rayleigh) profile computed from T/p. The quality of
the whole calibration is hostage to that window choice. This report makes the detector **pluggable**
(seven interchangeable strategies), tests them across six instruments at three sites (Payerne,
Amsterdam, EDT/Edmonton), checks what EarthCARE does, and selects the best configuration.

---

## 1. Summary

- **The molecular detector is now a selectable strategy** (`molecular_method` in `options.json`),
  with six implementations sharing one grid search. Default is `improved` (the production fix),
  so existing behaviour is unchanged unless you opt in.
- **Seven methods:** `main` (legacy), `improved` (production), `matlab` (Auto_Calib_25), `calipso`
  (CALIOP-type), `earlinet` (EARLINET/SCC-type), `optimal` (best-of), `bellini` (ALICENET / Bellini
  et al. 2024, AMT 17, 6119). **`optimal` is the only method that uses the time dimension** — it flags
  and excludes contaminated time-altitude cells (aerosol present only part of the night, a passing
  cloud) and fits on the time-cleaned profile (see the flowchart in §3).
- **Tested across six instruments at three sites** — Payerne (CHM15k, CL31, CL61), Amsterdam (CHM15k),
  EDT/Edmonton (CL51, CL61), ~30 clear nights each — ranked by night-to-night stability, usable-night
  count and cleanliness.
- **Selected best: `improved`** (best balance, already the default → no change), with **`optimal`
  recommended where aerosol/cloud is significant** (equal stability, cleaner, flags contaminated
  layers). `matlab` gets the most nights but admits aerosol (temporal CV 1.7); **`calipso` and `main`
  are unusable** (night-to-night CV 80–200 %).
- **`bellini` (ALICENET)** — its Breusch–Godfrey residual-autocorrelation test + E_CL gate give the
  cleanest fits (lowest rel-error, 2–3 % on CHM15k) and it does well on high-SNR CHM15k (Amsterdam
  20 nights), but its 3–7 km window costs SNR on the 910 nm units (fewer nights, CV ≈ 14 %): mid-pack
  overall, below the `improved`/`earlinet`/`optimal` group.
- **CL31/CL51 Rayleigh is marginal** for every method (no clean molecular window) — confirms they
  should use the liquid-cloud calibration; CHM15k and CL61 calibrate well.
- **EarthCARE/ATLID needs no window search**: it is a High-Spectral-Resolution Lidar with a direct
  molecular channel; it only uses a high-altitude pure-Rayleigh band to calibrate channel cross-talk.

---

## 2. Architecture — how to choose a method

All logic lives in `ALC_rayleigh_calibration/rayleigh_calibration/molecular_methods.py`:

- `compute_window_grid(signal, p_mol, range, …, signal_stack=…)` runs the grid search **once** and
  returns every candidate window's statistics (free-fit slope/intercept/R²; forced-intercept-0
  slope/RMSE/R²; median signal/molecular ratio; in-window SNR; Rayleigh-shape residual; scattering
  ratio; and — if the per-profile `signal_stack` is given — the **temporal** variability).
- Six **selectors** apply different eligibility masks and selection rules to that grid.
- `select_molecular_window(method, …)` is the dispatcher.

**To switch method**, set `"molecular_method"` in `options.json` to one of
`improved | main | matlab | calipso | earlinet | optimal`. The production pipeline
(`calibrate_rayleigh` → `find_optimal_molecular_window`) reads it; `improved` keeps the exact,
already-validated in-line path, the others are dispatched to `molecular_methods`. The
`min_window_*`/`max_window_*` gates in `options.json` tune `improved`; the other methods use their
own documented defaults.

---

## 3. The seven methods

| method | fit | eligibility (gate) | selection rule | basis |
|---|---|---|---|---|
| **main** | free intercept | none (any finite fit) | min Σ\|intercept\| centre, then max R² | the legacy main-branch code — **degenerate** (selects high-altitude noise where signal→0 ⇒ intercept→0 *and* R²→0) |
| **improved** | free intercept | start ≥2 km, R²≥0.5, slope>0, \|b\|<a, slope≈median ratio | **max R²** among eligible | production fix; Mattis 2016 shape gate + above-aerosol |
| **matlab** | **intercept forced = 0** | R²₀≥0.5, slope>0; centres ≤5 km | min Σ RMSE centre, then max R² | Auto_Calib_25/rayleigh_fit.m (Hervo & Poltera 2014) |
| **calipso** | free | SNR gate, scattering ratio ≈ min (least aerosol) | **highest** clean window | CALIOP molecular normalization (R≈1, highest clean layer) |
| **earlinet** | forced 0 | Rayleigh-shape residual ≤10 %, SNR gate, scattering ratio ≤1.1 | **lowest** qualifying window | EARLINET SCC (Mattis 2016: "take the lowest"), Freudenthaler 2018 |
| **optimal** | free | all of the above **+ temporal-variability ≤ threshold** | max composite quality (R² + shape + purity + SNR + **temporal steadiness**) | best-of; the temporal idea is novel here |
| **bellini** | free | 3–7 km, width 600–3000 m; residuals not autocorrelated (Breusch–Godfrey); slope>0, intercept≈0; border-residual sign<0; E_CL≤40% | max **M_Ray = (adjR²+(1−\|b\|))/std(b)** | ALICENET (Bellini et al. 2024, AMT 17, 6119) |

The three gated alternatives (`calipso`, `earlinet`, `optimal`) also require the fit slope to match
the median signal/molecular ratio to within the pipeline's `threshold_quality` (`rel_error ≤ 15 %`),
so a window they select also passes the downstream proportionality QC — keeping selection and
pipeline consistent. (`main` and `matlab` are left faithful to their originals, so they can select a
window the pipeline's QC then rejects — see §6.)

### The `optimal` method's temporal-variability aerosol rejection

All other methods collapse the night to **one mean profile** and must guess whether a smooth,
linear-looking layer is molecular or aerosol. `optimal` additionally uses the **full profile time
series**: molecular scattering is **steady in time** (it only tracks T/p), whereas **aerosol advects
and fluctuates**. For each candidate window it computes the temporal coefficient of variation of the
window-mean signal/molecular ratio across the night (averaging over the window's range bins first
suppresses photon noise, leaving mostly atmospheric variability). Temporally variable windows are
**rejected as aerosol** and steadier windows are rewarded. This catches aerosol that *looks* linear
in the mean (high R²) but betrays itself by fluctuating — information every single-profile method
discards. In the test it gives `optimal` consistently among the **lowest** temporal CVs (see §6).

**Time-resolved layer flagging.** Beyond the per-window temporal CV, `optimal` flags individual
*time-altitude cells* contaminated by aerosol/cloud: per altitude it takes the temporal median and
MAD of signal/molecular and flags cells exceeding `median + 4·MAD` (an upper-tail outlier test, so
clean molecular noise is left intact and the cleaned mean stays unbiased — an earlier percentile
threshold biased it low and had to be replaced). It then fits on the *time-cleaned* mean, using the
clean part of an otherwise-contaminated night (aerosol only at the start, a cloud only at the end).
The flagged cells are drawn hatched on the pcolor figures in §6 (`flag_contaminated_cells`).

**Flowchart of `optimal`** (the red box is the distinctive time-resolved flagging step):

![Flowchart of the optimal molecular-window detection method](figs_paper_validation/molecular_methods/optimal_flowchart.png)

Step by step: (1) take the night's per-profile range-normalized signal `signal(t,z)` and the
molecular profile `p_mol(z)`; (2) **flag and remove** aerosol/cloud cells per altitude
(`signal/p_mol > median + 4·MAD`), which excises an aerosol layer present only part of the night or a
passing cloud without biasing the clean molecular noise; (3) average the un-flagged cells into a
**time-cleaned mean** profile; (4) grid-search windows (centre × half-length) and fit
`signal = a·p_mol + b` in each; (5) keep only windows passing **all** eligibility gates — start above
the boundary-layer aerosol, R²≥0.5, slope>0, |b|<a, Rayleigh-shape residual ≤12 %, scattering ratio
≤1.1, adequate in-window SNR, slope ≈ median ratio (≤15 %), and **temporal CV ≤0.5** (steady = molecular);
(6) if none pass, emit **no calibration** (flag −2) and let the Kalman skip the night; (7) otherwise
choose the window maximising the composite quality `Q = R² − 0.25|R−1| − 0.20·resid − 0.10·SNR −
0.35·tCV − 0.20·rel + 0.10·n`; (8) the lidar constant is `C_L = median(signal/p_mol)` in that window,
passed to the Klett β_att step and the downstream pipeline QC.

### The `bellini` method (ALICENET)

`bellini` implements the ALICENET Rayleigh calibration (Bellini et al. 2024, AMT 17, 6119, Supplement
S3). It searches windows within **3–7 km** (widths 600–3000 m) and, crucially, rejects any window
whose linear-fit **residuals are autocorrelated** — a **Breusch–Godfrey** test (here a lag-1
autocorrelation proxy): coherent residual structure betrays an *undetected aerosol layer* even when
R² looks good. Surviving windows must also have slope>0, intercept≈0, and predominantly negative
border residuals (±200 m); among them the one maximising **M_Ray = (adjR² + (1−|b|))/std(b)** is
selected, then rejected if the relative calibration uncertainty **E_CL = err(slope)/slope +
std(C_L)/median(C_L) > 40 %**. (Its "negative-AOD" QC is enforced by the pipeline's Klett step.) It is
the most quality-stringent method; its higher search band is well suited to high-SNR CHM15k but costs
nights on the weaker 910 nm units (see §7).

---

## 4. Does EarthCARE (ATLID) do molecular detection? — No, and here's why

ATLID is a **355 nm High-Spectral-Resolution Lidar (HSRL)** with three channels (co-polar Mie,
co-polar Rayleigh/molecular, cross-polar). A Fabry-Pérot **high-spectral-resolution etalon (HSRE)**
separates the spectrally **broad molecular** return (→ Rayleigh channel) from the **narrow
particulate** peak (→ Mie channel). Because it **measures the molecular return directly at every
range gate**, ATLID does **not** need to search for an aerosol-free window — the molecular channel
*is* the calibration reference throughout the profile (Donovan et al. 2024; Wehr et al. 2023).

It still uses a **high-altitude (~30–40 km) pure-Rayleigh band**, but for a *different* purpose than
an elastic lidar: to determine the **spectral cross-talk coefficients** (χ = Mie-into-Rayleigh,
ε = Rayleigh-into-Mie; a height-dependent 3×3 channel matrix) and the inter-channel gain — not to
fix an absolute backscatter constant. (The CALIPSO "X-factor" term does not appear in the ATLID
literature; the analogous quantities are χ/ε.) A weekly **fine spectral calibration** keeps the
laser centred on the etalon.

**What an elastic ceilometer can borrow:** the discipline of *defining* clean air (volcanic-aerosol
awareness, averaging to beat noise, treating the clean-air assumption as testable); **per-shot
background by interpolating a pre- and post-echo estimate**; proper **T/p-dependent** molecular term;
and using **surface and cloud echoes as auxiliary references** (parallels our liquid-cloud
calibration). **What cannot transfer:** the direct, continuous molecular reference (needs spectral
separation), the χ/ε cross-talk correction, and extinction without a lidar-ratio assumption — these
are exactly what a single-channel elastic system lacks, which is *why* the molecular-window search
remains necessary. *(Sources: Donovan et al. 2024 AMT 17 5301; Wehr et al. 2023 AMT 16 3581;
Eisinger et al. 2024 AMT 17 839; Irbah et al. 2023 AMT 16 3631.)*

---

## 5. Literature basis for the selection criteria

- **Mattis, D'Amico, Baars et al. 2016** (EARLINET SCC, AMT 9, 3009): the minimum-signal/
  minimum-background search is **degenerate** — it "would find a minimum also in the case that there
  are fewer particles than in other altitude regions only … large errors." Remedy: SNR/std gate +
  Rayleigh-shape test + **take the lowest qualifying window**. → basis of `earlinet`, and the reason
  `main` fails.
- **Freudenthaler et al. 2018** (EARLINET QA): Rayleigh-fit **relative residual ≤ 1 %** + deviation
  plot. → the shape-residual gate.
- **Wiegner & Geiß 2012** (AMT 5, 1953): scattering ratio **R ≤ 1.1** as an acceptance criterion and
  systematic-error term; slope matching. → the scattering-ratio gate.
- **Baars et al. 2016** (PollyNET, ACP 16, 5111): automated Rayleigh-shape test. **CALIOP** (Powell
  et al. 2009): scattering ratio 1.01 ± 0.01, normalize as **high** as possible. → basis of `calipso`.
- **Bellini et al. 2024** (ALICENET, AMT 17, 6119, Suppl. S3): two-step E-PROFILE-based Rayleigh fit
  in 3–7 km with a **Breusch–Godfrey residual-autocorrelation test** to reject undetected aerosol, the
  **M_Ray** window metric, and an **E_CL** relative-uncertainty gate. → basis of `bellini`.

---

## 6. Detailed examples — profiles, selected windows, and time-resolved flagging

Each method runs on the **identical prepared profile** per night (`compare_molecular_methods.py`).
In the profile figures the brackets are each method's selected window; in the pcolor figures the
colour is the signal/molecular ratio over the night (molecular = vertically and temporally uniform;
aerosol/cloud = enhanced and variable), the dashed lines are the selected window centres, and the
**hatched cells are those `optimal` flagged as aerosol/cloud and excluded** from its fit. The
**bottom row of each figure shows the resulting lidar calibration constant C_L ± uncertainty per
method** (only methods that calibrate through the pipeline; same colour code) — so the agreement
(or disagreement) of the constants is read directly beneath the windows that produced them.

**Payerne CL61 (910 nm)** — 2026-03-12 (aerosol-laden), -16 and -28 (cleaner). As you noted, -16 has
aerosol only at the *start* of the night and -28 a thin cloud near the *end* — both appear as hatched
flagged cells:

![Payerne CL61 vertical profiles, methods bracketed](figs_paper_validation/molecular_methods/profiles_Payerne_CL61.png)

![Payerne CL61 signal/molecular ratio + window centres + optimal flagged cells](figs_paper_validation/molecular_methods/pcolor_Payerne_CL61.png)

**Amsterdam CHM15k (1064 nm)** — high-SNR reference; the signal tracks the molecular line cleanly from
~2 km to 7 km on clear nights:

![Amsterdam CHM15k vertical profiles](figs_paper_validation/molecular_methods/profiles_Amsterdam_CHM15k.png)

![Amsterdam CHM15k ratio + window centres + optimal flagged cells](figs_paper_validation/molecular_methods/pcolor_Amsterdam_CHM15k.png)

**EDT CL61 (910 nm, Edmonton 53.5°N)** — high-latitude site (also exercises the darkness-adaptive
night window and the WV correction):

![EDT CL61 vertical profiles](figs_paper_validation/molecular_methods/profiles_EDT_CL61.png)

![EDT CL61 ratio + window centres + optimal flagged cells](figs_paper_validation/molecular_methods/pcolor_EDT_CL61.png)

### What the detailed examples show

1. **Aerosol night (Payerne CL61 2026-03-12) — no method calibrates.** The gated methods select no
   eligible window (no clean molecular layer above the aerosol); `main` picks a degenerate high
   window (R²=0.06, scattering ratio **−190**, **negative constant**) that the pipeline's
   proportionality QC then rejects. The profile shows usable signal only below ~2 km (the aerosol).
2. **Method differentiation (clear nights):** `calipso` always sits highest (→ noisy), `earlinet`
   favours the lowest qualifying window, `improved`/`optimal`/`matlab` land mid-profile; the gated
   methods agree on the constant to a few %.
3. **`optimal`'s time-resolved flagging works.** The hatched cells mark the aerosol present only at the
   *start* of 2026-03-16, the thin cloud near the *end* of 2026-03-28, and the persistent
   boundary-layer aerosol at low altitude. `optimal` excludes those cells and fits on the time-cleaned
   mean, so the clean part of a partly-contaminated night stays usable.
4. **The x-axis now spans the full profile** (floor set to the smallest positive signal, ≤5 decades),
   so the signal–vs–molecular agreement is visible all the way to 7 km instead of being clipped.

---

## 7. Multi-site comparison and method selection

All seven methods were run across **six instruments at three sites** — Payerne (CHM15k, CL31, CL61),
Amsterdam (CHM15k), and EDT/Edmonton (CL51, CL61) — over the clear nights of a Feb–Mar 2026 every-2nd-
day sample (production WV correction on the 910 nm units). For each (instrument, method) we record the
number of nights that calibrate through the full pipeline, the **night-to-night lidar-constant CV**
(stability — the dominant quality for a calibration constant), and the median R² / temporal-CV /
rel-error. Full numbers:
[`method_comparison_multisite.md`](figs_paper_validation/molecular_methods/method_comparison_multisite.md).

![Methods across sites: usable nights and night-to-night stability](figs_paper_validation/molecular_methods/summary_methods_multisite.png)

**Calibration-constant time series per method** (same colour code; per instrument). The unstable
methods (`calipso` orange, `main` grey) swing wildly while the gated methods (`improved`, `optimal`,
`earlinet`, `matlab`, `bellini`) track each other tightly — the y-range is set from the stable group,
so the `calipso`/`main` excursions clip at the top:

![Calibration-constant time series per method, per instrument](figs_paper_validation/molecular_methods/timeseries_methods_multisite.png)

**Ranking (mean over the six instruments):**

| method | calibrated-fraction | mean CL_CV % | mean temporal_cv | verdict |
|---|---|---|---|---|
| **improved** | 0.21 | 11.3 | 0.24 | **best balance** — most usable of the clean group |
| **optimal** | 0.15 | 8.9 | 0.21 | cleanest + stable; flags & excludes contaminated layers |
| **earlinet** | 0.16 | 8.7 | 0.18 | cleanest, lowest yield |
| matlab | 0.30 | 11.4 | **1.71** | most nights, but admits aerosol (no above-aerosol gate) |
| bellini | 0.14 | 14.4 | 0.36 | strictest QC (lowest rel-error); 3–7 km costs SNR on 910 nm → fewer nights, higher CV |
| calipso | 0.59 | **112.6** | 0.45 | unusable scatter (chases the highest layer into noise) |
| main | 0.36 | **87.3** | 0.43 | unusable scatter (degenerate) |

- **`improved`, `optimal`, `earlinet` form a "clean" group** at CV 9–11 % and temporal CV ~0.2. They
  differ mainly in yield: `improved` calibrates the most nights; `earlinet`/`optimal` slightly fewer
  but a touch cleaner/steadier.
- **`matlab` gives the most usable-looking nights (0.30) but temporal CV 1.71** — lacking an
  above-aerosol gate, it admits temporally-variable (aerosol) windows. Its CV looks fine, but the
  aerosol contamination is a *systematic* accuracy risk, so it is not the default.
- **`calipso` and `main` are unusable** (CV 80–200 %): they calibrate many nights but the constants
  scatter wildly — `calipso` chases the highest "clean" layer into noise, `main` is degenerate. The
  time-series figure shows this directly: their orange/grey lines spike off-scale while the gated
  methods overlap in a tight band.
- **`bellini` (ALICENET)** sits mid-pack: its Breusch–Godfrey + E_CL QC give the **lowest fit
  rel-error** (2–3 % on CHM15k) and it does well on high-SNR CHM15k (Amsterdam 20 nights), but the
  3–7 km search band loses SNR on the 910 nm units, so it yields fewer nights and a higher CV (14 %)
  than the `improved`/`earlinet`/`optimal` group. A solid, literature-grounded method, but not the
  best fit for the lower-SNR ALCs in this network.
- **Per instrument:** CHM15k (high SNR) and CL61 calibrate well with the gated methods. **CL31 and
  CL51 (noisy 910 nm) rarely yield a clean molecular window** — the gated methods correctly reject
  almost all their nights (only the ungated `main`/`calipso` "calibrate", at CV > 100 %), confirming
  that **CL31/CL51 should be calibrated by the liquid-cloud method, not Rayleigh.**

**Selection:** `improved` is the best general default (and is already the shipped default — no change
needed). `optimal` is the best method where aerosol/cloud contamination matters (equal stability,
cleaner, and it flags & excludes contaminated layers). `matlab` only when maximum night count is
needed and aerosol contamination is acceptable; never `calipso`/`main` for ceilometers.

---

## 8. Recommendations

- **Keep `improved` as the production default** — the multi-site test confirms it as the best balance
  of usable nights and night-to-night stability across all instrument types (no `options.json`
  change needed).
- **Switch to `optimal` where aerosol/cloud is significant** (polluted sites, summer, smoke episodes):
  it is equally stable, cleaner, and uniquely flags & excludes contaminated time-altitude layers — set
  `"molecular_method": "optimal"`.
- **CHM15k (1064 nm) and CL61 (910 nm):** `improved` or `optimal` — both calibrate reliably.
- **CL31 / CL51 (910 nm, noisy):** Rayleigh rarely finds a clean molecular window — use the
  **liquid-cloud calibration** as primary, Rayleigh only opportunistically.
- **`bellini` (ALICENET)** is a sound, citable alternative — best where SNR is high (CHM15k) and the
  strictest fit quality is wanted; on the lower-SNR 910 nm ALCs its 3–7 km band yields fewer points.
- **Never `calipso` or `main`** for ceilometers (CV 80–200 %); `matlab` only with the aerosol caveat.
- **Keep rejecting bad nights** — "no calibration" (flag −2) that the Kalman skips beats a spurious
  constant from an aerosol/noise window.
- **Borrow from ATLID where cheap:** pre/post-echo background interpolation, a proper T/p molecular
  term, and treating the clean-air assumption as testable (the temporal flag is one such test).

### Files

- Core: `ALC_rayleigh_calibration/rayleigh_calibration/molecular_methods.py` (7 methods +
  `flag_contaminated_cells`), `rayleigh_fit.py` (dispatch + `signal_stack`),
  `config.py` + `options.json` (`molecular_method`), `calibration.py` (pass-through +
  `signal_stack` + `fit_inputs_out` capture hook).
- Test/plots: `ALC_rayleigh_calibration/compare_molecular_methods.py` (multi-site harness),
  `make_optimal_flowchart.py`, `test_methods_smoke.py`. Outputs in
  `figs_paper_validation/molecular_methods/`: `profiles_*` / `pcolor_*` (Payerne_CL61,
  Amsterdam_CHM15k, EDT_CL61), `summary_methods_multisite.png`,
  `timeseries_methods_multisite.png`, `optimal_flowchart.png`, `method_comparison_multisite.md`.
- Not yet committed (changes live in the `ALC_rayleigh_calibration` repo). `options.json` default
  stays `molecular_method="improved"`, `apply_wv_correction=1`, `molecular_source="standard"`.
