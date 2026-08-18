# Clear-sky-night overlap retrieval (network CHM15k) — M2 estimator + paper figure

Reconstruct the **applied (onboard/manufacturer) overlap function** `O(z)` of every network CHM15k
**without a hood measurement**, from the noise of routine clear-sky night L1, and build the paper figure.

## Methodology

The CHM15k firmware outputs `rcs_0 = raw · z² / O(z)`, so in the range-uncorrected view `P = rcs_0/z²`
the detector **electronic noise is amplified by `1/O(z)`** exactly as it would be under a dark hood, while
the atmospheric contributions are not. On a clear night the per-gate white (successive-difference)
variance splits into a signal-**independent** electronic term and signal-**dependent** shot and turbulence
terms, `σ²(z) = β(z) + α(z)·P + γ(z)·P²`; a per-gate non-negative regression of the white variance against
the mean signal across many clear nights isolates the electronic variance as the **intercept `β(z)`**
(estimator **"M2"**), and the overlap follows from `O(z) ∝ 1/√β(z)`, normalised to 1 in the full-overlap
far band (1500–3500 m). Retrievals are restricted to **vertically and temporally homogeneous** nights (a
layeredness and a first-third-vs-last-third drift screen reject moving/elevated aerosol), the raw `O(z)` is
de-spiked (an under-sampled gate can drive `β→0` and blow `1/√β` up) and fitted with a **floor-anchored
double-logistic** (`K=2`, `O(0)=0`, `O(∞)=1`) — the same two-inflection form the manufacturer uses in its
`.cfg`. Validated against the Payerne terminal-hood measurement the method reproduces the applied overlap
to **~2 %** (150–1500 m) and, applied blind, matches four independent factory `.cfg` files to 0.7–5.7 %, so
it is a passive, hood-free network QC for a drifted or mis-configured onboard overlap.

## Usage

```bash
# 1. compute the network overlap product (needs the L1 archive; set L1 in clearsky_overlap.py)
python drive_network.py            # -> network_overlaps.npz   (curves, wmos, nnights, spikes)

# 2. build the 25 °C temperature-model correction panel (needs D:/TEMP_MODELS/202606)
python part_c_correction.py        # -> correction25.npz

# 3. render the paper figure (uses the two .npz above; no archive needed)
python fig_paper_network.py        # -> ../../doc/reports/figs_overlap_from_noise/08_network_overlap_vs_correction.png
```

The two `.npz` products (Apr–Jun 2026 run) are committed, so **step 3 reproduces the paper figure
standalone**; steps 1–2 only need re-running to refresh the data (they read machine-local archives:
`L1 = D:/E-PROFILE_L1_2026` in `clearsky_overlap.py`, and `D:/TEMP_MODELS/202606` in `part_c_correction.py`).

## Files

| file | role |
|---|---|
| `clearsky_overlap.py` | core module — night loading, homogeneity screen, M1/M2/M3 noise-split estimators, night selection |
| `overlap_fit_general.py` | K-logistic overlap model + fitter (the floor-anchored `K=2` used for the network curves) |
| `drive_network.py` | run M2 (anchored `K=2`) over the CHM15k list → `network_overlaps.npz` (parallel) |
| `part_c_correction.py` | stack the temperature-overlap models, evaluate the correction applied at 25 °C → `correction25.npz` |
| `fig_paper_network.py` | paper figure: (1) network applied-overlap binscatter + median/IQR + earliest/latest to complete (QC'd); (2) 25 °C correction binscatter |
| `chm15k_stations.json` | the network CHM15k WMO list scanned by `drive_network.py` |
| `network_overlaps.npz`, `correction25.npz` | committed data products (Apr–Jun 2026) |

Full write-up (principle, validation, factory-cfg reproduction, tilt link, and the spectral/PSD variant used
only for QC): `doc/reports/overlap_from_clearsky_noise.md`.
