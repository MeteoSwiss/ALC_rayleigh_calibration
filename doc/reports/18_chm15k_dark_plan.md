# 18 — CHM15k: identify the electronic distortion, deliver an unbiased Rayleigh constant

Status: PLAN (2026-08-23), built from the remote-dark campaign's measured verdicts (report 16).
The stake, measured at Payerne: the dark is worth **+24 % on C_L**, most of the profile
gradient (−10.5 → −2.5 %/km), and +24 usable nights/yr — and the CHM15k has no cloud-
calibration fallback (it saturates in liquid clouds).

## The evidence this plan stands on

| # | Fact (all measured, fig refs in report 16) | Consequence for the plan |
|---|---|---|
| 1 | Single-station sky retrieval of the dark amplitude is impossible (injection gain 0.48, contamination ×4; fig13/14) | absolute amplitude must be MEASURED, not fitted |
| 2 | The near-range dark SHAPE is sky-retrievable (corr 0.98) | shape checks and triage are remote |
| 3 | Clean-night selection removes ⅔ of the aerosol pollution, then saturates (fig14) | evidence-based upper bounds are usable for triage |
| 4 | No clean fit window exists: low = aerosol-biased, high = dark-biased (fig15) | window retuning mitigates nothing by itself |
| 5 | Opaque-cloud natural hood fails for CHM15k: saturation artefact ×320–3000, non-monotonic in charge, no Q→0 rescue (fig16/18/19) | clouds are a *health probe*, never a dark source |
| 6 | Co-location gives exact per-unit differences; one anchored unit calibrates a cluster (M2 POC) | multi-unit sites need ONE covered unit |
| 7 | Covered-telescope nights work but drift for hours; T-resolved reduction needed; dark is per-optical-module (D swap; fig9/12/20) | protocol: ≥24 h, T_int-binned, era-keyed by module serial |
| 8 | Noise floor per night detects hardware changes (×0.55 step, burn-in +16 %/30 d); cloud y(Q) knee tracks detector state (fig12/19) | monitoring layer says WHEN to re-measure |
| 9 | ALC_DARK_PROFILE subtraction is already plumbed into the operational pipeline | application is a registry + config problem, not new physics |

Principle: **identification is remote, the absolute amplitude is measured once per unit-era,
monitoring keeps it valid.**

## P0 — Stakes and triage (desk, existing data + one balfrin job)

1. Extend the v4 night caches to every network CHM15k (balfrin, same machinery as the
   Schiphol extension).
2. Per unit: upper bound on the dark's calibration impact = rerun the Rayleigh fit with the
   (contaminated, hence conservative) evidence intercept subtracted → ΔC/C bound. The
   contamination only inflates the bound, never hides a large dark.
3. Classify: **< 2 %** ignore; **2–5 %** monitor only; **> 5 %** campaign list, geographic
   spread respected. Expected outcome: a minority of units justify a covered night at all.
4. Byproduct: near-range shape stability across units (fact 2) — if the gexp family shape is
   uniform, future amplitude-only reductions get cheaper.

## P1 — Covered-night reduction pipeline (code, validated on Payerne before any campaign)

- **Auto-detection of covered periods** in the whole archive first: a covered telescope is
  unmistakable (no atmospheric structure, no clouds, profile ≈ dark family, noise floor
  unchanged). Maintenance covers happen — free darks may already exist network-wide.
- Reduction: discard settle, then **T_int-binned median profiles** (Pinstrument(r, T) à la
  Le & O'Connor — fig20 proved a single session median under-resolves the hour-scale drift),
  era-keyed by `optical_module_id`, uncertainty from within-session spread.
- QC gates: session ≥ 24 h target (diurnal T span), repeat-session agreement, and the
  validation: reproduce Payerne's +24 % / gradient repair from its hood sessions.
- Output: per-unit-era dark npz in the exact `ALC_DARK_PROFILE` format.

## P2 — Campaign wave 1 (operator email, triaged units only)

One-page instruction sheet: cover the **receiver telescope** (not the laser) for ≥ 24 h,
note begin/end, done — the data arrives through the normal L1 flow, weather-independent.
Reduction is passive as files land. Wave 1 = the P0 ">5 %" list.

## P3 — Co-location multiplier

At multi-unit sites (Schiphol, Lindenberg, …): one covered unit anchors the cluster; the
M2 same-sky regression transfers the absolute dark to siblings; triangle closure is the QC.

## P4 — Operational application (v2.3)

- Registry keyed `(wmo, ident, optical_module_id, era)` feeding `ALC_DARK_PROFILE`.
- **Kalman restart at era boundaries** (the CL31 whole-year lesson: smoothing across a swap
  is 2.5× wrong).
- Acceptance: Payerne canaries still rejected, availability/gradient indicators improve or
  hold, network ΔC distribution matches the P0 bounds.

## Monitoring layer (already built — productize into the daily/weekly run)

- Nightly noise floor + step alarm per unit (fig12 machinery) → detects swaps even when
  metadata is silent, triggers re-measurement.
- Winter cloud-fingerprint y(Q) knee + kernel (fig19) → detector ageing between campaigns.
- Evidence near-range shape check (fact 2) → flags shape changes.

## Cheap asks and long shots

- **One email to Lufft**: does factory QC keep per-serial termination/dark profiles? A yes
  shortcuts the whole campaign.
- The archive sweep for accidental covers (in P1) — potentially large free coverage.

## Explicitly out

Sky-only absolute retrieval (proven impossible), cloud-derived darks for CHM15k (proven
impossible), network-wide window retuning as a dark fix (fig15), and any shape transfer
between units without a per-unit measurement (the CL61 lesson generalized).

## P0 execution log

**2026-08-23 — jobs launched**: caches+bounds for all 154 census CHM15k (balfrin postproc,
restartable), archive cover-sweep back to 2015 (no CAMS needed).

**First flagged unit inspected — MESSINA (0-20000-0-00203_A, bound +225 %), fig21.** Not a
monster dark in absolute terms (evidence amplitude only x1.1 Payerne's) but a SICK unit:

- noise floor DOUBLED smoothly over 20 months (no step -> not a swap; progressive laser-power
  decay, the Le & O'Connor early-unit pattern);
- the measured signal flattens to a positive pedestal above ~5 km where the molecular should
  keep decaying -> the 4.5-6.5 km "Rayleigh" fit is mostly fitting OFFSET, nightly constants
  scatter over 3 decades, and v2.2's validity gates already reject most recent nights;
- the two are one mechanism: rcs_0 is laser-normalised, so a FIXED raw electronic offset
  inflates as 1/laser exactly like the noise floor - the bound and the floor rise together.

Triage consequences: (i) cross-reference every bound with the unit's noise-floor TREND -
rising-floor units have 1/laser-inflated bounds and are maintenance cases first, campaign
cases second; (ii) Messina goes to the operator-notification list (laser end-of-life), and
any covered night there should wait until after servicing.

![Messina inspection](figs_remote_dark/fig21_MESSINA.png)

**Second flagged unit — PLOVDIV_UNI (0-100-20000-000627_A, bound +402 %), fig21b.** The
opposite pathology to Messina, and the pair validates the triage cross-check by example:

- noise floor FLAT over 20 months (healthy laser) - not a degradation case;
- one era step in Jan-2025: the first days sit x1e11 away in units (firmware/processing
  change) - era-key any dark work from ~2025-01-15;
- the measured signal saturates on a constant positive PEDESTAL above ~7 km where the
  molecular should keep decaying; the pedestal/molecular crossing sits right at the 4.5-6.5 km
  fit window, so the fit splits between pedestal and atmosphere -> bound +402 % with an
  otherwise stable constant series.

A fixed pedestal on a healthy unit is exactly what a covered night measures directly ->
Plovdiv goes on the P2 campaign list (era-keyed), and its stability makes it a good first
target. Contrast card: Messina = maintenance first; Plovdiv = campaign first.

![Plovdiv inspection](figs_remote_dark/fig21_PLOVDIV.png)
