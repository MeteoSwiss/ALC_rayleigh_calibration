# Cloudnet droplet size at the calibration cloud base across Cloudnet sites

Draft of 2026-11-02. Extension of `04_multiple_scattering_v2.md` section 3.1 (eta study, step 1d), which took
the droplet equivalent-area radius a_G = 3.3 um of the multiple-scattering (eta) tables from three Cloudnet sites
(Payerne, Palaiseau, Lindenberg). This note asks whether that value holds at every Cloudnet site that produces
the droplet effective radius product (`der`). It reports measurements only: **no new value is proposed and nothing
is applied to the operational code.** Every number is read from a CSV or markdown table in
`C:/DATA/Projects/202611_eta_study/cloudnet_sites/` (named with each table); if the run is repeated, take the numbers
from those files, not from this text.

## Summary

- **3.3 um lies inside the network range but on its low side.** Over the 30 Cloudnet sites inside the E-PROFILE
  domain (the operational CAMS box N27-80, W30-E45) with at least 20 days of data, the site medians of a_G (lowest
  100 m above the liquid base, all liquid scenes with base 0.25-2.4 km above the site) have a median of
  **3.66 um** (P25-P75 3.42-3.94). The day-weighted value is 3.61 um. The season-balanced value is 3.48 um, over the
  20 sites with all four seasons. Only 6 of the 30 site medians are below 3.3 um. 26 of 30 lie inside the 2.5-5 um
  range that `04_multiple_scattering_v2.md` uses as the a_G uncertainty. 11 of 30 (37 %) are within +-10 % of 3.3 um
  and 22 of 30 (73 %) within +-20 %.
- **Why 3.3 sits low:** two of the three step-1d sites are among the lowest-a_G sites of the network: Lindenberg
  3.00 um and Palaiseau 3.21 um (Payerne 3.65 um). With the section-6 sensitivity of the v2 report (a_G +-0.5 um ->
  C_L level at 1 km -+2.0-2.2 %), the median site (+0.18 to +0.36 um above 3.3) corresponds to about -0.7 to -1.6 %
  in C_L at 1 km. The between-site standard deviation (0.7 um) corresponds to about +-3 %.
- **The calibration profiles see the same droplets as all scenes at their site.** At the seven co-located Vaisala
  units the calibration-profile median is 0.97-1.10 times the all-scene median of the same days. This includes the
  three new co-located units: Chilbolton CL51 4.29 um (a weak match, see below), L'Aquila CL51 5.47 um, Norunda CL51
  3.42 um. The two CL31 units 9-11 km away give Middle Wallop 4.62 um (1.13) and Akrotiri 3.82 um (1.07).
- **No latitude dependence.** Spearman rank correlation of the site medians with |latitude|: -0.05 (p 0.78)
  scene-weighted, +0.01 (p 0.95) day-weighted, -0.35 (p 0.14) season-balanced. Arctic and boreal sites (Ny-Ålesund
  3.88, Kenttärova 3.61, Hyytiälä 3.37, Norunda 3.14 um) and Mediterranean ones (Lampedusa 4.38, Limassol 3.83,
  Granada 3.82, Potenza 2.90 um) span the same range.
- **What drives the spread.** The largest a_G go with the lowest N_scaled (L'Aquila 5.59 um at 31 cm-3, Cabauw
  5.52 at 119, Lampedusa 4.38 at 82): figure panel (c), close to a_G ~ N^(-1/3) at a given water content. This is
  partly physics (clean, marine air) and partly the retrieval: r_e and N trade against each other through the
  liquid water path (LWP) scaling (see caveats). The two mountain-top stations sit low (Schneefernerhaus 1.86 um at
  2653 m, Maïdo 2.40 um at 2200 m). At many sites the seasonal cycle is as large as the spread between sites:
  the winter (DJF) / summer (JJA) ratio is 0.62-1.63, median 1.05, and DJF is larger at 13 of 21 sites.
- **The step-1d inputs moved.** Lindenberg was reprocessed on the portal (CloudnetPy 1.92 -> 1.97, 2026-10-17): its
  CL61 calibration-profile a_G drops 5.5 % (3.04 -> 2.88 um) and its all-scene value 2.7 %. Pooled over the four
  step-1d units the median is now 3.24 um (step 1d: 3.31 um).

This note does not propose a new value. Which a_G the eta tables should use, and whether the per-site differences
matter for individual units (L'Aquila, Cabauw-like sites), is left to the eta-study decision.

![Cloudnet a_G across sites: per-site summary, latitude, a_G vs N_scaled](figs_eta_cloudnet_sites/fig_aG_summary.png)

![Cloudnet a_G distributions, one panel per site](figs_eta_cloudnet_sites/fig_aG_sites.png)

## Results

### Per site (all liquid scenes)

Lowest 100 m above the lowest liquid base, base 0.25-2.4 km above the site, `der_retrieval_status == 1` pixels,
2025-01-01..2026-10-31. Sites north to south; * = outside the E-PROFILE CAMS box; "thin" = fewer than 20 days.
Source: `cloudnet_sites/site_table.csv` (all columns: 50 / 150 m layers, fixed-N r_e, N quartiles, status-1 fraction,
effective layer depth) and `note_tables.md`.

| site | lat | alt (m) | radar (main), band | LWP from | days | scenes | a_G median [P25-P75] (um) | day-weighted | season-balanced | N_scaled median (cm-3) |
|---|---|---|---|---|---|---|---|---|---|---|
| Ny-Ålesund (SJ) | 78.9 | 19 | RPG-FMCW-94, W | HATPRO | 534 | 96622 | **3.88** [2.8-5.1] | 4.16 | 4.28 | 266 |
| North Slope of Alaska (US) * | 71.3 | 8 | KAZR, Ka | WVR-1100 | 231 | 51309 | **4.06** [3.2-5.0] | 4.07 | 3.99 | 429 |
| Pallasjärvi (FI) | 68.0 | 275 | RPG-FMCW-94, W | HATPRO | 24 | 10090 | **3.94** [2.5-5.8] | 4.49 |  | 91 |
| Kenttärova (FI) | 68.0 | 345 | RPG-FMCW-94, W | radar 89 GHz channel | 320 | 80971 | **3.61** [2.3-4.7] | 3.07 | 3.05 | 164 |
| Sammaltunturinlaatas (FI) | 68.0 | 460 | RPG-FMCW-94, W | radar 89 GHz channel | 25 | 5700 | **5.02** [3.0-7.6] | 4.88 |  | 50 |
| Snæfellsnes (IS) | 64.8 | 14 | RPG-FMCW-94, W | HATPRO | 21 | 3109 | **4.05** [3.1-5.2] | 3.92 |  | 123 |
| Hyytiälä (FI) | 61.8 | 150 | RPG-FMCW-94, W | HATPRO | 464 | 113238 | **3.37** [2.2-4.5] | 3.39 | 3.43 | 322 |
| Norunda (SE) | 60.1 | 46 | RPG-FMCW-94, W | radar 89 GHz channel | 450 | 117948 | **3.14** [1.8-4.6] | 2.79 | 2.80 | 574 |
| Rzecin (PL) | 52.8 | 57 | BASTA, W | HATPRO | 51 | 2449 | **3.68** [2.4-5.7] | 4.20 |  | 738 |
| Lindenberg (DE) | 52.2 | 104 | MIRA-35, Ka | HATPRO | 552 | 170296 | **3.00** [2.3-4.0] | 3.05 | 3.13 | 241 |
| Falkenberg (DE) | 52.2 | 73 | RPG-FMCW-94, W | HATPRO | 94 | 14768 | **3.47** [2.3-4.7] | 3.37 |  | 526 |
| Cabauw (NL) | 52.0 | -1 | RPG-FMCW-94, W | HATPRO | 485 | 104212 | **5.52** [4.0-7.6] | 5.80 | 5.39 | 119 |
| Melpitz (DE) | 51.5 | 86 | MIRA-35, Ka | HATPRO | 205 | 57798 | **3.45** [2.7-4.6] | 3.60 | 3.49 | 209 |
| Leipzig (DE) | 51.4 | 126 | MIRA-35, Ka | HATPRO | 363 | 104375 | **3.58** [2.6-4.8] | 3.50 | 3.47 | 194 |
| Leipzig LIM (DE) | 51.3 | 126 | RPG-FMCW-94, W | HATPRO | 75 | 18779 | **3.48** [2.7-4.7] | 3.56 |  | 406 |
| Chilbolton (GB) | 51.1 | 85 | RAL Space Copernicus, Ka | HATPRO | 49 | 14923 | **3.98** [3.1-5.0] | 4.16 |  | 335 |
| Jülich (DE) | 50.9 | 111 | MIRA-35, Ka | HATPRO | 537 | 180298 | **3.23** [2.5-4.3] | 3.22 | 3.26 | 268 |
| Palaiseau (FR) | 48.7 | 156 | BASTA, W | HATPRO | 489 | 108291 | **3.21** [2.3-4.4] | 3.32 | 3.36 | 426 |
| Kappelrodeck (DE) | 48.6 | 199 | RPG-FMCW-94, W | HATPRO | 360 | 81843 | **3.76** [2.4-5.0] | 3.64 | 3.64 | 310 |
| Munich (DE) | 48.1 | 538 | MIRA-35, Ka | HATPRO | 494 | 118834 | **3.67** [2.7-4.9] | 3.62 | 3.73 | 145 |
| Schneefernerhaus (DE) | 47.4 | 2653 | MIRA-35, Ka | HATPRO | 136 | 15629 | **1.86** [1.3-2.6] | 1.89 |  | 1020 |
| Payerne (CH) | 46.8 | 491 | RPG-FMCW-94, W | HATPRO | 490 | 135862 | **3.65** [2.8-4.7] | 3.53 | 3.40 | 274 |
| Cluj-Napoca (RO) | 46.8 | 459 | RPG-FMCW-94, W | HATPRO | 403 | 72431 | **3.64** [2.6-4.7] | 3.31 | 3.41 | 265 |
| Galați (RO) | 45.4 | 40 | RPG-FMCW-94, W | HATPRO | 359 | 69906 | **4.13** [2.9-5.8] | 3.79 | 3.84 | 168 |
| Bucharest (RO) | 44.3 | 77 | MIRA-35, Ka | HATPRO | 427 | 79604 | **3.42** [2.4-4.7] | 3.47 | 3.57 | 270 |
| L'Aquila (IT) | 42.4 | 683 | RPG-FMCW-94, W | radar 89 GHz channel | 489 | 109273 | **5.59** [3.9-7.5] | 5.18 | 5.28 | 31 |
| Potenza (IT) | 40.6 | 760 | MIRA-35, Ka | HATPRO | 197 | 46654 | **2.90** [2.2-3.9] | 3.05 | 3.14 | 323 |
| Graciosa (PT) | 39.1 | 30 | KAZR, Ka | MWR-3C | 240 | 154976 | **3.93** [2.6-6.2] | 3.82 |  | 184 |
| Granada (ES) | 37.2 | 680 | RPG-FMCW-94, W | HATPRO | 315 | 66450 | **3.82** [2.5-5.1] | 3.70 | 3.79 | 224 |
| Southern Great Plains (US) * | 36.6 | 316 | KAZR, Ka | WVR-1100 | 344 | 59201 | **3.87** [2.8-5.3] | 4.38 | 4.39 | 540 |
| Lampedusa (IT) | 35.5 | 45 | MIRA-35, Ka | HATPRO | 376 | 99665 | **4.38** [3.3-6.0] | 4.41 | 4.42 | 82 |
| Limassol (CY) | 34.7 | 10 | MIRA-35, Ka | HATPRO | 221 | 18755 | **3.83** [2.8-5.5] | 4.01 |  | 97 |
| Bankhead National Forest (US) * | 34.3 | 293 | KAZR, Ka | MWR-3C | 451 | 71411 | **3.59** [2.6-5.1] | 4.16 | 4.34 | 654 |
| Phoenix (US) * | 33.6 | 379 | KAZR, Ka | MWR-3C | 33 | 1653 | **4.94** [3.6-6.7] | 4.92 |  | 739 |
| Mindelo (CV) * | 16.9 | 13 | MIRA-35, Ka | HATPRO | 414 | 204539 | **3.64** [2.5-5.1] | 3.52 | 3.46 | 109 |
| ATTO-Campina (BR) * | -2.2 | 30 | MIRA-35, Ka | Radiometrics | 475 | 71612 | **1.58** [1.1-2.3] | 1.52 | 1.52 | 143603 |
| Maïdo Observatory (RE) * | -21.1 | 2200 | BASTA, W | HATPRO | 229 | 17888 | **2.40** [1.6-3.5] | 2.23 | 2.17 | 1924 |
| Tāwhaki (NZ) * | -43.8 | 5 | RPG-FMCW-94, W | radar 89 GHz channel | 81 | 23733 | **4.48** [3.1-5.9] | 4.38 |  | 118 |
| Invercargill (NZ) * | -46.4 | 5 | MIRA-35, Ka | HATPRO | 404 | 134296 | **4.36** [3.2-5.7] | 4.30 | 4.33 | 86 |
| Neumayer Station (AQ) * | -70.7 | 43 | RPG-FMCW-94, W | HATPRO | 422 | 45110 | **3.23** [2.5-4.1] | 3.51 | 3.60 | 244 |
| Troll Station (AQ) * | -72.0 | 1320 | RPG-FMCW-35, Ka | LHATPRO | 166 | 2417 | **3.64** [1.2-7.6] | 4.51 | 3.71 | 169 |

- "radar (main)" is the radar of most days; Bucharest, Leipzig, Lindenberg, Munich, Jülich, Invercargill and others
  changed radar during the period (`site_radar_table.csv`). "LWP from": the radiometer named in the file, or the
  89-GHz passive channel of the RPG radar where no radiometer is listed.
- Only 27-38 % of the liquid-base scenes in the window carry status-1 data in the lowest 100 m at the continental
  sites, 45-59 % at the marine ones (Graciosa, Lampedusa, Mindelo) and 2-14 % at the polar ones (mixed phase).
- The ATTO-Campina retrieval is not usable here: N_scaled ~1.4 x 10^5 cm-3 is unphysical, so its 1.58 um reflects an
  LWP / reflectivity mismatch, not droplets.

### Calibration profiles of the Vaisala units

Source: `cloudnet_sites/calib_table.csv`, `calib_pooled.csv`. "site, same days" = all-scene median of the site on
the days on which the unit has calibration profiles with status-1 data.

| unit | distance (km) | accepted profiles (days) | with status-1 data (days) | a_G median [P25-P75] (um) | site, all scenes | site, same days | calibration / same days |
|---|---|---|---|---|---|---|---|
| Payerne CL61 | 0.1 | 5035 (49) | 2160 (41) | **3.37** [2.6-4.2] | 3.65 | 3.32 | 1.02 |
| Payerne CL31 (distorted unit) | 0.1 | 39772 (262) | 16689 (208) | **3.56** [2.9-4.4] | 3.65 | 3.57 | 1.00 |
| Palaiseau CL31 | 0.2 | 34893 (277) | 10459 (227) | **3.10** [2.5-3.9] | 3.21 | 3.16 | 0.98 |
| Lindenberg CL61 (>=2025-05) | 0.3 | 27814 (232) | 12980 (192) | **2.88** [2.3-3.7] | 3.00 | 2.95 | 0.98 |
| Chilbolton CL51 | 0.1 | 5844 (34) | 745 (25) | **4.29** [3.4-5.0] | 3.98 | 3.97 | 1.08 |
| L'Aquila CL51 | 0.1 | 48165 (325) | 20737 (267) | **5.47** [4.1-7.2] | 5.59 | 5.63 | 0.97 |
| Norunda CL51 | 0.1 | 30892 (263) | 11564 (203) | **3.42** [2.0-4.9] | 3.14 | 3.12 | 1.10 |
| Middle Wallop CL31 (9 km from Chilbolton) | 9.2 | 11014 (39) | 1432 (28) | **4.62** [3.8-5.6] | 3.98 | 4.08 | 1.13 |
| Akrotiri CL31 (11 km from Limassol) | 10.9 | 1767 (31) | 110 (14) | **3.82** [2.5-5.2] | 3.83 | 3.56 | 1.07 |

Pooled over units, every accepted profile counting: the four step-1d units 3.24 um (P25-P75 2.5-4.1; 42 288
profiles); all seven co-located units 3.66 um (2.6-5.1; 75 334 profiles; L'Aquila contributes 28 % of them), median
of the seven unit medians 3.42 um. The calibration selection does not bias the droplet size: the ratio to the same
days' all-scene median stays within 0.97-1.10 at the co-located units.

The Chilbolton matches are the least certain. Only 28 % of the CL51 profiles, and 27 % of the Middle Wallop
profiles, find Cloudnet liquid within +-150 m of their CBH; the other co-located units reach 63-94 %. A profile
without liquid near the CBH almost never has liquid elsewhere in the Cloudnet profile. That share is 92-98 % at the
co-located units and about 100 % for the Chilbolton CL51 (`site_table.md`, "Profiles without Cloudnet liquid"). So
the match is limited by where Cloudnet diagnoses liquid at all, not by a height offset between the two cloud bases.
Both Chilbolton values therefore rest on a small and possibly unrepresentative part of the profiles (745 and 1 432).

### Spread between sites, latitude, climate, season, radar

Source: `cloudnet_sites/site_table.md` (section "Across sites"), `across_sites.csv`, `season_table.csv`,
`site_radar_table.csv`.

| E-PROFILE-domain sites (>= 20 days) | sites | median of site medians [P25-P75] (um) | min - max (um) | SD (um) | below 3.3 | Spearman with abs(lat) (p) |
|---|---|---|---|---|---|---|
| scene-weighted | 30 | 3.66 [3.42-3.94] | 1.86 (Schneefernerhaus) - 5.59 (L'Aquila) | 0.73 | 6 | -0.05 (0.78) |
| day-weighted | 30 | 3.61 [3.33-4.12] | 1.89 (Schneefernerhaus) - 5.80 (Cabauw) | 0.75 | 6 | +0.01 (0.95) |
| season-balanced (all four seasons) | 20 | 3.48 [3.33-3.80] | 2.80 (Norunda) - 5.39 (Cabauw) | 0.68 | 5 | -0.35 (0.14) |

- **Latitude:** none detectable (table). Altitude: weak and not significant (Spearman -0.24 to -0.33, p 0.08-0.31),
  carried by the mountain-top stations. L'Aquila (683 m) and Granada (680 m) are not low.
- **Within-site variability dominates:** the median within-site IQR of the daily medians is 1.7 um, more than twice
  the between-site SD.
- **Season:** DJF / JJA day-weighted ratio 1.05 (median over 21 sites, range 0.62-1.63). Winter droplets are larger
  at 13 of 21 sites, strongly so at L'Aquila (5.70 vs 3.50 um), Hyytiälä (4.32 vs 2.85), Ny-Ålesund (5.56 vs 3.64),
  Cluj (4.05 vs 2.82) and Payerne (4.14 vs 3.14). The cycle is flat or reversed at Lindenberg, Palaiseau, Jülich,
  Bucharest and Lampedusa. A site's record therefore matters: Chilbolton's 3.98 um is a single winter.
- **Radar band:** no systematic difference across sites. Day-weighted medians: W-band (94 GHz) 3.67 um over 18
  sites, Ka-band (35 GHz) 3.55 um over 12. Where a site changed radar, the W-band period gave larger a_G (Bucharest
  3.99 vs 3.21 um with the MIRA-35, Lindenberg 3.39 vs 3.03, Leipzig 3.58 vs 3.46, Invercargill 5.17 vs 4.31).
  These periods differ in season and length, so a radar effect cannot be separated from the seasonal cycle with
  this data.

### Robustness

- **Portal reprocessing** (`cloudnet_sites/version_effect.md`, same code on the July files): Lindenberg -2.7 %
  (all scenes, same days) and -5.5 % (CL61 calibration profiles, 3.04 -> 2.88 um); Palaiseau, Payerne and Chilbolton
  unchanged.
- **Layer quantization.** The "lowest 100 m" holds whole gates (gate-centre rule), so it really spans 89-109 m
  depending on the grid (Troll 60 m; column "'100 m' spans" of `site_table.md`). On Payerne's 22.36 m chirp, the step-1d base
  edge (median spacing, 21.04 m) put the 5th gate centre at 99.96 m, inside the layer. The local spacing puts it at
  100.6 m, outside. Payerne's calibration-profile a_G thus moves from 3.48 to 3.37 um (CL61) and from 3.60 to 3.56 um
  (CL31) between step 1d and this note. Layer depth is defined only to +-half a gate, i.e. a few % in a_G.
- **Layer depth:** a_G grows with height in every site: lowest 50 / 100 / 150 m, e.g. Lindenberg 2.79 / 3.00 /
  3.34, Payerne 3.35 / 3.65 / 3.89, Hyytiälä 3.04 / 3.37 / 3.63 um (`site_table.csv`).
- **Matching tolerance:** 15 min instead of 6 min adds 0-0.8 % of the profiles with status-1 data; per scene,
  Palaiseau and Lindenberg are identical to step 1d.

## Retrieval caveats (why not all of the spread is microphysics)

- **Retrieval scaling.** The LWP-scaled radius follows from the retrieval's own assumptions: lognormal of fixed
  width, N constant with height, reflectivity Z ~ N D^6, liquid water content ~ N D^3. Under these, r_e scales as
  (radar calibration)^(1/3) x LWP^(-1/3), and N_scaled as LWP^2 / (calibration). So 1 dB of radar calibration error
  is about 8 % in a_G, and a 30 % LWP error is about -8 / +13 % in a_G. A site with a_G high and N_scaled low (L'Aquila,
  Cabauw, and the Norunda N_scaled tail with P75 ~7000 cm-3) may have an LWP or calibration issue as well as clean
  air. L'Aquila, Kenttärova, Norunda, Sammaltunturinlaatas and Tāwhaki take their LWP from the single 89-GHz
  channel of the RPG radar, not from a multi-channel radiometer.
- **Fixed width:** sigma = 0.35 enters both the Cloudnet retrieval and a_G = 0.832 r_e; it is not measured.
- **Geometry:** mountain-top stations (Schneefernerhaus, Maïdo) and the polar sites (status-1 in 2-14 % of the
  scenes, mixed phase) do not sample the cloud regime that an E-PROFILE ceilometer on a plain calibrates on.
- **Sampling:** records differ in length and season (Chilbolton one winter; Pallasjärvi, Sammaltunturinlaatas,
  Snæfellsnes, Rzecin and Falkenberg are 1-4-month campaigns). The season-balanced column removes most of this
  for the 20 long records.
- **Liquid detection** depends on the site's lidar (CHM15k, CL51, CL61, PollyXT, CS135, ...), which sets where the
  Cloudnet liquid base lies (60-95 m below the ceilometer CBH at the step-1d sites).

## Data

**Cloudnet `der`, 2025-01-01..2026-10-31.** The Cloudnet portal API lists 41 sites with `der` (17 084 daily files,
5.57 GB; `der_availability.csv`, file list `der_filelist.csv`). All of them were downloaded with SHA-256 verification
into `D:/Cloudnet_Droplet_radius/<site>/` (manifest `download_manifest_20261102.csv` there). ACTRIS Cloudnet data,
CC-BY 4.0.

- **Chilbolton** has `der` only from 2025-01-07 to 2025-03-07 (58 days, the days step 1d already had): the portal holds
  no Chilbolton cloud-radar data after 2025-03-07 (lidar continues to 2026-10, microwave radiometer to 2026-01), so no
  later `der` exists. Its sample is one winter.
- **Warsaw** has only radiometer data (no radar, hence no `der`); **Mace Head** has no Cloudnet product in 2025-2026.
- **Reprocessing.** 847 of the 1 682 files of the July collection used by step 1d (all 546 Lindenberg days, about 300
  Palaiseau days, 2 Payerne days) have been reprocessed on the portal since (versions of 2026-07 and 2026-10). The
  latest versions are used here; the July files are kept unchanged, and the effect is sized below.

**E-PROFILE Vaisala ceilometers at or near the Cloudnet sites** (`validation/scope_l1_2026_census.json`). The CHM15k
saturates in liquid clouds and has no cloud calibration, so only CL31/CL51/CL61 units are relevant.

| Cloudnet site | Vaisala unit (distance) | per-profile calibration scenes |
|---|---|---|
| Payerne | CL61 0-20000-0-06610_C, CL31 0-20000-0-06610_B (0.1 km; B = the distorted unit) | July ms_fov run (step 1d) |
| Palaiseau | CL31 0-250-1001-07151_A (0.2 km) | July ms_fov run (step 1d) |
| Lindenberg | CL61 0-20000-0-10393_C (0.3 km; from 2025-05, after its hardware step) | July ms_fov run (step 1d) |
| Chilbolton | CL51 0-826-300-3_A (0.1 km); CL31 Middle Wallop 0-20000-0-03749_A (9.2 km) | new dump (below) |
| L'Aquila | CL51 0-20008-0-CEO_A (0.1 km) | new dump |
| Norunda | CL51 0-20008-0-NOR_A (0.1 km) | new dump |
| Limassol | CL31 Akrotiri 0-20000-0-17601_A (10.9 km) | new dump |
| Falkenberg (campaign, 4.6 km from Lindenberg) | Lindenberg CL61 (4.8 km) | matched to Lindenberg's own Cloudnet only |
| Cabauw, Granada, Bucharest, Jülich, Leipzig, Leipzig LIM, Melpitz, Munich, Potenza | CHM15k only (0-5 km) | none (no cloud calibration) |
| Kappelrodeck, Jülich (Vaisala), Limassol (second unit) | Strasbourg CL31 35 km, Köln CL51 37 km, Elsenborn CL31 50 km, CAO-AMX CL51 40 km | not used (too far) |
| Hyytiälä, Kenttärova, Galați, Cluj, Lampedusa, Rzecin, Ny-Ålesund, Schneefernerhaus, the Finnish and Icelandic campaign sites, all sites outside Europe | no E-PROFILE ALC within 40 km | none |

The new per-profile scenes come from `validation/eta_study/cloudnet_sites_calib_dump.py`. It runs the operational
liquid-cloud calibration (`rayleigh_availability/cloud_profile_dump.process_day`, the `_do_cloud` replica, 0.4 deg CAMS
water-vapour correction, 2026-07 eta tables) on each unit's known cloud-calibration days
(`rayleigh_availability/cbh_native/cbh_scenes.csv`) that have a Cloudnet `der` file. Per-day profile counts equal those
of `cbh_scenes.csv` on 95-100 % of the days; a few days differ by at most 47 profiles. Result: Chilbolton CL51 5 844 accepted profiles (34 days), L'Aquila CL51 48 165 (325 days),
Norunda CL51 30 892 (263 days), Middle Wallop CL31 11 014 (39 days), Akrotiri CL31 1 767 (31 days).

## Method (as step 1d)

- Liquid pixel = `der_retrieval_status` 1-3 (status 4, "surrounding ice", never carries a radius).
- **All scenes (per site):** every 30-s Cloudnet profile whose *lowest* contiguous liquid layer (the layer a
  ceilometer reports as its first cloud base) has its base 0.25-2.4 km above the site, the calibration cloud-base
  window. The Cloudnet liquid base sits 60-95 m below the ceilometer CBH at the step-1d sites.
- **Calibration profiles (per Vaisala unit):** nearest Cloudnet profile within 15 min (step 1d: 6 min), the lowest
  liquid gate within +-150 m of the ceilometer CBH (cloud base MSL = CBH + Cloudnet site altitude; the ceilometer's
  own altitude for the two units 9-11 km away), followed down to the base of its contiguous liquid run.
- Medians over the `der_retrieval_status == 1` (reliable) pixels of the lowest 50 / 100 / 150 m above the liquid
  base. r_e = `der_scaled` (the LWP-scaled radar retrieval, method 2 of the CloudnetPy `der` product). a_G = r_e
  exp(-1.5 sigma^2) = 0.832 r_e for the lognormal width sigma = 0.35 that the CloudnetPy retrieval assumes. `der`
  (fixed N = 200 cm-3) and `N_scaled` are carried along.
- Weightings: per scene (every 30-s profile counts), per day (median of the daily medians) and season-balanced
  (median of the four seasonal day-medians, for sites with at least 10 days in every season).
- The base edge uses the local gate spacing. Step 1d used the median spacing, which is identical on regular grids,
  but several RPG-FMCW-94 grids change spacing with the chirp (Payerne 25.6 m below 1 km and 22.4 / 9.6 m above;
  Hyytiälä 9.6 / 34 m). Palaiseau and Lindenberg reproduce step 1d per scene; Payerne moves for the reason given in the
  robustness section.

Code: `validation/eta_study/cloudnet_der_download.py` (inventory, download),
`validation/eta_study/cloudnet_sites_calib_dump.py` (calibration scenes),
`validation/eta_study/cloudnet_sites_droplet_size.py` (analysis, tables, figures; `--files july` re-runs it on the
July collection).
