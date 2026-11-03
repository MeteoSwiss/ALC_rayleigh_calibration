# Cloudnet cloud extinction at the calibration cloud base across Cloudnet sites, and what explains the site differences

Draft of 2026-11-03. Companion of `eta_cloudnet_sites_2026-11.md` (droplet size a_G at the 41 Cloudnet sites that
produce the `der` product) for the second input of the generation-2 multiple-scattering (eta) tables, the in-cloud
extinction alpha (adopted 20 /km, `eta_switch_decision_2026-11-02.md` section 9), followed by an investigation of what
makes the sites differ. It reports measurements only: **no new value is proposed and nothing is applied to the
operational code.** Every number is read from a CSV or markdown table in
`C:/DATA/Projects/202611_eta_study/cloudnet_sites/` (named with each table); if the run is repeated, take the numbers
from those files, not from this text.

## Summary

- **The network median of the Cloudnet extinction equals the adopted 20 /km.** Over the 30 Cloudnet sites inside the
  E-PROFILE domain (the operational CAMS box N27-80, W30-E45; all 30 have at least 20 days), the site medians of
  alpha = 2 pi N_scaled a_G^2 in the lowest 100 m above the liquid base have a median of **19.1 /km** (P25-P75
  13.8-25.6). Day-weighted: 21.0 /km; season-balanced (20 sites): 20.9 /km. Half of the sites (15 of 30) lie below
  20 /km, and 20 of 30 lie inside 10-30 /km, the range that section 9 gives for the backscatter-forming layer. The
  extremes are L'Aquila (6.3 /km) and Rzecin (63.1 /km, a 51-day campaign).
- **Calibration profiles.** At the four step-1d units (the match behind the paper figure), alpha is 33 (Payerne CL61),
  27 (Payerne CL31), 35 (Palaiseau CL31) and 13 /km (Lindenberg CL61), 24 /km pooled over 41 589 profiles. Taken in
  the 41-site framework, the calibration profiles of the nine Vaisala units have 0.91-1.26 times the all-scene alpha
  of their site on the same days. The calibration selection therefore changes alpha by a quarter at most.
- **Height above the base.** The site medians grow from 16.5 /km (lowest 50 m) to 19.1 (100 m) and 22.5 /km (150 m).
  At the median site, alpha grows by a factor 1.31 from the 50-m to the 150-m layer, much less than the factor 2.1 of
  adiabatic growth. The network value lies between the ceilometer-derived 13 /km at the signal peak and 25-35 /km
  100 m higher (section 9).
- **No geographic pattern** (maps): Spearman correlation with |latitude| +0.18 (p 0.36) for alpha, -0.05 for a_G.
  Sites 4-5 km apart can differ more than regions do (section 4.4).
- **What explains the site differences: the retrieval, not the aerosol.**
  - None of 8 aerosol measures from 3 independent sources correlates with a_G, alpha or N_scaled at the 30 sites:
    AERONET AOD, fine-mode AOD, Angstrom exponent and aerosol index; EEA surface PM2.5 and PM10; EBAS in-situ N50 and
    N100. All |rho| are 0.31 or less and all p are 0.12 or more. The retrieval's droplet number N_scaled is 0.18-0.80
    times the in-situ N100 at nine lowland sites, but 4.6 times N100 at Ny-Ålesund and 5.4 times at
    Schneefernerhaus, the two cleanest sites.
  - Cloud-base temperature (CAMS) does not correlate either (rho +0.04 for a_G, +0.06 for alpha).
  - Inside the retrieval, alpha follows the liquid water content that the LWP scaling puts in the lowest 100 m
    (rho = +0.93). a_G follows the split of that water between droplet number and size (rho = -0.70 with N_scaled),
    and is uncorrelated with the water content itself (rho = -0.12).
  - The clearest evidence comes from neighbouring sites observing the same clouds on the same days. Falkenberg, 4.5 km
    from Lindenberg, retrieves the same a_G (same-day ratio 1.09) but 4.3 times the N_scaled and **4.5 times the
    alpha** (67 common days, larger on 99 % of them). Leipzig LIM, 4.1 km from Leipzig, retrieves 2.0 times the
    alpha.
  - Eight of the 30 site medians put more liquid water into the lowest 100 m than adiabatic growth from the Cloudnet
    base allows.
  - The site medians of alpha (and of N_scaled) therefore carry instrument and retrieval differences (radar, LWP,
    lidar-defined layer) of a factor 2-4.5. Those of a_G carry differences of up to about 30 %. These differences are as
    large as the whole between-site spread.
- **Recommendation:** use only the network central values from Cloudnet (alpha about 20 /km, consistent with the
  adopted value). Do not derive per-site or per-unit eta inputs from Cloudnet site medians. Treat the between-site
  spread as retrieval uncertainty, not as regional microphysics.

![Cloudnet extinction: (a) calibration profiles of the four step-1d units, (b) site medians of the 30 E-PROFILE-domain sites](figs_04v2/fig_cloudnet_sites_extinction.png)

## 1. Method

**Same scenes as the droplet note.** `validation/eta_study/cloudnet_sites_extinction.py` re-reads the 17 084 latest
`der` files (41 sites, 2025-01-01..2026-10-31). It imports from `cloudnet_sites_droplet_size.py` the file list, the scene
selection, the liquid-run and layer rules, the reliability flag and the validity masks:

- every 30-s profile whose lowest liquid layer (status 1-3) has its base 0.25-2.4 km above the site;
- the gates whose centre lies in [base edge, base edge + D), D = 50 / 100 / 150 m, truncated at the top of the liquid run;
- only `der_retrieval_status == 1` pixels, with 0 < r_e < 1 mm and N_scaled > 0.

The calibration profiles of the nine Vaisala units are matched as in that note (15 min, liquid within +-150 m of the
CBH). The re-read reproduces the droplet note's per-scene a_G and N_scaled bit for bit at all 41 sites (2 956 918
scenes; `run_extinction.log`). The two parameters are therefore measured on exactly the same scenes.

**Extinction.** Per pixel, the extinction of a lognormal size distribution of width sigma = 0.35 with Q_ext = 2 is

alpha = 2 pi N_scaled r_e^2 exp(-3 sigma^2) = 2 pi N_scaled a_G^2,

with r_e = `der_scaled`. The scene value is the median over the layer's pixels. This is the step-1d definition
(`cloudnet_microphysics_match._layer_stats`, columns `alpha_km_50/100/150` of `step1d/scenes_msfov.csv`). The product
of the site medians, 2 pi N_med a_G,med^2, is within -11..+2 % of the site median of alpha at the 30 E-PROFILE-domain
sites (median -2 %; column `alpha_prodmed`). Site statistics use the same weightings as the droplet note:
scene-weighted median, day-weighted median, and season-balanced value (the median of the four seasonal day-medians,
at least 10 days per season).

**Panel (a)** of the figure uses the step-1d match of the four units, as panel (a) of the paper's droplet figure does.
Panel (b) uses the site medians of `site_table_extinction.csv`.

## 2. The extinction at the calibration cloud base

### 2.1 Across sites

Source: `site_table_extinction.csv`, `across_sites_extinction.csv` (generated table: `note_tables_extinction.md`).

| E-PROFILE-domain sites | sites | median of site medians [P25-P75] (/km) | min - max (/km) | below 10 / 10-30 / above 30 | Spearman with abs(lat) (p) | with altitude (p) |
|---|---|---|---|---|---|---|
| scene-weighted, lowest 100 m | 30 | **19.1** [13.8-25.6] | 6.3 (L'Aquila) - 63.1 (Rzecin) | 5 / 20 / 5 | +0.18 (0.36) | -0.15 (0.44) |
| day-weighted | 30 | 21.0 [13.6-28.1] | 7.1 (Limassol) - 69.0 (Rzecin) | 4 / 20 / 6 | +0.16 (0.40) | -0.06 (0.76) |
| season-balanced (4 seasons) | 20 | 20.9 [15.5-28.4] | 8.8 (Lampedusa) - 55.0 (Norunda) | 1 / 17 / 2 | +0.21 (0.37) | -0.02 (0.95) |
| lowest 50 m | 30 | 16.5 [11.5-22.3] | 5.5 - 60.8 | 6 / 20 / 4 | | |
| lowest 150 m | 30 | 22.5 [17.1-29.0] | 7.1 - 66.6 | 2 / 21 / 7 | | |

- The between-site geometric standard deviation is a factor 1.69 (scene-weighted). Within a site, the scenes spread
  much more: P75 / P25 of a single site is 4.1-21 (e.g. Lindenberg 6.3-26.3, Norunda 12.1-152.3 /km, Graciosa
  3.4-70.8 /km).
- Over all 41 sites the median of the site medians is 20.2 /km (P25-P75 14.4-31.1). The US ARM sites are high: North
  Slope 44.8, Southern Great Plains 49.2, Bankhead 48.8, Phoenix 109 /km (33 days). The marine Southern-Hemisphere
  sites are low: Invercargill 10.7, Tāwhaki 15.7. ATTO-Campina's 2313 /km is not usable: its N_scaled is unphysical,
  as noted for a_G.
- Season: winter (DJF) over summer (JJA), day-weighted, at the 21 E-PROFILE-domain sites with at least 10 days in both
  seasons. alpha: median ratio 0.68, range 0.11 (L'Aquila) to 3.02 (Lampedusa); summer is higher at 16 of 21 sites.
  N_scaled: median ratio 0.65. a_G: median ratio 1.05, winter larger at 13 of 21 sites. So within a site, alpha is
  mostly higher in summer, through a higher N_scaled (`season_ratio_extinction.csv`).

Per site (sites north to south; `note_tables_extinction.md`):

| site | lat | alt (m) | days | scenes | alpha median [P25-P75] (/km) | day-weighted | season-balanced | lowest 50 / 150 m | 2 pi N_med a_G,med^2 | a_G (um) | N_scaled (cm-3) |
|---|---|---|---|---|---|---|---|---|---|---|---|
| Ny-Ålesund (SJ) | 78.9 | 19 | 534 | 96622 | **25.8** [13.1-66.9] | 20.9 | 17.8 | 22.5 / 28.8 | 25.2 | 3.88 | 266 |
| Pallasjärvi (FI) | 68.0 | 275 | 24 | 10090 | **9.2** [4.3-29.6] | 8.8 |  | 5.6 / 11.0 | 8.9 | 3.94 | 91 |
| Kenttärova (FI) | 68.0 | 345 | 320 | 80971 | **14.4** [6.3-50.7] | 23.0 | 21.8 | 11.0 / 17.1 | 13.4 | 3.61 | 164 |
| Sammaltunturinlaatas (FI) | 68.0 | 460 | 25 | 5700 | **8.3** [2.6-27.8] | 10.0 |  | 6.2 / 10.2 | 7.9 | 5.02 | 50 |
| Snæfellsnes (IS) | 64.8 | 14 | 21 | 3109 | **12.5** [6.0-24.5] | 13.5 |  | 11.4 / 14.6 | 12.6 | 4.05 | 123 |
| Hyytiälä (FI) | 61.8 | 150 | 464 | 113238 | **24.8** [10.6-85.8] | 24.1 | 25.0 | 20.2 / 29.0 | 23.0 | 3.37 | 322 |
| Norunda (SE) | 60.1 | 46 | 450 | 117948 | **37.6** [12.1-152.3] | 66.3 | 55.0 | 31.5 / 42.8 | 35.5 | 3.14 | 574 |
| Rzecin (PL) | 52.8 | 57 | 51 | 2449 | **63.1** [27.3-293.8] | 69.0 |  | 60.8 / 66.6 | 62.9 | 3.68 | 738 |
| Lindenberg (DE) | 52.2 | 104 | 552 | 170296 | **13.7** [6.3-26.3] | 13.5 | 13.8 | 11.9 / 17.2 | 13.6 | 3.00 | 241 |
| Falkenberg (DE) | 52.2 | 73 | 94 | 14768 | **41.8** [18.2-151.4] | 55.2 |  | 38.3 / 44.8 | 39.8 | 3.47 | 526 |
| Cabauw (NL) | 52.0 | -1 | 485 | 104212 | **23.5** [10.4-54.2] | 25.4 | 28.2 | 21.1 / 25.6 | 22.7 | 5.52 | 119 |
| Melpitz (DE) | 51.5 | 86 | 205 | 57798 | **15.8** [6.6-33.1] | 13.8 | 14.6 | 14.5 / 18.7 | 15.6 | 3.45 | 209 |
| Leipzig (DE) | 51.4 | 126 | 363 | 104375 | **16.2** [6.7-41.6] | 16.0 | 15.8 | 14.6 / 19.3 | 15.5 | 3.58 | 194 |
| Leipzig LIM (DE) | 51.3 | 126 | 75 | 18779 | **31.1** [14.4-68.2] | 32.5 |  | 28.9 / 34.6 | 31.0 | 3.48 | 406 |
| Chilbolton (GB) | 51.1 | 85 | 49 | 14923 | **37.3** [18.0-80.5] | 56.9 |  | 35.6 / 39.9 | 33.4 | 3.98 | 335 |
| Jülich (DE) | 50.9 | 111 | 537 | 180298 | **18.0** [8.5-36.6] | 17.7 | 17.6 | 13.9 / 20.0 | 17.5 | 3.23 | 268 |
| Palaiseau (FR) | 48.7 | 156 | 489 | 108291 | **28.0** [13.5-57.5] | 29.8 | 29.4 | 24.5 / 31.2 | 27.6 | 3.21 | 426 |
| Kappelrodeck (DE) | 48.6 | 199 | 360 | 81843 | **28.9** [13.3-111.3] | 28.2 | 29.3 | 25.9 / 32.0 | 27.5 | 3.76 | 310 |
| Munich (DE) | 48.1 | 538 | 494 | 118834 | **12.5** [5.3-26.7] | 13.2 | 12.3 | 9.4 / 13.9 | 12.3 | 3.67 | 145 |
| Schneefernerhaus (DE) | 47.4 | 2653 | 136 | 15629 | **23.6** [9.0-66.5] | 21.0 |  | 21.7 / 26.9 | 22.2 | 1.86 | 1020 |
| Payerne (CH) | 46.8 | 491 | 490 | 135862 | **23.4** [11.6-50.8] | 27.7 | 31.6 | 19.7 / 26.8 | 22.9 | 3.65 | 274 |
| Cluj-Napoca (RO) | 46.8 | 459 | 403 | 72431 | **22.8** [11.2-58.9] | 30.2 | 29.1 | 17.5 / 25.0 | 22.1 | 3.64 | 265 |
| Galați (RO) | 45.4 | 40 | 359 | 69906 | **18.0** [4.9-44.3] | 22.1 | 23.5 | 15.4 / 20.3 | 18.0 | 4.13 | 168 |
| Bucharest (RO) | 44.3 | 77 | 427 | 79604 | **20.2** [8.5-52.1] | 20.5 | 20.7 | 17.9 / 23.6 | 19.8 | 3.42 | 270 |
| L'Aquila (IT) | 42.4 | 683 | 489 | 109273 | **6.3** [2.2-22.4] | 10.9 | 10.7 | 5.5 / 7.1 | 6.1 | 5.59 | 31 |
| Potenza (IT) | 40.6 | 760 | 197 | 46654 | **17.0** [6.8-35.9] | 16.7 | 19.0 | 15.5 / 19.7 | 17.1 | 2.90 | 323 |
| Graciosa (PT) | 39.1 | 30 | 240 | 154976 | **17.7** [3.4-70.8] | 18.5 |  | 15.4 / 21.7 | 17.9 | 3.93 | 184 |
| Granada (ES) | 37.2 | 680 | 315 | 66450 | **21.0** [8.5-67.1] | 21.5 | 21.1 | 17.5 / 23.4 | 20.5 | 3.82 | 224 |
| Lampedusa (IT) | 35.5 | 45 | 376 | 99665 | **9.9** [3.3-23.7] | 9.5 | 8.8 | 9.3 / 10.9 | 9.9 | 4.38 | 82 |
| Limassol (CY) | 34.7 | 10 | 221 | 18755 | **8.9** [2.4-23.2] | 7.1 |  | 8.4 / 9.6 | 8.9 | 3.83 | 97 |

### 2.2 Calibration profiles of the Vaisala units

Source: `calib_table_extinction.csv` (41-site framework, 15-min match) and `step1d/scenes_msfov.csv` (step-1d match:
6 min, median-spacing base edge, July file collection; the numbers of section 9). "Site, same days" is the all-scene
median of the site on the days on which the unit has calibration profiles with status-1 data.

| unit | profiles with status-1 data (days) | alpha median [P25-P75] (/km) | lowest 50 / 150 m | site, all scenes | site, same days | calibration / same days | step 1d: 50 / 100 / 150 m (profiles) |
|---|---|---|---|---|---|---|---|
| Payerne CL61 | 2160 (41) | **31.7** [19.3-56.7] | 26.0 / 35.8 | 23.4 | 29.3 | 1.08 | 26.4 / 32.8 / 35.8 (2142) |
| Payerne CL31 (distorted unit) | 16689 (208) | **26.1** [15.2-47.0] | 22.4 / 29.5 | 23.4 | 25.0 | 1.05 | 22.5 / 26.6 / 29.4 (16606) |
| Palaiseau CL31 | 10459 (227) | **35.5** [21.5-54.1] | 30.8 / 39.7 | 28.0 | 29.6 | 1.20 | 30.8 / 35.4 / 39.6 (10401) |
| Lindenberg CL61 (>=2025-05) | 12980 (192) | **15.1** [8.0-24.7] | 12.9 / 19.8 | 13.7 | 13.3 | 1.13 | 10.6 / 12.8 / 16.7 (12440) |
| Chilbolton CL51 | 745 (25) | **57.2** [37.7-102.1] | 55.0 / 60.9 | 37.3 | 45.5 | 1.26 |  |
| L'Aquila CL51 | 20737 (267) | **5.7** [2.3-15.7] | 4.8 / 6.6 | 6.3 | 5.8 | 0.98 |  |
| Norunda CL51 | 11564 (203) | **37.1** [13.9-126.0] | 29.9 / 43.8 | 37.6 | 40.7 | 0.91 |  |
| Middle Wallop CL31 (9 km from Chilbolton) | 1432 (28) | **42.1** [29.0-69.2] | 40.5 / 45.1 | 37.3 | 45.7 | 0.92 |  |
| Akrotiri CL31 (11 km from Limassol) | 110 (14) | **13.6** [5.9-32.8] | 11.9 / 14.9 | 8.9 | 15.0 | 0.91 |  |
| step 1d, four units pooled | | | | | | | 20.2 / 24.0 / 27.8 (41589) |

The step-1d values reproduce section 9 (26 / 33 / 36, 22 / 27 / 29, 31 / 35 / 40, 11 / 13 / 17 /km). The P10-P90 of the
pooled step-1d profiles spans 5.2-108 /km. The two matches differ by less than 2.5 /km at the four units: they differ
in the base-edge rule, the file versions and the time window. The calibration selection keeps alpha within 0.91-1.26 of
the same days' scenes. The two Chilbolton values rest on few profiles (745 and 1 432, see the caveats of the droplet
note).

### 2.3 Against the adopted 20 /km and the section-9 profile

- The adopted 20 /km equals the network median of the site medians (19.1 scene-weighted, 21.0 day-weighted, 20.9
  season-balanced /km) to within 5 %. The calibration units span 5.7-57 /km, and the four step-1d units give 24 /km
  pooled.
- Section 9 measured the extinction profile of the calibration clouds with the ceilometer: about 5 /km at the
  onset, 13 /km at the signal peak 50 m higher, 25-35 /km another 100 m up. The Cloudnet lowest-100-m value of the
  median site (19 /km) lies between the signal-peak and the deeper value, as expected for a 100-m layer starting near
  the base.
- The registration of the Cloudnet layer on the ceilometer profile is uncertain by several tens of metres. Section 9
  places the Cloudnet liquid base at the Vaisala cloud base, 30-60 m above the signal peak. The droplet note measured
  the Cloudnet liquid base 60-95 m below the ceilometer CBH. Over such a distance the extinction changes by a factor
  2-3 (section 9), so the Cloudnet value cannot pin down where in the 10-30 /km range the table value should lie.
- The Cloudnet profile near the base is flatter than adiabatic. The site medians grow from the lowest 50 m to the
  lowest 150 m by a factor 1.31 (range 1.10-1.96 over the 30 sites; `site_table_extinction.csv`, columns `alpha50_med`,
  `alpha150_med`). Adiabatic growth (alpha ~ h^(2/3) at constant N) gives a factor 2.1 between the layer mid-points
  (25 and 75 m). The steeper profile measured by the ceilometer forward inversion in section 9 is not reproduced by
  the radar retrieval, which assumes N constant with height and resolves the base with gates 7.5-36 m deep near 1 km
  (median 30 m; column `dz_1km_m` of `site_table.csv`).

## 3. Maps

![Cloudnet site medians of a_G (a) and alpha (b) on the E-PROFILE box, each with a zoom on central Europe](figs_04v2/fig_cloudnet_sites_maps.png)

The maps show no regional organisation. Low-extinction sites (Lampedusa, Limassol, L'Aquila, the Pallas cluster, Munich,
Lindenberg) and high-extinction sites (Rzecin, Falkenberg, Norunda, Chilbolton, Leipzig LIM) occur side by side; some
of these are campaigns (sampling, section 4.5). The a_G map is dominated by four outliers (Cabauw, L'Aquila,
Sammaltunturinlaatas high; Schneefernerhaus low) on a nearly uniform 3-4 um background. Markers of sites less than 10 km
apart are displaced slightly so that each colour is visible (figure note). The maps are drawn with plotly `scattergeo`
(Natural Earth 1:50 m coastlines and borders); cartopy is not installed.

## 4. What explains the site differences

### 4.1 Candidate explanations and data

`validation/eta_study/cloudnet_sites_explain.py` assembles one row per site (`site_explain_table.csv`).
`validation/eta_study/cloudnet_sites_ebas_n100.py` adds the in-situ particle number. Every external file was
downloaded without login and is cached under `cloudnet_sites/external/`.

- **AERONET** version 3 direct-sun daily averages, level 1.5 (cloud-screened and quality-controlled; level 2.0 is not
  yet available for much of 2025-2026), 2025-01-01..2026-10-31. Source: the AERONET web service
  `https://aeronet.gsfc.nasa.gov/cgi-bin/print_web_data_v3` (`AOD15`, `AOD20`, `SDA15`, `AVG=20`).
  - Site choice: the AERONET site nearest to the Cloudnet site within 100 km. Where several lie within 10 km, the
    longest record is used, and a photometer within 500 m of the site altitude is preferred (so Schneefernerhaus uses
    Zugspitze_KIT, not the valley photometer at Garmisch).
  - Quantities: AOD at 500 nm, Angstrom exponent 440-870 nm, aerosol index AOD x AE (Nakajima et al. 2001), and the
    SDA fine-mode AOD at 500 nm (O'Neill et al. 2003).
  - Coverage: 28 of the 30 E-PROFILE-domain sites have a photometer within 7 km. There is none within 100 km of
    Norunda or Snæfellsnes.
  - Level check: where at least 30 level-2.0 days exist (27 sites), the level-2.0 site medians are 0.78-1.34 times the
    level-1.5 medians (median 1.00, Spearman 0.93; `aeronet_sites.csv`).
  - The header line naming the site principal investigator is not stored.
- **Surface PM2.5 and PM10**: mean of the daily means 2025-01-01..2026-10-31, from the EEA air-quality download service
  (`https://eeadmz1-downloads-api-appservice.azurewebsites.net/`, `ParquetFile/urls`). Verified E1a data cover 2025
  and are completed with the up-to-date E2a data for 2026. Station metadata come from
  `https://discomap.eea.europa.eu/map/fme/metadata/PanEuropean_metadata.csv`.
  - Station choice: the nearest background station within 100 km, within 500 m of the site altitude first, then rural
    or suburban before urban. A station must have at least 40 % of the days.
  - Coverage: 26 sites for PM2.5, 24 for PM10 (`eea_pm_sites.csv`). Schneefernerhaus is excluded: no background
    station within 100 km lies within 500 m of its altitude, and the one selected (AT72710) lies 1 668 m lower.
  - Not available: Ny-Ålesund, Lampedusa and Chilbolton. For Chilbolton, the UK files of the service end in 2019
    (verified) and start in 2026 (up-to-date).
  - The OpenAQ API needs an account key and was not used.
- **In-situ N50 and N100** (dry particles larger than 50 / 100 nm, the usual CCN proxies at stratocumulus
  supersaturations): from the public EBAS THREDDS server (`https://thredds.nilu.no/thredds/catalog/ebas/catalog.html`).
  - Files: level-2 DMPS / SMPS / TDMPS / TSMPS / MPSS size distributions. For each station, the most recent file of
    at least 26 weeks, at most its last 3 years. The integral is taken per hourly spectrum, then the median of the
    hours is used.
  - Matching: the nearest station within 100 km and within 1 000 m of altitude.
  - Data age: the records on that server end between 2010 and 2024, so this is a **station climatology, not the
    2025-2026 period**.
  - Coverage: 11 sites (`ebas_n100_sites.csv`).
  - Pallas (FI0096G, used for the three Finnish Pallas sites) is flagged and not used. Its median hourly dN/dlogD
    maximum is 58 cm-3, below the 214 cm-3 of the Arctic Zeppelin station, which points to another storage convention.
- **Cloud-base temperature**: CAMS temperature at the liquid-base altitude of every scene, from the nearest 3-hourly
  step at the site's 0.4 deg grid point. Sources: the calibration's CAMS model-level files `D:/CAMS_Monthly_04`
  (2025-01..2026-06) and `D:/CAMS_daily` (2026-07..2026-10). The adiabatic liquid-water lapse rate Gamma_ad(T, p)
  follows from it (`cams_tbase_sites.csv`).
- **Cloud and retrieval quantities from the der files** (`site_table_extinction.csv`, `retrieval_consistency_sites.csv`):
  - liquid-layer depth and base height;
  - the share of window scenes whose lowest-100-m pixels carry status 3 (rain flag of the profile) or status 2;
  - the near-base gradient a_G(50 m) / a_G(150 m);
  - the retrieved liquid water content in the lowest 100 m, LWC = 4/3 pi rho_w N r_v^3 with r_v = r_e exp(-sigma^2),
    and its ratio to the adiabatic value at the layer's median height, LWC_ad = Gamma_ad x min(50 m, depth / 2).
  - What status 2 means, from the CloudnetPy source (version 1.95.5 inspected; the files were produced with versions up
    to 1.97): it is not specific to ice. It marks droplet pixels that also carry the categorize "falling" bit
    (`der.py`, `DropletClassification._find_mixed`). That bit is set on every radar echo that is neither clutter nor
    insects. It is removed inside a liquid layer only when Z increases with height and no echo exists directly above
    the lidar-detected liquid top (`falling.py`, `_fix_liquid_dominated_radar`). Status 2 therefore flags drizzle,
    ice, and radar cloud extending above the lidar-detected liquid top alike.
- **Instruments and sampling**: radar band, LWP source (multi-channel radiometer or the 89-GHz channel of the RPG
  radar), site altitude, |latitude|, record length, share of winter scenes.

Statistics: Spearman rank correlations of the site medians of a_G, alpha and N_scaled with every candidate
(`correlations_explain.csv`). The main set is the 30 E-PROFILE-domain sites. Three sensitivity sets drop
Schneefernerhaus; drop L'Aquila and Cabauw; or use all 40 sites with at least 20 days except ATTO-Campina. Two
instrument tests complete the correlations: the same-day comparison of neighbouring sites
(`neighbour_pairs_extinction.csv`) and the sites whose radar changed (`radar_change_extinction.csv`).

### 4.2 Comparison table

Sites sorted by a_G; E-PROFILE domain. Source: `site_explain_table.csv`, generated as `explain_tables.md`. The PM
and AERONET distances are in km. The status columns give the share of window scenes with status-2 / status-3 pixels in
the lowest 100 m. The N100 column is an EBAS climatology of the given years.

| site | alt (m) | a_G (um) | alpha (/km) | N (cm-3) | T_base (C) | depth (m) | status 2 / 3 | LWC/LWC_ad | AERONET site (km) | AOD 500 | AE | fine AOD | PM2.5 (station, km) | PM10 | N100 in situ (station, years) |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| Schneefernerhaus | 2653 | 1.86 | 23.6 | 1020 | -0.4 | 218 | 70 / 2 % | 0.63 | Zugspitze_KIT (0) | 0.027 | 1.10 | 0.020 | 7.8 (AT72710, 18), station -1668 m: not used | 10.9 | 189 (DE0054R, 2016-01..2017-06) |
| Potenza | 760 | 2.90 | 17.0 | 323 | 3.3 | 218 | 60 / 7 % | 0.52 | IMAA_Potenza (0) | 0.099 | 1.30 | 0.067 | 5.4 (IT2283A, 38) | 16.8 |  |
| Lindenberg | 104 | 3.00 | 13.7 | 241 | 2.5 | 249 | 61 / 7 % | 0.42 | MetObs_Lindenberg (0) | 0.105 | 1.27 | 0.076 | 10.0 (DEBB092, 31) | 14.9 |  |
| Norunda | 46 | 3.14 | 37.6 | 574 | 2.8 | 211 | 64 / 5 % | 1.23 |  |  |  |  | 3.4 (E159388, 1) | 7.5 |  |
| Palaiseau | 156 | 3.21 | 28.0 | 426 | 5.5 | 275 | 57 / 11 % | 0.89 | Palaiseau (1) | 0.107 | 1.11 | 0.072 | 8.3 (FR04181, 29) | 12.7 | 753 (FR0020R, 2019-01..2019-12) |
| Jülich | 111 | 3.23 | 18.0 | 268 | 3.7 | 252 | 55 / 7 % | 0.58 | FZJ-JOYCE (0) | 0.122 | 1.21 | 0.088 | 8.7 (DENW096, 27) | 15.9 | 822 (DE0069B, 2018-01..2020-12) |
| Hyytiälä | 150 | 3.37 | 24.8 | 322 | 2.6 | 195 | 64 / 5 % | 0.93 | Hyytiala (1) | 0.091 | 1.31 | 0.072 | 4.4 (FI00801, 46) |  | 403 (FI0050R, 2008-12..2011-12) |
| Bucharest | 77 | 3.42 | 20.2 | 270 | 4.4 | 218 | 70 / 3 % | 0.73 | Magurele_Inoe (2) | 0.152 | 1.47 | 0.121 | 12.5 (RO0071A, 2) | 22.2 |  |
| Melpitz | 86 | 3.45 | 15.8 | 209 | 3.6 | 249 | 60 / 5 % | 0.55 | Melpitz (1) | 0.118 | 1.26 | 0.085 | 7.7 (DESN130, 24) | 11.7 | 1124 (DE0044R, 2010-07..2013-06) |
| Falkenberg | 73 | 3.47 | 41.8 | 526 | 8.5 | 151 | 60 / 7 % | 1.56 | MetObs_Lindenberg (4) | 0.105 | 1.27 | 0.076 | 10.0 (DEBB092, 33) | 14.5 |  |
| Leipzig LIM | 126 | 3.48 | 31.1 | 406 | -1.4 | 238 | 65 / 7 % | 1.30 | Leipzig (4) | 0.119 | 1.23 | 0.079 | 9.5 (DESN059, 7) | 13.9 | 1447 (DE0067K, 2014-11..2017-11) |
| Leipzig | 126 | 3.58 | 16.2 | 194 | 0.9 | 249 | 59 / 8 % | 0.62 | Leipzig (0) | 0.119 | 1.23 | 0.079 | 9.5 (DESN059, 11) | 13.9 | 944 (DE0055B, 2015-01..2017-12) |
| Kenttärova | 345 | 3.61 | 14.4 | 164 | 0.3 | 212 | 67 / 2 % | 0.58 | Matorova_FMI (1) | 0.054 | 1.41 | 0.045 | 1.7 (FI00427, 1) | 3.4 |  |
| Cluj-Napoca | 459 | 3.64 | 22.8 | 265 | 1.5 | 215 | 69 / 5 % | 0.94 | CLUJ_UBB (1) | 0.140 | 1.47 | 0.108 | 15.1 (RO0075A, 1) | 18.6 |  |
| Payerne | 491 | 3.65 | 23.4 | 274 | 2.9 | 187 | 56 / 13 % | 0.93 | Payerne (0) | 0.104 | 1.29 | 0.076 | 7.7 (CH0002R, 0) | 11.7 |  |
| Munich | 538 | 3.67 | 12.5 | 145 | 2.3 | 252 | 62 / 9 % | 0.49 | Munich_University (0) | 0.108 | 1.46 | 0.084 | 8.6 (DEBY089, 6) | 15.2 | 802 (DE0043G, 2010-01..2011-12) |
| Rzecin | 57 | 3.68 | 63.1 | 738 | 8.2 | 250 | 77 / 9 % | 2.24 | POLWET_Rzecin (0) | 0.124 | 1.35 | 0.096 | 16.5 (PL0245A, 54) | 23.0 |  |
| Kappelrodeck | 199 | 3.76 | 28.9 | 310 | 3.1 | 151 | 56 / 12 % | 1.39 | KLOCX_Kappelrodeck (0) | 0.085 | 1.01 | 0.054 | 7.1 (DEBW076, 21) | 11.5 |  |
| Granada | 680 | 3.82 | 21.0 | 224 | 4.4 | 185 | 54 / 12 % | 0.87 | Granada (0) | 0.096 | 1.01 | 0.053 | 10.3 (ES1656A, 71) | 30.1 |  |
| Limassol | 10 | 3.83 | 8.9 | 97 | 12.9 | 187 | 55 / 7 % | 0.28 | CUT-TEPAK (1) | 0.154 | 1.15 | 0.089 | 13.7 (CY0008A, 56) | 28.1 |  |
| Ny-Ålesund | 19 | 3.88 | 25.8 | 266 | -3.0 | 138 | 82 / 4 % | 1.67 | Ny_Alesund_AWI (0) | 0.071 | 1.34 | 0.062 |  |  | 58 (NO0042G, 2008-01..2010-11) |
| Graciosa | 30 | 3.93 | 17.7 | 184 | 10.6 | 210 | 53 / 3 % | 0.52 | ARM_Graciosa (0) | 0.088 | 0.57 | 0.031 | 3.0 (PT07001, 75) | 8.0 |  |
| Pallasjärvi | 275 | 3.94 | 9.2 | 91 | 0.9 | 288 | 56 / 4 % | 0.41 | Matorova_FMI (4) | 0.054 | 1.41 | 0.045 | 1.7 (FI00427, 4) | 3.4 |  |
| Chilbolton | 85 | 3.98 | 37.3 | 335 | 1.7 | 210 | 51 / 11 % | 1.66 | Chilbolton (1) | 0.111 | 1.12 | 0.081 |  |  | 662 (GB0036R, 2010-01..2010-06) |
| Snæfellsnes | 14 | 4.05 | 12.5 | 123 | 4.1 | 238 | 56 / 10 % | 0.48 |  |  |  |  | 3.4 (IS0056A, 99) |  |  |
| Galați | 40 | 4.13 | 18.0 | 168 | 2.6 | 213 | 63 / 7 % | 0.76 | Galati_Rexdan (1) | 0.144 | 1.45 | 0.100 | 10.3 (RO0142A, 1) | 17.6 |  |
| Lampedusa | 45 | 4.38 | 9.9 | 82 | 10.9 | 218 | 39 / 7 % | 0.35 | Lampedusa (0) | 0.143 | 0.70 | 0.066 |  |  |  |
| Sammaltunturinlaatas | 460 | 5.02 | 8.3 | 50 | -0.6 | 224 | 55 / 2 % | 0.43 | Matorova_FMI (7) | 0.054 | 1.41 | 0.045 | 1.5 (FI00356, 1) | 3.3 |  |
| Cabauw | -1 | 5.52 | 23.5 | 119 | 5.0 | 185 | 60 / 9 % | 1.28 | Cabauw (0) | 0.122 | 1.14 | 0.082 | 9.2 (NL00644, 1) | 15.5 | 592 (NL0011R, 2022-02..2024-12) |
| L'Aquila | 683 | 5.59 | 6.3 | 31 | 2.4 | 179 | 61 / 11 % | 0.39 | LAQUILA_Coppito (4) | 0.087 | 1.34 | 0.062 | 6.9 (IT0989A, 36) | 11.1 |  |

### 4.3 Correlations

Spearman rho (p), n; the 30 E-PROFILE-domain sites. Source: `correlations_explain.csv` (`explain_tables.md` also gives
the three sensitivity sets).

| candidate | a_G | alpha | N_scaled |
|---|---|---|---|
| N_scaled (cm-3) | -0.70 (1.6e-05), 30 | +0.87 (5.1e-10), 30 |  |
| AERONET AOD 500 nm | +0.02 (0.91), 28 | +0.07 (0.74), 28 | -0.01 (0.96), 28 |
| AERONET fine-mode AOD 500 nm (SDA) | -0.11 (0.56), 28 | +0.15 (0.43), 28 | +0.12 (0.56), 28 |
| AERONET Angstrom exponent 440-870 | +0.01 (0.98), 28 | -0.21 (0.28), 28 | -0.17 (0.4), 28 |
| AERONET aerosol index AOD x AE | -0.17 (0.4), 28 | +0.12 (0.53), 28 | +0.13 (0.52), 28 |
| surface PM2.5 (ug m-3) | -0.10 (0.64), 26 | +0.31 (0.12), 26 | +0.26 (0.19), 26 |
| surface PM10 (ug m-3) | -0.11 (0.61), 24 | +0.23 (0.28), 24 | +0.22 (0.3), 24 |
| in-situ N100, EBAS climatology (cm-3) | -0.07 (0.83), 11 | -0.27 (0.42), 11 | -0.15 (0.67), 11 |
| in-situ N50, EBAS climatology (cm-3) | -0.13 (0.71), 11 | -0.10 (0.77), 11 | -0.05 (0.89), 11 |
| cloud-base temperature (C, CAMS) | +0.04 (0.83), 30 | +0.06 (0.76), 30 | +0.01 (0.96), 30 |
| adiabatic LWC lapse rate (g m-3 km-1) | +0.10 (0.59), 30 | +0.10 (0.6), 30 | +0.01 (0.97), 30 |
| site altitude (m) | -0.29 (0.12), 30 | -0.15 (0.44), 30 | +0.10 (0.58), 30 |
| absolute latitude | -0.05 (0.78), 30 | +0.17 (0.36), 30 | +0.09 (0.65), 30 |
| liquid-layer depth (m) | -0.29 (0.12), 30 | -0.27 (0.15), 30 | -0.03 (0.88), 30 |
| liquid-base height above site (m) | -0.26 (0.17), 30 | -0.16 (0.4), 30 | -0.10 (0.59), 30 |
| share of scenes with status 3 (rain in profile) | +0.14 (0.46), 30 | +0.21 (0.26), 30 | +0.12 (0.54), 30 |
| share of scenes with status 2 (droplets + falling bit) | -0.39 (0.034), 30 | +0.34 (0.064), 30 | +0.38 (0.037), 30 |
| share of scenes with reliable pixels | +0.24 (0.21), 30 | -0.44 (0.014), 30 | -0.41 (0.025), 30 |
| a_G(50 m) / a_G(150 m) | +0.09 (0.64), 30 | +0.42 (0.02), 30 | +0.35 (0.062), 30 |
| retrieved LWC in the lowest 100 m (g m-3) | -0.12 (0.54), 30 | +0.93 (5.6e-14), 30 | +0.67 (5e-05), 30 |
| retrieved / adiabatic LWC in the lowest 100 m | -0.16 (0.41), 30 | +0.94 (1.1e-14), 30 | +0.71 (1.1e-05), 30 |
| W-band radar (1) vs Ka (0) | +0.28 (0.13), 30 | +0.31 (0.099), 30 | +0.04 (0.84), 30 |
| LWP from radiometer (1) vs radar 89-GHz channel (0) | -0.15 (0.44), 30 | +0.25 (0.18), 30 | +0.25 (0.18), 30 |
| days of data | -0.21 (0.26), 30 | +0.00 (0.98), 30 | -0.02 (0.91), 30 |
| share of scenes in DJF | +0.06 (0.76), 30 | +0.11 (0.56), 30 | +0.04 (0.84), 30 |

The table holds 72 tests, so about 3-4 of them would fall below p = 0.05 by chance. The status-share correlations
(p 0.014-0.037) are at that level. The N_scaled and LWC correlations (p < 1e-4) are not. No aerosol measure reaches
p < 0.1 in any of the four sets: in the sensitivity sets the largest |rho| of an aerosol measure is 0.41 (N50, 10
sites, p 0.24), and the smallest p is 0.12.

![Site medians of a_G (top) and alpha (bottom) against N_scaled, the AERONET fine-mode AOD and the retrieved-to-adiabatic LWC ratio](figs_04v2/fig_cloudnet_sites_vs_aerosol.png)

### 4.4 Same clouds, different instruments

Source: `neighbour_pairs_extinction.csv`. Common days are days with at least 20 scenes at both sites. The table gives
the median over the common days of the same-day ratio (second site over first) and the share of those days on which
the second site is larger.

| pair (distance) | common days | instruments (radar / lidar / LWP) | a_G ratio | N_scaled ratio | alpha ratio (days larger) |
|---|---|---|---|---|---|
| Lindenberg -> Falkenberg (4.5 km) | 67 (2026-05-06..09-01) | MIRA-35 (Ka) -> RPG-FMCW-94 (W); CHM15k at both; HATPRO at both | 1.09 | 4.30 | **4.46** (99 %) |
| Leipzig -> Leipzig LIM (4.1 km) | 60 (2025-01-01..2026-04-20) | MIRA-35 or RPG-FMCW-94 -> RPG-FMCW-94; PollyXT and other lidars -> CHM15k; HATPRO or the 89-GHz channel at both | 0.96 | 2.33 | **2.02** (77 %) |
| Kenttärova -> Pallasjärvi (4.7 km) | 21 (2026-09-20..10-21) | RPG-FMCW-94 at both; CL61 -> CL31; 89 GHz -> HATPRO | 1.07 | 0.84 | 0.73 (43 %) |
| Kenttärova -> Sammaltunturinlaatas (6.3 km) | 17 (2026-09-21..10-16) | RPG-FMCW-94 at both; CL61 at both; 89 GHz at both | 1.18 | 0.70 | 0.90 (47 %) |
| Pallasjärvi -> Sammaltunturinlaatas (6.1 km) | 18 (2026-09-21..10-19) | RPG-FMCW-94 at both; CL31 -> CL61; HATPRO -> 89 GHz | 1.16 | 0.52 | 0.75 (44 %) |

Two sites 4-5 km apart see the same air mass and the same cloud field on the same days, so a systematic same-day
difference cannot come from the aerosol or the cloud climate. Falkenberg retrieves 4.5 times Lindenberg's extinction
on 99 % of 67 common days, with the same droplet size. A factor 4.5 is the spread between the low and the high sites of
the whole network (e.g. Lampedusa 9.9 and Norunda 37.6 /km). In the Pallas cluster (same radar type, autumn 2026)
a_G agrees within 7-18 % and alpha within 10-27 %, while N_scaled differs by up to a factor 2.

Where a site changed its radar (`radar_change_extinction.csv`), a pure reflectivity offset k would multiply a_G by
k^(1/3), N_scaled by 1/k and alpha by k^(-1/3).

- Two changes follow that signature closely. At Bucharest the W-band period gives a_G x1.30, N_scaled x0.48 (offset
  signature 0.46) and alpha x0.80 (0.77), equivalent to 3.4 dB. At Invercargill the W-band period gives a_G x1.22,
  N_scaled x0.60 (0.56) and alpha x0.87 (0.82), equivalent to 2.5 dB.
- Others do not. Leipzig's W-band period gives N_scaled x1.43 and alpha x1.65 against 0.82 / 0.94 for an offset;
  Lindenberg's (January-February only) gives 1.01 / 1.37 against 0.66 / 0.87.
- The radar periods differ in season, so this is suggestive, not conclusive.

### 4.5 Reading, hypothesis by hypothesis

- **Aerosol (Twomey): not supported.**
  - Expectation: more CCN should give a higher N, a smaller a_G and a somewhat higher alpha.
  - The site medians of N_scaled, a_G and alpha do not follow any of the eight aerosol measures.
  - The aerosol contrast between the sites is large: fine-mode AOD 0.020 (above Zugspitze) to 0.121 (Bucharest),
    PM2.5 1.5-16.5 ug m-3, in-situ N100 58-1 447 cm-3.
  - The retrieval's N_scaled does not track the in-situ particle number, and at the clean end it exceeds it.
    Ny-Ålesund retrieves 266 cm-3 droplets against 58 cm-3 particles above 100 nm (83 above 50 nm). Schneefernerhaus
    retrieves 1 020 cm-3 against 189 (396 above 50 nm). At the polluted continental sites N_scaled is 0.18-0.28 times
    N100 (Melpitz, Leipzig, Munich, Cabauw, Leipzig LIM). A droplet number above the number of particles larger than
    50 nm is not physical.
  - The EBAS records are a climatology of other years (2008-2024) and of the station level. That can shift N100 by
    tens of percent, not by the factor 16-30 that separates N_scaled / N100 at the clean and the polluted sites.
- **Cloud-base temperature and adiabatic water content: not supported.**
  - The site medians of the cloud-base temperature range from -3.0 C (Ny-Ålesund) to 12.9 C (Limassol), and Gamma_ad
    from 1.39 to 2.16 g m-3 km-1. Neither correlates with a_G, alpha or N_scaled.
  - The warm Mediterranean and marine sites (Limassol, Lampedusa, Graciosa) have low retrieved-to-adiabatic LWC
    ratios (0.28, 0.35, 0.52, in a network range of 0.28-2.24), not high ones.
- **Mountain and cold sites: supported only at the mountain-top station, and not through temperature or thin clouds.**
  - Schneefernerhaus is the a_G outlier (1.86 um). Its clouds are neither colder (T_base -0.4 C; Ny-Ålesund -3.0 C
    has 3.88 um) nor thinner (218 m liquid depth; the site medians span 138-288 m, median 216 m) than elsewhere.
  - Its small a_G comes with the highest N_scaled of the domain (1 020 cm-3), 5.4 times the in-situ N100, at an
    ordinary water content (LWC/LWC_ad 0.63). This means a low reflectivity for the retrieved water. That can be
    genuinely small droplets (the in-situ N100 does not support this), a reflectivity bias of the radar, or the
    mixed-phase environment (status 2 in 70 % of the scenes). These data cannot separate them.
  - Altitude alone correlates weakly: rho -0.29 (p 0.12) for a_G, -0.22 (p 0.26) without Schneefernerhaus.
- **Drizzle inflating r_e: not supported in the reliable sample.**
  - The rain flag (status 3; 2-13 % of the scenes) does not correlate with a_G (+0.14).
  - Status 2 (droplets plus falling bit; 39-82 % of the scenes) correlates negatively with a_G (-0.39, p 0.034). This
    is the opposite of drizzle-inflated radii.
  - Pixels with a falling bit are excluded from the reliable sample, so sites with frequent drizzle keep their less
    drizzly scenes.
- **Retrieval and instruments: supported; this is what the data show.**
  - How the scaled Frisch retrieval works: CloudnetPy's `der_scaled` distributes the column LWP of the microwave
    radiometer over each lidar-detected liquid layer in proportion to sqrt(Z). Hence N_scaled ~ LWP^2 / (integral of
    sqrt(Z) dz)^2 and r_e ~ Z^(1/6) (integral of sqrt(Z) dz)^(1/3) / LWP^(1/3) (`der.py`; Frisch et al. 2002).
  - Water content: the retrieved LWC is set by the LWP and the lidar-defined layer, not by the radar calibration. alpha
    follows it: rho = +0.93 with LWC and +0.94 with LWC/LWC_ad. alpha ~ LWP^(4/3), so a 20 % LWP error moves alpha by
    about 27 %.
  - The retrieved water content spans 0.027-0.203 g m-3, against an adiabatic 0.061-0.105 g m-3. Eight site medians
    are super-adiabatic: Rzecin 2.24, Ny-Ålesund 1.67, Chilbolton 1.66, Falkenberg 1.56, Kappelrodeck 1.39, Leipzig LIM
    1.30, Cabauw 1.28, Norunda 1.23. Making them adiabatic would need a condensation level 12-62 m below the Cloudnet
    liquid base (`retrieval_consistency_sites.csv`). At these sites 56-77 % of the scenes are super-adiabatic.
  - Partition between number and size: at a given LWC, the split is set by the reflectivity. A reflectivity offset k
    moves a_G as k^(1/3) and N_scaled as 1/k, with LWC unchanged. 1 dB is therefore 8 % in a_G, the size of the
    interquartile half-width of the site medians of a_G (3.42-3.94 um around 3.66, droplet note).
  - The a_G-N_scaled anticorrelation (rho -0.70; panel (a), close to a_G ~ N^(-1/3) at constant water content) and the
    absence of any a_G-LWC relation (rho -0.12) are what such offsets, or real number differences, produce. The aerosol
    data give no support to real number differences.
  - The neighbour pairs (4.4) show directly that changing instruments changes N_scaled and alpha by factors 2-4 under
    the same clouds.
- **Radar band and LWP source: no significant site-level difference.** W-band sites have higher alpha (rho +0.31,
  p 0.10). The LWP source shows rho +0.25 (p 0.18) for alpha. Instrument effects act through each unit's calibration
  and set-up, not through the band as a class.
- **Sampling: no significant effect on the site medians.** Record length gives rho -0.21 for a_G and 0.00 for alpha;
  the share of winter scenes +0.06 and +0.11. But the within-site seasonal cycle is as large as the between-site
  spread (alpha DJF/JJA 0.11-3.02). This matters for the short records: Pallasjärvi and Sammaltunturinlaatas
  (autumn 2026), Snæfellsnes and Rzecin (July-October 2026), Falkenberg (May-August 2026), Chilbolton (one winter);
  Leipzig LIM has 75 days spread over 16 months.

### 4.6 Conclusion

The between-site differences of the Cloudnet droplet radius and extinction at the calibration cloud base are not
explained by the aerosol load, by the cloud-base temperature, by latitude, or (apart from the mountain-top station) by
altitude. They are explained by the retrieval and its instruments:

- alpha is set by the liquid water that the LWP scaling places in the lowest 100 m;
- a_G is set by how the radar reflectivity splits that water into number and size;
- sites 4-5 km apart that see the same clouds differ by factors 2-4.5 in alpha and N_scaled.

The site medians are therefore not regional microphysics. Only the network central values carry physical
information: alpha about 20 /km, a_G about 3.7 um. The between-site spread measures the retrieval and calibration
uncertainty of the Cloudnet products (a factor about 1.7 in alpha, about 20 % in a_G).

**Recommendation:** keep one network value for each table input, and do not derive per-site or per-unit eta inputs from
Cloudnet site medians. The adopted alpha = 20 /km is consistent with the network median (19-21 /km by any weighting).

## 5. Caveats

- **Extinction from the radar retrieval.** alpha is the 2 pi N a_G^2 of the Frisch scaled retrieval: lognormal of fixed
  width 0.35, N constant with height, Z ~ N D^6. It is not an independent optical measurement. Any error of N_scaled
  enters linearly.
- **Reliable-pixel selection.** Status 1 requires the falling bit to be removed (radar echo ending at the lidar-detected
  liquid top, Z increasing with height). Only 14-55 % of the window scenes have such pixels in the lowest 100 m, and the
  share correlates with alpha (-0.44). The sample is not all liquid clouds.
- **Adiabatic index.** LWC/LWC_ad uses the layer median height (50 m), the Cloudnet liquid base (lidar-detected) and
  CAMS temperatures at a 0.4 deg grid point. A true condensation level 12-33 m below the Cloudnet base would bring the
  1.2-1.7 values to 1, and 62 m the 2.24 of Rzecin (`base_offset_for_adiabatic_m`). The ranking between sites is the
  robust part, not the threshold.
- **Aerosol proxies.**
  - AOD is a column quantity and includes elevated dust and smoke. The fine-mode AOD and the aerosol index reduce, but
    do not remove, that problem: AOD is known to be a weak proxy of cloud-base CCN (Stier 2016, ACP 16, 6595-6607).
  - Surface PM is mass, not number, and some stations are urban or up to 100 km away (`eea_pm_sites.csv`).
  - The EBAS N100 is a climatology of other years.
  - None of the three is a measurement at cloud base. A weak aerosol effect of the size of the Twomey sensitivity
    could be hidden by the retrieval spread. The data exclude only that aerosol drives the observed differences.
  - Published satellite droplet-number and CCN climatologies (Bennartz and Rausch 2017, ACP; Grosvenor et al. 2018, Rev.
    Geophys.; Fanourgakis et al. 2019, ACP) and the European size-distribution climatology of Asmi et al. (2011, ACP) were
    not digitised. No number in this note comes from them.
- **Neighbour pairs.** The Lindenberg-Falkenberg and Leipzig pairs use different radars and lidars and different
  radiometer units. Which component causes the factor 2-4.5 cannot be isolated with the der files alone: that needs the
  categorize files (LWP, Z).
- **Small samples:** Rzecin (51 days), Falkenberg (94 days in 4 months), Pallasjärvi, Sammaltunturinlaatas, Snæfellsnes
  (21-25 days) and Chilbolton (one winter) are campaigns. The neighbour pairs in the Pallas cluster have 17-21 common days.

## 6. Files and data sources

Code (repository):

- `validation/eta_study/cloudnet_sites_extinction.py`: re-read of the der files (16 processes, about 4 min), per-scene
  extinction, site and unit tables, the extinction figure and the maps.
- `validation/eta_study/cloudnet_sites_explain.py`: AERONET and EEA downloads (cached), CAMS cloud-base temperature,
  retrieval-consistency index, neighbour pairs, radar changes, correlations, comparison figure; `--offline` re-runs
  from the cache.
- `validation/eta_study/cloudnet_sites_ebas_n100.py`: EBAS size-distribution inventory and N50 / N100.
- The scene selection and constants are imported from `validation/eta_study/cloudnet_sites_droplet_size.py`
  (unchanged).

Outputs (`C:/DATA/Projects/202611_eta_study/cloudnet_sites/`):

- extinction: `extinction_scenes/<site>.npz` and `<site>_days.csv` (per scene alpha 50 / 100 / 150 m, a_G, N, depth,
  base; per-day status counts), `site_table_extinction.csv`, `season_table_extinction.csv`,
  `calib_scenes_extinction.csv`, `calib_table_extinction.csv`, `across_sites_extinction.csv`,
  `note_tables_extinction.md`, `run_extinction.log` (identity check against the droplet cache);
- explanation: `site_explain_table.csv`, `correlations_explain.csv`, `explain_tables.md`, `aeronet_sites.csv`,
  `eea_pm_sites.csv`, `ebas_n100_sites.csv`, `cams_tbase_sites.csv`, `cams_tbase/<site>.npz`,
  `retrieval_consistency_sites.csv`, `season_ratio_extinction.csv`, `radar_change_extinction.csv`,
  `neighbour_pairs_extinction.csv`;
- raw external data: `external/aeronet/` (daily AOD15 / AOD20 / SDA15 files, location list), `external/eea/` (station
  metadata, URL lists, parquet files), `external/ebas/` (catalogue, station list, hourly N50 / N100 per station).

Figures: `doc/reports/figs_04v2/fig_cloudnet_sites_extinction.png`, `fig_cloudnet_sites_maps.png`,
`fig_cloudnet_sites_vs_aerosol.png`.

Data sources and references:

- ACTRIS Cloudnet `der` product (CC-BY 4.0), CloudnetPy (Tukiainen et al. 2020, JOSS 5, 2123,
  https://doi.org/10.21105/joss.02123). Retrieval: Frisch et al. 2002, J. Atmos. Oceanic Technol. 19, 835-842,
  https://doi.org/10.1175/1520-0426(2002)019%3C0835:TROSCD%3E2.0.CO;2. The status definitions were read in the CloudnetPy
  source, `cloudnetpy/products/der.py` and `cloudnetpy/categorize/falling.py`.
- AERONET version 3, https://aeronet.gsfc.nasa.gov/ (web service `cgi-bin/print_web_data_v3`, site list
  `aeronet_locations_v3.txt`); Giles et al. 2019, Atmos. Meas. Tech. 12, 169-209 (version 3); O'Neill et al. 2003, J.
  Geophys. Res. 108(D17), 4559 (spectral deconvolution); Nakajima et al. 2001, Geophys. Res. Lett. 28, 1171-1174
  (aerosol index).
- European Environment Agency air-quality download service, https://eeadmz1-downloads-api-appservice.azurewebsites.net/
  (E1a verified, E2a up-to-date) and station metadata
  https://discomap.eea.europa.eu/map/fme/metadata/PanEuropean_metadata.csv.
- EBAS (NILU), public THREDDS catalogue https://thredds.nilu.no/thredds/catalog/ebas/catalog.html (level-2 particle
  number size distributions).
- CAMS global atmospheric composition model levels (the calibration's water-vapour input, local copies).
- Twomey 1977, J. Atmos. Sci. 34, 1149-1152 (aerosol effect on droplet number and size); Bolton 1980, Mon. Weather Rev.
  108, 1046-1053 (saturation vapour pressure in Gamma_ad); Stier 2016, Atmos. Chem. Phys. 16, 6595-6607 (AOD as a CCN
  proxy).
- Published climatologies named for completeness, not used quantitatively: Asmi et al. 2011 (ACP), Bennartz and Rausch
  2017 (ACP), Grosvenor et al. 2018 (Rev. Geophys.), Fanourgakis et al. 2019 (ACP).
