# In-cloud extinction α and the other inputs of the multiple-scattering factor η in the literature (2026-11)

Literature check of 2026-11-02, written for the η study (`04_multiple_scattering_v2.md`,
`eta_switch_decision_2026-11-02.md`). It establishes what the published liquid-cloud ceilometer
calibrations assumed when they computed η (B = 1/(2ηS)): in-cloud extinction α, liquid water content and
droplet number, droplet size, receiver field of view (FOV) and laser divergence and their convention, lidar
ratio S, the model, and how η depends on cloud-base height (CBH) and instrument. It then compares them with
the inputs now in the code (`calibration/cloud/_filters.py`, generation v2: PVC, a_G = 3.3 µm, α = 30 km⁻¹,
S = 18.8 sr). **It does not propose a value for α.** Section 5 runs the literature inputs through our PVC port;
that part is an inference, not a literature statement, and is labelled as such.

Quotations: one short verbatim quotation is given (§2.1). Every other statement is a close paraphrase with its
page, section, equation, table or figure, so that it can be checked in the documents listed in §8 (most are
in `doc/`).

## 1. Answer in brief

| question | answer |
|---|---|
| α in O'Connor et al. (2004) | **15 and 20 km⁻¹, homogeneous (constant with depth), assumed**: chosen to span typical stratocumulus (in-situ climatology of Miles et al. 2000) and to ensure full attenuation. Not measured on the calibration scenes, not an adiabatic profile. |
| α in Hopkin et al. (2019) | **not stated**. η is computed per gate with the Hogan (2006) code; no cloud or optics inputs are given anywhere in the paper or its review. Their operational code applies one η ladder (0.829 → 0.779 over 0.25–1.875 km) to all Vaisala types, taken from the University of Reading calibration page. |
| α in Le et al. (2026), CL61 | **not stated**. Only droplet effective diameters 8–20 µm and the CL61 optics (divergence ±0.2 × ±0.35 mrad, FOV ±0.56 mrad) are given; the MS model and the numerical S are not named. |
| α in Hogan (2006) | none prescribed: α(r) and a_G(r) are model inputs. Test cases are ice clouds; Hogan (2008) tests a 3 km⁻¹ liquid layer over a wide range of receiver footprints. |
| Cloudnet / CloudnetPy | no liquid-cloud calibration, no η, no α, no S in CloudnetPy: it multiplies by a calibration factor supplied per instrument and date. The Reading web tool (Hogan's code, CT75K defaults 905 nm, divergence 0.75, FOV 0.66 mrad) shows the O'Connor cloud set (D 8 / 20 µm × α 15 / 20 km⁻¹); its calculator now returns a server error (HTTP 500), so its η cannot be retrieved. |
| adiabatic profiles | **no source uses one for η.** O'Connor's Fig. 1 uses a linearly increasing LWC only to illustrate single-scattering profiles. |
| independent observations and climatologies (§2.9) | no long-term climatology of the extinction near the base of boundary-layer stratocumulus exists. Dual-FOV Raman lidar over Leipzig (29 thin altocumulus): **11.5 / 19.4 / 25.5 km⁻¹ at 0–30 / 30–70 / 70–120 m above base**. CALIPSO + MODIS cloud tops over the oceans: mostly 15–30 km⁻¹ (up to 50 off polluted coasts). Ground-based τ/H for continental stratus (ARM SGP): 21–30 km⁻¹, a whole-cloud mean. Adiabatic theory with typical continental droplet numbers (150–300 cm⁻³): mean over the layer a ground lidar sees (τ ≤ 3) **≈ 24–32 km⁻¹**. |
| our choices vs the literature | α 30 km⁻¹ is **measured** and 1.5–2× O'Connor's assumed 15–20 km⁻¹; a_G 3.3 µm is below O'Connor's smallest size (4 µm) and at the small end of Le's band; S 18.8 sr is the value of every source that gives one; our FOV values are the documented half-angles (Wiegner et al. 2014 Table 1, CL61 User Guide), the same numbers as Le et al. Table 1. |
| consequence (our model) | using 15–20 instead of 30 km⁻¹ (a_G 3.3 µm, our optics) would raise η, hence C_L, by +3.3…+6.7 % at CBH 0.5 km, +5.0…+9.7 % at 1 km, +6.7…+12.2 % at 2 km (§4.1). |
| cross-check (inference) | the published η of O'Connor (2004) and the ηS envelope of Le et al. (2026) are reproduced by our PVC port with their own cloud inputs **only if the FOV and divergence are halved**; with the documented optics they correspond to α ≈ 5–10 km⁻¹ (Le) or remain 0.07–0.12 above our η (O'Connor). Our v2 CL61 table runs 0.008–0.010 below the lower edge of the Le envelope, parallel to it (§5). |

## 2. What each source states

### 2.1 O'Connor et al. (2004), J. Atmos. Oceanic Technol. 21, 777–786

| item | where | statement (paraphrase unless quoted) |
|---|---|---|
| definition of η | Sect. 2, Eqs. (3)–(6), p. 778 | Platt's factor η (between 0.5 and 1) multiplies the optical-depth increment; for a totally attenuating cloud with S and η constant along the path, B = 1/(2ηS). η assumed constant with penetration depth, supported by Kunkel and Weinman (1976) and Platt (1981) for optically thick water clouds. |
| model | Sect. 2, p. 780 | Eloranta's (1998) multiple-scattering model is used to compute the increase of B for typical water clouds as a function of cloud height for a given lidar. The code is credited in the acknowledgements (p. 786). Scattering order not stated. |
| **α** | Sect. 2, p. 780 | "extinction coefficients of 15 and 20 km−1 were chosen to describe the range" (O'Connor et al. 2004, p. 780); the range is the one of typical stratocumulus (Miles et al. 2000) that also guarantees full attenuation. |
| α constant with depth | Sect. 3, Fig. 5 text, p. 781 | the theoretical lines are for clouds of 15 and 20 km⁻¹, which reach optical depth 1 at 66 and 50 m: homogeneous clouds. |
| droplet size | Sect. 2, p. 780; p. 781 | mean droplet effective diameters 8 and 20 µm, the extremes of Miles et al. (2000) for clouds thick enough to attenuate fully (8 µm continental / supercooled, ~20 µm marine); monodisperse in Fig. 5. In PVC terms a_G = 4 and 10 µm. Droplet size is named as the main uncertainty of η. |
| LWC, N | Sect. 2, Fig. 1, p. 778 | not inputs of η. Fig. 1 (no multiple scattering) illustrates B with LWC increasing linearly with height: 1 g m⁻³ km⁻¹ with N = 2000 cm⁻³, and 0.25 g m⁻³ km⁻¹ with N = 50 cm⁻³, monodisperse. |
| optics | Sect. 2, p. 780; Sect. 3, p. 780; abstract, p. 777 | Vaisala CT75K, 905 nm. Sect. 2: half-angle divergence 0.75 mrad and receiver FOV 0.66 mrad. Sect. 3: divergence 0.66 and FOV 0.75 mrad, both half-angles (swapped). Abstract: beam divergence of 1–1.5 mrad (full angle). Assumptions: Gaussian beam, full overlap, parameters known (p. 780). |
| S | abstract, p. 777; Sect. 2, Eq. (7), Figs. 2–3, pp. 778–779 | 18.8 ± 0.8 sr at 905 nm: Mie, water n = 1.327 + 0.672×10⁻⁶ i (25 °C), normalized gamma distributions with median volume diameter D₀ 8–20 µm and shape µ 2–10 (Miles et al. 2000); mean 18.8 sr over D₀ 10–50 µm. Below D₀ ≈ 8 µm S begins to rise at 905 nm for broad distributions (Sect. 5, p. 783). |
| CBH and instrument dependence | Sect. 2, p. 780 | MS grows with range as the beam widens; η must be recomputed for each lidar (wavelength and FOV in particular). |
| η values | p. 780; p. 782 | CT75K: η = 0.83 ± 0.09 at 1 km (ηS 14.5–17 sr) and 0.73 ± 0.06 at 4 km (ηS 13–14.5 sr); 0.79 at 1.5 km (ηS ≈ 15 sr). At 2.7 km: 0.74 printed as ± 0.4, evidently ± 0.04 given the quoted ηS 13.1–14.7 sr; the 354.7 nm Chilbolton UV lidar (divergence 0.05, FOV 0.2 mrad, half-angles) 0.86 ± 0.06. |
| how η is used | Sect. 3, Fig. 5, p. 781 | the calibration scales β so that the per-height median of the apparent ηS lies inside the theoretical envelope of the four cloud cases (medians within 10 % of 15.5 sr): an envelope, not one η curve. |

### 2.2 University of Reading lidar-calibration page (companion tool of O'Connor et al. 2004)

`http://www.met.reading.ac.uk/~swr99ejo/lidar_calibration/index.html` (page still online on 2026-11-02; HTML and
example figure read). **The calculator no longer works**: on 2026-11-02 its script
(`/~swr99ejo/cgi-bin/lidar_calibration.cgi`) returned HTTP 500 (Internal Server Error) for the CT75K defaults,
submitted both from the app's browser and from the command line (HTTP and HTTPS), and also for a request without
any parameters, so the fault is in the server script, not in the inputs. The CL61 and CL31 inputs were therefore
not submitted. The Internet Archive holds only the page and `example.png` (2024-07-03), no calculator output.

- Computes ηS versus height with **Hogan's multiscatter code** (not Eloranta's). Inputs: wavelength, divergence
  (form field `rho_l`), telescope FOV (`rho_t`), lidar altitude. Defaults, given as a CT75K example at sea level:
  905 nm, `rho_l` 0.75 mrad, `rho_t` 0.66 mrad, altitude 0 m. The field names are Hogan's model symbols (§2.4),
  and the defaults match Sect. 2 of the paper, not Sect. 3.
- The example figure overlays CT75K observations with five theoretical lines: D = 20 µm or 8 µm, each with
  σ = 20 or 15 km⁻¹, plus single scattering: the cloud set of the paper.
- S = 18.8 ± 0.8 sr at 905 nm for stratocumulus without drizzle; for a 1–1.5 mrad beam MS adds ~10 % calibration
  uncertainty.
- It is the source cited by the E-PROFILE code of Hopkin et al. (2019) for its single η ladder (§2.3).

### 2.3 Hopkin et al. (2019), Atmos. Meas. Tech. 12, 4131–4147

| item | where | statement |
|---|---|---|
| definition and model | Sect. 3, Eq. (1), p. 4134 | B = 1/(2ηS); η depends on wavelength, beam divergence, telescope FOV and altitude; computed for each gate with the fast method and code of Hogan (2006); usually 0.7–0.85 for 905–1064 nm in liquid clouds. Target B = 0.0266 sr⁻¹ (printed as m⁻¹) for S = 18.8 sr. |
| S | Sect. 3, p. 4134 | 18.8 sr for 905/910 and 1064 nm, following O'Connor et al. (2004). |
| α, droplet size, LWC, N | — | **not stated** anywhere (text, tables, figures). |
| FOV, divergence | Table 1, p. 4133 | not listed (Table 1 gives wavelength, optical design, pulse energy, PRF, resolution, overlap, range only). |
| review | amt-2018-427 RC1, RC2, AC1, AC2 | no reviewer asked for the η inputs; the replies only add that Hogan's code gives negligible MS in drizzle and that near-range η may be wrong under imperfect telescope alignment. |
| check of η | Sect. 4, Fig. 6, p. 4138 | after the η and water-vapour corrections, B at Middle Wallop (CL31, >100 000 profiles) is 0.021 sr⁻¹ in 16 of 21 100-m height bins; below 500 m B is 9.5 % lower (saturation or near-range η suspected); clouds at 200–500 m are rejected. |
| operational code | E-PROFILE Python of the study, `VaisCeil_DayCalibration_EProf_EH.py`, function `scatter_correct_Vais` (copy in the MeteoSwiss `alc_calib_py` repository) | one ladder for **every** Vaisala type, 0.25-km steps: 0.82881 … 0.77877 from 0.25 to 1.875 km (nothing above); source given as the Reading page (§2.2); a code comment notes the values are instrument dependent. These are the legacy "CL51" rows of `_filters.py`. |

The cloud and optics inputs behind the ladder are unknown. Our port does not reproduce it with the CT75K
defaults of the Reading page either (mid of the O'Connor cloud set: 0.885 → 0.688 over the same heights,
against 0.829 → 0.779), consistent with `04_multiple_scattering.md` §3b. A PhD thesis cited in that report was
not located in this check.

### 2.4 Hogan (2006), Appl. Opt. 45, 5984–5992, and Hogan (2008), J. Atmos. Sci. 65, 3621–3635

- **Conventions** (Hogan 2006, Sect. 2): the laser beam is Gaussian with a **1/e angular half-width ρ_l** of the
  energy density (Eq. 1, p. 5985); the telescope has a **half-angle FOV ρ_t** and accepts photons whose lateral
  offset is at most ρ_t R, i.e. a sharp-edged cone (p. 5986); the forward diffraction lobe is a Gaussian of
  1/e half-width Θ = λ/(π a_G), with a_G the **equivalent-area radius** of the size distribution (Eq. 5, p. 5986).
  Eloranta's model chose the lobe width to preserve peak and energy instead of energy and variance (p. 5986).
- **α**: an input profile; no default and no adiabatic profile for liquid water. Validation: a 4-km ice cloud
  (α 2 / 0.5 / 0.1 / 0.5 km⁻¹ in 1-km layers, a_G 100 µm, S 20 sr) seen by a ground-based 532 nm lidar with
  ρ_l 0.5 and ρ_t 0.75 mrad: agreement with high-order Eloranta within 4 % (p. 5989), equally good for divergence
  and FOV from 0.005 to 50 mrad (p. 5990). Optically thick liquid clouds with large footprints are outside the
  validity of both models (Sect. 6, p. 5991).
- Hogan (2008): the PVC family; a test liquid layer at 1–1.5 km with constant α = 3 km⁻¹ and a_G = 10 µm, with an
  ice cloud at 3.5–6 km above it, evaluated over a wide range of receiver footprints (Sect. 4, p. 3630); a **Gaussian receiver of 1/e half-width** as an
  alternative to the top-hat for diffraction-limited receivers (Eq. 28, p. 3628); in liquid cloud the anisotropic
  near-180° phase function lowers the apparent backscatter by up to 50 % (p. 3633), the error exceeding 50 % for
  footprints above ~20 m (p. 3632). A ceilometer footprint at 1 km is 2 R ρ_t ≈ 1–2 m (our arithmetic).

### 2.5 Cloudnet / CloudnetPy

- CloudnetPy (`actris-cloudnet/cloudnetpy`, `main` at a359dc4, 2026-10-22) holds **no liquid-cloud calibration and
  no multiple-scattering code** (repository code search: no hits). `instruments/ceilo.py::ceilo2nc` multiplies β by
  a `calibration_factor` from `site_meta` (its docstring recommends supplying one because the default is probably
  wrong); `cli.py` fetches it from the Cloudnet API calibration endpoint per instrument PID and date. No S and no α.
- Microphysics products that feed our study, not a calibration: `products/der.py` (Frisch et al. 2002; defaults
  N = 200 cm⁻³, lognormal width σ_x = 0.35; `der_scaled` uses the LWP) gives our a_G; `products/lwc.py` is
  scaled-adiabatic LWC.
- The Cloudnet-era calibration tool that documents its inputs is the Reading page (§2.2).

### 2.6 Le et al. (2026), Atmos. Meas. Tech. 19, 4923–4941 (CL61)

(Preprint egusphere-2025-6331 in `doc/`; final paper, supplement and discussion read online.)

| item | where | statement |
|---|---|---|
| method | Sect. 3.2, p. 4929 | the theoretical MS contribution is computed for droplet diameters 8–20 µm with the CL61 beam divergence and receiver FOV; the integrated backscatter is scaled into the expected theoretical values; fits are monthly. |
| optics | Table 1, p. 4925 | divergence ±0.2 × ±0.35 mrad, receiver FOV ±0.56 mrad, 910.55 nm, pulse 160 ns FWHM, from the CL61 User Guide M212475EN. The FOV row was added at the request of referee 2 (AC2). |
| envelope | supplement Figs. S8–S11 | monthly C with the theoretical ηS for droplet effective diameters 8–20 µm (two curves, 0.5–5 km). |
| α, LWC, N, S value, model | — | **not stated**: no extinction in text, figures or replies; Hogan (2006) is not cited; S is described as known (citing O'Connor 2004 and Hopkin 2019) without a number; the edge or 1/e meaning of the ± values is not given. Cloud-calibration uncertainty ~10 % (citing Hopkin et al. 2019). |

### 2.7 Receiver FOV values: Wiegner et al. (2014), Atmos. Meas. Tech. 7, 1979–1997

- Table 1, p. 1982: receiver FOV as a **half angle**, from the Vaisala and Jenoptik manuals: CT25K 0.66,
  CL31 0.83, CL51 0.56 mrad. The CT25K value equals O'Connor's CT75K FOV of Sect. 2. No divergences, no edge
  convention.
- Sect. 4.2.1, p. 1985: the unknown multiple-scattering factors further limit the accuracy of the cloud method.
- The brief and a code comment (`_filters.py`, line 25) call this table "Wiegner & Gasteiger 2014". It is
  Wiegner et al. (2014) (Gasteiger is a co-author); Wiegner and Gasteiger (2015) is the water-vapour paper.

### 2.8 Other documents checked

- PROBE ALC chapter (2024, `doc/PROBE_jnd_chapter2a_ALC-V2024_05_31-1.pdf`), Appendix A: generic B = 1/(2ηS),
  S = 18.8 ± 0.8 sr, no α or droplet size. It adds that CL61 factory calibration is traced to a reference
  instrument calibrated with the cloud method.
- Kotthaus et al. (2016): CL31 signal artefacts; nothing on η.

### 2.9 Independent observations and climatologies of the extinction near cloud base

Added after the operator asked whether satellite data or a stratocumulus climatology (or optical depth plus
cloud depth) constrain α. What matters for η is the extinction in the layer that carries B: the lowest
~100–150 m above base, i.e. the first ~3 optical depths seen from below. The cloud mean and the cloud top are
different quantities.

| source | what is measured | layer | extinction | context |
|---|---|---|---|---|
| dual-FOV Raman lidar, Leipzig (Schmidt et al. 2015, Table 3, p. 10692) | α (532 nm), r_e, LWC and N from the lidar's own multiple scattering | 0–30 / 30–70 / 70–120 m above base | **11.5 ± 5.7 / 19.4 ± 7.0 / 25.5 ± 11.4 km⁻¹** (ranges 2.6–28, 3.9–36, 5.1–44) | 29 thin liquid altocumulus layers over polluted central Europe, 2010–2012: base 1.1–4.4 km (mean 2.9), thickness 95–300 m, τ 1.5–5.9 (Table 2, p. 10691); N 112 / 92 / 72 cm⁻³, r_e 5.8 / 9.0 / 10 µm, LWC 0.05 / 0.12 / 0.19 g m⁻³ in the three layers; α uncertainty 10–20 % (noise) plus 10–20 % (retrieval) (Table 1, p. 10690) |
| dual-FOV polarization lidar, Punta Arenas (Jimenez et al. 2020b, Sect. 4, p. 15272; Table 2, p. 15276) | the same quantities from the depolarization at two FOVs | 75 m above base | ≈ 20 km⁻¹ (altocumulus field, N 30–80 cm⁻³, r_e 7–10 µm, LWC ≈ 0.1 g m⁻³); 10.0 ± 2.3 and 9.3 ± 1.8 km⁻¹ (cloud at the top of the boundary layer, N 42 and 21 cm⁻³) | pristine Southern Ocean air |
| CALIPSO depolarization + MODIS r_e (Hu et al. 2007, Fig. 6, p. 3356) | monthly mean α of low water clouds in 5° boxes | top ~100 m of the cloud, oceans only | mostly 15–30 km⁻¹; 30–50 km⁻¹ downwind of polluted coasts (read from the colour scale); 10–20 km⁻¹ over the remote southern oceans | July and October 2006, January 2007; relation between depolarization and extinction from Monte Carlo simulations |
| ground-based τ and cloud depth, ARM SGP continental stratus (Dong et al. 1999, Table 1; Dong et al. 2000) | τ from surface transmission and LWP; boundaries from radar and lidar | whole cloud | τ/H ≈ 21–30 km⁻¹ (τ 21–26 and H 0.8–1.0 km by season; annual mean τ 25.1) | 25 months, Oklahoma; N 131–244 cm⁻³, r_e 8.1–9.7 µm, LWP ≈ 130 g m⁻² |
| droplet-number climatologies (in situ, MODIS) | N | whole cloud or cloud top | — | stratocumulus N from < 10 to > 500 cm⁻³; ≤ 50 cm⁻³ over remote oceans; > 200 cm⁻³ over and downwind of industrialised Northern Hemisphere continents and over the North Sea (Wood 2012, Sect. 5a and Fig. 23, pp. 2399–2400); continental mean 288 cm⁻³ in the compilation of Miles et al. (2000), as quoted by Dong et al. (1999) |
| lidar + radar + MWR synergy, Mace Head (Martucci and O'Dowd 2011) | Klett-inverted lidar extinction of the lowest 100–200 m, used to retrieve N | — | no α statistics given | two cases: continental N 382 cm⁻³, r_e 4.3 µm, depth 246 m; marine N 25 cm⁻³ |

The closest to what we need are the Leipzig statistics, which concern thin altocumulus rather than the thick
boundary-layer clouds used for calibration, and the CALIPSO statistics, which concern cloud tops over the oceans.

**τ/H is not the near-base extinction.** The whole-cloud mean weights the upper cloud, where LWC is largest. With
a linear LWC profile (α ∝ z^(2/3)), τ = 25 over 0.85 km reaches τ = 3 only ≈ 240 m above base (mean α ≈ 13 km⁻¹
below that), whereas an adiabatic cloud with the same N (≈ 235 cm⁻³) reaches τ = 3 at ≈ 100–110 m. Thick continental
stratus are therefore not adiabatic over their whole depth, or their radar–lidar boundaries overstate the liquid
layer, and τ/H cannot be turned into a near-base value without the profile.

**Adiabatic conversion** (our arithmetic; script table E). For an adiabatic cloud with constant N,
α(z) = (3/2)(4πk/3)^(1/3) ρ_w^(−2/3) N^(1/3) (f_ad Γ_l z)^(2/3), with Q_ext = 2, k = (r_v/r_e)³ = 0.8 and Γ_l the
adiabatic increase of LWC with height. At 900 hPa and f_ad = 1:

| cloud base (Γ_l) | N (cm⁻³) | α at 25 / 50 / 100 m above base (km⁻¹) | τ = 1 / τ = 3 reached at | mean α over τ ≤ 3 (km⁻¹) |
|---|---|---|---|---|
| 0 °C (1.61 g m⁻³ km⁻¹) | 50 | 9.7 / 15.4 / 24.5 | 79 / 154 m | 19.5 |
| 0 °C | 100 | 12.2 / 19.4 / 30.8 | 69 / 134 m | 22.4 |
| 0 °C | 200 | 15.4 / 24.5 / 38.8 | 60 / 116 m | 25.8 |
| 10 °C (2.09) | 100 | 14.6 / 23.1 / 36.7 | 62 / 120 m | 24.9 |
| 10 °C | 200 | 18.4 / 29.1 / 46.3 | 54 / 105 m | 28.6 |
| 10 °C | 300 | 21.0 / 33.4 / 53.0 | 50 / 97 m | 31.1 |
| 20 °C (2.45) | 200 | 20.4 / 32.4 / 51.4 | 51 / 98 m | 30.5 |
| 20 °C | 500 | 27.7 / 43.9 / 69.7 | 42 / 82 m | 36.6 |

The τ ≤ 3 mean depends only weakly on the inputs, as N^(1/5) (f_ad Γ_l)^(2/5) k^(1/5): over 50–500 cm⁻³ and
0–20 °C it spans 19.5–36.6 km⁻¹, and f_ad = 0.7 lowers it by 13 %.

What this says:

- The measured Leipzig profile follows the model. With N ≈ 100 cm⁻³ and a cold (≈ 0 °C) base, the model gives
  ≈ 9, 19.4 and ≈ 30 km⁻¹ at 15, 50 and 95 m; the lidar measured 11.5, 19.4 and 25.5 km⁻¹ in the three layers.
- For European boundary-layer stratocumulus (N ≈ 150–300 cm⁻³, base 0–15 °C) the τ ≤ 3 mean is ≈ 24–32 km⁻¹
  if adiabatic, ≈ 21–28 km⁻¹ with f_ad = 0.7. Our measured α_k (28–31 km⁻¹ by type) lies in this range.
  O'Connor's 15–20 km⁻¹ corresponds to clean, cold, sub-adiabatic clouds (e.g. N 50–100 cm⁻³ at 0 °C with
  f_ad 0.7: 17–19 km⁻¹) or to the local extinction 25–50 m above base.
- The regional and seasonal pattern of α_k found in the v2 study (subarctic maritime units −5.6 km⁻¹ below the
  type median, higher values in August) has the direction and roughly the size this scaling predicts for cleaner,
  colder clouds: N 100 cm⁻³ at 0 °C against 200 cm⁻³ at 10 °C gives 22.4 against 28.6 km⁻¹, and 10 → 20 °C at
  fixed N gives +7 %. This is consistent, not a proof.
- The extinction rises by a factor 2–3 between 25 and 100 m above base, across the very layer that carries B. A
  constant α is an equivalent value. Which constant reproduces the η of an adiabatic cloud, and whether the
  tail-slope estimator (fitted 20–150 m above the peak) returns that constant, is not addressed by any source
  (open point, §7).

## 3. Comparison table

CloudnetPy is left out: it has no η, α or S (§2.5).

| | O'Connor et al. 2004 | Reading calibration page | Hopkin et al. 2019 | Le et al. 2026 (CL61) | Hogan 2006 / 2008 | this project, v2 (2026-11) |
|---|---|---|---|---|---|---|
| η model, how obtained | Eloranta (1998); 4 homogeneous clouds → envelope; per-height median ηS fitted into it | Hogan's multiscatter code; ηS vs height | Hogan (2006) code, η per gate; code: one ladder for all Vaisala | not named; ηS envelope for 8–20 µm; monthly fit | PVC forward model (inputs α(r), a_G(r)) | PVC port of multiscatter 1.2.11 'original' (verified < 1e-4); deep homogeneous cloud to τ = 6; table at 13 CBH nodes per type; β·η(range) per gate |
| in-cloud α | **15 and 20 km⁻¹** | 15 and 20 km⁻¹ (example) | not stated | not stated | none for liquid (tests: ice 0.1–2; liquid 3 km⁻¹) | **30 km⁻¹**, one value for CL31/CL51/CL61 |
| nature of α | assumed, homogeneous; Miles et al. (2000) Sc range + full attenuation | as the paper | unknown | unknown | — | measured: PVC-inverted log-slope of the in-cloud decay at +20…+150 m above the peak (unit medians 27.8 / 30.2 / 31.0 km⁻¹); homogeneous in the model |
| droplet size | monodisperse effective diameter 8, 20 µm (a_G 4, 10 µm) | D 8, 20 µm | not stated | effective diameter 8–20 µm (r_e 4–10 µm) | a_G = equivalent-area radius | a_G 3.3 µm (Cloudnet `der_scaled`, lowest 100 m; r_e ≈ 4.0 µm, σ 0.35) |
| LWC, N implied (Q_ext = 2; our arithmetic) | N 24–200 cm⁻³, LWC 0.04–0.13 g m⁻³ | same | — | — | — | N ≈ 440 cm⁻³, LWC ≈ 0.08 g m⁻³ |
| receiver FOV | CT75K 0.66 mrad (Sect. 2; 0.75 in Sect. 3) | 0.66 mrad (default) | not stated | ±0.56 mrad | ρ_t | CL31 0.83, CL51 0.56, CL61 0.56 mrad |
| divergence | 0.75 mrad (Sect. 2; 0.66 in Sect. 3); abstract 1–1.5 mrad full angle | 0.75 mrad (default) | not stated | ±0.2 × ±0.35 mrad | ρ_l | CL31 0.57, CL51 0.21, CL61 0.28 mrad (RMS of the two axes) |
| convention | half-angles; Gaussian beam assumed; FOV edge unspecified | form fields are Hogan's ρ_l, ρ_t | — | manual "±", unspecified | ρ_l 1/e half-width (Gaussian beam), ρ_t sharp half-angle; Gaussian receiver option (2008) | ρ_t as top-hat edge, ρ_l as 1/e half-width; ×/÷ √2 in the uncertainty band |
| S | 18.8 ± 0.8 sr (905 nm; Mie, gamma, D₀ 8–20 µm, µ 2–10) | 18.8 ± 0.8 sr | 18.8 sr (905–1064 nm) | not given numerically | ice tests 20 sr | 18.8 sr (Mie at r_e 4.0 µm: 18.35 sr, not adopted) |
| CBH dependence | η falls with range | ηS vs height | per gate; code ladder 0.25–1.875 km | envelope 0.5–5 km | through the geometry | 0.25–3.75 km nodes |
| instrument dependence | recompute per lidar | from the inputs | per optics (text); none (code) | CL61-specific | — | per type |
| η at CBH ≈ 1 km | 0.83 ± 0.09 (CT75K) | not retrieved | 0.7–0.85 (all heights, 905–1064 nm) | ≈ 0.83–0.93 (digitised, S 18.8 sr) | — | CL31 0.773, CL51 0.816, CL61 0.820 |

## 4. Where our choices differ

### 4.1 α

- The only α stated in this literature is O'Connor's **assumed, homogeneous 15–20 km⁻¹**. Hopkin et al. and
  Le et al. state none; Hogan prescribes none for liquid water; CloudnetPy has none. No source measured α on its
  calibration scenes, and **no source used an adiabatic (height-increasing) extinction profile for η**.
- Ours, 30 km⁻¹, is 1.5–2× O'Connor's value and 3× the generic 10 km⁻¹ of the 2026-07 tables. It is the only
  measured value: the PVC-equivalent homogeneous α that reproduces the observed tail decay over the first
  ~150 m. It is still a constant in the forward model.
- Independent observations (§2.9) do not contradict it: dual-FOV lidar over Leipzig measures 19–26 km⁻¹ at
  30–120 m above the base of thin altocumulus, and the adiabatic model with typical continental droplet numbers
  gives ≈ 24–32 km⁻¹ over the layer a ground lidar sees, against 17–22 km⁻¹ for clean or cold clouds.
- Orientation for the constant-versus-adiabatic question (our arithmetic, not from the sources): a
  monodisperse adiabatic cloud with N = 200 cm⁻³ and LWC rising by 2 g m⁻³ km⁻¹ has α ≈ 17, 30 and 48 km⁻¹ at 20,
  50 and 100 m above base (α ∝ N^(1/3) z^(2/3)); it reaches optical depth 1 at ≈ 53 m and 3 at ≈ 100 m, so its mean α
  over the layer that carries B is ≈ 30 km⁻¹. A homogeneous 15–20 km⁻¹ cloud reaches optical depth 3 only at
  150–200 m. How the profile shape (low α at the base, high α deeper) changes η relative to a constant of the
  same layer mean is not addressed by any of the sources.
- Size of the difference in our model (a_G 3.3 µm, our optics, S fixed; `eta_literature_crosscheck.py`, table D).
  C_L change if α were 20 / 15 / 10 km⁻¹ instead of 30:

| type | CBH 0.5 km | CBH 1 km | CBH 2 km |
|---|---|---|---|
| CL31 | +4.2 / +6.7 / +9.6 % | +5.9 / +9.7 / +14.3 % | +7.2 / +12.2 / +18.8 % |
| CL51 | +3.4 / +5.3 / +7.5 % | +5.1 / +8.2 / +11.9 % | +6.8 / +11.3 / +17.0 % |
| CL61 | +3.3 / +5.2 / +7.3 % | +5.0 / +8.0 / +11.6 % | +6.7 / +11.1 / +16.6 % |

### 4.2 Droplet size

- Ours, a_G = 3.3 µm (equivalent-area radius, r_e ≈ 4.0 µm), is **below O'Connor's smallest case** (monodisperse
  D = 8 µm, a_G = 4 µm) and at the small end of Le's 8–20 µm effective-diameter band. It equals Le's lower bound only
  if that bound is r_e = 4 µm in a lognormal of width 0.35 (a_G = 0.832 r_e); entered as a monodisperse size it is
  4 µm. Le et al. do not say which.
- For S, O'Connor's 18.8 ± 0.8 sr holds for D₀ 8–20 µm. Our scenes have D₀ ≈ 8.4 µm (r_e e^(σ²/2), σ 0.35): at
  the lower edge, where O'Connor note that S starts to rise for broad distributions (p. 783). The operator kept
  18.8 sr (decision of 2026-11-02); our Mie gives 18.35 sr at r_e 4.0 µm.

### 4.3 FOV and divergence

- The numbers agree: our CL31 / CL51 FOV are those of Wiegner et al. (2014) Table 1, the CL61 values are those of
  Le et al. (2026) Table 1 (both from the Vaisala manuals).
- The convention is unspecified in every source that gives instrument values. O'Connor call theirs half-angles and
  assume a Gaussian beam; Wiegner calls the FOV a half angle; Le et al. reproduce the manual's "±". Only the model
  defines its inputs (ρ_l 1/e Gaussian half-width, ρ_t sharp half-angle), and our mapping follows it.
- O'Connor's swap of FOV and divergence between Sect. 2 and Sect. 3 is immaterial: in our port the two readings
  differ by ≤ 0.001 in η (table A).

### 4.4 S

All calibration sources use 18.8 sr (Le et al. by citation, not by number); so do we.

### 4.5 How η is used

O'Connor and Le et al. fit the observations into an **envelope** spanned by the droplet-size (and α) extremes.
Hopkin et al. and this project apply **one η curve** per gate. With an envelope, α and droplet size act as an
uncertainty; with one curve they set the level.

## 5. Cross-check: the literature η through our PVC port (inference)

Script `validation/eta_study/eta_literature_crosscheck.py` (a few seconds; tables A–D and the figure). The PVC
port is `validation/multiple_scattering_eta.py`, verified against Hogan's multiscatter 1.2.11 to < 1e-4
(`04_multiple_scattering.md` §3b).

![Literature eta vs our PVC port](figs_eta_literature/fig_eta_literature.png)

*(a) O'Connor's cloud set with the CT75K optics: the published η (black) sit above the PVC band computed with the
documented half-angles (blue) and inside the band computed with the FOV and divergence halved (orange). (b) The
same for the CL61 with Le et al.'s droplet sizes and O'Connor's α: the digitised Le et al. curves (black) follow
the halved-optics band. (c) The CL61 tables of `_filters.py` against the Le et al. envelope: the operational v2
table runs just below the lower edge and parallel to it, the 2026-07 table inside it, the legacy ladder across
it. PVC: homogeneous cloud integrated to optical depth 6.*

### 5.1 O'Connor et al. (2004): published η vs PVC with the same inputs

η range over D 8 / 20 µm × α 15 / 20 km⁻¹ (mid in brackets):

| CBH | published | PVC, documented optics (FOV 0.66, div. 0.75 mrad) | optics × 0.71 | optics × 0.5 |
|---|---|---|---|---|
| 1.0 km | 0.83 ± 0.09 | 0.686–0.829 (0.757) | 0.729–0.864 (0.796) | 0.772–0.894 (0.833) |
| 1.5 km | 0.79 | 0.640–0.785 (0.713) | 0.681–0.824 (0.753) | 0.724–0.860 (0.792) |
| 2.7 km | 0.74 ± 0.04 | 0.584–0.715 (0.650) | 0.617–0.758 (0.687) | 0.655–0.799 (0.727) |
| 4.0 km | 0.73 ± 0.06 | 0.555–0.669 (0.612) | 0.581–0.711 (0.646) | 0.613–0.753 (0.683) |

With the documented optics, PVC gives 0.07–0.12 less than the published centres; the published centres match the
weakest-scattering corner (D 8 µm, α 15 km⁻¹: 0.829 and 0.785 at 1 and 1.5 km). With the optics halved, the band
centres match the published centres within 0.015 up to 2.7 km (0.05 at 4 km).

### 5.2 Le et al. (2026): which α reproduces the ηS envelope?

Digitised from supplement Fig. S8, panel 2023-12 (reading error ~0.05 sr; η = ηS / 18.8 sr):

| CBH | 0.6 km | 1.0 km | 1.5 km | 2.0 km | 3.0 km | 4.0 km |
|---|---|---|---|---|---|---|
| ηS, D 20 µm curve (sr) | 16.55 | 15.57 | 14.70 | 14.04 | 13.09 | 12.45 |
| ηS, D 8 µm curve (sr) | 18.00 | 17.51 | 16.98 | 16.53 | 15.79 | 15.19 |
| η, D 20 µm | 0.880 | 0.828 | 0.782 | 0.747 | 0.696 | 0.662 |
| η, D 8 µm | 0.957 | 0.931 | 0.903 | 0.879 | 0.840 | 0.808 |

Joint rms difference in η (both curves, 0.6–4 km) between PVC with the CL61 optics and the digitised curves:

| optics | α 5 | α 10 | α 15 | α 20 | α 30 km⁻¹ |
|---|---|---|---|---|---|
| documented (FOV 0.56, div. 0.28 mrad) | 0.064 | **0.038** | 0.097 | 0.141 | 0.204 |
| × 0.71 (0.40 / 0.20 mrad) | 0.106 | **0.027** | 0.044 | 0.086 | 0.149 |
| × 0.5 (0.28 / 0.14 mrad) | 0.142 | 0.069 | **0.023** | 0.036 | 0.094 |

With the documented optics, the D 20 µm curve is matched at α ≈ 10 km⁻¹ (rms 0.014) and the D 8 µm curve at
α ≈ 5 km⁻¹ (rms 0.022). With the optics halved, O'Connor's 15–20 km⁻¹ fits both. α and the optics convention are
degenerate in these curves, so **they do not identify the α used by Le et al.**

### 5.3 Operational CL61 tables against the Le et al. envelope

η of the table minus the lower (D 20 µm) edge of the envelope:

| table | 0.6 km | 1.0 km | 1.5 km | 2.0 km | 3.0 km | 4.0 km |
|---|---|---|---|---|---|---|
| v2, operational (a_G 3.3 µm, α 30 km⁻¹) | −0.010 | −0.008 | −0.008 | −0.008 | −0.008 | −0.010 |
| 2026-07 (a_G 5.5 µm, α 10 km⁻¹) | +0.032 | +0.047 | +0.056 | +0.060 | +0.059 | +0.043 |
| legacy ladder (CL51 rows, extrapolated above 2.375 km) | −0.062 | −0.022 | +0.008 | +0.028 | +0.050 | +0.056 |

The v2 table has the shape of the envelope's lower edge, about 1 % lower.

### 5.4 What this does and does not show

- Both published η sets (O'Connor's Eloranta values; Le et al.'s curves) are reproduced by our port **with their
  stated droplet sizes and O'Connor's α only if the FOV and divergence are about half our documented
  half-angles**. With our optics they need less extinction (Le: 5–10 km⁻¹) or remain 0.07–0.12 higher (O'Connor).
- The cause cannot be identified from the papers. Candidates: O'Connor's Eloranta sum truncated at a low order
  (the explicit sum needs fifth or sixth order to converge, Hogan 2006 p. 5989; truncation underestimates MS); the
  lobe-width convention (p. 5986); an optics convention (full vs half angle, edge vs 1/e, Gaussian receiver); the
  integration depth of B; or, for Le et al., a lower α.
- Consequences for this project: (i) the literature η values cannot serve as an independent test of α = 30 km⁻¹;
  (ii) the optics convention, already the largest term of our budget, may involve a factor up to 2 between
  implementations, wider than the ×/÷ √2 now carried; (iii) per type, the v2 CL61 table sits ~1 % below the Le et al.
  envelope, i.e. it applies slightly more MS than the strongest case of the published CL61 calibration.

## 6. Statements in our documents not supported by the literature (not edited here)

| document, place | statement | what the sources say |
|---|---|---|
| `04_multiple_scattering.md` §3a, bullets after the fit | the legacy tables bracket O'Connor's constant η = 0.7 | O'Connor give no constant 0.7: 0.83 ± 0.09 (1 km), 0.79 (1.5 km), 0.74 (2.7 km), 0.73 ± 0.06 (4 km); 0.5–1 is only the possible range. |
| `04_multiple_scattering.md` §10, Literature paragraph | all published methods feed a fixed droplet size to Hogan's code; O'Connor 2004 used a single representative r_e | O'Connor used Eloranta's model with two sizes × two extinctions (an envelope); only the Reading page uses Hogan's code. |
| `04_multiple_scattering.md` §3c | Le & O'Connor (2026) still use the same Hogan code | the paper names no MS model; plausible (the Reading page uses Hogan's code) but not documented. |
| `eta_switch_decision_2026-11-02.md` §7 | S = 18.8 sr attributed to O'Connor, Hopkin and Cloudnet | CloudnetPy contains no S; the Cloudnet-side source is the Reading page (18.8 ± 0.8 sr). Le et al. do not give the number. |
| `calibration/cloud/_filters.py` line 25 (comment) | type optics as in Wiegner & Gasteiger 2014 | the FOV table is Wiegner et al. (2014), Table 1. |

## 7. Open points

| point | why it matters | what would settle it |
|---|---|---|
| optics convention (edge / 1/e, half / full angle) | largest term of the η budget; §5.4 suggests up to a factor 2 between implementations | written answer from Vaisala (already pending); the inputs and code options behind Le et al. Figs. S8–S11 (FMI) |
| what the Reading tool returns | it is the documented source of the legacy ladder and probably of the Le et al. envelope | tried 2026-11-02 with the operator's go: the calculator returns HTTP 500 for every request and no output is archived. Remaining route: ask the authors (FMI) for the tool's cloud cases, model options and optics convention, together with the inputs of Le et al. Figs. S8–S11 |
| constant vs adiabatic α in η | no source addresses it; the extinction rises ×2–3 across the layer that carries B (§2.9) | a PVC run with an adiabatic α(z), a_G(z) profile (N, Γ_l as in §2.9) against constants, and the tail-slope estimator applied to the same simulated profile, to see which constant it returns |

## 8. References

- Dong, X., Minnis, P., Ackerman, T. P., Clothiaux, E. E., Mace, G. G., Long, C. N., and Liljegren, J. C.: A
  25-month database of stratus cloud properties generated from ground-based measurements at the Atmospheric
  Radiation Measurement Southern Great Plains site, J. Geophys. Res., 105, 4529–4537, 2000; the seasonal numbers
  used here are from Table 1 of the corresponding extended abstract, Proc. Ninth ARM Science Team Meeting, San
  Antonio, 1999, https://www.arm.gov/publications/proceedings/conf09/extended_abs/dong_x.pdf (cited as Dong et al.
  1999).
- Eloranta, E. W.: Practical model for the calculation of multiply scattered lidar returns, Appl. Opt., 37,
  2464–2472, 1998 (cited, not read).
- Frisch, S., Shupe, M., Djalalova, I., Feingold, G., and Poellot, M.: The retrieval of stratus cloud droplet
  effective radius with cloud radars, J. Atmos. Oceanic Technol., 19, 835–842,
  https://doi.org/10.1175/1520-0426(2002)019<0835:TROSCD>2.0.CO;2, 2002 (via CloudnetPy `der.py`).
- Hogan, R. J.: Fast, approximate calculation of multiply scattered lidar returns, Appl. Opt., 45, 5984–5992,
  https://doi.org/10.1364/AO.45.005984, 2006; author copy
  https://www.met.reading.ac.uk/~swrhgnrj/publications/multiscatter.pdf.
- Hogan, R. J.: Fast lidar and radar multiple-scattering models. Part I: Small-angle scattering using the photon
  variance–covariance method, J. Atmos. Sci., 65, 3621–3635, https://doi.org/10.1175/2008JAS2642.1, 2008; author copy
  https://www.met.reading.ac.uk/~swrhgnrj/publications/multiscat_small_angle.pdf.
- Hogan, R. J.: multiscatter code, version 1.2.11 (2015), https://www.met.reading.ac.uk/clouds/multiscatter/.
- Hopkin, E., Illingworth, A. J., Charlton-Perez, C., Westbrook, C. D., and Ballard, S.: A robust automated technique
  for operational calibration of ceilometers using the integrated backscatter from totally attenuating liquid
  clouds, Atmos. Meas. Tech., 12, 4131–4147, https://doi.org/10.5194/amt-12-4131-2019, 2019
  (`doc/Hopkins-amt-12-4131-2019.pdf`); discussion amt-2018-427 RC1, RC2, AC1, AC2,
  https://amt.copernicus.org/articles/12/4131/2019/amt-12-4131-2019-discussion.html.
- Hu, Y., Vaughan, M., McClain, C., Behrenfeld, M., Maring, H., Anderson, D., Sun-Mack, S., Flittner, D., Huang, J.,
  Wielicki, B., Minnis, P., Weimer, C., Trepte, C., and Kuehn, R.: Global statistics of liquid water content and
  effective number concentration of water clouds over ocean derived from combined CALIPSO and MODIS measurements,
  Atmos. Chem. Phys., 7, 3353–3359, https://doi.org/10.5194/acp-7-3353-2007, 2007.
- Jimenez, C., Ansmann, A., Engelmann, R., Donovan, D., Malinka, A., Schmidt, J., Seifert, P., and Wandinger, U.: The
  dual-field-of-view polarization lidar technique: a new concept in monitoring aerosol effects in liquid-water
  clouds – theoretical framework, Atmos. Chem. Phys., 20, 15247–15263, https://doi.org/10.5194/acp-20-15247-2020,
  2020a; – case studies (with Wiesen, R., Radenz, M., Yin, Z., Bühl, J., and Barja, B.), Atmos. Chem. Phys., 20,
  15265–15284, https://doi.org/10.5194/acp-20-15265-2020, 2020b.
- Kotthaus, S., O'Connor, E., Münkel, C., Charlton-Perez, C., Haeffelin, M., Gabey, A. M., and Grimmond, C. S. B.:
  Recommendations for processing atmospheric attenuated backscatter profiles from Vaisala CL31 ceilometers,
  Atmos. Meas. Tech., 9, 3769–3791, https://doi.org/10.5194/amt-9-3769-2016, 2016 (`doc/Kotthaus-amt-9-3769-2016.pdf`).
- Le, V., O'Connor, E. J., Filioglou, M., and Vakkari, V.: Operational performance of the Vaisala CL61 ceilometer
  for atmospheric profiling, Atmos. Meas. Tech., 19, 4923–4941, https://doi.org/10.5194/amt-19-4923-2026, 2026;
  supplement https://amt.copernicus.org/articles/19/4923/2026/amt-19-4923-2026-supplement.pdf; preprint
  egusphere-2025-6331 ("Long-term performance of the Vaisala CL61 ceilometer for atmospheric profiling",
  `doc/LE-OConnor_egusphere-2025-6331.pdf` and `-supplement.pdf`) with discussion RC1, RC2, AC1, AC2,
  https://egusphere.copernicus.org/preprints/2026/egusphere-2025-6331/.
- Martucci, G. and O'Dowd, C. D.: Ground-based retrieval of continental and marine warm cloud microphysics, Atmos.
  Meas. Tech., 4, 2749–2765, https://doi.org/10.5194/amt-4-2749-2011, 2011.
- Miles, N. L., Verlinde, J., and Clothiaux, E. E.: Cloud droplet size distributions in low-level stratiform clouds,
  J. Atmos. Sci., 57, 295–311, 2000 (cited by O'Connor et al. and Dong et al.; not read).
- O'Connor, E. J., Illingworth, A. J., and Hogan, R. J.: A technique for autocalibration of cloud lidar,
  J. Atmos. Oceanic Technol., 21, 777–786, https://doi.org/10.1175/1520-0426(2004)021<0777:ATFAOC>2.0.CO;2, 2004
  (`doc/OConnor-2004-atot-1520-0426_2004_021_0777_atfaoc_2_0_co_2.pdf`).
- PROBE: ALC chapter 2a, version 2024-05-31 (`doc/PROBE_jnd_chapter2a_ALC-V2024_05_31-1.pdf`).
- Schmidt, J., Ansmann, A., Bühl, J., and Wandinger, U.: Strong aerosol–cloud interaction in altocumulus during
  updraft periods: lidar observations over central Europe, Atmos. Chem. Phys., 15, 10687–10700,
  https://doi.org/10.5194/acp-15-10687-2015, 2015.
- University of Reading lidar-calibration page, http://www.met.reading.ac.uk/~swr99ejo/lidar_calibration/index.html
  (HTML and `example.png` read 2026-11-02; its calculator script returned HTTP 500 that day).
- Wiegner, M., Madonna, F., Binietoglou, I., Forkel, R., Gasteiger, J., Geiß, A., Pappalardo, G., Schäfer, K., and
  Thomas, W.: What is the benefit of ceilometers for aerosol remote sensing? An answer from EARLINET, Atmos. Meas.
  Tech., 7, 1979–1997, https://doi.org/10.5194/amt-7-1979-2014, 2014 (`doc/Wiegner-amt-7-1979-2014.pdf`).
- Wood, R.: Stratocumulus clouds, Mon. Weather Rev., 140, 2373–2423, https://doi.org/10.1175/MWR-D-11-00121.1, 2012.
- Wiegner, M. and Gasteiger, J.: Correction of water vapor absorption for aerosol remote sensing with ceilometers,
  Atmos. Meas. Tech., 8, 3971–3984, https://doi.org/10.5194/amt-8-3971-2015, 2015 (named here only to separate it
  from Wiegner et al. 2014).
- Code: CloudnetPy, https://github.com/actris-cloudnet/cloudnetpy (`instruments/ceilo.py`, `cli.py`,
  `products/der.py`, `products/lwc.py`; `main` at a359dc4, 2026-10-22). E-PROFILE Python of Hopkin et al. (2019):
  `VaisCeil_DayCalibration_EProf_EH.py` (MeteoSwiss `alc_calib_py`). This repository:
  `validation/multiple_scattering_eta.py`, `calibration/cloud/_filters.py`,
  `validation/eta_study/eta_literature_crosscheck.py`.
