# Dashboard d'inter-comparaison — v3 (page unique, multi-sites)

**Entrée : [`index.html`](index.html)** — autonome (double-clic), ou servi via l'entrée
`intercomp` de `.claude/launch.json`. Une seule page porte **Payerne, Amsterdam et Lindenberg**
(sélecteur en haut), ~10 Mo tout compris.

## Architecture v3 (pourquoi la page est libre ET légère)

La v2 précalculait chaque combinaison (variante × WV × λ × période) côté serveur : 18 combos et
39 Mo pour Payerne seul, et aucune liberté par instrument. La v3 expédie des **ingrédients par
canal** et recombine dans le navigateur :

```
sortie(h, z) = A(jour, z) · [brut(h, z) − dark(z)] / C_L(t_h) + B(jour, z)
```

- `brut` : bloc horaire × porte, quantifié log-uint16 (~0,007 % de pas), non calibré ;
- `A/B` : facteurs de correction par (jour, porte) — division WV, conversion λ affine ;
- `dark(z)` : profil capot mesuré (Payerne), soustraction **exacte par porte** ;
- `C_L(t)` : une série par variante d'étalonnage, interpolée sur l'axe horaire.

Filtrage de période, choix de référence, médianes, IQR, histogrammes, table de stats : de
l'arithmétique par heure côté client — **exacte**. `check_v3.py` rejoue **tous les combos v2**
(jusqu'à 16 par site : 4 variantes × WV on/off × molecular/none/angstrom, dark compris) contre le
câblage v3 et **sort en erreur** au-delà de 0,5 pt d'écart absolu (pire écart sain : 0,28 pt ;
les lignes L1 des variantes dark sont informatives — l'instantané v2 précède le refit des darks
du 2026-08-16 — leurs lignes L2, sans constante, restent verrouillées).

## Les trois fichiers

| fichier | rôle |
|---|---|
| `variants_v3.py` | **la matrice de contrôle** — méthodes par type, échelles de variantes, chemins des runs, étiquettes (par type), jumeaux dark, motifs de grisage. *C'est le fichier à éditer quand un nouveau run arrive.* |
| `build_v3.py` | calcule les payloads par site (ingrédients + Hopkin binnés + rideaux + PWV) |
| `render_v3.py` | assemble la page (contrôles, panneaux, JS de recombinaison) |

Reconstruction : `python build_v3.py` puis `python render_v3.py` (payloads sur
`C:/DATA/Projects/202606_E-PROFILE_calibration/`).

## Doctrine des contrôles

- **Omis vs grisé** : un menu n'offre que ce qui a un sens physique pour l'instrument (pas
  d'entrée WV pour un CHM15k à 1064 nm) ; le grisé + infobulle est réservé aux variantes
  sensées dont le run n'existe pas encore — le motif est toujours affiché.
- **Échelle λ CL61** (la question du spectre d'émission) : 910,74/1,0 (modèle historique) →
  910,55/1,0 → 910,55/0,1 → **910,55 σ0,08 ★ (spectre constructeur**, Le & O'Connor réponse à
  RC1, comm. pers. Vaisala ; corroboré DA10 M212895EN-E + Mariani 2021**)** → sans WV. Produite
  par `ALC_WV_SPECTRUM` (JSON, ex. `{"CL61": [910.55, 0.188]}`) qui pilote les DEUX méthodes via
  l'unique table `water_vapor.LASER_SPECTRUM` ; `ALC_WV_DISABLE=1` coupe la correction.
- **Case « soustraire le dark mesuré (capot) »** : profils corrigés par porte + bascule
  automatique vers la série de constantes jumelle (`DARK_TWIN`) — jamais de profil corrigé
  divisé par une constante non corrigée. En mode nuage, la constante est immune au dark
  (+0,08 % CL31, ~0,001 % CL61, mesuré) : l'info ⓘ le dit, ce demi-état est l'état correct.
- **Nombres dans la prose** : tout nombre dépendant de l'état des données est injecté au build
  depuis le payload — jamais écrit en dur (règle anti-fossile, après audit complet 2026-08-16).
- **Cohérence partout, pas seulement au clic** (revue adverse 2026-08-16) : le panneau PWV
  soustrait le dark de sa référence quand ses constantes viennent d'un run dark (l'hybride
  gonflait chaque pente de +0,4–0,8 %/mm) ; le tableau récapitulatif balaye chaque variante
  Rayleigh dans SON appariement dark cohérent ; l'aller-retour nuage→Rayleigh→nuage restaure le
  dark du canal nuage ; hors Payerne, les runs « dark » sont étiquetés **dark ESTIMÉ ciel clair**
  (une seule id de variante, deux sources très différentes).

## Filtre bruit (SNR) — 2026-08-24

Sélecteur « **Filtre bruit** » : retire ou masque les données non détectées à **SNR ≥ 3**
(SNR = médiane de fenêtre / (σ robuste/√n), le critère « SNR3 » opérationnel, évalué sur le
signal **dark-soustrait** là où le b(z) mesuré existe — la case « dark mesuré » ne change PAS
les masques). Deux menus + un seuil :

- **Fenêtre d'évaluation** : 5 min (défaut) / 30 min / 60 min / 3 h, alignée horloge. Le payload
  étant horaire, les décisions sous-horaires sont **précalculées au build** sur les flux L1
  natifs (`nf_v3.py`, lu via `l1_l2_io.read_l1_native`) : l'heure est admise si ≥ 50 % de ses
  données (`NF_FMIN`) tombent dans des fenêtres détectées — la valeur horaire reste la médiane
  inconditionnelle (l'alternative « médiane des survivants » coûterait un bloc de valeurs par
  mode × fenêtre ; l'approximation est chiffrée par `check_v3.py`, non bloquant).
- **Mode** — car un filtrage SNR conditionne l'échantillon au signal : le CL31 (bruyant) ne
  garderait que ses scènes chargées et sa médiane serait biaisée haute vs CL61/CHM15k :
  - `par instrument` : le filtre littéral, chaque canal sous son propre SNR. Le **biais
    d'échantillonnage est affiché** sous l'histogramme : médiane du référent sous le masque du
    CL31 vs sans restriction — de l'échantillonnage pur, mesuré sur l'instrument silencieux.
  - `masque commun (∩)` : un pixel n'entre que si TOUS les instruments le détectent —
    échantillon identique, couverture dictée par le plus bruité.
  - `moyenne d'abord, SNR après` : rien n'est retiré (un bruit à moyenne nulle ne biaise pas une
    médiane) ; le profil de période est masqué là où le SNR de l'**agrégat** < 3 (sommes
    mensuelles `sn`/`ss` par porte). Statistiquement le plus propre.
  - `masque scène (réf. CHM15k)` : masque construit depuis le **référent** (`nf_ref`, le CHM15k
    du site — décision opérateur 2026-08-24) et appliqué à tous ; la sélection dépend de
    l'atmosphère, pas du bruit de chacun. Trois seuils **β_att absolus** (décision opérateur) :
    **0,1 / 0,25 / 0,5 Mm⁻¹ sr⁻¹** à 1064 nm — la fenêtre passe quand le référent mesure
    β_att ≥ max(seuil, sa propre limite 3σ/√n), « β_att détecté par la référence ».
    Défaut : 0,25.
- Les masques sont des **bits par (heure appariée × porte)** : u1 par instrument (bit =
  fenêtre), u1 intersection, u2 scène (bit = seuil×4+fenêtre) ; constantes **gelées aux
  variantes par défaut** du site — changer de variante ne change jamais l'échantillon. Panneaux
  L2 : masques hérités du L1 (mêmes photons). La courbe pointillée (axe du haut des profils) =
  % d'heures retenues par porte. Rideaux et panneau PWV : statiques, non filtrés.
- Vérification (`check_v3.py`, bloquant à 0,1 %) : admission 60 min reconstruite via
  `retime_hourly`/`snr_mask` (l'implémentation opérationnelle), intersection et scène
  reconstruites, sommes p2 directes, force brute 5 min sur 2 jours.

## Correction de recouvrement CHM15k (Hervo 2016) — 2026-08-25

Case « **recouvrement T (Hervo 2016)** », colonne *Recouvrement*, offerte **uniquement pour les
CHM15k** (le modèle et la physique sont propres à cet instrument ; les Vaisala affichent « — »).

- Le L1 `rcs_0` porte déjà le recouvrement de **référence statique** de l'instrument ; ce qui
  reste est un résidu qui dépend de la **température interne** : `Dif(z) = a(z)·T_int + b(z)` [%],
  d'où `rcs_corrigé = rcs₀ · (1 + Dif/100)` — le recouvrement de référence s'annule.
- Noyau d'application : [`validation/paper/overlap.py`](../validation/paper/overlap.py) (vérifié
  contre la chaîne MATLAB / `overlap_probe_eprofile`). Modèles : un NetCDF par instrument dans
  `ALC_OVERLAP_DIR` (défaut `D:/TEMP_MODELS/202606`), indexés par (`wigos_station_id`,
  `instrument_id`) — les 6 CHM15k des trois sites en ont un.
- Le payload n'embarque que `a(z)`, `b(z)` sur la grille d'affichage (6 chiffres significatifs) et
  la **température interne par heure** (`temp_int`, portée par le même retimage que le signal) : le
  navigateur recalcule le facteur, borné comme le noyau (`fmin`/`fmax`).
- Ordre : **après** le dark (offset électronique additif), **avant** la division par C_L. Le
  facteur étant **sans dimension**, il s'applique aussi au **L2** (même instrument, même optique) —
  contrairement au dark, qui est en unités `rcs₀`.
- **Portée honnête, chiffrée au build et affichée en infobulle** : `a` et `b` s'annulent une fois
  le recouvrement complet (~690 m à Payerne, 750–840 m ailleurs), donc la correction **remodèle le
  proche champ uniquement** et est un **no-op exact au-dessus** ; sur la bande de statistiques
  500–3000 m elle déplace la médiane de **+0,01 %** (Payerne). C'est aussi pourquoi l'appairage
  avec des constantes non corrigées reste cohérent : les runs d'étalonnage ne l'ont jamais vue et
  elle ne touche pas leur fenêtre d'ajustement. Les masques du filtre bruit, gelés au build, ne
  changent pas non plus (même doctrine que le dark).
- Vérification (`check_v3.py`, bloquant à 1e-5) : le facteur reconstruit depuis le payload est
  comparé à `overlap.correct_rcs` sur la grille d'affichage — écart 2,1e-06 (précision de stockage
  de `a`/`b`), et le facteur est exactement 1 au-dessus de `ztop`.

## Provenance des données

Runs par variante sous `C:/DATA/Projects/202606_E-PROFILE_calibration/` (`calout_v22_04` réseau,
`diag_v22_*` locaux — voir `variants_v3.py`), extraits Hopkin par configuration dans
`cloud_profile_dump_<tag>/`. L1 : `D:/E-PROFILE_L1_2026` ; CAMS **0,4°** uniquement
(`A:/CAMS_Monthly_04;D:/CAMS_daily` — jamais `D:/CAMS` ni `D:/CAMS_run_v20`, monolithes 1°).
L2 as distributed: the CHM15k **and the CL61** carry the operational Rayleigh v1.0 calibration;
the CL31 stays uncalibrated (default constant). The CL61's applied constant is the default 1.0
until 2026-07-28, then 2.1257 from 07-29 to 10-13 and 0.7932 from 10-14; the CHM15k's changes on
the same two dates (5.51 → 5.88 → 5.83 → 5.30 → 4.87e11). The page dates these runs itself
from the L2 files (L2-row tooltip and pill); no date is written in its prose (corrected
2026-11-02: earlier text said the CL61 switch was mid-June).

**Hood dark switched from v1 to v3 (2026-11-02).** Payerne's measured-dark state now uses the
hood profile **v3** on both sides of the pairing: profile side `dark_profiles_payerne_v3.npz`
(`sites.py` `dark_npz`; 10 sessions May–Sep 2026 for the CHM15k and the CL61, the CL31's 6
post-optical-block-swap sessions; covered windows detected in the data, pooled median smoothed
in P = rcs₀/r²), constants side the reruns `diag_v22_dark_v3` (`variants_v3.DARK_MEAS`, also
`L55S008DARK`; eprof_v2.2, CL61 on the operational 910.55/0.188 spectrum, Rayleigh only).
Before, both sides used the 6-session logbook-window profile **v1** (`dark_profiles_payerne.npz`,
runs `diag_v22_dark` / `diag_v22_l55s008dark`). Reason: v1 produced a spurious CL61-Rayleigh-vs-
CHM15k gap; with v3 the two Rayleigh calibrations agree every month
([`doc/reports/eta_payerne_gap_2026-11-02.md`](../doc/reports/eta_payerne_gap_2026-11-02.md)).
Effect on the summary table (default state, CL61 switched to Rayleigh, 2026-06-11..08-14, CL61
★ + dark vs CHM15k v2.2 + dark, 500–3000 m median relative difference): **+7.8 % → −0.7 %**;
per month all hours Jun/Jul/Aug +18.6/+6.2/−2.8 % → +2.7/−0.1/−4.3 %, night only (sun < −10°)
+21.1/+7.2/−0.9 % → +1.9/+0.7/−2.3 %. The v3 rerun carries only the operational CL61 spectrum,
so the CL61 `v2.2dark` entry (910.74 + dark) is now greyed out at Payerne rather than shown on
another dark version.

**λ910,74 WV-comparison mode pinned to the retired CL61 line (2026-11-02).** `build_v3` took the
"nominal" water-vapour array from `validation/paper/intercompare.WV_PARAMS`, whose CL61 entry moved
to 910.55/0.188 on 2026-10-20 (commit a25977e). Every rebuild since then would have turned the
page's « λ910,74 » mode (and the PWV panel's 910,74 states) into a copy of the ★ mode; `check_v3`
caught it (all WV-on CL61 rows off by up to 38 points). `variants_v3.NOMINAL_WV_SPECTRUM` now pins
that mode to 910.74/1.0 for the CL61; the other types keep `WV_PARAMS`.

**Window extended to 2026-10-31 (2026-11-02).** The Payerne L2 for Aug–Oct and the October CAMS
daily files were synced from `balfrin:/scratch/mch/mhrvo/{E-PROFILE_L2,CAMS}` into
`D:/E-PROFILE_L2_2026` and `D:/CAMS_daily`; the stream cache now spans 2026-06-01..10-31 (the
2026-06-01..08-14 one is kept as `_streams_payerne_to20260814.npz`). Constants: every run-based
variant is carried past the end of its original run by a rerun with its own settings
(`run_payerne_ext.ps1` → `dashboard_ext_20261031/<config>/`, mapped in `variants_v3.EXT_TAG`;
2026-08-14..10-31, the hood-dark v3 pair 2026-10-09..10-31). The builder reads [original series,
extension] and appends only the extension nights dated after the original series' last night, so
the original series is never altered. « v2.0 » (published NetCDFs) cannot be rerun: its Kalman is
held after its last night, as the page prose states. The reruns reproduce the eta-study extension
runs of the same settings bit for bit on every common night; the only new nights are the October
ones that had no CAMS before.
`check_v3` keeps its v2 regression on the snapshot window (2026-06-11..08-14) with snapshot-era
constants: each series' Kalman is re-run on its nights up to the snapshot end date, because later
nights re-estimate the filter's series-wide noise and IQR screen and move C_L inside the old window
(up to ~7 % for the CHM15k v2.2; leaving out its single 08-14 night alone moves it up to 16 %).
That recipe is gated (1e-4) against the payload on the series that have no later night.

L'ancienne chaîne v2 (`sites.py`, `l1_l2_calib.py`, `build_l1_l2_dashboard.py`,
`render_l1_l2_dashboard.py`, pages `index_<site>.html`) reste fonctionnelle mais n'est plus le
point d'entrée.
