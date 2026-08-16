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

## Provenance des données

Runs par variante sous `C:/DATA/Projects/202606_E-PROFILE_calibration/` (`calout_v22_04` réseau,
`diag_v22_*` locaux — voir `variants_v3.py`), extraits Hopkin par configuration dans
`cloud_profile_dump_<tag>/`. L1 : `D:/E-PROFILE_L1_2026` ; CAMS **0,4°** uniquement
(`A:/CAMS_Monthly_04;D:/CAMS_daily` — jamais `D:/CAMS` ni `D:/CAMS_run_v20`, monolithes 1°).
L2 tel que distribué : CHM15k **et CL61** = calibration Rayleigh v1.0 opérationnelle (activée
mi-2026 pour le CL61 — d'où ses constantes observées 1,00 sur 48 j puis 2,1257 dès la mi-juin) ;
CL31 = non calibré (constante par défaut).

L'ancienne chaîne v2 (`sites.py`, `l1_l2_calib.py`, `build_l1_l2_dashboard.py`,
`render_l1_l2_dashboard.py`, pages `index_<site>.html`) reste fonctionnelle mais n'est plus le
point d'entrée.
