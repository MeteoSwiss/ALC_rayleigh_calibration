# Bascule v2.2 sur zueub434 — instructions pour une session Claude Code

À exécuter **sur `zueub434.meteoswiss.ch`**, dans un clone du dépôt, branche `rayleigh-availability`
(commit `c3386ba` ou plus récent). Le calcul est déjà fait : l'archive v2.2 complète existe sur
balfrin et le dashboard est **déjà publié** sur EWC. Ce document ne fait que ramener le serveur
opérationnel au même millésime et faire en sorte que le quotidien continue tout seul.

> **Contexte à ne pas perdre.** Ce qui est publié aujourd'hui est un *instantané* : les 263 574
> payloads du panneau interactif sont dans le bucket, mais le cron du serveur tourne encore
> `eprof_v2` et ne génère aucun payload. Tant que l'étape 5 n'est pas faite, le site se fige
> progressivement : les courbes avanceront (elles viennent des CSV) mais le panneau journalier
> restera bloqué à la date de la release.

---

## 0. Orientation (toujours en premier)

```bash
source ops/config.sh                 # exporte tous les ALC_* + creds S3 + ALC_VENV
source "$ALC_VENV/bin/activate"
git -C "$ALC_REPO" rev-parse --abbrev-ref HEAD    # QUELLE branche est opérationnelle ?
```

⚠️ **Avant tout `git checkout`/`pull` : sauvegarder le census.** `ALC_CENSUS` pointe *dans* le
checkout et est réécrit chaque jour par `refresh_census.py` ; un checkout ferait disparaître toute
station mise en service depuis le dernier commit.

```bash
cp -a "$ALC_CENSUS" ~/census_backup_$(date +%Y%m%d).json
```

Pour devenir `rem` : **uniquement** `sudo su - rem`, et les commandes passent par **STDIN**
(`echo 'cmd' | sudo su - rem`), jamais `-c`. Le shell de `rem` a `noclobber` → `rm -f fichier`
avant toute redirection sur un fichier existant.

---

## 1. Mettre à jour le code (il est DÉJÀ sur zueub434)

Le dépôt opérationnel vit sur le serveur : **rien ne se transfère depuis balfrin côté code**, c'est
un simple `git pull` dans le clone existant. Balfrin n'a servi qu'au calcul, et seules les
**données** en viennent (étape 2).

```bash
cd "$ALC_REPO"
git fetch origin
git rev-parse --abbrev-ref HEAD          # branche réellement déployée ?
git checkout rayleigh-availability && git pull
git log --oneline -1                     # doit être >= 955efb4
```

⚠️ Si le clone est sur une autre branche (`main` est gelée et très en retard — chemins `/data/pay`,
`eprof_v1.2`, pas de `publish.sh`), **ne pas basculer de branche à l'aveugle** : vérifier d'abord ce
que le cron lance réellement, et décider avec l'opérateur si on bascule la branche ou si on
rapatrie les commits. Une bascule de branche change aussi le census (voir §0).

Ce que ce code apporte et que l'ancien n'a pas :
* la colonne `version` dans `CSV_FIELDS` (sans elle le dashboard ne sait pas nommer le millésime) ;
* le correctif « jour vide » de la sensibilité (sans lui, 15 flux sur 434 produisent le produit) ;
* le panneau journalier interactif + la carte C vs base de nuage ;
* `window.__payloadBase` (le panneau va chercher ses données dans le bucket).

---

## 2. Transférer l'archive v2.2 depuis balfrin

Depuis zueub434 (`ssh -n -o BatchMode=yes balfrin` fonctionne dans ce sens) :

| | |
|---|---|
| source | `balfrin:/scratch/mch/mhrvo/E_PROFILE_calout_v22_rel` |
| taille | **11 Go** (dont 2,8 Go de PNG OmB/sens, 957 Mo de CSV+npz+NetCDF) |
| contenu | 434 flux, 429 `_omb_cache.npz`, 429 `_sens_cache.npz`, 1 686 NetCDF annuels |

```bash
DEST=/data/zue/E_PROFILE/ALC/Calibration/ALC_calibration_v2.2
mkdir -p "$DEST"
rsync -a --partial --inplace --info=progress2 \
  balfrin:/scratch/mch/mhrvo/E_PROFILE_calout_v22_rel/ "$DEST/"
```

**Les deux `.npz` par flux DOIVENT voyager.** Le garde-fou de régression du runner compare le cache
au CSV et **saute** le produit au lieu de le reconstruire : une archive transférée sans ses caches
fige OmB et la sensibilité *définitivement*.

**`classification/` n'est pas dans la nouvelle archive** (0 répertoire — c'est normal, ce produit ne
dépend pas de la méthode et n'a pas été régénéré). Deux options, dans cet ordre de préférence :

1. ne rien faire : le build liste les images de classification **dans le bucket** via
   `_bucket_keys()` (vérifié : 570 objets pour la seule station Payerne A) ;
2. si tu veux les avoir en local, les copier depuis l'ancienne archive :
   `rsync -a <ancienne>/*/classification/ "$DEST"/…` — jamais les régénérer.

Contrôle de complétude avant de continuer :

```bash
for f in _cal.csv _kalman.csv _hk.csv _status.csv _omb.csv _sens.csv; do
  echo -n "$f: "; ls "$DEST"/*/*$f 2>/dev/null | wc -l
done
echo -n "_omb_cache.npz: "; ls "$DEST"/*/_omb_cache.npz | wc -l
echo -n "_sens_cache.npz: "; ls "$DEST"/*/_sens_cache.npz | wc -l   # doit valoir 429, pas 15
```

---

## 3. La bascule : **une seule** édition de `ops/config.sh`

```bash
cp -a ops/config.sh ops/config.sh.pre_v22_$(date +%Y%m%d)
```

Puis, dans **le même commit / la même édition** (jamais l'un sans l'autre) :

```bash
export ALC_MOLECULAR_METHOD=eprof_v2.2
export ALC_WV_SPECTRUM='{"CL61": [910.55, 0.188]}'
export ALC_FULLCAL_DIR=/data/zue/E_PROFILE/ALC/Calibration/ALC_calibration_v2.2
```

⚠️ **Méthode et arbre de sortie basculent ensemble.** Avec la méthode en v2.2 et l'ancien arbre, le
cron fusionne des millésimes dans les mêmes CSV et `_preserve_existing_rows` conserve les deux par
(méthode, fenêtre) — **irréversible**.

Ne pas toucher `options.json` (c'est `ops/config.sh` qui pilote l'opérationnel).
Laisser `ALC_DARK_PROFILE` **non défini** (décision : réseau homogène).

---

## 4. Combler le trou D_END → hier

L'archive s'arrête au **2026-08-13**. Il faut rattraper les jours suivants **sous le même `flock`
que le cron**, sinon deux passes 433-flux se télescopent :

```bash
flock /var/lock/alc_daily.lock \
  python scripts/run_network_calibration.py \
    --start 20260814 --end $(date -u -d yesterday +%Y%m%d) \
    --sens --omb --per-type 0 --ignore-coverage --workers 8
```

(adapter le chemin du lock à celui qu'utilise `ops/run_daily.sh` — le vérifier d'abord).

---

## 5. Ajouter l'étape « payloads » au flux quotidien — **c'est l'étape qui rend le site vivant**

Sans elle, le panneau interactif reste figé à la date de la release. Dans `ops/ops_daily.py`, entre
l'étape *calibrate* et l'étape *build dashboard*, pour chaque jour cible :

```bash
python scripts/build_station_dashboard.py \
  $(for d in "$ALC_FULLCAL_DIR"/*/; do echo --key $(basename "$d"); done) \
  --start <D> --end <D> \
  --cal-dir "$ALC_FULLCAL_DIR" --l1-root "$ALC_L1_ROOT" \
  --out "$ALC_DASH_DIR" --payloads --no-pages --workers 8
```

Ordre de grandeur mesuré : **~456 payloads/jour, ~40 Mo/jour** (~15 Go/an), ~10,5 s CPU par payload.

Puis les pousser dans le bucket, avec le reste des images :

```bash
aws --profile ewc --endpoint-url "$ALC_S3_ENDPOINT" \
    s3 sync "$ALC_DASH_DIR/data" "s3://$ALC_S3_BUCKET/data" --only-show-errors --size-only
```

À décider avec l'opérateur : **horizon de rétention**. Aujourd'hui le store grossit sans limite
(~15 Go/an sur une base de 102 Go). Si un horizon est adopté, ajouter l'élagage des objets plus
anciens que la fenêtre, sur le modèle de l'élagage des PNG diag déjà présent dans `publish.sh`.

---

## 6. Le timer (cron **ou** systemd — vérifier lequel est réellement actif)

```bash
crontab -l | grep -i alc
sudo systemctl list-timers | grep -i alc
```

Le dépôt documente `cron 0 15 * * *` → `ops/run_daily.sh`. Si c'est un **timer systemd** :

```bash
systemctl cat alc-daily.timer alc-daily.service    # lire AVANT d'éditer
sudoedit /etc/systemd/system/alc-daily.service     # si l'unité doit pointer ailleurs
sudo systemctl daemon-reload
sudo systemctl restart alc-daily.timer
systemctl status alc-daily.timer --no-pager
```

Contraintes à respecter quel que soit le mécanisme :
* **après 03:30 UTC** — le L1 d'un jour n'arrive que le lendemain matin ; un déclenchement avant
  l'aube ne trouve aucune donnée pour « hier » ;
* **un seul run à la fois** (`flock`), sinon deux passes réseau se télescopent ;
* garder `ALC_DAY_LAG=1` et `ALC_BACKFILL_DAYS=5` (auto-réparation).

Il y a aussi un cron `18 18 * * *` de `hem` qui produisait l'overlay v13, **retiré du code en
2026-07** : le supprimer de la crontab s'il existe encore.

---

## 7. Publier sur EWC

Le bucket est **déjà** à jour (payloads + images OmB/sens v2.2) et **CORS est déjà configuré**
(règle GET/HEAD, `ops/cscs/bucket_cors.json`) : ne pas la retirer, le panneau cesserait de charger
ses données.

Le seul piège du déploiement HTML, appris à la dure : **les permissions**.

```bash
bash ops/publish.sh          # images -> bucket, HTML -> VM
```

Si le HTML est déployé par extraction d'une archive construite ailleurs, **normaliser les droits**,
sinon nginx renvoie 403/404 (une archive faite sur balfrin porte du 640) :

```bash
ssh -i ~/.ssh/EWC hem@136.156.139.31 '
  find /var/www/alc -type d -exec chmod 755 {} + ;
  find /var/www/alc -type f -exec chmod 644 {} +'
```

`ops/publish.sh` renvoie parfois **rc=2** sur le rsync HTML alors que le HTML est bien arrivé : si
le site public est en retard, relancer `bash ops/publish.sh`, ne pas déboguer rc=2.

---

## 8. Vérification (le lendemain, puis à J+7)

```bash
grep -c REGRESSION-GUARD <log du build>        # doit valoir 0
ls "$ALC_DASH_DIR"/stations/*.html | wc -l     # >= 434
```

Dans un navigateur, sur trois stations de types différents (CHM15k / CL31 / CL61) :
* la console affiche `[ALC] rangesync ready` ;
* le panneau journalier se dessine **et affiche la date d'hier** (c'est le test de l'étape 5) ;
* changer la période affiche `unchanged: none` ;
* les trois cartes du bas (OmB, classification, sensibilité) sont présentes.

---

## 9. Retour arrière

```bash
cp -a ops/config.sh.pre_v22_<date> ops/config.sh   # restaure méthode + spectre + arbre
# rebuild COMPLET (jamais partiel avec rsync --delete), puis publish
```

Rien dans cette release ne supprime d'objet du bucket : l'ancien site redevient cohérent dès que
son HTML est restauré, les payloads v2.2 restant simplement non référencés. Un instantané du
docroot d'avant bascule existe déjà sur la VM : `/home/hem/alc_pre_v22` (à supprimer une fois la
release acceptée).

---

## Ce qu'il ne faut PAS faire

* ne pas régénérer `classification/` (indépendant de la méthode, coûteux, déjà dans le bucket) ;
* ne pas transférer l'archive **sans** les `.npz` (fige OmB et la sensibilité pour toujours) ;
* ne pas basculer la méthode sans basculer l'arbre dans la même édition ;
* ne pas écrire dans `/mnt/amaroc_data/alc_calib` (ancien emplacement, plafond d'inodes) ;
* ne pas retirer la règle CORS du bucket ;
* ne pas livrer les constantes au hub L1→L2 : c'est une étape coordonnée **séparée**, avec préavis
  aux utilisateurs, et c'est elle — pas cette release — qui fera quitter aux CL31/CL51 leur
  défaut 1e8.
