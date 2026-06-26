# E-PROFILE ALC dashboard — European Weather Cloud deployment

**Deployed 2026-06-26.** Live: **https://alc-calib.ch-meteoswiss-emermet.f.ewcloud.host/**

Publish-only port of the calibration monitoring dashboard onto the EWC: the site is built where the
calibration runs (on-prem / CSCS balfrin); only *publishing* lives on the EWC. The bulky per-night
diagnostic images go to a public S3 bucket; the static HTML is served by a small VM. The HTML points at
the images by absolute URL, so the browser fetches them straight from object storage.

```
 build host (balfrin / on-prem)                EUROPEAN WEATHER CLOUD (cci2, meteoswiss-emermet)
  build_dashboard.py --img-base-url ┌── rsync HTML ──▶ VM alc-calib (136.156.139.31)
   -> dashboard_ewc/ (HTML+images) ─┤                  nginx /var/www/alc, Let's Encrypt TLS (443)
                                     └── aws s3 sync ─▶ bucket eprofile-alc-dashboard (public-read)
                                        (diag/ombsens)   object-store.os-api.cci2.ecmwf.int
 browser → HTML from the VM, images straight from the bucket (no CORS needed for <img>)
```

## Components

- **Bucket** `eprofile-alc-dashboard` on cci2 (`object-store.os-api.cci2.ecmwf.int`), public-read
  (anonymous `s3:GetObject`, no listing — see `ops/ewc_bucket_policy.json`). Public base URL =
  `https://object-store.os-api.cci2.ecmwf.int/eprofile-alc-dashboard/` → this is `ALC_IMG_BASE_URL`.
  S3 keys live in Morpheus → Tools → Cypher; **region-bound to cci2**; use **path-style** addressing.
- **VM** `alc-calib` (Ubuntu 22.04, public IP `136.156.139.31`, security group `ssh-https` = 22+443).
  Login user **`hem`** (passwordless sudo via `/etc/sudoers.d/90-hem`). nginx docroot `/var/www/alc`.
  DNS `alc-calib.ch-meteoswiss-emermet.f.ewcloud.host`.
- **TLS**: Let's Encrypt ECC cert via **acme.sh** TLS-ALPN-01 on 443 (certbot's standalone does *not*
  support tls-alpn-01; only http-01/dns-01 — and port 80 is closed). Auto-renews via acme.sh cron with
  stop/start-nginx hooks; next ~2026-08-25. nginx redirects 80→443.
- **Code** (this change): `ALC_IMG_BASE_URL` / `--img-base-url` makes the diagnostic, OmB/sensitivity and
  flag-example image references absolute under the bucket base (`monitoring/config.py`, `render.py`,
  `static/diag.js`, `scripts/build_dashboard.py`); empty = unchanged relative behaviour.
- **Publish step**: `ops/publish.sh` (images→bucket via rclone/aws, HTML→VM via rsync `--chmod=D755,F644`)
  + `ops/ops_daily.py` `publish()` (gated on `ALC_PUBLISH=1`, non-fatal) + `ops/config.sh` vars.

## Build + first publish (run on balfrin, 2026-06-26)

```bash
source ~/ALC_rayleigh_calibration/cscs_env_2025_2026.sh
source ~/miniforge3/etc/profile.d/conda.sh && conda activate alc
python scripts/build_dashboard.py --fullcal "$ALC_FULLCAL_DIR" --out /scratch/mch/mhrvo/dashboard_ewc \
  --img-base-url https://object-store.os-api.cci2.ecmwf.int/eprofile-alc-dashboard/ \
  --opcoeff /scratch/mch/mhrvo/opcoeff.csv --oldray "$ALC_OLDRAY_DIR" \
  --l2dir /scratch/mch/mhrvo/E_PROFILE_L2 --workers 8           # 427 stations, ~5 min
# HTML -> VM (aws/rsync profile 'ewc' uses ca_bundle=/etc/ssl/ca-bundle.pem, path-style, region us-east-1)
rsync -az --delete --chmod=D755,F644 -e "ssh -i ~/.ssh/ewc_vm_key" \
  --exclude 'diag/' --exclude 'ombsens/' --exclude 'flagex/' --exclude 'calib_index.sqlite' \
  /scratch/mch/mhrvo/dashboard_ewc/ hem@136.156.139.31:/var/www/alc/
# images -> bucket (recent subset first; full set is ~212 GB)
aws --profile ewc --endpoint-url https://object-store.os-api.cci2.ecmwf.int s3 sync \
  /scratch/mch/mhrvo/dashboard_ewc/diag s3://eprofile-alc-dashboard/diag --follow-symlinks \
  --exclude '*' --include '*_202604*' --include '*_202605*' --include '*_202606*'
aws ... s3 sync /scratch/mch/mhrvo/dashboard_ewc/ombsens s3://eprofile-alc-dashboard/ombsens --follow-symlinks
```

## Maintenance / follow-ups

- **Daily publishing**: set `ALC_IMG_BASE_URL`, `ALC_S3_REMOTE`/`ALC_S3_BUCKET` (or the aws profile),
  `ALC_VM_RSYNC_TARGET`, `ALC_PUBLISH=1` in `ops/config.sh`; `ops_daily.py` then publishes after each
  incremental build. (A one-time full rebuild with `--img-base-url` is needed at cutover so every page
  carries absolute URLs.)
- **Image backfill**: the first push uploaded only **recent diag images (Apr–Jun 2026, ~35 GB of 212 GB)**.
  Backfill older months by re-running the `aws s3 sync` of `diag/` without the date `--include` filters.
- **TLS** renews automatically (acme.sh); no action unless the DNS name changes.
- **Rotate** the bucket S3 keys in Cypher (they passed through a chat session during setup).
- **Trust note**: this Windows dev box's curl/schannel distrusts LE's ECDSA root (stale local store) —
  the cert is valid (verified from the VM); browsers are fine.
