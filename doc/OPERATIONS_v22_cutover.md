# v2.2 cutover runbook (host `zueub434`)

Take the operational network from **`eprof_v2`** to **`eprof_v2.2` + the CL61 constructor
water-vapour spectrum**, by transferring a recomputed archive from CSCS and flipping one config
block. Written to be executed step by step on `zueub434` (by a human or a Claude Code session
there); every step says what it changes and how to undo it.

**Companion documents.** `doc/OPERATIONS.md` (the standing ops guide — §10 Gotchas is assumed
knowledge here), `doc/reports/phase4_network_validation.md` §5.1 (why v2.2 is deployable),
`ops/cscs/alc_v22_release.sbatch` + `alc_v22_addons.sbatch` (how the archive was produced),
`scripts/check_v22_archive.py` (the gates that authorise the flip).

## What actually changes

| | before | after |
|---|---|---|
| Rayleigh method | `eprof_v2` (version code 200) | `eprof_v2.2` (220) |
| CL61 WV spectrum | 910.74 nm / FWHM 1.0 | **910.55 nm / FWHM 0.188** (constructor) |
| CL31 / CL51 / CHM15k spectra | unchanged | unchanged |
| dark correction | none | none (network kept homogeneous — operator decision) |
| output tree | `ALC_calibration_v2.0/` | `ALC_calibration_v2.2/` |
| `_status.csv` | 9 columns | +`mean_cloud_cover`, `cloud_cover_n`, `cloud_src` |

**Not** in this release: what E-PROFILE distributes in L2. `operational_coefficients.csv` only
*reads* the constants out of the distributed L2 files — nothing here writes back to the L1→L2 hub.
Handing the new constants to the hub is a separate, coordinated step, and it is that step (not this
one) that would take CL31/CL51 off their uncalibrated 1e8 default. No user notice is needed for the
cutover itself.

## Host gotchas that bite in this procedure

- `sudo` is **stdin-only** to `rem`: `echo "cmd" | sudo /bin/su - rem`. There is no `sudo -u rem`.
- `rem` runs with **`noclobber`**: inside its shell `>` refuses to overwrite. Use `>>`, `>|`, or a
  fresh path.
- `publish.sh` returning **rc=2** is the known rsync-protocol quirk and is **non-fatal**; the HTML
  usually lands anyway. Re-run `bash ops/publish.sh` if the live site lags.
- Internet egress (CAMS) needs the **8080 HTTP proxy**; the HTML rsync to the VM uses a **1080
  SOCKS** ProxyCommand. Different ports — do not conflate them.

---

## Step 0 — BLOCKING: establish what is actually deployed

The repository does not record which branch the server runs, and `main` is stale (2 commits, wrong
paths, `eprof_v1.2`). Find out before touching anything.

```bash
source /data/zue/E_PROFILE/ALC/Calibration/ALC_calibration_v2.0_code/ops/config.sh
cd "$ALC_REPO" && git rev-parse --abbrev-ref HEAD && git log --oneline -3 && git status --porcelain
echo "$ALC_FULLCAL_DIR"; ls -d "$ALC_FULLCAL_DIR" | head
```

Record the branch and commit in the change log. If the working tree is dirty in tracked files, stop
and find out why before continuing.

## Step 1 — BLOCKING: rescue the live census

`ALC_CENSUS` points **inside the git checkout** and `refresh_census.py` rewrites it every day, so a
`git checkout`/`reset` would silently discard every station commissioned since the last commit.

```bash
cp "$ALC_CENSUS" /data/zue/E_PROFILE/ALC/Calibration/census_live_$(date -u +%Y%m%d).json
python - <<'PY'
import json, os
c = json.load(open(os.environ["ALC_CENSUS"]))
rows = c if isinstance(c, list) else c.get("streams", [])
print("live census streams:", len({r["wmo"] + "_" + r["ident"] for r in rows}))
PY
ls -d "$ALC_FULLCAL_DIR"/*/ | wc -l      # streams that actually have history
```

Keep both numbers. Step 4 refuses to proceed if the new archive does not cover them.

## Step 2 — Transfer the archive from CSCS

The recompute lives on balfrin at `/scratch/mch/mhrvo/E_PROFILE_calout_v22_rel`. **Scratch is
periodically purged** — do this promptly after the run.

```bash
NEW=/data/zue/E_PROFILE/ALC/Calibration/ALC_calibration_v2.2
mkdir -p "$NEW"
df -h /data/zue | tail -1        # ABORT if free space < 1.5x the size of ALC_calibration_v2.0
rsync -a --info=progress2 balfrin:/scratch/mch/mhrvo/E_PROFILE_calout_v22_rel/ "$NEW/"
```

The `.npz` caches must come across. They are not optional: the daily runner compares cache against
CSV and **skips** OmB/sensitivity rather than rebuilding, so a tree without them freezes those
products permanently.

```bash
ls "$NEW"/*/_omb_cache.npz  | wc -l
ls "$NEW"/*/_sens_cache.npz | wc -l
ls -d "$NEW"/*/ | wc -l          # the three counts must agree
```

## Step 3 — Carry over what the recompute deliberately did not produce

Both are method-independent, so copying them is correct rather than lazy.

```bash
OLD="$ALC_FULLCAL_DIR"
# Cloudnet classification curtains: absent from the new tree; without this the whole gallery
# disappears from the dashboard network-wide at the flip.
for d in "$OLD"/*/classification; do
  k=$(basename "$(dirname "$d")"); [ -d "$d" ] && cp -a "$d" "$NEW/$k/"
done
# Per-night diagnostic PNGs: NOT copied on purpose (they are ~100-200 GB and already served from
# the S3 bucket). The daily run regenerates them from the flip date onward.
```

## Step 4 — Run the gates

```bash
cd "$ALC_REPO"
python scripts/check_v22_archive.py \
    --new "$NEW" \
    --ref  balfrin:/scratch/... (copy the reference locally, or run this ON balfrin) \
    --old "$OLD" \
    --census "$ALC_CENSUS" \
    --json /tmp/v22_gates.json
```

Exit 0 is the authorisation to flip. Anything else stops the procedure. `SKIPPED` is not a pass —
supply the missing tree and re-run. In particular G3 must confirm that every stream in the **live**
census and in the **current** tree exists in the new one; streams commissioned after the CSCS run
will be missing and must be calibrated locally before the flip:

```bash
# example, for the stragglers G3 lists
ALC_FULLCAL_DIR="$NEW" ALC_MOLECULAR_METHOD=eprof_v2.2 \
ALC_WV_SPECTRUM='{"CL61": [910.55, 0.188]}' \
python scripts/run_network_calibration.py --stream <key> --start 20250101 --end <D_END> \
       --methods rayleigh,cloud --sens --omb --force
```

## Step 5 — Fill the seam (D_END → yesterday)

The CSCS archive ends at the date printed in `RUN_PROVENANCE.json`. Fill the remaining days on the
server **under the same lock as the cron**, or two 433-stream runs will collide.

```bash
flock -n /tmp/alc_daily.lock -c '
  ALC_FULLCAL_DIR="'"$NEW"'" ALC_MOLECULAR_METHOD=eprof_v2.2 \
  ALC_WV_SPECTRUM='"'"'{"CL61": [910.55, 0.188]}'"'"' \
  python '"$ALC_REPO"'/ops/ops_daily.py --start <D_END+1> --end <yesterday> \
         --no-dashboard --no-publish'
```

Then re-run Step 4's gates. G4 (seam) in the plan is this check: pick 5 stations × 3 days near
D_END, re-run them through the daily driver and confirm the flags match and |ΔC/C| < 0.1 %. Any
larger difference must be attributed (the CSCS run used the 0.4° monthly CAMS plus dailies; the
server uses its own daily cache) and written down here.

## Step 6 — Dry-run the dashboard OFF the live site

```bash
ALC_FULLCAL_DIR="$NEW" ALC_DASHBOARD_DIR=/data/zue/E_PROFILE/ALC/Calibration/dashboard_v22_dryrun \
python scripts/build_dashboard.py --fullcal "$NEW" --out /data/zue/.../dashboard_v22_dryrun
```

Expect ≥ 434 station pages and **zero** `REGRESSION-GUARD` lines. Open a few pages: the
availability card should now carry the cloud-cover row, and the constants should look like the
station's history, not like a step change on every stream.

## Step 7 — THE FLIP

```bash
cp "$ALC_REPO/ops/config.sh" "$ALC_REPO/ops/config.sh.pre_v22_$(date -u +%Y%m%d)"
# uncomment the four lines of the "v2.2 CUTOVER BLOCK" -- all of them, together
$EDITOR "$ALC_REPO/ops/config.sh"
source "$ALC_REPO/ops/config.sh"
echo "$ALC_FULLCAL_DIR"; echo "$ALC_MOLECULAR_METHOD"; echo "$ALC_WV_SPECTRUM"
```

Then a **full** rebuild (never `--changed-only` across a tree change — the page set is regenerated
from a different archive) and publish:

```bash
rm -f "$ALC_DASHBOARD_DIR/.last_build"
python scripts/build_dashboard.py --changed-only=false 2>/dev/null || \
python scripts/build_dashboard.py
bash ops/publish.sh          # rc=2 is the known non-fatal rsync quirk
```

## Step 8 — Rollback (rehearse this BEFORE you need it)

The repository had no rollback procedure; this is it. It restores config, not data, because the
v2.0 tree was never written to.

```bash
cp "$ALC_REPO/ops/config.sh.pre_v22_<date>" "$ALC_REPO/ops/config.sh"
source "$ALC_REPO/ops/config.sh"
rm -f "$ALC_DASHBOARD_DIR/.last_build"
python scripts/build_dashboard.py && bash ops/publish.sh
```

The live site must show v2.0 numbers again within ~15 minutes. Gate G5 asks you to *execute* this
once on the dry-run directory and flip forward again, so the procedure is proven rather than
believed.

## Step 9 — Verify at D+1 and D+7

- the 15:00 cron ran green, and its log shows `ALC_FULLCAL_DIR=.../ALC_calibration_v2.2`;
- **zero** `REGRESSION-GUARD` lines (a single one means a `.npz` did not travel — copy it from
  balfrin, never delete the CSV);
- station page count stable, no station lost its history;
- the new days' constants are continuous with the CSCS archive across the seam;
- `operational_coefficients.csv` still updates (it reads the distributed L2 and is unaffected by
  the flip — if it changes, something else did).

## Deliberately out of scope

Uncertainty-weighted Kalman, the two-pass altitude correction, any network dC/dCBH correction,
regenerating 12 months of diagnostic images, and handing the constants to the L1→L2 hub.
