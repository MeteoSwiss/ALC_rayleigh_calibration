"""
homogenize_cl61_daily.py — bring the existing Cloudnet CL61 raw into ONE canonical daily
file per day:  <site>/<YYYYMMDD>.nc  (xarray.open_mfdataset(combine='by_coords'), per the
user's notebook download_raw_cl61_cloudnet_hyy.ipynb).

The current tree is a mix of intermediates per day:
  - a per-day FOLDER  <site>/<YYYYMMDD>/*.nc        (the raw 'live' files)
  - 24 HOURLY files   <site>/<YYYYMMDD>_HHMM.nc     (hourly concatenations; hyy)
  - sometimes already  <site>/<YYYYMMDD>.nc         (the daily file we want)

For each day we build <YYYYMMDD>.nc from the cheapest complete source (hourly files if
present -> far fewer opens; else the folder's live files), verify it, then DELETE the
redundant sources (the hourly files and the folder). Resumable; parallel across days.

Usage:  python homogenize_cl61_daily.py [site ...]      (default: lindenberg hyy)
        KEEP_FOLDERS=1  keep the source folders/hourly files (no deletion)
        WORKERS=N       parallelism (default 8)
"""
from __future__ import annotations
import os, sys, glob, shutil
from pathlib import Path
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed

ROOT = Path("R:/CL61/RAW_cloudnet_dl")
SITES = ["lindenberg", "hyy"]


def _verify(path) -> int:
    import xarray as xr
    try:
        ds = xr.open_dataset(path)
        n = int(ds.sizes.get("time", 0))
        ds.close()
        return n
    except Exception:
        return 0


def concat_day(args):
    """Build <day>.nc for one day from hourly files (preferred) or the folder; clean up."""
    site, day = args
    import warnings
    warnings.filterwarnings("ignore")
    import xarray as xr
    import dask
    dask.config.set(scheduler="synchronous")

    base = ROOT / site
    out = base / f"{day}.nc"
    hourly = sorted(glob.glob(str(base / f"{day}_*.nc")))          # YYYYMMDD_HHMM.nc
    folder = base / day
    ffiles = sorted(glob.glob(str(folder / "*.nc"))) if folder.is_dir() else []
    keep = os.environ.get("KEEP_FOLDERS", "") == "1"

    def cleanup():
        if keep:
            return
        for h in hourly:
            try:
                os.remove(h)
            except OSError:
                pass
        if folder.is_dir():
            shutil.rmtree(folder, ignore_errors=True)

    # already have a valid daily file -> just drop the redundant sources
    if out.exists() and out.stat().st_size > 0:
        if _verify(out) > 0:
            had_src = bool(hourly or ffiles)
            cleanup()
            return (site, day, "exists+cleaned" if (had_src and not keep) else "exists")
        out.unlink(missing_ok=True)  # corrupt -> rebuild below

    src = hourly if hourly else ffiles
    if not src:
        return (site, day, "nofiles")
    try:
        if len(src) == 1:
            ds = xr.open_dataset(src[0])
        else:
            ds = xr.open_mfdataset(src, combine="by_coords", data_vars="minimal",
                                   coords="minimal", compat="override")
        tmp = str(out) + ".tmp"
        ds.to_netcdf(tmp)
        ds.close()
        if _verify(tmp) > 0:
            os.replace(tmp, out)
            cleanup()
            return (site, day, f"ok({len(src)} {'hourly' if hourly else 'live'})")
        Path(tmp).unlink(missing_ok=True)
        return (site, day, "empty")
    except Exception as e:  # noqa: BLE001
        return (site, day, f"ERR {type(e).__name__}: {str(e)[:70]}")


def days_for_site(site):
    """All days that have a folder and/or hourly files (i.e. still need consolidating)."""
    days = set()
    for d in glob.glob(str(ROOT / site / "2*")):
        name = Path(d).name
        if Path(d).is_dir() and len(name) == 8 and name.isdigit():
            days.add(name)
        elif name.endswith(".nc") and "_" in name:          # hourly YYYYMMDD_HHMM.nc
            pre = name.split("_")[0]
            if len(pre) == 8 and pre.isdigit():
                days.add(pre)
    return sorted(days)


def main():
    sites = [a for a in sys.argv[1:] if not a.startswith("-")] or SITES
    jobs = []
    for site in sites:
        for day in days_for_site(site):
            jobs.append((site, day))
    workers = int(os.environ.get("WORKERS", "8"))
    print(f"homogenize: {len(jobs)} days across {sites}, {workers} workers, "
          f"delete_sources={os.environ.get('KEEP_FOLDERS','')!='1'}", flush=True)
    nok = 0
    counts = defaultdict(int)
    with ProcessPoolExecutor(max_workers=workers) as ex:
        futs = {ex.submit(concat_day, j): j for j in jobs}
        for f in as_completed(futs):
            site, day, status = f.result()
            counts[status.split("(")[0].split(":")[0]] += 1
            if status.startswith(("ok", "exists")):
                nok += 1
            print(f"  {site}/{day}: {status}", flush=True)
    print(f"HOMOGENIZE_DONE: {nok}/{len(jobs)} ok  | {dict(counts)}")


if __name__ == "__main__":
    main()
