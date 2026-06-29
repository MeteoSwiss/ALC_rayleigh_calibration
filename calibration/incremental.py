"""Per-day cached intermediates for the OmB and sensitivity products, so a daily
run can update a FULL-PERIOD snapshot by adding only yesterday's contribution --
never re-reading the raw L1 archive.

Design
------
Both sensitivity and OmB are *windowed* products: their headline figures are an
aggregate over many days.  Re-reading the whole archive every day is exactly what
we want to avoid.  Instead each station keeps a compact per-day cache:

* sensitivity ``<key>_sens_cache.npz`` : one ``bmin_night`` / ``bmin_day`` column
  (min detectable backscatter vs altitude) per day.  Aggregation is the existing
  :func:`calibration.sensitivity.network.combine_sens_results` (concatenate the
  daily columns, take the per-altitude nanmedian).  EXACT -- the cached-and-
  aggregated result equals a single direct multi-day pass.

* OmB ``<key>_omb_cache.npz`` : per-CAMS-time columns of the bias and the
  interpolated observation/CAMS profiles.  Aggregation concatenates along the
  CAMS-time axis and recomputes the scalar/profile statistics exactly as
  :func:`calibration.omb.omb.compute_omb` does over a single window.

A daily run computes yesterday's column(s) from yesterday's L1 only, appends them
to the cache (de-duplicating by date), and re-derives the snapshot from the cache.
"""

from __future__ import annotations

from pathlib import Path
from typing import Optional

import numpy as np

# ----------------------------------------------------------------------------
# sensitivity
# ----------------------------------------------------------------------------

def sens_cache_path(key_dir: Path) -> Path:
    return Path(key_dir) / "_sens_cache.npz"


def sens_cache_update(key_dir: Path, day_result) -> None:
    """Merge a (1-day-or-more) ``SensResult`` into the station's sensitivity cache.

    ``day_result.dates`` are datetime64[D]; existing dates are overwritten so a
    forced re-run of a day replaces (not duplicates) its column.
    """
    p = sens_cache_path(key_dir)
    z = np.asarray(day_result.z_ctr, dtype="float64")
    dates = np.asarray(day_result.dates, dtype="datetime64[D]")
    bn = np.asarray(day_result.bmin_night, dtype="float32")
    bd = np.asarray(day_result.bmin_day, dtype="float32")
    wl = float(getattr(day_result, "wavelength", float("nan")) or float("nan"))
    if p.exists():
        c = np.load(p, allow_pickle=False)
        if "wavelength" in c and not np.isfinite(wl):
            wl = float(c["wavelength"])
        if c["z_ctr"].size == z.size and np.allclose(c["z_ctr"], z, equal_nan=True):
            old_dates = c["dates"].astype("datetime64[D]")
            keep = ~np.isin(old_dates, dates)            # drop days we are replacing
            dates = np.concatenate([old_dates[keep], dates])
            bn = np.concatenate([c["bmin_night"][:, keep], bn], axis=1)
            bd = np.concatenate([c["bmin_day"][:, keep], bd], axis=1)
        # else: z grid changed (range-grid change) -> start a fresh cache
    order = np.argsort(dates)
    Path(key_dir).mkdir(parents=True, exist_ok=True)
    np.savez(p, z_ctr=z, dates=dates[order].astype("datetime64[D]"),
             bmin_night=bn[:, order], bmin_day=bd[:, order], wavelength=wl)


def sens_cache_aggregate(key_dir: Path):
    """Full-period ``SensResult`` rebuilt from the cache (ext/mass/ICAO/sigma
    recomputed from the per-day median), or ``None`` if the cache is empty."""
    from calibration.sensitivity.network import SensResult, combine_sens_results
    p = sens_cache_path(key_dir)
    if not p.exists():
        return None
    c = np.load(p, allow_pickle=False)
    if c["dates"].size == 0:
        return None
    seed = SensResult(wavelength=float(c["wavelength"]) if "wavelength" in c else float("nan"),
                      z_ctr=c["z_ctr"], dates=c["dates"].astype("datetime64[D]"),
                      bmin_night=c["bmin_night"], bmin_day=c["bmin_day"])
    return combine_sens_results([seed])     # recomputes the derived fields


# ----------------------------------------------------------------------------
# OmB
# ----------------------------------------------------------------------------

def omb_cache_path(key_dir: Path) -> Path:
    return Path(key_dir) / "_omb_cache.npz"


def _omb_part_from_result(res) -> dict:
    """Extract the per-time columns we need to aggregate later from one OmBResult.

    We keep, per CAMS time step: the CAMS backscatter profile, the per-source
    interpolated observation and bias (all (n_lev, n_cams)), the source obs_mean
    (n_cams, n_r) for the pcolor panel, plus the time and the (n_lev, n_cams)
    altitude grid so the vertical axis can be re-meaned over the full period.
    """
    srcs = list(res.bias.keys())
    # obs_mean only carries the BASE sources (e.g. 'op','ours'); the water-vapor
    # variants ('op_wv','ours_wv') exist for bias/obs_interp/prof/scalar but not
    # for obs_mean -- so key the obsmean columns by obs_mean's own keys.
    mean_srcs = list(res.obs_mean.keys())
    full_srcs = list(getattr(res, "obs_full", {}).keys())            # unscreened obs for the pcolor
    n = len(np.asarray(res.time_cams))
    cloud_base = np.asarray(getattr(res, "cloud_base", np.full(n, np.nan)), dtype="float32")
    if cloud_base.size != n:
        cloud_base = np.full(n, np.nan, dtype="float32")
    return dict(
        wavelength=float(res.wavelength),
        time_cams=np.asarray(res.time_cams),
        range_mean=np.asarray(res.range_mean, dtype="float64"),
        z_cams=np.asarray(res.z_cams, dtype="float64"),               # (n_lev,) mean for THIS part
        cams_beta=np.asarray(res.cams_beta, dtype="float32"),         # (n_lev, n_cams)
        srcs=np.array(srcs),
        mean_srcs=np.array(mean_srcs),
        full_srcs=np.array(full_srcs),
        cloud_base=cloud_base,                                        # (n_cams,) cloud base AGL [m]
        **{f"bias__{k}": np.asarray(res.bias[k], dtype="float32") for k in srcs},
        **{f"obsint__{k}": np.asarray(res.obs_interp[k], dtype="float32") for k in srcs},
        **{f"obsmean__{k}": np.asarray(res.obs_mean[k], dtype="float32") for k in mean_srcs},
        **{f"obsfull__{k}": np.asarray(res.obs_full[k], dtype="float32") for k in full_srcs},
    )


def omb_cache_update(key_dir: Path, res) -> None:
    """Append one window's OmB columns to the station cache, de-duplicating by the
    CAMS timestamp so a forced re-run replaces rather than doubles those steps."""
    p = omb_cache_path(key_dir)
    part = _omb_part_from_result(res)
    if p.exists():
        c = dict(np.load(p, allow_pickle=True))
        if (c.get("range_mean") is not None and c["range_mean"].size == part["range_mean"].size
                and c["cams_beta"].shape[0] == part["cams_beta"].shape[0]
                and list(c["srcs"]) == list(part["srcs"])
                and list(c.get("mean_srcs", c["srcs"])) == list(part["mean_srcs"])):
            old_t = c["time_cams"]
            keep = ~np.isin(old_t, part["time_cams"])
            part["time_cams"] = np.concatenate([old_t[keep], part["time_cams"]])
            part["cams_beta"] = np.concatenate([c["cams_beta"][:, keep], part["cams_beta"]], axis=1)
            # z_cams kept as a running per-level mean weighted by column count
            for k in list(part["srcs"]):
                part[f"bias__{k}"] = np.concatenate([c[f"bias__{k}"][:, keep], part[f"bias__{k}"]], axis=1)
                part[f"obsint__{k}"] = np.concatenate([c[f"obsint__{k}"][:, keep], part[f"obsint__{k}"]], axis=1)
            for k in list(part["mean_srcs"]):
                part[f"obsmean__{k}"] = np.concatenate([c[f"obsmean__{k}"][keep], part[f"obsmean__{k}"]], axis=0)
            for k in list(part.get("full_srcs", [])):
                if f"obsfull__{k}" in c:
                    part[f"obsfull__{k}"] = np.concatenate([c[f"obsfull__{k}"][keep], part[f"obsfull__{k}"]], axis=0)
            if "cloud_base" in c and "cloud_base" in part:
                part["cloud_base"] = np.concatenate([c["cloud_base"][keep], part["cloud_base"]])
    order = np.argsort(part["time_cams"])
    part["time_cams"] = part["time_cams"][order]
    part["cams_beta"] = part["cams_beta"][:, order]
    for k in list(part["srcs"]):
        part[f"bias__{k}"] = part[f"bias__{k}"][:, order]
        part[f"obsint__{k}"] = part[f"obsint__{k}"][:, order]
    for k in list(part["mean_srcs"]):
        part[f"obsmean__{k}"] = part[f"obsmean__{k}"][order]
    for k in list(part.get("full_srcs", [])):
        part[f"obsfull__{k}"] = part[f"obsfull__{k}"][order]
    if "cloud_base" in part:
        part["cloud_base"] = part["cloud_base"][order]
    Path(key_dir).mkdir(parents=True, exist_ok=True)
    np.savez(p, **part)


def omb_cache_aggregate(key_dir: Path):
    """Full-period ``OmBResult`` rebuilt from the cache, with scalar/prof recomputed
    over all cached CAMS times exactly as compute_omb does, or ``None`` if empty."""
    from calibration.omb.omb import OmBResult
    p = omb_cache_path(key_dir)
    if not p.exists():
        return None
    c = dict(np.load(p, allow_pickle=True))
    if c["time_cams"].size == 0:
        return None
    srcs = list(c["srcs"])
    mean_srcs = list(c["mean_srcs"]) if "mean_srcs" in c else srcs
    cams_beta = c["cams_beta"]
    prof, scalar, obs_interp, obs_mean, bias = {}, {}, {}, {}, {}
    for k in srcs:
        b = c[f"bias__{k}"].astype("float64")
        oi = c[f"obsint__{k}"].astype("float64")
        n_cams = b.shape[1]
        # period filter: keep altitudes valid in >= 25% of the CAMS steps (matches compute_omb
        # valid_frac_min default). valid_frac is returned for the figure's second x-axis.
        valid_frac = (np.sum(np.isfinite(b), axis=1) / n_cams
                      if n_cams else np.zeros(b.shape[0]))
        keep = valid_frac >= 0.25
        b_filt = np.where(keep[:, None], b, np.nan)
        oi_filt = np.where(keep[:, None], oi, np.nan)
        with np.errstate(invalid="ignore"):
            rms_profile = np.sqrt(np.nanmean(b_filt**2, axis=1))
            prof[k] = dict(obs_med=np.nanmedian(oi_filt, axis=1),
                           obs_p25=np.nanpercentile(oi_filt, 25, axis=1),
                           obs_p75=np.nanpercentile(oi_filt, 75, axis=1),
                           bias_med=np.nanmedian(b_filt, axis=1),
                           bias_p25=np.nanpercentile(b_filt, 25, axis=1),
                           bias_p75=np.nanpercentile(b_filt, 75, axis=1),
                           valid_frac=valid_frac)
            scalar[k] = dict(mean_bias=float(np.nanmean(b_filt)),
                             median_bias=float(np.nanmedian(b_filt)),
                             rms=float(np.nanmean(rms_profile)),
                             n_obs=int(np.sum(np.isfinite(b_filt))))
        obs_interp[k] = oi
        bias[k] = b
    for k in mean_srcs:
        obs_mean[k] = c[f"obsmean__{k}"].astype("float64")
    obs_full = {k: c[f"obsfull__{k}"].astype("float64")
                for k in list(c.get("full_srcs", [])) if f"obsfull__{k}" in c}
    cloud_base = (c["cloud_base"].astype("float64") if "cloud_base" in c
                  else np.full(c["time_cams"].size, np.nan))
    with np.errstate(invalid="ignore"):
        cams_med = np.nanmedian(cams_beta, axis=1)
    return OmBResult(wavelength=float(c["wavelength"]) if "wavelength" in c else float("nan"),
                     range_mean=c["range_mean"], z_cams=c["z_cams"], time_cams=c["time_cams"],
                     cams_beta=cams_beta, cams_med=cams_med, obs_mean=obs_mean,
                     obs_interp=obs_interp, bias=bias, prof=prof, scalar=scalar,
                     obs_full=obs_full, cloud_base=cloud_base)
