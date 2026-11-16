"""Droplet effective-radius climatology from the Cloudnet `der` product for Payerne,
Palaiseau and Lindenberg (A:/Cloudnet_Droplet_radius, 2025-). Histograms of the measured
droplet effective radius, for (a) ALL valid liquid-cloud pixels and (b) the CLOUD-BASE region
(the MS-relevant part the ceilometer beam meets first), compared to:
  - our fitted a_G = 8 um (equivalent-area RADIUS = 16 um diameter),
  - Le & O'Connor (2026) band: 8-20 um effective DIAMETER = 4-10 um effective radius,
  - literature continental/marine stratocumulus effective radius (~5-9 / ~8-13 um).

I/O-lean + threaded (the A: drive is the bottleneck; netCDF releases the GIL on read).
"""
from __future__ import annotations

import sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

CLOUDNET = Path("D:/Cloudnet_Droplet_radius")   # moved from A: (2026-10)
OUT = Path("C:/DATA/Projects/202606_E-PROFILE_calibration/old/ms_fov")   # moved from .../ms_fov (2026-10)
REPO_FIGS = Path(__file__).resolve().parents[1] / "doc" / "reports" / "figs_multiple_scattering"
SITES = ["payerne", "palaiseau", "lindenberg"]
WORKERS = 30


def read_one(fp):
    """Return (all_re_um, base_re_um) for one day's der file, or None."""
    from netCDF4 import Dataset
    try:
        with Dataset(fp) as nc:
            der = np.asarray(nc.variables["der"][:], float)   # (time,height) m, effective RADIUS
            h = np.asarray(nc.variables["height"][:], float)
    except Exception:
        return None
    der[~np.isfinite(der) | (der <= 0) | (der > 40e-6)] = np.nan
    allre = der[np.isfinite(der)] * 1e6
    # cloud-base re: lowest valid gate per profile (time row)
    base = []
    for ti in range(der.shape[0]):
        row = der[ti]
        idx = np.where(np.isfinite(row))[0]
        if idx.size:
            base.append(row[idx[0]] * 1e6)
    return allre, np.array(base)


# Cloudnet re at the VALID-CALIBRATION scenes: the match cache holds the re matched to the
# in-cloud profiles that a cloud calibration actually used (median over cloud base .. +300 m).
# We take the Vaisala instruments (CL61/CL31 — the ones we cloud-calibrate; CHM15k saturates).
CAL_VAIS = {"payerne": ["Payerne__CL61__re", "Payerne__CL31__re"],
            "palaiseau": ["Palaiseau__CL31__re"],
            "lindenberg": ["Lindenberg__CL61__re"]}


def _calib_scene_re(site):
    fp = OUT / "match_cache.npz"
    if not fp.exists():
        return np.array([])
    z = np.load(fp, allow_pickle=False)
    vals = []
    for k in CAL_VAIS.get(site, []):
        if k in z.files:
            vals.append(z[k])
    x = np.concatenate(vals) if vals else np.array([])
    return x[np.isfinite(x) & (x > 0) & (x < 40)]


def main():
    REPO_FIGS.mkdir(parents=True, exist_ok=True)
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    per = {s: dict(all=[], base=[], cal=_calib_scene_re(s)) for s in SITES}
    for s in SITES:
        files = sorted(CLOUDNET.glob(f"*_{s}_der.nc"))
        print(f"{s}: {len(files)} daily files", flush=True)
        # ProcessPool (not threads): HDF5/netCDF is NOT thread-safe -> concurrent Dataset opens
        # segfault. This module is light (numpy/netCDF4 only, no calibration stack), so process
        # spawn is cheap and each process has its own HDF5 state.
        with ProcessPoolExecutor(max_workers=WORKERS) as ex:
            for res in ex.map(read_one, files, chunksize=4):
                if res is None:
                    continue
                a, b = res
                if a.size:
                    per[s]["all"].append(a)
                if b.size:
                    per[s]["base"].append(b)
        per[s]["all"] = np.concatenate(per[s]["all"]) if per[s]["all"] else np.array([])
        per[s]["base"] = np.concatenate(per[s]["base"]) if per[s]["base"] else np.array([])

    # ---- stats printout ----
    print("\nCloudnet droplet EFFECTIVE RADIUS climatology (um):")
    print(f"{'site':11s} {'sample':6s} {'n':>9s} {'median':>7s} {'mean':>6s} {'P10':>6s} {'P90':>6s}  "
          f"{'=> median eff.DIAM':>18s}")
    stats = {}
    for s in SITES:
        for lab in ("all", "base", "cal"):
            x = per[s][lab]; x = x[np.isfinite(x)]
            if x.size < 10:
                continue
            md = np.median(x)
            print(f"{s:11s} {lab:6s} {x.size:9d} {md:7.1f} {np.mean(x):6.1f} "
                  f"{np.percentile(x,10):6.1f} {np.percentile(x,90):6.1f}  {2*md:16.0f} um")
            stats[(s, lab)] = md

    # ---- figure: histograms per site (radius axis; diameter on top) ----
    bins = np.arange(0, 22, 0.75)
    fig, axes = plt.subplots(1, 3, figsize=(16, 5.4), sharey=True)
    disp = {"payerne": "Payerne", "palaiseau": "Palaiseau", "lindenberg": "Lindenberg"}
    labmap = {"all": "all cloud pixels", "base": "cloud base", "cal": "VALID-CALIBRATION scenes"}
    for ax, s in zip(axes, SITES):
        for lab, c, al, ht in (("all", "C0", 0.35, "bar"), ("base", "C3", 0.45, "bar"),
                               ("cal", "gold", 1.0, "step")):
            x = per[s][lab]; x = x[np.isfinite(x) & (x < 22)]
            if x.size:
                ax.hist(x, bins=bins, density=True, color=c, alpha=al, histtype=ht,
                        lw=(2.4 if ht == "step" else 0),
                        label=f"{labmap[lab]} (med {np.median(x):.1f} $\\mu$m, n={x.size})")
        # reference markers (RADIUS). Omit the fitted-value line in --paper mode.
        if "--paper" not in sys.argv:
            ax.axvline(8.0, color="k", ls="--", lw=1.8, label="our fitted a$_G$=8 $\\mu$m (=16 $\\mu$m D)")
        ax.axvspan(4.0, 10.0, color="green", alpha=0.10,
                   label="Le et al. 2026 band 8-20 $\\mu$m diam. (=4-10 $\\mu$m radius)")
        ax.axvspan(5.0, 9.0, color="orange", alpha=0.0)  # continental Sc note in legend only
        ax.set_title(disp[s]); ax.set_xlabel("droplet effective radius [$\\mu$m]")
        ax.set_xlim(0, 20); ax.grid(alpha=0.3); ax.legend(fontsize=7.5, loc="upper right")
        secax = ax.secondary_xaxis("top", functions=(lambda r: 2 * r, lambda d: d / 2))
        secax.set_xlabel("effective diameter [$\\mu$m]", fontsize=8)
    axes[0].set_ylabel("probability density")
    fig.suptitle("Cloudnet droplet effective-radius climatology (2025-2026) vs the MS-correction "
                 "droplet-size assumptions\ncontinental stratocumulus literature: eff. radius "
                 "~5-9 $\\mu$m (Miles 2000; Frisch 2002); marine ~8-13 $\\mu$m", fontweight="bold")
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    name = "fig_droplet_climatology_paper.png" if "--paper" in sys.argv else "fig_droplet_climatology.png"
    for d in (OUT, REPO_FIGS):
        fig.savefig(d / name, dpi=150)
    plt.close(fig)
    print(f"\n  {name}")


if __name__ == "__main__":
    main()
