# -*- coding: utf-8 -*-
"""For each CORRECTABLE station (network_offset scan: Uccle CL51, Diepenbeek CL51, Akrotiri CL31),
show the effect of subtracting the fixed digitizer-ripple correction b_phys from L1 rcs_0:
  * pcolor (time-height) of the uncorrected and corrected attenuated backscatter, plus the removed
    ripple pattern (uncorrected - corrected);
  * time-averaged profiles (full range + a free-troposphere zoom) where the coherent 40 m ripple is
    visible in the uncorrected mean and flattened in the corrected mean.
A clear day is auto-selected per station (max cloud-free fraction) so the weak free-tropospheric signal
- where the ripple is a large fractional error - is visible."""
import sys, warnings
warnings.filterwarnings("ignore")
sys.path.insert(0, r"C:/Users/hervo/OneDrive/Documents/ALC_rayleigh_calibration")
try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import matplotlib.dates as mdates
from validation.paper._offset_lib import read_file, list_files, L1ROOT

D = Path("C:/DATA/Projects/202606_E-PROFILE_calibration/figs_paper_validation/paper_python")
NPZ = D / "network_offset"
REPORT_FIG = Path("C:/Users/hervo/OneDrive/Documents/ALC_rayleigh_calibration/doc/reports/figs_paper_report")
BETA_SCALE = 1e-2                       # rcs_0 (V*m^2) -> attenuated backscatter [Mm^-1 sr^-1], nominal C=1e8
ZTOP = 6000.0                          # pcolor top
STATIONS = [("CL51", "0-20000-0-06447", "A", "Uccle"),
            ("CL51", "0-20000-0-06477", "A", "Diepenbeek"),
            ("CL31", "0-20000-0-17601", "A", "Akrotiri")]
DAY_OVERRIDE = {"Uccle": ("2026", "05", "20260526"), "Diepenbeek": ("2026", "05", "20260526"),
                "Akrotiri": ("2026", "05", "20260527")}    # pre-found clearest days (skip the scan)


def load_day(wmo, ident, ymd):
    yr, mm, ds = ymd
    f = L1ROOT / wmo / yr / mm / f"L1_{wmo}_{ident}{ds}.nc"
    if not f.exists():
        return None
    tt, r, x, cbh, _ = read_file(f)
    return tt, r, x, cbh


def pick_clear_day(wmo, ident, months=("05", "04", "06")):
    """Return (datetimes, range, rcs_0, cbh) for the day with the most cloud-free profiles."""
    best = None
    for f in list_files(wmo, ident, months):
        try:
            tt, r, x, cbh, _ = read_file(f)
        except Exception:
            continue
        if x.ndim != 2 or tt.size < 200:
            continue
        nocloud = np.nanmax(cbh, axis=1) < 0 if cbh is not None else np.ones(tt.size, bool)
        score = int(nocloud.sum())
        if best is None or score > best[0]:
            best = (score, tt, r, x, cbh)
    return best[1:] if best else None


for itype, wmo, ident, site in STATIONS:
    z = np.load(NPZ / f"{wmo}_{ident}_{itype}.npz", allow_pickle=True)
    rng_b, b_phys, Lam = z["rng"], z["b_phys"], float(z["Lam_dom"])
    got = load_day(wmo, ident, DAY_OVERRIDE[site]) if site in DAY_OVERRIDE else None
    if got is None:
        got = pick_clear_day(wmo, ident)
    if got is None:
        print(f"{site}: no data"); continue
    tt, rng, X, cbh = got
    dr = float(np.median(np.diff(rng)))
    day = tt[len(tt) // 2].strftime("%Y-%m-%d")
    corr = np.interp(rng, rng_b, b_phys, left=0.0, right=0.0)      # rcs_0 units
    Xc = X - corr[None, :]
    beta_u = X * BETA_SCALE
    beta_c = Xc * BETA_SCALE
    # profiles: average over the cloud-free profiles of the day (coherent ripple survives, turbulence averages down)
    nocloud = np.nanmax(cbh, axis=1) < 0 if cbh is not None else np.ones(tt.size, bool)
    prof_u = np.nanmedian(beta_u[nocloud], axis=0)
    prof_c = np.nanmedian(beta_c[nocloud], axis=0)
    zk = rng / 1000.0
    print(f"{site} {itype}: day {day}, {tt.size} profiles ({int(nocloud.sum())} clear), Λ={Lam:.0f}m", flush=True)

    fig = plt.figure(figsize=(17, 9.2))
    gs = fig.add_gridspec(2, 3, height_ratios=[1.0, 1.0], hspace=0.30, wspace=0.30)
    zm = rng <= ZTOP
    tnum = mdates.date2num(tt)
    vmin = max(np.nanpercentile(beta_u[:, zm][beta_u[:, zm] > 0], 5), 1e-3)
    vmax = np.nanpercentile(beta_u[:, zm], 99.5)
    ext = [tnum[0], tnum[-1], rng[zm].min() / 1000, rng[zm].max() / 1000]

    def draw_pcolor(ax, beta, title):
        im = ax.imshow(np.clip(beta[:, zm].T, vmin, vmax), aspect="auto", origin="lower", extent=ext,
                       cmap="turbo", norm=LogNorm(vmin=vmin, vmax=vmax))
        ax.xaxis_date(); ax.xaxis.set_major_formatter(mdates.DateFormatter("%H"))
        ax.set_ylabel("height [km]"); ax.set_title(title, fontsize=10.5)
        ax.axhline(1.8, color="w", lw=0.8, ls=":")            # correction valid above here
        return im

    # --- row 1: full context (0-6 km, log) ---
    ax0 = fig.add_subplot(gs[0, 0]); im = draw_pcolor(ax0, beta_u, f"(a) uncorrected  β_att   ({site} {itype}, {day})")
    ax1 = fig.add_subplot(gs[0, 1]); draw_pcolor(ax1, beta_c, "(b) corrected  β_att")
    fig.colorbar(im, ax=[ax0, ax1], label=r"β_att [Mm$^{-1}$sr$^{-1}$]", pad=0.01, fraction=0.03)
    # (c) day-mean profile, full range (log)
    axc = fig.add_subplot(gs[0, 2])
    axc.plot(prof_u, zk, color="0.55", lw=1.5, label="uncorrected")
    axc.plot(prof_c, zk, color="#d62728", lw=1.0, label="corrected")
    axc.set_xscale("log"); axc.set_ylim(0, ZTOP / 1000); axc.set_xlabel(r"β_att [Mm$^{-1}$sr$^{-1}$]")
    axc.set_ylabel("height [km]"); axc.set_title("(c) day-mean profile (clear)"); axc.legend(fontsize=9); axc.grid(alpha=0.3)

    # --- row 2: free-troposphere zoom where the ripple is resolvable ---
    zlo, zhi = 2.5, 4.5
    zmz = (rng >= zlo * 1000) & (rng <= zhi * 1000)
    extz = [tnum[0], tnum[-1], zlo, zhi]
    bz = beta_u[:, zmz].T; bzc = beta_c[:, zmz].T
    vlo, vhi = np.nanpercentile(bz, 3), np.nanpercentile(bz, 97)

    def zoom_pcolor(ax, b, title):
        im = ax.imshow(b, aspect="auto", origin="lower", extent=extz, cmap="turbo", vmin=vlo, vmax=vhi)
        ax.xaxis_date(); ax.xaxis.set_major_formatter(mdates.DateFormatter("%H"))
        ax.set_ylabel("height [km]"); ax.set_xlabel("hour [UTC]"); ax.set_title(title, fontsize=10.5)
        return im

    ax3 = fig.add_subplot(gs[1, 0]); imz = zoom_pcolor(ax3, bz, f"(d) uncorrected zoom — {Lam:.0f} m banding")
    ax4 = fig.add_subplot(gs[1, 1]); zoom_pcolor(ax4, bzc, "(e) corrected zoom — banding removed")
    fig.colorbar(imz, ax=[ax3, ax4], label=r"β_att [Mm$^{-1}$sr$^{-1}$]", pad=0.01, fraction=0.03)
    # (f) zoom mean profile (linear) — the ripple removed
    axf = fig.add_subplot(gs[1, 2])
    zz = (rng >= zlo * 1000) & (rng <= zhi * 1000)
    axf.plot(prof_u[zz], zk[zz], color="0.55", lw=1.7, label="uncorrected")
    axf.plot(prof_c[zz], zk[zz], color="#d62728", lw=1.7, label="corrected")
    axf.set_ylim(zlo, zhi); axf.set_xlabel(r"β_att [Mm$^{-1}$sr$^{-1}$]"); axf.set_ylabel("height [km]")
    axf.set_title("(f) zoom day-mean — ripple removed"); axf.legend(fontsize=9); axf.grid(alpha=0.3)

    fig.suptitle(f"{site} {itype} — digitizer-ripple correction on the aerosol backscatter "
                 f"(subtract b_phys from L1 rcs₀; Λ={Lam:.0f} m)", fontweight="bold", fontsize=13)
    out = D / f"fig_corrected_signal_{site.lower()}_{itype}.png"
    fig.savefig(out, dpi=145, bbox_inches="tight")
    fig.savefig(REPORT_FIG / out.name, dpi=145, bbox_inches="tight")
    print("  saved", out.name, flush=True)
print("CORRECTABLE_SIGNAL_FIGURES_DONE")
