"""PAPER figure (2 subplots, English):
 (1) applied overlap of every network CHM15k (clear-sky noise, M2 K=2 anchored): binscatter density
     + median + 25-75% band + the latest- and earliest-completing instrument (quality-controlled).
 (2) temperature-model correction applied at 25 deg C (a*25+b) from D:/TEMP_MODELS/202606: binscatter
     + median + 25-75% + the largest and smallest correction."""
import numpy as np, os
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from scipy.signal import medfilt
S = os.path.dirname(__file__)
sm = lambda y, k=5: medfilt(np.nan_to_num(y, nan=0), k)

def binscatter(ax, curves, rng, vlo, vhi, nvbins=60, cmap='magma', vmin=0.006, vmax=1.0):
    vbins = np.linspace(vlo, vhi, nvbins + 1); vc = 0.5 * (vbins[:-1] + vbins[1:])
    H = np.full((len(rng), nvbins), np.nan)
    for i in range(len(rng)):
        col = curves[:, i]; col = col[np.isfinite(col)]
        if len(col) < 5: continue
        h, _ = np.histogram(col, bins=vbins); s = h.sum()
        if s > 0:
            d = h / s; d[d == 0] = np.nan; H[i, :] = d
    return ax.pcolormesh(vc, rng, np.ma.masked_invalid(H), cmap=cmap, shading='auto',
                         norm=LogNorm(vmin=vmin, vmax=vmax))

def first_cross(O, thr, rng):
    """First range where O reaches thr (robust to plateaus/non-monotonicity); NaN if never."""
    idx = np.argmax(O >= thr)
    return rng[idx] if O[idx] >= thr else np.nan

def quality_mask(curves, rng):
    """Keep physically-valid retrievals: completes at far range (not blown up), low near-ground floor,
    a real mid-range rise, a PHYSICAL maximum gradient (rejects near-vertical artefacts), sensible z50,
    reaches 0.9, and near-monotonic."""
    far = np.nanmedian(curves[:, (rng >= 1200) & (rng <= 1600)], axis=1)
    near = np.nanmedian(curves[:, (rng >= 30) & (rng <= 90)], axis=1)
    rise = (np.nanmedian(curves[:, (rng >= 400) & (rng <= 600)], axis=1)
            - np.nanmedian(curves[:, (rng >= 80) & (rng <= 150)], axis=1))
    band = (rng >= 80) & (rng <= 1500)
    mg = np.array([np.nanmax(np.gradient(medfilt(np.clip(c, 0, 1.2), 5), rng)[band]) * 1000 for c in curves])
    z50 = np.array([first_cross(np.clip(c, 0, 1), 0.5, rng) for c in curves])   # first-crossing (robust)
    z90 = np.array([first_cross(np.clip(c, 0, 1), 0.9, rng) for c in curves])
    nonmono = np.array([np.mean(np.diff(c[(rng >= 100) & (rng <= 1400)]) < -0.005) for c in curves])
    ok = ((far > 0.85) & (far < 1.05) & (near < 0.30) & (rise > 0.20) & (mg < 18.0)
          & (z50 > 150) & (z50 < 1100) & np.isfinite(z90) & (z90 < 1600) & (nonmono < 0.08)
          & np.isfinite(curves).all(axis=1))
    return ok, z90, mg

fig, ax = plt.subplots(1, 2, figsize=(17, 8))

# ---------- subplot 1: network overlaps ----------
net = np.load(S + "/network_overlaps.npz", allow_pickle=True)
rng = net['rng']; curves_all = np.clip(net['curves'], 0.0, 1.0); wmos_all = net['wmos']  # overlap is physically <= 1
spk_all = net['spikes'] if 'spikes' in net.files else np.zeros(len(curves_all))
qm, z90a, mga = quality_mask(curves_all, rng)
qm = qm & (spk_all < 0.25)                                    # reject badly spike-corrupted retrievals (beta->0)
curves = curves_all[qm]; wmos = wmos_all[qm]; z90 = z90a[qm]; nrej = int((~qm).sum())
a = ax[0]
im = binscatter(a, curves, rng, 0.0, 1.02, cmap='magma')
med = np.array([np.nanmedian(curves[:, i]) for i in range(curves.shape[1])])
p25 = np.array([np.nanpercentile(curves[:, i], 25) for i in range(curves.shape[1])])
p75 = np.array([np.nanpercentile(curves[:, i], 75) for i in range(curves.shape[1])])
a.plot(sm(med), rng, 'k-', lw=3.2); a.plot(sm(med), rng, color='cyan', lw=1.8, label='median')
a.plot(sm(p25), rng, 'w--', lw=1.5); a.plot(sm(p75), rng, 'w--', lw=1.5, label='25–75 %')
khi = int(np.argmax(z90)); klo = int(np.argmin(z90))
a.plot(sm(curves[khi]), rng, color='lime', lw=2, label=f'latest to complete ({wmos[khi].split("-")[-1]}, z₉₀={z90[khi]:.0f} m)')
a.plot(sm(curves[klo]), rng, color='deepskyblue', lw=2, label=f'earliest to complete ({wmos[klo].split("-")[-1]}, z₉₀={z90[klo]:.0f} m)')
a.set_xlabel("overlap  O(z)"); a.set_ylabel("Range AGL [m]"); a.set_xlim(0, 1.02); a.set_ylim(0, 1600)
a.set_title(f"(1) Applied overlap — {len(curves)} network CHM15k ({nrej} rejected by QC)\n(clear-sky night noise, M2 anchored)")
a.legend(fontsize=8.5, loc='lower right', framealpha=0.85)
cb = fig.colorbar(im, ax=a, pad=0.01); cb.set_label("station density (per altitude)", fontsize=8)

# ---------- subplot 2: 25C correction ----------
cor = np.load(S + "/correction25.npz", allow_pickle=True)
rc = cor['rng']; C = cor['C']; names = cor['names']
b = ax[1]
imc = binscatter(b, C, rc, -50, 100, cmap='magma')
medc = np.array([np.nanmedian(C[:, i]) for i in range(C.shape[1])])
p25c = np.array([np.nanpercentile(C[:, i], 25) for i in range(C.shape[1])])
p75c = np.array([np.nanpercentile(C[:, i], 75) for i in range(C.shape[1])])
b.plot(sm(medc), rc, 'k-', lw=3.2); b.plot(sm(medc), rc, color='cyan', lw=1.8, label='median')
b.plot(sm(p25c), rc, 'w--', lw=1.5); b.plot(sm(p75c), rc, 'w--', lw=1.5, label='25–75 %')
metric = np.nanmean(np.abs(C[:, (rc >= 200) & (rc <= 500)]), axis=1)
kbig = int(np.nanargmax(metric)); ksml = int(np.nanargmin(metric))
b.plot(sm(C[kbig]), rc, color='lime', lw=2, label=f'largest correction ({names[kbig].split("_")[0]})')
b.plot(sm(C[ksml]), rc, color='deepskyblue', lw=2, label=f'smallest correction ({names[ksml].split("_")[0]})')
b.axvline(0, color='w', lw=0.8, alpha=0.6)
b.set_xlabel("correction applied at 25 °C  [%]"); b.set_ylabel("Range AGL [m]"); b.set_xlim(-50, 100); b.set_ylim(0, 1600)
b.set_title(f"(2) Temperature-model correction at 25 °C\n({len(C)} models, a·25+b ; ≈0 above 500 m)")
b.legend(fontsize=8.5, loc='upper right', framealpha=0.85)
cbc = fig.colorbar(imc, ax=b, pad=0.01); cbc.set_label("station density (per altitude)", fontsize=8)

plt.suptitle("Network CHM15k applied overlap (clear-sky noise) vs temperature-model correction at 25 °C  —  April–June 2026", fontsize=13)
plt.tight_layout(rect=[0, 0, 1, 0.96])
out = os.path.normpath(os.path.join(S, "..", "..", "doc", "reports", "figs_overlap_from_noise",
                                    "08_network_overlap_vs_correction.png"))            # canonical paper figure
os.makedirs(os.path.dirname(out), exist_ok=True)
plt.savefig(out, dpi=125, bbox_inches='tight'); print("saved", out)
print(f"overlap: {len(curves)} kept, {nrej} rejected | z90 median={np.median(z90):.0f}m, {np.min(z90):.0f}-{np.max(z90):.0f}m")
