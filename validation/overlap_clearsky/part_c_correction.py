"""Load all temperature-overlap models (D:/TEMP_MODELS/202606) and compute the correction applied at
25 deg C:  corr%(z) = a(z)*25 + b(z).  Stack -> data for the binscatter subplot."""
import netCDF4, numpy as np, glob, os
S = os.path.dirname(__file__)
fs = sorted(glob.glob("D:/TEMP_MODELS/202606/Overlap_correction_model_*.nc"))
rng = None; C = []; names = []
for f in fs:
    nc = netCDF4.Dataset(f); r = np.asarray(nc['range'][:], float); a = np.asarray(nc['a'][:], float); b = np.asarray(nc['b'][:], float); nc.close()
    if rng is None: rng = r
    c = a * 25.0 + b
    C.append(np.interp(rng, r, c)); names.append(os.path.basename(f).split('model_')[1].replace('.nc', ''))
C = np.array(C)                                              # (nstation, nrange)
np.savez(S + "/correction25.npz", rng=rng, C=C.astype('f4'), names=np.array(names))
print(f"{len(C)} temperature models loaded")
# distribution over range
for h in [150, 200, 300, 400, 500, 600, 800, 1000]:
    col = C[:, np.argmin(np.abs(rng - h))]
    print(f"  z={h:4d}m: median={np.nanmedian(col):+6.1f}%  p25={np.nanpercentile(col,25):+6.1f}  p75={np.nanpercentile(col,75):+6.1f}  min={np.nanmin(col):+7.1f}  max={np.nanmax(col):+7.1f}")
# ranking metric: mean |corr| over 200-500 m
metric = np.nanmean(np.abs(C[:, (rng >= 200) & (rng <= 500)]), axis=1)
print(f"\nlargest correction: {names[np.argmax(metric)]} (mean|corr|200-500m={metric.max():.1f}%)")
print(f"smallest correction: {names[np.argmin(metric)]} (mean|corr|200-500m={metric.min():.1f}%)")
