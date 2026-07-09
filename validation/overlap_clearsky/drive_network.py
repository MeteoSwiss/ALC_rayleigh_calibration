"""Quantify the applied overlap of every network CHM15k from clear-sky night noise (M2, K=2 anchored).
~12 clear nights per station (convergence: ~10 nights -> RMSE 0.02 vs hood)."""
import numpy as np, json, os, sys, warnings
warnings.filterwarnings("ignore")
import clearsky_overlap as CS
from overlap_fit_general import fit_K, model
from scipy.signal import medfilt

MONTHS = ["2026-04", "2026-05", "2026-06"]                   # contiguous recent 3-month window
COMMON = np.arange(15.0, 2501.0, 10.0)                       # common range grid (m) for mixed 1024/1536 CHM grids

def despike(rng, O2):
    """Continuity: interpolate over unphysical spikes (beta->0), then smooth."""
    bad = ~np.isfinite(O2) | (O2 > 1.15) | (O2 < -0.05)
    Od = np.interp(rng, rng[~bad], O2[~bad]) if (~bad).sum() > 20 else O2.copy()
    return medfilt(np.clip(Od, 0, 1.15), 5)

def anchored_k2(rng, O2):
    p, _ = fit_K(rng, despike(rng, O2), 2, 150, 1850)        # despike BEFORE the fit (continuity/robustness)
    f = model(rng / 1000.0, p); f0 = model(np.array([0.0]), p)[0]; finf = model(np.array([5.0]), p)[0]
    if f0 > 0 and finf > f0:                                  # guard: only anchor when floor>0
        out = (f - f0) / (finf - f0)
    else:
        out = f
    return np.clip(out, 0.0, 1.0).astype('f4')               # overlap is physically in [0, 1]

def process_station(wmo):
    try:
        files = CS.select_clear_nights(wmo, MONTHS, max_nights=25, min_clear=0.80)
        if len(files) < 6:
            return (wmo, None, len(files), 1.0)
        R = CS.derive(files, verbose=False); rng = R['rng']; O2 = R['O2']
        spike = float(np.mean(O2[(rng >= 150) & (rng <= 2000)] > 1.15))   # unphysical-spike fraction (retrieval quality)
        Oa = anchored_k2(rng, O2)
        Oc = np.interp(COMMON, rng, Oa).astype('f4')         # onto common grid (handles 1024/1535/1536)
        return (wmo, Oc, int(R['ndays']), spike)
    except Exception:
        return (wmo, None, -1, 1.0)

if __name__ == "__main__":
    from concurrent.futures import ProcessPoolExecutor, as_completed
    wmos = json.load(open(os.path.dirname(__file__) + "/chm15k_stations.json"))
    rng = COMMON
    results = {}; nnights = {}; spikes = {}
    print(f"{len(wmos)} CHM15k stations, 10 workers", flush=True)
    with ProcessPoolExecutor(max_workers=10) as ex:
        futs = {ex.submit(process_station, w): w for w in wmos}
        done = 0
        for fu in as_completed(futs):
            w, O, n, sp = fu.result(); done += 1
            if O is not None:
                results[w] = O; nnights[w] = n; spikes[w] = sp
            if done % 10 == 0 or O is None:
                print(f"[{done}/{len(wmos)}] {w}: {'OK '+str(n)+'n spk'+format(sp,'.2f') if O is not None else 'skip('+str(n)+')'}", flush=True)
    curves = np.array([results[w] for w in results])
    np.savez(os.path.dirname(__file__) + "/network_overlaps.npz",
             rng=rng, curves=curves, wmos=np.array(list(results.keys())),
             nnights=np.array([nnights[w] for w in results]),
             spikes=np.array([spikes[w] for w in results]))
    print(f"\n{len(results)}/{len(wmos)} stations retrieved -> network_overlaps.npz", flush=True)
