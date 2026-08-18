import sys, json, os
from pathlib import Path
import numpy as np
REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO)); sys.path.insert(0, str(REPO / "rayleigh_availability"))
from _adv_within_night import job_band
B = Path("C:/DATA/Projects/202606_E-PROFILE_calibration/rayleigh_availability")
OUT = Path("C:/Users/hervo/AppData/Local/Temp/claude/C--Users-hervo-OneDrive-Documents-ALC-rayleigh-calibration/265f36ab-1a40-4592-bae6-6c3550f2fac1/scratchpad/adv/bands2.json")

LABELS = ["PAYERNE_CHM15k_A", "LINDENBERG_CHM15k_0", "SAINT-CHRISTOPHE_A_CHM15k_0",
          "PALAISEAU_CHM15k_B", "SAINT-CHRISTOPHE_A_CL61_B", "CAMBORNE_CL61_C",
          "LANZHOT_CL61_A", "HOHENPEISSENBERG_CHM15k_0"]
NPER = int(sys.argv[1]) if len(sys.argv) > 1 else 30

def kept(lab):
    d = json.loads((B / "baselines" / ("base_eprof_v2_%s.json" % lab)).read_text(encoding="utf-8"))
    return sorted(k for k, v in d.items() if v[0] in (1.0, 0.5))

if __name__ == "__main__":
    from multiprocessing import Pool
    tasks = []
    rng = np.random.default_rng(7)
    for lab in LABELS:
        ds = kept(lab)
        if len(ds) > NPER:
            ds = [ds[i] for i in sorted(rng.choice(len(ds), NPER, replace=False))]
        tasks += [(lab, d, "eprof_v2", None) for d in ds]
    print("tasks:", len(tasks), flush=True)
    with Pool(8) as p:
        rows = []
        for i, r in enumerate(p.imap_unordered(job_band, tasks, chunksize=2)):
            rows.append(r)
            if (i + 1) % 25 == 0: print("  %d/%d" % (i + 1, len(tasks)), flush=True)
    OUT.write_text(json.dumps(rows))
    print("wrote", OUT, len(rows))
