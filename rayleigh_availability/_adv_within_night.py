"""ADVERSARIAL CHECK: WITHIN-night altitude dependence of C_L, read off cl_profile(z).

cl_profile is the retrieval's own C_L evaluated at every range gate, with the aerosol
two-way transmission correction applied. If the constant were altitude independent, its
slope inside (and just around) the chosen molecular window would be zero.
"""
import sys, json, logging, warnings, os
from pathlib import Path
import numpy as np
REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))
import calibration.rayleigh.calibration as CAL
from calibration import calibrate_rayleigh, CalibrationOptions, InstrumentInfo, DataLevel
from calibration.config import InstrumentType

L1_ROOT = Path("D:/E-PROFILE_L1_2026")
MANIFEST = json.loads((REPO / "rayleigh_availability" / "scope_availability.json").read_text())
SCRATCH = Path("C:/Users/hervo/AppData/Local/Temp/claude/C--Users-hervo-OneDrive-Documents-ALC-rayleigh-calibration/265f36ab-1a40-4592-bae6-6c3550f2fac1/scratchpad/adv")

_STASH = {}
_ORIG = CAL._compute_cl_for_perturbation
def _patched(*a, **kw):
    r = _ORIG(*a, **kw)
    if kw.get("return_diagnostics"):
        _STASH["cl_profile"] = r.cl_profile
    return r
CAL._compute_cl_for_perturbation = _patched

def opts(tag, method, params):
    o = CalibrationOptions.from_json(REPO / "options.json")
    o.folder_root = L1_ROOT; o.data_level = DataLevel.L1
    o.cams_folder = Path("D:/CAMS_run_v20")
    o.molecular_method = method
    if params: o.molecular_params = params
    o.plot_main = o.plot_all = False
    o.folder_output = SCRATCH / ("out_" + tag)
    o.folder_output.mkdir(parents=True, exist_ok=True)
    return o

def info_for(label):
    for i in MANIFEST:
        if i["label"] == label:
            return InstrumentInfo(site_name=i["label"], wmo_id=i["wmo"], identifier=i["ident"],
                                  instrument_type=InstrumentType(i["type"]),
                                  latitude=i["lat"], longitude=i["lon"], altitude=i["alt"]), i
    raise SystemExit("no " + label)

def job(a):
    label, ds, method, params = a
    warnings.filterwarnings("ignore"); logging.getLogger().setLevel(logging.ERROR)
    inf, meta = info_for(label)
    tag = "%s_%d" % (label.replace('/','_')[:10], os.getpid())
    row = dict(label=label, date=ds)
    _STASH.clear(); fit = {}
    try:
        r = calibrate_rayleigh(ds, inf, opts(tag, method, params), fit_inputs_out=fit)
    except Exception as e:
        row["err"] = f"{type(e).__name__}: {e}"[:100]; return row
    row["flag"] = float(r.flag)
    cp = _STASH.get("cl_profile")
    if cp is None or r.flag not in (1.0, 0.5): return row
    rng = np.asarray(fit["range_alc"], float)
    cp = np.asarray(cp, float)
    row.update(C=float(r.lidar_constant), alt=float(meta["alt"]),
               zb=float(r.calibration_bottom_height), zt=float(r.calibration_top_height),
               rng=rng.tolist(), cl=cp.tolist())
    return row

def slopes(row):
    """within-night d ln C_L / dz in %/km, inside the window and over the 2-6 km search band."""
    if "cl" not in row: return None
    rng = np.asarray(row["rng"]); cl = np.asarray(row["cl"]); alt = row["alt"]
    zasl = rng + alt
    out = {}
    for name, lo, hi in (("win", row["zb"], row["zt"]), ("band", alt+2000., alt+6000.)):
        m = (zasl >= lo) & (zasl <= hi) & np.isfinite(cl) & (cl > 0)
        if m.sum() < 15: out[name] = (np.nan, np.nan, m.sum()); continue
        x = zasl[m]/1000.; y = np.log(cl[m])
        A = np.column_stack([np.ones_like(x), x])
        b, res, *_ = np.linalg.lstsq(A, y, rcond=None)
        r = y - A@b; s2 = r@r/(len(x)-2)
        cov = s2*np.linalg.pinv(A.T@A)
        out[name] = (100*b[1], 100*np.sqrt(cov[1,1]), int(m.sum()))
    return out


# ---------------------------------------------------------------------------
# FORCED-WINDOW within-night experiment: same night, low vs high fit window.
# ---------------------------------------------------------------------------
def job_band(a):
    """Run one night twice with the molecular-window SEARCH BAND forced low / high."""
    label, ds, method, params = a
    warnings.filterwarnings("ignore"); logging.getLogger().setLevel(logging.ERROR)
    inf, meta = info_for(label)
    tag = "%s_%d" % (label.replace('/', '_')[:10], os.getpid())
    row = dict(label=label, date=ds)
    for name, lo, hi in (("low", 2000., 3400.), ("high", 4600., 6000.)):
        o = opts(tag + "_" + name, method, params)
        o.range_start_m = lo; o.range_end_m = hi; o.min_window_start_m = lo
        try:
            r = calibrate_rayleigh(ds, inf, o)
        except Exception as e:
            row[name + "_err"] = f"{type(e).__name__}: {e}"[:80]; continue
        row[name + "_flag"] = float(r.flag)
        if r.flag in (1.0, 0.5):
            row[name + "_C"] = float(r.lidar_constant)
            row[name + "_zb"] = float(r.calibration_bottom_height)
            row[name + "_zt"] = float(r.calibration_top_height)
            row[name + "_u"] = float(r.uncertainty)
    return row
