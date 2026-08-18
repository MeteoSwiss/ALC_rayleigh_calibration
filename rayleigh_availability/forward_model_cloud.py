"""Forward model of a liquid-water cloud for the O'Connor/Hopkin cloud calibration.

Purpose
-------
Closure test of the SHIPPED cloud estimator (``calibration.cloud``): synthesise the
attenuated backscatter an ALC *would* measure for a KNOWN lidar constant ``C_true`` and a
KNOWN multiple-scattering factor ``eta_true(z)``, push it through the unmodified pipeline
(:func:`calibration.cloud.calibration.liquid_cloud_calibration_from_data`) and compare the
retrieved Wiegner constant with ``C_true``.  Nothing of the estimator is re-implemented
here: only the *measurement* is simulated.

The estimator being tested (read from calibration/cloud/calibration.py + _filters.py)
-------------------------------------------------------------------------------------
For every profile i, on the operational working grid (30 s x 10 m), with
``beta`` = physical attenuated backscatter [m^-1 sr^-1] = rcs / C_applied:

  1. 910 nm only:  beta /= trans2_wv(z, t)                      (water vapour, two-way)
  2. beta_corr(z) = eta_table(z) * beta(z)                      (multiple scattering)
     eta_table is a PER-INSTRUMENT-TYPE lookup in RANGE (km), linearly interpolated and
     linearly EXTRAPOLATED (``_interp1_linear_extrap``) from the PVC/Hogan-2006 tables
     _ETA_CL31/_ETA_CL51/_ETA_CL61/_ETA_CHM15K/... in calibration/cloud/_filters.py.
     It is applied POINTWISE in range, i.e. as a weight inside the integral below.
  3. S_apparent(i) = 1 / (2 * trapz_z(beta_corr(z, i) dz))      [sr]
     integrated on the FIXED gate window  cal_minheight = 100 m .. cal_maxheight = 2400 m
     AGL (config, not derived from the cloud base). Profiles with >10 % NaN or <3 valid
     gates in that window are dropped.
  4. cloud filters: (a) beta 300 m ABOVE the in-gate peak must be < peak/20, (b) same
     300 m BELOW, (c) below-cloud fraction  sum(beta[100 m .. peak-50 m]) /
     sum(beta[100..2400 m]) <= ratio_filter = 0.05, (d) 500 m <= CBH <= 2400 m.
  5. temporal consistency: >= 5 consecutive profiles within +-10 % of their mean.
  6. C(i) = S_consistent(i) / S_THEORETICAL,  S_THEORETICAL = 18.8 sr (O'Connor 2004).
  7. optional aerosol transmission correction (operational: ON, LR_aer = 50 sr):
     B = sum(beta[gate 4 .. peak-5]) * dr ;  C <- C * exp(-2 * LR_aer * C * B).
     NB the integrand is the FULL below-cloud attenuated backscatter, molecular included.
  8. C_L = calibration_constant_applied / median_i(C)           (Wiegner constant, headline)

So the whole estimator is "one number per profile = 1 / (2 * integral of eta*beta over a
fixed 100-2400 m window)", and every altitude dependence must come from (2), (3) or (7).

The forward model implemented here
----------------------------------
Range grid r [m AGL] (r == height AGL; the instrument is assumed vertical).  Units:
range m, backscatter m^-1 sr^-1, extinction m^-1, LWC g m^-3.

  beta_att(r) = [beta_mol(r) + beta_aer(r) + beta_cld(r)]
                * exp(-2 * INT_0^r [ alpha_mol + alpha_aer + eta_true * alpha_cld ] dr')

  * molecular: US Std 1976 (repo ``load_standard_atmosphere`` + Bucholtz-1995
    ``calculate_molecular_properties``), alpha_mol = beta_mol * 8*pi/3.
  * aerosol: exponential layer alpha_aer(r) = alpha0 exp(-r/H) below the cloud base,
    scaled to a prescribed below-cloud AOD; beta_aer = alpha_aer / LR_aer_true.
  * liquid cloud: adiabatic LWC(r) = Gamma_l * (r - CBH) with a fixed droplet number N,
    r_eff = (3 LWC / (4 pi rho_w N))^(1/3), alpha_cld = 3 LWC / (2 rho_w r_eff)
    (so alpha ~ (r-CBH)^(2/3), the adiabatic ramp), beta_cld = alpha_cld / S_cloud
    with S_cloud = 18.8 sr = the value the estimator divides by.
  * MULTIPLE SCATTERING IS PHYSICS HERE, not a correction: following Platt (1981) the
    forward-scattered photons reduce the APPARENT cloud extinction, so eta_true(r)
    multiplies alpha_cld inside the exponent (and only there -- molecular and aerosol
    single-scatter).  Analytically, with u = 2 INT eta alpha_cld, du = 2 eta alpha dr,
        INT eta_true beta_att dr = (1/(2 S)) INT e^-u du -> 1/(2 S)   for a fully
    attenuating cloud, EXACTLY, for any eta(r) profile -- which is why the pipeline's
    pointwise eta multiplication is the correct inverse of this physics and the closure
    test below must return C_true.

  The synthetic instrument reports rcs = C_true * beta_att, and the pipeline is fed
  beta_in = rcs / C_applied, so the retrieval must give C_L = C_true * (18.8 / S_cloud).

Usage
-----
    python rayleigh_availability/forward_model_cloud.py            # full experiment set
    python rayleigh_availability/forward_model_cloud.py --quick    # coarse CBH grid
"""

from __future__ import annotations

import argparse
import sys
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, Optional, Tuple

import numpy as np
from numpy.typing import NDArray

_REPO = Path(__file__).resolve().parent.parent
if str(_REPO) not in sys.path:
    sys.path.insert(0, str(_REPO))

from calibration.rayleigh.atmosphere import (           # noqa: E402
    load_standard_atmosphere, calculate_molecular_properties, MOLECULAR_LIDAR_RATIO)
from calibration.cloud import _filters as _CF           # noqa: E402
from calibration.cloud.calibration import (             # noqa: E402
    CloudCalConfig, CeiloData, set_defaults, liquid_cloud_calibration_from_data,
    INSTRUMENT_CAL_DEFAULT)

# --- physical constants (SI unless stated) ---------------------------------------------
RHO_W = 1.0e6            # g m^-3, density of liquid water
S_CLOUD_TRUE = 18.8      # sr, liquid-water lidar ratio assumed by O'Connor (2004)
LR_AER_TRUE = 50.0       # sr, aerosol lidar ratio used to build the synthetic aerosol
GAMMA_L = 2.0e-3         # g m^-3 per m, adiabatic LWC lapse (~2 g m^-3 km^-1, T ~ 5-10 C)
N_DROP = 2.0e8           # m^-3, cloud droplet number concentration (200 cm^-3)
WAVELENGTH_M = {"CL31": 909.7e-9, "CL51": 910.0e-9, "CL61": 910.74e-9,
                "CHM15k": 1064.47e-9, "CHM15K": 1064.47e-9}

# --- LEGACY (pre-2026-07) multiple-scattering tables, Hopkin et al. 2019 / a_G = 8 um ---
# Recovered from the pre-split module (calibration/cloud/calibration.py in the
# hardcore-mccarthy worktree). CL61 borrowed the CL51 table back then.
LEGACY_ETA_CL31 = np.array([
    [0.250, 0.82854], [0.375, 0.82371], [0.625, 0.81608], [0.875, 0.80811],
    [1.125, 0.79969], [1.375, 0.79027], [1.625, 0.78227], [1.875, 0.77480],
    [2.125, 0.76710], [2.375, 0.76088]])
LEGACY_ETA_CL51 = np.array([
    [0.250, 0.82881], [0.375, 0.82445], [0.625, 0.81752], [0.875, 0.81021],
    [1.125, 0.80241], [1.375, 0.79356], [1.625, 0.78595], [1.875, 0.77877],
    [2.125, 0.77100], [2.375, 0.76400]])
LEGACY_ETA_CL61 = LEGACY_ETA_CL51.copy()      # legacy: CL61 used the CL51 table


def shipped_eta_table(instrument: str) -> NDArray:
    """The operational (PVC / Hogan 2006, a_G = 5.5 um) eta table for an instrument type."""
    return {"CL31": _CF._ETA_CL31, "CL51": _CF._ETA_CL51, "CL61": _CF._ETA_CL61,
            "CHM15K": _CF._ETA_CHM15K, "CHM15k": _CF._ETA_CHM15K}[instrument]


def legacy_eta_table(instrument: str) -> NDArray:
    return {"CL31": LEGACY_ETA_CL31, "CL51": LEGACY_ETA_CL51,
            "CL61": LEGACY_ETA_CL61}[instrument]


def eta_profile(table: NDArray, range_m: NDArray) -> NDArray:
    """eta(range) exactly as the pipeline evaluates it (linear interp + linear extrap)."""
    return _CF._interp1_linear_extrap(table[:, 0], table[:, 1], range_m / 1000.0)


# =======================================================================================
#  Forward model
# =======================================================================================
@dataclass
class Scene:
    """A synthetic measurement scene (one homogeneous cloud, repeated in time)."""
    range_m: NDArray                 # (n_range,) m AGL, uniform grid
    beta_att: NDArray                # (n_range,) m^-1 sr^-1, what the instrument measures
    alpha_mol: NDArray               # m^-1
    alpha_aer: NDArray               # m^-1
    alpha_cld: NDArray               # m^-1
    beta_mol: NDArray                # m^-1 sr^-1
    beta_aer: NDArray
    beta_cld: NDArray
    eta_true: NDArray                # (n_range,) the eta used IN THE PHYSICS
    cbh_m: float
    aod_below: float                 # aerosol optical depth 0 -> CBH (one-way)
    tau_cloud: float                 # cloud optical depth (geometric, no eta)


def make_range_grid(dr: float = 10.0, top_m: float = 5000.0) -> NDArray:
    """Operational working-grid range vector: dr = 10 m, first gate at dr."""
    return np.arange(1, int(round(top_m / dr)) + 1) * dr


def molecular_columns(range_m: NDArray, station_alt_m: float,
                      wavelength_m: float) -> Tuple[NDArray, NDArray]:
    """(beta_mol, alpha_mol) from the repo's US-Std-1976 + Bucholtz-1995 code."""
    prof = load_standard_atmosphere(None, range_m + station_alt_m)
    mol = calculate_molecular_properties(prof.temperature, prof.pressure, range_m, wavelength_m)
    return mol.beta_mol, mol.beta_mol * MOLECULAR_LIDAR_RATIO


def adiabatic_cloud_extinction(range_m: NDArray, cbh_m: float, thickness_m: float,
                               gamma_l: float = GAMMA_L, n_drop: float = N_DROP,
                               alpha_const: Optional[float] = None) -> NDArray:
    """Adiabatic-ish liquid cloud extinction [m^-1].

    LWC(z) = gamma_l * (z - CBH)  [g m^-3];  r_eff = (3 LWC / (4 pi rho_w N))^(1/3) [m];
    alpha = 3 LWC / (2 rho_w r_eff)  ->  alpha ~ (z - CBH)^(2/3).
    ``alpha_const`` (m^-1) overrides the adiabatic ramp with a slab of constant extinction
    (used to check the sensitivity of the closure to the in-cloud shape).
    """
    alpha = np.zeros_like(range_m)
    inside = (range_m >= cbh_m) & (range_m <= cbh_m + thickness_m)
    if alpha_const is not None:
        alpha[inside] = alpha_const
        return alpha
    dz = np.clip(range_m[inside] - cbh_m, 0.0, None)
    lwc = gamma_l * dz                                    # g m^-3
    with np.errstate(invalid="ignore", divide="ignore"):
        r_eff = (3.0 * lwc / (4.0 * np.pi * RHO_W * n_drop)) ** (1.0 / 3.0)   # m
        a = np.where(r_eff > 0, 3.0 * lwc / (2.0 * RHO_W * r_eff), 0.0)
    alpha[inside] = np.nan_to_num(a)
    return alpha


def exponential_aerosol(range_m: NDArray, cbh_m: float, aod_below: float,
                        scale_height_m: float = 1000.0) -> NDArray:
    """Aerosol extinction [m^-1]: exponential layer confined below the cloud base,
    normalised so that INT_0^CBH alpha_aer dz == aod_below."""
    alpha = np.zeros_like(range_m)
    if aod_below <= 0:
        return alpha
    dr = float(range_m[1] - range_m[0])
    below = range_m < cbh_m
    shape = np.exp(-range_m / scale_height_m) * below
    norm = np.sum(shape) * dr
    if norm <= 0:
        return alpha
    return shape * (aod_below / norm)


def simulate_scene(range_m: NDArray, cbh_m: float, *,
                   instrument: str = "CL61",
                   eta_table_true: Optional[NDArray] = None,
                   station_alt_m: float = 491.0,
                   thickness_m: float = 500.0,
                   aod_below: float = 0.0,
                   aer_scale_height_m: float = 1000.0,
                   s_cloud: float = S_CLOUD_TRUE,
                   lr_aer_true: float = LR_AER_TRUE,
                   alpha_const: Optional[float] = None,
                   multiple_scattering: bool = True,
                   include_molecular: bool = True,
                   n_sub: int = 20) -> Scene:
    """Build one synthetic attenuated-backscatter profile (the forward equation of the
    module docstring).  ``eta_table_true`` defaults to the instrument's shipped PVC table,
    i.e. the truth EQUALS what the pipeline will assume (the closure configuration).

    The profile is built on a sub-grid of ``n_sub`` points per output gate and then
    GATE-AVERAGED, which is both the honest detector model (the receiver integrates the
    returned power over the gate) and what keeps the numerical closure error negligible:
    an adiabatic cloud attenuates within ~30 m, i.e. ~3 gates of 10 m, so evaluating
    beta_att only at gate centres costs ~0.5 % on the integral.
    """
    dr = float(range_m[1] - range_m[0])
    wl = WAVELENGTH_M[instrument]
    # --- fine sub-grid: n_sub points per gate, centred on the gate ---
    off = (np.arange(n_sub) + 0.5) / n_sub - 0.5                     # in units of dr
    r_f = (range_m[:, None] + off[None, :] * dr).ravel()             # (n_range*n_sub,)
    dr_f = dr / n_sub

    beta_mol_f, alpha_mol_f = molecular_columns(r_f, station_alt_m, wl)
    if not include_molecular:      # diagnostic: isolate the pure numerical closure error
        beta_mol_f = np.zeros_like(beta_mol_f)
        alpha_mol_f = np.zeros_like(alpha_mol_f)
    alpha_aer_f = exponential_aerosol(r_f, cbh_m, aod_below, aer_scale_height_m)
    beta_aer_f = alpha_aer_f / lr_aer_true
    alpha_cld_f = adiabatic_cloud_extinction(r_f, cbh_m, thickness_m,
                                             alpha_const=alpha_const)
    beta_cld_f = alpha_cld_f / s_cloud

    if eta_table_true is None:
        eta_table_true = shipped_eta_table(instrument)
    eta_f = (eta_profile(eta_table_true, r_f) if multiple_scattering
             else np.ones_like(r_f))

    # one-way optical depth at the CENTRE of each sub-gate
    ext_eff = alpha_mol_f + alpha_aer_f + eta_f * alpha_cld_f
    tau_f = np.cumsum(ext_eff) * dr_f - 0.5 * ext_eff * dr_f
    beta_att_f = (beta_mol_f + beta_aer_f + beta_cld_f) * np.exp(-2.0 * tau_f)

    def _gate_mean(x):
        return x.reshape(len(range_m), n_sub).mean(axis=1)

    return Scene(range_m=range_m, beta_att=_gate_mean(beta_att_f),
                 alpha_mol=_gate_mean(alpha_mol_f), alpha_aer=_gate_mean(alpha_aer_f),
                 alpha_cld=_gate_mean(alpha_cld_f), beta_mol=_gate_mean(beta_mol_f),
                 beta_aer=_gate_mean(beta_aer_f), beta_cld=_gate_mean(beta_cld_f),
                 eta_true=_gate_mean(eta_f), cbh_m=cbh_m,
                 aod_below=float(np.sum(alpha_aer_f) * dr_f),
                 tau_cloud=float(np.sum(alpha_cld_f) * dr_f))


# =======================================================================================
#  Retrieval through the SHIPPED pipeline
# =======================================================================================
def _ceilo_data_from_scene(scene: Scene, c_true: float, c_applied: float,
                           n_time: int = 12, station_alt_m: float = 491.0,
                           jitter_pct: float = 0.0, seed: int = 0) -> CeiloData:
    """Wrap a scene as the CeiloData the pipeline consumes.

    rcs = C_true * beta_att  ->  the pipeline is fed beta = rcs / C_applied.
    ``jitter_pct`` adds a small multiplicative profile-to-profile jitter (still well
    inside the +-10 % temporal-consistency window) so the consistency filter is exercised
    on a non-degenerate series.
    """
    rng = np.random.default_rng(seed)
    scale = c_true / c_applied
    beta = np.repeat((scene.beta_att * scale)[:, None], n_time, axis=1)
    if jitter_pct > 0:
        beta = beta * (1.0 + jitter_pct / 100.0 * rng.standard_normal(n_time))[None, :]
    t0 = np.datetime64("2026-01-15T22:00:00")
    time = t0 + np.arange(n_time) * np.timedelta64(30, "s")
    time_num = 19738.0 + np.arange(n_time) * 30.0 / 86400.0
    return CeiloData(
        time=time, time_num=time_num, station_altitude=station_alt_m,
        station_latitude=46.81, station_longitude=6.94,
        range=scene.range_m.copy(), range_resol=float(scene.range_m[1] - scene.range_m[0]),
        beta=beta, cbh=np.full(n_time, scene.cbh_m), quality_flag=None,
        window_transmission=None, laser_energy=None, trans2_wv=None,
        calibration_constant_applied=c_applied)


def _make_config(instrument: str, *, transmission_correction: bool,
                 cal_maxheight: float = 2400.0, cbh_maxheight: float = 2400.0,
                 cal_minheight: float = 100.0, cbh_minheight: float = 500.0,
                 ratio_filter: float = 0.05) -> CloudCalConfig:
    return set_defaults(CloudCalConfig(
        instrument=instrument,
        apply_wv_correction=False,            # WV is a separate, multiplicative effect
        apply_transmission_correction=transmission_correction,
        aerosol_lidar_ratio=50.0,             # operational runner value
        cal_minheight=cal_minheight, cal_maxheight=cal_maxheight,
        cbh_minheight=cbh_minheight, cbh_maxheight=cbh_maxheight,
        ratio_filter=ratio_filter,
        average_time_s=0.0, average_range_m=0.0))


class eta_table_override:
    """Context manager patching the eta table the PIPELINE uses (the tables live in
    calibration.cloud._filters; apply_multiple_scattering_correction reads them by name
    at call time, so patching the module attribute is enough)."""

    _ATTR = {"CL31": "_ETA_CL31", "CL51": "_ETA_CL51", "CL61": "_ETA_CL61",
             "CHM15K": "_ETA_CHM15K", "CHM15k": "_ETA_CHM15K"}

    def __init__(self, instrument: str, table: Optional[NDArray]):
        self.attr = self._ATTR[instrument]
        self.table = table

    def __enter__(self):
        self.saved = getattr(_CF, self.attr).copy()
        if self.table is not None:
            setattr(_CF, self.attr, np.asarray(self.table, dtype="float64"))
        return self

    def __exit__(self, *exc):
        setattr(_CF, self.attr, self.saved)
        return False


def retrieve(scene: Scene, *, instrument: str = "CL61", c_true: float = 1.25,
             c_applied: Optional[float] = None, transmission_correction: bool = False,
             eta_table_retrieval: Optional[NDArray] = None,
             config_kw: Optional[dict] = None, n_time: int = 12,
             jitter_pct: float = 0.0) -> Dict[str, float]:
    """Push a scene through the unmodified shipped cloud calibration.

    Returns a dict with the retrieved Wiegner constant, the O'Connor multiplier, the
    apparent lidar ratio, the relative error vs C_true and the rejection stats.
    """
    if c_applied is None:
        c_applied = INSTRUMENT_CAL_DEFAULT.get(instrument, 1.0)
    data = _ceilo_data_from_scene(scene, c_true, c_applied, n_time=n_time,
                                  jitter_pct=jitter_pct)
    cfg = _make_config(instrument, transmission_correction=transmission_correction,
                       **(config_kw or {}))
    with eta_table_override(instrument, eta_table_retrieval):
        res = liquid_cloud_calibration_from_data(data, cfg)
    c_l = float(res.lidar_constant)
    out = {
        "cbh_m": scene.cbh_m,
        "C_L": c_l,
        "C_oconnor": float(res.calibration_coefficient),
        "S_apparent": float(np.nanmedian(res.S_apparent)) if res.S_apparent is not None else np.nan,
        "n_profiles": int(res.n_profiles),
        "err_pct": 100.0 * (c_l / c_true - 1.0) if np.isfinite(c_l) else np.nan,
        "aod_below": scene.aod_below,
    }
    stats = dict(res.cloud_stats or {})
    stats.update(res.consistency_stats or {})
    out["rejects"] = stats
    return out


def sweep_cbh(cbh_grid: NDArray, *, instrument: str = "CL61", c_true: float = 1.25,
              eta_table_true: Optional[NDArray] = None,
              eta_table_retrieval: Optional[NDArray] = None,
              transmission_correction: bool = False, aod_below: float = 0.0,
              thickness_m: float = 500.0, config_kw: Optional[dict] = None,
              range_m: Optional[NDArray] = None,
              alpha_const: Optional[float] = None,
              multiple_scattering: bool = True,
              include_molecular: bool = True) -> Dict[str, NDArray]:
    """Retrieve C_L over a grid of cloud-base heights. Returns arrays keyed cbh/C_L/err_pct."""
    if range_m is None:
        range_m = make_range_grid()
    rows = []
    for cbh in cbh_grid:
        sc = simulate_scene(range_m, float(cbh), instrument=instrument,
                            eta_table_true=eta_table_true, thickness_m=thickness_m,
                            aod_below=aod_below, alpha_const=alpha_const,
                            multiple_scattering=multiple_scattering,
                            include_molecular=include_molecular)
        rows.append(retrieve(sc, instrument=instrument, c_true=c_true,
                             transmission_correction=transmission_correction,
                             eta_table_retrieval=eta_table_retrieval,
                             config_kw=config_kw))
    return {
        "cbh_m": np.array([r["cbh_m"] for r in rows]),
        "C_L": np.array([r["C_L"] for r in rows]),
        "err_pct": np.array([r["err_pct"] for r in rows]),
        "S_apparent": np.array([r["S_apparent"] for r in rows]),
        "n_profiles": np.array([r["n_profiles"] for r in rows]),
        "rows": rows,
    }


def slope_pct_per_km(cbh_m: NDArray, err_pct: NDArray,
                     lo: float = 500.0, hi: float = 2400.0) -> Tuple[float, float]:
    """Least-squares dC/dCBH in %/km over [lo, hi] m, plus the peak-to-peak span in %."""
    m = np.isfinite(err_pct) & (cbh_m >= lo) & (cbh_m <= hi)
    if m.sum() < 2:
        return float("nan"), float("nan")
    p = np.polyfit(cbh_m[m] / 1000.0, err_pct[m], 1)
    return float(p[0]), float(np.nanmax(err_pct[m]) - np.nanmin(err_pct[m]))


def beta_weighted_height(scene: Scene, lo: float = 100.0, hi: float = 2400.0) -> float:
    """Centroid height [m AGL] of the integrand the estimator actually sees, i.e. the
    altitude at which the eta table is EFFECTIVELY evaluated for that profile."""
    m = (scene.range_m >= lo) & (scene.range_m <= hi)
    w = scene.beta_att[m]
    return float(np.sum(w * scene.range_m[m]) / np.sum(w))


# =======================================================================================
#  Experiments
# =======================================================================================
FIG_DIR = _REPO / "doc" / "reports" / "figs_altitude_audit"
C_TRUE = 1.25            # arbitrary "true" Wiegner constant (C_applied for CL61 is 1.0)


def _fmt(sweep, lo=500.0, hi=2400.0) -> str:
    s, span = slope_pct_per_km(sweep["cbh_m"], sweep["err_pct"], lo, hi)
    e = sweep["err_pct"]
    m = np.isfinite(e) & (sweep["cbh_m"] >= lo) & (sweep["cbh_m"] <= hi)
    if not m.any():
        return f"NO CBH CALIBRATED (all {e.size} rejected by the pipeline filters)"
    return (f"max|dev| = {np.nanmax(np.abs(e[m])):6.3f} %   "
            f"dC/dCBH = {s:+7.3f} %/km   span = {span:6.3f} %   "
            f"(n = {m.sum()}/{e.size} CBH)")


def run_experiments(quick: bool = False) -> dict:
    out = {}
    dcbh = 200.0 if quick else 100.0
    cbh_op = np.arange(500.0, 2400.0 + 1, dcbh)          # operational CBH gate 500-2400 m
    cbh_ext = np.arange(500.0, 3000.0 + 1, dcbh)         # extended, for the 0.5-3 km ask
    ext_cfg = dict(cal_maxheight=5000.0, cbh_maxheight=3200.0)
    rg = make_range_grid(10.0, 6000.0)

    print("=" * 100)
    print("FORWARD MODEL -- LIQUID-CLOUD (O'Connor/Hopkin) CLOSURE")
    print(f"C_true = {C_TRUE}, S_cloud_true = {S_CLOUD_TRUE} sr, adiabatic cloud "
          f"(Gamma_l = {GAMMA_L*1e3:.1f} g m^-3 km^-1, N = {N_DROP/1e6:.0f} cm^-3, 500 m thick),")
    print("US Std 1976 molecular atmosphere, station 491 m ASL, 10 m x 30 s grid, "
          "eta_true = eta_table (PVC).")
    print("=" * 100)

    # --- E1: closure -------------------------------------------------------------------
    print("\n[E1] CLOSURE  (retrieval eta table == true eta; transmission correction OFF)")
    for inst in ("CL61", "CL31", "CHM15k"):
        num = sweep_cbh(cbh_ext, instrument=inst, c_true=C_TRUE, range_m=rg,
                        config_kw=ext_cfg, include_molecular=False)
        full = sweep_cbh(cbh_ext, instrument=inst, c_true=C_TRUE, range_m=rg,
                         config_kw=ext_cfg)
        opr = sweep_cbh(cbh_op, instrument=inst, c_true=C_TRUE, range_m=rg)
        print(f"  {inst:7s} cloud+eta only (no molecular), CBH 0.5-3.0 km : {_fmt(num, 500, 3000)}")
        print(f"  {inst:7s} + molecular atmosphere,        CBH 0.5-3.0 km : {_fmt(full, 500, 3000)}")
        print(f"  {inst:7s} operational gate 100-2400 m,   CBH 0.5-2.4 km : {_fmt(opr)}")
        print(f"  {inst:7s} operational gate, UNTRUNCATED, CBH 0.5-2.2 km : {_fmt(opr, 500, 2200)}")
        out[f"closure_num_{inst}"] = num
        out[f"closure_{inst}"] = full
        out[f"closure_op_{inst}"] = opr

    # --- E2: eta-table mismatch --------------------------------------------------------
    print("\n[E2] ETA-TABLE MISMATCH  (legacy Hopkin a_G=8um  vs  operational PVC a_G=5.5um)")
    for inst in ("CL61", "CL31"):
        pvc, leg = shipped_eta_table(inst), legacy_eta_table(inst)
        e_pvc, e_leg = eta_profile(pvc, cbh_op), eta_profile(leg, cbh_op)
        print(f"  {inst}: eta_PVC {e_pvc[0]:.4f}->{e_pvc[-1]:.4f}, "
              f"eta_legacy {e_leg[0]:.4f}->{e_leg[-1]:.4f} over CBH 0.5->2.4 km; "
              f"table error {100*(e_leg[0]/e_pvc[0]-1):+.2f} % at 0.5 km, "
              f"{100*(e_leg[-1]/e_pvc[-1]-1):+.2f} % at 2.4 km")
        a = sweep_cbh(cbh_op, instrument=inst, c_true=C_TRUE, range_m=rg,
                      eta_table_true=pvc, eta_table_retrieval=leg)
        b = sweep_cbh(cbh_op, instrument=inst, c_true=C_TRUE, range_m=rg,
                      eta_table_true=leg, eta_table_retrieval=pvc)
        print(f"    truth PVC   / retrieval LEGACY : {_fmt(a)}")
        print(f"    truth LEGACY/ retrieval PVC    : {_fmt(b)}")
        out[f"mismatch_legacy_{inst}"] = a
        out[f"mismatch_pvc_{inst}"] = b

    # --- E2b: pure eta scale error -> mapping x% eta -> y%/km --------------------------
    print("\n[E2b] MAPPING  eta table error -> C_L error (constant scale and pure slope error)")
    inst = "CL61"
    base = shipped_eta_table(inst)
    for e_pct in (-5.0, -2.0, 2.0, 5.0):
        tab = base.copy()
        tab[:, 1] = tab[:, 1] * (1.0 + e_pct / 100.0)
        s = sweep_cbh(cbh_op, instrument=inst, c_true=C_TRUE, range_m=rg,
                      eta_table_true=base, eta_table_retrieval=tab)
        sl, _ = slope_pct_per_km(s["cbh_m"], s["err_pct"])
        print(f"  eta scaled by {e_pct:+5.1f} % at ALL heights -> C_L error "
              f"{np.nanmean(s['err_pct']):+6.3f} % (mean), dC/dCBH {sl:+6.3f} %/km")
    for tilt in (-5.0, 5.0):   # eta error growing linearly from 0 at 0.5 km to tilt at 2.4 km
        tab = base.copy()
        f = 1.0 + (tilt / 100.0) * (tab[:, 0] - 0.5) / (2.4 - 0.5)
        tab[:, 1] = tab[:, 1] * f
        s = sweep_cbh(cbh_op, instrument=inst, c_true=C_TRUE, range_m=rg,
                      eta_table_true=base, eta_table_retrieval=tab)
        sl, _ = slope_pct_per_km(s["cbh_m"], s["err_pct"])
        print(f"  eta error ramping 0 -> {tilt:+5.1f} % between 0.5 and 2.4 km -> "
              f"dC/dCBH {sl:+6.3f} %/km")

    # --- E3: fixed 100-2400 m gate: truncation ----------------------------------------
    print("\n[E3] FIXED 100-2400 m INTEGRATION GATE")
    print("  (a) GATE TOP: the cloud is cut off. What matters is not the cloud THICKNESS")
    print("      (an adiabatic cloud extinguishes the beam in ~50 m) but the in-cloud")
    print("      extinction, i.e. how many metres of cloud the gate must contain.")
    cbh_hi = np.arange(2100.0, 2400.0 + 1, 10.0)
    for a_c, lbl in ((None, "adiabatic ramp"), (0.050, "slab 50 /km"),
                     (0.010, "slab 10 /km")):
        s_ = sweep_cbh(cbh_hi, instrument="CL61", c_true=C_TRUE, range_m=rg,
                       alpha_const=a_c, thickness_m=1000.0)
        e = s_["err_pct"]
        ok = np.isfinite(e)
        def _first(th):
            bad = ok & (np.abs(e) > th)
            return s_["cbh_m"][bad].min() if bad.any() else np.nan
        print(f"    {lbl:16s}: |err| > 1 % above CBH {_first(1.0):6.0f} m, "
              f">5 % above {_first(5.0):6.0f} m, "
              f"err at CBH 2350 m = {e[np.argmin(np.abs(s_['cbh_m']-2350))]:+8.2f} %, "
              f"at 2390 m = {e[np.argmin(np.abs(s_['cbh_m']-2390))]:+8.2f} %")
        out[f"trunc_{lbl}"] = s_
    sc = simulate_scene(rg, 2300.0, instrument="CL61")
    m = (rg >= 100) & (rg <= 2400)
    frac = np.sum(sc.beta_att[m]) / np.sum(sc.beta_att[rg >= 100])
    print(f"    adiabatic cloud, CBH 2300 m: 100 m of cloud inside the gate -> "
          f"{100*frac:.2f} % of the column integrated (err "
          f"{retrieve(sc, c_true=C_TRUE)['err_pct']:+.2f} %)")

    print("  (b) GATE BOTTOM (100 m): reached only if the CBH gate is opened below ~500 m")
    cbh_lo = np.arange(120.0, 620.0 + 1, 20.0)
    s_lo = sweep_cbh(cbh_lo, instrument="CL61", c_true=C_TRUE, range_m=rg,
                     config_kw=dict(cbh_minheight=100.0))
    txt = "  ".join(f"{c:.0f}m:{e:+.2f}%" for c, e in zip(s_lo["cbh_m"], s_lo["err_pct"])
                    if c <= 300 or c % 100 == 0)
    print(f"    {txt}")
    out["trunc_low"] = s_lo

    # --- E4: molecular signal treated as aerosol by the transmission correction --------
    print("\n[E4] AEROSOL TRANSMISSION CORRECTION (operational: ON, LR_aer = 50 sr)")
    for inst in ("CL61", "CHM15k"):
        s = sweep_cbh(cbh_op, instrument=inst, c_true=C_TRUE, range_m=rg,
                      transmission_correction=True)
        print(f"  {inst:7s} CLEAN air, correction ON : {_fmt(s)}   "
              f"<- pure molecular column mis-read as aerosol")
        out[f"transm_{inst}"] = s

    # --- E5: below-cloud aerosol -------------------------------------------------------
    print("\n[E5] BELOW-CLOUD AEROSOL (exponential layer, H = 1000 m, LR_true = 50 sr)")
    aods = [0.0, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.12, 0.20]
    for tc in (False, True):
        print(f"  transmission correction {'ON ' if tc else 'OFF'}:")
        for aod in aods:
            s = sweep_cbh(cbh_op, instrument="CL61", c_true=C_TRUE, range_m=rg,
                          aod_below=aod, transmission_correction=tc)
            ok = np.isfinite(s["err_pct"])
            at1 = s["err_pct"][np.argmin(np.abs(s["cbh_m"] - 1000))]
            at2 = s["err_pct"][np.argmin(np.abs(s["cbh_m"] - 2000))]
            print(f"    AOD(0->CBH) = {aod:4.2f} : err {at1:+6.2f} % @CBH 1.0 km, "
                  f"{at2:+6.2f} % @CBH 2.0 km, calibrated {ok.sum():2d}/{ok.size} CBH "
                  f"(rest rejected by the 5 % below-cloud ratio filter)")
            out[f"aer_{'on' if tc else 'off'}_{aod:.2f}"] = s
        # AOD needed for a +-5 % bias at CBH 1 km, and the AOD the 5 % ratio filter allows
        e1 = np.array([out[f"aer_{'on' if tc else 'off'}_{a:.2f}"]["err_pct"][
            np.argmin(np.abs(cbh_op - 1000))] for a in aods])
        keep = np.array([np.isfinite(out[f"aer_{'on' if tc else 'off'}_{a:.2f}"]["err_pct"]).sum()
                         for a in aods])
        ok = np.isfinite(e1)
        d = e1[ok] - e1[0]
        aod_5 = np.interp(5.0, np.abs(d), np.array(aods)[ok]) if np.nanmax(np.abs(d)) >= 5 else np.nan
        aod_max = np.array(aods)[keep > 0].max()
        print(f"    -> AOD for a 5 %% shift vs clean air @CBH 1 km: "
              f"{aod_5 if np.isfinite(aod_5) else float('nan'):.3f}"
              f"{'' if np.isfinite(aod_5) else '  (NEVER reached: filter rejects first)'}; "
              f"largest AOD surviving the ratio filter = {aod_max:.2f}, "
              f"max achievable |shift| = {np.nanmax(np.abs(d)):.2f} %%")

    # --- E6: sanity on the in-cloud shape ---------------------------------------------
    print("\n[E6] SENSITIVITY TO THE IN-CLOUD EXTINCTION SHAPE (closure must not care)")
    for a_c, lbl in ((0.005, "slab alpha =  5 /km"), (0.010, "slab alpha = 10 /km"),
                     (0.050, "slab alpha = 50 /km"), (None, "adiabatic ramp")):
        s = sweep_cbh(cbh_op, instrument="CL61", c_true=C_TRUE, range_m=rg,
                      alpha_const=a_c, thickness_m=1000.0)
        rej = s["rows"][len(s["rows"]) // 2]["rejects"]
        print(f"  {lbl:22s}: {_fmt(s, 500, 2200)}"
              + ("" if np.isfinite(s["err_pct"]).any() else f"   rejects={rej}"))
        out[f"shape_{lbl}"] = s

    print("\n[E7] EFFECTIVE HEIGHT AT WHICH THE ETA TABLE IS EVALUATED "
          "(backscatter-weighted centroid of the integrand)")
    for cbh in (600.0, 1200.0, 2000.0):
        sc = simulate_scene(rg, cbh, instrument="CL61")
        zb = beta_weighted_height(sc)
        print(f"  CBH {cbh:6.0f} m -> centroid {zb:7.1f} m  (CBH + {zb-cbh:5.1f} m); "
              f"eta_table(centroid) = {eta_profile(shipped_eta_table('CL61'), np.array([zb]))[0]:.5f}, "
              f"eta_table(CBH) = {eta_profile(shipped_eta_table('CL61'), np.array([cbh]))[0]:.5f}")

    out["_grids"] = dict(cbh_op=cbh_op, cbh_ext=cbh_ext, rg=rg, ext_cfg=ext_cfg)
    return out


# =======================================================================================
#  Figure
# =======================================================================================
def make_figure(res: dict, path=None) -> Path:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    path = Path(path or (FIG_DIR / "fwd_cloud_closure.png"))
    path.parent.mkdir(parents=True, exist_ok=True)
    rg = res["_grids"]["rg"]

    fig, axes = plt.subplots(2, 4, figsize=(23.0, 11.0))
    ax = axes.ravel()

    # (a) extinction profiles
    a = ax[0]
    for cbh, c in ((800.0, "tab:blue"), (1800.0, "tab:red")):
        sc = simulate_scene(rg, cbh, instrument="CL61", aod_below=0.08)
        a.plot(sc.alpha_cld * 1e3, sc.range_m, color=c, lw=1.8,
               label=f"cloud, CBH {cbh:.0f} m")
        a.plot(sc.alpha_aer * 1e3, sc.range_m, color=c, lw=1.2, ls="--",
               label=f"aerosol (AOD 0.08), CBH {cbh:.0f} m")
    a.plot(sc.alpha_mol * 1e3, sc.range_m, color="k", lw=1.2, ls=":", label="molecular")
    a.set_xscale("log")
    a.set_xlim(1e-3, 3e2)
    a.set_ylim(0, 3000)
    a.set_xlabel(r"extinction $\alpha$  [km$^{-1}$]")
    a.set_ylabel("altitude AGL  [m]")
    a.set_title("(a) forward model: extinction")

    # (b) simulated attenuated backscatter
    b = ax[1]
    for cbh, c in ((800.0, "tab:blue"), (1800.0, "tab:red")):
        sc = simulate_scene(rg, cbh, instrument="CL61", aod_below=0.08)
        b.plot(sc.beta_att, sc.range_m, color=c, lw=1.8, label=f"CBH {cbh:.0f} m")
    b.axhspan(100, 2400, color="0.85", zorder=0, label="integration gate 100-2400 m")
    b.set_xscale("log")
    b.set_xlim(1e-10, 1e-2)
    b.set_ylim(0, 3000)
    b.set_xlabel(r"$\beta_{att}$ measured  [m$^{-1}$ sr$^{-1}$]")
    b.set_ylabel("altitude AGL  [m]")
    b.set_title(r"(b) simulated $\beta_{att}$")

    # (c) eta tables
    c_ax = ax[2]
    zz = np.linspace(200.0, 3000.0, 200)
    for inst, col in (("CL61", "tab:green"), ("CL31", "tab:purple")):
        c_ax.plot(eta_profile(shipped_eta_table(inst), zz), zz, color=col, lw=2.0,
                  label=f"{inst} PVC 5.5 um (operational)")
        c_ax.plot(eta_profile(legacy_eta_table(inst), zz), zz, color=col, lw=1.5,
                  ls="--", label=f"{inst} legacy Hopkin 8 um")
    c_ax.plot(eta_profile(shipped_eta_table("CHM15k"), zz), zz, color="tab:orange",
              lw=2.0, label="CHM15k PVC 5.5 um")
    c_ax.set_xlabel(r"multiple-scattering factor $\eta$")
    c_ax.set_ylabel("altitude AGL  [m]")
    c_ax.set_title(r"(c) $\eta$ tables: truth vs correction")

    # (d) closure
    d = ax[3]
    for inst, col in (("CL61", "tab:green"), ("CL31", "tab:purple"),
                      ("CHM15k", "tab:orange")):
        s = res[f"closure_{inst}"]
        d.plot(s["err_pct"], s["cbh_m"], color=col, lw=1.8, marker="o", ms=3,
               label=f"{inst}, full physics")
        n = res[f"closure_num_{inst}"]
        d.plot(n["err_pct"], n["cbh_m"], color=col, lw=1.0, ls=":",
               label=f"{inst}, no molecular")
    d.axvline(0, color="k", lw=1.0)
    d.set_xlim(-1.0, 1.0)
    d.set_xlabel(r"retrieved $C_L$ error  [%]")
    d.set_ylabel("cloud base height  [m AGL]")
    d.set_title(r"(d) CLOSURE: table $\eta$ = true $\eta$")

    # (e) eta mismatch
    e = ax[4]
    for inst, col in (("CL61", "tab:green"), ("CL31", "tab:purple")):
        s = res[f"mismatch_legacy_{inst}"]
        e.plot(s["err_pct"], s["cbh_m"], color=col, lw=2.0,
               label=f"{inst}: truth PVC, retrieved LEGACY")
        s = res[f"mismatch_pvc_{inst}"]
        e.plot(s["err_pct"], s["cbh_m"], color=col, lw=1.5, ls="--",
               label=f"{inst}: truth LEGACY, retrieved PVC")
    e.axvline(0, color="k", lw=1.0)
    e.set_xlabel(r"retrieved $C_L$ error  [%]")
    e.set_ylabel("cloud base height  [m AGL]")
    e.set_title(r"(e) $\eta$-table error $\rightarrow$ altitude dependence")

    # (f) gate truncation
    f = ax[5]
    for lbl, ls, col in (("adiabatic ramp", "-", "tab:red"),
                         ("slab 50 /km", "--", "tab:brown"),
                         ("slab 10 /km", ":", "tab:gray")):
        s = res[f"trunc_{lbl}"]
        f.plot(s["err_pct"], s["cbh_m"], color=col, lw=1.8, ls=ls,
               marker="o", ms=3, label=lbl)
    s = res["trunc_low"]
    f.plot(s["err_pct"], s["cbh_m"], color="tab:blue", lw=1.8, marker="s", ms=3,
           label="gate bottom (CBH gate opened to 100 m)")
    f.set_xlim(-30, 5)
    f.axhline(2400, color="k", lw=1.0, ls=":", label="gate top / CBH gate 2400 m")
    f.axvline(0, color="k", lw=1.0)
    f.set_xlabel(r"retrieved $C_L$ error  [%]")
    f.set_ylabel("cloud base height  [m AGL]")
    f.set_title("(f) fixed gate top: cloud truncation")

    # (g) transmission correction on clean air
    g = ax[6]
    for inst, col in (("CL61", "tab:green"), ("CHM15k", "tab:orange")):
        s = res[f"transm_{inst}"]
        g.plot(s["err_pct"], s["cbh_m"], color=col, lw=2.0,
               label=f"{inst}, correction ON")
        s0 = res[f"closure_op_{inst}"]
        g.plot(s0["err_pct"], s0["cbh_m"], color=col, lw=1.2, ls=":",
               label=f"{inst}, correction OFF")
    g.axvline(0, color="k", lw=1.0)
    g.set_xlabel(r"retrieved $C_L$ error  [%]")
    g.set_ylabel("cloud base height  [m AGL]")
    g.set_title("(g) transmission correction, CLEAN air\n(molecular column read as aerosol)")

    # (h) below-cloud aerosol
    h = ax[7]
    cmap = plt.get_cmap("viridis")
    aods = [0.0, 0.02, 0.04, 0.05, 0.06, 0.08]
    for i, aod in enumerate(aods):
        s = res[f"aer_on_{aod:.2f}"]
        h.plot(s["err_pct"], s["cbh_m"], color=cmap(i / (len(aods) - 1)), lw=1.8,
               marker="o", ms=2.5, label=f"AOD {aod:.2f}")
    h.axvline(0, color="k", lw=1.0)
    h.set_xlabel(r"retrieved $C_L$ error  [%]")
    h.set_ylabel("cloud base height  [m AGL]")
    h.set_title("(h) below-cloud aerosol\n(operational config, correction ON)")

    for a_ in ax:
        a_.grid(alpha=0.3)
        a_.legend(fontsize=7, loc="best")
    fig.suptitle("Forward model of the O'Connor/Hopkin liquid-cloud calibration: "
                 f"closure and altitude dependence (C_true = {C_TRUE}, S_cloud = 18.8 sr)",
                 fontsize=13)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    fig.savefig(path, dpi=140)
    plt.close(fig)
    return path


def main() -> None:
    ap = argparse.ArgumentParser(description="Cloud forward model / closure test")
    ap.add_argument("--quick", action="store_true", help="coarse CBH grid (200 m)")
    ap.add_argument("--no-fig", action="store_true")
    args = ap.parse_args()
    res = run_experiments(quick=args.quick)
    if not args.no_fig:
        p = make_figure(res)
        print(f"\nfigure -> {p}")


if __name__ == "__main__":
    main()
