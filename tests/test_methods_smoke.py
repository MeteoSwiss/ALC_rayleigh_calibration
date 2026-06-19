"""Synthetic smoke test for molecular_methods: all 6 methods + temporal aerosol rejection."""
import numpy as np
from calibration.rayleigh.molecular_methods import (
    METHODS, compute_window_grid, select_molecular_window,
)

rng = np.random.default_rng(0)
dz = 30.0
z = np.arange(0, 8000, dz)                  # range AGL
H = 8000.0
p_mol = np.exp(-z / H)                       # molecular power ~ density (relative)
CL = 2.0

# Aerosol layer 0-1.8 km (decays with height), molecular elsewhere.
aer = 1.5 * np.exp(-z / 700.0) * (z < 1800)

nT = 80
stack = np.empty((nT, z.size))
for t in range(nT):
    # molecular is steady; aerosol fluctuates in time (advection); photon noise grows w/ range
    aer_t = aer * (1.0 + 0.6 * np.sin(t / 5.0) + 0.3 * rng.standard_normal())
    noise = rng.standard_normal(z.size) * (0.002 + 0.02 * (z / H) ** 2)
    stack[t] = CL * p_mol + aer_t + noise
signal = np.nanmean(stack, axis=0)

half = tuple(range(250, 2000, 240))
grid = compute_window_grid(signal, p_mol, z, half, range_start_m=2000, range_end_m=6000,
                           increment_bins=8, signal_stack=stack)
print(f"grid {grid.r2.shape}  temporal_cv finite: {np.isfinite(grid.temporal_cv).sum()}")

for m in METHODS:
    w = select_molecular_window(m, signal, p_mol, z, half, range_start_m=2000,
                                range_end_m=6000, increment_bins=8, grid=grid)
    if w.ok:
        print(f"  {m:9s} ok  start={w.start_m:5.0f}-{w.end_m:5.0f} m  R2={w.r2:.3f} "
              f"CL={w.cl:.3f} Rscat={w.scattering_ratio:.2f} tcv={w.temporal_cv:.2f}  {w.message}")
    else:
        print(f"  {m:9s} FAIL  {w.message}")
print("SMOKE_OK")
