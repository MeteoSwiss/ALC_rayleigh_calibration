import cProfile, pstats, io, time, warnings, logging
warnings.filterwarnings("ignore")
logging.getLogger().setLevel(logging.ERROR)

from netCDF4 import Dataset
import numpy as np
import run_lindenberg_cl61_cal as m
from rayleigh_calibration.cloud_calibration import liquid_cloud_calibration, CloudCalConfig

nc = str(m.DATA_ROOT / m.WMO / "20241001.nc")
with Dataset(nc) as ds:
    t = ds.variables["time"]
    r = ds.variables["range"]
    tv = t[:]
    rv = r[:]
    dt = np.median(np.diff(tv.astype("float64")))
    print(f"file: {nc}")
    print(f"  n_time={tv.size}  n_range={rv.size}")
    print(f"  time units: {getattr(t,'units','?')}")
    print(f"  median dt (raw units) = {dt:.4f}")
    print(f"  range res = {float(rv[1]-rv[0]):.2f} m, top = {float(rv[-1]):.0f} m")

cfg = CloudCalConfig(nc_file=nc, instrument="CL61", apply_wv_correction=True,
                     cams_folder=m.CAMS_FOLDER, abs_cs_lookup_table=m.ABS_CS_LUT,
                     station_latitude=m.SITE["lat"], station_longitude=m.SITE["lon"],
                     aerosol_lidar_ratio=50.0)

t0 = time.time()
pr = cProfile.Profile()
pr.enable()
res = liquid_cloud_calibration(cfg)
pr.disable()
dt_run = time.time() - t0
print(f"\nWV ON: elapsed = {dt_run:.1f} s, n_profiles={res.n_profiles}, cal_mean={res.cal_mean}")

s = io.StringIO()
ps = pstats.Stats(pr, stream=s).sort_stats("cumulative")
ps.print_stats(20)
print(s.getvalue())
