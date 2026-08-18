from netCDF4 import Dataset
import numpy as np

nc = r"A:/CL61_Cloudnet/Lindenberg/20241001.nc"
with Dataset(nc) as ds:
    t = ds.variables["time"]
    r = ds.variables["range"]
    tv = np.asarray(t[:], dtype="float64")
    rv = np.asarray(r[:], dtype="float64")
    units = getattr(t, "units", "?")
    dt = float(np.median(np.diff(tv)))
    # convert dt to seconds
    if "second" in units:
        dt_s = dt
    elif "hour" in units:
        dt_s = dt * 3600
    elif "day" in units:
        dt_s = dt * 86400
    else:
        dt_s = dt
print("n_time =", tv.size)
print("n_range =", rv.size)
print("time units =", units)
print("median dt =", dt, "->", round(dt_s, 2), "s")
print("range res =", round(float(rv[1] - rv[0]), 3), "m")
print("range top =", round(float(rv[-1]), 1), "m")
print("beta matrix =", tv.size, "x", rv.size, "=", tv.size * rv.size, "elements")
