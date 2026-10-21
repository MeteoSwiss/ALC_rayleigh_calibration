# -*- coding: utf-8 -*-
"""q16e check 4: the night selector's solar zenith angle (data_loader.filter_time_range) vs NOAA.

The selector uses mean solar time (UTC + lon/15, no equation of time) and a cosine declination with
the CURRENT date's day of year for the whole night. Compare the SZA>100 deg window it selects with
the NOAA/Spencer algorithm for three latitudes over a year: shift of the window edges [min] and the
dark duration [h].
"""
import numpy as np

LATS = (45.0, 52.0, 56.0)
LON = 10.0


def code_sza(doy, solar_h, lat):
    decl = np.radians(-23.44) * np.cos(2 * np.pi * (doy + 10) / 365.25)
    ha = np.radians((solar_h - 12) * 15)
    la = np.radians(lat)
    return np.degrees(np.arccos(np.clip(np.sin(la) * np.sin(decl) + np.cos(la) * np.cos(decl) * np.cos(ha), -1, 1)))


def noaa_sza(doy, utc_h, lat, lon):
    g = 2 * np.pi / 365 * (doy - 1 + (utc_h - 12) / 24)
    eot = 229.18 * (0.000075 + 0.001868 * np.cos(g) - 0.032077 * np.sin(g) - 0.014615 * np.cos(2 * g) - 0.040849 * np.sin(2 * g))
    decl = (0.006918 - 0.399912 * np.cos(g) + 0.070257 * np.sin(g) - 0.006758 * np.cos(2 * g) + 0.000907 * np.sin(2 * g)
            - 0.002697 * np.cos(3 * g) + 0.00148 * np.sin(3 * g))
    tst = utc_h * 60 + eot + 4 * lon
    ha = np.radians(tst / 4 - 180)
    la = np.radians(lat)
    return np.degrees(np.arccos(np.clip(np.sin(la) * np.sin(decl) + np.cos(la) * np.cos(decl) * np.cos(ha), -1, 1)))


def main():
    # night of date D: UTC from D-1 12h to D 12h in 1-min steps (relative hours -12..12 around 00 UTC)
    rel = np.arange(-12 * 60, 12 * 60) / 60.0
    print("lat  doy  dark_h(code)  dark_h(NOAA)  start_shift[min]  end_shift[min]")
    for lat in LATS:
        for doy in (15, 80, 172, 266, 355):
            utc = np.mod(rel, 24)
            day = np.where(rel < 0, doy - 1, doy)
            solar = np.mod(utc + LON / 15, 24)
            c = code_sza(doy, solar, lat) > 100
            n = noaa_sza(day, utc, lat, LON) > 100
            def edges(m):
                if not m.any():
                    return np.nan, np.nan
                i = np.where(m)[0]
                return rel[i[0]] * 60, rel[i[-1]] * 60
            cs, ce = edges(c)
            ns, ne = edges(n)
            print(f"{lat:4.0f} {doy:4d} {c.sum() / 60:12.2f} {n.sum() / 60:13.2f} {cs - ns:16.1f} {ce - ne:14.1f}")


if __name__ == "__main__":
    main()
