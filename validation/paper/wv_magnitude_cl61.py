"""wv_magnitude_cl61.py - Section 5.1.1: size of the water-vapour bias of the CL61 lidar constant.

Paper text: "left uncorrected, [water-vapour absorption] biases the retrieved constant low, by about 3-4 % for
the narrow emission line of the CL61". Ratio C_L(no WV correction) / C_L(WV correction with the operational
910.55 nm / 0.188 nm spectrum) on the nights calibrated by both runs:

    Payerne CL61     Rayleigh 0.962 (-3.8 %, 43 nights, Mar-Aug 2026)    cloud 0.973 (69 nights)
    Lindenberg CL61  Rayleigh 0.973 (-2.7 %, 190 nights, Jan 2025-Jun 2026)  cloud 0.980 (304 nights)

Inputs (calibration runs, <key>/<key>_cal.csv):
    DATA/diag_v22_l55s008  WV on, CL61 spectrum 910.55/0.188  (rayleigh_availability/run_lambda_correct.ps1)
    DATA/diag_v22_nowv     identical but ALC_WV_DISABLE=1      (rayleigh_availability/run_diag_v22_nowv.ps1)

Run:  python -m validation.paper.wv_magnitude_cl61
History: inline check of 2026-10-10, transcribed unchanged except for the period printout.
"""
import pandas as pd

R = "C:/DATA/Projects/202606_E-PROFILE_calibration"


def main():
    for key in ("0-20000-0-06610_C", "0-20000-0-10393_C"):
        for m in ("rayleigh", "cloud"):
            try:
                a = pd.read_csv(f"{R}/diag_v22_l55s008/{key}/{key}_cal.csv")
                b = pd.read_csv(f"{R}/diag_v22_nowv/{key}/{key}_cal.csv")
            except FileNotFoundError as e:
                print(key, m, "missing", e.filename)
                continue
            a = a[(a.method == m) & a.flag.isin([1, 0.5])][["date", "cal_value"]]
            b = b[(b.method == m) & b.flag.isin([1, 0.5])][["date", "cal_value"]]
            j = a.merge(b, on="date", suffixes=("_wv", "_nowv"))
            if len(j):
                r = j.cal_value_nowv / j.cal_value_wv
                print(f"{key} {m}: {len(j)} common nights; C_L(no WV)/C_L(WV, 910.55/0.188) median "
                      f"{r.median():.3f} (IQR {r.quantile(.25):.3f}-{r.quantile(.75):.3f}); "
                      f"nights {j.date.min()} .. {j.date.max()}")


if __name__ == "__main__":
    main()
