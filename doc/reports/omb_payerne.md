# OmB spot-check - Payerne L1 vs CAMS (2026-05)

Per-station Observation-minus-Background against the CAMS aerosol forecast, comparing the operational L2 constant and our Kalman best-estimate C_L. 910 nm streams (CL31, CL61) use the Angstrom interpolation (532/1064) plus the reused water-vapour correction; CHM15k (1064 nm) is native. Bias = observation - background, in Mm^-1 sr^-1.

> CAMS is on a 1 deg grid; the nearest point to Payerne sits in an Alpine cell whose surface is ~1400 m ASL, so the comparison starts ~900 m above the Payerne plateau (490 m). This is inherent to CAMS resolution over complex terrain, as in the MATLAB reference. Because the near-surface (where 910 nm water-vapour absorption is largest) is below this floor, the reported WV-correction shift is a **lower bound** on the full-column effect.

## Summary (median bias, Mm^-1 sr^-1)

| Instrument | lambda [nm] | median C_L ours | median C_op | bias ours | bias op | bias ours (WV) | RMS ours |
| --- | --- | --- | --- | --- | --- | --- | --- |
| CHM15k | 1064 | 7.267e+11 | 6.240e+11 | -0.0710 | -0.0605 | nan | 1.4814 |
| CL31 | 910 | 3.434e+07 | 1.000e+08 | -0.0474 | -0.2719 | 0.0421 | 4.9019 |
| CL61 | 911 | 1.237e+00 | 1.000e+00 | -0.1150 | -0.1038 | -0.1029 | 1.7010 |

## CHM15k (1064 nm)

![omb CHM15k](omb_payerne/figs/omb_CHM15k.png)

## CL31 (910 nm)

![omb CL31](omb_payerne/figs/omb_CL31.png)

## CL61 (911 nm)

![omb CL61](omb_payerne/figs/omb_CL61.png)
