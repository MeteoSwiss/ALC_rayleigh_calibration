# Molecular-window methods — multi-site comparison (Payerne, Amsterdam, EDT)

Sampled 35 nights (every 2nd day, Mar+Apr 2026) per instrument. `n_cal` = nights that calibrate through the full pipeline (rel_error<=15%). `CL_CV` = night-to-night lidar-constant scatter (lower = more stable).

inst | type | method | nights | n_cal | CL_CV% | med_R2 | med_tcv | med_rel%
---|---|---|---|---|---|---|---|---
Payerne_CHM15k | CHM15k | main | 518 | 215 | 1454.8 | 0.908 | 0.34 | 4.9
Payerne_CHM15k | CHM15k | improved | 518 | 135 | 1149.2 | 0.975 | 0.17 | 7.3
Payerne_CHM15k | CHM15k | matlab | 518 | 184 | 574.3 | 0.940 | 0.24 | 5.5
Payerne_CHM15k | CHM15k | calipso | 518 | 389 | 1534.2 | 0.650 | 0.33 | 8.1
Payerne_CHM15k | CHM15k | earlinet | 518 | 133 | 38.2 | 0.887 | 0.15 | 8.0
Payerne_CHM15k | CHM15k | optimal | 518 | 148 | 1207.1 | 0.967 | 0.18 | 5.4
Payerne_CHM15k | CHM15k | bellini | 518 | 129 | 54.3 | 0.888 | 0.20 | 4.6
Lindenberg_CHM15k | CHM15k | main | 495 | 258 | 647.7 | 0.901 | 0.33 | 5.1
Lindenberg_CHM15k | CHM15k | improved | 495 | 181 | 243.6 | 0.965 | 0.20 | 6.9
Lindenberg_CHM15k | CHM15k | matlab | 495 | 223 | 627.0 | 0.942 | 0.28 | 4.4
Lindenberg_CHM15k | CHM15k | calipso | 495 | 400 | 297.9 | 0.639 | 0.32 | 7.6
Lindenberg_CHM15k | CHM15k | earlinet | 495 | 158 | 29.7 | 0.903 | 0.15 | 6.9
Lindenberg_CHM15k | CHM15k | optimal | 495 | 182 | 86.6 | 0.963 | 0.18 | 5.0
Lindenberg_CHM15k | CHM15k | bellini | 495 | 168 | 36.9 | 0.815 | 0.28 | 4.2
Aosta_CHM15k | CHM15k | main | 29 | 17 | 13.3 | 0.894 | 0.44 | 3.5
Aosta_CHM15k | CHM15k | improved | 29 | 13 | 10.8 | 0.947 | 0.29 | 6.4
Aosta_CHM15k | CHM15k | matlab | 29 | 13 | 10.3 | 0.947 | 0.36 | 6.4
Aosta_CHM15k | CHM15k | calipso | 29 | 23 | 327.1 | 0.686 | 0.41 | 5.7
Aosta_CHM15k | CHM15k | earlinet | 29 | 10 | 8.8 | 0.842 | 0.24 | 9.3
Aosta_CHM15k | CHM15k | optimal | 29 | 10 | 9.4 | 0.987 | 0.23 | 4.0
Aosta_CHM15k | CHM15k | bellini | 29 | 15 | 27.9 | 0.829 | 0.32 | 8.4
Palaiseau_CHM15k | CHM15k | main | 524 | 292 | 574.2 | 0.896 | 0.33 | 4.6
Palaiseau_CHM15k | CHM15k | improved | 524 | 209 | 616.3 | 0.965 | 0.20 | 6.8
Palaiseau_CHM15k | CHM15k | matlab | 524 | 267 | 571.7 | 0.937 | 0.29 | 4.6
Palaiseau_CHM15k | CHM15k | calipso | 524 | 417 | 623.4 | 0.605 | 0.31 | 7.0
Palaiseau_CHM15k | CHM15k | earlinet | 524 | 167 | 31.9 | 0.895 | 0.16 | 6.5
Palaiseau_CHM15k | CHM15k | optimal | 524 | 189 | 34.4 | 0.961 | 0.19 | 5.0
Palaiseau_CHM15k | CHM15k | bellini | 524 | 155 | 32.3 | 0.830 | 0.27 | 4.4
Granada_CHM15k | CHM15k | main | 465 | 306 | 496.0 | 0.953 | 0.30 | 3.8
Granada_CHM15k | CHM15k | improved | 465 | 257 | 150.5 | 0.988 | 0.17 | 5.9
Granada_CHM15k | CHM15k | matlab | 465 | 288 | 469.3 | 0.971 | 0.25 | 4.9
Granada_CHM15k | CHM15k | calipso | 465 | 437 | 189.3 | 0.730 | 0.37 | 6.2
Granada_CHM15k | CHM15k | earlinet | 465 | 253 | 87.5 | 0.917 | 0.14 | 7.2
Granada_CHM15k | CHM15k | optimal | 465 | 263 | 86.7 | 0.980 | 0.19 | 5.8
Granada_CHM15k | CHM15k | bellini | 465 | 279 | 83.9 | 0.899 | 0.24 | 3.7
Magurele_CHM15k | CHM15k | main | 368 | 191 | 565.6 | 0.956 | 0.23 | 4.4
Magurele_CHM15k | CHM15k | improved | 368 | 151 | 136.5 | 0.985 | 0.15 | 5.1
Magurele_CHM15k | CHM15k | matlab | 368 | 178 | 589.9 | 0.968 | 0.18 | 4.6
Magurele_CHM15k | CHM15k | calipso | 368 | 319 | 375.0 | 0.769 | 0.25 | 7.6
Magurele_CHM15k | CHM15k | earlinet | 368 | 188 | 50.4 | 0.904 | 0.13 | 7.0
Magurele_CHM15k | CHM15k | optimal | 368 | 204 | 47.2 | 0.970 | 0.16 | 4.1
Magurele_CHM15k | CHM15k | bellini | 368 | 169 | 23.7 | 0.893 | 0.20 | 3.9
Bergen_CHM15k | CHM15k | main | 411 | 139 | 600.8 | 0.885 | 0.39 | 5.6
Bergen_CHM15k | CHM15k | improved | 411 | 100 | 358.1 | 0.961 | 0.21 | 6.2
Bergen_CHM15k | CHM15k | matlab | 411 | 130 | 590.2 | 0.898 | 0.37 | 6.2
Bergen_CHM15k | CHM15k | calipso | 411 | 231 | 460.0 | 0.634 | 0.35 | 8.2
Bergen_CHM15k | CHM15k | earlinet | 411 | 74 | 38.8 | 0.913 | 0.14 | 7.5
Bergen_CHM15k | CHM15k | optimal | 411 | 86 | 38.9 | 0.966 | 0.18 | 6.8
Bergen_CHM15k | CHM15k | bellini | 411 | 65 | 38.2 | 0.756 | 0.30 | 4.3
Oslo_CHM15k | CHM15k | main | 400 | 197 | 569.9 | 0.885 | 0.32 | 4.7
Oslo_CHM15k | CHM15k | improved | 400 | 141 | 809.1 | 0.952 | 0.21 | 6.3
Oslo_CHM15k | CHM15k | matlab | 400 | 169 | 541.3 | 0.944 | 0.28 | 5.4
Oslo_CHM15k | CHM15k | calipso | 400 | 279 | 386.6 | 0.647 | 0.30 | 7.1
Oslo_CHM15k | CHM15k | earlinet | 400 | 146 | 43.0 | 0.878 | 0.17 | 7.7
Oslo_CHM15k | CHM15k | optimal | 400 | 164 | 43.5 | 0.963 | 0.19 | 4.1
Oslo_CHM15k | CHM15k | bellini | 400 | 113 | 34.6 | 0.807 | 0.28 | 5.1
Hamburg_CHM15k | CHM15k | main | 495 | 216 | 492.0 | 0.935 | 0.26 | 4.9
Hamburg_CHM15k | CHM15k | improved | 495 | 175 | 294.0 | 0.976 | 0.15 | 6.4
Hamburg_CHM15k | CHM15k | matlab | 495 | 200 | 499.9 | 0.955 | 0.21 | 5.3
Hamburg_CHM15k | CHM15k | calipso | 495 | 375 | 403.4 | 0.667 | 0.27 | 7.8
Hamburg_CHM15k | CHM15k | earlinet | 495 | 172 | 17.5 | 0.906 | 0.13 | 6.1
Hamburg_CHM15k | CHM15k | optimal | 495 | 199 | 17.8 | 0.969 | 0.16 | 4.4
Hamburg_CHM15k | CHM15k | bellini | 495 | 166 | 33.5 | 0.862 | 0.24 | 5.2
Hohenpeiss_CHM15k | CHM15k | main | 482 | 246 | 572.3 | 0.913 | 0.32 | 4.9
Hohenpeiss_CHM15k | CHM15k | improved | 482 | 178 | 121.5 | 0.979 | 0.18 | 5.6
Hohenpeiss_CHM15k | CHM15k | matlab | 482 | 218 | 612.8 | 0.947 | 0.28 | 5.1
Hohenpeiss_CHM15k | CHM15k | calipso | 482 | 388 | 549.3 | 0.672 | 0.33 | 6.7
Hohenpeiss_CHM15k | CHM15k | earlinet | 482 | 170 | 25.5 | 0.921 | 0.13 | 7.3
Hohenpeiss_CHM15k | CHM15k | optimal | 482 | 198 | 54.5 | 0.977 | 0.17 | 4.5
Hohenpeiss_CHM15k | CHM15k | bellini | 482 | 164 | 40.1 | 0.887 | 0.28 | 3.9
Brest-MPL_MPL | Mini-MPL | main | 436 | 268 | 78.2 | 0.998 | 0.10 | 1.7
Brest-MPL_MPL | Mini-MPL | improved | 436 | 245 | 61.3 | 0.999 | 0.11 | 2.3
Brest-MPL_MPL | Mini-MPL | matlab | 436 | 257 | 74.1 | 0.997 | 0.09 | 2.3
Brest-MPL_MPL | Mini-MPL | calipso | 436 | 308 | 71.0 | 0.995 | 0.24 | 3.9
Brest-MPL_MPL | Mini-MPL | earlinet | 436 | 246 | 55.9 | 0.990 | 0.09 | 4.6
Brest-MPL_MPL | Mini-MPL | optimal | 436 | 195 | 45.7 | 0.999 | 0.06 | 1.5
Brest-MPL_MPL | Mini-MPL | bellini | 436 | 224 | 58.0 | 0.993 | 0.10 | 1.2
Toulouse-MPL_MPL | Mini-MPL | main | 436 | 338 | 66.7 | 0.999 | 0.07 | 1.9
Toulouse-MPL_MPL | Mini-MPL | improved | 436 | 305 | 58.3 | 0.999 | 0.06 | 3.0
Toulouse-MPL_MPL | Mini-MPL | matlab | 436 | 322 | 62.7 | 0.998 | 0.06 | 2.3
Toulouse-MPL_MPL | Mini-MPL | calipso | 436 | 369 | 66.7 | 0.996 | 0.08 | 3.8
Toulouse-MPL_MPL | Mini-MPL | earlinet | 436 | 318 | 54.8 | 0.990 | 0.06 | 4.8
Toulouse-MPL_MPL | Mini-MPL | optimal | 436 | 290 | 49.2 | 0.999 | 0.05 | 1.1
Toulouse-MPL_MPL | Mini-MPL | bellini | 436 | 255 | 55.6 | 0.997 | 0.06 | 1.5
Corsica-MPL_MPL | Mini-MPL | main | 431 | 353 | 48.0 | 0.999 | 0.05 | 1.8
Corsica-MPL_MPL | Mini-MPL | improved | 431 | 324 | 41.2 | 1.000 | 0.05 | 2.2
Corsica-MPL_MPL | Mini-MPL | matlab | 431 | 336 | 45.6 | 0.999 | 0.05 | 1.9
Corsica-MPL_MPL | Mini-MPL | calipso | 431 | 370 | 48.8 | 0.999 | 0.06 | 3.3
Corsica-MPL_MPL | Mini-MPL | earlinet | 431 | 328 | 41.5 | 0.994 | 0.05 | 5.4
Corsica-MPL_MPL | Mini-MPL | optimal | 431 | 294 | 31.4 | 0.999 | 0.04 | 1.1
Corsica-MPL_MPL | Mini-MPL | bellini | 431 | 261 | 40.4 | 0.999 | 0.05 | 1.0
SIRTA-MPL_MPL | Mini-MPL | main | 432 | 307 | 53.6 | 0.999 | 0.10 | 1.7
SIRTA-MPL_MPL | Mini-MPL | improved | 432 | 300 | 47.8 | 0.999 | 0.10 | 2.1
SIRTA-MPL_MPL | Mini-MPL | matlab | 432 | 300 | 50.7 | 0.998 | 0.10 | 1.6
SIRTA-MPL_MPL | Mini-MPL | calipso | 432 | 372 | 62.9 | 0.997 | 0.13 | 3.3
SIRTA-MPL_MPL | Mini-MPL | earlinet | 432 | 306 | 46.7 | 0.983 | 0.10 | 4.1
SIRTA-MPL_MPL | Mini-MPL | optimal | 432 | 281 | 40.1 | 0.999 | 0.08 | 1.1
SIRTA-MPL_MPL | Mini-MPL | bellini | 432 | 260 | 48.5 | 0.997 | 0.11 | 0.7

## Ranking (mean over instruments)

method | calibrated-fraction | mean CL_CV% | mean temporal_cv | score
---|---|---|---|---
earlinet | 0.45 | 40.7 | 0.13 | 0.196
bellini | 0.42 | 43.4 | 0.21 | 0.160
optimal | 0.45 | 128.0 | 0.15 | 0.062
calipso | 0.79 | 385.4 | 0.27 | 0.032
improved | 0.46 | 292.7 | 0.16 | 0.027
matlab | 0.52 | 380.0 | 0.22 | 0.022
main | 0.57 | 445.2 | 0.25 | 0.020

**Best overall: `earlinet`** (highest usable-night fraction at competitive night-to-night stability).
