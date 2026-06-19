# Long-run (full archive) — molecular-window methods, ROBUST aggregation

14 instruments (10 CHM15k + 4 Mini-MPL), 5 nights/month over the full E-PROFILE L2 archive (~80-113 months each). `rob_CV` = robust night-to-night scatter (1.4826·MAD/median, %); `std_CV` = classic std/mean (outlier-sensitive, shown for reference). n_cal = nights calibrating through the pipeline.

inst | type | method | nights | n_cal | rob_CV% | std_CV% | med_R2 | med_tcv
---|---|---|---|---|---|---|---|---
Aosta_CHM15k | CHM15k | eprof_v10 | 29 | 13 | 16.1 | 15 | - | -
Aosta_CHM15k | CHM15k | main | 29 | 17 | 14.8 | 13 | 0.894 | 0.44
Aosta_CHM15k | CHM15k | improved | 29 | 13 | 14.7 | 11 | 0.947 | 0.29
Aosta_CHM15k | CHM15k | optimal | 29 | 10 | 10.2 | 9 | 0.987 | 0.23
Aosta_CHM15k | CHM15k | matlab | 29 | 13 | 13.7 | 10 | 0.947 | 0.36
Aosta_CHM15k | CHM15k | calipso | 29 | 23 | 21.3 | 327 | 0.686 | 0.41
Aosta_CHM15k | CHM15k | earlinet | 29 | 10 | 8.3 | 9 | 0.842 | 0.24
Aosta_CHM15k | CHM15k | bellini | 29 | 15 | 15.4 | 28 | 0.829 | 0.32
Bergen_CHM15k | CHM15k | eprof_v10 | 411 | 89 | 35.9 | 596 | - | -
Bergen_CHM15k | CHM15k | main | 411 | 139 | 48.6 | 592 | 0.885 | 0.39
Bergen_CHM15k | CHM15k | improved | 411 | 100 | 48.8 | 358 | 0.961 | 0.21
Bergen_CHM15k | CHM15k | optimal | 411 | 86 | 53.1 | 39 | 0.966 | 0.18
Bergen_CHM15k | CHM15k | matlab | 411 | 130 | 50.3 | 590 | 0.898 | 0.37
Bergen_CHM15k | CHM15k | calipso | 411 | 231 | 49.3 | 460 | 0.634 | 0.35
Bergen_CHM15k | CHM15k | earlinet | 411 | 74 | 53.6 | 39 | 0.913 | 0.14
Bergen_CHM15k | CHM15k | bellini | 411 | 65 | 35.0 | 38 | 0.756 | 0.30
Brest-MPL_MPL | Mini-MPL | eprof_v10 | 436 | 197 | 64.5 | 66 | - | -
Brest-MPL_MPL | Mini-MPL | main | 436 | 268 | 75.0 | 77 | 0.998 | 0.10
Brest-MPL_MPL | Mini-MPL | improved | 436 | 245 | 77.9 | 61 | 0.999 | 0.11
Brest-MPL_MPL | Mini-MPL | optimal | 436 | 195 | 63.3 | 46 | 0.999 | 0.06
Brest-MPL_MPL | Mini-MPL | matlab | 436 | 257 | 72.0 | 74 | 0.997 | 0.09
Brest-MPL_MPL | Mini-MPL | calipso | 436 | 308 | 90.3 | 71 | 0.995 | 0.24
Brest-MPL_MPL | Mini-MPL | earlinet | 436 | 246 | 68.0 | 56 | 0.990 | 0.09
Brest-MPL_MPL | Mini-MPL | bellini | 436 | 224 | 64.7 | 58 | 0.993 | 0.10
Corsica-MPL_MPL | Mini-MPL | eprof_v10 | 431 | 317 | 35.4 | 41 | - | -
Corsica-MPL_MPL | Mini-MPL | main | 431 | 353 | 37.7 | 47 | 0.999 | 0.05
Corsica-MPL_MPL | Mini-MPL | improved | 431 | 324 | 34.6 | 41 | 1.000 | 0.05
Corsica-MPL_MPL | Mini-MPL | optimal | 431 | 294 | 30.3 | 31 | 0.999 | 0.04
Corsica-MPL_MPL | Mini-MPL | matlab | 431 | 336 | 36.6 | 46 | 0.999 | 0.05
Corsica-MPL_MPL | Mini-MPL | calipso | 431 | 370 | 42.0 | 49 | 0.999 | 0.06
Corsica-MPL_MPL | Mini-MPL | earlinet | 431 | 328 | 33.4 | 42 | 0.994 | 0.05
Corsica-MPL_MPL | Mini-MPL | bellini | 431 | 261 | 34.1 | 40 | 0.999 | 0.05
Granada_CHM15k | CHM15k | eprof_v10 | 465 | 275 | 17.6 | 123 | - | -
Granada_CHM15k | CHM15k | main | 465 | 306 | 16.5 | 496 | 0.953 | 0.30
Granada_CHM15k | CHM15k | improved | 465 | 257 | 14.8 | 151 | 0.988 | 0.17
Granada_CHM15k | CHM15k | optimal | 465 | 263 | 13.1 | 87 | 0.980 | 0.19
Granada_CHM15k | CHM15k | matlab | 465 | 288 | 15.9 | 469 | 0.971 | 0.25
Granada_CHM15k | CHM15k | calipso | 465 | 437 | 19.7 | 189 | 0.730 | 0.37
Granada_CHM15k | CHM15k | earlinet | 465 | 253 | 13.2 | 87 | 0.917 | 0.14
Granada_CHM15k | CHM15k | bellini | 465 | 279 | 15.2 | 84 | 0.899 | 0.24
Hamburg_CHM15k | CHM15k | eprof_v10 | 495 | 168 | 22.9 | 486 | - | -
Hamburg_CHM15k | CHM15k | main | 495 | 216 | 24.4 | 492 | 0.935 | 0.26
Hamburg_CHM15k | CHM15k | improved | 495 | 175 | 21.7 | 294 | 0.976 | 0.15
Hamburg_CHM15k | CHM15k | optimal | 495 | 199 | 16.3 | 18 | 0.969 | 0.16
Hamburg_CHM15k | CHM15k | matlab | 495 | 200 | 22.8 | 500 | 0.955 | 0.21
Hamburg_CHM15k | CHM15k | calipso | 495 | 375 | 23.6 | 403 | 0.667 | 0.27
Hamburg_CHM15k | CHM15k | earlinet | 495 | 172 | 14.9 | 17 | 0.906 | 0.13
Hamburg_CHM15k | CHM15k | bellini | 495 | 166 | 24.7 | 33 | 0.862 | 0.24
Hohenpeiss_CHM15k | CHM15k | eprof_v10 | 482 | 193 | 28.5 | 572 | - | -
Hohenpeiss_CHM15k | CHM15k | main | 482 | 246 | 29.0 | 572 | 0.913 | 0.32
Hohenpeiss_CHM15k | CHM15k | improved | 482 | 178 | 32.3 | 121 | 0.979 | 0.18
Hohenpeiss_CHM15k | CHM15k | optimal | 482 | 198 | 26.5 | 54 | 0.977 | 0.17
Hohenpeiss_CHM15k | CHM15k | matlab | 482 | 218 | 27.9 | 613 | 0.947 | 0.28
Hohenpeiss_CHM15k | CHM15k | calipso | 482 | 388 | 35.5 | 549 | 0.672 | 0.33
Hohenpeiss_CHM15k | CHM15k | earlinet | 482 | 170 | 21.7 | 25 | 0.921 | 0.13
Hohenpeiss_CHM15k | CHM15k | bellini | 482 | 164 | 26.6 | 40 | 0.887 | 0.28
Lindenberg_CHM15k | CHM15k | eprof_v10 | 495 | 209 | 29.9 | 681 | - | -
Lindenberg_CHM15k | CHM15k | main | 495 | 258 | 31.7 | 648 | 0.901 | 0.33
Lindenberg_CHM15k | CHM15k | improved | 495 | 181 | 28.3 | 244 | 0.965 | 0.20
Lindenberg_CHM15k | CHM15k | optimal | 495 | 182 | 29.9 | 87 | 0.963 | 0.18
Lindenberg_CHM15k | CHM15k | matlab | 495 | 223 | 32.3 | 627 | 0.942 | 0.28
Lindenberg_CHM15k | CHM15k | calipso | 495 | 400 | 37.9 | 298 | 0.639 | 0.32
Lindenberg_CHM15k | CHM15k | earlinet | 495 | 158 | 30.0 | 30 | 0.903 | 0.15
Lindenberg_CHM15k | CHM15k | bellini | 495 | 168 | 32.6 | 37 | 0.815 | 0.28
Magurele_CHM15k | CHM15k | eprof_v10 | 368 | 172 | 22.8 | 590 | - | -
Magurele_CHM15k | CHM15k | main | 368 | 191 | 24.5 | 566 | 0.956 | 0.23
Magurele_CHM15k | CHM15k | improved | 368 | 151 | 28.9 | 137 | 0.985 | 0.15
Magurele_CHM15k | CHM15k | optimal | 368 | 204 | 21.7 | 47 | 0.970 | 0.16
Magurele_CHM15k | CHM15k | matlab | 368 | 178 | 24.9 | 590 | 0.968 | 0.18
Magurele_CHM15k | CHM15k | calipso | 368 | 319 | 28.7 | 375 | 0.769 | 0.25
Magurele_CHM15k | CHM15k | earlinet | 368 | 188 | 22.1 | 50 | 0.904 | 0.13
Magurele_CHM15k | CHM15k | bellini | 368 | 169 | 22.2 | 24 | 0.893 | 0.20
Oslo_CHM15k | CHM15k | eprof_v10 | 400 | 150 | 28.6 | 623 | - | -
Oslo_CHM15k | CHM15k | main | 400 | 197 | 40.5 | 563 | 0.885 | 0.32
Oslo_CHM15k | CHM15k | improved | 400 | 141 | 34.7 | 809 | 0.952 | 0.21
Oslo_CHM15k | CHM15k | optimal | 400 | 164 | 21.3 | 43 | 0.963 | 0.19
Oslo_CHM15k | CHM15k | matlab | 400 | 169 | 27.9 | 541 | 0.944 | 0.28
Oslo_CHM15k | CHM15k | calipso | 400 | 279 | 40.0 | 387 | 0.647 | 0.30
Oslo_CHM15k | CHM15k | earlinet | 400 | 146 | 18.9 | 43 | 0.878 | 0.17
Oslo_CHM15k | CHM15k | bellini | 400 | 113 | 14.9 | 35 | 0.807 | 0.28
Palaiseau_CHM15k | CHM15k | eprof_v10 | 524 | 232 | 23.0 | 605 | - | -
Palaiseau_CHM15k | CHM15k | main | 524 | 292 | 25.8 | 574 | 0.896 | 0.33
Palaiseau_CHM15k | CHM15k | improved | 524 | 209 | 26.7 | 616 | 0.965 | 0.20
Palaiseau_CHM15k | CHM15k | optimal | 524 | 189 | 22.0 | 34 | 0.961 | 0.19
Palaiseau_CHM15k | CHM15k | matlab | 524 | 267 | 22.7 | 572 | 0.937 | 0.29
Palaiseau_CHM15k | CHM15k | calipso | 524 | 417 | 26.5 | 623 | 0.605 | 0.31
Palaiseau_CHM15k | CHM15k | earlinet | 524 | 167 | 22.0 | 32 | 0.895 | 0.16
Palaiseau_CHM15k | CHM15k | bellini | 524 | 155 | 23.8 | 32 | 0.830 | 0.27
Payerne_CHM15k | CHM15k | eprof_v10 | 518 | 163 | 43.6 | 1267 | - | -
Payerne_CHM15k | CHM15k | main | 518 | 215 | 47.6 | 1455 | 0.908 | 0.34
Payerne_CHM15k | CHM15k | improved | 518 | 135 | 56.5 | 1149 | 0.975 | 0.17
Payerne_CHM15k | CHM15k | optimal | 518 | 148 | 40.3 | 1207 | 0.967 | 0.18
Payerne_CHM15k | CHM15k | matlab | 518 | 184 | 44.2 | 574 | 0.940 | 0.24
Payerne_CHM15k | CHM15k | calipso | 518 | 389 | 62.5 | 1534 | 0.650 | 0.33
Payerne_CHM15k | CHM15k | earlinet | 518 | 133 | 41.2 | 38 | 0.887 | 0.15
Payerne_CHM15k | CHM15k | bellini | 518 | 129 | 42.3 | 54 | 0.888 | 0.20
SIRTA-MPL_MPL | Mini-MPL | eprof_v10 | 432 | 274 | 49.0 | 46 | - | -
SIRTA-MPL_MPL | Mini-MPL | main | 432 | 307 | 51.8 | 53 | 0.999 | 0.10
SIRTA-MPL_MPL | Mini-MPL | improved | 432 | 300 | 53.0 | 48 | 0.999 | 0.10
SIRTA-MPL_MPL | Mini-MPL | optimal | 432 | 281 | 45.7 | 40 | 0.999 | 0.08
SIRTA-MPL_MPL | Mini-MPL | matlab | 432 | 300 | 49.8 | 51 | 0.998 | 0.10
SIRTA-MPL_MPL | Mini-MPL | calipso | 432 | 372 | 61.4 | 63 | 0.997 | 0.13
SIRTA-MPL_MPL | Mini-MPL | earlinet | 432 | 306 | 52.0 | 47 | 0.983 | 0.10
SIRTA-MPL_MPL | Mini-MPL | bellini | 432 | 260 | 53.8 | 49 | 0.997 | 0.11
Toulouse-MPL_MPL | Mini-MPL | eprof_v10 | 436 | 285 | 47.5 | 61 | - | -
Toulouse-MPL_MPL | Mini-MPL | main | 436 | 338 | 56.0 | 66 | 0.999 | 0.07
Toulouse-MPL_MPL | Mini-MPL | improved | 436 | 305 | 54.1 | 58 | 0.999 | 0.06
Toulouse-MPL_MPL | Mini-MPL | optimal | 436 | 290 | 44.1 | 49 | 0.999 | 0.05
Toulouse-MPL_MPL | Mini-MPL | matlab | 436 | 322 | 52.7 | 63 | 0.998 | 0.06
Toulouse-MPL_MPL | Mini-MPL | calipso | 436 | 369 | 56.2 | 67 | 0.996 | 0.08
Toulouse-MPL_MPL | Mini-MPL | earlinet | 436 | 318 | 50.3 | 55 | 0.990 | 0.06
Toulouse-MPL_MPL | Mini-MPL | bellini | 436 | 255 | 53.5 | 56 | 0.997 | 0.06

## Ranking (mean over instruments, ROBUST CV)

method | calibrated-fraction | mean rob_CV% | mean temporal_cv | score
---|---|---|---|---
calipso | 0.79 | 42.5 | 0.27 | 0.293
optimal | 0.45 | 31.3 | 0.15 | 0.253
earlinet | 0.45 | 32.1 | 0.13 | 0.248
main | 0.57 | 37.4 | 0.25 | 0.242
matlab | 0.52 | 35.2 | 0.22 | 0.242
bellini | 0.42 | 32.8 | 0.21 | 0.212
improved | 0.46 | 37.7 | 0.16 | 0.211
eprof_v10 | 0.46 | 33.2 | - | 0.140

**Best overall (robust, full archive): `calipso`.**
