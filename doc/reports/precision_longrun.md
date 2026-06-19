# Long-run precision — drift-insensitive metrics (full archive, 14 sites)

CV mixes precision with seasonal + long-term drift; the metrics below remove slow drift. All scatters are % of the median CL. `sigma_SD` (successive-difference), `sigma_dt` (detrended), `sigma_im` (within-month) are the precision metrics; CV is kept only for reference.

## Per method (mean over the 14 instruments)

method | valid_% | sigma_night% | **sigma_SD%** | sigma_dt% | sigma_im% | CV% (ref)
---|---|---|---|---|---|---
optimal | 45 | 8.6 | **13.6** | 12.5 | 8.4 | 128
earlinet | 45 | 4.2 | **15.9** | 13.6 | 10.0 | 41
eprof_v10 | 46 | 4.6 | **16.5** | 13.4 | 10.3 | 412
bellini | 42 | 10.7 | **18.3** | 15.4 | 11.0 | 43
matlab | 52 | 16.7 | **18.5** | 15.3 | 12.2 | 380
improved | 46 | 9.4 | **19.2** | 16.1 | 12.3 | 293
main | 57 | 20.0 | **20.4** | 16.0 | 12.1 | 444
calipso | 79 | 14.9 | **24.9** | 19.9 | 16.4 | 385

**Most precise (lowest sigma_SD): `optimal`.** Note CV ≫ sigma_SD for every method — most of the CV was real drift, not noise.

## Per instrument × method

inst | type | method | valid_% | sigma_night% | sigma_SD% | sigma_dt% | sigma_im% | CV%
---|---|---|---|---|---|---|---|---
Aosta_CHM15k | CHM15k | eprof_v10 | 45 | 4.8 | 10.6 | 10.0 | 10.2 | 15
Aosta_CHM15k | CHM15k | main | 59 | 28.3 | 14.9 | 9.9 | 4.9 | 13
Aosta_CHM15k | CHM15k | improved | 45 | 18.4 | 9.9 | 13.7 | 10.0 | 11
Aosta_CHM15k | CHM15k | optimal | 34 | 9.7 | 8.6 | 9.6 | 7.2 | 9
Aosta_CHM15k | CHM15k | matlab | 45 | 19.1 | 10.1 | 8.1 | 14.4 | 10
Aosta_CHM15k | CHM15k | calipso | 79 | 24.4 | 11.7 | 16.9 | 14.8 | 327
Aosta_CHM15k | CHM15k | earlinet | 34 | 5.7 | 12.0 | 5.2 | 4.9 | 9
Aosta_CHM15k | CHM15k | bellini | 52 | 13.5 | 14.4 | 6.0 | 4.3 | 28
Bergen_CHM15k | CHM15k | eprof_v10 | 22 | 6.1 | 14.3 | 11.9 | 8.7 | 596
Bergen_CHM15k | CHM15k | main | 33 | 34.5 | 19.8 | 13.8 | 9.0 | 592
Bergen_CHM15k | CHM15k | improved | 24 | 14.0 | 12.9 | 11.9 | 8.2 | 358
Bergen_CHM15k | CHM15k | optimal | 21 | 10.9 | 9.3 | 9.8 | 5.8 | 39
Bergen_CHM15k | CHM15k | matlab | 32 | 33.2 | 12.8 | 10.0 | 8.3 | 590
Bergen_CHM15k | CHM15k | calipso | 56 | 21.2 | 22.9 | 14.9 | 13.0 | 460
Bergen_CHM15k | CHM15k | earlinet | 18 | 5.7 | 11.8 | 8.7 | 5.3 | 39
Bergen_CHM15k | CHM15k | bellini | 16 | 16.9 | 17.8 | 12.9 | 6.6 | 38
Brest-MPL_MPL | Mini-MPL | eprof_v10 | 45 | 3.5 | 28.0 | 22.2 | 16.9 | 66
Brest-MPL_MPL | Mini-MPL | main | 61 | 2.9 | 46.9 | 37.9 | 27.7 | 77
Brest-MPL_MPL | Mini-MPL | improved | 56 | 2.6 | 41.9 | 38.3 | 26.7 | 61
Brest-MPL_MPL | Mini-MPL | optimal | 45 | 2.6 | 22.2 | 25.5 | 14.4 | 46
Brest-MPL_MPL | Mini-MPL | matlab | 59 | 3.0 | 38.5 | 33.4 | 24.1 | 74
Brest-MPL_MPL | Mini-MPL | calipso | 71 | 4.0 | 56.4 | 50.5 | 39.2 | 71
Brest-MPL_MPL | Mini-MPL | earlinet | 56 | 1.6 | 34.3 | 31.9 | 24.1 | 56
Brest-MPL_MPL | Mini-MPL | bellini | 51 | 2.2 | 39.8 | 34.2 | 26.2 | 58
Corsica-MPL_MPL | Mini-MPL | eprof_v10 | 74 | 1.2 | 17.8 | 14.5 | 10.4 | 41
Corsica-MPL_MPL | Mini-MPL | main | 81 | 1.3 | 23.0 | 16.9 | 13.1 | 47
Corsica-MPL_MPL | Mini-MPL | improved | 75 | 1.2 | 20.0 | 14.9 | 14.1 | 41
Corsica-MPL_MPL | Mini-MPL | optimal | 68 | 1.7 | 15.3 | 13.0 | 9.4 | 31
Corsica-MPL_MPL | Mini-MPL | matlab | 78 | 1.3 | 23.7 | 18.1 | 14.3 | 46
Corsica-MPL_MPL | Mini-MPL | calipso | 86 | 2.1 | 21.7 | 19.8 | 13.4 | 49
Corsica-MPL_MPL | Mini-MPL | earlinet | 76 | 1.2 | 18.7 | 14.8 | 11.9 | 42
Corsica-MPL_MPL | Mini-MPL | bellini | 61 | 0.7 | 18.9 | 15.2 | 11.8 | 40
Granada_CHM15k | CHM15k | eprof_v10 | 59 | 5.0 | 11.3 | 10.1 | 7.6 | 123
Granada_CHM15k | CHM15k | main | 66 | 20.0 | 12.3 | 10.2 | 8.5 | 496
Granada_CHM15k | CHM15k | improved | 55 | 9.9 | 13.9 | 10.9 | 9.2 | 151
Granada_CHM15k | CHM15k | optimal | 57 | 11.3 | 10.9 | 9.4 | 6.2 | 87
Granada_CHM15k | CHM15k | matlab | 62 | 15.9 | 12.4 | 10.6 | 8.7 | 469
Granada_CHM15k | CHM15k | calipso | 94 | 19.1 | 17.9 | 13.4 | 11.7 | 189
Granada_CHM15k | CHM15k | earlinet | 54 | 4.8 | 11.1 | 9.0 | 7.1 | 87
Granada_CHM15k | CHM15k | bellini | 60 | 12.4 | 12.6 | 11.1 | 8.0 | 84
Hamburg_CHM15k | CHM15k | eprof_v10 | 34 | 5.7 | 13.3 | 11.3 | 7.6 | 486
Hamburg_CHM15k | CHM15k | main | 44 | 24.1 | 15.0 | 11.8 | 8.4 | 492
Hamburg_CHM15k | CHM15k | improved | 35 | 12.4 | 18.7 | 12.8 | 10.5 | 294
Hamburg_CHM15k | CHM15k | optimal | 40 | 11.2 | 10.4 | 9.6 | 6.6 | 18
Hamburg_CHM15k | CHM15k | matlab | 40 | 20.6 | 15.0 | 12.4 | 8.4 | 500
Hamburg_CHM15k | CHM15k | calipso | 76 | 18.5 | 18.2 | 13.6 | 10.4 | 403
Hamburg_CHM15k | CHM15k | earlinet | 35 | 4.9 | 11.1 | 9.5 | 7.1 | 17
Hamburg_CHM15k | CHM15k | bellini | 34 | 13.2 | 16.3 | 12.9 | 10.2 | 33
Hohenpeiss_CHM15k | CHM15k | eprof_v10 | 40 | 5.1 | 15.3 | 12.6 | 8.7 | 572
Hohenpeiss_CHM15k | CHM15k | main | 51 | 25.4 | 13.5 | 12.2 | 8.1 | 572
Hohenpeiss_CHM15k | CHM15k | improved | 37 | 10.9 | 12.9 | 8.3 | 7.1 | 121
Hohenpeiss_CHM15k | CHM15k | optimal | 41 | 11.2 | 9.9 | 8.2 | 5.2 | 54
Hohenpeiss_CHM15k | CHM15k | matlab | 45 | 21.0 | 14.1 | 10.1 | 7.2 | 613
Hohenpeiss_CHM15k | CHM15k | calipso | 80 | 19.0 | 22.3 | 13.5 | 13.2 | 549
Hohenpeiss_CHM15k | CHM15k | earlinet | 35 | 4.5 | 9.2 | 8.1 | 5.1 | 25
Hohenpeiss_CHM15k | CHM15k | bellini | 34 | 14.4 | 10.8 | 8.7 | 7.3 | 40
Lindenberg_CHM15k | CHM15k | eprof_v10 | 42 | 5.9 | 11.2 | 9.0 | 7.6 | 681
Lindenberg_CHM15k | CHM15k | main | 52 | 28.4 | 14.2 | 10.8 | 8.9 | 648
Lindenberg_CHM15k | CHM15k | improved | 37 | 12.8 | 15.5 | 12.1 | 8.6 | 244
Lindenberg_CHM15k | CHM15k | optimal | 37 | 11.3 | 11.0 | 8.1 | 5.0 | 87
Lindenberg_CHM15k | CHM15k | matlab | 45 | 25.9 | 14.3 | 11.8 | 9.6 | 627
Lindenberg_CHM15k | CHM15k | calipso | 81 | 20.6 | 22.6 | 14.8 | 13.6 | 298
Lindenberg_CHM15k | CHM15k | earlinet | 32 | 5.6 | 10.8 | 8.7 | 5.6 | 30
Lindenberg_CHM15k | CHM15k | bellini | 34 | 16.0 | 14.7 | 13.3 | 9.3 | 37
Magurele_CHM15k | CHM15k | eprof_v10 | 47 | 5.1 | 14.7 | 10.3 | 9.1 | 590
Magurele_CHM15k | CHM15k | main | 52 | 16.4 | 16.5 | 10.3 | 9.2 | 566
Magurele_CHM15k | CHM15k | improved | 41 | 8.3 | 19.9 | 14.9 | 10.5 | 137
Magurele_CHM15k | CHM15k | optimal | 55 | 10.6 | 11.9 | 10.6 | 7.3 | 47
Magurele_CHM15k | CHM15k | matlab | 48 | 13.1 | 15.6 | 11.8 | 9.6 | 590
Magurele_CHM15k | CHM15k | calipso | 87 | 14.7 | 18.8 | 14.3 | 11.1 | 375
Magurele_CHM15k | CHM15k | earlinet | 51 | 4.7 | 12.2 | 11.1 | 7.6 | 50
Magurele_CHM15k | CHM15k | bellini | 46 | 11.4 | 12.6 | 8.8 | 6.9 | 24
Oslo_CHM15k | CHM15k | eprof_v10 | 38 | 5.1 | 12.8 | 10.1 | 7.6 | 623
Oslo_CHM15k | CHM15k | main | 48 | 31.1 | 14.6 | 10.4 | 6.9 | 563
Oslo_CHM15k | CHM15k | improved | 35 | 13.2 | 9.3 | 8.0 | 6.4 | 809
Oslo_CHM15k | CHM15k | optimal | 41 | 12.0 | 8.2 | 9.0 | 4.0 | 43
Oslo_CHM15k | CHM15k | matlab | 42 | 26.5 | 13.5 | 10.0 | 6.0 | 541
Oslo_CHM15k | CHM15k | calipso | 70 | 18.7 | 15.5 | 11.6 | 9.5 | 387
Oslo_CHM15k | CHM15k | earlinet | 36 | 6.2 | 8.5 | 7.7 | 4.4 | 43
Oslo_CHM15k | CHM15k | bellini | 28 | 18.3 | 13.1 | 11.8 | 6.9 | 35
Palaiseau_CHM15k | CHM15k | eprof_v10 | 44 | 7.1 | 11.1 | 9.0 | 6.7 | 605
Palaiseau_CHM15k | CHM15k | main | 56 | 31.1 | 13.1 | 12.1 | 8.0 | 574
Palaiseau_CHM15k | CHM15k | improved | 40 | 14.5 | 13.0 | 10.5 | 7.8 | 616
Palaiseau_CHM15k | CHM15k | optimal | 36 | 11.8 | 12.5 | 9.8 | 8.2 | 34
Palaiseau_CHM15k | CHM15k | matlab | 51 | 26.2 | 13.1 | 11.7 | 9.5 | 572
Palaiseau_CHM15k | CHM15k | calipso | 80 | 21.6 | 18.5 | 13.6 | 11.7 | 623
Palaiseau_CHM15k | CHM15k | earlinet | 32 | 6.2 | 12.9 | 10.9 | 7.9 | 32
Palaiseau_CHM15k | CHM15k | bellini | 30 | 16.8 | 12.8 | 12.7 | 7.9 | 32
Payerne_CHM15k | CHM15k | eprof_v10 | 31 | 5.7 | 17.6 | 15.8 | 10.1 | 1267
Payerne_CHM15k | CHM15k | main | 42 | 31.7 | 18.2 | 12.3 | 9.9 | 1455
Payerne_CHM15k | CHM15k | improved | 26 | 9.1 | 25.7 | 20.4 | 11.5 | 1149
Payerne_CHM15k | CHM15k | optimal | 29 | 11.6 | 13.8 | 12.1 | 9.6 | 1207
Payerne_CHM15k | CHM15k | matlab | 36 | 24.6 | 17.0 | 14.1 | 9.0 | 574
Payerne_CHM15k | CHM15k | calipso | 75 | 19.0 | 31.1 | 20.1 | 17.8 | 1534
Payerne_CHM15k | CHM15k | earlinet | 26 | 5.3 | 12.4 | 11.3 | 8.3 | 38
Payerne_CHM15k | CHM15k | bellini | 25 | 11.9 | 14.8 | 14.1 | 9.0 | 54
SIRTA-MPL_MPL | Mini-MPL | eprof_v10 | 63 | 2.6 | 24.8 | 20.3 | 15.9 | 46
SIRTA-MPL_MPL | Mini-MPL | main | 71 | 1.9 | 27.2 | 26.3 | 22.1 | 53
SIRTA-MPL_MPL | Mini-MPL | improved | 69 | 1.7 | 26.9 | 25.7 | 20.5 | 48
SIRTA-MPL_MPL | Mini-MPL | optimal | 65 | 2.2 | 22.8 | 19.4 | 15.3 | 40
SIRTA-MPL_MPL | Mini-MPL | matlab | 69 | 1.8 | 27.1 | 27.3 | 21.0 | 51
SIRTA-MPL_MPL | Mini-MPL | calipso | 86 | 3.1 | 36.2 | 32.0 | 25.8 | 63
SIRTA-MPL_MPL | Mini-MPL | earlinet | 71 | 1.7 | 27.6 | 22.2 | 20.4 | 47
SIRTA-MPL_MPL | Mini-MPL | bellini | 60 | 1.4 | 28.8 | 27.6 | 20.5 | 49
Toulouse-MPL_MPL | Mini-MPL | eprof_v10 | 65 | 2.1 | 28.5 | 20.9 | 16.5 | 61
Toulouse-MPL_MPL | Mini-MPL | main | 77 | 2.2 | 36.7 | 29.0 | 24.4 | 66
Toulouse-MPL_MPL | Mini-MPL | improved | 70 | 2.2 | 29.0 | 22.8 | 20.5 | 58
Toulouse-MPL_MPL | Mini-MPL | optimal | 67 | 2.3 | 23.9 | 21.2 | 12.8 | 49
Toulouse-MPL_MPL | Mini-MPL | matlab | 74 | 2.3 | 32.1 | 25.1 | 20.5 | 63
Toulouse-MPL_MPL | Mini-MPL | calipso | 85 | 2.9 | 35.1 | 29.7 | 24.8 | 67
Toulouse-MPL_MPL | Mini-MPL | earlinet | 73 | 1.4 | 29.5 | 30.8 | 19.7 | 55
Toulouse-MPL_MPL | Mini-MPL | bellini | 58 | 1.4 | 28.3 | 27.0 | 18.5 | 56
