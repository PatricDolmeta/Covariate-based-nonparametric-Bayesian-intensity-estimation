# One-chain cap validation against long two-chain fits

**NOT VALIDATED: 36 of 40 comparisons passing (1 missing, 3 failed).**

Generated 2026-10-16 10:50 from C:\Users\patri\Dropbox\Point Processes with Covariates\Data Analysis\Simulation Study\Results\Replicate_Study and C:\Users\patri\Dropbox\Point Processes with Covariates\Data Analysis\Simulation Study\Results\Long_Fits. No sampling.

Tolerances: mean difference < 1% of the truth norm and absolute relative-L2 error change < 0.002.
Expected: 40 comparisons (2 scenarios x 4 sample sizes x datasets 1-5).

Missing comparisons:
- 2d_hole n=1000 r=3: long reference fit not finished

Failed comparisons:
- 1d_skn n=50 r=4
- 2d_hole n=50 r=3
- 2d_hole n=500 r=5

| Scenario | n | Dataset | Reference | Ref. error | Replicate-study error | Mean diff. (%) | Error change | Pass |
|---|---:|---:|---|---:|---:|---:|---:|---|
| 1d_skn | 50 | 1 | 2 chains x 100000 sweeps | 0.16041 | 0.16140 | 0.224 | 0.00099 | TRUE |
| 1d_skn | 50 | 2 | 2 chains x 100000 sweeps | 0.10757 | 0.10749 | 0.263 | 0.00008 | TRUE |
| 1d_skn | 50 | 3 | 2 chains x 100000 sweeps | 0.13537 | 0.13598 | 0.170 | 0.00061 | TRUE |
| 1d_skn | 50 | 4 | 2 chains x 100000 sweeps | 0.17411 | 0.17879 | 0.594 | 0.00468 | FALSE |
| 1d_skn | 50 | 5 | 2 chains x 100000 sweeps | 0.22040 | 0.22105 | 0.244 | 0.00066 | TRUE |
| 1d_skn | 250 | 1 | 2 chains x 100000 sweeps | 0.07251 | 0.07287 | 0.204 | 0.00036 | TRUE |
| 1d_skn | 250 | 2 | 2 chains x 100000 sweeps | 0.05741 | 0.05675 | 0.189 | 0.00066 | TRUE |
| 1d_skn | 250 | 3 | 2 chains x 100000 sweeps | 0.04597 | 0.04607 | 0.080 | 0.00010 | TRUE |
| 1d_skn | 250 | 4 | 2 chains x 100000 sweeps | 0.05265 | 0.05304 | 0.097 | 0.00039 | TRUE |
| 1d_skn | 250 | 5 | 2 chains x 100000 sweeps | 0.04931 | 0.04984 | 0.141 | 0.00052 | TRUE |
| 1d_skn | 500 | 1 | 2 chains x 150000 sweeps | 0.03539 | 0.03568 | 0.121 | 0.00029 | TRUE |
| 1d_skn | 500 | 2 | 2 chains x 100000 sweeps | 0.03858 | 0.03890 | 0.118 | 0.00032 | TRUE |
| 1d_skn | 500 | 3 | 2 chains x 100000 sweeps | 0.04571 | 0.04571 | 0.133 | 0.00001 | TRUE |
| 1d_skn | 500 | 4 | 2 chains x 250000 sweeps | 0.05090 | 0.05015 | 0.156 | 0.00074 | TRUE |
| 1d_skn | 500 | 5 | 2 chains x 150000 sweeps | 0.04253 | 0.04223 | 0.163 | 0.00031 | TRUE |
| 1d_skn | 1000 | 1 | 2 chains x 200000 sweeps | 0.04256 | 0.04285 | 0.144 | 0.00029 | TRUE |
| 1d_skn | 1000 | 2 | 2 chains x 200000 sweeps | 0.04318 | 0.04366 | 0.059 | 0.00048 | TRUE |
| 1d_skn | 1000 | 3 | 2 chains x 200000 sweeps | 0.03283 | 0.03232 | 0.134 | 0.00051 | TRUE |
| 1d_skn | 1000 | 4 | 2 chains x 250000 sweeps | 0.02785 | 0.02857 | 0.190 | 0.00071 | TRUE |
| 1d_skn | 1000 | 5 | 2 chains x 150000 sweeps | 0.03748 | 0.03820 | 0.131 | 0.00072 | TRUE |
| 2d_hole | 50 | 1 | 2 chains x 200000 sweeps | 0.18373 | 0.18406 | 0.192 | 0.00032 | TRUE |
| 2d_hole | 50 | 2 | 2 chains x 200000 sweeps | 0.19384 | 0.19374 | 0.437 | 0.00010 | TRUE |
| 2d_hole | 50 | 3 | 2 chains x 200000 sweeps | 0.17825 | 0.18087 | 0.520 | 0.00262 | FALSE |
| 2d_hole | 50 | 4 | 2 chains x 200000 sweeps | 0.20034 | 0.19984 | 0.219 | 0.00050 | TRUE |
| 2d_hole | 50 | 5 | 2 chains x 200000 sweeps | 0.17693 | 0.17771 | 0.365 | 0.00078 | TRUE |
| 2d_hole | 250 | 1 | 2 chains x 200000 sweeps | 0.07969 | 0.07981 | 0.349 | 0.00012 | TRUE |
| 2d_hole | 250 | 2 | 2 chains x 200000 sweeps | 0.10600 | 0.10554 | 0.243 | 0.00046 | TRUE |
| 2d_hole | 250 | 3 | 2 chains x 200000 sweeps | 0.09036 | 0.09010 | 0.208 | 0.00027 | TRUE |
| 2d_hole | 250 | 4 | 2 chains x 250000 sweeps | 0.09584 | 0.09581 | 0.218 | 0.00003 | TRUE |
| 2d_hole | 250 | 5 | 2 chains x 200000 sweeps | 0.08813 | 0.08867 | 0.321 | 0.00054 | TRUE |
| 2d_hole | 500 | 1 | 2 chains x 250000 sweeps | 0.05713 | 0.05724 | 0.277 | 0.00011 | TRUE |
| 2d_hole | 500 | 2 | 2 chains x 200000 sweeps | 0.06603 | 0.06621 | 0.178 | 0.00017 | TRUE |
| 2d_hole | 500 | 3 | 2 chains x 200000 sweeps | 0.05911 | 0.06082 | 0.381 | 0.00171 | TRUE |
| 2d_hole | 500 | 4 | 2 chains x 200000 sweeps | 0.07837 | 0.07895 | 0.308 | 0.00057 | TRUE |
| 2d_hole | 500 | 5 | 2 chains x 250000 sweeps | 0.06361 | 0.06593 | 0.658 | 0.00233 | FALSE |
| 2d_hole | 1000 | 1 | 2 chains x 300000 sweeps | 0.05923 | 0.05941 | 0.381 | 0.00018 | TRUE |
| 2d_hole | 1000 | 2 | 2 chains x 300000 sweeps | 0.06177 | 0.06231 | 0.247 | 0.00054 | TRUE |
| 2d_hole | 1000 | 4 | 2 chains x 250000 sweeps | 0.05784 | 0.05778 | 0.261 | 0.00007 | TRUE |
| 2d_hole | 1000 | 5 | 2 chains x 250000 sweeps | 0.05112 | 0.05107 | 0.346 | 0.00005 | TRUE |

Datasets compared: 39; passing: 36; largest mean difference 0.658%; largest error change 0.00468.
