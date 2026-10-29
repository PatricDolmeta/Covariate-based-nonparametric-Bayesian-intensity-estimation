# Replicate study: fixed-cap one-chain fits

Updated 2026-10-16 23:55. Protocol replicate_study: one chain, 60000 sweeps, 20000 burn-in, thinning 10. 

| Scenario | n | Datasets | Flagged | GP | Kernel oracle | Kernel CV | Log-poly 1 | Log-poly 2 | Log-poly 3 | Log-poly 4 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 1d_skn | 1000 | 50 | 0 | 0.0342 (0.0085) | 0.0396 (0.0087) | 0.0423 (0.0090) | 0.6235 (0.0015) | 0.2402 (0.0059) | 0.1405 (0.0064) | 0.0639 (0.0066) |
| 1d_skn | 250 | 50 | 0 | 0.0547 (0.0132) | 0.0637 (0.0143) | 0.0698 (0.0177) | 0.6242 (0.0026) | 0.2401 (0.0104) | 0.1431 (0.0116) | 0.0725 (0.0108) |
| 1d_skn | 50 | 50 | 0 | 0.1453 (0.0406) | 0.1195 (0.0347) | 0.1288 (0.0363) | 0.6249 (0.0057) | 0.2557 (0.0231) | 0.1657 (0.0314) | 0.1144 (0.0319) |
| 1d_skn | 500 | 50 | 0 | 0.0455 (0.0151) | 0.0526 (0.0117) | 0.0587 (0.0135) | 0.6235 (0.0020) | 0.2423 (0.0074) | 0.1429 (0.0104) | 0.0696 (0.0151) |
| 2d_hole | 1000 | 50 | 2 | 0.0523 (0.0060) | 0.0636 (0.0066) | 0.0662 (0.0067) | 0.3163 (0.0010) | 0.2622 (0.0028) | 0.1755 (0.0017) | 0.1509 (0.0023) |
| 2d_hole | 250 | 50 | 0 | 0.0900 (0.0115) | 0.0945 (0.0111) | 0.1017 (0.0138) | 0.3178 (0.0030) | 0.2665 (0.0054) | 0.1832 (0.0050) | 0.1639 (0.0059) |
| 2d_hole | 50 | 50 | 0 | 0.1780 (0.0221) | 0.1456 (0.0214) | 0.1616 (0.0274) | 0.3240 (0.0088) | 0.2838 (0.0142) | 0.2226 (0.0188) | 0.2245 (0.0224) |
| 2d_hole | 500 | 50 | 2 | 0.0694 (0.0084) | 0.0781 (0.0080) | 0.0836 (0.0090) | 0.3170 (0.0017) | 0.2644 (0.0040) | 0.1787 (0.0029) | 0.1558 (0.0044) |

Cells: mean relative L2 error over datasets (between-dataset SD).
Flags (relative RMS posterior-mean MC noise > 1%, or within-chain split R-hat > 1.05 on rho(z) probes/mean intensity) mark datasets for inspection; they are never excluded.
Unresolved technical failures: 0

Tables: error_table.csv, paired_differences_vs_GP.csv, computation_table.csv, hyperparameter_descriptive.csv, dataset_summary.csv, method_errors.csv.
