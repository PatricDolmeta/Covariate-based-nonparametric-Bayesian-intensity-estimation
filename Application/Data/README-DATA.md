# Canadian wildfire data (June, 2004-2022)

Data of the wildfire application (paper Section 4, Supplement Section S5). They were extracted from the Canadian Wildland Fire Information System (CWFIS, http://cwfis.cfs.nrcan.gc.ca/home): for each year from 2004 to 2022 (n = 19), the locations of the wildfire hotspots detected in June, the month of peak activity, together with the June averages of three meteorological covariates on a pixel grid.

| File | Province | Object prefix | Valid pixels | Events (2004-2022) |
|---|---|---|---|---|
| `ontario_data_june_orig_scale.RData` | Ontario (main text) | `ontario` | 21615 | 648 |
| `saska_data_june_orig_scale.RData` | Saskatchewan (Supplement) | `saska` | 14300 | 2610 |
| `bt_columbia_data_june_orig_scale.RData` | British Columbia (Supplement) | `bt_columbia` | 20795 | 994 |

Each file contains, for its prefix, lists of `spatstat` objects named by year:

- `<prefix>_temp_june`: average June temperature (degrees Celsius), one pixel image (`im`) per year;
- `<prefix>_precip_june`: average June precipitation (mm/m^2), one `im` per year;
- `<prefix>_wind_june`: average June wind speed (km/h), one `im` per year;
- `<prefix>_ppp_june`: wildfire locations, one point pattern (`ppp`) per year.

"orig_scale" means that the covariates are in their original units; the rescaling to [0, 1] is done in the code. Coordinates are in pixel units of the grid; pixels outside the province are NA. The conversion to longitude and latitude for the axis labels is in `CANADA_GEO` in `../../R/canada_data.R`.

The files are read by `../../R/canada_data.R` without loading `spatstat` (the objects are read as plain lists). The pre-processing (selection of the years, empirical-c.d.f. transform of the covariates, pixel quadrature, assignment of border events to the nearest valid pixel) is described in `../README-APPLICATION.md`. `run_canada.R` stores MD5 checksums of these files in each run's `config.rds`, so the files must not be modified. The original files for Alberta and Manitoba and a natural-scale Ontario file, not used in the paper, are in `archive.zip` (`canada_original/`).
