# Future range analysis inputs

`Future_outside_historical_range.R` compares adjusted daily CMIP6 VPD and SM against the pixel-specific historical observational ranges used by the final projection workflow.

## Required local structure

```text
Future_range/
├── CMIP6_adjusted/
│   └── <model>_<Mid-or-End>_<ssp245-or-ssp585>.rds
└── Reference/
    ├── Combined_0.25D_Intercepts.rds
    └── Hist_SM_range_ls.rds
```

## Server sources

Adjusted future projections:

`/fs/ess/PAS2204/Results/CONUS_Threshold_Final/CMIP6_Projection_Final3/CMIP6_adjusted/`

Copy the 60 model files corresponding to 15 models, two periods (`Mid` and `End`), and two scenarios (`ssp245` and `ssp585`). The required pattern is `<model>_<Mid-or-End>_<ssp245-or-ssp585>.rds`. Files beginning with `Mean_models_` are not required.

Reference raster:

`/fs/ess/PAS2204/Results/CONUS_Threshold_Final/CMIP6_Projection_Final2/0.25D_intercepts/Combined_0.25D_Intercepts.rds`

Historical SM empirical distributions:

`/fs/ess/PAS2204/Results/CONUS_Threshold_Final/CMIP6_Projection_Final3/Hist_obs/Hist_SM_range_ls.rds`

The historical monthly observation files, original NetCDF files, `CMIP6_raw` files, frequency stacks, and existing projection figures are not required.

Because the 60 adjusted daily projection files may be large, an alternative is to run `Future_outside_historical_range.R` on the server after replacing its four local paths, then download only the compact output folder.
