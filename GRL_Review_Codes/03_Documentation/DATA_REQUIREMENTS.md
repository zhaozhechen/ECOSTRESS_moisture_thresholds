# Data requirements

This file records the inputs needed for each reviewer analysis before implementation begins.

## Clear-sky net-radiation comparison

### Raw AmeriFlux inputs

- Sites: `US-A32` and `US-CF3`
- Temporal resolution: half-hourly (`SUBSET_HH`) or hourly (`SUBSET_HR`)
- Required fields: `TIMESTAMP_START`, `NETRAD`, and the associated `NETRAD` QC field when available
- Current OSC root: `/fs/ess/PAS2204/SharedData/AmeriFlux_All_Sites/`

### Processed manuscript inputs

- `AMF_Full_range_df_US-A32.csv`
- `AMF_Full_range_df_US-CF3.csv`
- Required fields: `time` and `ESI_daily`
- Current OSC root: `/fs/ess/PAS2204/Results/Validation_ET_ESI_All_AMF/Full_range_df/AMF/`

### Software dependencies

R packages: `dplyr`, `ggplot2`, `RColorBrewer`, `cowplot`, and `scales`.
