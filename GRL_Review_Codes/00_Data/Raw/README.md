# Raw data

Raw datasets needed for reviewer analyses are documented here but are not committed because many files are large.

For `Clear_Sky_comparison.R`, the raw inputs are AmeriFlux BASE data for `US-A32` and `US-CF3`. Each site requires its half-hourly or hourly SUBSET CSV containing:

- `TIMESTAMP_START`
- a net-radiation variable matching `NETRAD`
- the corresponding `NETRAD` quality-control variable, when available

The script currently expects these site folders under:

`/fs/ess/PAS2204/SharedData/AmeriFlux_All_Sites/`
