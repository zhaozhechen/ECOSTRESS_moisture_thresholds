# AmeriFlux alpha sensitivity to excluding low net radiation

## Results

| Site | Complete observations, all / filtered | Alpha, all | Alpha, filtered | Change |
|---|---:|---:|---:|---:|
| US-A32 | 210 / 157 | 52.2122 degrees | 55.3196 degrees | +3.1074 degrees |
| US-CF3 | 124 / 92 | 83.5511 degrees | 79.5899 degrees | -3.9613 degrees |

The all-observation boundary slopes reproduce the original saved values of 1.2897607 and 8.8470461 within 1e-6. Complete observations retained are 74.76% and 74.19%, respectively. No complete alpha observations lack daily net radiation.

## Exact calculation

1. Read the existing AMF_Full_range_df files for US-A32 and US-CF3. ESI is reused without recalculation or z-score transformation. These files retain the original preprocessing: the most recent four years at most, May-September, exclusion of precipitation above 1 mm/day and temperature below 0 degrees C, and the original variable quality controls. Thus, 'all' means all observations eligible under the original analysis, not all calendar days.
2. Read each site's raw half-hourly or hourly AmeriFlux SUBSET file. Retain the most recent 4*365 days if longer, as in Clear_Sky_comparison.R. Select NETRAD with the existing Var_QC function (QC=0 where available), replace -9999 by NA, set negative radiation to zero, then average within each day using the original Daily_mean function. This is the original nonnegative-radiation convention, not a signed 24-hour net-radiation mean.
3. Join daily radiation to the existing processed observations by date. Calculate the 25th percentile separately for every year-month, among rows with nonmissing ESI, exactly as in the previous radiation comparison. R's default quantile type 7 is used. The cutoff population is selected before requiring complete SM and VPD. Save every cutoff and its sample count.
4. Fit the original full processed dataframe for the all-observation result. Fit rows with daily radiation greater than or equal to their year-month cutoff for the filtered result. Ties at the cutoff are included. Rows lacking radiation cannot enter the filtered dataset; none of the complete baseline observations here have that issue.
5. Apply the functions from CONUS_Threshold/Validation_Final_2/Validation_thresholds_20250202.R. The separate helper file preserves the original code and comments, with only summarise changed to transmute in Normalize_df for compatibility with current dplyr. This preserves the original multiple-row numerical output.
6. Within EACH input dataset, independently rescale ESI, SM, and VPD to 0-100 using their observed minimum and maximum. Subtract the rescaled median ESI to obtain Z, then remove incomplete ESI/SM/VPD rows. These are min-max coordinates, not empirical percentile ranks, despite the original 'quantile' terminology. Normalization occurs before complete-case removal, as in the original code.
7. Average Z in the 7x7 SM-VPD grid. Apply the original 3x3 Gaussian kernel with sigma=1 using raster::focal, including its original mean and edge handling. Construct weights from smoothed observation counts. Fit weighted linear regressions to interior row and column slices, requiring at least four valid cells per slice.
8. Retain slices with R-squared >=0.8 and exclude statistically significant slopes with the incorrect direction (negative ESI-SM or positive ESI-VPD). Calculate each retained slice's zero crossing, discard crossings outside observed normalized SM-VPD support, and require at least three crossing points. Fit VPD_coordinate ~ SM_coordinate by ordinary least squares to those points.
9. Calculate alpha = atan(slope)*180/pi. As in the original 20250202 script, no additional final contour filter is imposed in SM_VPD. Both comparisons also have contour R-squared above 0.5 and angles above -20 degrees. The summary includes slopes, intercepts, contour R-squared, normalization ranges, sample counts, and original threshold statistics.
10. Preserve the original uncertainty routine: 100 bootstrap resamples of the boundary points, seed=1, with the original 1.5-IQR outlier removal before taking SD. This is not a bootstrap of independent observation datasets or a statistical test of the difference between the nested datasets. The figure shows point estimates only.

## Interpretation and limitations

The point estimates change by approximately 3-4 degrees after excluding the lower-radiation observations. This describes sensitivity at these two sites; it does not establish statistical equivalence or that clear-sky sampling has no effect across CONUS.

Independent normalization is retained to satisfy the request to rerun exactly the original method on two datasets. Therefore, the comparison includes both changes in retained observations and changes in their reference ranges and median ESI. At US-CF3, for example, the SM normalization range changes from 0.03404-0.32780 to 0.11775-0.30221. The two fitted lines should not be overlaid on a shared normalized coordinate system as if their axes were identical.

The original uncertainty calculation produces a large alpha SD for the US-CF3 baseline (76.40 degrees, based on only four boundary points). This result is retained in the table; the similar point estimates should not be described as a formal demonstration of robustness or a nonsignificant difference.

## Files and reproduction

- Run 01_Codes/AMF_alpha_netradiation.R; it sources 01_Codes/AMF_alpha_functions.R and General_functions.R.
- Figures: 02_Results/Figures/AMF_alpha_netradiation.png and .pdf.
- Summary: 02_Results/Tables/AMF_alpha_netradiation_summary.csv.
- Monthly cutoffs: 02_Results/Tables/AMF_alpha_netradiation_monthly_cutoffs.csv.
- Local audit files: 02_Results/Intermediate/AMF_alpha_netradiation, including daily selection flags, every slice regression, boundary points, input hashes, and R session information. This directory and the source data remain ignored by Git.
- The original server-code archive and existing review analyses were not modified. No commit or push was performed.
