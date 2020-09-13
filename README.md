# Meta-analysis_COPD

## AIM: Meta-analysis of gene expression experiments on COPD patients

This repository was created to perform a meta-analysis from COPD patients using public data.

All data was downloaded into the computer `DNA.lavis.unam.mx` and `10.200.0.42`(restricted access) and raw data downloaded and normalized (e.i. RMA).

We followed this steps, and script can be found in `/vignettes`:

- 0: Data selection `/0-Data-selection.RMD`

- 1: Dowload RAW data `/1-download_raw-data.RMD`

- 2: Normalizing data `/2-normalizing-data.RMD`

- 3: Curating meta information `/3-curatig-data.RMD`

- 4: Differential gene expression analysis `/4-DE.RMD`

- 5: Joing data `/5-merge-data.RMD`

- 6: Meta-analysis `/6-meta-analysis.RMD`

- 7: Enrichment-analysis `/7-enrichment-analysis.RMD`

- 8: Dash-board `/8-dashboard.RMD`


Additional files that can be found in:
https://drive.google.com/drive/folders/1Lw57iquTxWMonUanZPo5G-Yb0QZDDVdP?usp=sharing

For running the dashboard app you can:
```{r}
rmarkdown::run("vignettes/2020-06-13_Meta-analysis_7GSE/Meta-analysis_Dashboard.Rmd")
```



