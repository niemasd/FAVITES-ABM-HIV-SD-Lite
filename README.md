# FAVITES-ABM-HIV-SD-Lite
A simplified simulation pipeline specifically for ABM-HIV SD

Assumes that [`abm_hiv-HRSA_SD`](https://github.com/mathematica-pub/abm_hiv/tree/HRSA_SD) is in `/usr/local/bin`, so the [`abm_hiv_commandline.R`](https://github.com/mathematica-pub/abm_hiv/blob/HRSA_SD/abm_hiv_commandline.R) script is located at `/usr/local/bin/abm_hiv-HRSA_SD/abm_hiv_commandline.R` and the [`modules`](https://github.com/mathematica-pub/abm_hiv/tree/HRSA_SD/modules) folder is located at `/usr/local/bin/abm_hiv-HRSA_SD/modules` (but the user can override via command-line arguments).

Example files can be found in [`examples`](examples). Example usage:

```bash
rm -rf tmp && ./run_favites_lite.py -o tmp --abm_hiv_params_xlsx example/data.xlsx --abm_hiv_trans_start .25 --abm_hiv_trans_end .5 --abm_hiv_trans_time 25 --sample_time_probs_csv example/time_probs.csv --abm_hiv_sd_demographics_csv example/demographics.csv --coatran_eff_pop_size 3000 --time_tree_seed example/seed_tree.nwk --mutation_rate_loc 0.0000667 --mutation_rate_scale 0.000041667 --time_tree_tmrca 1938.34
```
* Effective Population Size (~3000) from [Seo *et al*. (Genetics 2002)](https://doi.org/10.1093/genetics/160.4.1283)
