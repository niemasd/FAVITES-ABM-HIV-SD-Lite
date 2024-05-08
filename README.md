# FAVITES-ABM-HIV-SD-Lite
A simplified simulation pipeline specifically for ABM-HIV SD

Assumes that [`abm_hiv-HRSA_SD`](https://github.com/mathematica-pub/abm_hiv/tree/HRSA_SD) is in `/usr/local/bin`, so the [`abm_hiv_commandline.R`](https://github.com/mathematica-pub/abm_hiv/blob/HRSA_SD/abm_hiv_commandline.R) script is located at `/usr/local/bin/abm_hiv-HRSA_SD/abm_hiv_commandline.R` and the [`modules`](https://github.com/mathematica-pub/abm_hiv/tree/HRSA_SD/modules) folder is located at `/usr/local/bin/abm_hiv-HRSA_SD/modules` (but the user can override via command-line arguments).

Example files can be found in [`examples`](examples). Example usage:

```bash
rm -rf tmp && ./run_favites_lite.py -o tmp --abm_hiv_params_xlsx example/data.xlsx --abm_hiv_trans_start .25 --abm_hiv_trans_end .5 --abm_hiv_trans_time 25 --sample_time_probs_csv example/sample_time_probs.csv --abm_hiv_sd_demographics_csv example/demographics.csv --coatran_eff_pop_size 3000 --time_tree_seed example/time_tree_seed.nex --mutation_rate_mean -8.670092204790605 --mutation_rate_sigma 1.855881729127663 --time_tree_tmrca 1938.34 --sim_start_time 2015
```
* Effective Population Size (~3000) from [Seo *et al*. (Genetics 2002)](https://doi.org/10.1093/genetics/160.4.1283)

Running on real dataset:

```bash
rm -rf tmp && ./run_favites_lite.py -o tmp --abm_hiv_params_xlsx real_data/data.xlsx --abm_hiv_trans_start .25 --abm_hiv_trans_end .5 --abm_hiv_trans_time 25 --sample_time_probs_csv real_data/sample_time_probs.csv --abm_hiv_sd_demographics_csv real_data/demographics.csv --coatran_eff_pop_size 3000 --time_tree_seed real_data/time_tree_seed.nex --mutation_rate_mean -8.670092204790605 --mutation_rate_sigma 1.855881729127663 --time_tree_tmrca 1938.34 --sim_start_time 2019
```

Running calibration:

```bash
rm -rf tmp && ./run_calibration.py -o tmp --abm_hiv_params_xlsx real_data/data.xlsx --abm_hiv_trans_start .25 --abm_hiv_trans_end .5 --abm_hiv_trans_time 25 --sample_time_probs_csv real_data/sample_time_probs.csv --abm_hiv_sd_demographics_csv real_data/demographics.csv --coatran_eff_pop_size 3000 --time_tree_seed real_data/time_tree_seed.nex --mutation_rate_mean -8.670092204790605 --mutation_rate_sigma 1.855881729127663 --time_tree_tmrca 1938.34 --sim_start_time 2019 --calibration_csv real_data/sd_calibration.csv --calibration_mode epi+genetic --zip_output
```
