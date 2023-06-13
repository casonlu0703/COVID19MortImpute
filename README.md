This README describes the `Rscripts` used to impute county-level age-specific provisional COVID-19 deaths in 2020. The order of the following scripts is the order of conducting the analysis. 

1. `create_data.R`: creates and saves data sets in the `data/` directory. 
2. `data_explore.R`: conducts data exploration to investigate the prevalence and pattern of suppressed data in the county-level provisional COVID-19 deaths. 
3. `impute/impute_model.R`: interfaces between `RStan` and `Stan` to conduct Bayesian imputation. There are a total of 26 models created and examined over time. The final model selected are model \#8 (M1), model \#16 (M2), and model \#18 (M3). The corresponding `stan` scripts specified in this `Rscript` is in the directory `impute/stan/`
4. `impute/gather_impute_result.R`: obtains diagnostics for each model to assess model performance; creates sampled datasets from the Bayesian model to reduce computation burden for post-analyses. 
5. `impute/aggregate_death.R`: creates post simulation analyses, reports, and figures. 

