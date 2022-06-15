This README describes the directory `impute/` and the Stan models. Here is the outline of 
the directory structure. 
   
```
/impute/
  |
  |-- impute_model.R
  |-- impute_model_retired.R
  |-- gather_impute_result.R
  |-- aggregate_death.R
  |-- bayes_impute_aggXX.RDS
  |-- impute1.sh
  |
  |-- stan 
  |     |-- impute_hurdle_agg_nat.stan
  |     |-- impute_hurdle_agg_nat_prior.stan
  |     |-- impute_hurdle_agg_nat_prior_bifur.stan
  |     |-- impute_hurdle_structure_prior.stan
  |     |-- retired stan 
  |
  |-- retired scripts
```

Based on this directory structure, we will describe the top level files (`R` scripts, `.RDS` dataset, 
and a shell file) and the Stan models in the `stan/` folder. 
\
\
At the top level, there are three `R` scripts that perform the entire Bayesian sampling process and data summary:`impute_model.R`/`impute_model_retired.R`, `gather_impute_result.R`, `aggregate_death.R`, following the order of execution. The `.RDS` datasets are the 1,000 posterior sample sets for the suppressed data by each 
model specification. The shell file `impute1.sh` is an example shell file used on the Aspen cluster computers. 
\
\
In the `stan/` folder, it contains the Stan model used in different model specifications. 
\
\
In this document, we will focus only on the description of the three process `R` scripts and their relation 
with each other and with other files. 


**1. `impute_model.R`/`impute_model_retired.R`**

  Both scripts were used for estimating various Bayesian models.
  Originally, `impute_model_retired.R` was used, including
  all model specifications from model 1 through 22. In the beginning, the evolution of the model specifications
  focused more on the combination of variables, interaction terms, and adding more targets in the likelihood. 
  Later, the model specification focused more on the manipulation of prior distributions of the suppressed data. 
  To simplify and focus more on the current few model changes, `impute_model.R` was created to fit Bayesian 
  models with different prior assumptions as variations of model 8. In `impute_model.R`, 
  it includes models 15-24. 
  \
  \
  In a script, each model specification will call the corresponding stan model in the `stan/` folder and 
  conduct HMC to get the model fits. At the end of the script, a stan object will be generated for further
  model diagnosis and calculation. Due to the size of the fitted stan models, it is impossible for us to 
  share these on GitHub. 

**2. `gather_impute_result.R`**
  
  After a fitted Bayesian model was generated, we used this script to perform model diagnosis, 
  summarize data, and perform posterior predictive check. The model diagnosis includes the R-hat, effective
  sample size, and leave-one-out cross validation. In addition, we examined the trend of the distribution
  that combines imputed suppressed and unsuppressed COVID-19 deaths by county, age, and quarter. Finally, 
  we conducted the posterior predictive check to compared the county-level COVID-19 deaths between the 
  simulations from the fitted Bayesian model and observed county-level data. 

**3. `aggregate_death.R`**

  This script sampled 1,000 posterior sets from the fitted Bayesian model and used the 1,000 sample sets
  to validate the data reported form the state level and the national level. The 1,000 sample sets are all saved 
  as `bayes_impute_aggXX.RDS` at the top level directory. 
  



  
  

