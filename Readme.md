# Sliding window

Bayesian sliding window for R(t) work in 1 dimension

use in order:
* `01_simulate_data.R`
* `02_run_stan.R`
* `03_extract_output.R`
* `04_plot.R`

Essentially this is one step beyond the Cori 2013 method. But doing it this way means you can get a posterior distribution for M(t)!

### TO-DO
* re-implement back imputation and negative indexing

### next steps
* adding the J dimension back in -- confirm that this works if P is diagonal ones
* do `Flu2009`


![Alt text](plot.png)


