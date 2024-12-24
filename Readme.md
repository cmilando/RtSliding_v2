# Sliding window

Bayesian sliding window for R(t) work in 1 dimension

use in order:
* `01_simulate_data.R`
* `02_run_stan.R`
* `03_extract_output.R`
* `04_plot.R`

Doing it this way means you can get a posterior distribution for M(t)!

### Back-imputation
Right now the way this is fixed is just imputing what the M[0] value was. I tried implementing the exponential growth model from EpiEstim for back-calculation but it didn't quite work. I think its because its being applied in the serial interval part of the calculation. I think w would have to come into play but imputing the original value seems to work well enough. 

### next steps
* do `Flu2009`
* adding the J dimension back in -- confirm that this works if P is diagonal ones. This should be "pretty easy" -- just any time you see a n, add a j. 

![Alt text](plot.png)


