## Test environments
* local R installation, version 4.3.1 (2023-06-16)
* x86_64-apple-darwin15.6.0 (64-bit)
* macOS Big Sur 11.6.8

## R CMD check results 

0 errors | 0 warnings | 0 notes

## 0.1.0 Submission

After manual feedback from CRAN, I reduced length of title, and reduced
information messaging to console. I now only have a reduced number of if(verbose) cat(), as suggested by Jelena Saf at CRAN.



### 0.2.0

Added new function rPLSIM 


## 0.2.1

* Added  vignette

## 0.2.2

* Updated citation

## 1.0.0

* Updated vignette
* Updated with reference to new JSS article

The DOI in CITATION is for a new JSS publication that will be registered after publication on CRAN. 

## 1.1.0

* Improved speed for rPLSIM
* Changed default A type in rIG from symmetric to triangular
* Added check of correct margin names in vita. 

## 1.2.0

* New calibrate/simulate API: `calibrate_ig()`, `calibrate_plsim()`,
  `calibrate_vita()` return reusable calibration objects; sampling is done
  via S3 methods on `stats::simulate()`. Existing `rIG()`, `rPLSIM()`,
  `vita()` are unchanged in signature and return shape and now delegate
  to the new API.
* `vita()` is roughly 3-4x faster. The pair-copula parameter search now
  uses `rvinecopulib::inverse_rosenblatt()` with cached uniform inputs
  so the root-finding objective is deterministic, replacing the previous
  multi-stage stochastic search with a single `uniroot()` call.
* Pair-copulas within each tree are now calibrated in parallel via
  `parallel::mclapply` when `cores > 1`.
* New vignette: `calibrate-simulate.Rmd`.

(Note: update the "Test environments" section above with your actual
local R version and the results of `devtools::check_win_devel()` and
`rhub::rhub_check()` before submitting.)



