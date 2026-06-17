# covsim 1.2.0

New calibrate/simulate API: `calibrate_ig()`, `calibrate_plsim()`, and
`calibrate_vita()` return reusable calibration objects. Sampling is done via
S3 methods on `stats::simulate()`, so a single (potentially expensive)
calibration can drive many simulations without recomputation. The existing
`rIG()`, `rPLSIM()`, and `vita()` functions are unchanged in signature and
return shape; they now delegate to the new API.

`vita()` is roughly 3–4× faster. The pair-copula parameter search now uses
`rvinecopulib::inverse_rosenblatt()` with cached uniform inputs so the
root-finding objective is deterministic, replacing the previous multi-stage
stochastic search with a single `uniroot()` call. The `numrootpoints`,
`conflevel`, and `numpoints` arguments are kept for backward compatibility
but no longer have any effect. Seeded outputs are not bit-identical with
covsim 1.1.0; the matched target covariance is reached at similar precision.

# covsim 1.0.0

Third release. Accompanies journal publication in JSS. 

# covsim 1.1.0

rIG: The A matrix is now triangular by default

rPLSIM: Added tests for whether an identical pair has been calibrated. 
This speeds up rPLSIM if sigma.target has many identical values, and skewness and kurtosis values are identical. 

vita: Now checks margins are correctly named. 
