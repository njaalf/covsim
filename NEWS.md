# covsim 1.2.0

New calibrate/simulate API: `calibrate_ig()`, `calibrate_plsim()`, and
`calibrate_vita()` return reusable calibration objects. Sampling is done via
S3 methods on `stats::simulate()`, so a single (potentially expensive)
calibration can drive many simulations without recomputation. The existing
`rIG()`, `rPLSIM()`, and `vita()` functions are unchanged in signature and
return shape; they now delegate to the new API.

# covsim 1.0.0

Third release. Accompanies journal publication in JSS. 

# covsim 1.1.0

rIG: The A matrix is now triangular by default

rPLSIM: Added tests for whether an identical pair has been calibrated. 
This speeds up rPLSIM if sigma.target has many identical values, and skewness and kurtosis values are identical. 

vita: Now checks margins are correctly named. 
