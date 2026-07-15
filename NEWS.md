# truncpois 0.1.0

* Initial release.
* Implements `dtruncpois`, `ptruncpois`, `qtruncpois`, and `rtruncpois` for the left-truncated, right-truncated, and doubly-truncated Poisson distribution.
* Implements `extruncpois`, `vartruncpois`, `medtruncpois`, and `modtruncpois` for distributional moments using closed-form formulae.
* All computations use log-scale arithmetic for numerical stability.
* Three random-sampling methods available in `rtruncpois`: `"direct"`, `"inversion"`, and `"bounded"`.
* Implements `mletruncpois` for maximum likelihood estimation of the rate parameter from observed counts, including a standard error for the estimate.
* Implements `plottruncpois` for visualizing the PMF, CDF, or quantile function of a truncated Poisson distribution, with an optional overlay of the corresponding non-truncated distribution.
