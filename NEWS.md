# truncpois 0.1.0

* Initial release.
* Implements `dtruncpois`, `ptruncpois`, `qtruncpois`, and `rtruncpois` for the left-truncated, right-truncated, and doubly-truncated Poisson distribution.
* Implements `extruncpois`, `vartruncpois`, `medtruncpois`, and `modtruncpois` for distributional moments using closed-form formulae.
* All computations use log-scale arithmetic for numerical stability.
* Three random-sampling methods available in `rtruncpois`: `"direct"`, `"inversion"`, and `"bounded"`.
