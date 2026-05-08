# truncpois 0.1.0

* Initial release.
* Implements `dtruncpois`, `ptruncpois`, `qtruncpois`, `rtruncpois` for the truncated Poisson distribution.
* Implements `extruncpois`, `vartruncpois`, `mtruncpois`, `mode_truncpois` for distributional moments.
* All computations use log-scale arithmetic for numerical stability.
* Three random-sampling methods available in `rtruncpois`: `"direct"`, `"inversion"`, `"bounded"`.
