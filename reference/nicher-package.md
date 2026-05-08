# nicher: Estimates Ecological Niche Models Using Ellipses

Ecological niche model estimation using ellipsoidal geometry under an M
hypothesis. Fits seven likelihood families (presence-only, weighted,
inverse-probability-weighted, skew-normal, skew-normal-weighted, skew-t,
and skew-t-weighted) via multi-start optimisation with Sobol sequences.
Methods for the optimisation of ellipses parameters are as described in
Jimenez et al. (2022)
[doi:10.1016/j.ecolmodel.2021.109823](https://doi.org/10.1016/j.ecolmodel.2021.109823)
.

Estimates ecological niche models by fitting multivariate distributions
(Gaussian, skew-normal, or skew-t) to species occurrence data in
environmental space.

## Likelihood families

Seven likelihood formulations are available via \[optimize_niche()\]:

- `"presence_only"`:

  Multivariate normal density estimated from presence records alone
  (Jimenez et al. 2019, Eq. 2).

- `"weighted"`:

  Paper-faithful weighted-normal model with KDE weights accounting for
  the density of available environments in M (Jimenez & Soberon 2022,
  Eq. 5/8), plus a ridge prior on \\\log\sigma\\. Default and
  recommended.

- `"ip_weighted"`:

  Inverse-probability-weighted (IPW) normal. Assumes the Jimenez et
  al. (2019) DGP (presences drawn from \\f\\ restricted to \\M\\) and
  uses \\w = 1/\hat{g}\\ (Horvitz–Thompson correction) to estimate the
  fundamental niche on uniform environmental space (Lebesgue measure),
  removing the distortion caused by non-uniform environmental density in
  \\M\\. See
  [`loglik_niche_math_ip_weighted`](https://alrobles.github.io/nicher/reference/loglik_niche_math_ip_weighted.md)
  for the full derivation, model comparison, and known limitations.
  Formerly `"kde_bias_corrected"`.

- `"skew_normal"`:

  Presence-only multivariate skew-normal (Azzalini & Capitanio 1999).

- `"skew_normal_weighted"`:

  Weighted skew-normal (Eq. 8 with skew-normal density).

- `"skew_t"`:

  Presence-only non-central skew-t (Branco & Dey 2001) via 32-node
  Gauss-Laguerre quadrature.

- `"skew_t_weighted"`:

  Weighted skew-t (Eq. 8 with NCST density).

See the package website for the full mathematical specification of each
model.

## References

Jimenez, L., Soberon, J., Christen, J. A., & Soto, D. (2019). On the
problem of modeling a fundamental niche from occurrence data.
*Ecological Modelling*, 397, 109823.
[doi:10.1016/j.ecolmodel.2019.01.009](https://doi.org/10.1016/j.ecolmodel.2019.01.009)

Jimenez, L., & Soberon, J. (2022). Weighted-normal model for the
fundamental niche. *Ecological Modelling*, 438, 109982.

Azzalini, A., & Capitanio, A. (1999). Statistical applications of the
multivariate skew normal distribution. *J. R. Stat. Soc. B*, 61(3),
579–602.

Branco, M. D., & Dey, D. K. (2001). A general class of multivariate
skew-elliptical distributions. *J. Multivariate Anal.*, 79(1), 99–113.

## See also

Useful links:

- <https://alrobles.github.io/nicher/>

- <https://github.com/alrobles/nicher-devel>

- Report bugs at <https://github.com/alrobles/nicher-devel/issues>

## Author

**Maintainer**: Angel Robles <a.l.robles.fernandez@gmail.com>

Authors:

- Laura Jimenez <ljimenez@cmm.uchile.cl>
  ([ORCID](https://orcid.org/0000-0002-4674-4270)) \[contributor\]
