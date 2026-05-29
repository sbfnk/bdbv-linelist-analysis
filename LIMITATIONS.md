# Limitations

## Data

- n = 52 cases is small. Posteriors are tight for the well-anchored
  delays (onset → admit n = 40; onset → death via convolution
  n ≥ 22) and wide for sparse or heavy-tailed ones (admit → discharge
  n = 15; onset → notification has max observed 86 d).
- No transmission pairs in the deposit. Serial interval, generation
  interval and R₀ cannot be estimated. The Rosello deposit is
  canonical for Isiro and the underlying source data is not available.
- No exposure dates. Incubation period cannot be fitted. Use
  MacNeil 2010 (mean 6.3 d, 95% CI 5.2–7.3, n = 24, 2007 Uganda)
  instead.
- 5 admission-date encoding outliers (−89, −5, −4, −1, 328720 days
  from onset) and 1 notification-date outlier (−62 days) are set to
  missing during loading. These look like onset/secondary date swaps
  but we cannot verify the direction with the custodian, so they are
  dropped rather than swapped.
- 3 of 52 cases have missing age, imputed at the sample mean before
  standardisation.
- Notification has a long right tail (47, 50, 53, 60, 65, 75, 86 days
  after onset) — likely retrospective notifications of community
  deaths investigated post hoc. Gamma handles this better than
  LogNormal; the mean (≈ 20 d) is still well above the median
  (≈ 11 d).

## Model

- The HCW-stratified delay model (`bdbv_model_stratified`) still uses
  the marginalised `double_interval_censored` formulation and so
  treats per-case censoring as independent across delays. The main
  `bdbv_model` was converted to per-case latents (see "Model" in
  `MODEL.md`); the stratified variant is a smaller analysis used only
  for the HCW-vs-non-HCW odds-ratio table and has not yet been
  ported.
- No latent severity links delays and mortality. The delay fits and
  the CFR fits do not share information about case severity.
- Age effect is linear in standardised age. Nonlinearities
  (e.g. U-shape) are not captured.
- Sex is not in the CFR model. Rosello records sex (40 F / 12 M in
  Isiro 2012) but Kratz 2015's community-vs-ETC stratification is
  not reproducible from the deposit.
- 5 community deaths (died without recorded admission) contribute to
  the CFR block but not to the delay block; the model has no pathway
  for non-admitted deaths.
- The stratified delay model is restricted to the Gamma family. The
  LogNormal parametrisation in the stratified model would need a
  different SD prior to be consistent with the unstratified fit, and
  the Weibull shape prior needed extra truncation that complicated
  the stratified parametrisation. Gamma is the WAIC winner anyway.

## Inference

- InitFromPrior works here. The bounded-primary reparametrisation
  on the per-case latents means NUTS sees a smoothly-bounded support
  (no hard `t_s > t_p` constraint), so warm-up is well-behaved.
- Single seed (20260519). The convergence diagnostics are strong but
  a multi-seed run is a useful robustness check.

## Generalisability

BDBV is the third *Ebolavirus* species, with ~250 cases globally
across three outbreaks (Uganda 2007–08, DRC 2012, DRC/Uganda 2026).
Isiro is the only delay parametrisation we have. The Isiro outbreak
was rural, mining-influenced, with MSF ETC support. The 2026 Ituri
outbreak has urban (Bunia) and conflict-zone components that may
push delays longer.

## Out of scope

- Whether the HCW protective effect is causal (faster access to care)
  or selection (HCWs were preferentially confirmed). The
  probable-vs-confirmed adjustment partially controls for the latter
  but does not eliminate it.
- Transmission. Use Zaire-EVD generation interval priors for
  downstream Rt or nowcasting work.
- Reservoir or spillover dynamics.

## Downstream priors

For downstream re-applications anchored on the 2012 Isiro outbreak
(e.g. `epiforecasts/BVDOutbreakSize`), the exploratory exponential
fit in [Analysis walkthrough §
Early-phase growth rate](analysis.md#Early-phase-growth-rate) is the
recommended source of a growth-rate prior. As a one-line summary,
use **`Normal(0.008, 0.005)` on `r` per day** (posterior mean and SD
from a Poisson regression on weekly onset counts from week 1 through
the peak week). The 95% CrI on `r` covers zero — Isiro was a slow,
noisy rise — so treat this as a weakly-informative prior rather than
a tight constraint, and re-run `fit_growth_rate(ll)` if you want the
live posterior summary.
