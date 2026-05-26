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

- Per-case censoring is treated as independent across delays. The
  day-window uncertainty on `T_onset` is shared between the
  `onset → admission` and `onset → notification` observations of the
  same case; treating them as independent inflates posterior variance
  slightly but is unbiased.
- A latent-time joint model was prototyped (v3 with per-case
  `T_onset`, `T_admit` latents and ordering constraints). It hit
  wedge-shaped boundary geometry NUTS handles poorly at the 6
  same-day-admission cases (308 divergent transitions). The marginal
  formulation used here matches what CensoredDistributions.jl
  supports and what Charniga *et al.* 2024 recommend for retrospective
  complete-outbreak data. The natural-history identity is enforced
  via convolution post-processing.
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

- InitFromPrior works here. Per-case latent-time variants would need
  stricter initialisation.
- AutoForwardDiff is used rather than Enzyme —
  `CensoredDistributions`' integral-based likelihood is not yet
  uniformly Enzyme-friendly. ForwardDiff is fast enough at this size.
- Single seed (20260519). The convergence diagnostics are strong but
  a multi-seed run is a useful robustness check.

## Generalisability

BDBV is the third *Ebolavirus* species, with ~250 cases globally
across three outbreaks (Uganda 2007–08, DRC 2012, DRC/Uganda 2026).
Isiro is the only delay parametrisation we have. The Isiro outbreak
was rural, mining-influenced, with MSF ETC support. The 2026 Ituri
outbreak has urban (Bunia) and conflict-zone components that may
push delays longer.

## Charniga checklist compliance

Self-audit against the [Charniga *et al.* 2024](https://doi.org/10.1371/journal.pcbi.1012520)
reporting checklist (Tables 2 and 3). Items grouped by status; the
checklist *table* itself is issue #2.

**Estimation:**

- Double interval censoring — `double_interval_censored(...; interval = 1.0)`
  on all four atomic delays.
- Right truncation — *not applied; not applicable.* The 2012 Isiro
  outbreak ended over a decade before fitting; the deposit is a
  retrospective, complete-outbreak line list.
- Dynamical bias — *not applied; not applicable* for the same reason.
- Multiple distributions fitted — LogNormal, Gamma, Weibull compared
  via WAIC. Gamma wins; reported as canonical.
- Distributions visualised — per-family posterior predictive panels
  and a family-comparison figure in the analysis walkthrough.
- Parameter conversion — mean, SD, shape, scale reported per draw and
  summarised with 95% CrI for each atomic delay.
- Stratified estimates — CFR is stratified by HCW status, case
  definition, and standardised age. A `bdbv_model_stratified`
  HCW-stratified delay model exists in the package but is not surfaced
  in the walkthrough (the n = 52 sample size keeps the HCW × delay
  estimates very wide; treated as a follow-up).
- Model diagnostics — R̂, bulk ESS, and divergent transition counts
  reported for every family in the comparison table.

**Reporting:**

- Variability — SD reported alongside mean and median for atomic and
  convolved marginals, with 95% CrI throughout.
- Quantiles and distribution parameters — median, mean, SD, shape and
  scale for each atomic delay; median, mean, SD and 95th percentile
  for the convolved marginals.
- Credible intervals — 95% CrI throughout.
- Contextual information — sample size, observation window,
  age/sex/HCW breakdown, case-definition mix and outcome counts in the
  walkthrough's "Outbreak context" block; epidemic curve adjacent; ETC
  and outbreak setting in the [Generalisability](#generalisability)
  section.
- Code and data — full source on
  [GitHub](https://github.com/epiforecasts/bdbv-linelist-analysis);
  data CSV bundled; outputs published as the `main-latest` rolling
  release; Zenodo metadata in `.zenodo.json`.
- Posterior samples — `output/posterior_<family>.csv` exposes
  per-draw mean, median, SD, shape and scale for each atomic delay and
  the convolved marginals; the raw MCMC chain is not serialised by
  default (downstream users typically only need the per-draw
  parameters; the CSV is the canonical deliverable).

**Incubation period / serial interval (Charniga Table 3):**

- Incubation period, serial interval and generation interval — not
  fitted. No exposure dates and no transmission pairs in the deposit
  (see [Data](#data) above). Items on multiple exposures, transmission
  pair confirmation, transmission direction and negative-interval
  distributions are all not applicable for the same reason.

**Other recommendations:**

- Time-varying delays — not implemented. n = 52 is too small for a
  pre-/post-ETC cohort split on the forward delays. The death-pathway
  mixture in `fit_death_mixture` partially addresses this for the
  onset → death marginal.
- Priors — weakly informative Normal priors on log-mean / log-shape;
  three-point sensitivity sweep (`prior_scale ∈ {0.5, 1.0, 2.0}`).
- Pooled vs meta-analysis — individual-level reanalysis of the
  original line list rather than a meta-analysis, so the meta-analysis
  sensitivity item is not applicable.

## Out of scope

- Whether the HCW protective effect is causal (faster access to care)
  or selection (HCWs were preferentially confirmed). The
  probable-vs-confirmed adjustment partially controls for the latter
  but does not eliminate it.
- Transmission. Use Zaire-EVD generation interval priors for
  downstream Rt or nowcasting work.
- Reservoir or spillover dynamics.
