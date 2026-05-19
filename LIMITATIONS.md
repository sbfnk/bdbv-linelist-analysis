# Limitations

## Data

- **n = 52 cases** is small. Posteriors are tight where the delay
  component is well-anchored (onset → admit n = 40, onset → death
  via convolution effectively n ≥ 22) and wider where data are sparse
  (admit → discharge n = 15) or heavy-tailed (onset → notification
  with max observed 86 d).
- **No transmission pairs** in the public deposit. Serial interval,
  generation interval and R₀ cannot be estimated from this data set.
  The Rosello deposit is canonical for Isiro — the underlying source
  data is not available, so this is not fixable upstream.
- **No exposure dates.** Incubation period cannot be fitted from
  this data. Borrow MacNeil 2010 (mean 6.3 d, 95% CI 5.2–7.3, n = 24,
  2007 Uganda) for downstream uses.
- **5 admission-date encoding outliers** (−89, −5, −4, −1, 328720
  days from onset) are programmatically set to missing during
  loading. **1 notification-date outlier** (−62 days) is similarly
  dropped. Both are encoding errors in the deposit (likely
  onset/secondary date swaps) — we cannot verify the direction with
  the custodian, so they are conservatively dropped rather than
  swapped.
- **Three cases have missing age** (3 / 52). Imputed at the sample
  mean before standardisation. With only 3 imputations the effect on
  the age coefficient is small but non-zero.
- **Notification has a long right tail** (47, 50, 53, 60, 65, 75, 86
  days after onset) — likely retrospective notifications of
  community deaths investigated post hoc. The Gamma fit handles this
  better than LogNormal (the basis of the WAIC win); even so the
  mean (≈ 20 d) is markedly above the median (≈ 11 d).

## Model

- **Per-case censoring is treated as independent across delays.**
  In reality the day-window uncertainty on `T_onset` is shared
  between the `onset → admission` and `onset → notification`
  observations of the same case. Treating these as independent in
  the likelihood inflates posterior variance slightly but is
  unbiased.
- **A full latent-time joint model was prototyped** (v3 with per-case
  `T_onset`, `T_admit` latents and ordering constraints). It hit
  wedge-shaped boundary geometry NUTS handles poorly at the 6
  same-day-admission cases (308 divergent transitions). The marginal
  formulation used here matches what CensoredDistributions.jl
  naturally supports and what Charniga *et al.* 2024 recommends for
  retrospective complete-outbreak data. The natural-history identity
  is still enforced via convolution post-processing.
- **No latent severity** linking delays and mortality. The delay
  fits and the CFR fits cannot directly share information about
  case severity. A v3 with per-case latent severity is a defensible
  follow-up but would not change the marginal headline numbers
  appreciably at this sample size.
- **Age effect is linear in standardised age.** Strong
  nonlinearities (e.g. U-shape with high mortality at extremes)
  cannot be captured.
- **Sex is not in the CFR model.** Rosello records sex (40 F / 12 M
  in Isiro 2012) but Kratz 2015's community-vs-ETC stratification
  is not directly reproducible from the deposit. Adding sex would
  be one extra coefficient and is a defensible follow-up.
- **5 community deaths** (died without recorded admission) contribute
  to the CFR block but not to the delay model — the model has no
  pathway for non-admitted deaths.
- **Stratified delay model is fitted for Gamma and LogNormal only**
  (Weibull's shape prior needed extra truncation that complicated
  the stratified parametrisation).

## Inference

- **InitFromPrior** initialisation works for this model. Larger
  variants with per-case latent event times would need stricter
  initialisation.
- **AutoForwardDiff** is used rather than Enzyme —
  `CensoredDistributions`' integral-based likelihood is not yet
  uniformly Enzyme-friendly. ForwardDiff is fast enough at this
  size (~1–2 min per family on a laptop).
- **Single seed** (20260519) reported. Sensitivity to seed should
  be small given the strong convergence diagnostics, but a
  multi-seed run is a sensible robustness check.

## Generalisability

- BDBV is the third species in the *Ebolavirus* genus; sample size
  globally is ~250 cases across 3 outbreaks (Uganda 2007–08, DRC
  2012, DRC/Uganda 2026 ongoing). The Isiro-derived delay
  distributions are the only delay parametrisation we have for BDBV.
  Whether they generalise to other settings depends on health-
  system context: the Isiro outbreak was rural, mining-influenced,
  with MSF ETC support; the current 2026 Ituri outbreak has urban
  (Bunia) and conflict-zone components that may push delays longer.

## Things this analysis cannot tell you

- Whether the HCW protective effect is causal (faster access to
  care) or selection (HCWs were preferentially confirmed). The
  probable-vs-confirmed adjustment partially controls for the
  latter but does not eliminate it.
- Anything about transmission. Use Zaire-EVD generation interval
  priors when running `EpiNow2` or similar downstream.
- Reservoir or spillover dynamics.
