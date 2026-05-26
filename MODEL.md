# Model

Two blocks fitted to the 2012 Isiro BDBV line list (n = 52):

1. Four atomic delay components, each a doubly-censored distribution
   wrapped in `CensoredDistributions.double_interval_censored` with
   `interval = 1.0` (1-day reporting bins on both endpoints, default
   `primary_event = Uniform(0, 1)`).
2. Stratified case-fatality: Bernoulli outcome with a logistic link
   on HCW status, case definition (Probable vs Confirmed), and
   standardised age.

Three parametric families are supported (`:lognormal`, `:gamma`,
`:weibull`), selected via the `family` keyword to `bdbv_model`. The
canonical run uses Gamma — lowest WAIC.

`bdbv_model_stratified(d; family)` adds a `β_*_hcw` log-mean shift to
each delay so HCW and non-HCW cases share the shape parameter but can
differ on the central tendency.

## Common parametrisation across families

All three families share a log-mean / log-shape parametrisation so the
priors are comparable:

| Family   | Latent                | Distribution constructed              |
|---       |---                    |---                                    |
| LogNormal| (log_median, log_sd)  | `LogNormal(log_median, log_sd)`       |
| Gamma    | (log_mean, log_shape) | `Gamma(k, mean/k)` with `k = exp(log_shape)` |
| Weibull  | (log_mean, log_shape) | `Weibull(α, mean/Γ(1+1/α))` with `α = exp(log_shape)` |

For Weibull, `log_shape` is truncated to `(−1, 1)` so `Γ(1 + 1/α)`
doesn't blow up. The Weibull stratified model is not supported.

## Compound delays via convolution post-processing

Fitting `onset → death` and `onset → discharge` separately alongside
their atomic components would violate the per-case identity
`onset → death = onset → admit + admit → death`. So the model fits
only the four atomic components:

| Symbol | Component                | n  |
|---     |---                       |--- |
| `d_oa` | onset → admission        | 40 |
| `d_ad` | admission → death        | 22 |
| `d_ac` | admission → discharge    | 15 |
| `d_on` | onset → notification     | 38 |

and derives `d_od` and `d_oc` in post-processing as sample-level
convolutions. The mean of the convolved marginal equals the sum of
the atomic means by linearity of expectation; medians and quantiles
come from Monte Carlo (500 samples per posterior draw).

Per-case censoring noise across the four delays of the same case is
treated as independent — see [LIMITATIONS.md](LIMITATIONS.md).

## Priors

Weakly informative, centred on a plausible day count for each delay
(`prior_scale` keyword scales the SDs jointly; default 1.0):

| Parameter | Prior |
|---|---|
| `log_median` / `log_mean` for d_oa | Normal(log 3, `prior_scale`) |
| `log_median` / `log_mean` for d_ad | Normal(log 6, `prior_scale`) |
| `log_median` / `log_mean` for d_ac | Normal(log 13, `prior_scale`) |
| `log_median` / `log_mean` for d_on | Normal(log 7, `prior_scale`) |
| `log_sd` (LogNormal) | half-Normal(0, `prior_scale`) |
| `log_shape` (Gamma)  | Normal(0, `prior_scale`) |
| `log_shape` (Weibull)| Normal(0, `prior_scale`), truncated to (−1, 1) |
| `β_0` (CFR baseline)        | Normal(0, 2) |
| `β_hcw`, `β_def`, `β_age`   | Normal(0, 1) |
| `β_*_hcw` (stratified delay shifts) | Normal(0, 0.5) |

The 1.0 default SD on the log-medians gives ≈ ×3-fold prior latitude
either way on the central tendency. A 3-scale sensitivity sweep
(0.5, 1.0, 2.0) shifts posterior means by < 5%.

## Likelihood

For each atomic delay `j ∈ {oa, ad, ac, on}`:

```julia
dic_j = double_interval_censored(dist_j; interval = 1.0)
observed_j ~ Turing.filldist(dic_j, length(observed_j))
```

CFR block:

```julia
η_i = β_0 + β_hcw · 1[HCW_i] + β_def · 1[Probable_i] + β_age · age_z_i
outcome_i ~ Bernoulli(logistic(η_i))
```

`age_z` is standardised age (mean 0, SD 1) with the population mean
imputed for the 3 of 52 missing ages. `logistic(η)` is clamped to
`(1e-10, 1 - 1e-10)` so `Bernoulli`'s domain check doesn't fail when
NUTS proposes large `|η|` during warmup.

## Inference

NUTS via `Turing.jl`. Defaults:

- 1000 post-warmup samples per chain
- 4 chains in parallel (`MCMCThreads`)
- `target_accept = 0.95`
- `AutoForwardDiff()` AD backend
- `DynamicPPL.InitFromPrior()` with a different RNG per chain
- Seed 20260519 (parent), per-chain seeds derived from a child
  `MersenneTwister` so `Random.seed!` is never called on the global
  RNG

About 1–2 minutes per family on a laptop. max R̂ ≤ 1.002, min bulk
ESS ≥ 3,200, 0 divergent transitions across families.

The diagnostics are: R̂ (the ratio of between-chain to within-chain
variance — values near 1 indicate the four chains agree on the
posterior), bulk ESS (effective sample size in the body of the
distribution — values in the thousands indicate the chain has
explored the posterior efficiently), and divergent transitions
(HMC steps that hit numerical instability — zero divergences across
all fits means the geometry is well-behaved).

## Model comparison

WAIC is computed pointwise: for each observation, the per-draw
log-likelihood is evaluated from the per-draw distribution parameters,
and WAIC = −2(lppd − p_waic) aggregated across the 167 observations
(40 + 22 + 15 + 38 atomic delays + 52 outcomes).

PSIS-LOO would be preferable but adds dependencies; WAIC is adequate
here.

## Post-processing

For each atomic delay the script reports the posterior median, mean,
and (for LogNormal) log-SD with 95% CrI. Gamma and Weibull medians
come from per-draw `quantile(dist, 0.5)`.

For each derived marginal (`d_od = d_oa ⊛ d_ad`,
`d_oc = d_oa ⊛ d_ac`) the script samples 500 realisations per
posterior draw from the convolution and reports median, mean, SD,
and 95th percentile with 95% CrI across draws.

For the CFR block the script reports the log-OR coefficients and the
implied marginal CFR at each combination of HCW status and case
definition, with age held at its sample mean.
