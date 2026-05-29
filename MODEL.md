# Model

Two blocks fitted to the 2012 Isiro BDBV line list (n = 52):

1. Four atomic delay components, fitted via per-case latent event
   times (`T_onset`, `T_admit`, `T_death`, `T_disch`, `T_notif`,
   sampled where observed) within their respective day windows.
   Latents are shared across that case's delays, so the
   natural-history identity `D_oa + D_ad = D_od` holds *per case*
   for every posterior draw.
2. Stratified case-fatality: Bernoulli outcome with a logistic link
   on HCW status, case definition (Probable vs Confirmed), and
   standardised age.

Three parametric families are supported (`:lognormal`, `:gamma`,
`:weibull`), selected via the `family` keyword to `bdbv_model`. The
canonical run uses Gamma — lowest WAIC.

`bdbv_model_stratified(d; family)` adds a `β_*_hcw` log-mean shift to
each delay so HCW and non-HCW cases share the shape parameter but can
differ on the central tendency. (This variant still uses the
marginalised `double_interval_censored` formulation; see
[LIMITATIONS.md](LIMITATIONS.md).)

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

The model fits only the four atomic components — the natural-history
identity (`onset → death = onset → admit + admit → death`) is
satisfied per case by construction, so fitting `onset → death`
separately would just duplicate the information already in the
shared latents.

| Symbol | Component                | n  |
|---     |---                       |--- |
| `d_oa` | onset → admission        | 40 |
| `d_ad` | admission → death        | 22 |
| `d_ac` | admission → discharge    | 15 |
| `d_on` | onset → notification     | 38 |

`d_od` and `d_oc` are derived in post-processing as sample-level
convolutions of the population atomic distributions. The mean of the
convolved marginal equals the sum of the atomic means by linearity
of expectation; medians and quantiles come from Monte Carlo (500
samples per posterior draw).

## In-hospital length of stay via mixture post-processing

The length of stay in hospital (admission → departure) is the time an
admitted case occupies a bed before leaving, by either death or
discharge. It is derived in post-processing as a two-component mixture
of the fitted `d_ad` (admission → death) and `d_ac`
(admission → discharge), with the mixing weight set by the in-hospital
fatality among admitted cases:

```
d_los = w · d_ad + (1 − w) · d_ac,   w ~ Beta(1 + n_died, 1 + n_discharged)
```

`w` is the conjugate posterior on the admitted-case fatality fraction
(`n_died = 22`, `n_discharged = 15`), drawn once per posterior draw so
the marginal propagates uncertainty in both the delay distributions and
the fatality share. Sampled per draw (500 realisations) like the
convolved marginals. The fatal and survivor pathways are reported
separately alongside the overall mixture.

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

For each case `i`, sample the per-case event-time latents (where
observed) in reverse-chain order with **Sam Abbott's
bounded-primary trick** — the secondary event time bounds the
primary's upper window edge directly, so the support is smoothly
parametrised and NUTS doesn't see the wedge-shaped corner that the
naive `T_admit ≤ T_death` constraint would produce at same-day
cases:

```julia
T_death ~ Uniform(day, day + 1)                   # leaf event
T_disch ~ Uniform(day, day + 1)
T_admit ~ Uniform(day, min(day + 1, T_death, T_disch))
T_notif ~ Uniform(day, day + 1)
T_onset ~ Uniform(day, min(day + 1, T_admit, T_notif))
```

For each ordered pair both observed in the case, add the delay
likelihood `logpdf(dist_*, t_later − t_earlier)`.

Each bounded prior also carries a `log(upper − L)` Jacobian, which
restores the implicit independent-uniform-over-day-window prior of
the equivalent marginalised double-interval-censoring model. The
Jacobian vanishes for multi-day cases (where the upper bound is
just the natural day-window edge) and only contributes when the
ordering constraint binds (6 same-day cases here).

See Park *et al.* 2024 (medRxiv
[2024.01.12.24301247](https://doi.org/10.1101/2024.01.12.24301247))
§2.3.3 for the latent-variable formulation that this builds on.

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

For the in-hospital length of stay (`d_los`) the script reports the
fatal pathway (`d_ad`), the survivor pathway (`d_ac`), and the overall
mixture (median, mean, SD, 95th percentile with 95% CrI), along with the
in-hospital fatality weight. Supplied only when `summarise` is given the
model input `d`.

For the CFR block the script reports the log-OR coefficients and the
implied marginal CFR at each combination of HCW status and case
definition, with age held at its sample mean.
