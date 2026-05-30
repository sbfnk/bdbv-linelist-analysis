# Model

Two blocks fitted to the 2012 Isiro BDBV line list (n = 52):

1. Four atomic delay components, fitted via per-case latent event
   times ($T_\mathrm{onset}$, $T_\mathrm{admit}$, $T_\mathrm{death}$,
   $T_\mathrm{disch}$, $T_\mathrm{notif}$, sampled where observed)
   within their respective day windows. Latents are shared across
   that case's delays, so the natural-history identity
   $D_\mathrm{oa} + D_\mathrm{ad} = D_\mathrm{od}$ holds *per case*
   for every posterior draw.
2. Stratified case-fatality: Bernoulli outcome with a logistic link
   on HCW status, case definition (Probable vs Confirmed), and
   standardised age.

Three parametric families are supported (`:lognormal`, `:gamma`,
`:weibull`), selected via the `family` keyword to `bdbv_model`. The
canonical run uses Gamma — lowest WAIC.

`bdbv_model_stratified(d; family)` adds a per-delay log-mean shift
$\beta^{\mathrm{hcw}}_*$ so HCW and non-HCW cases share the shape
parameter but can differ on the central tendency. (This variant
still uses the marginalised `double_interval_censored` formulation;
see [LIMITATIONS.md](LIMITATIONS.md).)

## Common parametrisation across families

All three families share a log-location / log-shape parametrisation
so the priors are comparable. The latent log-location $\log\mu$
(LogNormal: log-median $\log m$) and log-shape $\log\alpha$ map to
the natural parameters $\mu = \exp(\log\mu)$ (or $m$) and $\alpha =
\exp(\log\alpha)$:

| Family    | Latent pair              | Distribution constructed                                           |
|---        |---                       |---                                                                 |
| LogNormal | $(\log m,\ \log\sigma)$  | $\mathrm{LogNormal}(\log m,\ \sigma)$                              |
| Gamma     | $(\log\mu,\ \log\alpha)$ | $\mathrm{Gamma}(\alpha,\ \mu/\alpha)$                              |
| Weibull   | $(\log\mu,\ \log\alpha)$ | $\mathrm{Weibull}\!\left(\alpha,\ \mu / \Gamma(1 + 1/\alpha)\right)$ |

Here $\mathrm{Gamma}(\alpha,\theta)$ uses the shape-scale
parametrisation (mean $\alpha\theta$); $\mathrm{Weibull}(\alpha,
\theta)$ uses shape $\alpha$ and scale $\theta$ (mean $\theta\,
\Gamma(1 + 1/\alpha)$). For LogNormal, $m$ is the median and
$\sigma$ the log-scale standard deviation.

For Weibull, $\log\alpha$ is truncated to $(-1.5,\ 1.5)$ so $\Gamma(1 +
1/\alpha)$ stays well-defined. The Weibull stratified model is not
supported.

## Compound delays via convolution post-processing

The model fits only the four atomic components — the natural-history
identity $D_\mathrm{od} = D_\mathrm{oa} + D_\mathrm{ad}$ (and
likewise $D_\mathrm{oc} = D_\mathrm{oa} + D_\mathrm{ac}$) holds per
case by construction, so fitting onset → death separately would just
duplicate the information already in the shared latents.

| Symbol            | Component                | $n$ |
|---                |---                       |---  |
| $D_\mathrm{oa}$   | onset → admission        | 40  |
| $D_\mathrm{ad}$   | admission → death        | 22  |
| $D_\mathrm{ac}$   | admission → discharge    | 15  |
| $D_\mathrm{on}$   | onset → notification     | 38  |

The marginal compound distributions
$D_\mathrm{od} = D_\mathrm{oa} * D_\mathrm{ad}$ and
$D_\mathrm{oc} = D_\mathrm{oa} * D_\mathrm{ac}$ ($*$ denoting
convolution) are derived in post-processing as sample-level
convolutions of the population atomic distributions. The mean of
the convolved marginal equals the sum of the atomic means by
linearity of expectation; medians and quantiles come from Monte
Carlo (500 samples per posterior draw).

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

Weakly informative, centred on a plausible day count for each delay.
The `prior_scale` keyword (default $s = 1$) scales the SDs jointly:

```math
\begin{aligned}
\log m_{\mathrm{oa}} \ \text{or}\ \log\mu_{\mathrm{oa}} &\sim \mathrm{Normal}(\log 3,\ s) \\
\log m_{\mathrm{ad}} \ \text{or}\ \log\mu_{\mathrm{ad}} &\sim \mathrm{Normal}(\log 6,\ s) \\
\log m_{\mathrm{ac}} \ \text{or}\ \log\mu_{\mathrm{ac}} &\sim \mathrm{Normal}(\log 13,\ s) \\
\log m_{\mathrm{on}} \ \text{or}\ \log\mu_{\mathrm{on}} &\sim \mathrm{Normal}(\log 7,\ s)
\end{aligned}
```

with one scale-parameter prior per family,

```math
\begin{aligned}
\log\sigma                       &\sim \mathrm{Normal}_{+}(0,\ s)                 && \text{(LogNormal)} \\
\log\alpha                       &\sim \mathrm{Normal}(0,\ s)                     && \text{(Gamma)} \\
\log\alpha\,|\,(-1.5,\,1.5)      &\sim \mathrm{Normal}(0,\ s)\ \text{truncated}   && \text{(Weibull)}
\end{aligned}
```

and CFR coefficients

```math
\beta_0 \sim \mathrm{Normal}(0,\ 2), \quad
\beta_{\mathrm{hcw}},\ \beta_{\mathrm{def}},\ \beta_{\mathrm{age}}
   \sim \mathrm{Normal}(0,\ 1), \quad
\beta^{\mathrm{hcw}}_{*} \sim \mathrm{Normal}(0,\ 0.5),
```

where the last is the per-delay HCW log-mean shift in the stratified
model only.

The $s = 1$ default on the log-medians gives roughly threefold prior
latitude either way on the central tendency. A 3-scale sensitivity
sweep ($s \in \{0.5,\ 1.0,\ 2.0\}$) shifts posterior means by less
than 5%.

## Likelihood

For each case $i$ with reported day $d^{(i)}_e$ for event $e \in
\{\mathrm{onset}, \mathrm{admit}, \mathrm{death}, \mathrm{disch},
\mathrm{notif}\}$, sample the within-day latent times $T^{(i)}_e$ in
reverse-chain order with **Sam Abbott's bounded-primary trick** —
the secondary event time bounds the primary's upper window edge
directly, so the support is smoothly parametrised and NUTS doesn't
see the wedge-shaped corner that the naive ordering constraint
$T^{(i)}_{\mathrm{admit}} \le T^{(i)}_{\mathrm{death}}$ would produce
at same-day cases. For each case, with the (possibly-tightened)
upper edge $U^{(i)}_e$,

```math
\begin{aligned}
T^{(i)}_{\mathrm{death}} &\sim \mathrm{Uniform}\!\left(d^{(i)}_{\mathrm{death}},\ d^{(i)}_{\mathrm{death}} + 1\right) \\
T^{(i)}_{\mathrm{disch}} &\sim \mathrm{Uniform}\!\left(d^{(i)}_{\mathrm{disch}},\ d^{(i)}_{\mathrm{disch}} + 1\right) \\
T^{(i)}_{\mathrm{admit}} &\sim \mathrm{Uniform}\!\left(d^{(i)}_{\mathrm{admit}},\ U^{(i)}_{\mathrm{admit}}\right),
   &U^{(i)}_{\mathrm{admit}} &= \min\!\left(d^{(i)}_{\mathrm{admit}} + 1,\ T^{(i)}_{\mathrm{death}},\ T^{(i)}_{\mathrm{disch}}\right) \\
T^{(i)}_{\mathrm{notif}} &\sim \mathrm{Uniform}\!\left(d^{(i)}_{\mathrm{notif}},\ d^{(i)}_{\mathrm{notif}} + 1\right) \\
T^{(i)}_{\mathrm{onset}} &\sim \mathrm{Uniform}\!\left(d^{(i)}_{\mathrm{onset}},\ U^{(i)}_{\mathrm{onset}}\right),
   &U^{(i)}_{\mathrm{onset}} &= \min\!\left(d^{(i)}_{\mathrm{onset}} + 1,\ T^{(i)}_{\mathrm{admit}},\ T^{(i)}_{\mathrm{notif}}\right),
\end{aligned}
```

where any latent whose event is unobserved for case $i$ is dropped.

For each ordered pair $(e_1, e_2)$ both observed in case $i$, the
delay likelihood contributes

```math
\ell^{(i)}_{e_1 e_2} \;=\; \log f_{e_1 e_2}\!\left(T^{(i)}_{e_2} - T^{(i)}_{e_1}\right),
```

where $f_{e_1 e_2}$ is the density of the fitted delay distribution
$D_{e_1 e_2}$ for that ordered pair (the four atomic components
$D_\mathrm{oa}, D_\mathrm{ad}, D_\mathrm{ac}, D_\mathrm{on}$).

Each bounded prior carries a Jacobian $\log\!\left(U^{(i)}_e -
d^{(i)}_e\right)$, which restores the implicit
independent-uniform-over-day-window prior of the equivalent
marginalised double-interval-censoring model. The Jacobian vanishes
for multi-day cases (where $U^{(i)}_e = d^{(i)}_e + 1$) and only
contributes when the ordering constraint binds (6 same-day cases
here).

See Park *et al.* 2024 (medRxiv
[2024.01.12.24301247](https://doi.org/10.1101/2024.01.12.24301247))
§2.3.3 for the latent-variable formulation that this builds on.

**CFR block.** For each case $i$ with outcome $y_i \in \{0, 1\}$,
HCW indicator $h_i$, probable-case indicator $p_i$, and standardised
age $z_i$,

```math
\eta_i \;=\; \beta_0 + \beta_{\mathrm{hcw}}\,h_i + \beta_{\mathrm{def}}\,p_i + \beta_{\mathrm{age}}\,z_i,
\qquad
y_i \;\sim\; \mathrm{Bernoulli}\!\left(\operatorname{logistic}(\eta_i)\right).
```

The standardised age $z_i$ has sample mean 0 and SD 1, with the
population mean imputed for the 3 of 52 missing ages.
$\operatorname{logistic}(\eta)$ is clamped to $(10^{-10},\ 1 -
10^{-10})$ so the Bernoulli's domain check doesn't fail when NUTS
proposes large $|\eta|$ during warmup.

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
and $\mathrm{WAIC} = -2\,(\mathrm{lppd} - p_{\mathrm{waic}})$
aggregated across the 167 observations ($40 + 22 + 15 + 38$ atomic
delays $+\ 52$ outcomes).

PSIS-LOO would be preferable but adds dependencies; WAIC is adequate
here.

## Post-processing

For each atomic delay the script reports the posterior median, mean,
and (for LogNormal) log-SD with 95% CrI. Gamma and Weibull medians
come from the per-draw inverse-CDF at 0.5.

For each derived marginal
($D_\mathrm{od} = D_\mathrm{oa} * D_\mathrm{ad}$ and
$D_\mathrm{oc} = D_\mathrm{oa} * D_\mathrm{ac}$) the script samples
500 realisations per posterior draw from the convolution and
reports median, mean, SD, and 95th percentile with 95% CrI across
draws.

For the in-hospital length of stay ($D_\mathrm{los}$) the script
reports the fatal pathway ($D_\mathrm{ad}$), the survivor pathway
($D_\mathrm{ac}$), and the overall mixture (median, mean, SD, 95th
percentile with 95% CrI), along with the in-hospital fatality
weight. Supplied only when `summarise` is given the model input `d`.

For the CFR block the script reports the log-OR coefficients and the
implied marginal CFR at each combination of HCW status and case
definition, with age held at its sample mean.
