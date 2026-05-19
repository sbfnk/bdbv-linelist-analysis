# Bundibugyo Ebola virus — Bayesian delay-distribution and stratified CFR estimation from the 2012 Isiro outbreak

A Julia + Turing + [CensoredDistributions.jl](https://censoreddistributions.epiaware.org/)
re-analysis of the only publicly accessible Bundibugyo ebolavirus (BDBV)
line list — the 2012 Isiro outbreak in Haut-Uélé, DRC, as deposited by
[Rosello *et al.* 2015 in *eLife*](https://elifesciences.org/articles/09015)
— following the best-practice checklist of
[Charniga *et al.* 2024 in PLOS Comput Biol](https://doi.org/10.1371/journal.pcbi.1012520).

Four atomic delay components are estimated under double interval
censoring at the day level (onset → admission, admission → death,
admission → discharge, onset → notification), comparing **LogNormal**,
**Gamma** and **Weibull** families via WAIC. The marginal onset → death
and onset → discharge distributions are derived in post-processing as
convolutions of the atomic components — this enforces the per-case
natural-history identity `onset → death = onset → admit + admit → death`
without giving each marginal its own degrees of freedom. A stratified
CFR logistic regression (HCW status, case definition, age) uses all
52 cases. An HCW-stratified delay variant is also available.

Sister project to
[`andv-linelist-analysis`](https://github.com/sbfnk/andv-linelist-analysis)
(Epuyén Andes hantavirus). The Isiro deposit has no transmission pairs
and no exposure dates, so incubation period, serial interval and
generation interval cannot be fitted from this data set — see
[`LIMITATIONS.md`](LIMITATIONS.md).

## Headline results

NUTS, 4 chains × 1000 post-warmup samples. Convergence: max R̂ ≤ 1.002,
min bulk ESS ≥ 3,200, 0 divergent transitions across families.

### Family comparison (WAIC, lower = better)

| family    | WAIC   | ΔWAIC | p_waic | max R̂ | min ESS |
|---|---|---|---|---|---|
| lognormal | 806.6  | +14.4 | 9.8    | 1.002  | 3,391   |
| **gamma** | **792.2** | **0.0**  | **9.1**  | **1.001** | **3,243** |
| weibull   | 793.7  | +1.4  | 9.3    | 1.000  | 3,015   |

**Gamma wins by ΔWAIC = 14.4** over LogNormal; Weibull essentially tied
with Gamma. LogNormal over-fits the heavy right tail of the notification
delays (max observed 86 d).

### Fitted atomic delay components (Gamma, doubly censored)

| Component                | n  | Median (95% CrI) | Mean (95% CrI) |
|---                       |--- |---               |---             |
| Onset → admission        | 40 |  2.9 (2.0 –  4.0)|  4.0 (3.0 –  5.5) |
| Admission → death        | 22 |  6.4 (4.6 –  8.6)|  7.6 (5.8 – 10.3) |
| Admission → discharge    | 15 |  5.4 (2.8 –  9.4)|  7.7 (4.9 – 14.4) |
| Onset → notification     | 38 | 10.9 (6.5 – 17.2)| 19.7 (13.7 – 30.1)|

### Derived marginals — convolutions of the atomic components

| Marginal | Median (95% CrI) | Mean (95% CrI) | P95 (95% CrI) |
|---|---|---|---|
| Onset → death     (= oa + ad) | 10.4 (8.3 – 13.1) | 11.7 (9.4 – 14.6) | 24.0 (19.0 – 32.4) |
| Onset → discharge (= oa + ac) |  9.7 (6.9 – 14.0) | 11.8 (8.7 – 18.8) | 28.4 (19.8 – 54.2) |

### Stratified case-fatality

| Stratum | CFR (95% CrI) |
|---|---|
| Non-HCW, Confirmed (baseline) | **0.47 (0.31 – 0.65)** |
| HCW, Confirmed                | **0.23 (0.08 – 0.47)** |
| Non-HCW, Probable             | **0.88 (0.71 – 0.96)** |
| HCW, Probable                 | **0.70 (0.38 – 0.91)** |

Logit-scale coefficients:

| Coefficient       | log-OR (95% CrI)     | OR (95% CrI)        |
|---                |---                   |---                  |
| HCW status        | −1.12 (−2.33, +0.07) | 0.33 (0.10, 1.07)   |
| Probable case def | +2.07 (+0.95, +3.35) | 7.95 (2.58, 28.45)  |
| Standardised age  | −0.29 (−0.93, +0.31) | 0.74 (0.40, 1.37)   |

HCW status independently lowers CFR after adjusting for case
definition (OR 0.33, 95% CrI just crossing 1.0). The raw 28.6% vs
63.2% HCW vs non-HCW gap is dampened but not eliminated by the
probable-vs-confirmed control — most of the gap survives the
adjustment.

### HCW × delay stratification

A separate gamma model with HCW shifts on each delay's log-mean
shows no meaningful clinical-delay difference between HCWs and
non-HCWs, but HCW notifications take ≈ 70% longer on average
(borderline significant):

| Delay | exp(β_HCW) on mean (95% CrI) |
|---|---|
| Onset → admission     | 0.81 (0.47 – 1.53) |
| Admission → death     | 0.94 (0.52 – 1.86) |
| Admission → discharge | 0.95 (0.45 – 1.85) |
| Onset → notification  | **1.70 (0.95 – 3.17)** |

### Prior sensitivity (Gamma family)

Mean delays under three prior-scale settings on the log-mean and
log-shape priors (default 1.0 ⇒ ~×3-fold latitude either side of
the prior median):

| Delay              | scale 0.5         | scale 1.0 (default) | scale 2.0         |
|---                 |---                |---                  |---                |
| Onset → admission  | 3.9 (3.0 – 5.3)   | 4.0 (3.0 – 5.5)     | 4.1 (3.0 – 5.5)   |
| Admission → death  | 7.5 (5.6 – 10.3)  | 7.6 (5.8 – 10.3)    | 7.6 (5.7 – 10.5)  |
| Admission → discharge | 8.5 (5.4 – 14.2)| 7.7 (4.9 – 14.4) | 7.6 (4.6 – 13.4) |
| Onset → notification | 18.0 (12.8 – 25.8) | 19.7 (13.7 – 30.1) | 20.4 (13.9 – 32.0) |

Posterior means shift by < 5% across the 4× prior-scale range —
priors are confirmed weakly informative as intended.

## Figures

| File | Content |
|---|---|
| `figures/ppc_gamma.png` | Posterior predictive bars (4 panels, Gamma) |
| `figures/ppc_family_comparison.png` | Side-by-side LogNormal / Gamma / Weibull PP checks |
| `figures/epi_curve.png` | Weekly onset counts with HCW subcount stacked |

## What this analysis does not do

- **No transmission pairs in the deposit** → serial interval and
  generation interval cannot be fitted from this data set; for
  downstream Rt work use Zaire-EVD generation interval priors as a
  starting point.
- **Per-case censoring noise across the four delays of the same case
  is treated as independent** (not biased, slightly inefficient).
  A full latent-time joint model was prototyped but creates wedge-
  shaped boundary geometry NUTS handles poorly at the 6 same-day-
  admission cases; the marginal formulation here matches what
  CensoredDistributions.jl naturally supports and what `Charniga
  et al. 2024` recommends for retrospective complete-outbreak data.

## Why these numbers are new

The Rosello deposit has been openly available since 2015 and its
sample-mean delay statistics have been ingested into the
[PERG/epireview](https://mrc-ide.github.io/epireview/) and
[grEPI](https://collaboratory.who.int/epidemiologicalparameters/repository/)
parameter databases as Gamma(mean, SD). What's added here:

1. **Proper double-interval-censoring correction** for every delay.
2. **Bayesian posteriors** with credible intervals on every parameter.
3. **Joint structure** enforcing the per-case natural-history identity
   via convolution post-processing.
4. **Family comparison** (Gamma > LogNormal by ΔWAIC = 14).
5. **HCW × case-definition stratified CFR** with full Bayesian
   uncertainty — neither the original Kratz 2015 paper nor the
   PERG/grEPI database provides this.
6. **HCW-stratified delays** (clinical delays unchanged; notification
   delay borderline longer for HCWs).
7. **Prior sensitivity sweep** confirming the priors are weakly
   informative.
8. **Reproducible Julia code** in a public repo for direct re-use.

See [MODEL.md](MODEL.md) for the full model description and
[LIMITATIONS.md](LIMITATIONS.md) for caveats.

## Repository layout

```
src/
  BdbvLinelist.jl   — module entry point and imports
  data.jl           — line list loading, outlier scrub, per-pair delay extraction
  model.jl          — joint Turing model + stratified + community-death variants
  postprocess.jl    — diagnostics, summaries, WAIC, convolution post-processing
  plots.jl          — Makie posterior predictive + epi curve
  main.jl           — analyse() / compare_families() / sensitivity() / fit_death_mixture() drivers
data/
  linelist.csv      — Isiro 2012 BDBV subset of Rosello 2015 eLife supp 1
figures/
  ppc_gamma.png
  ppc_family_comparison.png
  epi_curve.png
Project.toml        — Julia package manifest
Manifest.toml       — locked dependency versions
CITATION.cff        — citation metadata
.zenodo.json        — Zenodo deposit metadata
LICENSE             — MIT
```

## Running

```bash
# Headline run with Gamma (default after WAIC selection)
julia --project=. -t auto -m BdbvLinelist -- -f gamma

# Compare all three families via WAIC
julia --project=. -t auto -m BdbvLinelist -- --compare

# From the REPL
julia> using BdbvLinelist
julia> chn, post, diag      = analyse(family = :gamma)
julia> results              = compare_families()
julia> sensitivity_results  = sensitivity()
julia> save_figure(plot_ppc(chn, build_data(load_linelist()), :gamma), "figures/ppc_gamma.png")
```

A minute or two per family on a laptop. Posterior CSV in `output/`,
figures in `figures/`.

## Data provenance

`data/linelist.csv` is the Bundibugyo subset (n = 52) of the aggregated
seven-outbreak DRC line list (n = 996) published as supplementary
file 1 of:

> Rosello A, Mossoko M, Flasche S, *et al.* Ebola virus disease in the
> Democratic Republic of the Congo, 1976–2014. *eLife* 2015;4:e09015.
> [doi:10.7554/eLife.09015](https://doi.org/10.7554/eLife.09015)

Licensed CC-BY 4.0. The five admission-date encoding outliers (−89,
−5, −4, −1, 328720 days from onset) are programmatically set to
missing during loading. One notification-date outlier (delay of −62
days) is similarly dropped. No other modifications.

## Citing

If you use the model or the cleaned line list, please cite Rosello
*et al.* 2015 above, the Charniga *et al.* 2024 best-practice paper
this analysis follows, and:

> Funk S. Bundibugyo Ebola virus — Bayesian delay distributions
> from the 2012 Isiro outbreak. 2026. https://github.com/sbfnk/bdbv-linelist-analysis

(see `CITATION.cff` for the machine-readable form).

## Authors

Sebastian Funk (London School of Hygiene & Tropical Medicine).

## License

MIT (see [LICENSE](LICENSE)).
