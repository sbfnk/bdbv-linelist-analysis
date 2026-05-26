# Charniga *et al.* 2024 checklist

[Charniga K, Park SW, Akhmetzhanov AR *et al.* "Best practices for
estimating and reporting epidemiological delay distributions of
infectious diseases." *PLoS Comput Biol* 2024;20(10):e1012520.](https://doi.org/10.1371/journal.pcbi.1012520)

The table below maps each Table 2 checklist item from Charniga *et al.*
2024 to where (and how) it is addressed in this analysis. ✅ = done,
🟡 = partial / done with documented caveats, ❌ = not done. We have
tried to be honest about partial coverage rather than claim more than
the deposit supports.

## Estimation

| # | Charniga 2024 item | Status | Where in this repo | Note |
|---|---|---|---|---|
| 1 | Adjust for censoring (always), right truncation (when needed), and dynamical bias (when needed); state which adjustments were made | 🟡 | [`src/model.jl`](https://github.com/epiforecasts/bdbv-linelist-analysis/blob/main/src/model.jl) and [Model](model.md#Likelihood) | Double interval censoring at 1-day bins on both endpoints via `CensoredDistributions.double_interval_censored`. No right-truncation or dynamical-bias adjustment because the 2012 Isiro outbreak is fitted retrospectively in 2026 on the closed Rosello deposit — neither bias applies. |
| 2 | Fit more than one probability distribution and use model-selection criteria; visualise the fit | ✅ | [Analysis walkthrough — Family comparison](analysis.md#Family-comparison) | LogNormal, Gamma, and Weibull all fitted; WAIC ranking + side-by-side posterior-predictive panels per family. Gamma wins. |
| 3 | Convert distribution parameters to summary statistics correctly | ✅ | [Model — Common parametrisation](model.md#Common-parametrisation-across-families) | All families share a log-mean / log-shape parametrisation; conversion uses `Distributions.jl` built-ins (`Gamma(k, μ/k)`, `Weibull(α, μ/Γ(1+1/α))`, `LogNormal(log_median, log_sd)`). Documented in a single table. |
| 4 | Stratify or add subgroups when sample size allows and group differences are hypothesised | 🟡 | [`bdbv_model_stratified`](api.md) and [Model](model.md) | HCW stratification implemented for the delay block (`β_*_hcw` log-mean shifts on Gamma); HCW, case definition (Probable vs Confirmed), and standardised age in the CFR block. Sex is recorded but not stratified — [Limitations](limitations.md#Model) explains why (Kratz 2015's community-vs-ETC stratum is not reproducible from the Rosello deposit). |
| 5 | Check model diagnostics (R̂, divergent transitions, ESS) | ✅ | [Analysis walkthrough — Family comparison](analysis.md#Family-comparison) | Per-family table reports max R̂, min bulk ESS, and divergent transition count. Across all three families: max R̂ ≤ 1.002, min ESS ≥ 3,200, 0 divergent transitions. |

## Reporting

| # | Charniga 2024 item | Status | Where in this repo | Note |
|---|---|---|---|---|
| 6 | Report measures of central tendency and variability (mean **and** SD/variance/dispersion) | ✅ | [Analysis walkthrough — Headline estimates](analysis.md#Headline-Gamma-estimates) and [— Derived (convolved) marginals](analysis.md#Derived-(convolved)-marginals) | Median, mean **and SD** reported for every atomic and convolved delay (the convolved-marginals table gained an `sd (95% CrI)` column once `_convolve_delays` started returning per-draw SDs). Per-draw shape/scale/SD also exposed in `posterior_<family>.csv`. |
| 7 | Report key quantiles (e.g. 2.5, 5, 25, 50, 75, 95, 97.5, 99) of the distribution | 🟡 | [Analysis walkthrough](analysis.md) | The 50th percentile is reported as the median; the 95th percentile is reported for the convolved onset → death and onset → discharge marginals. A full quantile grid (5/25/75/95/99) is not produced — `output/` posteriors allow it to be computed downstream. |
| 8 | Report the fitted-distribution parameters | 🟡 | [`output/`](https://github.com/epiforecasts/bdbv-linelist-analysis/tree/main/output) posterior CSV; parametrisation table in [Model](model.md#Common-parametrisation-across-families) | The on-page tables emphasise summary statistics (median, mean, 95% CrI). Latent log-mean / log-shape draws are written to the `main-latest` rolling-release CSV bundle so the underlying Gamma `(shape, scale)`, LogNormal `(logmean, logsd)`, or Weibull `(shape, scale)` can be reconstructed exactly. |
| 9 | Report uncertainty (90% or 95% intervals) on every estimate | ✅ | [Analysis walkthrough](analysis.md) throughout | All headline numbers, derived marginals, CFR strata, and log-OR coefficients are reported as `median (2.5% – 97.5%)`. |
| 10 | Report study-sample characteristics (sample size, age, sex, location, route of exposure, vaccination status if any) | ✅ | [Analysis walkthrough — Outbreak context](analysis.md#Outbreak-context); [README — Data provenance](https://github.com/epiforecasts/bdbv-linelist-analysis#data-provenance); [Limitations — Data](limitations.md#Data) | The walkthrough opens with an "Outbreak context" block tabulating n, onset window, age (median + IQR), sex split, HCW split, case-definition mix, and outcome counts directly from `data/linelist.csv`. Route of exposure is not recorded in the deposit; no vaccine was available in 2012 (vaccination status not applicable). |
| 11 | Report the epidemic curve and any control measures in place | 🟡 | [Analysis walkthrough — Epidemic curve](analysis.md#Epidemic-curve) | Weekly onset epi curve with HCW subcounts is rendered live on the docs site. Control measures are mentioned narratively in [Limitations](limitations.md#Generalisability) (MSF ETC support; rural mining-influenced setting) rather than as a structured timeline. |
| 12 | Provide anonymised / de-identified line-list data and documented code | ✅ | [`data/linelist.csv`](https://github.com/epiforecasts/bdbv-linelist-analysis/blob/main/data/linelist.csv); [`src/`](https://github.com/epiforecasts/bdbv-linelist-analysis/tree/main/src); [`output/` rolling release](https://github.com/epiforecasts/bdbv-linelist-analysis/releases/tag/main-latest) | The bundled CSV is Rosello *et al.* 2015's supplementary file 1 (CC-BY 4.0) subset to BDBV — already anonymised at source. All Julia code (model, post-processing, plotting) is MIT-licensed in this repo. The full executable walkthrough is published as the docs site, and posterior CSVs + figures ship as a `main-latest` GitHub release. |

## Items 13+ (Table 3 — incubation period & serial interval)

Charniga *et al.* 2024 Table 3 adds recommendations specific to
incubation-period and serial-interval estimation (exposure-window
encoding, transmission-pair confidence, transmission-direction
ordering, support for negative serial intervals).

**Not applicable here.** The Rosello deposit has no exposure dates and
no transmission pairs, so neither the incubation period nor the serial
interval is estimated in this work. See
[Limitations — Data](limitations.md#Data) for the rationale and the
MacNeil 2010 / Zaire-EVD priors recommended for downstream nowcasting
or R<sub>t</sub> work that needs these distributions.

## How to update this page

When the model or reporting changes, edit the table directly. Status
should err on the side of 🟡 over ✅ whenever a documented caveat
exists. The LLM-friendly column layout (Charniga item ↔ status ↔
location ↔ note) is deliberate — it lets future automated reviews
re-verify each row by following the link in the third column.
