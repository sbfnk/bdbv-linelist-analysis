```@meta
CurrentModule = BdbvLinelist
```

# Bundibugyo Ebola virus — Bayesian delay-distribution and stratified CFR estimation from the 2012 Isiro outbreak

A Julia + Turing + [CensoredDistributions.jl](https://censoreddistributions.epiaware.org/) re-analysis of the only publicly accessible Bundibugyo ebolavirus (BDBV) line list — the 2012 Isiro outbreak in Haut-Uélé, DRC, as deposited by [Rosello *et al.* 2015 in *eLife*](https://elifesciences.org/articles/09015) — following the best-practice checklist of [Charniga *et al.* 2024 in PLOS Comput Biol](https://doi.org/10.1371/journal.pcbi.1012520).

## Key finding

The onset-to-notification mean of 8.83 d listed in [Rosello *et al.* 2015 Table 5](https://elifesciences.org/articles/09015) (and consequently propagated through PERG, WHO grEPI and other downstream parameter repositories) derives from a 30-day cap on the underlying dates. The uncapped Gamma fit gives a mean of 19.7 d (95% CrI 13.7 – 30.1) — a factor of ≈ 2.2 longer. The same cap is applied across at least four of the Rosello outbreaks; the published Isiro point estimate is downward-biased by roughly the same factor in all of them.

## Pages

- [Model](model.md) — likelihood, priors, parametrisation, inference, post-processing.
- [Limitations](limitations.md) — known caveats around data, model, inference, and generalisability.
- [Analysis walkthrough](analysis.md) — executable walkthrough that regenerates the headline tables and figures from the fitted model at build time.
- [API Reference](api.md) — exported functions.

Raw artefacts (posterior CSV and figures) from the most recent main build live on the [`main-latest` release](https://github.com/sbfnk/bdbv-linelist-analysis/releases/tag/main-latest).

## Citing

> Rosello A, Mossoko M, Flasche S, *et al.* Ebola virus disease in the Democratic Republic of the Congo, 1976–2014. *eLife* 2015;4:e09015. [doi:10.7554/eLife.09015](https://doi.org/10.7554/eLife.09015)

The reporting follows:

> Charniga K, Park SW, Akhmetzhanov AR, *et al.* Best practices for estimating and reporting epidemiological delay distributions of infectious diseases. *PLoS Comput Biol* 2024;20(10):e1012520. [doi:10.1371/journal.pcbi.1012520](https://doi.org/10.1371/journal.pcbi.1012520)
