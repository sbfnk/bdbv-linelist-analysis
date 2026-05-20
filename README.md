# Bundibugyo Ebola virus — Bayesian delay-distribution and stratified CFR estimation from the 2012 Isiro outbreak

[![Docs (dev)](https://img.shields.io/badge/docs-dev-blue.svg)](https://sbfnk.github.io/bdbv-linelist-analysis/dev)

Full documentation — write-up with all tables, model description, API
reference — at <https://sbfnk.github.io/bdbv-linelist-analysis/dev>.

A Julia + Turing + [CensoredDistributions.jl](https://censoreddistributions.epiaware.org/)
re-analysis of the only publicly accessible Bundibugyo ebolavirus (BDBV)
line list — the 2012 Isiro outbreak in Haut-Uélé, DRC, as deposited by
[Rosello *et al.* 2015 in *eLife*](https://elifesciences.org/articles/09015)
— following the best-practice checklist of
[Charniga *et al.* 2024 in PLOS Comput Biol](https://doi.org/10.1371/journal.pcbi.1012520).

## Key finding

The onset-to-notification mean of 8.83 d listed in
[Rosello *et al.* 2015 Table 5](https://elifesciences.org/articles/09015)
(and consequently propagated through PERG, WHO grEPI and other
downstream parameter repositories) derives from a 30-day cap on the
underlying dates. The uncapped Gamma fit gives a mean of 19.7 d
(95% CrI 13.7 – 30.1) — a factor of ≈ 2.2 longer.

## Results

Rendered walkthrough with all tables and figures regenerated from the
current model on every push to `main`:
<https://sbfnk.github.io/bdbv-linelist-analysis/dev/analysis>.

Raw artefacts from the most recent main build (posterior CSV and all figures):
<https://github.com/sbfnk/bdbv-linelist-analysis/releases/tag/main-latest>.

## Methods and limitations

Model description and priors are in [MODEL.md](MODEL.md).
Known caveats are in [LIMITATIONS.md](LIMITATIONS.md).

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

## Authors

Sebastian Funk (London School of Hygiene & Tropical Medicine).

## License

MIT (see [LICENSE](LICENSE)).
