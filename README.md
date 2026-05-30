# Bundibugyo Ebola virus — Bayesian delay-distribution and stratified CFR estimation from the 2012 Isiro outbreak

[![Docs](https://img.shields.io/badge/docs-dev-blue.svg)](https://epiforecasts.io/bdbv-linelist-analysis/dev)

A Julia + Turing + [CensoredDistributions.jl](https://censoreddistributions.epiaware.org/)
re-analysis of the only publicly accessible Bundibugyo ebolavirus (BDBV)
line list — the 2012 Isiro outbreak in Haut-Uélé, DRC, as deposited by
[Rosello *et al.* 2015 in *eLife*](https://elifesciences.org/articles/09015)
— following the best-practice checklist of
[Charniga *et al.* 2024 in PLOS Comput Biol](https://doi.org/10.1371/journal.pcbi.1012520).

🌐 **[Analysis walkthrough](https://epiforecasts.io/bdbv-linelist-analysis/dev/analysis)** — full tables, figures, and diagnostics regenerated from the current model on every push to `main` (HTML)

⚠️ **[Limitations](https://epiforecasts.io/bdbv-linelist-analysis/dev/limitations)** — read before using these estimates: data, model, inference, and generalisability caveats (HTML)

✅ **[Charniga 2024 checklist](https://epiforecasts.io/bdbv-linelist-analysis/dev/charniga-checklist)** — best-practice item-by-item compliance table (HTML)

📦 **[Posterior CSV and figures](https://github.com/epiforecasts/bdbv-linelist-analysis/releases/tag/main-latest)** — `main-latest` rolling release bundle

📖 **[Model description](https://epiforecasts.io/bdbv-linelist-analysis/dev/model)** — priors, likelihood, inference (HTML)

📑 **[API reference](https://epiforecasts.io/bdbv-linelist-analysis/dev/api)** — exported functions (HTML)

## Headline estimates

<!-- HEADLINE:START -->
Live headline estimates (Gamma fit, four atomic delay components, with Rosello 2015 means side by side) are rendered on the [docs home page](https://epiforecasts.io/bdbv-linelist-analysis/dev) — regenerated from the current model on every push to `main`.
<!-- HEADLINE:END -->

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

A minute or two per family on a laptop. Posterior CSV is written
to `output/` and figures to `figures/` locally; both directories
are gitignored. The rolling
[`main-latest` release](https://github.com/epiforecasts/bdbv-linelist-analysis/releases/tag/main-latest)
publishes the canonical posterior CSV and figure bundle regenerated
on every push to `main`, and the [docs site](https://epiforecasts.io/bdbv-linelist-analysis/dev)
renders the same figures inline.

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

Sebastian Funk and Sam Abbott (London School of Hygiene & Tropical Medicine).

## Repository administration

Branch protection on `main` is configured via
[`scripts/setup-branch-protection.sh`](scripts/setup-branch-protection.sh).
Run it once with an admin-scoped `gh` token to require pull requests,
passing status checks, conversation resolution, and to block force-pushes
and deletion on `main`.

## License

MIT (see [LICENSE](https://github.com/epiforecasts/bdbv-linelist-analysis/blob/main/LICENSE)).
