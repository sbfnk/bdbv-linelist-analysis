## Regenerate the canonical posterior + all figures for the rolling
## `main-latest` GitHub release.
##
## Fits all three families (so plot_family_comparison has data),
## saves the canonical Gamma posterior to output/posterior_gamma.csv,
## and writes ppc_gamma.png, ppc_family_comparison.png and
## epi_curve.png into figures/.

using BdbvLinelist

ll = load_linelist()
d  = build_data(ll)

results = compare_families()
chains_by_family = Dict(f => results[f].chain for f in (:lognormal, :gamma, :weibull))

post_gamma = summarise(results[:gamma].chain, :gamma)
save_posterior(post_gamma, joinpath("output", "posterior_gamma.csv"))

save_figure(plot_epi_curve(ll), joinpath("figures", "epi_curve.png"))
save_figure(plot_ppc(results[:gamma].chain, d, :gamma),
            joinpath("figures", "ppc_gamma.png"))
save_figure(plot_family_comparison(chains_by_family, d),
            joinpath("figures", "ppc_family_comparison.png"))
