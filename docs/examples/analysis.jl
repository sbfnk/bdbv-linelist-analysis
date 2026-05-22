# # Analysis walkthrough
#
# This page fits the BDBV delay + CFR model in `BdbvLinelist.jl` to the
# Rosello *et al.* 2015 Isiro deposit (n = 52) and renders the headline
# outputs.
# Three parametric families (LogNormal, Gamma, Weibull) are compared by
# WAIC; the Gamma fit is the canonical run.
# The marginal onset → death and onset → discharge distributions are
# derived in post-processing as Monte Carlo convolutions of the atomic
# components, enforcing the per-case natural-history identity at the
# population level.
# A logistic regression stratifies case-fatality by HCW status, case
# definition (Probable vs Confirmed), and standardised age.
#
# Priors, the doubly-censored likelihood, and post-processing are
# detailed on the [Model](model.md) page.
# Caveats are on the [Limitations](limitations.md) page.

using BdbvLinelist
using DataFrames
using DataFramesMeta
using Chain
using FlexiChains
using Statistics
using Printf
using Random
using CairoMakie
using Turing: @varname

# ## Load the line list
#
# `load_linelist` parses the bundled CSV from `data/linelist.csv` and
# scrubs the five admission-date encoding outliers and the one
# notification-date outlier described in [Limitations](limitations.md).
# `build_data` re-encodes the per-pair delays into the vectors used by
# the doubly-censored likelihood.

ll = load_linelist()
d  = build_data(ll)

@chain ll begin
    @rsubset(!ismissing(:Date_of_onset_symp))
    @select(:Person_ID, :Date_of_onset_symp,
        :Date_of_Hospitalisation, :Date_of_Death,
        :Date_hospital_discharge, :Date_of_notification,
        :is_hcw, :Case_definition)
    first(8)
end

# ## Headline Gamma estimates
#
# Fit all three families (Gamma wins on WAIC; comparison table and
# plot are in the [Family comparison](#Family-comparison) section
# below) and tabulate the four atomic delay components alongside
# Rosello *et al.* 2015 Table 5 means. The Rosello fits applied a
# 30-day cap on the underlying dates — this is binding on
# onset → notification (where 9 of 38 cases exceed 30 d) and
# accounts for most of the difference there.

## Suppress the package's stdout printing — we build clean DataFrame
## tables from the returned posterior vectors below.
results = redirect_stdout(devnull) do
    compare_families()
end
nothing #hide

chn_gamma = results[:gamma].chain
post = redirect_stdout(devnull) do
    summarise(chn_gamma, :gamma)
end
nothing #hide

# Local helper: format a vector as `median (2.5% – 97.5%)`.

qci(x) = (quantile(x, 0.025), quantile(x, 0.5), quantile(x, 0.975))

function fmt(x)
    lo, med, hi = qci(x)
    return @sprintf("%.2f (%.2f – %.2f)", med, lo, hi)
end

headline_estimates = DataFrame(
    "Delay" => [
        "Onset → admission",
        "Admission → death",
        "Admission → discharge",
        "Onset → notification",
    ],
    "n" => [
        length(d.onset_to_admit),
        length(d.admit_to_death),
        length(d.admit_to_discharge),
        length(d.onset_to_notif),
    ],
    "Gamma median (95% CrI), days" =>
        [fmt(post.median_oa), fmt(post.median_ad),
         fmt(post.median_ac), fmt(post.median_on)],
    "Gamma mean (95% CrI), days" =>
        [fmt(post.mean_oa),   fmt(post.mean_ad),
         fmt(post.mean_ac),   fmt(post.mean_on)],
    "Rosello mean" => [4.00, 7.59, 8.00, 8.83],
)

# The same numbers are written to disk as a Markdown snippet so the docs
# build can splice them into the home page — keeps the two views in sync.

include(joinpath(pkgdir(BdbvLinelist), "docs", "examples", "_helpers.jl"))    #hide
let snippet_path = joinpath(BdbvLinelist.OUTPUT_DIR, "cache", "headline.md")  #hide
    mkpath(dirname(snippet_path))                                             #hide
    write(snippet_path, headline_snippet(post, d))                            #hide
    nothing                                                                   #hide
end                                                                           #hide

# ## Epidemic curve
#
# Weekly onset counts with HCW subcounts stacked.

plot_epi_curve(ll)

# ## Family comparison
#
# WAIC ranking and convergence diagnostics for each fit.
# Lower WAIC is better; ΔWAIC is relative to the best family.

families = (:lognormal, :gamma, :weibull)
waics = [results[f].waic.waic for f in families]
best  = minimum(waics)

family_comparison = DataFrame(
    family      = collect(families),
    WAIC        = round.(waics, digits = 1),
    ΔWAIC       = round.(waics .- best, digits = 1),
    p_waic      = round.([results[f].waic.p_waic for f in families], digits = 1),
    max_Rhat    = round.([results[f].diag.rhat for f in families], digits = 3),
    min_ESS     = [round(Int, results[f].diag.ess) for f in families],
    n_divergent = [results[f].diag.ndiv for f in families],
)

# Side-by-side posterior-predictive comparison across families. One
# column per family, four rows (one per delay). The LogNormal
# tail-over-fit on the onset → notification panel is what drives the
# WAIC penalty. The black dashed line in each panel marks the
# Rosello *et al.* 2015 Table 5 mean for that delay (capped at 30
# days in their fits) — the onset → notification gap is the most
# visible consequence of that cap.

chains_by_family = Dict(f => results[f].chain for f in families)
plot_family_comparison(chains_by_family, d)

# ## Derived (convolved) marginals
#
# Onset → death = (Onset → admission) ⊛ (Admission → death);
# Onset → discharge = (Onset → admission) ⊛ (Admission → discharge).
# Sampled per posterior draw (500 realisations per draw).

convolved_marginals = DataFrame(
    marginal = [
        "Onset → death (oa ⊛ ad)",
        "Onset → discharge (oa ⊛ ac)",
    ],
    median = [fmt(post.od_median), fmt(post.oc_median)],
    mean   = [fmt(post.od_mean),   fmt(post.oc_mean)],
    P95    = [fmt(post.od_p95),    fmt(post.oc_p95)],
)

# ## Stratified case-fatality

cfr_table = DataFrame(
    stratum = [
        "Non-HCW, Confirmed (baseline)",
        "HCW, Confirmed",
        "Non-HCW, Probable",
        "HCW, Probable",
    ],
    CFR = [
        fmt(post.cfr_baseline),
        fmt(post.cfr_hcw_conf),
        fmt(post.cfr_nonhcw_prob),
        fmt(post.cfr_hcw_prob),
    ],
)

# Logit-scale coefficients and odds ratios.

function fmt_or(β)
    lo, med, hi = qci(β)
    log_or = @sprintf("%+.2f (%+.2f, %+.2f)", med, lo, hi)
    or     = @sprintf("%.2f (%.2f, %.2f)", exp(med), exp(lo), exp(hi))
    return log_or, or
end

logit_coefficients = DataFrame(
    map(((label, β),) -> begin
            log_or, or = fmt_or(β)
            (; coefficient = label, log_OR = log_or, OR = or)
        end,
        [("HCW status", post.β_hcw),
         ("Probable case definition", post.β_def),
         ("Standardised age", post.β_age)]),
)

# ## Posterior predictive check (Gamma)
#
# Four panels — one per atomic delay — overlaying the observed integer
# day histogram with the simulated double-interval-censored posterior
# predictive (median + 95% band). The dashed black line in each
# panel marks the Rosello *et al.* 2015 Table 5 mean (4.00, 7.59,
# 8.00, 8.83 d) for visual comparison.

plot_ppc(chn_gamma, d, :gamma)

# ## Prior sensitivity
#
# Refit the Gamma model under three prior-scale settings (0.5, 1.0
# default, 2.0).
# The default 1.0 prior gives ≈ ×3-fold prior latitude on the central
# tendency of each delay; a 4× span across the sweep shifts posterior
# means by < 5%.

sens = redirect_stdout(devnull) do
    sensitivity()
end
nothing #hide

# Pull the per-draw posterior mean for each atomic delay out of each
# sensitivity chain. The Gamma submodel parametrises each delay as
# `(log_mean, log_shape)`, so `exp(log_mean)` is the population mean.

const DELAY_LM_VARS = (
    onset_to_admit     = @varname(dist_oa.log_mean),
    admit_to_death     = @varname(dist_ad.log_mean),
    admit_to_discharge = @varname(dist_ac.log_mean),
    onset_to_notif     = @varname(dist_on.log_mean),
)

chain_mean(chn, var) = exp.(vec(collect(chn[var])))

prior_sensitivity = DataFrame(
    delay = ["Onset → admission", "Admission → death",
             "Admission → discharge", "Onset → notification"],
    scale_05 = [fmt(chain_mean(sens[0.5].chain, v))
                for v in DELAY_LM_VARS],
    scale_10 = [fmt(chain_mean(sens[1.0].chain, v))
                for v in DELAY_LM_VARS],
    scale_20 = [fmt(chain_mean(sens[2.0].chain, v))
                for v in DELAY_LM_VARS],
)

# ## See also
#
# - [Model](model.md) — likelihood, priors and post-processing details.
# - [Limitations](limitations.md) — known caveats.
# - [API Reference](api.md) — exported functions.
