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
# Per-draw posterior CSV and figure bundle from this build is in the
# rolling [`main-latest` release](https://github.com/epiforecasts/bdbv-linelist-analysis/releases/tag/main-latest).

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

# ## Outbreak context
#
# Per the Charniga *et al.* 2024 reporting checklist (item 12: provide
# the contextual information needed to interpret the delays). Sample
# size, observation window, demographics, case-definition mix and
# care-setting are summarised here; epidemic curve and control measures
# are below ([Epidemic curve](#Epidemic-curve)) and on the
# [Limitations](limitations.md) page (outbreak setting and ETC support).
let onsets = collect(skipmissing(ll.Date_of_onset_symp)),
    ages   = collect(skipmissing(ll.Age)),
    sexes  = collect(skipmissing(ll.Sex))
    DataFrame(
        "Quantity" => [
            "Total cases", "Cases with onset date",
            "Onset window",
            "Median age (IQR), years",
            "Sex (Female / Male)",
            "HCW (yes / no)",
            "Case definition (Confirmed / Probable)",
            "Outcome (Dead / Alive / Unknown)",
        ],
        "Value" => [
            string(nrow(ll)),
            string(length(onsets)),
            @sprintf("%s to %s", minimum(onsets), maximum(onsets)),
            @sprintf("%.0f (%.0f – %.0f)",
                     quantile(ages, 0.5),
                     quantile(ages, 0.25),
                     quantile(ages, 0.75)),
            @sprintf("%d / %d",
                     count(==("Female"), sexes), count(==("Male"), sexes)),
            @sprintf("%d / %d", sum(ll.is_hcw), nrow(ll) - sum(ll.is_hcw)),
            @sprintf("%d / %d",
                     count(==("Confirmed"), skipmissing(ll.Case_definition)),
                     count(==("Probable"),  skipmissing(ll.Case_definition))),
            @sprintf("%d / %d / %d",
                     count(==("Dead"),  skipmissing(ll.Outcome)),
                     count(==("Alive"), skipmissing(ll.Outcome)),
                     count(o -> !(o in ("Dead", "Alive")),
                           skipmissing(ll.Outcome))),
        ],
    )
end

# ## Headline Gamma estimates
#
# Fit all three families (Gamma wins on WAIC; comparison table and
# plot are in the [Family comparison](#Family-comparison) section
# below) and tabulate the four atomic delay components alongside
# Rosello *et al.* 2015 Table 5 empirical means. Rosello capped the
# raw delays at 30 days before computing those summary statistics —
# this is binding on onset → notification (where 9 of 38 cases
# exceed 30 d) and accounts for most of the difference there.

## Suppress the package's stdout printing — we build clean DataFrame
## tables from the returned posterior vectors below. The fit itself
## emits a fair amount of Turing diagnostics (initial step size,
## sampling progress, the occasional divergent-transition warning)
## that aren't part of the narrative, so the call and its captured
## output are tucked into the dropdown below.

#md # ```@raw html
#md # <details><summary>Fit: <code>compare_families()</code> — Turing diagnostics</summary>
#md # ```

results = redirect_stdout(devnull) do
    compare_families()
end
nothing #hide

#md # ```@raw html
#md # </details>
#md # ```

#md # ```@raw html
#md # <details><summary>Post-processing: <code>summarise()</code> on the Gamma fit</summary>
#md # ```

chn_gamma = results[:gamma].chain
post = redirect_stdout(devnull) do
    summarise(chn_gamma, :gamma; d = d)
end
nothing #hide

#md # ```@raw html
#md # </details>
#md # ```

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

# ## Posterior predictive check (Gamma)
#
# Four panels — one per atomic delay — overlaying the observed integer
# day histogram with the simulated double-interval-censored posterior
# predictive (median + 95% band). The dashed black line in each
# panel marks the Rosello *et al.* 2015 Table 5 mean (4.00, 7.59,
# 8.00, 8.83 d) for visual comparison.

plot_ppc(chn_gamma, d, :gamma)

# ## Prior versus posterior

# Pairwise distributions of the four atomic Gamma `log_mean` parameters
# under the prior used in fitting and the resulting posterior. Diagonal
# panels show the marginals; off-diagonal panels show the pairwise
# contours. The shrinkage from prior (grey) to posterior (red) is the
# visual analogue of the prior-sensitivity table immediately below.

using PairPlots

prior_vs_posterior = let
    posterior_log_means = (;
        oa = vec(collect(chn_gamma[@varname(dist_oa.log_mean)])),
        ad = vec(collect(chn_gamma[@varname(dist_ad.log_mean)])),
        ac = vec(collect(chn_gamma[@varname(dist_ac.log_mean)])),
        on = vec(collect(chn_gamma[@varname(dist_on.log_mean)])),
    )
    ## Sample the prior directly: each delay's log_mean is
    ## Normal(log(plausible_median_d), 1.0) at the default prior scale —
    ## see `bdbv_model` in `src/model.jl`.
    S = length(posterior_log_means.oa)
    rng = Random.MersenneTwister(20260519)
    prior_log_means = (;
        oa = log(3.0)  .+ randn(rng, S),
        ad = log(6.0)  .+ randn(rng, S),
        ac = log(13.0) .+ randn(rng, S),
        on = log(7.0)  .+ randn(rng, S),
    )
    pairplot(
        PairPlots.Series(prior_log_means;     label = "prior",     color = :grey),
        PairPlots.Series(posterior_log_means; label = "posterior", color = :firebrick),
    )
end

# ## Prior sensitivity
#
# Refit the Gamma model under three prior-scale settings (0.5, 1.0
# default, 2.0).
# The default 1.0 prior gives ≈ ×3-fold prior latitude on the central
# tendency of each delay; a 4× span across the sweep shifts posterior
# means by < 5%.

#md # ```@raw html
#md # <details><summary>Fit: <code>sensitivity()</code> — three prior-scale refits</summary>
#md # ```

sens = redirect_stdout(devnull) do
    sensitivity()
end
nothing #hide

#md # ```@raw html
#md # </details>
#md # ```

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

# ## Epidemic curve
#
# Weekly onset counts with HCW subcounts stacked.

plot_epi_curve(ll)

# ## Early-phase growth rate
#
# An exploratory exponential-growth fit to the rising phase of the
# weekly onset curve — week 1 (2012-05-28) through the peak week
# (2012-09-10), inclusive. The model is a Poisson regression
# `log(λ_t) = α + r·t` with weakly-informative `Normal(0, 5)` and
# `Normal(0, 1)` priors on `α` and `r_week`. Intended as a prior
# source for downstream re-applications (e.g. the outbreak-size work
# in `epiforecasts/BVDOutbreakSize`) rather than as a primary
# headline of this analysis. The CrI on `r` covers zero — Isiro was
# a slow, noisy rise — so use the posterior as a weakly-informative
# prior, not a tight constraint.

growth = redirect_stdout(devnull) do
    fit_growth_rate(ll)
end
nothing #hide

# Doubling time `log(2)/r` is reported as a median only; its
# distribution is heavy-tailed because the posterior on `r` includes
# values close to zero. The credibility statement worth quoting is
# the CrI on `r` itself.

growth_estimates = DataFrame(
    "Quantity" => [
        "Growth rate r (per week)",
        "Growth rate r (per day)",
        "Doubling time (days, median)",
        "P(r > 0)",
    ],
    "Posterior summary" => [
        fmt(growth.r_week),
        fmt(growth.r_day),
        @sprintf("%.1f", quantile(growth.doubling_time, 0.5)),
        @sprintf("%.2f", mean(growth.r_day .> 0)),
    ],
)

# **Recommended downstream prior on `r` (per day):**
# `Normal(mean(r_day), sd(r_day))` from the posterior above. A
# single-line summary intended for `epiforecasts/BVDOutbreakSize`
# and any other downstream model that needs an Isiro-anchored growth
# prior. See [Limitations](limitations.md#downstream-priors) for
# caveats.

prior_summary = @sprintf("Normal(%.4f, %.4f)",
                         mean(growth.r_day),
                         std(growth.r_day))

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
# Rosello *et al.* 2015 Table 5 empirical mean for that delay
# (computed after capping the raw delays at 30 days) — the
# onset → notification gap is the most visible consequence of that cap.

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
    sd     = [fmt(post.od_sd),     fmt(post.oc_sd)],
    P95    = [fmt(post.od_p95),    fmt(post.oc_p95)],
)

# ## Length of stay in hospital
#
# Time from admission to leaving hospital (admission → departure). The
# overall length of stay is a mixture of the fatal pathway
# (admission → death) and the survivor pathway (admission → discharge),
# weighted per posterior draw by the in-hospital fatality among admitted
# cases — `Beta(1 + n_died, 1 + n_discharged)` with 22 deaths and 15
# discharges. The fatal and survivor rows are the corresponding atomic
# components; the overall row is the bed-occupancy-relevant marginal
# across both outcomes.

length_of_stay = DataFrame(
    pathway = [
        "Fatal (admission → death)",
        "Survivor (admission → discharge)",
        "Overall (mixture)",
    ],
    median = [fmt(post.median_ad), fmt(post.median_ac), fmt(post.los_median)],
    mean   = [fmt(post.mean_ad),   fmt(post.mean_ac),   fmt(post.los_mean)],
    sd     = ["—",                 "—",                 fmt(post.los_sd)],
    P95    = ["—",                 "—",                 fmt(post.los_p95)],
)

# ## Gamma shape, scale and SD per atomic delay
#
# Underlying-distribution parameters for downstream consumers that need
# to reconstruct each atomic Gamma delay rather than just its central
# tendency.

gamma_parameters = DataFrame(
    "delay" => [
        "Onset → admission",
        "Admission → death",
        "Admission → discharge",
        "Onset → notification",
    ],
    "shape (95% CrI)" =>
        [fmt(post.shape_oa), fmt(post.shape_ad),
         fmt(post.shape_ac), fmt(post.shape_on)],
    "scale (95% CrI)" =>
        [fmt(post.scale_oa), fmt(post.scale_ad),
         fmt(post.scale_ac), fmt(post.scale_on)],
    "sd (95% CrI)" =>
        [fmt(post.sd_oa), fmt(post.sd_ad),
         fmt(post.sd_ac), fmt(post.sd_on)],
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

# ## See also
#
# - [Model](model.md) — likelihood, priors and post-processing details.
# - [Limitations](limitations.md) — known caveats.
# - [API Reference](api.md) — exported functions.
