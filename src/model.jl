## Turing model for the 2012 Isiro BDBV outbreak.
##
## Per-case latent event times (T_onset, T_admit, T_death/T_disch,
## T_notif — sampled where observed) are shared across that case's
## delays, so the four atomic components are fitted jointly and the
## natural-history identity (onset→death = onset→admit ⊕ admit→death)
## holds per case for every posterior draw:
##
##   d_oa  onset → admission     (n = 40)
##   d_ad  admission → death     (n = 22)
##   d_ac  admission → discharge (n = 15)
##   d_on  onset → notification  (n = 38)
##
## Within-day positions use the bounded-primary reparametrisation
## (the secondary event time bounds the primary's upper window edge
## directly) so NUTS doesn't see the wedge-shaped corner that a
## naive `T_admit ≤ T_death` ordering constraint produces at same-day
## cases. Each bounded prior carries a `log(upper − L)` Jacobian to
## restore the implicit independent-uniform-over-day-window prior of
## the equivalent marginalised double-interval-censoring model.
##
## Three parametric families are supported (LogNormal, Gamma, Weibull),
## selected via `family ∈ (:lognormal, :gamma, :weibull)`. All three
## share a common log-mean / log-shape parametrisation so the priors
## are comparable:
##
##   log_mean  ~ Normal(log(plausible_median_d), 1.0)
##   log_shape ~ Normal(0, 1.0)
##
## The distribution constructor maps (mean, shape) onto the canonical
## (μ, σ) for LogNormal, (k, θ = mean/k) for Gamma, and
## (α, θ = mean/Γ(1+1/α)) for Weibull.
##
## The marginal onset → death and onset → discharge population
## distributions are derived in post-processing as Monte-Carlo
## convolutions of the fitted atomic components.
##
## A separate CFR block (Bernoulli with logistic link on HCW, case
## definition, standardised age) uses all 52 cases.

# Numerically safe logistic — pin away from {0, 1} to keep
# Bernoulli's domain check happy when η drifts during NUTS warmup.
_logistic(x) = clamp(inv(1 + exp(-x)), 1e-10, 1.0 - 1e-10)

# Family singletons used for dispatch. The symbol-keyed public API
# (`:lognormal`, `:gamma`, `:weibull`) is converted to one of these
# at a single boundary, `delay_family`; everything internal dispatches
# on the type.
abstract type DelayFamily end
struct LogNormalDelay <: DelayFamily end
struct GammaDelay     <: DelayFamily end
struct WeibullDelay   <: DelayFamily end

delay_family(s::Symbol) =
    s === :lognormal ? LogNormalDelay() :
    s === :gamma     ? GammaDelay()     :
    s === :weibull   ? WeibullDelay()   :
    throw(ArgumentError("unknown family $s"))

family_symbol(::LogNormalDelay) = :lognormal
family_symbol(::GammaDelay)     = :gamma
family_symbol(::WeibullDelay)   = :weibull

# Construct a delay distribution from already-sampled location and
# scale parameters. For LogNormal the pair is (log_median, log_sd);
# for Gamma and Weibull it is (log_mean, log_shape).
build_delay_dist(::LogNormalDelay, log_median, log_sd) =
    LogNormal(log_median, exp(log_sd))

function build_delay_dist(::GammaDelay, log_mean, log_shape)
    shape = exp(log_shape)
    return Gamma(shape, exp(log_mean) / shape)
end

function build_delay_dist(::WeibullDelay, log_mean, log_shape)
    shape = exp(log_shape)
    return Weibull(shape, exp(log_mean) / SpecialFunctions.gamma(1 + 1 / shape))
end

# Prior submodels — one method per family. Each samples the
# family's natural location/scale pair and returns the constructed
# distribution. Used with `to_submodel()` so the sample names are
# prefixed by the LHS variable.
@model function delay_prior(::LogNormalDelay, loc_loc, loc_scale, sd_scale)
    log_median ~ Normal(loc_loc, loc_scale)
    log_sd     ~ truncated(Normal(0.0, sd_scale); lower = 0.0)
    return build_delay_dist(LogNormalDelay(), log_median, log_sd)
end

@model function delay_prior(::GammaDelay, mean_loc, mean_scale, shape_scale)
    log_mean  ~ Normal(mean_loc, mean_scale)
    log_shape ~ Normal(0.0, shape_scale)
    return build_delay_dist(GammaDelay(), log_mean, log_shape)
end

# Weibull's log-shape is truncated to keep `α ∈ (≈0.37, ≈2.7)` so
# `Γ(1 + 1/α)` stays well-defined and the scale stays positive under NUTS.
@model function delay_prior(::WeibullDelay, mean_loc, mean_scale, shape_scale)
    log_mean  ~ Normal(mean_loc, mean_scale)
    log_shape ~ truncated(Normal(0.0, shape_scale); lower = -1.0, upper = 1.0)
    return build_delay_dist(WeibullDelay(), log_mean, log_shape)
end

"""
$(TYPEDSIGNATURES)

Turing model for the BDBV Isiro 2012 line list, parametrised by
the choice of delay-distribution family (`:lognormal`, `:gamma`,
or `:weibull`). Estimates four delay components and a stratified
CFR logistic regression.

Per-case latent event times (`T_onset`, `T_admit`, `T_death`,
`T_disch`, `T_notif`, sampled where observed) are shared across
that case's atomic delays, so the natural-history identity
(`D_oa + D_ad = D_od` per case) is enforced *in* the model rather
than recovered in post-processing. Within-day positions are
sampled with Sam Abbott's bounded-primary reparametrisation —

    T_death/T_disch ~ Uniform(day, day + 1)                       # leaf
    T_admit         ~ Uniform(day, min(day + 1, T_death, T_disch))
    T_notif         ~ Uniform(day, day + 1)
    T_onset         ~ Uniform(day, min(day + 1, T_admit, T_notif))

— which absorbs the ordering constraint into the support and
avoids the wedge-shaped boundary geometry NUTS handles poorly at
same-day cases. A `log(upper − L)` Jacobian on each bounded prior
restores the implicit independent-uniform-over-day-window prior of
the equivalent marginalised model.

Reference: Park *et al.* 2024 (medRxiv
[2024.01.12.24301247](https://doi.org/10.1101/2024.01.12.24301247))
§2.3.3 for the latent-variable formulation; the bounded-primary
trick is from Sam Abbott.
"""
@model function bdbv_model(d; family::Symbol = :gamma,
        prior_scale::Float64 = 1.0)
    fam = delay_family(family)

    dist_oa ~ to_submodel(delay_prior(fam, log(3.0),  prior_scale, prior_scale))
    dist_ad ~ to_submodel(delay_prior(fam, log(6.0),  prior_scale, prior_scale))
    dist_ac ~ to_submodel(delay_prior(fam, log(13.0), prior_scale, prior_scale))
    dist_on ~ to_submodel(delay_prior(fam, log(7.0),  prior_scale, prior_scale))

    # Vectorised latent sampling — each event type is one `arraydist`
    # over the cases that observe it, so the underlying storage is a
    # concrete-eltype Vector (fast under ForwardDiff) rather than the
    # abstract `Vector{Real}` a per-case loop would need to hold both
    # `Float64` and `Dual`.
    cases = d.case_events
    death_idx = findall(c -> !ismissing(c.death), cases)
    disch_idx = findall(c -> !ismissing(c.disch), cases)
    notif_idx = findall(c -> !ismissing(c.notif), cases)
    admit_idx = findall(c -> !ismissing(c.admit), cases)
    onset_idx = findall(c -> !ismissing(c.onset), cases)

    # Leaf events (no upper bound from chain ordering).
    T_death ~ Turing.arraydist([Uniform(cases[i].death, cases[i].death + 1.0)
                                for i in death_idx])
    T_disch ~ Turing.arraydist([Uniform(cases[i].disch, cases[i].disch + 1.0)
                                for i in disch_idx])
    T_notif ~ Turing.arraydist([Uniform(cases[i].notif, cases[i].notif + 1.0)
                                for i in notif_idx])

    # Per-case lookup of leaf samples (case index → sampled time).
    T_death_at = Dict(death_idx .=> T_death)
    T_disch_at = Dict(disch_idx .=> T_disch)
    T_notif_at = Dict(notif_idx .=> T_notif)

    # T_admit's upper bound is shrunk by the secondary event time when
    # the case also has a recorded death and/or discharge (Sam Abbott's
    # bounded-primary trick — see docstring).
    admit_uppers = [let c = cases[i]
        u = c.admit + 1.0
        haskey(T_death_at, i) && (u = min(u, T_death_at[i]))
        haskey(T_disch_at, i) && (u = min(u, T_disch_at[i]))
        u
    end for i in admit_idx]
    T_admit ~ Turing.arraydist([Uniform(cases[admit_idx[k]].admit, admit_uppers[k])
                                for k in eachindex(admit_idx)])
    # Jacobian to match the marginalised model's implicit
    # independent-uniform-over-day-window prior. Vanishes
    # (`log(1) = 0`) on the multi-day cases where the ordering
    # constraint doesn't bind.
    Turing.@addlogprob!(sum(log(admit_uppers[k] - cases[admit_idx[k]].admit)
                            for k in eachindex(admit_idx); init = 0.0))

    T_admit_at = Dict(admit_idx .=> T_admit)

    onset_uppers = [let c = cases[i]
        u = c.onset + 1.0
        haskey(T_admit_at, i) && (u = min(u, T_admit_at[i]))
        haskey(T_notif_at, i) && (u = min(u, T_notif_at[i]))
        u
    end for i in onset_idx]
    T_onset ~ Turing.arraydist([Uniform(cases[onset_idx[k]].onset, onset_uppers[k])
                                for k in eachindex(onset_idx)])
    Turing.@addlogprob!(sum(log(onset_uppers[k] - cases[onset_idx[k]].onset)
                            for k in eachindex(onset_idx); init = 0.0))

    T_onset_at = Dict(onset_idx .=> T_onset)

    # Delay likelihoods: one term per ordered pair both observed in
    # the case, using the shared per-case latents.
    Turing.@addlogprob!(sum(
        logpdf(dist_oa, T_admit[k] - T_onset_at[admit_idx[k]])
        for k in eachindex(admit_idx) if haskey(T_onset_at, admit_idx[k]);
        init = 0.0))
    Turing.@addlogprob!(sum(
        logpdf(dist_ad, T_death[k] - T_admit_at[death_idx[k]])
        for k in eachindex(death_idx) if haskey(T_admit_at, death_idx[k]);
        init = 0.0))
    Turing.@addlogprob!(sum(
        logpdf(dist_ac, T_disch[k] - T_admit_at[disch_idx[k]])
        for k in eachindex(disch_idx) if haskey(T_admit_at, disch_idx[k]);
        init = 0.0))
    Turing.@addlogprob!(sum(
        logpdf(dist_on, T_notif[k] - T_onset_at[notif_idx[k]])
        for k in eachindex(notif_idx) if haskey(T_onset_at, notif_idx[k]);
        init = 0.0))

    # CFR logistic regression: intercept, HCW indicator, case-definition
    # indicator (probable vs confirmed), and standardised age.
    β_0   ~ Normal(0.0, 2.0)
    β_hcw ~ Normal(0.0, 1.0)
    β_def ~ Normal(0.0, 1.0)
    β_age ~ Normal(0.0, 1.0)

    for i in eachindex(d.outcome)
        η = β_0 +
            β_hcw * (d.hcw[i]      ? 1.0 : 0.0) +
            β_def * (d.probable[i] ? 1.0 : 0.0) +
            β_age * d.age_z[i]
        d.outcome[i] ~ Bernoulli(_logistic(η))
    end
end

"""
$(TYPEDSIGNATURES)

Small standalone model: fits a single doubly-censored delay
distribution for the community-died onset → death pathway, plus a
Beta-Binomial posterior for `p_admit` (the population share of
fatal cases that pass through the admit-and-died pathway). Used by
[`fit_death_mixture`](@ref) together with the main `bdbv_model` to
build the death-pathway mixture marginal.
"""
@model function community_death_model(delays; family::Symbol = :gamma,
        n_admit_died::Int = 0, n_comm_died::Int = 0)
    fam = delay_family(family)
    dist_cd ~ to_submodel(delay_prior(fam, log(8.0), 1.0, 1.0))
    if !isempty(delays)
        delays ~ Turing.filldist(
            double_interval_censored(dist_cd; interval = 1.0),
            length(delays),
        )
    end
    # p_admit ~ Beta(1+n_admit, 1+n_comm); independent of the delay fit.
    p_admit ~ Beta(1 + n_admit_died, 1 + n_comm_died)
end

"""
$(TYPEDSIGNATURES)

Stratified-by-HCW model. Adds an `β_*_hcw` log-mean shift to each
of the four atomic delay components, so HCW vs non-HCW cases share
the shape parameter but can differ on the central tendency.

Reports per-delay HCW odds-ratio-like shifts (`exp(β_*_hcw)`) =
multiplicative effect on the delay mean for HCWs vs non-HCWs.
"""
@model function bdbv_model_stratified(d; family::Symbol = :gamma)
    family === :gamma ||
        throw(ArgumentError("stratified model is only supported for the :gamma family. Use bdbv_model for other families."))
    fam = delay_family(family)

    # Shared shape and baseline (non-HCW) log-mean per delay component.
    # Suffixes match the unstratified model: oa, ad, ac, on
    # (onset→admission, admission→death, admission→discharge, onset→notification).
    log_shape_oa ~ Normal(0.0, 1.0)
    log_shape_ad ~ Normal(0.0, 1.0)
    log_shape_ac ~ Normal(0.0, 1.0)
    log_shape_on ~ Normal(0.0, 1.0)

    log_mean_oa ~ Normal(log(3.0),  1.0)
    log_mean_ad ~ Normal(log(6.0),  1.0)
    log_mean_ac ~ Normal(log(13.0), 1.0)
    log_mean_on ~ Normal(log(7.0),  1.0)

    # HCW shifts on log-mean. Mildly informative prior — exp(β) ∈ (~0.4, ~2.7).
    β_oa_hcw ~ Normal(0.0, 0.5)
    β_ad_hcw ~ Normal(0.0, 0.5)
    β_ac_hcw ~ Normal(0.0, 0.5)
    β_on_hcw ~ Normal(0.0, 0.5)

    # Per-stratum likelihoods per delay. Field-name suffix convention:
    # `_h` = HCW subset, `_n` = non-HCW subset (e.g. `d.ac_n` is the
    # non-HCW admission→discharge delays). These subsets are pre-split
    # fields of the data tuple so Turing treats them as observations
    # rather than parameters.
    if !isempty(d.oa_h)
        d.oa_h ~ Turing.filldist(double_interval_censored(
            build_delay_dist(fam, log_mean_oa + β_oa_hcw, log_shape_oa);
            interval = 1.0), length(d.oa_h))
    end
    if !isempty(d.oa_n)
        d.oa_n ~ Turing.filldist(double_interval_censored(
            build_delay_dist(fam, log_mean_oa, log_shape_oa);
            interval = 1.0), length(d.oa_n))
    end

    if !isempty(d.ad_h)
        d.ad_h ~ Turing.filldist(double_interval_censored(
            build_delay_dist(fam, log_mean_ad + β_ad_hcw, log_shape_ad);
            interval = 1.0), length(d.ad_h))
    end
    if !isempty(d.ad_n)
        d.ad_n ~ Turing.filldist(double_interval_censored(
            build_delay_dist(fam, log_mean_ad, log_shape_ad);
            interval = 1.0), length(d.ad_n))
    end

    if !isempty(d.ac_h)
        d.ac_h ~ Turing.filldist(double_interval_censored(
            build_delay_dist(fam, log_mean_ac + β_ac_hcw, log_shape_ac);
            interval = 1.0), length(d.ac_h))
    end
    if !isempty(d.ac_n)
        d.ac_n ~ Turing.filldist(double_interval_censored(
            build_delay_dist(fam, log_mean_ac, log_shape_ac);
            interval = 1.0), length(d.ac_n))
    end

    if !isempty(d.on_h)
        d.on_h ~ Turing.filldist(double_interval_censored(
            build_delay_dist(fam, log_mean_on + β_on_hcw, log_shape_on);
            interval = 1.0), length(d.on_h))
    end
    if !isempty(d.on_n)
        d.on_n ~ Turing.filldist(double_interval_censored(
            build_delay_dist(fam, log_mean_on, log_shape_on);
            interval = 1.0), length(d.on_n))
    end

    # CFR block — same as the unstratified model.
    β_0   ~ Normal(0.0, 2.0)
    β_hcw ~ Normal(0.0, 1.0)
    β_def ~ Normal(0.0, 1.0)
    β_age ~ Normal(0.0, 1.0)

    for i in eachindex(d.outcome)
        η = β_0 +
            β_hcw * (d.hcw[i]      ? 1.0 : 0.0) +
            β_def * (d.probable[i] ? 1.0 : 0.0) +
            β_age * d.age_z[i]
        d.outcome[i] ~ Bernoulli(_logistic(η))
    end
end
