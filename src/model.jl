## Turing model for the 2012 Isiro BDBV outbreak. The four atomic
## delay components and the CFR block share no latent variables —
## the natural-history identity (onset→death = onset→admit ⊕
## admit→death) is recovered in post-processing via convolution
## (see LIMITATIONS.md for why the latent-time joint variant was
## abandoned).
##
## Four atomic delay components are fitted with double interval
## censoring at the day level (`CensoredDistributions.double_interval_censored`,
## default `primary_event = Uniform(0, 1)`, `interval = 1.0`):
##
##   d_oa  onset → admission     (n = 40)
##   d_ad  admission → death     (n = 22)
##   d_ac  admission → discharge (n = 15)
##   d_on  onset → notification  (n = 38)
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
## The marginal onset → death and onset → discharge distributions are
## derived in post-processing as convolutions of d_oa with d_ad and
## d_ac (only for LogNormal, since the convolution is closed-form-
## adjacent there; Gamma/Weibull convolutions are sampled).
##
## A separate CFR block (Bernoulli with logistic link on HCW, case
## definition, standardised age) uses all 52 cases.
##
## Per-case censoring noise across the four delays of the same case is
## treated as independent — see `LIMITATIONS.md`.

# Numerically safe logistic — pin away from {0, 1} to keep
# Bernoulli's domain check happy when η drifts during NUTS warmup.
_logistic(x) = clamp(inv(1 + exp(-x)), 1e-10, 1.0 - 1e-10)

# Compress a delay vector into the unique values and their integer
# multiplicities. Used to weight the censored-likelihood contributions
# in the marginalised models so each unique (lower, upper) pair only
# needs one `logpdf` call per evaluation. At day-level censoring and
# small N this collapses 38 observations down to ~10-15 unique values.
function _unique_counts(v::AbstractVector)
    counts = Dict{eltype(v), Int}()
    for x in v
        counts[x] = get(counts, x, 0) + 1
    end
    uniques = collect(keys(counts))
    return uniques, [counts[u] for u in uniques]
end

# Weighted sum of `logpdf(dist, u)` over unique observations `uniques`
# with integer multiplicities `counts`. Equivalent to summing
# `logpdf(dist, x)` over the original (de-duplicated) observation
# vector, but with one `logpdf` call per unique value.
@inline function _weighted_loglik(dist, uniques::AbstractVector, counts::AbstractVector{Int})
    s = zero(logpdf(dist, first(uniques)))
    @inbounds for i in eachindex(uniques)
        s += counts[i] * logpdf(dist, uniques[i])
    end
    return s
end

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
or `:weibull`). Estimates four doubly-censored delay components via
the submodel pattern of CensoredDistributions.jl, and a stratified
CFR logistic regression. The delay components and CFR block share
no latent variables; the marginal onset→death and onset→discharge
distributions are derived in post-processing as convolutions of
the atomic components.
"""
@model function bdbv_model(d; family::Symbol = :gamma, prior_scale::Float64 = 1.0)
    fam = delay_family(family)

    # Delay-component suffixes used throughout: oa = onset→admission,
    # ad = admission→death, ac = admission→discharge, on = onset→notification.
    # `dist_*` is the underlying continuous distribution; `dic_*` is its
    # doubly-interval-censored counterpart used in the likelihood.
    dist_oa ~ to_submodel(delay_prior(fam, log(3.0),  prior_scale, prior_scale))
    dist_ad ~ to_submodel(delay_prior(fam, log(6.0),  prior_scale, prior_scale))
    dist_ac ~ to_submodel(delay_prior(fam, log(13.0), prior_scale, prior_scale))
    dist_on ~ to_submodel(delay_prior(fam, log(7.0),  prior_scale, prior_scale))

    dic_oa = double_interval_censored(dist_oa; interval = 1.0)
    dic_ad = double_interval_censored(dist_ad; interval = 1.0)
    dic_ac = double_interval_censored(dist_ac; interval = 1.0)
    dic_on = double_interval_censored(dist_on; interval = 1.0)

    # Day-level censoring means many cases share the same delay value;
    # compress each observation vector into uniques + counts and weight
    # the censored-likelihood contribution by the multiplicity. One
    # `logpdf` call per unique value rather than per case (issue #4).
    u_oa, c_oa = _unique_counts(d.onset_to_admit)
    u_ad, c_ad = _unique_counts(d.admit_to_death)
    u_ac, c_ac = _unique_counts(d.admit_to_discharge)
    u_on, c_on = _unique_counts(d.onset_to_notif)

    Turing.@addlogprob! _weighted_loglik(dic_oa, u_oa, c_oa)
    Turing.@addlogprob! _weighted_loglik(dic_ad, u_ad, c_ad)
    Turing.@addlogprob! _weighted_loglik(dic_ac, u_ac, c_ac)
    Turing.@addlogprob! _weighted_loglik(dic_on, u_on, c_on)

    # CFR logistic regression coefficients: intercept, HCW indicator,
    # case-definition indicator (probable vs confirmed), and standardised age.
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
        # Weighted-by-multiplicity likelihood — see `bdbv_model`.
        u_cd, c_cd = _unique_counts(delays)
        Turing.@addlogprob! _weighted_loglik(
            double_interval_censored(dist_cd; interval = 1.0), u_cd, c_cd,
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

    # Per-stratum weighted-by-multiplicity likelihood per delay. Field-name
    # suffix convention: `_h` = HCW subset, `_n` = non-HCW subset (e.g.
    # `d.ac_n` is the non-HCW admission→discharge delays). These subsets
    # are pre-split fields of the data tuple. Compression to uniques +
    # counts uses the same pattern as `bdbv_model` (see issue #4).
    if !isempty(d.oa_h)
        u, c = _unique_counts(d.oa_h)
        Turing.@addlogprob! _weighted_loglik(double_interval_censored(
            build_delay_dist(fam, log_mean_oa + β_oa_hcw, log_shape_oa);
            interval = 1.0), u, c)
    end
    if !isempty(d.oa_n)
        u, c = _unique_counts(d.oa_n)
        Turing.@addlogprob! _weighted_loglik(double_interval_censored(
            build_delay_dist(fam, log_mean_oa, log_shape_oa);
            interval = 1.0), u, c)
    end

    if !isempty(d.ad_h)
        u, c = _unique_counts(d.ad_h)
        Turing.@addlogprob! _weighted_loglik(double_interval_censored(
            build_delay_dist(fam, log_mean_ad + β_ad_hcw, log_shape_ad);
            interval = 1.0), u, c)
    end
    if !isempty(d.ad_n)
        u, c = _unique_counts(d.ad_n)
        Turing.@addlogprob! _weighted_loglik(double_interval_censored(
            build_delay_dist(fam, log_mean_ad, log_shape_ad);
            interval = 1.0), u, c)
    end

    if !isempty(d.ac_h)
        u, c = _unique_counts(d.ac_h)
        Turing.@addlogprob! _weighted_loglik(double_interval_censored(
            build_delay_dist(fam, log_mean_ac + β_ac_hcw, log_shape_ac);
            interval = 1.0), u, c)
    end
    if !isempty(d.ac_n)
        u, c = _unique_counts(d.ac_n)
        Turing.@addlogprob! _weighted_loglik(double_interval_censored(
            build_delay_dist(fam, log_mean_ac, log_shape_ac);
            interval = 1.0), u, c)
    end

    if !isempty(d.on_h)
        u, c = _unique_counts(d.on_h)
        Turing.@addlogprob! _weighted_loglik(double_interval_censored(
            build_delay_dist(fam, log_mean_on + β_on_hcw, log_shape_on);
            interval = 1.0), u, c)
    end
    if !isempty(d.on_n)
        u, c = _unique_counts(d.on_n)
        Turing.@addlogprob! _weighted_loglik(double_interval_censored(
            build_delay_dist(fam, log_mean_on, log_shape_on);
            interval = 1.0), u, c)
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
