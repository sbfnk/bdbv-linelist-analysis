## Joint Turing model for the 2012 Isiro BDBV outbreak.
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

"""
$(TYPEDSIGNATURES)

Submodel returning a LogNormal delay distribution, parametrised by
log-mean and log-SD (the canonical LogNormal `(μ, σ)`). Used with
`to_submodel()` so the (log_mean, log_sd) pair is prefixed by the
LHS name.
"""
@model function delay_lognormal(log_median_loc, log_median_scale, log_sd_scale)
    log_median ~ Normal(log_median_loc, log_median_scale)
    log_sd     ~ truncated(Normal(0.0, log_sd_scale); lower = 0.0)
    return LogNormal(log_median, log_sd)
end

"""
$(TYPEDSIGNATURES)

Submodel returning a Gamma delay distribution parametrised by log-mean
and log-shape (shape `k`, scale `mean / k`).
"""
@model function delay_gamma(log_mean_loc, log_mean_scale, log_shape_scale)
    log_mean  ~ Normal(log_mean_loc, log_mean_scale)
    log_shape ~ Normal(0.0, log_shape_scale)
    shape = exp(log_shape)
    scale = exp(log_mean) / shape
    return Gamma(shape, scale)
end

"""
$(TYPEDSIGNATURES)

Submodel returning a Weibull delay distribution parametrised by
log-mean and log-shape. Weibull `α` is the shape; scale `θ` is set
so that `θ Γ(1 + 1/α) = mean`. The log-shape prior is truncated
to keep `α ∈ (≈0.37, ≈2.7)` so `Γ(1 + 1/α)` doesn't blow up and
the scale stays positive under NUTS.
"""
@model function delay_weibull(log_mean_loc, log_mean_scale, log_shape_scale)
    log_mean  ~ Normal(log_mean_loc, log_mean_scale)
    log_shape ~ truncated(Normal(0.0, log_shape_scale); lower = -1.0, upper = 1.0)
    shape = exp(log_shape)
    scale = exp(log_mean) / SpecialFunctions.gamma(1 + 1 / shape)
    return Weibull(shape, scale)
end

# Map family symbol → submodel constructor. Used by `bdbv_model` to
# select the family without writing three near-identical model bodies.
const _DELAY_SUBMODEL = Dict(
    :lognormal => delay_lognormal,
    :gamma     => delay_gamma,
    :weibull   => delay_weibull,
)

"""
$(TYPEDSIGNATURES)

Joint Turing model for the BDBV Isiro 2012 line list, parametrised
by the choice of delay-distribution family (`:lognormal`, `:gamma`,
or `:weibull`). Estimates four doubly-censored delay components via
the submodel pattern of CensoredDistributions.jl, and a stratified
CFR logistic regression. The marginal onset→death and onset→discharge
distributions are derived in post-processing.
"""
@model function bdbv_model(d; family::Symbol = :gamma, prior_scale::Float64 = 1.0)
    submodel = _DELAY_SUBMODEL[family]

    dist_oa ~ to_submodel(submodel(log(3.0),  prior_scale, prior_scale))
    dist_ad ~ to_submodel(submodel(log(6.0),  prior_scale, prior_scale))
    dist_ac ~ to_submodel(submodel(log(13.0), prior_scale, prior_scale))
    dist_on ~ to_submodel(submodel(log(7.0),  prior_scale, prior_scale))

    dic_oa = double_interval_censored(dist_oa; interval = 1.0)
    dic_ad = double_interval_censored(dist_ad; interval = 1.0)
    dic_ac = double_interval_censored(dist_ac; interval = 1.0)
    dic_on = double_interval_censored(dist_on; interval = 1.0)

    d.onset_to_admit     ~ Turing.filldist(dic_oa, length(d.onset_to_admit))
    d.admit_to_death     ~ Turing.filldist(dic_ad, length(d.admit_to_death))
    d.admit_to_discharge ~ Turing.filldist(dic_ac, length(d.admit_to_discharge))
    d.onset_to_notif     ~ Turing.filldist(dic_on, length(d.onset_to_notif))

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
    submodel = _DELAY_SUBMODEL[family]
    dist_cd ~ to_submodel(submodel(log(8.0), 1.0, 1.0))
    if !isempty(delays)
        delays ~ Turing.filldist(
            double_interval_censored(dist_cd; interval = 1.0),
            length(delays),
        )
    end
    # p_admit ~ Beta(1+n_admit, 1+n_comm); independent of the delay fit.
    p_admit ~ Beta(1 + n_admit_died, 1 + n_comm_died)
end

# Inline construction of a delay distribution from (log_mean, log_shape)
# given a family. Matches the parametrisation of the submodels above.
function _build_delay_dist(family::Symbol, log_mean, log_shape)
    if family == :lognormal
        # For LogNormal, `log_mean` here is interpreted as log_median.
        return LogNormal(log_mean, exp(log_shape))
    elseif family == :gamma
        shape = exp(log_shape)
        return Gamma(shape, exp(log_mean) / shape)
    elseif family == :weibull
        shape = exp(log_shape)
        return Weibull(shape, exp(log_mean) / SpecialFunctions.gamma(1 + 1 / shape))
    else
        error("unknown family $family")
    end
end

"""
$(TYPEDSIGNATURES)

Stratified-by-HCW joint model. Adds an `β_*_hcw` log-mean shift to
each of the four atomic delay components, so HCW vs non-HCW cases
share the shape parameter but can differ on the central tendency.

Reports per-delay HCW odds-ratio-like shifts (`exp(β_*_hcw)`) =
multiplicative effect on the delay mean for HCWs vs non-HCWs.
"""
@model function bdbv_model_stratified(d; family::Symbol = :gamma)
    family === :gamma ||
        throw(ArgumentError("stratified model is only supported for the :gamma family. Use bdbv_model for other families."))

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

    # Per-stratum likelihoods per delay. The HCW/non-HCW subsets are
    # pre-split fields of the data tuple so Turing treats them as
    # observations rather than parameters.
    if !isempty(d.oa_h)
        d.oa_h ~ Turing.filldist(double_interval_censored(
            _build_delay_dist(family, log_mean_oa + β_oa_hcw, log_shape_oa);
            interval = 1.0), length(d.oa_h))
    end
    if !isempty(d.oa_n)
        d.oa_n ~ Turing.filldist(double_interval_censored(
            _build_delay_dist(family, log_mean_oa, log_shape_oa);
            interval = 1.0), length(d.oa_n))
    end

    if !isempty(d.ad_h)
        d.ad_h ~ Turing.filldist(double_interval_censored(
            _build_delay_dist(family, log_mean_ad + β_ad_hcw, log_shape_ad);
            interval = 1.0), length(d.ad_h))
    end
    if !isempty(d.ad_n)
        d.ad_n ~ Turing.filldist(double_interval_censored(
            _build_delay_dist(family, log_mean_ad, log_shape_ad);
            interval = 1.0), length(d.ad_n))
    end

    if !isempty(d.ac_h)
        d.ac_h ~ Turing.filldist(double_interval_censored(
            _build_delay_dist(family, log_mean_ac + β_ac_hcw, log_shape_ac);
            interval = 1.0), length(d.ac_h))
    end
    if !isempty(d.ac_n)
        d.ac_n ~ Turing.filldist(double_interval_censored(
            _build_delay_dist(family, log_mean_ac, log_shape_ac);
            interval = 1.0), length(d.ac_n))
    end

    if !isempty(d.on_h)
        d.on_h ~ Turing.filldist(double_interval_censored(
            _build_delay_dist(family, log_mean_on + β_on_hcw, log_shape_on);
            interval = 1.0), length(d.on_h))
    end
    if !isempty(d.on_n)
        d.on_n ~ Turing.filldist(double_interval_censored(
            _build_delay_dist(family, log_mean_on, log_shape_on);
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
