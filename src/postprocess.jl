## Posterior summaries, model comparison, and CSV output for the
## BDBV delay + CFR model.

# ---------------------------------------------------------------------------
# Diagnostics
# ---------------------------------------------------------------------------

function _scalar_stats(summary)
    out = Float64[]
    for p in FlexiChains.parameters(summary)
        v = summary[p]
        if v isa Number
            ismissing(v) && continue
            push!(out, Float64(v))
        else
            for x in skipmissing(vec(collect(v)))
                push!(out, Float64(x))
            end
        end
    end
    return out
end

function _num_divergences(chn)
    for e in FlexiChains.extras(chn)
        e.name === :numerical_error || continue
        return Int(sum(skipmissing(vec(chn[e]))))
    end
    return 0
end

"""
$(TYPEDSIGNATURES)

Convergence diagnostics for `chn`: `(; rhat, ess, ndiv)` — max R̂
across scalar parameter entries, min bulk ESS, and divergent
transition count.
"""
function diagnostics(chn)
    rhats = _scalar_stats(FlexiChains.rhat(chn))
    esses = _scalar_stats(FlexiChains.ess(chn; kind = :bulk))
    return (; rhat = maximum(rhats), ess = minimum(esses),
            ndiv = _num_divergences(chn))
end

# ---------------------------------------------------------------------------
# Per-family delay parameter extraction
# ---------------------------------------------------------------------------

# LogNormal stores (log_median, log_sd); mean = exp(log_median + log_sd²/2).
# Gamma   stores (log_mean, log_shape); mean = exp(log_mean).
# Weibull stores (log_mean, log_shape); mean = exp(log_mean).
function _delay_summary(chn, name::Symbol, family::Symbol)
    if family == :lognormal
        μ = vec(collect(chn[_delay_varname(name, :log_median)]))
        σ = vec(collect(chn[_delay_varname(name, :log_sd)]))
        median = exp.(μ)
        mean   = exp.(μ .+ σ.^2 ./ 2)
        return (; median, mean, μ, σ)
    else
        log_mean  = vec(collect(chn[_delay_varname(name, :log_mean)]))
        log_shape = vec(collect(chn[_delay_varname(name, :log_shape)]))
        mean   = exp.(log_mean)
        # Median for Gamma/Weibull has no closed form — compute per draw.
        median = Vector{Float64}(undef, length(mean))
        for i in eachindex(mean)
            s = exp(log_shape[i])
            if family == :gamma
                median[i] = quantile(Gamma(s, mean[i] / s), 0.5)
            else
                θ = mean[i] / SpecialFunctions.gamma(1 + 1 / s)
                median[i] = quantile(Weibull(s, θ), 0.5)
            end
        end
        return (; median, mean, log_mean, log_shape)
    end
end

# VarName lookup table — submodel-prefixed by the LHS distribution
# name. Hard-coded because `@varname` is a macro requiring a literal
# expression so we cannot synthesise it from `(name, p)` at runtime.
const _DELAY_VARNAMES = Dict{Tuple{Symbol,Symbol}, Any}(
    (:dist_oa, :log_median) => @varname(dist_oa.log_median),
    (:dist_oa, :log_sd)     => @varname(dist_oa.log_sd),
    (:dist_oa, :log_mean)   => @varname(dist_oa.log_mean),
    (:dist_oa, :log_shape)  => @varname(dist_oa.log_shape),
    (:dist_ad, :log_median) => @varname(dist_ad.log_median),
    (:dist_ad, :log_sd)     => @varname(dist_ad.log_sd),
    (:dist_ad, :log_mean)   => @varname(dist_ad.log_mean),
    (:dist_ad, :log_shape)  => @varname(dist_ad.log_shape),
    (:dist_ac, :log_median) => @varname(dist_ac.log_median),
    (:dist_ac, :log_sd)     => @varname(dist_ac.log_sd),
    (:dist_ac, :log_mean)   => @varname(dist_ac.log_mean),
    (:dist_ac, :log_shape)  => @varname(dist_ac.log_shape),
    (:dist_on, :log_median) => @varname(dist_on.log_median),
    (:dist_on, :log_sd)     => @varname(dist_on.log_sd),
    (:dist_on, :log_mean)   => @varname(dist_on.log_mean),
    (:dist_on, :log_shape)  => @varname(dist_on.log_shape),
)

_delay_varname(name::Symbol, p::Symbol) = _DELAY_VARNAMES[(name, p)]

# ---------------------------------------------------------------------------
# Per-parameter summaries
# ---------------------------------------------------------------------------

qci(x; q = (0.025, 0.5, 0.975)) =
    (quantile(x, q[1]), quantile(x, q[2]), quantile(x, q[3]))

# Sample one realisation of the convolution of two distributions per
# draw, returning summary statistics across draws. The family chooses
# which distribution to reconstruct from the per-draw parameters.
function _draw_dist(family::Symbol, params, idx::Int)
    if family == :lognormal
        return LogNormal(params.μ[idx], params.σ[idx])
    elseif family == :gamma
        m = exp(params.log_mean[idx]); s = exp(params.log_shape[idx])
        return Gamma(s, m / s)
    else
        m = exp(params.log_mean[idx]); s = exp(params.log_shape[idx])
        return Weibull(s, m / SpecialFunctions.gamma(1 + 1 / s))
    end
end

function _convolve_delays(rng, family::Symbol, p_a, p_b; n_per_draw = 500)
    K = length(p_a.mean)
    means   = Vector{Float64}(undef, K)
    medians = Vector{Float64}(undef, K)
    sds     = Vector{Float64}(undef, K)
    p95     = Vector{Float64}(undef, K)
    for k in 1:K
        s = rand(rng, _draw_dist(family, p_a, k), n_per_draw) .+
            rand(rng, _draw_dist(family, p_b, k), n_per_draw)
        means[k]   = mean(s)
        medians[k] = quantile(s, 0.5)
        sds[k]     = std(s)
        p95[k]     = quantile(s, 0.95)
    end
    return (; means, medians, sds, p95)
end

# Per-draw two-component mixture of fitted delay distributions. Each of
# the `n_per_draw` realisations is drawn from `p_a` with probability
# `w[k]` (the per-draw mixing weight), otherwise from `p_b`. Used for
# the in-hospital length-of-stay marginal, where `p_a` = admit→death,
# `p_b` = admit→discharge, and `w` = in-hospital fatality among admitted
# cases. Summarised across draws like `_convolve_delays`.
function _mixture_delays(rng, family::Symbol, p_a, p_b, w; n_per_draw = 500)
    K = length(p_a.mean)
    means   = Vector{Float64}(undef, K)
    medians = Vector{Float64}(undef, K)
    sds     = Vector{Float64}(undef, K)
    p95     = Vector{Float64}(undef, K)
    for k in 1:K
        a = rand(rng, _draw_dist(family, p_a, k), n_per_draw)
        b = rand(rng, _draw_dist(family, p_b, k), n_per_draw)
        s = ifelse.(rand(rng, n_per_draw) .< w[k], a, b)
        means[k]   = mean(s)
        medians[k] = quantile(s, 0.5)
        sds[k]     = std(s)
        p95[k]     = quantile(s, 0.95)
    end
    return (; means, medians, sds, p95)
end

"""
$(TYPEDSIGNATURES)

Build the named tuple of posterior draws and print a headline
summary table. Onset → death and onset → discharge are derived in
post-processing as convolutions of the fitted onset→admit with
admit→death and admit→discharge respectively.

When the model input `d` is supplied, the in-hospital length-of-stay
marginal (admission → departure) is also derived: a per-draw mixture
of admit→death and admit→discharge, weighted by the in-hospital
fatality among admitted cases (`Beta(1 + n_died, 1 + n_discharged)`).
"""
function summarise(chn, family::Symbol; d = nothing, seed = 20260519, n_per_draw = 500)
    p_oa = _delay_summary(chn, :dist_oa, family)
    p_ad = _delay_summary(chn, :dist_ad, family)
    p_ac = _delay_summary(chn, :dist_ac, family)
    p_on = _delay_summary(chn, :dist_on, family)

    β_0   = vec(collect(chn[@varname(β_0)]))
    β_hcw = vec(collect(chn[@varname(β_hcw)]))
    β_def = vec(collect(chn[@varname(β_def)]))
    β_age = vec(collect(chn[@varname(β_age)]))

    rng = Random.MersenneTwister(seed)
    od_conv = _convolve_delays(rng, family, p_oa, p_ad; n_per_draw)
    oc_conv = _convolve_delays(rng, family, p_oa, p_ac; n_per_draw)

    # In-hospital length of stay (admission → departure): mixture of the
    # fitted admit→death and admit→discharge delays, weighted per draw by
    # the in-hospital fatality among admitted cases. Drawn after the
    # convolutions so their RNG stream is unchanged. Skipped when `d` is
    # not supplied (the weight needs the admitted-case counts).
    los_conv = nothing
    p_die    = nothing
    if d !== nothing
        n_died       = length(d.admit_to_death)
        n_discharged = length(d.admit_to_discharge)
        p_die    = rand(rng, Beta(1 + n_died, 1 + n_discharged), length(p_ad.mean))
        los_conv = _mixture_delays(rng, family, p_ad, p_ac, p_die; n_per_draw)
    end

    cfr_baseline    = _logistic.(β_0)
    cfr_hcw_conf    = _logistic.(β_0 .+ β_hcw)
    cfr_nonhcw_prob = _logistic.(β_0 .+ β_def)
    cfr_hcw_prob    = _logistic.(β_0 .+ β_hcw .+ β_def)

    println("=== Fitted atomic delay components ($(family), doubly censored) ===\n")
    _print_delay(:onset_to_admit,        p_oa)
    _print_delay(:admit_to_death,        p_ad)
    _print_delay(:admit_to_discharge,    p_ac)
    _print_delay(:onset_to_notification, p_on)

    pp_oa = _delay_params_full(family, p_oa)
    pp_ad = _delay_params_full(family, p_ad)
    pp_ac = _delay_params_full(family, p_ac)
    pp_on = _delay_params_full(family, p_on)

    println("\n=== Drop-in distribution parameters (mean, SD, Gamma shape/scale) ===\n")
    _print_delay_params(:onset_to_admit,        family, pp_oa)
    _print_delay_params(:admit_to_death,        family, pp_ad)
    _print_delay_params(:admit_to_discharge,    family, pp_ac)
    _print_delay_params(:onset_to_notification, family, pp_on)

    println("\n=== Derived (convolved) marginals ===\n")
    _print_conv(:onset_to_death,     od_conv)
    _print_conv(:onset_to_discharge, oc_conv)

    if los_conv !== nothing
        println("\n=== Length of stay in hospital (admission → departure) ===\n")
        println("  Overall = mixture of admit→death (fatal) and admit→discharge")
        println("  (survivor), weighted by the in-hospital fatality among admitted")
        println("  cases. The fatal and survivor rows repeat the atomic components.\n")
        _print_cfr(:in_hospital_fatality, p_die)
        _print_delay(:los_fatal,    p_ad)
        _print_delay(:los_survivor, p_ac)
        _print_conv(:los_overall,   los_conv)
    end

    println("\n=== Stratified case-fatality ===\n")
    _print_cfr(:CFR_baseline_nonHCW_confirmed, cfr_baseline)
    _print_cfr(:CFR_HCW_confirmed,             cfr_hcw_conf)
    _print_cfr(:CFR_nonHCW_probable,           cfr_nonhcw_prob)
    _print_cfr(:CFR_HCW_probable,              cfr_hcw_prob)

    println("\n=== Logit-scale coefficients ===\n")
    _print_logit(:beta_HCW,      β_hcw)
    _print_logit(:beta_probable, β_def)
    _print_logit(:beta_age_z,    β_age)

    base = (;
        family,
        mean_oa = p_oa.mean, median_oa = p_oa.median,
        mean_ad = p_ad.mean, median_ad = p_ad.median,
        mean_ac = p_ac.mean, median_ac = p_ac.median,
        mean_on = p_on.mean, median_on = p_on.median,
        sd_oa    = pp_oa.sds,    sd_ad    = pp_ad.sds,
        sd_ac    = pp_ac.sds,    sd_on    = pp_on.sds,
        shape_oa = pp_oa.shapes, shape_ad = pp_ad.shapes,
        shape_ac = pp_ac.shapes, shape_on = pp_on.shapes,
        scale_oa = pp_oa.scales, scale_ad = pp_ad.scales,
        scale_ac = pp_ac.scales, scale_on = pp_on.scales,
        od_mean = od_conv.means, od_median = od_conv.medians,
        od_sd   = od_conv.sds,   od_p95    = od_conv.p95,
        oc_mean = oc_conv.means, oc_median = oc_conv.medians,
        oc_sd   = oc_conv.sds,   oc_p95    = oc_conv.p95,
        β_0, β_hcw, β_def, β_age,
        cfr_baseline, cfr_hcw_conf, cfr_nonhcw_prob, cfr_hcw_prob,
    )
    los_conv === nothing && return base
    return merge(base, (;
        in_hospital_fatality = p_die,
        los_mean = los_conv.means, los_median = los_conv.medians,
        los_sd   = los_conv.sds,   los_p95    = los_conv.p95,
    ))
end

function _print_delay(label, p)
    m_lo, m_med, m_hi = qci(p.median)
    e_lo, e_med, e_hi = qci(p.mean)
    @printf("  %-22s  median %5.2f (%4.2f – %5.2f)  mean %5.2f (%4.2f – %5.2f)\n",
            string(label), m_med, m_lo, m_hi, e_med, e_lo, e_hi)
end

# Charniga Table 2 also wants SD and underlying distribution
# parameters (with CrIs) so downstream users can reconstruct the
# distribution. For each delay component compute mean, SD, Gamma
# shape and scale per draw and summarise across draws.
function _delay_params_full(family::Symbol, p)
    S = length(p.mean)
    means  = collect(p.mean)
    sds    = Vector{Float64}(undef, S)
    shapes = Vector{Float64}(undef, S)
    scales = Vector{Float64}(undef, S)
    for s in 1:S
        dist = _draw_dist(family, p, s)
        sds[s]    = std(dist)
        shapes[s] = _shape_param(family, p, s)
        scales[s] = _scale_param(family, p, s)
    end
    return (; means, sds, shapes, scales)
end

_shape_param(::Symbol, p, s) = exp(haskey(p, :log_shape) ? p.log_shape[s] : p.σ[s])
_scale_param(family::Symbol, p, s) =
    family == :lognormal ? exp(p.μ[s]) :
    family == :gamma     ? exp(p.log_mean[s]) / exp(p.log_shape[s]) :
                           exp(p.log_mean[s]) / SpecialFunctions.gamma(1 + 1 / exp(p.log_shape[s]))

function _print_delay_params(label, family::Symbol, pp)
    m = qci(pp.means);   sd = qci(pp.sds)
    sh = qci(pp.shapes); sc = qci(pp.scales)
    param_label = family == :lognormal ? "(log_μ, log_σ)" :
                  family == :gamma     ? "(shape, scale)" :
                                         "(shape, scale)"
    @printf("  %-22s  mean %5.2f (%4.2f – %5.2f)  SD %5.2f (%4.2f – %5.2f)  shape %5.2f (%4.2f – %5.2f)  scale %5.2f (%4.2f – %5.2f)\n",
            string(label), m[2], m[1], m[3], sd[2], sd[1], sd[3],
            sh[2], sh[1], sh[3], sc[2], sc[1], sc[3])
end

function _print_conv(label, conv)
    m_lo, m_med, m_hi = qci(conv.medians)
    e_lo, e_med, e_hi = qci(conv.means)
    s_lo, s_med, s_hi = qci(conv.sds)
    p_lo, p_med, p_hi = qci(conv.p95)
    @printf("  %-22s  median %5.2f (%4.2f – %5.2f)  mean %5.2f (%4.2f – %5.2f)  SD %5.2f (%4.2f – %5.2f)  P95 %5.2f (%4.2f – %5.2f)\n",
            string(label), m_med, m_lo, m_hi, e_med, e_lo, e_hi,
            s_med, s_lo, s_hi, p_med, p_lo, p_hi)
end

function _print_cfr(label, p)
    lo, med, hi = qci(p)
    @printf("  %-35s  %.3f (%.3f – %.3f)\n", string(label), med, lo, hi)
end

function _print_logit(label, β)
    lo, med, hi = qci(β)
    or_lo, or_med, or_hi = exp(lo), exp(med), exp(hi)
    @printf("  %-15s  log-OR %+.2f (%+.2f – %+.2f)   OR %.2f (%.2f – %.2f)\n",
            string(label), med, lo, hi, or_med, or_lo, or_hi)
end

# ---------------------------------------------------------------------------
# WAIC
# ---------------------------------------------------------------------------

"""
$(TYPEDSIGNATURES)

Compute WAIC directly from the per-draw delay parameters and the
CFR coefficients in the fitted chain. Independent of any
`pointwise_loglikelihoods` machinery so it works uniformly across
FlexiChains and family choices.

Returns `(; waic, lppd, p_waic, n_obs)`. Lower WAIC is better.
WAIC = −2 (lppd − p_waic) by the standard scaling.
"""
function compute_waic(chn, d, family::Symbol)
    p_oa = _delay_summary(chn, :dist_oa, family)
    p_ad = _delay_summary(chn, :dist_ad, family)
    p_ac = _delay_summary(chn, :dist_ac, family)
    p_on = _delay_summary(chn, :dist_on, family)

    β_0   = vec(collect(chn[@varname(β_0)]))
    β_hcw = vec(collect(chn[@varname(β_hcw)]))
    β_def = vec(collect(chn[@varname(β_def)]))
    β_age = vec(collect(chn[@varname(β_age)]))

    S = length(β_0)

    # Per-observation log-likelihood: rows = draws, cols = observations.
    L_blocks = Vector{Matrix{Float64}}()

    for (delays, p) in (
        (d.onset_to_admit, p_oa),
        (d.admit_to_death, p_ad),
        (d.admit_to_discharge, p_ac),
        (d.onset_to_notif, p_on),
    )
        n = length(delays)
        L = Matrix{Float64}(undef, S, n)
        for s in 1:S
            dist = _draw_dist(family, p, s)
            dic  = double_interval_censored(dist; interval = 1.0)
            for i in 1:n
                L[s, i] = logpdf(dic, delays[i])
            end
        end
        push!(L_blocks, L)
    end

    # CFR observations.
    n_out = length(d.outcome)
    L_out = Matrix{Float64}(undef, S, n_out)
    for s in 1:S
        for i in 1:n_out
            η = β_0[s] +
                β_hcw[s] * (d.hcw[i]      ? 1.0 : 0.0) +
                β_def[s] * (d.probable[i] ? 1.0 : 0.0) +
                β_age[s] * d.age_z[i]
            p_d = _logistic(η)
            L_out[s, i] = d.outcome[i] ? log(p_d) : log(1 - p_d)
        end
    end
    push!(L_blocks, L_out)

    lppd_terms   = Float64[]
    p_waic_terms = Float64[]
    for L in L_blocks
        for i in axes(L, 2)
            col = @view L[:, i]
            Lmax = maximum(col)
            lppd_i = Lmax + log(mean(exp.(col .- Lmax)))
            push!(lppd_terms, lppd_i)
            push!(p_waic_terms, var(col; corrected = true))
        end
    end

    lppd   = sum(lppd_terms)
    p_waic = sum(p_waic_terms)
    waic   = -2 * (lppd - p_waic)
    return (; waic, lppd, p_waic, n_obs = length(lppd_terms))
end

# ---------------------------------------------------------------------------
# Persistence
# ---------------------------------------------------------------------------

"""
$(TYPEDSIGNATURES)

Save the posterior summary to a CSV.
"""
function save_posterior(post, path)
    # Drop the family symbol — DataFrame doesn't accept scalars
    # mixed with vectors.
    df = DataFrame(Base.structdiff(post, NamedTuple{(:family,)}))
    mkpath(dirname(path))
    CSV.write(path, df)
end
