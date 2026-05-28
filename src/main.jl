## CLI entry point, analyse() driver, and multi-family comparison.

const FAMILIES = (:lognormal, :gamma, :weibull)

"""
$(TYPEDSIGNATURES)

Run NUTS on `model` using ForwardDiff and per-chain prior init.

A parent `MersenneTwister(seed)` is used to draw `chains` child
seeds, which are passed as an explicit RNG vector to `sample`. This
avoids mutating the global RNG and keeps results reproducible
regardless of unrelated `Random` calls between runs.
"""
function sample_fit(model;
        samples = 1000,
        chains = 4,
        target_accept = 0.95,
        seed = 20260519,
        progress = false,
    )
    adtype = AutoForwardDiff()
    parent = MersenneTwister(seed)
    rng = MersenneTwister(rand(parent, UInt64))
    return sample(
        rng,
        model,
        NUTS(target_accept; adtype),
        MCMCThreads(),
        samples, chains;
        initial_params = fill(DynamicPPL.InitFromPrior(), chains),
        progress = progress,
    )
end

"""
$(TYPEDSIGNATURES)

Load the line list, fit the BDBV delay + CFR model under one
distribution `family`, print the summary, and save the posterior to
`output/posterior_<family>.csv`. Returns `(chain, post, diag)`.
"""
function analyse(;
        data = LINELIST_PATH,
        output = OUTPUT_DIR,
        family = :gamma,
        samples = 1000,
        chains = 4,
        seed = 20260519,
        target_accept = 0.95,
        progress = true,
    )
    family in FAMILIES ||
        throw(ArgumentError("family must be one of $(FAMILIES); got $family"))

    ll = load_linelist(data)
    d  = build_data(ll)
    @info "Loaded line list" n_cases=d.N n_onset_admit=length(d.onset_to_admit) n_admit_death=length(d.admit_to_death) n_admit_disch=length(d.admit_to_discharge) n_onset_notif=length(d.onset_to_notif) n_hcw=sum(d.hcw)

    chn = sample_fit(bdbv_model(d; family);
        samples = samples, chains = chains, seed = seed,
        target_accept = target_accept, progress = progress)

    diag = diagnostics(chn)
    @info "Diagnostics" family max_rhat=round(diag.rhat, digits = 3) min_ess=round(Int, diag.ess) n_divergent=diag.ndiv

    println("\n=== Posterior summary [$family] (median, 95% CrI) ===\n")
    post = summarise(chn, family; d = d)
    save_posterior(post, joinpath(output, "posterior_$(family).csv"))

    return chn, post, diag
end

"""
$(TYPEDSIGNATURES)

Death-pathway mixture analysis. Fits

1. the main `bdbv_model` (for `oa`, `ad` and the admit-pathway
   convolution `od_admit = oa ⊛ ad`),
2. a separate doubly-censored Gamma for the community-died
   onset → death distribution `od_comm`,
3. a Beta-Binomial posterior on `p_admit = N_admit_died /
   (N_admit_died + N_comm_died)`.

The population-mixture marginal is then

   od_mix = p_admit · od_admit  +  (1 − p_admit) · od_comm

sampled per-draw and summarised. Reports all three (community-only,
admit-pathway, mixture) side-by-side along with their interpretation
ranges in operational terms (early outbreak vs after ETC).
"""
function fit_death_mixture(;
        data = LINELIST_PATH,
        family = :gamma,
        samples = 1000,
        chains = 4,
        seed = 20260519,
        target_accept = 0.95,
        progress = false,
        n_per_draw = 500,
    )
    ll = load_linelist(data); d = build_data(ll)
    @info "Death mixture inputs" n_admit_died=d.n_admit_died n_comm_died=d.n_comm_died comm_delays=d.onset_to_comm_death

    chn_main = sample_fit(bdbv_model(d; family);
        samples = samples, chains = chains, seed = seed,
        target_accept = target_accept, progress = progress)
    chn_comm = sample_fit(community_death_model(d.onset_to_comm_death;
            family = family,
            n_admit_died = d.n_admit_died,
            n_comm_died  = d.n_comm_died);
        samples = samples, chains = chains, seed = seed + 7,
        target_accept = target_accept, progress = progress)
    # Drop-in single-Gamma fit to all 27 onset→death pairs — the
    # population-level equivalent of Rosellos direct empirical fit,
    # for downstream users who want one distribution.
    # No counts passed: p_admit defaults to Beta(1, 1) prior. The
    # share is already estimated in chn_comm; only the delay fit from
    # chn_all is used downstream.
    chn_all = sample_fit(community_death_model(d.onset_to_death_all;
            family = family);
        samples = samples, chains = chains, seed = seed + 17,
        target_accept = target_accept, progress = progress)

    # Atomic-delay parameters (admit-pathway components).
    p_oa = _delay_summary(chn_main, :dist_oa, family)
    p_ad = _delay_summary(chn_main, :dist_ad, family)

    # Community-pathway: only one delay distribution.
    p_cd = if family == :lognormal
        μ = vec(collect(chn_comm[@varname(dist_cd.log_median)]))
        σ = vec(collect(chn_comm[@varname(dist_cd.log_sd)]))
        m = exp.(μ .+ σ.^2 ./ 2)
        (; μ, σ, mean = m, log_mean = μ, log_shape = σ)
    else
        lm = vec(collect(chn_comm[@varname(dist_cd.log_mean)]))
        ls = vec(collect(chn_comm[@varname(dist_cd.log_shape)]))
        (; log_mean = lm, log_shape = ls, mean = exp.(lm))
    end
    p_admit = vec(collect(chn_comm[@varname(p_admit)]))

    S = min(length(p_oa.mean), length(p_cd.mean), length(p_admit))

    rng = Random.MersenneTwister(seed)
    means_admit = Vector{Float64}(undef, S)
    means_comm  = Vector{Float64}(undef, S)
    means_mix   = Vector{Float64}(undef, S)
    medians_mix = Vector{Float64}(undef, S)

    for s in 1:S
        oa_draws  = rand(rng, _draw_dist(family, p_oa, s), n_per_draw)
        ad_draws  = rand(rng, _draw_dist(family, p_ad, s), n_per_draw)
        admit_od  = oa_draws .+ ad_draws
        comm_od   = rand(rng, _draw_dist(family, p_cd, s), n_per_draw)
        means_admit[s] = mean(admit_od)
        means_comm[s]  = mean(comm_od)

        # Mixture: each of the n_per_draw realisations is drawn from
        # admit-pathway with prob p_admit[s], else community pathway.
        z = rand(rng, n_per_draw) .< p_admit[s]
        mix = ifelse.(z, admit_od, comm_od)
        means_mix[s]   = mean(mix)
        medians_mix[s] = quantile(mix, 0.5)
    end

    function _qci(x, label)
        q = quantile(x, [0.025, 0.5, 0.975])
        @printf("  %-30s  median %5.2f   mean %5.2f   95%% CrI [%5.2f – %5.2f]\n",
                label, q[2], mean(x), q[1], q[3])
    end

    println("\n=== Death-pathway mixture (Isiro 2012, $family) ===\n")
    @printf("  Pathway shares: N_admit_died = %d, N_community_died = %d\n",
            d.n_admit_died, d.n_comm_died)
    _qci(p_admit, "p_admit (admit pathway share)")
    println()
    _qci(means_admit, "admit-pathway od mean (oa+ad)")
    _qci(means_comm,  "community-pathway od mean")
    _qci(means_mix,   "mixture marginal od mean")
    _qci(medians_mix, "mixture marginal od median")

    # Drop-in single-Gamma fit summaries — mean, SD, shape, scale
    # with 95% CrIs.
    println("\n=== Drop-in single distributions (Gamma parameters) ===\n")
    println("  Use these for downstream EpiNow2 / scenario modelling. The")
    println("  appropriate fit depends on the assumed care-access regime:")
    println()
    for (label, chn) in (
        ("ALL DEATHS (population-avg, Rosello-equivalent)", chn_all),
        ("Community deaths only (early-phase / no ETC)",   chn_comm),
    )
        if family == :gamma
            lm = vec(collect(chn[@varname(dist_cd.log_mean)]))
            ls = vec(collect(chn[@varname(dist_cd.log_shape)]))
            pp = (; mean = exp.(lm), log_mean = lm, log_shape = ls)
            full = _delay_params_full(family, pp)
            println("  $label:")
            _print_delay_params(:onset_to_death, family, full)
        end
    end

    println("\n  Operational reading:")
    println("    Early outbreak / pre-ETC:     use community-pathway above")
    println("    After ETC operational:        use admit-pathway convolution (oa ⊛ ad)")
    println("    Population-average / steady:  use ALL DEATHS drop-in above")

    return (; chn_main, chn_comm, chn_all, p_admit,
            means_admit, means_comm, means_mix, medians_mix)
end

"""
$(TYPEDSIGNATURES)

Prior sensitivity sweep: fit the canonical (gamma) model under
three prior-scale settings (tight 0.5 / default 1.0 / wide 2.0)
and tabulate the atomic delay means with their 95% CrIs.
"""
function sensitivity(;
        data = LINELIST_PATH,
        family = :gamma,
        samples = 1000,
        chains = 4,
        seed = 20260519,
        target_accept = 0.95,
        progress = false,
        scales = (0.5, 1.0, 2.0),
    )
    ll = load_linelist(data); d = build_data(ll)
    results = Dict{Float64, Any}()
    for s in scales
        @info "Prior sensitivity" prior_scale=s
        chn = sample_fit(bdbv_model(d; family = family, prior_scale = s);
            samples = samples, chains = chains, seed = seed,
            target_accept = target_accept, progress = progress)
        results[s] = (; chain = chn, diag = diagnostics(chn))
    end

    println("\n=== Prior sensitivity (means of atomic delays, 95% CrI) ===\n")
    println("                    prior scale 0.5         prior scale 1.0         prior scale 2.0")
    for (sym, label) in (
        (:dist_oa, "onset→admit"),
        (:dist_ad, "admit→death"),
        (:dist_ac, "admit→discharge"),
        (:dist_on, "onset→notif"),
    )
        row = "  " * rpad(label, 18)
        for s in scales
            p = _delay_summary(results[s].chain, sym, family)
            q = quantile(p.mean, [0.025, 0.5, 0.975])
            row *= @sprintf("  %5.2f (%4.2f – %5.2f)   ", q[2], q[1], q[3])
        end
        println(row)
    end
    println()
    return results
end

"""
$(TYPEDSIGNATURES)

Fit all three families and print a side-by-side comparison of the
headline atomic delay means alongside model-comparison via WAIC
(per [`compute_waic`](@ref)).

Returns a Dict keyed by family symbol with `(; chain, post, diag,
waic)` for each.
"""
function compare_families(;
        data = LINELIST_PATH,
        output = OUTPUT_DIR,
        samples = 1000,
        chains = 4,
        seed = 20260519,
        target_accept = 0.95,
        progress = false,
    )
    ll = load_linelist(data)
    d  = build_data(ll)
    @info "Loaded line list" n_cases=d.N

    results = Dict{Symbol, Any}()
    for fam in FAMILIES
        @info "Fitting family" family=fam
        model = bdbv_model(d; family = fam)
        chn = sample_fit(model;
            samples = samples, chains = chains, seed = seed,
            target_accept = target_accept, progress = progress)
        diag = diagnostics(chn)
        waic = compute_waic(chn, d, fam)
        results[fam] = (; chain = chn, diag = diag, waic = waic)
        @info "Fit done" family=fam rhat=round(diag.rhat, digits=3) ess=round(Int, diag.ess) ndiv=diag.ndiv waic=round(waic.waic, digits=1) p_waic=round(waic.p_waic, digits=1)
    end

    println("\n=== Family comparison ===\n")
    @printf("  %-10s  %10s  %10s  %8s  %6s  %6s\n",
            "family", "WAIC", "ΔWAIC", "p_waic", "R̂", "ESS")
    waics = [results[f].waic.waic for f in FAMILIES]
    best  = minimum(waics)
    for (i, f) in enumerate(FAMILIES)
        r = results[f]
        @printf("  %-10s  %10.1f  %+10.1f  %8.1f  %6.3f  %6d\n",
                string(f), r.waic.waic, r.waic.waic - best, r.waic.p_waic,
                r.diag.rhat, round(Int, r.diag.ess))
    end

    return results
end

function (@main)(args)
    s = ArgParseSettings(; description = "Fit BDBV Isiro 2012 delay + CFR model")
    @add_arg_table! s begin
        "--data", "-d"
            help = "path to line-list CSV"
            default = LINELIST_PATH
        "--output", "-o"
            help = "output directory for posterior CSVs"
            default = OUTPUT_DIR
        "--family", "-f"
            help = "delay distribution family (lognormal | gamma | weibull)"
            default = "gamma"
        "--samples", "-n"
            help = "NUTS samples per chain"
            arg_type = Int
            default = 1000
        "--chains", "-c"
            help = "number of parallel chains"
            arg_type = Int
            default = 4
        "--seed", "-s"
            help = "random seed"
            arg_type = Int
            default = 20260519
        "--compare"
            help = "fit all three families and compare via WAIC"
            action = :store_true
    end
    p = parse_args(args, s)
    if p["compare"]
        compare_families(;
            data = p["data"], output = p["output"],
            samples = p["samples"], chains = p["chains"], seed = p["seed"],
            progress = false)
    else
        analyse(;
            data = p["data"], output = p["output"],
            family = Symbol(p["family"]),
            samples = p["samples"], chains = p["chains"], seed = p["seed"],
            progress = false)
    end
    return 0
end
