## Posterior predictive figures.
##
## For each delay component we produce a panel showing:
##  - the observed integer-delay histogram (bars),
##  - the posterior predictive day-bin counts under the fitted
##    double-censoring model (median + 95% band per day),
##  - the latent fitted continuous distribution (median + 95% band).
##
## Posterior predictive bars are simulated by drawing from the fitted
## distribution per posterior draw, applying primary + secondary
## interval censoring (the actual observation process), and binning
## to integer days — matching the `double_interval_censored` likelihood.

const _DELAY_LABELS = (;
    oa = "onset → admission",
    ad = "admission → death",
    ac = "admission → discharge",
    on = "onset → notification",
)

# Simulate a single double-interval-censored observation given a
# latent distribution: primary event ~ Uniform(0, 1), latent delay
# from `dist`, secondary observed in a 1-day bin. Returns the
# observed integer delay.
function _sim_one(rng, dist)
    primary = rand(rng)
    latent  = rand(rng, dist)
    secondary = primary + latent
    return floor(Int, secondary) - 0  # primary day = 0; observed bin = floor(secondary)
end

# Posterior predictive: for each draw, simulate `n_obs` observations
# and aggregate counts per integer day bin. Returns
# (; bin_lo, bin_hi, ppc_lo, ppc_med, ppc_hi) — vectors of length
# (max_bin - min_bin + 1) giving the 2.5/50/97.5% posterior bin
# counts at each day.
function _ppc_bins(rng, S::Int, draws_to_dist, n_obs::Int; max_bin = 100,
                   n_sim_per_draw = 5)
    bins_per_draw = Matrix{Float64}(undef, S, max_bin + 1)
    bins = Vector{Int}(undef, max_bin + 1)
    for s in 1:S
        dist = draws_to_dist(s)
        fill!(bins, 0)
        for _ in 1:n_sim_per_draw, _ in 1:n_obs
            x = _sim_one(rng, dist)
            0 <= x <= max_bin || continue
            bins[x + 1] += 1
        end
        @views bins_per_draw[s, :] .= bins ./ n_sim_per_draw
    end
    bin_centres = collect(0:max_bin)
    ppc_lo  = [quantile(@view(bins_per_draw[:, b + 1]), 0.025) for b in bin_centres]
    ppc_med = [quantile(@view(bins_per_draw[:, b + 1]), 0.5)   for b in bin_centres]
    ppc_hi  = [quantile(@view(bins_per_draw[:, b + 1]), 0.975) for b in bin_centres]
    return (; bin_centres, ppc_lo, ppc_med, ppc_hi)
end

# Posterior density of the latent continuous distribution at a grid
# of points. Returns (; xs, dens_lo, dens_med, dens_hi).
function _latent_density(draws_to_dist; xmax = 100.0, n_x = 200)
    S = length(draws_to_dist)
    xs = range(0.01, xmax; length = n_x)
    pdfs = Matrix{Float64}(undef, S, n_x)
    for s in 1:S
        dist = draws_to_dist(s)
        for (j, x) in pairs(xs)
            pdfs[s, j] = pdf(dist, x)
        end
    end
    dens_lo  = [quantile(pdfs[:, j], 0.025) for j in 1:n_x]
    dens_med = [quantile(pdfs[:, j], 0.5)   for j in 1:n_x]
    dens_hi  = [quantile(pdfs[:, j], 0.975) for j in 1:n_x]
    return (; xs = collect(xs), dens_lo, dens_med, dens_hi)
end

# Build a draws_to_dist closure for the named delay component.
# Returns (closure, S = number of posterior draws).
function _draws_to_dist(chn, name::Symbol, family::Symbol)
    p = _delay_summary(chn, name, family)
    return (s -> _draw_dist(family, p, s), length(p.mean))
end

"""
$(TYPEDSIGNATURES)

Posterior predictive check figure: four panels (one per atomic
delay), each showing the observed histogram, the simulated
day-binned posterior predictive (median + 95% band), and the
fitted continuous latent density (median + 95% band, on a secondary
y-axis).
"""
function plot_ppc(chn, d, family::Symbol; seed = 20260519)
    rng = Random.MersenneTwister(seed)

    panels = [
        (name = :dist_oa, label = _DELAY_LABELS.oa, obs = d.onset_to_admit,     xmax = 15.0),
        (name = :dist_ad, label = _DELAY_LABELS.ad, obs = d.admit_to_death,     xmax = 30.0),
        (name = :dist_ac, label = _DELAY_LABELS.ac, obs = d.admit_to_discharge, xmax = 30.0),
        (name = :dist_on, label = _DELAY_LABELS.on, obs = d.onset_to_notif,     xmax = 90.0),
    ]

    fig = Figure(; size = (1000, 700))
    for (k, p) in enumerate(panels)
        row, col = fldmod1(k, 2)
        ax = Axis(fig[row, col]; title = String(p.label),
                  xlabel = "delay (days)", ylabel = "count")
        ax.titlefont = :regular
        n_obs = length(p.obs)
        xmax_int = ceil(Int, p.xmax)

        # Observed histogram (counts per integer day).
        obs_counts = zeros(Int, xmax_int + 1)
        for x in p.obs
            xi = Int(round(x))
            0 <= xi <= xmax_int && (obs_counts[xi + 1] += 1)
        end
        barplot!(ax, 0:xmax_int, obs_counts;
                 color = (:steelblue, 0.5), strokecolor = :steelblue,
                 strokewidth = 1, label = "observed")

        # Posterior predictive bins.
        draws_to_dist, S = _draws_to_dist(chn, p.name, family)
        ppc = _ppc_bins(rng, S, draws_to_dist, n_obs; max_bin = xmax_int,
                        n_sim_per_draw = 5)
        idx = ppc.bin_centres .<= xmax_int
        band!(ax, ppc.bin_centres[idx], ppc.ppc_lo[idx], ppc.ppc_hi[idx];
              color = (:firebrick, 0.25))
        lines!(ax, ppc.bin_centres[idx], ppc.ppc_med[idx];
               color = :firebrick, linewidth = 2, label = "posterior predictive")

        xlims!(ax, -0.5, p.xmax)
        if k == 1
            axislegend(ax; position = :rt)
        end
    end
    rowsize!(fig.layout, 1, Relative(0.5))
    rowsize!(fig.layout, 2, Relative(0.5))

    return fig
end

"""
$(TYPEDSIGNATURES)

Side-by-side posterior predictive comparison across families. One
column per family, four rows (one per delay). Used to make the
WAIC comparison visually concrete.
"""
function plot_family_comparison(chns_by_family, d; seed = 20260519)
    rng = Random.MersenneTwister(seed)
    families = collect(keys(chns_by_family))
    panels = [
        (name = :dist_oa, label = _DELAY_LABELS.oa, obs = d.onset_to_admit,     xmax = 15.0),
        (name = :dist_ad, label = _DELAY_LABELS.ad, obs = d.admit_to_death,     xmax = 30.0),
        (name = :dist_ac, label = _DELAY_LABELS.ac, obs = d.admit_to_discharge, xmax = 30.0),
        (name = :dist_on, label = _DELAY_LABELS.on, obs = d.onset_to_notif,     xmax = 90.0),
    ]

    fig = Figure(; size = (300 * length(families), 800))
    family_colour = Dict(:lognormal => :navy, :gamma => :firebrick, :weibull => :seagreen)

    for (row, p) in enumerate(panels)
        xmax_int = ceil(Int, p.xmax)
        n_obs = length(p.obs)
        for (col, fam) in enumerate(families)
            ax = Axis(fig[row, col]; xlabel = "delay (days)", ylabel = "count")
            if row == 1
                ax.title = String(fam)
                ax.titlefont = :bold
            end
            if col == 1
                ax.ylabel = String(p.label) * "\ncount"
            end

            obs_counts = zeros(Int, xmax_int + 1)
            for x in p.obs
                xi = Int(round(x))
                0 <= xi <= xmax_int && (obs_counts[xi + 1] += 1)
            end
            barplot!(ax, 0:xmax_int, obs_counts;
                     color = (:steelblue, 0.5), strokecolor = :steelblue, strokewidth = 1)

            chn = chns_by_family[fam]
            draws_to_dist, S = _draws_to_dist(chn, p.name, fam)
            ppc = _ppc_bins(rng, S, draws_to_dist, n_obs; max_bin = xmax_int,
                            n_sim_per_draw = 5)
            idx = ppc.bin_centres .<= xmax_int
            band!(ax, ppc.bin_centres[idx], ppc.ppc_lo[idx], ppc.ppc_hi[idx];
                  color = (family_colour[fam], 0.25))
            lines!(ax, ppc.bin_centres[idx], ppc.ppc_med[idx];
                   color = family_colour[fam], linewidth = 2)
            xlims!(ax, -0.5, p.xmax)
        end
    end

    return fig
end

"""
$(TYPEDSIGNATURES)

Epidemic curve: weekly onset counts with HCW subcounts stacked.
"""
function plot_epi_curve(ll)
    wc = weekly_onset_counts(ll)
    fig = Figure(; size = (900, 360))
    ax = Axis(fig[1, 1];
              title = "BDBV Isiro 2012: weekly onset counts",
              xlabel = "Week (ISO Monday)", ylabel = "Cases by onset",
              xticklabelrotation = 0.5)
    ax.titlefont = :regular

    week_ints  = Dates.value.(wc.week_start)
    counts     = wc.count
    hcw_counts = wc.hcw_count
    nonhcw     = counts .- hcw_counts

    # Non-HCW bars (base layer), then HCW stacked on top.
    barplot!(ax, week_ints, nonhcw;
             color = (:steelblue, 0.7), strokecolor = :steelblue, strokewidth = 1,
             label = "non-HCW")
    barplot!(ax, week_ints, hcw_counts;
             color = (:firebrick, 0.7), strokecolor = :firebrick, strokewidth = 1,
             offset = nonhcw, label = "HCW")

    # Date tick labels.
    tick_idx = 1:2:length(wc.week_start)
    ax.xticks = (week_ints[tick_idx], string.(wc.week_start[tick_idx]))

    axislegend(ax; position = :rt)
    return fig
end

"""
$(TYPEDSIGNATURES)

Save a Makie figure to `path` (PNG by default).
"""
function save_figure(fig, path)
    mkpath(dirname(path))
    CairoMakie.save(path, fig)
    return path
end
