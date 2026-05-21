module BdbvLinelist

using ArgParse: ArgParseSettings, @add_arg_table!, parse_args
using CSV: CSV
using Chain: @chain
using DataFrames: DataFrame, nrow, eachrow, passmissing, sort!, dropmissing
using DataFramesMeta: @select, @transform, @subset, @combine, @by, @rtransform,
                      @rsubset, @orderby, @rename
using Dates: Dates, Date, Day
using Distributions: Normal, LogNormal, Gamma, Weibull, Bernoulli, Beta,
                     Uniform, truncated, logpdf, cdf
using SpecialFunctions: SpecialFunctions
using CensoredDistributions: double_interval_censored
using Makie: Figure, Axis, Relative, barplot!, lines!, band!, vlines!,
             xlims!, rowsize!, axislegend
using CairoMakie: CairoMakie
using DocStringExtensions: TYPEDSIGNATURES
using Printf: @printf, @sprintf
using Random: Random, MersenneTwister
using Statistics: quantile, mean, std, var
using Turing: Turing, @model, NUTS, MCMCThreads, MCMCSerial, sample, DynamicPPL,
              to_submodel, @varname
using ADTypes: AutoForwardDiff
import FlexiChains

include("data.jl")
include("model.jl")
include("postprocess.jl")
include("plots.jl")
include("main.jl")

export load_linelist, build_data, atomic_delays
export bdbv_model, bdbv_model_stratified, community_death_model
export DelayFamily, LogNormalDelay, GammaDelay, WeibullDelay,
       delay_family, family_symbol, delay_prior, build_delay_dist
export sample_fit, analyse, compare_families, sensitivity, fit_death_mixture
export summarise, save_posterior, compute_waic
export diagnostics
export plot_ppc, plot_family_comparison, plot_epi_curve, weekly_onset_counts, save_figure

end
