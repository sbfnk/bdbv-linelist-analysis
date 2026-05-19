## Line-list loading, outlier scrub, and per-pair delay extraction
## for the 2012 Isiro Bundibugyo outbreak (Rosello et al. 2015
## eLife supplement 1, n = 52).

const LINELIST_PATH = joinpath(pkgdir(@__MODULE__), "data", "linelist.csv")
const OUTPUT_DIR    = joinpath(pkgdir(@__MODULE__), "output")
const FIGURES_DIR   = joinpath(pkgdir(@__MODULE__), "figures")

# Admission-date outliers from the Rosello deposit (encoding errors:
# −89, −5, −4, −1, 328720 d offset from onset). The deposit is
# canonical so we set these admission dates to missing rather than
# chasing a fix.
const ADMIT_DELAY_OUTLIERS = (-89, -5, -4, -1, 328720)

"""
$(TYPEDSIGNATURES)

Load the 2012 Isiro BDBV line list cleaned from the Rosello 2015
eLife supplement. See top of `data.jl` for the outlier scrub.
"""
function load_linelist(path = LINELIST_PATH)
    ll = CSV.read(path, DataFrame; missingstring = ["NA", ""])

    for col in (:Date_of_onset_symp, :Date_hospital_discharge,
                :Date_of_notification, :Date_of_Hospitalisation,
                :Date_disease_ended, :Date_of_Death)
        ll[!, col] = _to_date.(ll[!, col])
    end

    raw_admit_delay = [
        (ismissing(a) || ismissing(o)) ? missing : Dates.value(a - o)
        for (a, o) in zip(ll.Date_of_Hospitalisation, ll.Date_of_onset_symp)
    ]
    for i in eachindex(raw_admit_delay)
        d = raw_admit_delay[i]
        if !ismissing(d) && d in ADMIT_DELAY_OUTLIERS
            ll.Date_of_Hospitalisation[i] = missing
        end
    end

    ll.is_hcw = [!ismissing(o) && (o == "HCW" || o == "possible HCW")
                 for o in ll.Occupation]

    return ll
end

_to_date(x::Date) = x
_to_date(x::Missing) = missing
_to_date(x::AbstractString) = Date(x)
_to_date(x) = ismissing(x) ? missing : Date(string(x))

"""
$(TYPEDSIGNATURES)

Tabulate onset-date counts per ISO week for the loaded line list.
Returns a `DataFrame` with columns `(week_start, count, hcw_count)`.
"""
function weekly_onset_counts(ll)
    onsets = collect(skipmissing(ll.Date_of_onset_symp))
    hcw_onsets = ll.Date_of_onset_symp[ll.is_hcw]
    hcw_onsets = collect(skipmissing(hcw_onsets))

    min_week = Dates.firstdayofweek(minimum(onsets))
    max_week = Dates.firstdayofweek(maximum(onsets))
    weeks = collect(min_week:Day(7):max_week)
    counts     = [count(o -> Dates.firstdayofweek(o) == w, onsets) for w in weeks]
    hcw_counts = [count(o -> Dates.firstdayofweek(o) == w, hcw_onsets) for w in weeks]
    return DataFrame(week_start = weeks, count = counts, hcw_count = hcw_counts)
end

# Integer day-difference for each case where both dates exist and the
# delay is non-negative. Used to scrub biologically-impossible
# encoding errors at the same time as collecting per-pair delays.
function _pair_delays(start_dates, end_dates)
    out = Float64[]
    for i in eachindex(start_dates)
        s = start_dates[i]; e = end_dates[i]
        (ismissing(s) || ismissing(e)) && continue
        δ = Float64(Dates.value(e - s))
        δ >= 0 && push!(out, δ)
    end
    return out
end

# Same as `_pair_delays` but also returns a parallel vector of HCW
# indicators (Bool) so we can stratify by HCW status downstream.
function _pair_delays_with_hcw(start_dates, end_dates, is_hcw)
    out = Float64[]
    hcw = Bool[]
    for i in eachindex(start_dates)
        s = start_dates[i]; e = end_dates[i]
        (ismissing(s) || ismissing(e)) && continue
        δ = Float64(Dates.value(e - s))
        δ >= 0 || continue
        push!(out, δ)
        push!(hcw, Bool(is_hcw[i]))
    end
    return out, hcw
end

"""
$(TYPEDSIGNATURES)

Build the model input from a cleaned line-list `DataFrame`.

Returns a named tuple with:

- four vectors of integer day-delays (one per atomic component):
  `onset_to_admit`, `admit_to_death`, `admit_to_discharge`,
  `onset_to_notif`. Each vector contains only cases for which both
  endpoint dates are observed and the delay is non-negative.
- CFR covariates over all `N = 52` cases: `outcome` (Bool died),
  `hcw`, `probable`, `age_z`.
- `N`: total number of cases.

The marginal onset → death and onset → discharge are NOT in the data
tuple — they are derived in post-processing from the fitted
component distributions.
"""
function build_data(ll)
    onset = ll.Date_of_onset_symp
    admit = ll.Date_of_Hospitalisation
    death = ll.Date_of_Death
    disch = ll.Date_hospital_discharge
    notif = ll.Date_of_notification

    onset_to_admit,     hcw_oa = _pair_delays_with_hcw(onset, admit, ll.is_hcw)
    admit_to_death,     hcw_ad = _pair_delays_with_hcw(admit, death, ll.is_hcw)
    admit_to_discharge, hcw_ac = _pair_delays_with_hcw(admit, disch, ll.is_hcw)
    onset_to_notif,     hcw_on = _pair_delays_with_hcw(onset, notif, ll.is_hcw)

    outcome  = [o == "Dead" for o in ll.Outcome]
    hcw      = ll.is_hcw
    probable = [d == "Probable" for d in ll.Case_definition]

    age_raw = [ismissing(a) ? missing : Float64(a) for a in ll.Age]
    age_obs = collect(skipmissing(age_raw))
    ā = mean(age_obs); s̄ = std(age_obs)
    age_z = [ismissing(a) ? 0.0 : (a - ā) / s̄ for a in age_raw]

    # Pre-split each delay vector by HCW status for the stratified
    # model — Turing's `~` treats *local* variables as parameters, so
    # the observation vectors need to be data-tuple fields.
    oa_h = onset_to_admit[hcw_oa];      oa_n = onset_to_admit[.!hcw_oa]
    ad_h = admit_to_death[hcw_ad];      ad_n = admit_to_death[.!hcw_ad]
    ac_h = admit_to_discharge[hcw_ac];  ac_n = admit_to_discharge[.!hcw_ac]
    on_h = onset_to_notif[hcw_on];      on_n = onset_to_notif[.!hcw_on]

    # Community-pathway deaths: died without recorded admission.
    # Used by the death-pathway mixture analysis. Some of these were
    # also scrubbed of an outlier admission date in load_linelist; we
    # treat them all as community-pathway because we can't distinguish
    # genuine community cases from admission-date-missing-but-admitted.
    is_community_died = [
        o == "Dead" && ismissing(a)
        for (o, a) in zip(ll.Outcome, ll.Date_of_Hospitalisation)
    ]
    onset_to_comm_death, _ = _pair_delays_with_hcw(
        ll.Date_of_onset_symp[is_community_died],
        ll.Date_of_Death[is_community_died],
        ll.is_hcw[is_community_died],
    )

    n_admit_died = sum(.!ismissing.(admit) .& .!ismissing.(death))
    n_comm_died  = length(onset_to_comm_death)

    # All onset → death pairs (admit-pathway + community), for the
    # drop-in Rosello-equivalent single-distribution fit.
    onset_to_death_all = _pair_delays(onset, death)

    return (;
        onset_to_admit, admit_to_death, admit_to_discharge, onset_to_notif,
        hcw_oa, hcw_ad, hcw_ac, hcw_on,
        oa_h, oa_n, ad_h, ad_n, ac_h, ac_n, on_h, on_n,
        onset_to_comm_death, n_admit_died, n_comm_died,
        onset_to_death_all,
        outcome, hcw, probable, age_z,
        N = nrow(ll),
    )
end
