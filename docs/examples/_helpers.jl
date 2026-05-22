## Build the Markdown headline-table snippet that the docs build
## splices into the README home page. Kept here (not in the
## package) because it is a docs-rendering concern, not part of the
## published API. Called once from a `#src` line in `analysis.jl`.

using Printf: @sprintf
using Statistics: quantile

function headline_snippet(post, d)
    fmt(v) = let q = quantile(v, [0.025, 0.5, 0.975])
        @sprintf("%.2f (%.2f – %.2f)", q[2], q[1], q[3])
    end
    # Note: separator cells must use **≥3** dash characters — Julia's
    # Markdown.jl rejects 2-char alignment markers (`-:`, `:-`) and
    # silently falls back to parsing the whole table as a paragraph.
    return """
    All four delays fitted as **Gamma** (WAIC-selected, doubly censored). Posterior median and mean with 95% credible intervals, alongside the Rosello *et al.* 2015 Table 5 means for direct comparison (their fits, also Gamma, applied a 30-day cap on the raw delays).

    | Delay                  |  n   | Gamma median (95% CrI), days | Gamma mean (95% CrI), days | Rosello mean |
    | :---                   | ---: | :---:                        | :---:                      | ---:         |
    | Onset → admission      | $(length(d.onset_to_admit))     | $(fmt(post.median_oa)) | $(fmt(post.mean_oa)) | 4.00 |
    | Admission → death      | $(length(d.admit_to_death))     | $(fmt(post.median_ad)) | $(fmt(post.mean_ad)) | 7.59 |
    | Admission → discharge  | $(length(d.admit_to_discharge)) | $(fmt(post.median_ac)) | $(fmt(post.mean_ac)) | 8.00 |
    | Onset → notification   | $(length(d.onset_to_notif))     | $(fmt(post.median_on)) | $(fmt(post.mean_on)) | 8.83 |

    Rosello's onset → notification fit applied a 30-day cap; without it our posterior mean runs roughly twice as long. See the [analysis walkthrough](https://epiforecasts.io/bdbv-linelist-analysis/dev/analysis) for the full Gamma summary, posterior predictive checks and convolved marginals.
    """
end
