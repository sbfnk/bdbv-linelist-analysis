using Pkg: Pkg
Pkg.instantiate()

using DocumenterVitepress
using Documenter
using Literate
using CairoMakie
using BdbvLinelist
using Printf: @sprintf
using Statistics: quantile

# Retina-quality figures in the rendered docs.
CairoMakie.activate!(; px_per_unit = 2.0)

DocMeta.setdocmeta!(
    BdbvLinelist, :DocTestSetup,
    :(using BdbvLinelist); recursive = true
)

const REPO_ROOT    = dirname(@__DIR__)
const LITERATE_SRC = joinpath(@__DIR__, "examples", "analysis.jl")
const LITERATE_OUT = joinpath(@__DIR__, "src")

# Pre-fit the canonical Gamma model so the README home page can show
# live headline estimates. The walkthrough re-fits (and runs the
# family comparison) inside its own @example blocks; this one extra
# Gamma fit adds ~1–2 min to the build but keeps the home-page table
# from drifting out of sync with the model.
let
    chn, post, _ = redirect_stdout(devnull) do
        BdbvLinelist.analyse(family = :gamma, progress = false)
    end

    fmt(v) = let q = quantile(v, [0.025, 0.5, 0.975])
        @sprintf("%.2f (%.2f – %.2f)", q[2], q[1], q[3])
    end

    # Sample sizes: regenerate from the data the fit just used so the
    # snippet doesn't depend on hard-coded counts.
    d = BdbvLinelist.build_data(BdbvLinelist.load_linelist())

    headline_md = """
    Gamma fit (WAIC-selected, doubly censored), posterior median and mean with 95% credible intervals. Rosello's 2015 Table 5 means included for direct comparison.

    | Delay                  |  n | Median, days (95% CrI) | Mean, days (95% CrI) | Rosello mean |
    | ---                    | -: | ---                    | ---                  | -:           |
    | Onset → admission      | $(length(d.onset_to_admit))     | $(fmt(post.median_oa)) | $(fmt(post.mean_oa)) | 4.00 |
    | Admission → death      | $(length(d.admit_to_death))     | $(fmt(post.median_ad)) | $(fmt(post.mean_ad)) | 7.59 |
    | Admission → discharge  | $(length(d.admit_to_discharge)) | $(fmt(post.median_ac)) | $(fmt(post.mean_ac)) | 8.00 |
    | Onset → notification   | $(length(d.onset_to_notif))     | $(fmt(post.median_on)) | $(fmt(post.mean_on)) | 8.83 |

    Rosello's onset → notification fit applied a 30-day cap; without it our posterior mean runs roughly twice as long. See the [analysis walkthrough](https://epiforecasts.io/bdbv-linelist-analysis/dev/analysis) for the full Gamma summary, posterior predictive checks and convolved marginals.
    """

    readme = read(joinpath(REPO_ROOT, "README.md"), String)
    # Replace the HEADLINE block in the README with the live table.
    # GitHub readers see the static fallback prose between the marker
    # comments; the docs home page gets the rebuilt-each-time table.
    rendered = replace(readme,
        r"<!-- HEADLINE:START -->.*?<!-- HEADLINE:END -->"s => headline_md)
    write(joinpath(LITERATE_OUT, "index.md"), rendered)
end

# Render the executable walkthrough through Literate so figures and
# tables are generated at build time from a single source script.
Literate.markdown(LITERATE_SRC, LITERATE_OUT;
    name = "analysis",
    flavor = Literate.DocumenterFlavor(),
    mdstrings = true, credit = false)

makedocs(;
    sitename = "BDBV linelist analysis",
    authors = "Sebastian Funk",
    clean = true,
    doctest = false,
    linkcheck = true,
    warnonly = [:docs_block, :missing_docs, :autodocs_block, :linkcheck],
    modules = [BdbvLinelist],
    pages = [
        "Home" => "index.md",
        "Model" => "model.md",
        "Limitations" => "limitations.md",
        "Analysis walkthrough" => "analysis.md",
        "API Reference" => "api.md"
    ],
    format = DocumenterVitepress.MarkdownVitepress(;
        repo = "github.com/epiforecasts/bdbv-linelist-analysis",
        devbranch = "main",
        devurl = "dev"
    )
)

DocumenterVitepress.deploydocs(;
    repo = "github.com/epiforecasts/bdbv-linelist-analysis",
    target = "build",
    branch = "gh-pages",
    devbranch = "main",
    push_preview = true
)
