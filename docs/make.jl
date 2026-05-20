using Pkg: Pkg
Pkg.instantiate()

using DocumenterVitepress
using Documenter
using Literate
using CairoMakie
using BdbvLinelist

# Retina-quality figures in the rendered docs.
CairoMakie.activate!(; px_per_unit = 2.0)

DocMeta.setdocmeta!(
    BdbvLinelist, :DocTestSetup,
    :(using BdbvLinelist); recursive = true
)

# Render the executable walkthrough through Literate so figures and
# tables are generated at build time from a single source script.
const LITERATE_SRC = joinpath(@__DIR__, "examples", "analysis.jl")
const LITERATE_OUT = joinpath(@__DIR__, "src")
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
