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

const REPO_ROOT    = dirname(@__DIR__)
const LITERATE_SRC = joinpath(@__DIR__, "examples", "analysis.jl")
const LITERATE_OUT = joinpath(@__DIR__, "src")
const SNIPPET_PATH = joinpath(REPO_ROOT, "output", "cache", "headline.md")

# Single-fit pipeline. `execute = true` makes Literate actually run
# the walkthrough at pre-processing time (instead of leaving
# `@example` blocks for Documenter to evaluate later), so the fit
# happens exactly once. A `#src` side effect inside the walkthrough
# writes the headline-table Markdown to `SNIPPET_PATH`, which the
# README-substitution step below splices into the docs home page.
Literate.markdown(LITERATE_SRC, LITERATE_OUT;
    name = "analysis",
    flavor = Literate.DocumenterFlavor(),
    execute = true,
    mdstrings = true, credit = false)

let
    readme  = read(joinpath(REPO_ROOT, "README.md"), String)
    snippet = read(SNIPPET_PATH, String)
    # GitHub readers see the static prose between the marker
    # comments; the deployed home page gets the rebuilt-each-time
    # table the walkthrough just produced.
    rendered = replace(readme,
        r"<!-- HEADLINE:START -->.*?<!-- HEADLINE:END -->"s => snippet)
    write(joinpath(LITERATE_OUT, "index.md"), rendered)
end

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
        "Charniga 2024 checklist" => "charniga-checklist.md",
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
