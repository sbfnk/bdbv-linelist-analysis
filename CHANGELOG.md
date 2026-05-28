# Changelog

All notable changes to this project are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project aims to follow [Semantic Versioning](https://semver.org/spec/v2.0.0.html)
once a first tagged release is cut. Until then, see the `main-latest` rolling
release at
<https://github.com/epiforecasts/bdbv-linelist-analysis/releases/tag/main-latest>
for the current state of `main`.

## [Unreleased]

### Added

- Per-case shared-latent joint model promoted to the default fit, using
  Sam Abbott's bounded-primary reparametrisation (`T_primary ~ Uniform(day,
  min(day + 1, T_secondary))`) with the matching Jacobian correction. The
  natural-history identity `D_oa + D_ad = D_od` is now enforced inside the
  model rather than recovered in post-processing (#11).
- Gamma shape, scale and standard deviation per atomic delay (`oa`, `ad`,
  `ac`, `on`) are exposed in the returned `post` NamedTuple and written to
  `output/posterior_<family>.csv`, so downstream consumers (e.g. grEPI) can
  import the fitted delay distributions directly with their 95 % credible
  intervals (#12).
- Dashed reference lines for the Rosello et al. 2015 Table 5 means on both
  `plot_ppc` and `plot_family_comparison`, making the 30-day-cap finding on
  onset-to-notification visible directly on the figures (#1).
- Live headline-estimates table on the documentation home page and at the
  top of the walkthrough, generated from the single-fit pipeline so the
  rendered numbers always match the latest model run.
- Documenter site published to <https://epiforecasts.io/bdbv-linelist-analysis>
  with a rolling `main-latest` results release bundling posterior summaries
  and figures.

### Changed

- Headline tables clarify that the Rosello column is the empirical sample
  mean (SD) under the 30-day cap, not their fitted Gamma summary.
- Delay-family dispatch refactored to use singleton types instead of a
  symbol switch.
- Documentation home page is now sourced from `README.md` and the
  limitations page is lifted into the navigation.
- All repository URLs relocated from `sbfnk` to the `epiforecasts`
  organisation.

### Fixed

- Vitepress dead-link error for the `LIMITATIONS.md` cross-reference on the
  model page.
- Thread-safe sampling RNG and Gamma-only restriction on the stratified
  model.

## Release process (proposed)

Once issue #7 is resolved this section will describe the tag-based release
flow. Until then, the rolling `main-latest` release published by
`.github/workflows/Release.yml` is the canonical artefact for `main`.
