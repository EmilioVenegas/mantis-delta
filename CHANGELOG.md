# Changelog

All notable changes to mantis-delta are documented here.
This file is generated from the git history by [git-cliff](https://git-cliff.org).
The format follows [Keep a Changelog](https://keepachangelog.com) and the project
uses [Semantic Versioning](https://semver.org).

## [Unreleased]

### Documentation

- Add deterministic ODE simulation documentation and update project metadata

## [0.2.0] - 2026-05-31

### Documentation

- Add bioRxiv preprint manuscript in Quarto (.qmd) format

    Add a research article showcasing pycrn for chemical reaction network
    analysis, formatted as a bioRxiv preprint with professional scientific
    typography. Includes structural analysis with LaTeX equations, the six
    example workflows from examples/, a bibliography spanning the foundational
    CRNT literature (Horn-Jackson, Feinberg) through to applied work
    (Goldbeter-Koshland, NUPACK, scientific Python stack), and bioRxiv
    submission metadata (article category, subject area, abstract, author
    approvals, competing interests, funding, distribution/reuse).

- Render pycrn_preprint.pdf via XeLaTeX / Quarto (16 pages, 458 KB)

    Fix two rendering issues on headless Ubuntu:
    - Replace Times New Roman with Liberation Serif (not installed in env)
    - Remove Zotero CSL URL (blocked); use pandoc default citation style

    PDF renders clean with all figures, numbered equations, callout boxes
    for DZT and D1T, 30-entry bibliography, and full bioRxiv front-matter.

- Replace em-dashes with commas for improved readability in pycrn_preprint.qmd
- Add JOSS paper to the project repository

### Features

- Add time-course simulation, integrate STIXGeneral font, and expose custom integration horizons.
- Implement Gillespie and τ-leap stochastic simulation methods with associated tests

### Refactor

- Improve steady-state solver to use algebraic methods for open systems and update stability classification logic.
- Improve steady-state filtering and add publication-quality plotting utilities for reaction graphs and Brusselator dynamics
- Rename project from mantis to mantis-delta across documentation, examples, and configuration files

## [0.1.0] - 2026-05-14

### Features

- Add support for chemostatted species with fixed concentration values in ODE and steady-state analysis

### Refactor

- Rename package from crnpy to pycrn


