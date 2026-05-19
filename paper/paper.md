---
title: 'Mantis-Delta: Mass-Action Network Theory and Steady-State Characterization for Chemical Reaction Networks'
tags:
  - Python
  - chemical reaction network theory
  - mass-action kinetics
  - deficiency theorems
  - steady-state analysis
  - bifurcation analysis
  - DNA computing
  - systems biology
  - scientific Python
authors:
  - name: Emilio A. Venegas
    orcid: 0000-0002-7689-9185
    equal-contrib: true
    corresponding: true
    affiliation: 1
affiliations:
 - name: Tecnológico de Monterrey, Escuela de Medicina y Ciencias de la Salud, Ave. Eugenio Garza Sada 2501 Sur, Col. Tecnológico, Monterrey, N.L. 64700, México
   index: 1
date: 18 May 2026
bibliography: references.bib
---

# Summary

Chemical Reaction Network Theory (CRNT), developed by Horn, Jackson, and Feinberg [@horn1972general; @feinberg1972complex; @feinberg1987chemical], provides parameter-free structural theorems that constrain the asymptotic dynamics of mass-action systems regardless of the numerical values of the rate constants. These theorems allow researchers to determine, from the network topology alone, whether a biochemical system can exhibit multiple steady states, sustained oscillations, or is guaranteed to converge to a unique equilibrium. Despite the maturity of the theory, modern open-source implementations that combine CRNT structural analysis with symbolic ordinary differential equation (ODE) construction and robust numerical steady-state solving remain scarce.

`mantis-delta` is a pure-Python library that fills this gap. It ingests human-readable reaction strings (e.g., `"A + B -> C"`), builds the complex reaction graph, computes the network deficiency $\delta = n - \ell - s$, determines weak reversibility, and decides the applicability of the Deficiency Zero Theorem (DZT) and Deficiency One Theorem (D1T). For systems satisfying these structural conditions, `mantis-delta` certifies---without any simulation---the existence, uniqueness, and (for DZT) asymptotic stability of the positive steady state in every stoichiometric compatibility class. When structural theorems do not apply, the library constructs symbolic mass-action ODEs and Jacobians via SymPy and provides a hybrid numerical solver that combines stiff implicit integration with bound-constrained algebraic least-squares to locate both stable and unstable fixed points, including Hopf bifurcation centres inaccessible to forward integration.

# Statement of need

The dynamics of a biochemical reaction network under mass-action kinetics are governed by polynomial ODEs whose qualitative behaviour depends jointly on the network topology and on rate constants that are seldom known to better than one order of magnitude in living systems. Experimentalists designing enzymatic cycles, DNA strand-displacement cascades [@yin2008catalytic; @li2011nonenzymatic], or metabolic motifs [@ang2013considerations] need to certify---before any wet-lab work---whether a proposed network can exhibit pathological bistability or oscillation. Simultaneously, educators and researchers in CRNT need a compact, fully open-source reference implementation against which to develop further structural algorithms.

`mantis-delta` targets both communities. By coupling Feinberg's structural theorems [@feinberg1987chemical; @feinberg1995existence; @feinberg2019foundations] to symbolic and numerical tools, the library allows practitioners to (i) state programmatically whether a network is structurally prohibited from exhibiting bistability or oscillation, (ii) extract symbolic ODE and Jacobian expressions for downstream pipelines such as parameter inference or sensitivity analysis, and (iii) locate both stable and unstable fixed points across the full parameter regime, including chemostatted (open) networks where conservation laws are broken and forward integration alone cannot reach Hopf-unstable equilibria [@polettini2014irreversible].

# State of the field

The practical bridge between formal CRNT theorems and computational pipelines remains under-served. The Chemical Reaction Network Toolbox [@ji2011crnreal], a MATLAB-based graphical tool developed by Feinberg's group, implements the core deficiency computations but requires a commercial MATLAB license, does not expose a programmatic API, and has not been updated since 2011. The constraint-programming approach of @soliman2012crnpy provides invariant detection for biochemical models but does not implement the deficiency theorems directly and does not provide ODE construction or numerical steady-state solving. General-purpose reaction network simulators such as COPASI, BioNetGen, and PySB focus on simulation and parameter estimation rather than structural certification; none decide the DZT or D1T.

`mantis-delta` was built rather than contributed to an existing tool because no available package provides the specific combination of capabilities required: (a) an idiomatic Python API compatible with the modern scientific stack (NumPy, SciPy, SymPy, NetworkX), (b) end-to-end coverage from reaction parsing through structural theorem certification to numerical steady-state finding, (c) first-class support for chemostatted networks including algebraic location of Hopf-unstable fixed points, and (d) symbolic ODE/Jacobian expressions ready for downstream analysis such as thermodynamic constraint checking [@ederer2007thermodynamically]. This combination constitutes the scholarly contribution of `mantis-delta`: it makes CRNT accessible as a composable building block within Python-based research workflows.

# Software design

`mantis-delta` is organised into seven modules totalling approximately 2,500 lines of Python. The `parsing` module converts human-readable reaction strings into canonical complex representations. The `network` module constructs the `CRNetwork` object, which maintains the complex reaction graph (via NetworkX [@hagberg2008exploring]), the stoichiometric matrix, and methods for computing deficiency, linkage classes, weak reversibility, and conservation laws. The `crnt` module implements the decision procedures for the DZT and D1T, including the per-linkage-class deficiency checks required by D1T [@feinberg1995existence; @feinberg2019foundations]. The `symbolic` module constructs mass-action ODE and Jacobian expressions through SymPy [@meurer2017sympy], supporting both symbolic and numeric rate constants. The `analysis` module provides the hybrid steady-state solver and eigenvalue-based stability classification. The `stoichiometry` module computes conservation laws as non-negative integer left null-space vectors via exact rational arithmetic. The `plot` module provides reaction-graph visualisation and bifurcation scan plotting.

The key design trade-off was between adopting an existing simulation framework (e.g., extending PySB or Tellurium) and building a standalone library. We chose the standalone approach because the structural CRNT analysis---exact rational rank computation, graph-theoretic theorem decision, and algebraic steady-state solving---is fundamentally different from the forward-simulation pipelines those frameworks provide. The standalone design also keeps the dependency footprint minimal (NumPy [@harris2020array], SciPy [@virtanen2020scipy], SymPy, NetworkX, Matplotlib) and allows the library to serve as a composable component in larger workflows.

The hybrid solver merits particular attention. For closed networks with conservation laws, it first attempts forward integration using SciPy's stiff implicit Radau IIA integrator [@hairer1996solving], then runs bound-constrained algebraic least-squares seeded at both the user-supplied initial condition and random points on the conservation manifold. For open (chemostatted) networks, where conservation laws are broken, only the algebraic strategy is used. This two-stage approach reliably finds unstable fixed points---such as the Hopf-unstable equilibrium of the chemostatted Brusselator [@prigogine1968symmetry]---that pure forward integration would orbit indefinitely without converging to.

# Research impact statement

`mantis-delta` has been validated on six benchmark networks of increasing complexity, all bundled as reproducible example scripts: a reversible isomerisation (DZT applies), the Michaelis--Menten enzyme mechanism [@michaelis1913kinetik] ($\delta = 0$, not weakly reversible), the closed and chemostatted Brusselator (Hopf bifurcation), a catalytic hairpin assembly (CHA) miR-21 biosensor [@yin2008catalytic] ($\delta = 1$, D1T applies), and the Goldbeter--Koshland zero-order ultrasensitivity switch [@goldbeter1981amplified]. In every case, the CRNT-predicted qualitative behaviour (monostability, oscillation, uniqueness) is recovered numerically with a residual below $10^{-6}\,\mathrm{M\,s^{-1}}$, and the Goldbeter--Koshland dose-response curve agrees with the closed-form quasi-steady-state approximation to within $1\%$ over a $400\times$ kinase/phosphatase activity scan. The CHA miR-21 biosensor analysis, in particular, demonstrates the practical value of structural certification for DNA nanotechnology circuit design, where the D1T guarantee of a unique equilibrium eliminates concerns about hysteresis in dose-response curves. The library is accompanied by a 108-test pytest suite covering parsing, CRNT structural analysis, symbolic ODE construction, chemostatted networks, and numerical steady-state solving, and comprehensive API documentation generated via MkDocs.

![Brusselator dynamics with chemostatted $A$ and $B$. Top row ($B = 3$, oscillatory): (a) time series of $X$ and $Y$; (b) phase portrait showing the limit cycle (green) orbiting the unstable focus at $(1, 3)$ located by `mantis-delta`'s algebraic solver (black marker). Bottom row ($B = 1.5$, damped): (c) time series; (d) phase portrait spiralling into the stable equilibrium $(1, 1.5)$. The Hopf bifurcation threshold $B = 1 + A^2$ separates the two regimes. Generated by `examples/05_brusselator_chemostatted.py`.](../examples/brusselator_combined.png){#fig-bruss width=95%}

# AI usage disclosure

The core algorithms, software architecture, and numerical experiments in `mantis-delta` were designed and implemented by the author without generative AI assistance. AI coding assistants were used in a limited capacity during test development: specifically, for generating synthetic reaction networks as test fixtures, suggesting edge-case scenarios (e.g., degenerate stoichiometric matrices, networks at the boundary of theorem applicability), and drafting boilerplate test scaffolding. All AI-generated test code and data were manually reviewed, validated against known analytical results, and revised as needed before inclusion in the test suite. No AI tools were used in the writing of this manuscript.

# Acknowledgements

The author acknowledges the open-source scientific Python ecosystem---NumPy [@harris2020array], SciPy [@virtanen2020scipy], SymPy [@meurer2017sympy], and NetworkX [@hagberg2008exploring]---without which this work would not have been possible. This research received no specific grant from any funding agency in the public, commercial, or not-for-profit sectors.

# References
