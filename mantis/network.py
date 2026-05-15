"""CRNetwork: unified interface for CRNT structural analysis, symbolic ODEs, ODE simulation, and steady-state finding."""
from __future__ import annotations

from functools import cached_property
from typing import TYPE_CHECKING

import numpy as np

from .parsing import Reaction, Complex, parse_reactions, normalize_rate_key
from .stoichiometry import (
    build_species_list,
    build_stoichiometry_matrix,
    matrix_rank,
    conservation_laws_sympy,
)
from .crnt import (
    build_complex_graph,
    get_linkage_classes,
    check_weak_reversibility,
    compute_deficiency,
    crnt_analysis,
    CRNTResult,
)

if TYPE_CHECKING:
    import sympy
    import networkx as nx
    import matplotlib.axes


class CRNetwork:
    """
    Unified interface for Chemical Reaction Network Theory (CRNT) analysis.

    Construct with `CRNetwork.from_string(reaction_strings, rates=...)`.
    Structural properties (deficiency, conservation laws, weak reversibility) are
    computed lazily and cached; they do not require rate constants.  Symbolic and
    numerical methods require rate constants to be supplied.

    Parameters
    ----------
    reactions : list[Reaction]
        Directed Reaction objects.  Each reversible string `A <-> B` produces two.
    rates : dict[str, float], optional
        Rate constants keyed by reaction strings, e.g. ``{"A -> B": 1.0}``.
        Keys are automatically normalized to canonical form (species sorted
        alphabetically); use ``rate_keys()`` to inspect the expected form.
    """

    def __init__(
        self,
        reactions: list[Reaction],
        rates: dict[str, float] | None = None,
        chemostatted: dict[str, float] | None = None,
    ) -> None:
        self._reactions = reactions
        self._rates: dict[str, float] = {}
        if rates:
            # Normalize user-supplied keys to canonical form
            for key, val in rates.items():
                try:
                    norm = normalize_rate_key(key)
                except ValueError:
                    norm = key
                self._rates[norm] = val
        self._chemostatted: dict[str, float] = dict(chemostatted) if chemostatted else {}

    @classmethod
    def from_string(
        cls,
        reaction_strings: list[str],
        rates: dict[str, float] | None = None,
        chemostatted: dict[str, float] | None = None,
    ) -> CRNetwork:
        """
        Construct from human-readable reaction strings.

        Parameters
        ----------
        reaction_strings : list[str]
            Each string is one reaction, e.g. ``"A + 2B <-> C"`` or ``"ES -> E + P"``.
        rates : dict[str, float], optional
            Maps reaction strings to rate constants.  Order of species within each
            side does not matter; keys are normalized before lookup.
        chemostatted : dict[str, float], optional
            Species held at fixed concentrations by an external reservoir.
            They are excluded from the ODE system and stoichiometry matrix rows,
            but their concentrations are folded into flux expressions.
        """
        reactions = parse_reactions(reaction_strings)
        return cls(reactions, rates, chemostatted)

    # ── Structural properties (no rates needed) ──────────────────────────────

    @cached_property
    def species(self) -> list[str]:
        return build_species_list(self._reactions, set(self._chemostatted.keys()))

    @cached_property
    def complexes(self) -> list[Complex]:
        return list(self._crnt_graph.nodes())

    @cached_property
    def stoichiometry_matrix(self) -> np.ndarray:
        return build_stoichiometry_matrix(self._reactions, self.species)

    @cached_property
    def n_species(self) -> int:
        return len(self.species)

    @cached_property
    def n_complexes(self) -> int:
        return len(self.complexes)

    @cached_property
    def n_reactions(self) -> int:
        return len(self._reactions)

    @cached_property
    def _crnt_graph(self) -> nx.DiGraph:
        return build_complex_graph(self._reactions, set(self._chemostatted.keys()))

    @cached_property
    def n_linkage_classes(self) -> int:
        return len(get_linkage_classes(self._crnt_graph))

    @cached_property
    def deficiency(self) -> int:
        return self._crnt_result.deficiency

    @cached_property
    def is_weakly_reversible(self) -> bool:
        return self._crnt_result.is_weakly_reversible

    @cached_property
    def conservation_laws(self) -> list[sympy.Expr]:
        return conservation_laws_sympy(self.stoichiometry_matrix, self.species)

    @cached_property
    def _crnt_result(self) -> CRNTResult:
        return crnt_analysis(self._reactions, self.stoichiometry_matrix, self.species, set(self._chemostatted.keys()))

    @property
    def chemostatted(self) -> dict[str, float]:
        """Chemostatted species and their fixed concentrations."""
        return dict(self._chemostatted)

    def rate_keys(self) -> list[str]:
        """Return canonical rate key strings for all reactions (for debugging)."""
        return [r.rate_key for r in self._reactions]

    def crnt_summary(self) -> str:
        """Return a human-readable CRNT analysis summary."""
        r = self._crnt_result
        lines = list(r.summary_lines)
        if self._chemostatted:
            chem_str = ", ".join(f"{s}={v}" for s, v in sorted(self._chemostatted.items()))
            lines.insert(0, f"Chemostatted species: {chem_str}")
            lines.insert(1, "")
        lines.append("")
        laws = self.conservation_laws
        lines.append(f"Conservation laws ({len(laws)}):")
        for i, law in enumerate(laws, 1):
            lines.append(f"  [{i}] {law} = const")
        return "\n".join(lines)

    # ── Symbolic ──────────────────────────────────────────────────────────────

    def odes(self, numeric_rates: bool = True) -> dict[str, sympy.Expr]:
        """
        Return mass-action ODEs as SymPy expressions keyed by species name.

        Parameters
        ----------
        numeric_rates : bool
            If True (default) and rates were supplied, substitute numerical values
            so that the returned expressions contain floats rather than ``k_i`` symbols.
        """
        from .symbolic import (
            make_species_symbols,
            make_rate_symbols,
            build_odes,
            substitute_rates,
        )
        sp_syms = make_species_symbols(self.species)
        rate_syms, key_to_sym = make_rate_symbols(self._reactions)
        odes_sym = build_odes(
            self._reactions, self.species, sp_syms, rate_syms, key_to_sym,
            chemostatted_values=self._chemostatted or None,
        )
        if numeric_rates and self._rates:
            odes_sym = {
                sp: substitute_rates(expr, key_to_sym, self._rates)
                for sp, expr in odes_sym.items()
            }
        return odes_sym

    def jacobian(self) -> sympy.Matrix:
        """Return the symbolic Jacobian matrix (n_species × n_species), cached after first call."""
        if "_jacobian_cache" not in self.__dict__:
            from .symbolic import make_species_symbols, build_jacobian
            sp_syms = make_species_symbols(self.species)
            J = build_jacobian(self.odes(numeric_rates=False), self.species, sp_syms)
            self.__dict__["_jacobian_cache"] = J
        return self.__dict__["_jacobian_cache"]

    def jacobian_at(self, ss: dict[str, float]) -> sympy.Matrix:
        from .symbolic import make_species_symbols, substitute_rates, substitute_steady_state
        J_sym = self.jacobian()
        sp_syms = make_species_symbols(self.species)
        from .symbolic import make_rate_symbols
        rate_syms, key_to_sym = make_rate_symbols(self._reactions)
        J_rates = substitute_rates(J_sym, key_to_sym, self._rates)
        return substitute_steady_state(J_rates, sp_syms, ss)

    # ── Numerical ─────────────────────────────────────────────────────────────

    def steady_states(
        self,
        initial_conditions: dict[str, float],
        n_attempts: int = 50,
        seed: int | None = None,
        t_end: float = 1e4,
    ) -> list:
        """
        Find steady states from the given initial conditions.

        Uses ODE integration (Radau) as the primary strategy, which exactly
        preserves conservation laws, then falls back to multi-start least_squares
        for additional steady states. Duplicate solutions (relative L2 distance
        < 10%) and non-physical solutions (any concentration < −tol) are discarded.
        Results are sorted by residual ascending and filtered by relative residual
        against the best solution found to remove spurious stuck-at-IC artifacts.

        Parameters
        ----------
        initial_conditions : dict[str, float]
            Initial concentrations; missing species default to 0.
        n_attempts : int
            Total number of solver attempts (including the integration from IC).
        seed : int, optional
            RNG seed for reproducibility.
        t_end : float
            ODE integration horizon in seconds (default 1e4).  Increase for
            systems with slow reactions whose timescale τ = 1/(k·[X]) >> 1e4 s.

        Returns
        -------
        list[SteadyState]
            Each entry has ``.concentrations``, ``.eigenvalues``, ``.is_stable``,
            ``.is_oscillatory``, and ``.residual``.
        """
        from .analysis import find_steady_states
        return find_steady_states(
            self._reactions,
            self.species,
            self._rates,
            initial_conditions,
            n_attempts=n_attempts,
            seed=seed,
            chemostatted_values=self._chemostatted or None,
            t_end=t_end,
        )

    def bifurcation(
        self,
        parameter: str,
        range: tuple[float, float],
        n_points: int = 100,
        initial_conditions: dict[str, float] | None = None,
        plot: bool = False,
        t_end: float = 1e4,
    ) -> object:
        """
        Scan one rate constant over a log-spaced range and collect steady states.

        Parameters
        ----------
        parameter : str
            Rate key to vary, e.g. ``"miR21 + H1 -> miR21_H1"``.
        range : tuple[float, float]
            (min, max) values for the parameter (log scale).
        n_points : int
            Number of points in the scan.
        initial_conditions : dict[str, float], optional
            Initial concentrations for each scan point (defaults to all-zero).
        plot : bool
            If True, display a bifurcation diagram for the first species.
        t_end : float
            ODE integration horizon in seconds (default 1e4).  Increase for
            systems with slow reactions whose timescale τ = 1/(k·[X]) >> 1e4 s.

        Returns
        -------
        BifurcationResult
        """
        from .analysis import scan_bifurcation
        result = scan_bifurcation(
            self._reactions,
            self.species,
            self._rates,
            parameter=parameter,
            param_range=range,
            n_points=n_points,
            initial_conditions=initial_conditions or {},
            n_attempts=20,
            chemostatted_values=self._chemostatted or None,
            t_end=t_end,
        )
        if plot:
            from .plot import plot_bifurcation
            plot_bifurcation(result, species=self.species[0])
        return result

    def simulate(
        self,
        initial_conditions: dict[str, float],
        t_span: tuple[float, float],
        t_eval=None,
        rtol: float = 1e-8,
        atol: float = 1e-12,
    ):
        """
        Integrate the ODE system forward in time and return the full trajectory.

        Parameters
        ----------
        initial_conditions : dict[str, float]
            Initial concentrations; missing species default to 0.
        t_span : (t0, tf)
            Start and end times in seconds.
        t_eval : array-like, optional
            Times at which to store the solution.  Defaults to 200 log-spaced
            points across t_span.
        rtol, atol : float
            Solver tolerances (passed to scipy Radau).

        Returns
        -------
        SimulationResult
            ``.times``, ``.concentrations`` (dict species → 1-D array),
            ``.success``.  Use ``.final()`` for the last time-point dict
            or ``.at(t)`` for a specific time.
        """
        from .analysis import simulate_ode
        return simulate_ode(
            self._reactions,
            self.species,
            self._rates,
            initial_conditions,
            t_span=t_span,
            t_eval=t_eval,
            chemostatted_values=self._chemostatted or None,
            rtol=rtol,
            atol=atol,
        )

    def draw(self, ax=None, layout: str = "spring") -> matplotlib.axes.Axes:
        from .plot import draw_reaction_graph
        return draw_reaction_graph(self._crnt_graph, ax=ax, layout=layout)
