"""Numerical steady-state finding, stability analysis, and bifurcation scanning."""
from dataclasses import dataclass, field
from typing import Any

import numpy as np
from scipy.optimize import least_squares
from scipy.linalg import eigvals

from .parsing import Reaction
from .stoichiometry import (
    build_stoichiometry_matrix,
    conservation_laws_sympy,
    matrix_rank,
)


@dataclass
class SteadyState:
    concentrations: dict[str, float]
    eigenvalues: np.ndarray
    is_stable: bool
    is_oscillatory: bool
    residual: float


@dataclass
class BifurcationResult:
    parameter_name: str
    parameter_values: np.ndarray
    steady_states: list[list[SteadyState]]


def build_ode_function(
    reactions: list[Reaction],
    species: list[str],
    rate_values: dict[str, float],
    chemostatted: set[str] | None = None,
) -> Any:
    """
    Return a callable f(t, y) → dydt using pure numpy.
    Pre-computes N matrix and ordered rate array for speed.

    Chemostatted species are excluded from reactants_info (their contribution
    has already been folded into rate_values via fold_chemostatted_into_rates).
    """
    chem = chemostatted or set()
    sp_idx = {s: i for i, s in enumerate(species)}
    N = build_stoichiometry_matrix(reactions, species)
    rates_arr = np.array([rate_values.get(r.rate_key, 0.0) for r in reactions])

    # Pre-compute reactant info: list of (species_index, coeff) per reaction
    # Skip chemostatted species — their concentrations are already in rates_arr
    reactants_info = [
        [(sp_idx[name], coeff) for name, coeff in rxn.reactants if name not in chem]
        for rxn in reactions
    ]

    def f(t, y):
        y = np.maximum(y, 0.0)
        fluxes = np.empty(len(reactions))
        for j, (rate, ri) in enumerate(zip(rates_arr, reactants_info)):
            flux = rate
            for idx, coeff in ri:
                flux *= y[idx] ** coeff
            fluxes[j] = flux
        return N @ fluxes

    return f


def _conservation_law_vectors(
    reactions: list[Reaction],
    species: list[str],
) -> list[np.ndarray]:
    """
    Return list of non-negative integer numpy vectors spanning the left null space of N.
    Uses the same sign-normalized basis as conservation_laws_sympy() so that
    each vector has all non-negative entries (physical moiety totals).
    """
    from .stoichiometry import conservation_laws_sympy
    import sympy
    N = build_stoichiometry_matrix(reactions, species)
    laws = conservation_laws_sympy(N, species)
    sp_sym_map = {s: sympy.Symbol(s) for s in species}
    result = []
    for law in laws:
        v = np.array([
            float(law.coeff(sp_sym_map[s])) for s in species
        ], dtype=float)
        result.append(v)
    return result


def _build_augmented_system(
    reactions: list[Reaction],
    species: list[str],
    rate_values: dict[str, float],
    cl_vectors: list[np.ndarray],
    cl_totals: list[float],
    chemostatted: set[str] | None = None,
) -> Any:
    """
    Build f(y) = 0 system combining independent ODE equations + conservation constraints.
    Selects rank(N) linearly independent rows of N using QR pivoting.
    """
    from scipy.linalg import qr
    chem = chemostatted or set()
    sp_idx = {s: i for i, s in enumerate(species)}
    N = build_stoichiometry_matrix(reactions, species)
    rates_arr = np.array([rate_values.get(r.rate_key, 0.0) for r in reactions])
    reactants_info = [
        [(sp_idx[name], coeff) for name, coeff in rxn.reactants if name not in chem]
        for rxn in reactions
    ]

    rank = matrix_rank(N)
    if rank > 0:
        _, _, perm = qr(N.T, pivoting=True)
        independent_rows = perm[:rank]
        N_ind = N[independent_rows, :]
    else:
        N_ind = np.zeros((0, len(reactions)))

    def system(y):
        y_pos = np.maximum(y, 0.0)
        fluxes = np.empty(len(reactions))
        for j, (rate, ri) in enumerate(zip(rates_arr, reactants_info)):
            flux = rate
            for idx, coeff in ri:
                flux *= y_pos[idx] ** coeff
            fluxes[j] = flux
        ode_eqs = N_ind @ fluxes
        cl_eqs = np.array([
            v @ y_pos - tot
            for v, tot in zip(cl_vectors, cl_totals)
        ])
        return np.concatenate([ode_eqs, cl_eqs])

    return system


def _make_feasible_initial(
    cl_vectors: list[np.ndarray],
    cl_totals: list[float],
    n_sp: int,
    rng: np.random.Generator,
) -> np.ndarray | None:
    """
    Generate a random initial point that approximately satisfies conservation constraints.
    Uses a simple proportional allocation strategy.
    """
    y0 = rng.uniform(1e-12, 1e-7, size=n_sp)
    # Scale each conservation law to match its total
    for _ in range(20):
        for v, tot in zip(cl_vectors, cl_totals):
            mass = v @ y0
            if mass > 1e-30:
                mask = v > 0.5  # species in this law
                y0[mask] *= tot / mass
        y0 = np.maximum(y0, 0.0)
    return y0


def _full_jacobian_fn(
    reactions: list[Reaction],
    species: list[str],
    rate_values: dict[str, float],
    chemostatted: set[str] | None = None,
):
    """Return a function y → J(y) for the full ODE Jacobian dF/dy."""
    from .stoichiometry import build_stoichiometry_matrix
    chem = chemostatted or set()
    sp_idx = {s: i for i, s in enumerate(species)}
    N_full = build_stoichiometry_matrix(reactions, species)
    n_sp = len(species)
    rates_arr = np.array([rate_values.get(r.rate_key, 0.0) for r in reactions])
    reactants_info = [
        [(sp_idx[name], coeff) for name, coeff in rxn.reactants if name not in chem]
        for rxn in reactions
    ]

    def full_jac(y):
        y_pos = np.maximum(y, 0.0)
        J = np.zeros((n_sp, n_sp))
        for j, (rate, ri) in enumerate(zip(rates_arr, reactants_info)):
            flux = rate
            for idx, coeff in ri:
                flux *= y_pos[idx] ** coeff
            for k_idx, k_coeff in ri:
                if y_pos[k_idx] > 1e-30:
                    d_flux = flux * k_coeff / y_pos[k_idx]
                else:
                    d_flux = rate * k_coeff * (y_pos[k_idx] ** max(k_coeff - 1, 0))
                    for other_idx, other_coeff in ri:
                        if other_idx != k_idx:
                            d_flux *= y_pos[other_idx] ** other_coeff
                J[:, k_idx] += N_full[:, j] * d_flux
        return J

    return full_jac


def _integrate_to_ss(
    f_ode: Any,
    y0: np.ndarray,
    t_end: float = 1e6,
    rtol: float = 1e-10,
    atol: float = 1e-16,
) -> tuple[np.ndarray | None, float]:
    """
    Integrate ODE from y0 to t_end and return the final state.
    Returns (y_final, residual_norm) or (None, inf) on failure.
    """
    from scipy.integrate import solve_ivp
    try:
        sol = solve_ivp(
            f_ode,
            [0, t_end],
            y0,
            method="Radau",
            rtol=rtol,
            atol=atol,
            dense_output=False,
        )
        if not sol.success:
            return None, float("inf")
        y_final = np.maximum(sol.y[:, -1], 0.0)
        residual = float(np.linalg.norm(f_ode(0, y_final)))
        return y_final, residual
    except Exception:
        return None, float("inf")


def find_steady_states(
    reactions: list[Reaction],
    species: list[str],
    rate_values: dict[str, float],
    initial_conditions: dict[str, float],
    n_attempts: int = 50,
    tol: float = 1e-8,
    seed: int | None = None,
    chemostatted_values: dict[str, float] | None = None,
) -> list[SteadyState]:
    """
    Steady-state finder using two strategies:
    1. Primary: integrate ODE to t=1e6 s (conserves all conservation laws exactly).
    2. Secondary: multi-start least_squares for additional steady states.

    Chemostatted species are folded into effective rate constants so the ODE
    system only tracks dynamic species.
    """
    from .stoichiometry import fold_chemostatted_into_rates
    chem_vals = chemostatted_values or {}
    chem_keys = set(chem_vals.keys())
    eff_rates = fold_chemostatted_into_rates(reactions, rate_values, chem_vals) if chem_vals else rate_values

    n_sp = len(species)
    rng = np.random.default_rng(seed)

    cl_vecs = _conservation_law_vectors(reactions, species)
    y_ref = np.array([initial_conditions.get(s, 0.0) for s in species])
    cl_totals = [float(v @ y_ref) for v in cl_vecs]

    f_ode = build_ode_function(reactions, species, eff_rates, chem_keys)
    full_jac = _full_jacobian_fn(reactions, species, eff_rates, chem_keys)
    system_fn = _build_augmented_system(reactions, species, eff_rates, cl_vecs, cl_totals, chem_keys)

    collected: list[SteadyState] = []

    def _try_add(y_sol: np.ndarray, residual: float) -> None:
        if residual > 1e-4:
            return
        if np.any(y_sol < -tol):
            return
        y_sol = np.maximum(y_sol, 0.0)
        for existing in collected:
            y_ex = np.array([existing.concentrations[s] for s in species])
            if np.linalg.norm(y_sol - y_ex) / (np.linalg.norm(y_ex) + 1e-30) < 0.01:
                return
        J_num = full_jac(y_sol)
        eigs = eigvals(J_num)
        is_stable, is_osc = classify_steady_state(eigs)
        ss = SteadyState(
            concentrations=dict(zip(species, y_sol.tolist())),
            eigenvalues=eigs,
            is_stable=is_stable,
            is_oscillatory=is_osc,
            residual=residual,
        )
        collected.append(ss)

    # Strategy 1: integrate from initial conditions (guaranteed to respect CLs)
    y_final, res = _integrate_to_ss(f_ode, y_ref.copy())
    if y_final is not None:
        _try_add(y_final, res)

    # Strategy 2: multi-start least_squares for finding additional steady states
    for _ in range(n_attempts - 1):
        y0 = _make_feasible_initial(cl_vecs, cl_totals, n_sp, rng)
        if y0 is None:
            continue
        # Try ODE integration first (respects CLs)
        y_ivp, res_ivp = _integrate_to_ss(f_ode, y0, t_end=1e5)
        if y_ivp is not None and res_ivp < 1e-4:
            _try_add(y_ivp, res_ivp)
            continue
        # Fallback: direct least_squares
        try:
            result = least_squares(
                system_fn,
                y0,
                bounds=(0.0, np.inf),
                method="trf",
                ftol=tol,
                xtol=tol,
                gtol=tol,
                max_nfev=5000,
            )
            _try_add(result.x, float(np.linalg.norm(result.fun)))
        except Exception:
            pass

    return collected


def classify_steady_state(
    eigenvalues: np.ndarray,
    tol_zero: float = 1e-4,
) -> tuple[bool, bool]:
    """
    Returns (is_stable, is_oscillatory).
    Filters near-zero eigenvalues (from conservation law dimensions) before classifying.
    """
    if len(eigenvalues) == 0:
        return True, False
    max_abs = np.max(np.abs(eigenvalues))
    if max_abs < 1e-30:
        return True, False
    # Filter eigenvalues that are effectively zero relative to the largest
    significant = eigenvalues[np.abs(eigenvalues) > tol_zero * max_abs]
    if len(significant) == 0:
        return True, False
    is_stable = bool(np.all(np.real(significant) < 0))
    is_oscillatory = bool(
        np.any((np.abs(np.imag(significant)) > 1e-10) & (np.real(significant) < 0))
    )
    return is_stable, is_oscillatory


def compute_eigenvalues(J_numeric: np.ndarray) -> np.ndarray:
    return eigvals(J_numeric)


def scan_bifurcation(
    reactions: list[Reaction],
    species: list[str],
    base_rate_values: dict[str, float],
    parameter: str,
    param_range: tuple[float, float],
    n_points: int,
    initial_conditions: dict[str, float],
    n_attempts: int = 20,
    chemostatted_values: dict[str, float] | None = None,
) -> BifurcationResult:
    """Vary one rate constant over a log-spaced range and collect steady states."""
    from .parsing import normalize_rate_key
    try:
        norm_param = normalize_rate_key(parameter)
    except ValueError:
        norm_param = parameter

    param_values = np.logspace(
        np.log10(param_range[0]), np.log10(param_range[1]), n_points
    )
    all_ss: list[list[SteadyState]] = []
    for pval in param_values:
        rates = dict(base_rate_values)
        rates[norm_param] = pval
        ss_list = find_steady_states(
            reactions, species, rates, initial_conditions,
            n_attempts=n_attempts, seed=42,
            chemostatted_values=chemostatted_values,
        )
        all_ss.append(ss_list)

    return BifurcationResult(
        parameter_name=norm_param,
        parameter_values=param_values,
        steady_states=all_ss,
    )
