from .network import CRNetwork
from .parsing import Reaction, Complex
from .analysis import (
    SteadyState, BifurcationResult, SimulationResult,
    StochasticResult, gillespie_simulate, tau_leap_simulate,
)
from .crnt import CRNTResult

__all__ = [
    "CRNetwork", "Reaction", "Complex",
    "SteadyState", "BifurcationResult", "SimulationResult",
    "StochasticResult", "gillespie_simulate", "tau_leap_simulate",
    "CRNTResult",
]
from importlib.metadata import PackageNotFoundError, version as _pkg_version

try:
    # Derived from git tags at build time via setuptools-scm.
    __version__ = _pkg_version("mantis-delta")
except PackageNotFoundError:  # source tree that was never installed
    __version__ = "0.0.0+unknown"
