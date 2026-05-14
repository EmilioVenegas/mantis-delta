"""Example 02 — Michaelis-Menten enzyme kinetics.

System:   E + S <-> ES -> E + P
Demonstrates:
- CRNT analysis (δ=0, not weakly reversible — DZT stability applies only for WR)
- Enzyme conservation law: E + ES = E_total
- Quasi-steady-state approximation validation via numerical SS
- Eigenvalue analysis showing stable steady state
"""
from crnpy import CRNetwork

rn = CRNetwork.from_string(
    ["E + S <-> ES", "ES -> E + P"],
    rates={
        "E + S -> ES": 1e6,
        "ES -> E + S": 1e3,
        "ES -> E + P": 100.0,
    },
)

print(rn.crnt_summary())
print()

E0 = 1e-6   # 1 µM enzyme
S0 = 1e-4   # 100 µM substrate
ic = {"E": E0, "S": S0, "ES": 0.0, "P": 0.0}

ss_list = rn.steady_states(ic, n_attempts=30, seed=0)
print(f"Found {len(ss_list)} steady state(s):")
ss = ss_list[0]
c = ss.concentrations
print(f"  E  = {c['E']:.4e} M")
print(f"  S  = {c['S']:.4e} M")
print(f"  ES = {c['ES']:.4e} M")
print(f"  P  = {c['P']:.4e} M")
print()

E_total = c["E"] + c["ES"]
print(f"E conservation: E + ES = {E_total:.4e} M  (E0 = {E0:.4e} M,  error = {abs(E_total-E0)/E0*100:.3f}%)")
print(f"Stable: {ss.is_stable}  |  Residual: {ss.residual:.2e}")

# QSSA: [ES]* = E0*[S]/(Km+[S]) where Km = (k_r + k_cat)/k_f
Km = (1e3 + 100.0) / 1e6
ES_qssa = E0 * c["S"] / (Km + c["S"])
print(f"\nQSSA check: [ES]_QSSA = {ES_qssa:.4e} M  |  computed [ES] = {c['ES']:.4e} M")
