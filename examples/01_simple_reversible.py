"""Example 01 — Simple reversible reaction A <-> B.

Demonstrates:
- CRNetwork creation from reaction strings
- CRNT structural analysis (δ=0, DZT applies)
- Conservation law extraction
- Steady-state finding with analytical verification
- Symbolic ODE display
"""
from mantis import CRNetwork

rn = CRNetwork.from_string(
    ["A <-> B"],
    rates={"A -> B": 1.0, "B -> A": 0.5},
)

print(rn.crnt_summary())
print()

# Symbolic ODEs
odes = rn.odes(numeric_rates=True)
print("Symbolic ODEs:")
for sp, expr in sorted(odes.items()):
    print(f"  d{sp}/dt = {expr}")
print()

# Steady states from [A]=2, [B]=0 → total = 2
# Analytical: A* = 2*(0.5/1.5) = 2/3, B* = 2*(1/1.5) = 4/3
ic = {"A": 2.0, "B": 0.0}
ss_list = rn.steady_states(ic, n_attempts=20, seed=0)

print(f"Found {len(ss_list)} steady state(s):")
for i, ss in enumerate(ss_list):
    print(f"  SS {i+1}: A={ss.concentrations['A']:.6f}  B={ss.concentrations['B']:.6f}")
    print(f"         stable={ss.is_stable}  residual={ss.residual:.2e}")
    print(f"         eigenvalues: {ss.eigenvalues}")

print()
A_star = ss_list[0].concentrations["A"]
B_star = ss_list[0].concentrations["B"]
print(f"Analytical A* = {2/3:.6f},  computed = {A_star:.6f},  error = {abs(A_star - 2/3):.2e}")
print(f"Analytical B* = {4/3:.6f},  computed = {B_star:.6f},  error = {abs(B_star - 4/3):.2e}")
