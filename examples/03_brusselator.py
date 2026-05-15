"""Example 03 — Brusselator oscillator (open system, all species dynamic).

Reactions:
  A -> X           (k1)
  B + X -> Y + D   (k2)
  2X + Y -> 3X     (k3)
  X -> E           (k4)

Classic Brusselator parameters: A=1, B=3, k1=k2=k3=k4=1.
With all species dynamic (A,B,D,E not chemostatted), δ=0 via DZT.

Note: the oscillatory behaviour observed in the classical Brusselator arises
when A and B are treated as fixed (chemostatted) parameters. Here we model the
full closed system; A and B deplete over time and the system relaxes to a fixed
point. This example shows CRNT analysis and steady-state finding for the
full open-system description.

Demonstrates:
- Deficiency calculation for the all-species open system (δ=0)
- Symbolic ODE generation for a 3rd-order network
- Steady-state finding for a nonlinear system
"""
from mantis import CRNetwork

rn = CRNetwork.from_string(
    [
        "A -> X",
        "B + X -> Y + D",
        "2X + Y -> 3X",
        "X -> E",
    ],
    rates={
        "A -> X": 1.0,
        "B + X -> Y + D": 1.0,
        "2X + Y -> 3X": 1.0,
        "X -> E": 1.0,
    },
)

print(rn.crnt_summary())
print()

print("Symbolic ODEs:")
odes = rn.odes(numeric_rates=True)
for sp, expr in sorted(odes.items()):
    print(f"  d{sp}/dt = {expr}")
print()

# Initial conditions: start with A=10, B=5, rest small
ic = {"A": 10.0, "B": 5.0, "X": 0.01, "Y": 0.01, "D": 0.0, "E": 0.0}
ss_list = rn.steady_states(ic, n_attempts=30, seed=0)
print(f"Found {len(ss_list)} steady state(s):")
for i, ss in enumerate(ss_list):
    print(f"  SS {i+1}:")
    for sp in sorted(ss.concentrations):
        print(f"    {sp} = {ss.concentrations[sp]:.6f}")
    print(f"    stable={ss.is_stable}  oscillatory={ss.is_oscillatory}  residual={ss.residual:.2e}")
