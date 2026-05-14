"""Example 04 — Catalytic Hairpin Assembly (CHA) miR-21 cascade.

Uses kinetic parameters derived from NUPACK thermodynamic simulations of the
CHA biosensor targeting miR-21 (22-nt microRNA biomarker).

Demonstrates:
- Full CRNT analysis of a multi-step catalytic cascade (δ=1, Deficiency One Theorem)
- 4 conservation laws (H1-total, H2-total, CP-total, miR21-total)
- Steady-state finding with conservation law constraints
- Eigenvalue stability analysis
- Bifurcation scan over miR21 concentration
- Reaction graph visualization
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")   # headless backend for non-interactive use
import matplotlib.pyplot as plt

from crnpy import CRNetwork

# ── Reaction strings ────────────────────────────────────────────────────────
CHA_STRINGS = [
    "miR21 + H1 <-> miR21_H1",
    "miR21_H1 + H2 <-> H1H2 + miR21",
    "H1H2 + CP <-> H1H2_CP",
    "H1 + H2 <-> H1H2",
]

# Rate constants from NUPACK thermodynamic simulations (37°C, physiological salt)
CHA_RATES = {
    "miR21 + H1 -> miR21_H1": 3.0e5,       # k_on (M⁻¹s⁻¹)
    "miR21_H1 -> miR21 + H1": 2.687e-3,    # k_off (s⁻¹)
    "miR21_H1 + H2 -> H1H2 + miR21": 1.0e6,
    "H1H2 + miR21 -> miR21_H1 + H2": 3.383e-3,
    "H1H2 + CP -> H1H2_CP": 2.0e5,
    "H1H2_CP -> H1H2 + CP": 5.79,
    "H1 + H2 -> H1H2": 2.379,
    "H1H2 -> H1 + H2": 7.208e-17,          # essentially irreversible
}

CHA_IC = {
    "H1":      100e-9,
    "H2":      100e-9,
    "CP":      100e-9,
    "miR21":   10e-9,
    "miR21_H1": 0.0,
    "H1H2":    0.0,
    "H1H2_CP": 0.0,
}

# ── 1. Build network ─────────────────────────────────────────────────────────
print("=" * 60)
print("CHA miR-21 Catalytic Hairpin Assembly")
print("=" * 60)
rn = CRNetwork.from_string(CHA_STRINGS, rates=CHA_RATES)

# ── 2. CRNT structural analysis ──────────────────────────────────────────────
print()
print(rn.crnt_summary())

# ── 3. Conservation laws ─────────────────────────────────────────────────────
print()
print("Conservation laws:")
for i, law in enumerate(rn.conservation_laws, 1):
    print(f"  [{i}] {law} = const")

# ── 4. Symbolic ODEs ─────────────────────────────────────────────────────────
print()
print("Symbolic ODEs (with numeric rates):")
odes = rn.odes(numeric_rates=True)
for sp in sorted(odes):
    print(f"  d{sp}/dt = {odes[sp]}")

# ── 5. Steady-state analysis ─────────────────────────────────────────────────
print()
print("Solving for steady state ([miR21]₀ = 10 nM) ...")
ss_list = rn.steady_states(CHA_IC, n_attempts=30, seed=0)
print(f"  Found {len(ss_list)} steady state(s)")
ss = ss_list[0]

print()
print("Steady-state concentrations:")
for sp in sorted(ss.concentrations):
    print(f"  {sp:12s} = {ss.concentrations[sp]:.4e} M  ({ss.concentrations[sp]*1e9:.3f} nM)")

print()
print(f"  Stable:      {ss.is_stable}")
print(f"  Oscillatory: {ss.is_oscillatory}")
print(f"  Residual:    {ss.residual:.2e}")

print()
print("Eigenvalues (significant only):")
eigs = ss.eigenvalues
max_abs = np.max(np.abs(eigs))
sig = eigs[np.abs(eigs) > 1e-4 * max_abs]
for e in sorted(sig, key=lambda x: x.real):
    print(f"  {e.real:+.4e}  {'+' if e.imag >= 0 else '-'}  {abs(e.imag):.4e}i")

# ── 6. Conservation law verification ────────────────────────────────────────
c = ss.concentrations
print()
print("Conservation law verification at SS:")
print(f"  miR21 total: {(c['miR21'] + c['miR21_H1'])*1e9:.3f} nM  (expected 10.0 nM)")
print(f"  CP total:    {(c['CP'] + c['H1H2_CP'])*1e9:.3f} nM  (expected 100.0 nM)")

# ── 7. Bifurcation scan ──────────────────────────────────────────────────────
print()
print("Running bifurcation scan over miR21 initial concentration ...")
print("(scanning initial [miR21] from 0.1 nM to 100 nM — this may take ~1 min)")

mir21_values = np.logspace(-10, -7, 20)  # 0.1 nM → 100 nM
signal_values = []

for mir21_conc in mir21_values:
    ic_scan = dict(CHA_IC)
    ic_scan["miR21"] = float(mir21_conc)
    ic_scan["miR21_H1"] = 0.0
    ss_scan = rn.steady_states(ic_scan, n_attempts=10, seed=42)
    if ss_scan:
        signal_values.append(ss_scan[0].concentrations.get("H1H2_CP", 0.0))
    else:
        signal_values.append(0.0)

fig, ax = plt.subplots(figsize=(7, 4))
ax.loglog(mir21_values * 1e9, np.maximum(signal_values, 1e-15) * 1e9, "o-", color="#2196F3")
ax.set_xlabel("[miR21]₀ (nM)")
ax.set_ylabel("[H1H2_CP]* (nM)")
ax.set_title("CHA Signal vs. miR21 Concentration")
ax.grid(True, which="both", alpha=0.3)
plt.tight_layout()
import pathlib
_here = pathlib.Path(__file__).parent

out_path = _here / "cha_bifurcation.png"
fig.savefig(out_path, dpi=150)
print(f"  Bifurcation plot saved to {out_path}")

# ── 8. Reaction graph ────────────────────────────────────────────────────────
fig2, ax2 = plt.subplots(figsize=(10, 6))
rn.draw(ax=ax2)
ax2.set_title("CHA Reaction Graph")
out_path2 = _here / "cha_reaction_graph.png"
fig2.savefig(out_path2, dpi=150, bbox_inches="tight")
print(f"  Reaction graph saved to {out_path2}")

print()
print("Done.")
