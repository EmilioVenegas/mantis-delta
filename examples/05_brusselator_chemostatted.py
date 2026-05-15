"""Example 05 — Brusselator with chemostatted species (A and B held fixed).

In the full open-system Brusselator (example 03), all species are dynamic: A and B
deplete over time and the system relaxes to a stable fixed point (δ=0, DZT predicts
this). The *classical* Brusselator, however, assumes A and B are continuously
replenished — i.e. they are chemostatted at fixed concentrations. In that setting,
the 2D (X, Y) subsystem can undergo a Hopf bifurcation and sustain limit-cycle
oscillations.

This example demonstrates how to use mantis-delta's ODE machinery with chemostatted
species by wrapping the library's ODE function: species A and B are pinned to their
initial values at every time step, effectively treating them as kinetic parameters.

Demonstrates:
- CRNT analysis of the full 6-species Brusselator (δ=0, DZT: not WR)
- Chemostatting via ODE function wrapping (no library modification needed)
- Limit-cycle oscillations for B > 1 + A²  (Hopf condition)
- Time series and phase portrait plotting
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
matplotlib.rcParams.update({"font.family": "STIXGeneral", "mathtext.fontset": "stix"})
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

from mantis import CRNetwork
from mantis.analysis import build_ode_function

# ── 1. Full network CRNT analysis ───────────────────────────────────────────
rn = CRNetwork.from_string(
    [
        "A -> X",
        "B + X -> Y + D",
        "2X + Y -> 3X",
        "X -> E",
    ],
    rates={
        "A -> X":       1.0,
        "B + X -> Y + D": 1.0,
        "2X + Y -> 3X": 1.0,
        "X -> E":       1.0,
    },
)

print("=" * 60)
print("Full Brusselator — all species dynamic")
print("=" * 60)
print(rn.crnt_summary())

print()
print("Note: DZT is not applicable because the network is not weakly reversible.")
print("With all species dynamic, A and B deplete and the system reaches a fixed point.")
print()

# ── 2. Chemostatted setup ───────────────────────────────────────────────────
# Classical Brusselator parameters: A=1, B=3, k=1.
# Hopf condition: B > 1 + A²  →  3 > 1 + 1 = 2  ✓  (oscillatory)
A_val = 1.0
B_val = 3.0

# Species sorted alphabetically: [A, B, D, E, X, Y]
print(f"Species order: {rn.species}")
A_idx = rn.species.index("A")
B_idx = rn.species.index("B")

f_full = build_ode_function(rn._reactions, rn.species, rn._rates)

def f_chemostatted(t, y):
    """ODE with A and B pinned — treats them as time-invariant parameters."""
    y_pinned = y.copy()
    y_pinned[A_idx] = A_val
    y_pinned[B_idx] = B_val
    dy = f_full(t, y_pinned)
    dy[A_idx] = 0.0   # A is replenished externally
    dy[B_idx] = 0.0   # B is replenished externally
    return dy

# Fixed point of the 2D subsystem: X* = A, Y* = B/A
X_star = A_val
Y_star = B_val / A_val
print(f"Chemostatted fixed point:  X* = {X_star:.2f},  Y* = {Y_star:.2f}")
print(f"Hopf condition B > 1 + A²: {B_val:.1f} > {1 + A_val**2:.1f}  →  oscillatory\n")

# ── 3. Integrate ─────────────────────────────────────────────────────────────
# Start slightly off the fixed point to seed the oscillation.
# y0 order: [A, B, D, E, X, Y]
y0 = np.zeros(len(rn.species))
y0[A_idx] = A_val
y0[B_idx] = B_val
y0[rn.species.index("X")] = X_star + 0.5   # perturb X
y0[rn.species.index("Y")] = Y_star          # Y at fixed point

t_span = (0, 40.0)
t_eval = np.linspace(*t_span, 4000)

sol = solve_ivp(f_chemostatted, t_span, y0, method="Radau",
                t_eval=t_eval, rtol=1e-10, atol=1e-12)

X_sol = sol.y[rn.species.index("X")]
Y_sol = sol.y[rn.species.index("Y")]

print(f"Integration {'succeeded' if sol.success else 'FAILED'}")
print(f"Final X = {X_sol[-1]:.4f},  Y = {Y_sol[-1]:.4f}")
print(f"X amplitude (last 10 s):  [{X_sol[sol.t > 30].min():.3f}, {X_sol[sol.t > 30].max():.3f}]")
print()

# ── 4. Plot ──────────────────────────────────────────────────────────────────
import pathlib
_here = pathlib.Path(__file__).parent

fig, axes = plt.subplots(1, 2, figsize=(11, 4))

# Left: time series
ax1 = axes[0]
ax1.plot(sol.t, X_sol, color="#2196F3", label="X")
ax1.plot(sol.t, Y_sol, color="#F44336", label="Y")
ax1.axhline(X_star, color="#2196F3", linestyle=":", linewidth=0.8, alpha=0.6)
ax1.axhline(Y_star, color="#F44336", linestyle=":", linewidth=0.8, alpha=0.6)
ax1.set_xlabel("Time (s)")
ax1.set_ylabel("Concentration")
ax1.set_title("Time series  (A=1, B=3 chemostatted)")
ax1.legend()
ax1.set_xlim(0, 40)

# Right: phase portrait
ax2 = axes[1]
# Show transient in grey, limit cycle in colour
transient_mask = sol.t < 10
lc_mask = sol.t >= 10
ax2.plot(X_sol[transient_mask], Y_sol[transient_mask], color="grey",
         linewidth=0.8, alpha=0.5, label="Transient")
ax2.plot(X_sol[lc_mask], Y_sol[lc_mask], color="#4CAF50",
         linewidth=1.5, label="Limit cycle")
ax2.scatter([X_star], [Y_star], color="black", s=50, zorder=5, label=f"Fixed point ({X_star},{Y_star})")
ax2.set_xlabel("X")
ax2.set_ylabel("Y")
ax2.set_title("Phase portrait")
ax2.legend(framealpha=0.7)

plt.tight_layout()
out_path = _here / "brusselator_oscillations.png"
fig.savefig(out_path, dpi=150)
print(f"Plot saved to {out_path}")

# ── 5. Contrast: B < Hopf threshold (stable spiral) ────────────────────────
print()
print("─" * 60)
print("Contrast: B=1.5 < 1+A²=2  →  stable spiral (no oscillations)")
print("─" * 60)

B_stable = 1.5
def f_stable(t, y):
    y_pinned = y.copy()
    y_pinned[A_idx] = A_val
    y_pinned[B_idx] = B_stable
    dy = f_full(t, y_pinned)
    dy[A_idx] = 0.0
    dy[B_idx] = 0.0
    return dy

y0_s = y0.copy()
y0_s[B_idx] = B_stable
sol_s = solve_ivp(f_stable, (0, 60), y0_s, method="Radau",
                  t_eval=np.linspace(0, 60, 6000), rtol=1e-10, atol=1e-12)

X_s = sol_s.y[rn.species.index("X")]
Y_s = sol_s.y[rn.species.index("Y")]
print(f"B=1.5 fixed point: X*=1, Y*={B_stable/A_val:.2f}")
print(f"Final X = {X_s[-1]:.4f}  (converges to {A_val:.2f})")
print(f"Final Y = {Y_s[-1]:.4f}  (converges to {B_stable/A_val:.2f})")

fig2, axes2 = plt.subplots(1, 2, figsize=(11, 4))
axes2[0].plot(sol_s.t, X_s, color="#2196F3", label="X")
axes2[0].plot(sol_s.t, Y_s, color="#F44336", label="Y")
axes2[0].set_xlabel("Time (s)")
axes2[0].set_ylabel("Concentration")
axes2[0].set_title(f"Time series  (A=1, B={B_stable} — stable spiral)")
axes2[0].legend()

axes2[1].plot(X_s, Y_s, color="#FF9800", linewidth=1.5)
axes2[1].scatter([A_val], [B_stable/A_val], color="black", s=50, zorder=5)
axes2[1].set_xlabel("X")
axes2[1].set_ylabel("Y")
axes2[1].set_title("Phase portrait  (spirals inward)")
plt.tight_layout()
out_path2 = _here / "brusselator_stable.png"
fig2.savefig(out_path2, dpi=150)
print(f"Plot saved to {out_path2}")

# ── 6. Combined publication figure (2x2: oscillatory + stable) ─────────────
fig3, ax = plt.subplots(2, 2, figsize=(11, 7.5))

# Top row: oscillatory (B=3)
ax[0, 0].plot(sol.t, X_sol, color="#2196F3", label="X")
ax[0, 0].plot(sol.t, Y_sol, color="#F44336", label="Y")
ax[0, 0].axhline(X_star, color="#2196F3", linestyle=":", linewidth=0.8, alpha=0.6)
ax[0, 0].axhline(Y_star, color="#F44336", linestyle=":", linewidth=0.8, alpha=0.6)
ax[0, 0].set_xlabel("Time (s)")
ax[0, 0].set_ylabel("Concentration")
ax[0, 0].set_title("(a) Time series — B = 3 (oscillatory)")
ax[0, 0].legend(loc="upper right", framealpha=0.85)
ax[0, 0].set_xlim(0, 40)

ax[0, 1].plot(X_sol[transient_mask], Y_sol[transient_mask], color="grey",
              linewidth=0.8, alpha=0.5, label="Transient")
ax[0, 1].plot(X_sol[lc_mask], Y_sol[lc_mask], color="#4CAF50",
              linewidth=1.5, label="Limit cycle")
ax[0, 1].scatter([X_star], [Y_star], color="black", s=50, zorder=5,
                 label=f"Unstable focus ({X_star:.0f}, {Y_star:.0f})")
ax[0, 1].set_xlabel("X")
ax[0, 1].set_ylabel("Y")
ax[0, 1].set_title("(b) Phase portrait — B = 3")
ax[0, 1].legend(loc="upper right", framealpha=0.85, fontsize=9)

# Bottom row: stable spiral (B=1.5)
ax[1, 0].plot(sol_s.t, X_s, color="#2196F3", label="X")
ax[1, 0].plot(sol_s.t, Y_s, color="#F44336", label="Y")
ax[1, 0].set_xlabel("Time (s)")
ax[1, 0].set_ylabel("Concentration")
ax[1, 0].set_title(f"(c) Time series — B = {B_stable} (damped)")
ax[1, 0].legend(loc="upper right", framealpha=0.85)

ax[1, 1].plot(X_s, Y_s, color="#FF9800", linewidth=1.5, label="Trajectory")
ax[1, 1].scatter([A_val], [B_stable/A_val], color="black", s=50, zorder=5,
                 label=f"Stable spiral ({A_val:.0f}, {B_stable:.1f})")
ax[1, 1].set_xlabel("X")
ax[1, 1].set_ylabel("Y")
ax[1, 1].set_title(f"(d) Phase portrait — B = {B_stable}")
ax[1, 1].legend(loc="upper right", framealpha=0.85, fontsize=9)

plt.tight_layout()
out_path3 = _here / "brusselator_combined.png"
fig3.savefig(out_path3, dpi=150)
print(f"Combined figure saved to {out_path3}")
