"""Example 02 — Michaelis-Menten enzyme kinetics.

System:   E + S <-> ES -> E + P
Demonstrates:
- CRNT analysis (δ=0, not weakly reversible — DZT stability applies only for WR)
- Enzyme conservation law: E + ES = E_total
- Quasi-steady-state approximation validation via numerical SS
- Eigenvalue analysis showing stable steady state
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

from pycrn import CRNetwork
from pycrn.analysis import build_ode_function

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

# ── Figure: mass-action time course + QSS-vs-mass-action dynamic comparison
import pathlib
_here = pathlib.Path(__file__).parent

species_order = list(rn.species)
S_idx, E_idx, ES_idx, P_idx = (species_order.index(s) for s in ("S", "E", "ES", "P"))

f_ode = build_ode_function(rn._reactions, species_order, rn._rates)
y0 = np.zeros(len(species_order))
y0[S_idx], y0[E_idx] = S0, E0
t_span = (0.0, 1.5)
t_eval = np.linspace(*t_span, 6000)
sol = solve_ivp(f_ode, t_span, y0, method="Radau", t_eval=t_eval,
                rtol=1e-10, atol=1e-14)

S_t  = sol.y[S_idx]
E_t  = sol.y[E_idx]
ES_t = sol.y[ES_idx]
P_t  = sol.y[P_idx]
ES_qss_t = E0 * S_t / (Km + S_t)

fig, ax = plt.subplots(1, 2, figsize=(11, 4))

# Panel (a): full time course with twin axis for ES (concentration scale separation)
ax_a = ax[0]
ax_a2 = ax_a.twinx()
ax_a.plot(sol.t, S_t * 1e6, label="S (left)",  color="#1976D2")
ax_a.plot(sol.t, P_t * 1e6, label="P (left)",  color="#388E3C")
ax_a2.plot(sol.t, ES_t * 1e6, label="ES (right)", color="#F57C00")
ax_a2.plot(sol.t, E_t  * 1e6, label="E (right)",  color="#7B1FA2", linestyle="--")
ax_a.set_xlabel("Time (s)")
ax_a.set_ylabel("S, P  (µM)")
ax_a2.set_ylabel("E, ES  (µM)")
ax_a.set_title(f"(a) Mass-action time course  ($E_0$ = {E0*1e6:.0f} µM, $S_0$ = {S0*1e6:.0f} µM)")
lines1, labels1 = ax_a.get_legend_handles_labels()
lines2, labels2 = ax_a2.get_legend_handles_labels()
ax_a.legend(lines1 + lines2, labels1 + labels2, loc="center right", framealpha=0.85)

# Panel (b): dynamic QSS validation — full mass-action ES(t) vs QSS prediction
ax_b = ax[1]
ax_b.plot(sol.t, ES_t * 1e6, label="ES(t)  mass-action", color="#1976D2", linewidth=2)
ax_b.plot(sol.t, ES_qss_t * 1e6,
          label=r"$E_0\,[S](t)\,/\,(K_M+[S](t))$  (QSS)",
          color="#D32F2F", linestyle="--", linewidth=2)
ax_b.set_xlabel("Time (s)")
ax_b.set_ylabel("[ES]  (µM)")
ax_b.set_title(r"(b) QSS approximation tracks full mass-action $[ES](t)$")
ax_b.legend(framealpha=0.85)
ax_b.grid(True, alpha=0.3)

# annotate the worst-case relative error after the brief transient (t > 5 ms)
mask = sol.t > 5e-3
rel_err = np.abs(ES_t[mask] - ES_qss_t[mask]) / np.maximum(ES_t[mask], 1e-15)
ax_b.text(0.55, 0.15,
          f"max |ES − QSS| / ES = {rel_err.max()*100:.2f}%  (t > 5 ms)",
          transform=ax_b.transAxes, fontsize=9,
          bbox=dict(facecolor="white", alpha=0.85, edgecolor="grey"))

plt.tight_layout()
out_path = _here / "michaelis_menten.png"
fig.savefig(out_path, dpi=150)
print(f"Figure saved to {out_path}")
print(f"Max QSS relative error (after transient): {rel_err.max()*100:.3f}%")
