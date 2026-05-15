"""Example 04 — Catalytic Hairpin Assembly (CHA) miR-21 biosensor cascade.

Models a four-reaction DNA strand-displacement circuit that detects miR-21
(22-nt microRNA) via catalytic hairpin assembly.  Rate constants are derived
from NUPACK MFE calculations (codesign_h1_h2_cp.py) at 37 °C in 137 mM NaCl,
10 mM MgCl₂, using the Zhang & Winfree 2009 toehold-length → k_f lookup and
detailed balance (k_r = k_f · exp(ΔΔG / RT)).

Reactions
---------
  R1: miR21 + H1     ⇌  miR21·H1           toehold binding     (D1 = 7 nt)
  R2: miR21·H1 + H2  ⇌  H1·H2 + miR21      strand displacement (D2 = 15 nt)
  R3: H1·H2 + CP     ⇌  H1·H2·CP           capture-probe signal (T = 8 nt)
  R4: H1 + H2        ⇌  H1·H2              spontaneous leakage

CRNT result: δ = 8 − 4 − 3 = 1, weakly reversible, D1T applies →
at most one positive steady state per stoichiometric compatibility class.

Outputs
-------
  cha_bifurcation.png       — kinetic dose-response scan (signal at t = 1 h)
  cha_reaction_graph.png    — complex reaction graph coloured by linkage class
  cha_kinetics_species.png  — full species traces at four miR-21 concentrations
  cha_kinetics_analysis.png — CP-bound signal comparison and detection t½ curve
"""
import pathlib
import numpy as np
import matplotlib
matplotlib.use("Agg")   # headless backend for non-interactive use
matplotlib.rcParams.update({"font.family": "STIXGeneral", "mathtext.fontset": "stix"})
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from scipy.integrate import solve_ivp

from mantis import CRNetwork

_here = pathlib.Path(__file__).parent

# ── Reaction strings ────────────────────────────────────────────────────────
CHA_STRINGS = [
    "miR21 + H1 <-> miR21_H1",
    "miR21_H1 + H2 <-> H1H2 + miR21",
    "H1H2 + CP <-> H1H2_CP",
    "H1 + H2 <-> H1H2",
]

# Rate constants from codesign_h1_h2_cp.py (MFE, 37°C, physiological salt)
# Domain split: D1=7 nt toehold (3' miR-21), D2=15 nt stem (5' miR-21), T=8 nt tail
# k_f from Zhang & Winfree 2009 (≥7 nt → saturated 1e6); k_r = k_f·exp(ΔΔG/RT)
# ΔΔG: R1=-8.84, R2=-4.45, R3=-10.93 kcal/mol; leakage Ea=18.29 kcal/mol (H1 stem)
CHA_RATES = {
    "miR21 + H1 -> miR21_H1":          1.0e6,      # k_on  (M⁻¹s⁻¹), 7-nt D1 toehold
    "miR21_H1 -> miR21 + H1":          5.8925e-1,  # k_off (s⁻¹),    ΔΔG=-8.84 kcal/mol
    "miR21_H1 + H2 -> H1H2 + miR21":   1.0e6,      # k_on  (M⁻¹s⁻¹), 15-nt D2
    "H1H2 + miR21 -> miR21_H1 + H2":   7.3116e2,   # k_rev (M⁻¹s⁻¹), ΔΔG=-4.45 kcal/mol
    "H1H2 + CP -> H1H2_CP":            1.0e6,      # k_on  (M⁻¹s⁻¹), 8-nt tail T
    "H1H2_CP -> H1H2 + CP":            1.9836e-2,  # k_off (s⁻¹),    ΔΔG=-10.93 kcal/mol
    "H1 + H2 -> H1H2":                 1.2904e-7,  # k_leak (M⁻¹s⁻¹), Ea=18.29 kcal/mol
    "H1H2 -> H1 + H2":                 5.5596e-17, # essentially irreversible
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

# ── 3. Symbolic ODEs ─────────────────────────────────────────────────────────
print()
print("Symbolic ODEs (with numeric rates):")
odes = rn.odes(numeric_rates=True)
for sp in sorted(odes):
    print(f"  d{sp}/dt = {odes[sp]}")

# ── 4. Steady-state analysis ─────────────────────────────────────────────────
print()
print("Solving for steady state ([miR21]₀ = 10 nM) ...")
ss_list = rn.steady_states(CHA_IC, n_attempts=30, seed=0)
# D1T guarantees ≤1 SS per stoichiometry class; extras are numerical duplicates.
n_stable = sum(1 for s in ss_list if s.is_stable)
if len(ss_list) > 1:
    print(f"  Found {len(ss_list)} candidate(s) ({n_stable} stable) — "
          f"D1T guarantees ≤1; reporting best by residual.")
else:
    print(f"  Found {len(ss_list)} steady state(s)")
ss = ss_list[0]  # sorted by residual ascending; always the best solution

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

# ── 7. Kinetic scan over miR21 concentration (signal at t = 1 h) ─────────────
# True thermodynamic SS is dominated by the slow leakage pathway (τ ~ 4×10⁶ s)
# and is nearly miR21-independent.  The biologically relevant readout is the
# kinetic signal at a fixed incubation time (here 1 h = 3600 s).
print()
print("Running kinetic scan over miR21 initial concentration (signal at t = 1 h) ...")
print("(scanning initial [miR21] from 0.1 nM to 100 nM)")

mir21_values = np.logspace(-10, -7, 20)  # 0.1 nM → 100 nM
signal_values = []

for mir21_conc in mir21_values:
    ic_scan = dict(CHA_IC)
    ic_scan["miR21"] = float(mir21_conc)
    ic_scan["miR21_H1"] = 0.0
    sim = rn.simulate(ic_scan, t_span=(0, 3600))
    signal_values.append(sim.final().get("H1H2_CP", 0.0))

fig, ax = plt.subplots(figsize=(7, 4))
ax.loglog(mir21_values * 1e9, np.maximum(signal_values, 1e-15) * 1e9, "o-", color="#2196F3")
ax.set_xlabel("[miR21]₀ (nM)")
ax.set_ylabel("[H1H2_CP]* (nM)")
ax.set_title("CHA Signal vs. miR21 Concentration (t = 1 h)")
ax.grid(True, which="both", alpha=0.3)
plt.tight_layout()
out_path = _here / "cha_bifurcation.png"
fig.savefig(out_path, dpi=150)
print(f"  Bifurcation plot saved to {out_path}")

# ── 8. Reaction graph ────────────────────────────────────────────────────────
fig2, ax2 = plt.subplots(figsize=(11, 7))
rn.draw(ax=ax2)
plt.tight_layout()
out_path2 = _here / "cha_reaction_graph.png"
fig2.savefig(out_path2, dpi=160)
print(f"  Reaction graph saved to {out_path2}")

# ── 9. ODE kinetic simulation ────────────────────────────────────────────────
# Uses scipy solve_ivp directly with rates from CHA_RATES — no NUPACK calls.
print()
print("Running ODE kinetic simulations …")

K1_F = CHA_RATES["miR21 + H1 -> miR21_H1"]
K1_R = CHA_RATES["miR21_H1 -> miR21 + H1"]
K2_F = CHA_RATES["miR21_H1 + H2 -> H1H2 + miR21"]
K2_R = CHA_RATES["H1H2 + miR21 -> miR21_H1 + H2"]
K3_F = CHA_RATES["H1H2 + CP -> H1H2_CP"]
K3_R = CHA_RATES["H1H2_CP -> H1H2 + CP"]
K4_F = CHA_RATES["H1 + H2 -> H1H2"]
K4_R = CHA_RATES["H1H2 -> H1 + H2"]

C_H1 = CHA_IC["H1"]
C_H2 = CHA_IC["H2"]
C_CP = CHA_IC["CP"]

T_END  = 7200                          # 2 h in seconds
T_EVAL = np.linspace(0, T_END, 10_000)


def _cha_odes(t, y):
    # species order: H1, H2, CP, miR21, miR21_H1, H1H2, H1H2_CP
    H1, H2, CP, miR21, miR21_H1, H1H2, H1H2_CP = y
    r1_f = K1_F * miR21    * H1
    r1_r = K1_R * miR21_H1
    r2_f = K2_F * miR21_H1 * H2
    r2_r = K2_R * H1H2     * miR21
    r3_f = K3_F * H1H2     * CP
    r3_r = K3_R * H1H2_CP
    r4_f = K4_F * H1       * H2
    r4_r = K4_R * H1H2
    return [
        -r1_f + r1_r - r4_f + r4_r,              # dH1
        -r2_f + r2_r - r4_f + r4_r,              # dH2
        -r3_f + r3_r,                             # dCP
        -r1_f + r1_r + r2_f - r2_r,              # dmiR21
         r1_f - r1_r - r2_f + r2_r,              # dmiR21_H1
         r2_f - r2_r - r3_f + r3_r + r4_f - r4_r,  # dH1H2
         r3_f - r3_r,                             # dH1H2_CP
    ]


def _run_ode(c_mir21):
    y0 = [C_H1, C_H2, C_CP, c_mir21, 0.0, 0.0, 0.0]
    return solve_ivp(_cha_odes, (0, T_END), y0,
                     t_eval=T_EVAL, method="Radau", rtol=1e-10, atol=1e-16)


scenarios = [
    (10e-9,   "10 nM miR-21",             "#2196F3"),
    ( 1e-9,   " 1 nM miR-21",             "#4CAF50"),
    (100e-12, "100 pM miR-21",            "#FF9800"),
    ( 0,      "No trigger (leakage only)", "#F44336"),
]

sim_results = [(c, lbl, col, _run_ode(c)) for c, lbl, col in scenarios]
print("  Simulations complete.")

# ── Figure: full species traces (2×2 grid) ────────────────────────────────────
t_min = T_EVAL / 60

SPECIES_STYLE = {
    "H₁ (free)":          ("#4C72B0", "-",   1.6),
    "H₂ (free)":          ("#DD8452", "-",   1.6),
    "miR-21 (free)":      ("#55A868", ":",   1.4),
    "miR-21·H₁":          ("#C44E52", "-.",  1.4),
    "H₁·H₂ dimer":        ("#8172B3", "-",   2.2),
    "H₁·H₂·CP (signal)":  ("#1a1a1a", "--",  2.5),
}

fig_sp, axes_sp = plt.subplots(2, 2, figsize=(14, 9), sharex=True)
fig_sp.suptitle("CHA Cascade — Full Kinetics  (37 °C, [H₁]=[H₂]=[CP]=100 nM)",
                fontsize=13, y=1.01)

for ax_sp, (c_mir21, label, _, sol) in zip(axes_sp.flat, sim_results):
    H1, H2, CP, miR21, miR21_H1, H1H2, H1H2_CP = sol.y
    nM = 1e9
    traces = [
        ("H₁ (free)",         H1),
        ("H₂ (free)",         H2),
        ("miR-21 (free)",     miR21),
        ("miR-21·H₁",         miR21_H1),
        ("H₁·H₂ dimer",       H1H2),
        ("H₁·H₂·CP (signal)", H1H2_CP),
    ]
    for name, conc in traces:
        color, ls, lw = SPECIES_STYLE[name]
        ax_sp.plot(t_min, conc * nM, label=name, color=color, ls=ls, lw=lw)
    ax_sp.set_title(label, fontsize=11)
    ax_sp.set_xlabel("Time (min)")
    ax_sp.set_ylabel("Concentration (nM)")
    ax_sp.legend(fontsize=7.5, loc="center right")
    ax_sp.set_xlim(0, 120)
    ax_sp.set_ylim(-2, 112)
    ax_sp.grid(True, alpha=0.25)

plt.tight_layout()
out_sp = _here / "cha_kinetics_species.png"
fig_sp.savefig(out_sp, dpi=150, bbox_inches="tight")
print(f"  Species traces saved to {out_sp}")
plt.close(fig_sp)

# ── Figure: signal comparison + dose-response ─────────────────────────────────
fig_an = plt.figure(figsize=(14, 5))
gs_an  = gridspec.GridSpec(1, 2, figure=fig_an, wspace=0.38)
ax_sig, ax_dose = (fig_an.add_subplot(gs_an[i]) for i in range(2))

for c_mir21, label, color, sol in sim_results:
    signal_pct = sol.y[6] / C_CP * 100
    ls = "--" if c_mir21 == 0 else "-"
    ax_sig.plot(t_min, signal_pct, label=label, color=color, lw=2, ls=ls)

ax_sig.set_xlabel("Time (min)")
ax_sig.set_ylabel("CP bound (%)")
ax_sig.set_title("Detection Signal vs Time")
ax_sig.legend(fontsize=8.5)
ax_sig.grid(True, alpha=0.25)
ax_sig.set_xlim(0, 120)
ax_sig.set_ylim(-1, 102)

print("  Computing dose–response curve …")
concs = np.logspace(-11.5, -7, 40)   # ~30 pM → 100 nM
t_half_arr = []
for c in concs:
    sol_d = _run_ode(c)
    sig   = sol_d.y[6]
    max_s = sig[-1]
    if max_s < 1e-15:
        t_half_arr.append(np.nan)
        continue
    idx = np.searchsorted(sig, max_s * 0.5)
    t_half_arr.append(T_EVAL[min(idx, len(T_EVAL) - 1)] / 60)

t_half_arr = np.array(t_half_arr, dtype=float)
valid = ~np.isnan(t_half_arr)
ax_dose.loglog(concs[valid] * 1e9, t_half_arr[valid], "o-", color="#2196F3", lw=2, ms=5)
ax_dose.set_xlabel("[miR-21] (nM)")
ax_dose.set_ylabel("Time to 50% max signal (min)")
ax_dose.set_title("Detection Speed vs [miR-21]")
ax_dose.grid(True, which="both", alpha=0.25)
ax_dose.invert_yaxis()

fig_an.suptitle("CHA Kinetics — Analysis", fontsize=13)
plt.tight_layout()
out_an = _here / "cha_kinetics_analysis.png"
fig_an.savefig(out_an, dpi=150, bbox_inches="tight")
print(f"  Analysis panels saved to {out_an}")
plt.close(fig_an)

print()
print("Done.")
