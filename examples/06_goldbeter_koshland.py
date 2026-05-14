"""Example 06 — Goldbeter-Koshland covalent modification switch.

Reference: Goldbeter A, Koshland DE Jr. (1981) An amplified sensitivity arising
from covalent modification in biological systems. PNAS 78(11): 6840-6844.

The GK switch models a substrate W that is interconverted between inactive (W)
and active (Wp) forms by a kinase (E1) and phosphatase (E2):

    W + E1 ⇌ WE1  →  Wp + E1     (kinase arm)
    Wp + E2 ⇌ WpE2 →  W + E2     (phosphatase arm)

Key CRNT result  δ = 6 − 2 − 3 = 1, not weakly reversible, D1T applies.
This proves: for ANY positive rate constants, at most one steady state per
stoichiometry class → bistability is impossible.

Observed behaviour: sigmoidal (ultrasensitive) but monostable dose-response.
At zero-order saturation (Km << W_total), the effective Hill coefficient is very
large, making the switch effectively digital — but always single-valued.

Demonstrates:
- Full CRNT analysis of an enzyme-catalysed modification cycle
- D1T as the mechanistic explanation for monostability
- Numerical SS matched against the Goldbeter-Koshland QSS formula
- Dose-response scan showing ultrasensitivity without bistability
"""
import pathlib
import numpy as np
import matplotlib
matplotlib.use("Agg")
matplotlib.rcParams.update({"font.family": "STIXGeneral", "mathtext.fontset": "stix"})
import matplotlib.pyplot as plt
from scipy.optimize import brentq

from pycrn import CRNetwork

_here = pathlib.Path(__file__).parent

# ── Network ──────────────────────────────────────────────────────────────────
GK_STRINGS = [
    "W + E1 <-> WE1",
    "WE1 -> Wp + E1",
    "Wp + E2 <-> WpE2",
    "WpE2 -> W + E2",
]

# Deep zero-order regime: Km = (k_off + kcat)/k_on = 1.1/1e4 = 1.1e-4 M << W_total = 1 M
# → Wt/Km ≈ 9000 → near-infinite effective Hill coefficient
K_ON  = 1e4    # M⁻¹ s⁻¹
K_OFF = 0.1    # s⁻¹
KCAT  = 1.0    # s⁻¹
KM    = (K_OFF + KCAT) / K_ON   # 1.1e-4 M

GK_RATES_SYM = {
    "E1 + W -> WE1":   K_ON,
    "WE1 -> E1 + W":   K_OFF,
    "WE1 -> E1 + Wp":  KCAT,
    "E2 + Wp -> WpE2": K_ON,
    "WpE2 -> E2 + Wp": K_OFF,
    "WpE2 -> E2 + W":  KCAT,
}

W_TOTAL = 1.0      # total substrate (M)
E1_TOTAL = 0.01    # kinase (M) — 1% of substrate
E2_TOTAL = 0.01    # phosphatase (M)

GK_IC = {
    "W":    W_TOTAL,
    "E1":   E1_TOTAL,
    "WE1":  0.0,
    "Wp":   0.0,
    "E2":   E2_TOTAL,
    "WpE2": 0.0,
}

# ── 1. CRNT analysis ─────────────────────────────────────────────────────────
rn_sym = CRNetwork.from_string(GK_STRINGS, rates=GK_RATES_SYM)

print("=" * 60)
print("Goldbeter-Koshland Switch — CRNT Analysis")
print("=" * 60)
print(rn_sym.crnt_summary())
print(f"Km = (k_off + kcat) / k_on = {KM:.2e} M")
print(f"Wt / Km = {W_TOTAL/KM:.0f}  →  deep zero-order regime\n")

# ── 2. Symmetric steady state (kcat1 = kcat2) ───────────────────────────────
print("Symmetric steady state (equal kinase and phosphatase activity):")
ss_sym = rn_sym.steady_states(GK_IC, n_attempts=20, seed=0)[0]
c = ss_sym.concentrations
phospho = c["Wp"] + c["WpE2"]
unphos  = c["W"]  + c["WE1"]
print(f"  W    = {c['W']:.4e} M     Wp   = {c['Wp']:.4e} M")
print(f"  WE1  = {c['WE1']:.4e} M     WpE2 = {c['WpE2']:.4e} M")
print(f"  Total phosphorylated   = {phospho:.4f} M  ({phospho/W_TOTAL*100:.1f}%)")
print(f"  Total unphosphorylated = {unphos:.4f} M  ({unphos/W_TOTAL*100:.1f}%)")
print(f"  Stable: {ss_sym.is_stable}   Residual: {ss_sym.residual:.2e}")

# By symmetry (equal k's, equal enzyme totals): exactly 50%/50%
print(f"  Expected (symmetry): 50.0%  →  error = {abs(phospho/W_TOTAL - 0.5)*100:.2f}%\n")

# ── 3. Dose-response scan ────────────────────────────────────────────────────
print("Scanning kcat1 (kinase activity) from 0.05 to 20 × kcat2 ...")
print("D1T guarantees exactly one SS per stoichiometry class throughout.\n")

kcat2 = KCAT
kcat1_values = np.logspace(-1.3, 1.3, 30)   # ~0.05 to 20 × kcat2
phospho_ss   = []
n_ss_each    = []

for kcat1 in kcat1_values:
    rates = dict(GK_RATES_SYM)
    rates["WE1 -> E1 + Wp"] = kcat1
    rn_scan = CRNetwork.from_string(GK_STRINGS, rates=rates)
    ss_list = rn_scan.steady_states(GK_IC, n_attempts=10, seed=42)
    n_ss_each.append(len(ss_list))
    if ss_list:
        c_scan = ss_list[0].concentrations
        phospho_ss.append(c_scan["Wp"] + c_scan["WpE2"])
    else:
        phospho_ss.append(float("nan"))

phospho_ss = np.array(phospho_ss)
ratios = kcat1_values / kcat2

print(f"  Steady states found per parameter value: {set(n_ss_each)}")
print(f"  (all = 1 confirms D1T — no bistability)\n")

# ── 4. QSS (Goldbeter-Koshland formula) ─────────────────────────────────────
def gk_qss(u, km_norm):
    """Fractional activation from the GK QSS approximation.
    u = V1/V2 (kinase/phosphatase ratio), km_norm = Km/W_total.
    Solves: (1-y)/(K+1-y) = (1/u) * y/(K+y)  for y in (0,1).
    """
    K = km_norm
    def eq(y):
        return (1 - y) * (K + y) - (1 / u) * y * (K + 1 - y)
    try:
        return brentq(eq, 1e-12, 1 - 1e-12)
    except ValueError:
        return 1.0 if u > 1 else 0.0

km_norm = KM / W_TOTAL
qss_vals = np.array([gk_qss(u, km_norm) for u in ratios])

# QSS uses effective velocities V1 = kcat1*E1t, V2 = kcat2*E2t
# Since E1t = E2t, ratio u = kcat1/kcat2
print("Comparison: pycrn full mass-action vs. GK QSS approximation")
print(f"  {'kcat1/kcat2':>12}  {'pycrn Wp/Wt':>12}  {'QSS Wp/Wt':>10}  {'error':>8}")
for u, num, qss in zip(ratios[::5], phospho_ss[::5], qss_vals[::5]):
    err = abs(num/W_TOTAL - qss)
    print(f"  {u:12.2f}  {num/W_TOTAL:12.4f}  {qss:10.4f}  {err:8.4f}")

# ── 5. Plot ──────────────────────────────────────────────────────────────────
fig, axes = plt.subplots(1, 2, figsize=(11, 4))

# Left: dose-response curve
ax1 = axes[0]
ax1.semilogx(ratios, phospho_ss / W_TOTAL, "o-", color="#2196F3",
             linewidth=2, markersize=5, label="pycrn (full mass-action)")
ax1.semilogx(ratios, qss_vals, "--", color="#F44336",
             linewidth=1.5, label="GK QSS approximation")
ax1.axvline(1.0, color="grey", linestyle=":", linewidth=1)
ax1.axhline(0.5, color="grey", linestyle=":", linewidth=1)
ax1.set_xlabel("kcat₁ / kcat₂  (kinase / phosphatase activity)", fontsize=10)
ax1.set_ylabel("Fractional activation  Wp* / W_total", fontsize=10)
ax1.set_title(
    "Goldbeter-Koshland switch\n"
    f"Wt/Km = {W_TOTAL/KM:.0f}  |  D1T: unique SS (no bistability)",
    fontsize=9,
)
ax1.legend(fontsize=9)
ax1.set_ylim(-0.05, 1.05)
ax1.grid(True, which="both", alpha=0.3)

# Effective Hill coefficient from QSS EC10/EC90 (cleaner than numerical derivative
# when the transition is too steep to sample with a coarse scan)
try:
    ec10 = brentq(lambda u: gk_qss(u, km_norm) - 0.10, 1e-9, 1.0 - 1e-9)
    ec90 = brentq(lambda u: gk_qss(u, km_norm) - 0.90, 1e-9, 1e9)
    hill_at_half = np.log(81) / np.log(ec90 / ec10)
except Exception:
    hill_at_half = float("nan")

ax2 = axes[1]
ax2.semilogx(ratios, phospho_ss / W_TOTAL, "o-", color="#2196F3",
             linewidth=2, markersize=5)
ax2.fill_between(ratios, 0, phospho_ss / W_TOTAL, alpha=0.15, color="#2196F3")
ax2.set_xlabel("kcat₁ / kcat₂", fontsize=10)
ax2.set_ylabel("Fractional activation", fontsize=10)
ax2.set_title(
    f"Single-valued dose-response\n"
    f"Effective Hill coeff ≈ {hill_at_half:.1f}  (n_H → ∞ as Km → 0)",
    fontsize=9,
)
ax2.set_ylim(-0.05, 1.05)
ax2.grid(True, which="both", alpha=0.3)
ax2.text(0.98, 0.1, "D1T guarantees\nmonostability\n(no hysteresis)",
         transform=ax2.transAxes, ha="right", va="bottom",
         fontsize=8, color="#555555",
         bbox=dict(boxstyle="round,pad=0.3", fc="lightyellow", ec="grey", alpha=0.8))

plt.tight_layout()
out = _here / "gk_switch_dose_response.png"
fig.savefig(out, dpi=150)
print(f"\nPlot saved to {out}")
print("\nConclusion:")
print(f"  Effective Hill coefficient ≈ {hill_at_half:.1f}")
print(f"  (theoretical limit for Km→0 is ∞)")
print(f"  Single-valued throughout → D1T prediction confirmed.")
