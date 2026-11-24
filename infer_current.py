"""
FIRST-PASS ESTIMATE (superseded). Assumes an ideal jet with no losses on the collector.
Stage 3 (sim/stage3_*) showed the tubes absorb ~60 % of the ion force, which raises the
current to ~130 uA and lowers the efficiency to ~0.1 %. Final values: RESULTS.md.

No current was measured. Invert the EHD model: from the measured jet velocity,
infer the corona current, thrust and efficiency, then check plausibility.
See theory.md sections 3-6.
"""
import json
import numpy as np
from pathlib import Path

EPS0 = 8.854e-12
g = json.loads((Path(__file__).parent / "data/geometry.json").read_text())
rho = g["air_density_moist_kg_m3"]
d = g["gap_mm"] * 1e-3
A_full = g["flow_area_mm2"] * 1e-6
A_open = g["open_flow_area_mm2"] * 1e-6
L_em = g["emitter_length_model_mm"] * 1e-3
V_nom = g["supply"]["spark_test"]["peak_voltage_kV_estimate"] * 1e3


def infer(v, mu, A, V):
    I = rho * mu * A * v**2 / (2 * d)      # from v_w = sqrt(2 I d / (rho mu A))
    T = I * d / mu                          # sec. 3
    P = V * I                               # corona power (peak V: upper bound)
    eta = 0.5 * T * v / P                   # jet kinetic power / corona power
    return I, T, P, eta


v = 1.21  # time-averaged steady reading, WA_23.06.38, t=10-48 s
print(f"rho={rho} kg/m3, d={d*1e3:.2f} mm, A_full={A_full*1e6:.0f} mm2, A_open={A_open*1e6:.0f} mm2\n")
print(f"Nominal (v={v} m/s, mu=2.0e-4, V=11 kV, full-area actuator disk):")
I, T, P, eta = infer(v, 2.0e-4, A_full, V_nom)
print(f"  current I     = {I*1e6:6.1f} uA   ({I/L_em*1e6:.0f} uA per metre of emitter)")
print(f"  thrust T      = {T*1e3:6.2f} mN")
print(f"  corona power  = {P:6.2f} W   (battery side ~7-18 W)")
print(f"  EHD efficiency= {eta*100:6.2f} %")
print(f"  T/P           = {T/P*1e3:6.2f} mN/W  (theory d/(mu V) = {d/(2e-4*V_nom)*1e3:.2f})")
print(f"  ion drift speed mu*E = {2e-4*V_nom/d:.0f} m/s vs air {v} m/s")

J_mg = 9 / 8 * EPS0 * 2.0e-4 * V_nom**2 / d**3
print(f"\nPlausibility: Mott-Gurney space-charge limit over the opening = {J_mg*A_full*1e6:.0f} uA"
      f" -> inferred current is {I/(J_mg*A_full)*100:.0f} % of the limit")

print("\nUncertainty bands:")
for lab, vv in [("v low (instrument)", v - (0.5 + 0.05 * v)), ("v nom", v),
                ("v high (instrument)", v + (0.5 + 0.05 * v))]:
    for A, al in [(A_full, "full"), (A_open, "open")]:
        for mu in (1.8e-4, 2.2e-4):
            I, T, P, eta = infer(vv, mu, A, V_nom)
            print(f"  {lab:20s} A={al:4s} mu={mu:.1e}: I={I*1e6:6.1f} uA  T={T*1e3:5.2f} mN  eta={eta*100:.2f} %")
