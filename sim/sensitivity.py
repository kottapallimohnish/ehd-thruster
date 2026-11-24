"""
Uncertainty of the headline numbers, and how robust the onset claim is.

1. Onset claim (smooth wires). Without space charge, Laplace's equation is linear, so every
   surface field scales with V: each wire's smooth-wire onset voltage is
   V_on,i = 11 kV * E_onset / E_max,i(11 kV). Below onset there are no ions, so onset is set
   exactly by the Laplace field (no space-charge correction applies). Checked at the nominal gap and Peek factor m = 1, and in a pessimistic
   corner (gap 9.86 mm, m = 0.9), because the array field is dominated by the 2*pi*d/s term.

2. Operating point (current, efficiency). The measured jet fixes the thrust delivered;
   T = I * d / mu then gives
        I   = T_meas(v) / r^2 * mu / d,      T_meas = 1/2 rho A v^2
        eta = P_jet / (V * I),               P_jet  = 1/2 T_meas v  (ideal-disk model of the jet)
   Inputs and their ranges (sampled uniformly and independently, Monte Carlo):
        v      anemometer reading 1.21 m/s, instrument accuracy +-(0.5 m/s + 5 %)
        mu     1.8 - 2.2e-4 m^2/Vs
        d      9.86 - 12.36 mm (thread-count uncertainty)
        r_cfd  jet survival ratio from the 2D CFD: both meshes, each +-2 standard errors of
               its time mean (batch means), times the force-level factor
        f_3d   2D -> 3D correction on r, 0.7 - 1.3 (not computed; assumed +-30 %; inviscid
               actuator theory at the anemometer distance suggests ~1.1, inside this range)
        V      9 - 14 kV peak (efficiency and space-charge share only)
   Reported: median and 5-95 % interval, and the interval from each input alone. The
   anemometer accuracy is the largest single term; the assumed 2D -> 3D factor is second
   and nearly as large.
"""
import csv
import json
from pathlib import Path

import numpy as np
import params as P

RES = Path(__file__).parent / "results"

# ------------------------------------------------------------------ 1. onset claim
rows = list(csv.DictReader(open(RES / "stage1_wire_fields.csv")))
Emax = np.array([float(r["E_max_MV_m"]) for r in rows]) * 1e6
V0 = abs(P.V_WIRE)
V_on = V0 * P.E_ONSET_SMOOTH / Emax
counts = {f"{v/1e3:.0f} kV": int((V_on <= v).sum()) for v in (9e3, 10e3, 11e3, 12e3, 13e3, 14e3)}
# pessimistic corner: shorter gap raises the array field (dominant 2*pi*d/s term), m = 0.9 lowers onset
s, r_w = P.WIRE_PITCH, P.R_WIRE
geo = lambda d: 2 * np.pi * d / s + np.log(s / (2 * np.pi * r_w))
gap_factor = geo(P.GAP) / geo(9.86)
E_on_m09 = 3.1e6 * 0.9 * 0.944 * (1 + 0.0308 / np.sqrt(0.944 * r_w * 1e-3))
V_on_corner = V0 * E_on_m09 / (Emax * gap_factor)
inner = slice(2, 9)

# ------------------------------------------------------------------ 2. operating point
fine = json.loads((RES / "stage3_summary.json").read_text())
coarse = json.loads((RES / "stage3_summary_coarse.json").read_text())
hi = RES / "stage3_summary_coarse_hiforce.json"
v_ideal = fine["ideal_2d_actuator_velocity_mps"]
r_fine = fine["survival_ratio_r"]
r_coarse = coarse["survival_ratio_r"]
se_f = fine["anemometer_mean_batch_se_mps"] / v_ideal
se_c = coarse["anemometer_mean_batch_se_mps"] / v_ideal
force_factor = 1.0
if hi.exists():                      # same coarse mesh, 3.2x force: effect of the higher force on r
    h = json.loads(hi.read_text())
    force_factor = (h["cfd_velocity_at_anemometer_mps"] / h["ideal_2d_actuator_velocity_mps"]) / r_coarse
r_nom = r_fine * force_factor
r_lo = min(r_fine - 2 * se_f, r_coarse - 2 * se_c) * force_factor
r_hi = max(r_fine + 2 * se_f, r_coarse + 2 * se_c) * force_factor

d_axial = json.loads((RES / "stage2_high_current_check.json").read_text())["d_axial_mm"] * 1e-3
A = P.A_FLOW_MM2 * 1e-6
rng = np.random.default_rng(1)
N = 200_000
v_acc = 0.5 + 0.05 * P.V_JET
v = rng.uniform(P.V_JET - v_acc, P.V_JET + v_acc, N)
mu = rng.uniform(1.8e-4, 2.2e-4, N)
gap = rng.uniform(9.86, 12.36, N)                      # mm
dd = d_axial * gap / P.GAP
r = rng.uniform(r_lo, r_hi, N) * rng.uniform(0.7, 1.3, N)
V = rng.uniform(9e3, 14e3, N)
T_meas = 0.5 * P.RHO_AIR * A * v**2
I = T_meas / r**2 * mu / dd
eta = 0.5 * T_meas * v / (V * I)
I_MG = 9 / 8 * P.EPS0 * mu * V**2 / (gap * 1e-3) ** 3 * A
share = I / I_MG
pct = lambda x, k=1.0: [round(float(np.percentile(x, q)) * k, 3) for q in (5, 50, 95)]

# one input at a time (others at nominal): 5-95 % multiplicative range of the current
def alone(name):
    vv = v if name == "anemometer" else np.full(N, P.V_JET)
    mm = mu if name == "mobility" else np.full(N, 2e-4)
    gg = dd if name == "gap" else np.full(N, d_axial)
    if name == "cfd_r":
        rr = rng.uniform(r_lo, r_hi, N)
    elif name == "2d_to_3d":
        rr = r_nom * rng.uniform(0.7, 1.3, N)
    else:
        rr = np.full(N, r_nom)
    Ii = 0.5 * P.RHO_AIR * A * vv**2 / rr**2 * mm / gg
    I0 = 0.5 * P.RHO_AIR * A * P.V_JET**2 / r_nom**2 * 2e-4 / d_axial
    return [round(float(np.percentile(Ii / I0, q)), 2) for q in (5, 95)]

I_nom = 0.5 * P.RHO_AIR * A * P.V_JET**2 / r_nom**2 * 2e-4 / d_axial
eta_nom = 0.5 * (0.5 * P.RHO_AIR * A * P.V_JET**2) * P.V_JET / (V0 * I_nom)
out = dict(
    onset=dict(smooth_wire_onset_kV_per_wire=[round(x / 1e3, 2) for x in V_on],
               wires_above_onset_at_voltage=counts,
               inner7_min_onset_kV_nominal=round(float(V_on[inner].min()) / 1e3, 2),
               inner7_min_onset_kV_pessimistic=round(float(V_on_corner[inner].min()) / 1e3, 2),
               pessimistic_corner="gap 9.86 mm and Peek factor m = 0.9",
               note="Laplace fields: below onset there is no space charge, so onset is set by these fields"),
    survival_ratio_r=dict(nominal=round(r_nom, 3), cfd_range=[round(r_lo, 3), round(r_hi, 3)],
                          force_effect_factor=round(force_factor, 3), plus_2d_to_3d="x 0.7-1.3 assumed"),
    nominal=dict(current_uA=round(I_nom * 1e6, 1), efficiency_pct_at_11kV=round(100 * eta_nom, 3),
                 power_W_at_11kV=round(V0 * I_nom, 2), space_charge_share_at_11kV=round(
                     float(I_nom / (9 / 8 * P.EPS0 * 2e-4 * V0**2 / (P.GAP * 1e-3) ** 3 * A)), 3)),
    current_uA_5_50_95=pct(I, 1e6), efficiency_pct_5_50_95=pct(eta, 100),
    power_W_5_50_95=pct(V * I), space_charge_share_5_50_95=pct(share),
    current_uA_from_anemometer_alone=[round(I_nom * ((P.V_JET + sgn * v_acc) / P.V_JET) ** 2 * 1e6, 1) for sgn in (-1, 1)],
    smooth_to_nicked_current_ratio=[round((P.V_JET / 0.2) ** 2), round((P.V_JET / 0.1) ** 2)],
    current_factor_5_95_each_input_alone={k: alone(k) for k in ("anemometer", "2d_to_3d", "cfd_r", "mobility", "gap")},
    thrust_per_power_mN_per_W=dict(
        ion_force_per_watt_nominal=round(1e3 * d_axial / (2e-4 * V0), 2),
        jet_thrust_per_watt_nominal=round(1e3 * (0.5 * P.RHO_AIR * A * P.V_JET**2) / (V0 * I_nom), 2),
        jet_thrust_per_watt_5_50_95=pct(T_meas / (V * I), 1e3)),
)
(RES / "sensitivity.json").write_text(json.dumps(out, indent=2))
print(json.dumps(out, indent=2))
