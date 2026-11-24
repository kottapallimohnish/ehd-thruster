"""
Validation: one wire above an infinite grounded plate, where the exact answer is known.

Physics: a cylinder of radius r whose centre is a distance h from a grounded plane has
capacitance per unit length C = 2*pi*eps0 / acosh(h/r). Its charge q = C*V, so the
circumference-averaged surface field is

    E_mean = q / (2*pi*eps0*r) = V / ( r * acosh(h/r) )   (~ V / (r ln(2h/r)) for h >> r)

The solver must reproduce this within 2 % before the thruster model is trusted.
The same wire-surface element size is then used in the main model.
"""
import json
from pathlib import Path

import numpy as np

import params as P
from fem_tools import build_mesh, classify_boundaries, solve_laplace, surface_field

OUT = Path(__file__).parent / "results"
OUT.mkdir(exist_ok=True)

h = P.GAP              # wire centre to plate, mm (same gap as the thruster)
r = P.R_WIRE
V = P.V_WIRE
L = 300.0              # half-width of the domain; far boundaries are natural (field ~ dipole, decays fast)

E_exact = abs(V) / (r * 1e-3 * np.arccosh(h / r))
q_exact = 2 * np.pi * P.EPS0 * abs(V) / np.arccosh(h / r)
print(f"analytic: E_mean = {E_exact/1e6:.4f} MV/m, |q| = {q_exact*1e9:.4f} nC/m")

rows = []
for h_wire in (0.012, 0.006, 0.003):
    rect = (-L, -L, L, h)                                        # plate = top edge
    circles = [dict(name="wire", x=0.0, z=0.0, r=r, h=h_wire, grow=40.0)]
    mesh = build_mesh(rect, circles, {"far": 10.0})
    b = classify_boundaries(mesh, circles, rect, {"top": "plate"})
    basis, phi = solve_laplace(mesh, [(b["wire"], V), (b["plate"], 0.0)])
    s = surface_field(mesh, phi, b["wire"], (0.0, 0.0))
    err = (s["E_mean"] - E_exact) / E_exact
    errq = (abs(s["charge_per_m"]) - q_exact) / q_exact
    rows.append(dict(h_wire_mm=h_wire, elements=int(mesh.t.shape[1]),
                     E_mean_MV_m=s["E_mean"] / 1e6, E_exact_MV_m=E_exact / 1e6,
                     err_E_pct=100 * err, err_q_pct=100 * errq,
                     E_max_MV_m=s["E_max"] / 1e6, angle_max_deg=s["angle_max_deg"]))
    print(f"h_wire={h_wire:.3f} mm  elements={mesh.t.shape[1]:7d}  "
          f"E_mean={s['E_mean']/1e6:.4f} MV/m  error={100*err:+.2f} %  charge error={100*errq:+.2f} %  "
          f"E_max={s['E_max']/1e6:.3f} at {s['angle_max_deg']:.0f} deg")

(OUT / "benchmark_wire_plane.json").write_text(json.dumps(rows, indent=2))
ok = abs(rows[1]["err_E_pct"]) < 2.0
print("PASS" if ok else "FAIL", f"(|error| < 2 % at h_wire = {P.H_WIRE} mm)")
