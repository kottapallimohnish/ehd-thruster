"""
Sensitivity check: does the size of the air box change the stage 2 answer?
Re-runs stage 2 at the calibrated current in a much larger box (x +-70 mm, z -60..100 mm)
and compares thrust, effective axial length and wire surface fields with the baseline.
"""
import json
from pathlib import Path
import numpy as np
import params as P
P.X_MIN, P.X_MAX, P.Z_MIN, P.Z_MAX = -70.0, 70.0, -60.0, 100.0
import stage2_space_charge as S
from skfem import Basis, ElementTriP2, asm, condense, solve
from skfem.models.poisson import laplace, mass
from fem_tools import build_mesh, classify_boundaries

base = json.loads((S.RES / "stage2_summary.json").read_text())
I = base["total_current_uA"] * 1e-6
rect = (P.X_MIN, P.Z_MIN, P.X_MAX, P.Z_MAX)
wires, tubes = S.thruster_circles()
mesh = build_mesh(rect, wires + tubes, {"far": 3.0})
b = classify_boundaries(mesh, wires + tubes, rect)
basis = Basis(mesh, ElementTriP2())
A, M = asm(laplace, basis), asm(mass, basis)
phi_bc = basis.zeros(); D = []
for w in wires:
    d = basis.get_dofs(b[w["name"]]).all(); phi_bc[d] = P.V_WIRE; D.append(d)
for t in tubes:
    d = basis.get_dofs(b[t["name"]]).all(); D.append(d)
D = np.unique(np.concatenate(D))
phi_lap = solve(*condense(A, np.zeros(A.shape[0]), x=phi_bc, D=D))
sol = S.run(I, mesh, b, basis, A, M, D, phi_bc, phi_lap, maxit=25)
dg = S.diagnostics(sol, mesh, b, basis, phi_lap)
out = dict(box="x +-70, z -60..100 mm", I_uA=I * 1e6,
           thrust_mN_large_box=dg["T"] * 1e3, thrust_mN_baseline=base["thrust_mN"],
           thrust_change_pct=100 * (dg["T"] * 1e3 / base["thrust_mN"] - 1),
           d_axial_mm_large_box=dg["d_axial_eff_mm"], d_axial_mm_baseline=base["effective_axial_length_mm"],
           surface_field_large_box=[round(c, 3) for _, c in dg["Es"]],
           surface_field_baseline=[c for _, c in base["wire_surface_field_MV_m_laplace_vs_spacecharge"]],
           max_dphi_large_box_V=dg["max_dphi_space_charge"])
(S.RES / "stage2_domain_check.json").write_text(json.dumps(out, indent=2))
print(json.dumps(out, indent=2))
