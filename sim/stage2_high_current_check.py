"""
Consistency check after stage 3: the flow losses imply roughly 3x the stage-2 current.
Re-run the space-charge model at that current and check that T = I*d_axial/mu still holds
(i.e. that the stronger space charge does not change the thrust per ampere).
"""
import json
import numpy as np
import params as P
import stage2_space_charge as S
from skfem import Basis, ElementTriP2, asm, condense, solve
from skfem.models.poisson import laplace, mass
from fem_tools import build_mesh, classify_boundaries

# headline nominal current from sensitivity.py (run sensitivity.py first; rerun it afterwards,
# since it reads d_axial from this file - the dependence is < 0.1 %)
I = json.loads((S.RES / "sensitivity.json").read_text())["nominal"]["current_uA"] * 1e-6
rect = (P.X_MIN, P.Z_MIN, P.X_MAX, P.Z_MAX)
wires, tubes = S.thruster_circles()
mesh = build_mesh(rect, wires + tubes, {"far": P.H_FAR})
b = classify_boundaries(mesh, wires + tubes, rect)
basis = Basis(mesh, ElementTriP2()); A, M = asm(laplace, basis), asm(mass, basis)
phi_bc = basis.zeros(); D = []
for w in wires:
    d = basis.get_dofs(b[w["name"]]).all(); phi_bc[d] = P.V_WIRE; D.append(d)
for t in tubes:
    D.append(basis.get_dofs(b[t["name"]]).all())
D = np.unique(np.concatenate(D))
phi_lap = solve(*condense(A, np.zeros(A.shape[0]), x=phi_bc, D=D))
sol = S.run(I, mesh, b, basis, A, M, D, phi_bc, phi_lap, maxit=25, verbose=False)
# note: at fixed emission, the wire surface fields drop well below onset here. In real steady
# corona the surface field stays near onset (Kaptzov), so this low field is a property of the
# imposed-emission model, implying the nicks sustain emission at a few times local field enhancement.
dg = S.diagnostics(sol, mesh, b, basis, phi_lap)
from matplotlib.tri import Triangulation, LinearTriInterpolator
nv = mesh.p.shape[1]
f = LinearTriInterpolator(Triangulation(mesh.p[0], mesh.p[1], mesh.t.T), (sol["phi"] - phi_lap)[:nv])
X, Z = np.meshgrid(np.linspace(-20, 20, 81), np.linspace(1, 10, 19))
out = dict(I_uA=I * 1e6, thrust_mN=dg["T"] * 1e3, d_axial_mm=dg["d_axial_eff_mm"],
           d_axial_mm_at_45uA=json.loads((S.RES / "stage2_summary.json").read_text())["effective_axial_length_mm"],
           gap_mean_potential_change_V=float(np.nanmean(np.ma.filled(f(X, Z), np.nan))),
           surface_field_inner_wire_MV_m=dg["Es"][5][1], surface_field_end_wire_MV_m=dg["Es"][0][1])
(S.RES / "stage2_high_current_check.json").write_text(json.dumps(out, indent=2))
print(json.dumps(out, indent=2))
