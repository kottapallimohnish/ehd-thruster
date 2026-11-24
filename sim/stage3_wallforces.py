"""
Momentum budget check for stage 3: time-averaged force of the air on the tubes and frames
(pressure + viscous), compared with the total ion force pushing the air.
Drag on the walls is momentum the ions gave the air that never reaches the jet.
"""
import json, re
import numpy as np
import params as P
import stage3_build_case as B
from stage3_post import read_U

CASE = B.CASE
pm = CASE / "constant/polyMesh"
pts = B._read_foam_list(pm / "points", "point")
faces = B._read_foam_list(pm / "faces", "face")
own = B._read_foam_list(pm / "owner", "label")
C, V = B.cell_geometry(); C = C * 1e-3
bnd = (pm / "boundary").read_text()

def patch(name):
    blk = bnd[bnd.index(f"    {name}\n"):]
    n = int(re.search(r"nFaces\s+(\d+);", blk).group(1)); s = int(re.search(r"startFace\s+(\d+);", blk).group(1))
    return range(s, s + n)

def read_p(tdir):
    txt = (tdir / "p").read_text(); i = txt.index("internalField"); j = txt.index("(", txt.index("List<scalar>", i))
    return np.array(txt[j + 1: txt.index("\n)\n", j)].split(), float)

times = [d for d in CASE.iterdir() if re.fullmatch(r"[0-9.]+", d.name) and float(d.name) >= B.AVG_START - 1e-9]
pbar = np.mean([read_p(d) for d in times], axis=0)
Ubar = np.mean([read_U(d) for d in times], axis=0)
out = {}
for name in ("tubes", "frames"):
    Fp = np.zeros(3); Fv = np.zeros(3)
    for f in patch(name):
        q = pts[faces[f]]
        Avec = 0.5 * np.cross(q[2] - q[0], q[3] - q[1]) if len(q) == 4 else 0.5 * np.cross(q[1] - q[0], q[2] - q[0])
        A = np.linalg.norm(Avec); n = Avec / A                           # points out of the fluid
        c = own[f]; fc = q.mean(0)
        dn = abs(np.dot(fc - C[c], n))
        Fp += P.RHO_AIR * pbar[c] * A * n                                # p is kinematic (p/rho)
        Ut = Ubar[c] - np.dot(Ubar[c], n) * n
        Fv += P.RHO_AIR * B.NU * Ut / dn * A                             # wall shear, drags wall along flow
    depth = B.DEPTH * 1e-3
    out[name] = dict(pressure_N_per_m=float(Fp[1] / depth), viscous_N_per_m=float(Fv[1] / depth),
                     total_N_per_m=float((Fp[1] + Fv[1]) / depth))
T = json.loads((CASE / "force_mapping.json").read_text())["thrust_per_m_in_zones"]
out["ion_force_N_per_m"] = T
out["wall_drag_fraction_of_ion_force"] = (out["tubes"]["total_N_per_m"] + out["frames"]["total_N_per_m"]) / T
(B.HERE / f"results/stage3_wallforces{('_' + B.VARIANT) if B.VARIANT else ''}.json").write_text(json.dumps(out, indent=2))
print(json.dumps(out, indent=2))
