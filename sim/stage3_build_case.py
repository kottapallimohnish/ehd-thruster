"""
Stage 3: airflow driven by the ion body force (OpenFOAM, 2D, laminar, transient).

Builds the OpenFOAM case in ./stage3_case:
  1. Mesh (gmsh, extruded one cell thick): open air box, 7 tubes, the two printed frames
     around the 65 mm opening. The 0.13 mm wires are left out as obstacles (they block
     3 % of the width and their drag is negligible next to the tubes).
  2. Body force f = rho_ion * E from stage 2 mapped onto the CFD cells CONSERVATIVELY:
     every stage-2 element's force integral goes to the CFD cell nearest to it, so the
     total force on the air is exactly the stage-2 total.
  3. The packaged OpenFOAM (v1912) cannot compile coded sources, so the cell forces are
     grouped into K cellZones by k-means on (fx, fz); each zone gets its volume-weighted
     mean force through a vectorSemiImplicitSource. Zone means preserve the total force
     exactly; the grouping only smooths the local distribution.
  4. Boundary conditions: far-field boundaries are open (total pressure 0, air can enter
     or leave), tubes and frames are no-slip walls, front/back are 'empty' (2D).
Coordinates: OpenFOAM x = across the opening, y = flow direction (stage-2 z), z = depth.
"""
import json
import shutil
from pathlib import Path

import gmsh
import numpy as np
from scipy.cluster.vq import kmeans2
from scipy.spatial import cKDTree
from skfem import MeshTri, Basis, ElementTriP2

import params as P

import os
HERE = Path(__file__).parent
VARIANT = os.environ.get("STAGE3_VARIANT", "")          # "" = baseline, "coarse" = mesh-sensitivity case
CASE = HERE / ("stage3_case" + (f"_{VARIANT}" if VARIANT else ""))
RES = HERE / "results"

# ---------------------------------------------------------------- settings
BOX = (-150.0, -120.0, 150.0, 300.0)          # x0, y0, x1, y1 in mm
DEPTH = 1.0                                   # mm, extrusion thickness (one cell)
FRAME_INNER, FRAME_OUTER = 32.5, 42.5         # opening radius, frame half-width (85 mm)
WIRE_FRAME_Y = (-5.0, 0.0)                    # wire frame plate (wires on its inner face, y = 0)
TUBE_FRAME_Y = (16.25, 21.25)                 # tube frame plate (inner face 16.25 mm from wires)
H_TUBE, H_FORCE, H_JET, H_FAR = 0.35, 0.6, 1.2, 8.0
if VARIANT.startswith("coarse"):                                  # every cell ~1.5x larger (about 2x fewer cells: 37 892 vs 74 716)
    H_TUBE, H_FORCE, H_JET, H_FAR = 0.525, 0.9, 1.8, 8.0
FORCE_SCALE = float(os.environ.get("STAGE3_FORCE_SCALE", "1"))   # >1: rerun at the loss-corrected force
N_ZONES = 60                                  # each zone = one source term; cost per step grows with zones
NU = 1.67e-5                                  # m^2/s, air at 33 C
END_TIME = 3.0                                # s of simulated time (first run to 2.0 s, extended to 3.0 s)
AVG_START = 1.0                               # s, start of averaging: before ~1 s the anemometer signal still trends upward (start-up)
ANEMO_Y = (P.TUBE_Z + P.R_TUBE) + 60.0        # 60 mm behind the back of the tubes = 77.1 mm
VANE_HALF = 12.5                              # mm, vane radius assumed for averaging


def mesh_case():
    gmsh.initialize()
    gmsh.option.setNumber("General.Verbosity", 1)
    occ = gmsh.model.occ
    x0, y0, x1, y1 = BOX
    box = occ.addRectangle(x0, y0, 0, x1 - x0, y1 - y0)
    tools = [occ.addDisk(x, P.TUBE_Z, 0, P.R_TUBE, P.R_TUBE) for x in P.TUBE_X]
    for ya, yb in (WIRE_FRAME_Y, TUBE_FRAME_Y):
        for sgn in (-1, 1):
            xa = FRAME_INNER if sgn > 0 else -FRAME_OUTER
            tools.append(occ.addRectangle(xa, ya, 0, FRAME_OUTER - FRAME_INNER, yb - ya))
    occ.cut([(2, box)], [(2, t) for t in tools])
    occ.synchronize()

    tube_curves, other = [], []
    for _, c in gmsh.model.getEntities(1):
        xa, ya, _, xb, yb, _ = gmsh.model.getBoundingBox(1, c)
        if (xb - xa) < 2 * P.R_TUBE + 1e-3 and abs((ya + yb) / 2 - P.TUBE_Z) < 1e-3:
            tube_curves.append(c)
        else:
            other.append(c)
    f = gmsh.model.mesh.field
    d_tube = f.add("Distance"); f.setNumbers(d_tube, "CurvesList", tube_curves); f.setNumber(d_tube, "Sampling", 200)
    t_tube = f.add("Threshold"); f.setNumber(t_tube, "InField", d_tube)
    f.setNumber(t_tube, "SizeMin", H_TUBE); f.setNumber(t_tube, "SizeMax", H_FAR)
    f.setNumber(t_tube, "DistMin", 0.5); f.setNumber(t_tube, "DistMax", 60)
    b_force = f.add("Box")
    f.setNumber(b_force, "VIn", H_FORCE); f.setNumber(b_force, "VOut", H_FAR)
    f.setNumber(b_force, "XMin", -45); f.setNumber(b_force, "XMax", 45)
    f.setNumber(b_force, "YMin", -22); f.setNumber(b_force, "YMax", 24); f.setNumber(b_force, "Thickness", 40)
    b_jet = f.add("Box")
    f.setNumber(b_jet, "VIn", H_JET); f.setNumber(b_jet, "VOut", H_FAR)
    f.setNumber(b_jet, "XMin", -60); f.setNumber(b_jet, "XMax", 60)
    f.setNumber(b_jet, "YMin", -40); f.setNumber(b_jet, "YMax", 160); f.setNumber(b_jet, "Thickness", 60)
    fm = f.add("Min"); f.setNumbers(fm, "FieldsList", [t_tube, b_force, b_jet]); f.setAsBackgroundMesh(fm)
    for o in ("Mesh.MeshSizeExtendFromBoundary", "Mesh.MeshSizeFromPoints", "Mesh.MeshSizeFromCurvature"):
        gmsh.option.setNumber(o, 0)
    gmsh.option.setNumber("Mesh.Algorithm", 6)

    surf = [s for _, s in gmsh.model.getEntities(2)]
    ext = occ.extrude([(2, s) for s in surf], 0, 0, DEPTH, numElements=[1], recombine=True)
    occ.synchronize()
    vols = [t for d, t in ext if d == 3]
    gmsh.model.addPhysicalGroup(3, vols, name="fluid")
    groups = {"front_back": [], "tubes": [], "frames": [], "farfield": []}
    for _, s in gmsh.model.getEntities(2):
        xa, ya, za, xb, yb, zb = gmsh.model.getBoundingBox(2, s)
        if zb - za < 1e-6:
            groups["front_back"].append(s)
        elif (xb - xa) < 2 * P.R_TUBE + 1e-3 and abs((ya + yb) / 2 - P.TUBE_Z) < 1e-3:
            groups["tubes"].append(s)
        elif min(abs(xa - x0), abs(xb - x1), abs(ya - y0), abs(yb - y1)) < 1e-6 and \
                (abs(xa - xb) < 1e-6 and abs(xa - x0) + abs(xb - x0) < 1e-6 or abs(xa - xb) < 1e-6 and abs(xa - x1) < 1e-6
                 or abs(ya - yb) < 1e-6 and (abs(ya - y0) < 1e-6 or abs(ya - y1) < 1e-6)):
            groups["farfield"].append(s)
        else:
            groups["frames"].append(s)
    for name, tags in groups.items():
        gmsh.model.addPhysicalGroup(2, tags, name=name)
    gmsh.model.mesh.generate(3)
    gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)
    gmsh.option.setNumber("Mesh.Binary", 0)
    gmsh.write(str(CASE / "mesh.msh"))
    gmsh.finalize()
    return {k: len(v) for k, v in groups.items()}


def foam_header(cls, obj, loc=None):
    l = f'    location    "{loc}";\n' if loc else ""
    return ("FoamFile\n{\n    version     2.0;\n    format      ascii;\n"
            f"    class       {cls};\n{l}    object      {obj};\n}}\n\n")


def write(rel, text):
    p = CASE / rel
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(text)


def write_dicts():
    write("constant/transportProperties", foam_header("dictionary", "transportProperties", "constant") +
          f"transportModel  Newtonian;\nnu              {NU};\n")
    write("constant/turbulenceProperties", foam_header("dictionary", "turbulenceProperties", "constant") +
          "simulationType  laminar;\n")
    write("system/controlDict", foam_header("dictionary", "controlDict", "system") + f"""
application     pimpleFoam;
startFrom       startTime;
startTime       0;
stopAt          endTime;
endTime         {END_TIME};
deltaT          1e-5;
adjustTimeStep  yes;
maxCo           0.9;
writeControl    adjustable;
writeInterval   0.02;
purgeWrite      0;
writeFormat     ascii;
writePrecision  7;
timeFormat      general;
runTimeModifiable true;

// Function objects (probes, fieldAverage) crash in this packaged v1912 build ('sha1' error),
// so fields are written every 0.02 s and sampled/averaged in Python (stage3_post.py).
""")
    write("system/fvSchemes", foam_header("dictionary", "fvSchemes", "system") + """
ddtSchemes { default Euler; }
gradSchemes { default Gauss linear; grad(U) cellLimited Gauss linear 1; }
divSchemes
{
    default none;
    div(phi,U) Gauss linearUpwind grad(U);
    div((nuEff*dev2(T(grad(U))))) Gauss linear;
}
laplacianSchemes { default Gauss linear limited corrected 0.5; }
interpolationSchemes { default linear; }
snGradSchemes { default limited corrected 0.5; }
""")
    write("system/fvSolution", foam_header("dictionary", "fvSolution", "system") + """
solvers
{
    p      { solver GAMG; smoother GaussSeidel; tolerance 1e-6; relTol 0.01; }
    pFinal { $p; relTol 0; }
    "(U|UFinal)" { solver smoothSolver; smoother symGaussSeidel; tolerance 1e-7; relTol 0; }
}
PIMPLE
{
    nOuterCorrectors 1;
    nCorrectors 2;
    nNonOrthogonalCorrectors 1;
}
""")
    common = "    front_back { type empty; }\n"
    write("0/U", foam_header("volVectorField", "U", "0") + """dimensions [0 1 -1 0 0 0 0];
internalField uniform (0 0 0);
boundaryField
{
    farfield { type pressureInletOutletVelocity; value uniform (0 0 0); }
    tubes    { type noSlip; }
    frames   { type noSlip; }
""" + common + "}\n")
    write("0/p", foam_header("volScalarField", "p", "0") + """dimensions [0 2 -2 0 0 0 0];
internalField uniform 0;
boundaryField
{
    farfield { type totalPressure; p0 uniform 0; value uniform 0; }
    tubes    { type zeroGradient; }
    frames   { type zeroGradient; }
""" + common + "}\n")


def read_internal_field(path, vector):
    txt = Path(path).read_text()
    i = txt.index("internalField")
    j = txt.index("(", txt.index("List<", i))
    n = int(txt[txt.index(">", i) + 1:j].split()[-1]) if False else None
    body = txt[j + 1:txt.index("\n)\n", j)]
    if vector:
        return np.array([list(map(float, l.strip("() ").split())) for l in body.strip().splitlines()])
    return np.array(list(map(float, body.split())))


def _read_foam_list(path, kind):
    txt = Path(path).read_text()
    txt = txt[txt.index("}", txt.index("FoamFile")) + 1:]
    lines = [l for l in txt.splitlines() if l.strip() and not l.strip().startswith("//")]
    n = int(lines[0]); body = lines[2:2 + n]
    if kind == "label":
        return np.array(list(map(int, body)))
    if kind == "point":
        return np.array([list(map(float, l.strip("()").split())) for l in body])
    return [list(map(int, l[l.index("(") + 1:l.rindex(")")].split())) for l in body]   # faces


def cell_geometry():
    """Cell centres (mm) and volumes (m^3) from the polyMesh files. Every cell is a triangular
    prism (one extruded layer), so the centre is the mean of its 6 points and the volume is
    triangle area x depth."""
    pm = CASE / "constant/polyMesh"
    pts = _read_foam_list(pm / "points", "point")
    faces = _read_foam_list(pm / "faces", "face")
    own = _read_foam_list(pm / "owner", "label")
    nei = _read_foam_list(pm / "neighbour", "label")
    ncell = own.max() + 1
    cell_pts = [set() for _ in range(ncell)]
    for f, c in enumerate(own):
        cell_pts[c].update(faces[f])
    for f, c in enumerate(nei):
        cell_pts[c].update(faces[f])
    C = np.zeros((ncell, 3)); V = np.zeros(ncell)
    for c, ps in enumerate(cell_pts):
        q = pts[sorted(ps)]
        assert q.shape[0] == 6, "expected prism cells"
        C[c] = q.mean(0)
        tri = q[q[:, 2] < q[:, 2].mean()][:, :2]
        a = 0.5 * abs((tri[1, 0] - tri[0, 0]) * (tri[2, 1] - tri[0, 1]) - (tri[2, 0] - tri[0, 0]) * (tri[1, 1] - tri[0, 1]))
        V[c] = a * (q[:, 2].max() - q[:, 2].min())
    return C * 1e3, V


def stage2_element_forces():
    """Force integral (N per metre depth) of rho*E over every stage-2 element, and its centroid (mm)."""
    d = np.load(RES / "stage2_solution.npz")
    mesh = MeshTri(d["p"], d["t"])
    basis = Basis(mesh, ElementTriP2())
    import stage2_space_charge as S2
    rho2 = S2.p1_to_p2(mesh, basis, d["rho_v"])
    E = -np.asarray(basis.interpolate(d["phi"]).grad) * 1e3            # V/m, (2, nel, nq)
    rq = np.asarray(basis.interpolate(rho2))                           # C/m^3
    dx = np.asarray(basis.dx) * 1e-6                                   # m^2
    F = (rq[None] * E * dx[None]).sum(axis=2)                          # N/m per element, (2, nel)
    cen = mesh.p[:, mesh.t].mean(axis=1)                               # mm
    return F, cen


def map_forces():
    F, cen = stage2_element_forces()
    C, V = cell_geometry()                                             # mm, m^3
    tree = cKDTree(C[:, :2])
    _, idx = tree.query(cen.T)
    Fc = np.zeros((C.shape[0], 2))
    np.add.at(Fc, idx, F.T)                                            # N/m per CFD cell
    area = V / (DEPTH * 1e-3)                                          # m^2
    f_spec = FORCE_SCALE * Fc / area[:, None] / P.RHO_AIR              # m/s^2 (force per unit mass)
    Fc = FORCE_SCALE * Fc

    active = np.where(np.hypot(*f_spec.T) > 1e-6 * np.abs(f_spec).max())[0]
    X = f_spec[active]
    scale = X.std(axis=0)
    rng = np.random.default_rng(0)
    centroids, lab = kmeans2(X / scale, N_ZONES, minit="++", seed=rng)
    zones = []
    for k in range(N_ZONES):
        cells = active[lab == k]
        if cells.size == 0:
            continue
        w = V[cells]
        mean = (f_spec[cells] * w[:, None]).sum(0) / w.sum()          # volume-weighted: keeps total force
        zones.append((cells, mean))

    lines = [foam_header("regIOobject", "cellZones", "constant/polyMesh"), f"{len(zones)}\n(\n"]
    opts = [foam_header("dictionary", "fvOptions", "constant")]
    for k, (cells, mean) in enumerate(zones):
        lines.append(f"ionForce{k}\n{{\n    type cellZone;\ncellLabels List<label> {cells.size}\n(\n")
        lines.append("\n".join(map(str, cells)))
        lines.append("\n)\n;\n}\n")
        opts.append(f"""ionForce{k}
{{
    type            vectorSemiImplicitSource;
    active          true;
    selectionMode   cellZone;
    cellZone        ionForce{k};
    volumeMode      specific;
    injectionRateSuSp {{ U (({mean[0]:.6e} {mean[1]:.6e} 0) 0); }}
}}
""")
    lines.append(")\n")
    write("constant/polyMesh/cellZones", "".join(lines))
    write("constant/fvOptions", "".join(opts))

    T_exact = FORCE_SCALE * F[1].sum()
    T_cells = (Fc[:, 1]).sum()
    T_zones = sum((m[1] * P.RHO_AIR * (V[c] / (DEPTH * 1e-3))).sum() for c, m in zones)
    rel_err_local = np.sqrt(np.mean([(np.hypot(*(f_spec[c] - m).T) ** 2 * V[c]).sum() for c, m in zones]))
    info = dict(thrust_per_m_stage2=T_exact, thrust_per_m_on_cells=T_cells, thrust_per_m_in_zones=T_zones,
                cells_total=int(C.shape[0]), cells_with_force=int(active.size), zones=len(zones),
                max_specific_force=float(np.abs(f_spec).max()))
    (CASE / "force_mapping.json").write_text(json.dumps(info, indent=2))
    return info


if __name__ == "__main__":
    import subprocess, sys
    if CASE.exists():
        shutil.rmtree(CASE)
    (CASE / "system").mkdir(parents=True)
    write_dicts()
    print("patches:", mesh_case())
    env = f"source {os.environ.get('FOAM_BASHRC', '/usr/share/openfoam/etc/bashrc')} >/dev/null 2>&1; "
    def run(cmd):
        r = subprocess.run(["bash", "-c", env + f"cd {CASE} && {cmd}"], capture_output=True, text=True)
        if r.returncode:
            print(r.stdout[-3000:], r.stderr[-3000:]); sys.exit(1)
        return r.stdout
    run("gmshToFoam mesh.msh > log.gmshToFoam 2>&1")
    run("transformPoints -scale '(0.001 0.001 0.001)' > log.transformPoints 2>&1")   # gmsh works in mm, OpenFOAM in m
    # patch types
    run("foamDictionary -entry entry0/front_back/type -set empty constant/polyMesh/boundary >/dev/null 2>&1 || true")
    bnd = (CASE / "constant/polyMesh/boundary").read_text()
    for name, typ in (("front_back", "empty"), ("tubes", "wall"), ("frames", "wall"), ("farfield", "patch")):
        i = bnd.index(f"    {name}\n")
        j = bnd.index("type", i)
        k = bnd.index(";", j)
        bnd = bnd[:j] + f"type            {typ}" + bnd[k:]
        if typ != "empty":
            pass
    bnd = bnd.replace("physicalType    patch;", "")
    (CASE / "constant/polyMesh/boundary").write_text(bnd)
    print(run("checkMesh > log.checkMesh 2>&1; grep -E 'cells:|Mesh OK|Failed|non-orthogonality' log.checkMesh"))
    print(json.dumps(map_forces(), indent=2))
    write("system/decomposeParDict", foam_header("dictionary", "decomposeParDict", "system") +
          "numberOfSubdomains 2;\nmethod scotch;\n")
