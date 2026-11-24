"""
Shared finite-element tools: mesh a rectangle with circular conductors cut out (gmsh),
solve Laplace's equation with Dirichlet conductors (scikit-fem, P2 elements), and
evaluate the electric field on conductor surfaces.

Units: mm and volts inside; fields are converted to V/m when reported.
"""
import gmsh
import numpy as np
from skfem import MeshTri, Basis, FacetBasis, ElementTriP2, asm, condense, solve
from skfem.models.poisson import laplace


# ----------------------------------------------------------------------------- mesh
def build_mesh(rect, circles, sizes, verbose=False):
    """
    rect    : (x0, z0, x1, z1) air domain.
    circles : list of dicts {name, x, z, r, h, grow}: conductor of radius r cut out of the air,
              target element size h on its surface, growing to sizes['far'] over distance grow.
    sizes   : {'far': h_far}
    Returns a scikit-fem MeshTri (straight-sided triangles; the P2 *field* is second order).
    """
    x0, z0, x1, z1 = rect
    gmsh.initialize()
    gmsh.option.setNumber("General.Verbosity", 2 if verbose else 0)
    gmsh.model.add("domain")
    occ = gmsh.model.occ
    box = occ.addRectangle(x0, z0, 0, x1 - x0, z1 - z0)
    disks = [occ.addDisk(c["x"], c["z"], 0, c["r"], c["r"]) for c in circles]
    occ.cut([(2, box)], [(2, d) for d in disks])
    occ.synchronize()

    # Identify each circle's boundary curve by its bounding box.
    curves = [tag for (_, tag) in gmsh.model.getEntities(1)]
    fields = []
    for c in circles:
        tags = []
        for tag in curves:
            xa, za, _, xb, zb, _ = gmsh.model.getBoundingBox(1, tag)
            if abs((xa + xb) / 2 - c["x"]) < 1e-6 * max(1, abs(c["x"])) + 1e-6 and \
               abs((xb - xa) / 2 - c["r"]) < 1e-3 * c["r"] + 1e-9 and \
               abs((za + zb) / 2 - c["z"]) < 1e-6 * max(1, abs(c["z"])) + 1e-6:
                tags.append(tag)
        assert tags, f"curve for circle {c['name']} not found"
        # Distance field from this circle, then a size that grows linearly with distance.
        fd = gmsh.model.mesh.field.add("Distance")
        gmsh.model.mesh.field.setNumbers(fd, "CurvesList", tags)
        gmsh.model.mesh.field.setNumber(fd, "Sampling", 400)
        ft = gmsh.model.mesh.field.add("Threshold")
        gmsh.model.mesh.field.setNumber(ft, "InField", fd)
        gmsh.model.mesh.field.setNumber(ft, "SizeMin", c["h"])
        gmsh.model.mesh.field.setNumber(ft, "SizeMax", sizes["far"])
        gmsh.model.mesh.field.setNumber(ft, "DistMin", 0.0)
        gmsh.model.mesh.field.setNumber(ft, "DistMax", c["grow"])
        fields.append(ft)
    fmin = gmsh.model.mesh.field.add("Min")
    gmsh.model.mesh.field.setNumbers(fmin, "FieldsList", fields)
    gmsh.model.mesh.field.setAsBackgroundMesh(fmin)
    for opt in ("Mesh.MeshSizeExtendFromBoundary", "Mesh.MeshSizeFromPoints",
                "Mesh.MeshSizeFromCurvature"):
        gmsh.option.setNumber(opt, 0)
    gmsh.option.setNumber("Mesh.Algorithm", 6)   # Frontal-Delaunay: good-quality triangles
    gmsh.model.mesh.generate(2)

    node_tags, coords, _ = gmsh.model.mesh.getNodes()
    xyz = coords.reshape(-1, 3)
    index = {t: i for i, t in enumerate(node_tags)}
    _, elem_nodes = gmsh.model.mesh.getElementsByType(2)   # 3-node triangles
    tri = np.vectorize(index.get)(elem_nodes.reshape(-1, 3)).T
    gmsh.finalize()

    used = np.unique(tri)                                   # drop unused nodes (curve-only)
    remap = -np.ones(xyz.shape[0], int)
    remap[used] = np.arange(used.size)
    return MeshTri(xyz[used, :2].T, remap[tri])


def classify_boundaries(mesh, circles, rect, edge_names=None):
    """Assign boundary facets to named conductors by geometry (distance from circle centres)."""
    bf = mesh.boundary_facets()
    mid = mesh.p[:, mesh.facets[:, bf]].mean(axis=1)
    out, taken = {}, np.zeros(bf.size, bool)
    for c in circles:
        d = np.hypot(mid[0] - c["x"], mid[1] - c["z"])
        sel = np.abs(d - c["r"]) < 0.02 * c["r"]
        out[c["name"]] = bf[sel]
        taken |= sel
    x0, z0, x1, z1 = rect
    tol = 1e-9 + 1e-6 * max(abs(v) for v in rect)
    edges = {"left": np.abs(mid[0] - x0) < tol, "right": np.abs(mid[0] - x1) < tol,
             "bottom": np.abs(mid[1] - z0) < tol, "top": np.abs(mid[1] - z1) < tol}
    for k, sel in edges.items():
        name = (edge_names or {}).get(k, k)
        out[name] = np.concatenate([out.get(name, np.array([], int)), bf[sel & ~taken]])
        taken |= sel
    assert taken.all(), "some boundary facets were not classified"
    return out


# ---------------------------------------------------------------------------- solve
def solve_laplace(mesh, dirichlet):
    """
    Laplace equation  -div(grad phi) = 0  with P2 (quadratic) elements.
    dirichlet: {facet_array: potential}. All other boundaries are natural
    (zero normal derivative = no field crossing them).
    """
    basis = Basis(mesh, ElementTriP2())
    A = asm(laplace, basis)
    phi = basis.zeros()
    D = []
    for facets, value in dirichlet:
        dofs = basis.get_dofs(facets).all()
        phi[dofs] = value
        D.append(dofs)
    D = np.unique(np.concatenate(D))
    phi = solve(*condense(A, x=phi, D=D))
    return basis, phi


# -------------------------------------------------------------- surface diagnostics
def surface_field(mesh, phi, facets, centre, eps0=8.854e-12):
    """
    Electric field on a conductor surface.
    Returns dict: mean |E| (V/m, length-weighted), max |E|, angle of max (deg, 0 = +z,
    i.e. pointing downstream towards the collector), and charge per unit length (C/m)
    from Gauss's law q = eps0 * closed-integral(E . n_out) ds.
    """
    fb = FacetBasis(mesh, ElementTriP2(), facets=facets, intorder=6)
    g = fb.interpolate(phi).grad                       # V/mm, shape (2, nf, nq)
    E = -np.asarray(g) * 1e3                            # V/m
    Emag = np.hypot(E[0], E[1])
    w = np.asarray(fb.dx)                               # quadrature weight * length (mm)
    n_out_of_conductor = -np.asarray(fb.normals)        # skfem normal points out of the air
    En = (E * n_out_of_conductor).sum(axis=0)
    xq = np.asarray(fb.global_coordinates())
    ang = np.degrees(np.arctan2(xq[0] - centre[0], xq[1] - centre[1]))
    imax = np.unravel_index(np.argmax(Emag), Emag.shape)
    return {
        "E_mean": float((Emag * w).sum() / w.sum()),
        "E_max": float(Emag[imax]),
        "angle_max_deg": float(ang[imax]),
        "charge_per_m": float(eps0 * (En * w).sum() * 1e-3),
        "perimeter_mm": float(w.sum()),
        "_angles": ang.ravel(), "_Emag": Emag.ravel(),
    }
