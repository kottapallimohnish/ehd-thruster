"""
Stage 2: the ions' own charge (space charge), coupled to the electric field.

Equations (steady state, negative ions with mobility mu, diffusion and air motion neglected):
    Poisson:              -div(grad phi) = rho / eps0
    charge conservation:  div(rho * mu * E) = 0,   E = -grad phi

Method of characteristics. Expanding the conservation law with Poisson's equation gives,
along an ion path (a field line),
    d|rho|/dt = -mu * rho^2 / eps0     =>     1/|rho| = 1/|rho0| + tau/eps0,
    tau = integral( ds / |E| )  from the wire surface to the point.
So the charge density anywhere follows from (a) which wire its field line starts on,
(b) the charge density rho0 injected at that wire, and (c) tau along the line.
tau is found by tracing every mesh vertex BACK along E to the wire it came from.

Emission model (nicked wires): the cuts, not the smooth-wire field, set where ions are
made, and they are spread along every wire. So each wire emits the SAME current per metre
of wire; rho0 for each wire is set to deliver it. The total current is the single
calibration parameter, chosen so the computed thrust equals the thrust inferred from the
measured jet speed (T = 1/2 rho_air A v^2 = 2.73 mN).

Loop: field -> back-trace -> rho -> Poisson -> field ... until the potential stops changing.
Then thrust T = integral(rho * E_z) dA per metre, times each wire's real active length.
"""
import json
import time
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
from matplotlib.tri import Triangulation
from skfem import Basis, ElementTriP2, asm, condense, solve, BilinearForm, LinearForm
from skfem.models.poisson import laplace, mass
from skfem.helpers import dot

import params as P
from fem_tools import build_mesh, classify_boundaries, surface_field
from solve_laplace import thruster_circles

HERE = Path(__file__).parent
RES, FIG = HERE / "results", HERE / "figures"
RES.mkdir(exist_ok=True)
FIG.mkdir(exist_ok=True)
EPS0 = P.EPS0


# ---------------------------------------------------------------- geometry helpers
def wire_index_and_gap(pts):
    """Nearest wire and distance from its surface (mm) for points (2, N)."""
    d = np.hypot(pts[0][:, None] - P.WIRE_X[None, :], pts[1][:, None] - P.WIRE_Z) - P.R_WIRE
    i = d.argmin(axis=1)
    return i, d[np.arange(d.shape[0]), i]


def in_tube(pts):
    d = np.hypot(pts[0][:, None] - P.TUBE_X[None, :], pts[1][:, None] - P.TUBE_Z)
    return d.min(axis=1) <= P.R_TUBE * 1.003


# ---------------------------------------------------------------- field sampling
class Field:
    """
    E-field (V/m) at arbitrary points. The P2 gradient is L2-projected to a continuous field,
    its vertex values are interpolated linearly inside each triangle (matplotlib's fast
    triangle finder). Points that fall in a hole (inside a wire or tube) return NaN.
    """
    def __init__(self, basis, phi):
        from matplotlib.tri import LinearTriInterpolator
        mesh = basis.mesh
        g = basis.interpolate(phi).grad
        nv = mesh.p.shape[1]
        gx, gz = basis.project(g[0])[:nv], basis.project(g[1])[:nv]      # grad phi, V/mm
        tri = Triangulation(mesh.p[0], mesh.p[1], mesh.t.T)
        self.ix = LinearTriInterpolator(tri, -gx * 1e3)
        self.iz = LinearTriInterpolator(tri, -gz * 1e3)

    def E(self, pts):
        ex = np.ma.filled(self.ix(pts[0], pts[1]), np.nan)
        ez = np.ma.filled(self.iz(pts[0], pts[1]), np.nan)
        return np.vstack([ex, ez])


def back_trace(field, pts, max_len=600.0, max_steps=6000):
    """
    Trace each point backwards along the ion path (along E, towards the wires).
    Returns origin wire index (-1 if none reached) and tau = integral ds/|E| in m^2/V.
    """
    n = pts.shape[1]
    p = pts.copy()
    origin = -np.ones(n, int)
    tau = np.zeros(n)
    length = np.zeros(n)
    w0, g0 = wire_index_and_gap(p)
    done = g0 <= 0.02 * P.R_WIRE                 # vertices on a wire surface: tau = 0
    origin[done] = w0[done]
    for _ in range(max_steps):
        run = np.where(~done)[0]
        if run.size == 0:
            break
        q = p[:, run]
        _, gap = wire_index_and_gap(q)
        h = np.clip(0.15 * (gap + P.R_WIRE), 0.003, 0.25)            # mm
        E1 = field.E(q)
        m1 = np.hypot(*E1)
        bad = ~np.isfinite(m1) | (m1 == 0)
        if bad.any():                                    # started inside a hole: snap to nearest wire
            wb, _ = wire_index_and_gap(q[:, bad])
            origin[run[bad]] = wb
            done[run[bad]] = True
            m1 = np.where(bad, 1.0, m1); E1 = np.where(bad, 0.0, E1)
        qm = q + 0.5 * h * E1 / m1
        qm[0] = np.clip(qm[0], P.X_MIN + 1e-4, P.X_MAX - 1e-4)
        qm[1] = np.clip(qm[1], P.Z_MIN + 1e-4, P.Z_MAX - 1e-4)
        wi, gm = wire_index_and_gap(qm)
        hit = (gm <= 0.0) | in_tube(qm)
        qn = qm.copy()
        Em = np.empty_like(E1)
        ok = ~hit
        if ok.any():
            Em[:, ok] = field.E(qm[:, ok])
            mm = np.hypot(*Em[:, ok])
            mm = np.where(np.isfinite(mm) & (mm > 0), mm, m1[ok])
            Em[:, ok] = np.where(np.isfinite(Em[:, ok]), Em[:, ok], E1[:, ok])
            qn[:, ok] = q[:, ok] + h[ok] * Em[:, ok] / mm
            tau[run[ok]] += h[ok] * 1e-3 / mm                         # ds [m] / |E| [V/m]
        tau[run[hit]] += 0.5 * h[hit] * 1e-3 / m1[hit]
        length[run] += np.where(hit, 0.5 * h, h)
        p[:, run] = qn
        wi_n, gn = wire_index_and_gap(qn)
        arrived = (gn <= 0.02 * P.R_WIRE) | (hit & (gm <= 0.0))
        origin[run[arrived]] = np.where(hit[arrived], wi[arrived], wi_n[arrived])
        lost = (~arrived) & ((hit & ~(gm <= 0.0)) | (length[run] > max_len))
        done[run[arrived | lost]] = True
    return origin, tau


# ---------------------------------------------------------------- Poisson solve
def p1_to_p2(mesh, basis, vals):
    """Embed vertex values (piecewise linear) into the P2 space: edge midpoints = mean of ends."""
    out = basis.zeros()
    nv = mesh.p.shape[1]
    out[:nv] = vals                                     # P2 vertex dofs come first in scikit-fem
    out[nv:nv + mesh.facets.shape[1]] = 0.5 * (vals[mesh.facets[0]] + vals[mesh.facets[1]])
    return out


def solve_poisson(basis, A, M, D, phi_bc, rho_p2):
    """-lap(phi) = rho/eps0, lengths in mm: rho/eps0 [V/m^2] * 1e-6 -> V/mm^2."""
    b = M @ (rho_p2 / EPS0 * 1e-6)
    return solve(*condense(A, b, x=phi_bc, D=D))


@BilinearForm
def _unused(u, v, w):
    return u * v


# ---------------------------------------------------------------- main solver
def run(I_total, mesh, b, basis, A, M, D, phi_bc, phi_init, verbose=True, tol=0.5, maxit=30):
    wires, tubes = thruster_circles()
    Lw = P.WIRE_LEN_MM * 1e-3
    Iprime = I_total / Lw.sum()                      # A per metre of wire, equal for every wire
    phi = phi_init.copy()
    nv = mesh.p.shape[1]
    hist = []
    for it in range(maxit):
        field = Field(basis, phi)
        origin, tau = back_trace(field, mesh.p)
        # rho0 per wire from the current it must emit: I' = rho0 * mu * closed-integral(E.n) ds
        rho0 = np.zeros(P.N_WIRES)
        for i, w in enumerate(wires):
            sf = surface_field(mesh, phi, b[w["name"]], (w["x"], w["z"]))
            flux = abs(sf["charge_per_m"]) / EPS0      # closed-integral E.n ds  [V], Gauss's law
            rho0[i] = Iprime / (P.MU_ION * flux)
        rho_v = np.zeros(nv)
        ok = origin >= 0
        rho_v[ok] = -1.0 / (1.0 / rho0[origin[ok]] + tau[ok] / EPS0)   # negative ions
        phi_new = solve_poisson(basis, A, M, D, phi_bc, p1_to_p2(mesh, basis, rho_v))
        # Under-relaxed fixed-point update. Pointwise changes stay noisy at the ~1 % level because
        # each vertex is re-traced every iteration; the RMS change is the convergence measure.
        change_max = np.abs(phi_new - phi).max()
        change_rms = np.sqrt(np.mean((phi_new - phi) ** 2))
        phi = 0.5 * phi + 0.5 * phi_new
        hist.append(change_rms)
        if verbose:
            print(f"  it {it:2d}: rms |dphi| = {change_rms:7.3f} V  max = {change_max:7.2f} V  "
                  f"unreached vertices = {(~ok).sum()}")
        if change_rms < tol:
            break
    return dict(phi=phi, rho_v=rho_v, origin=origin, tau=tau, rho0=rho0, Iprime=Iprime, hist=hist)


def diagnostics(sol, mesh, b, basis, phi_laplace):
    wires, tubes = thruster_circles()
    phi, rho_v, origin = sol["phi"], sol["rho_v"], sol["origin"]
    Lw = P.WIRE_LEN_MM * 1e-3
    g = basis.interpolate(phi).grad
    Ez = -np.asarray(g[1]) * 1e3                                       # V/m at quadrature points
    out = dict(thrust_per_m=np.zeros(P.N_WIRES))
    for i in range(P.N_WIRES):
        rv = np.where(origin == i, rho_v, 0.0)
        rq = np.asarray(basis.interpolate(p1_to_p2(mesh, basis, rv)))  # C/m^3 at quadrature points
        dx = np.asarray(basis.dx)                                       # mm^2
        out["thrust_per_m"][i] = float((rq * Ez * dx).sum() * 1e-6)     # N per metre of wire
    T = float((out["thrust_per_m"] * Lw).sum())
    I = float(sol["Iprime"] * Lw.sum())
    # current collected by the tubes: closed-integral rho*mu*E.n over tube surfaces
    rho_p2 = p1_to_p2(mesh, basis, rho_v)
    from skfem import FacetBasis
    I_tubes = 0.0
    for t in tubes:
        fb = FacetBasis(mesh, ElementTriP2(), facets=b[t["name"]], intorder=6)
        Eq = -np.asarray(fb.interpolate(phi).grad) * 1e3
        rq = np.asarray(fb.interpolate(rho_p2))
        n_into_tube = -np.asarray(fb.normals)
        I_tubes += float((rq * P.MU_ION * (Eq * n_into_tube).sum(0) * np.asarray(fb.dx)).sum() * 1e-3)
    # surface fields with and without space charge
    Es = []
    for w in wires:
        s1 = surface_field(mesh, phi, b[w["name"]], (w["x"], w["z"]))
        s0 = surface_field(mesh, phi_laplace, b[w["name"]], (w["x"], w["z"]))
        Es.append((s0["E_mean"] / 1e6, s1["E_mean"] / 1e6))
    out.update(T=T, I=I, I_tubes_per_m_times_mean_len=abs(I_tubes) * Lw.mean(),
               I_tubes_per_m=abs(I_tubes), I_emitted_per_m=sol["Iprime"] * P.N_WIRES,
               d_axial_eff_mm=T * P.MU_ION / I * 1e3 if I > 0 else np.nan,
               Es=Es, max_dphi_space_charge=float(np.abs(phi - phi_laplace).max()),
               rho_max=float(np.abs(rho_v).max()), rho0=sol["rho0"].tolist())
    return out


def main():
    t0 = time.time()
    rect = (P.X_MIN, P.Z_MIN, P.X_MAX, P.Z_MAX)
    wires, tubes = thruster_circles()
    mesh = build_mesh(rect, wires + tubes, {"far": P.H_FAR})
    b = classify_boundaries(mesh, wires + tubes, rect)
    basis = Basis(mesh, ElementTriP2())
    A, M = asm(laplace, basis), asm(mass, basis)
    phi_bc = basis.zeros()
    D = []
    for w in wires:
        d = basis.get_dofs(b[w["name"]]).all(); phi_bc[d] = P.V_WIRE; D.append(d)
    for t in tubes:
        d = basis.get_dofs(b[t["name"]]).all(); phi_bc[d] = P.V_TUBE; D.append(d)
    D = np.unique(np.concatenate(D))
    phi_lap = solve(*condense(A, np.zeros(A.shape[0]), x=phi_bc, D=D))
    print(f"mesh {mesh.p.shape[1]} vertices; Laplace solved ({time.time()-t0:.0f} s)")

    # Calibration: thrust is close to proportional to current (T = I d/mu), so use secant steps.
    I = 45.3e-6                       # start near the previous calibration
    log = []
    phi_start = phi_lap
    for k in range(4):
        print(f"calibration step {k}: I_total = {I*1e6:.2f} uA")
        sol = run(I, mesh, b, basis, A, M, D, phi_bc, phi_start)
        dg = diagnostics(sol, mesh, b, basis, phi_lap)
        log.append(dict(I_uA=I * 1e6, T_mN=dg["T"] * 1e3, iterations=len(sol["hist"])))
        print(f"  -> thrust {dg['T']*1e3:.4f} mN (target {P.T_MEAS*1e3:.4f})  ({time.time()-t0:.0f} s)")
        if abs(dg["T"] - P.T_MEAS) < 1e-3 * P.T_MEAS:
            break
        I = I * P.T_MEAS / dg["T"]
        phi_start = sol["phi"]

    V = abs(P.V_WIRE)
    eta = 0.5 * dg["T"] * P.V_JET / (V * dg["I"])
    summary = dict(
        calibration=log,
        total_current_uA=dg["I"] * 1e6,
        current_per_m_of_wire_uA=sol["Iprime"] * 1e6,
        thrust_mN=dg["T"] * 1e3, thrust_target_mN=P.T_MEAS * 1e3,
        effective_axial_length_mm=dg["d_axial_eff_mm"],
        corona_power_W=V * dg["I"],
        ehd_efficiency_pct=100 * eta,
        thrust_per_power_mN_per_W=dg["T"] / (V * dg["I"]) * 1e3,
        charge_conservation_tubes_vs_emitted_pct=100 * (dg["I_tubes_per_m"] / dg["I_emitted_per_m"] - 1),
        max_potential_change_from_space_charge_V=dg["max_dphi_space_charge"],
        max_potential_change_pct_of_V=100 * dg["max_dphi_space_charge"] / V,
        max_ion_density_C_m3=dg["rho_max"],
        rho0_per_wire_C_m3=dg["rho0"],
        thrust_per_wire_mN=(dg["thrust_per_m"] * P.WIRE_LEN_MM * 1e-3 * 1e3).tolist(),
        wire_surface_field_MV_m_laplace_vs_spacecharge=[[round(a, 4), round(c, 4)] for a, c in dg["Es"]],
        outer_iterations=sol["hist"],
        runtime_s=time.time() - t0)
    # Space-charge potential change inside the wire-tube gap (|x| <= 20 mm, 1 <= z <= 10 mm).
    # The largest change sits in the far field upstream and depends on the box size, so the
    # gap values are the physically meaningful ones (see stage2_domain_check.py).
    from matplotlib.tri import LinearTriInterpolator
    nv = mesh.p.shape[1]
    f_dphi = LinearTriInterpolator(Triangulation(mesh.p[0], mesh.p[1], mesh.t.T), (sol["phi"] - phi_lap)[:nv])
    Xg, Zg = np.meshgrid(np.linspace(-20, 20, 81), np.linspace(1, 10, 19))
    dgap = np.ma.filled(f_dphi(Xg, Zg), np.nan)
    summary["gap_region_mean_potential_change_V"] = float(np.nanmean(dgap))
    summary["gap_region_max_potential_change_V"] = float(np.nanmin(dgap))
    summary["note_max_potential_change"] = ("max change occurs in the upstream far field and depends on the box size; "
                                            "gap-region values are the physical ones")
    (RES / "stage2_summary.json").write_text(json.dumps(summary, indent=2))
    print(json.dumps({k: v for k, v in summary.items() if not isinstance(v, list)}, indent=2))
    np.savez(RES / "stage2_solution.npz", p=mesh.p, t=mesh.t, phi=sol["phi"], phi_lap=phi_lap,
             rho_v=sol["rho_v"], origin=sol["origin"])
    plots(mesh, sol, phi_lap)


def plots(mesh, sol, phi_lap):
    tri = Triangulation(mesh.p[0], mesh.p[1], mesh.t.T)
    nv = mesh.p.shape[1]
    rho = -sol["rho_v"] * 1e6                                # micro-C/m^3, magnitude
    fig, axs = plt.subplots(1, 2, figsize=(15, 6))
    lev = np.linspace(0, np.percentile(rho, 99.5), 30)
    cs = axs[0].tricontourf(tri, np.clip(rho, 0, lev[-1]), levels=lev, cmap="magma")
    fig.colorbar(cs, ax=axs[0], label="ion charge density |rho| [uC/m^3]")
    axs[0].set_title("Negative-ion density (space charge)")
    dphi = sol["phi"][:nv] - phi_lap[:nv]
    cs2 = axs[1].tricontourf(tri, dphi, levels=30, cmap="coolwarm")
    fig.colorbar(cs2, ax=axs[1], label="potential change due to space charge [V]")
    axs[1].set_title("Change in potential caused by the ions")
    for ax in axs:
        for x in P.WIRE_X:
            ax.add_patch(Circle((x, P.WIRE_Z), 0.15, color="orange", zorder=5))
        for x in P.TUBE_X:
            ax.add_patch(Circle((x, P.TUBE_Z), P.R_TUBE, color="#b87333", zorder=5))
        ax.set(xlim=(P.X_MIN, P.X_MAX), ylim=(-15, 40), aspect="equal",
               xlabel="x [mm]", ylabel="z [mm] (flow direction)")
    fig.suptitle("Stage 2: space charge at the calibrated current")
    fig.tight_layout()
    fig.savefig(FIG / "stage2_space_charge.png", dpi=150)
    plt.close(fig)


if __name__ == "__main__":
    main()
