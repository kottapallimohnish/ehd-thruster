"""
Stage 1b: ion paths and their lengths.

Without space charge, a negative ion leaving a wire drifts along the field line
(direction +grad(phi), i.e. opposite to E, because its charge is negative).
The force density on the air is f = j/mu (a vector). Integrated along one ion path carrying
current I_k, its axial (z) component is (I_k/mu) * (z_end - z_start): sideways detours cancel
by symmetry, and an ion that first drifts upstream pushes the air backwards for that part.
So the length that sets the THRUST is the net axial displacement of the ions, d_axial,
used in T = I * d_axial / mu. The full arc length is reported too (it sets the ion transit
time, which matters for space charge in stage 2), but it does not set the thrust.

Weighting: lines are seeded uniformly around each wire. Emission is concentrated
where the surface field is highest, so a field-weighted mean is reported as well.
"""
import csv
import json
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
from skfem import MeshTri, Basis, ElementTriP2

import params as P

HERE = Path(__file__).parent
RES = HERE / "results"
FIG = HERE / "figures"
FIG.mkdir(exist_ok=True)

N_SEED = 24          # field lines per wire
MAX_LEN = 400.0      # mm; give up beyond this
MAX_STEPS = 20000


def load():
    d = np.load(RES / "stage1_solution.npz")
    mesh = MeshTri(d["p"], d["t"])
    basis = Basis(mesh, ElementTriP2())
    phi = d["phi"]
    grad = basis.interpolate(phi).grad
    gx, gz = basis.project(grad[0]), basis.project(grad[1])   # L2 projection to P2, smooth field
    return mesh, basis, phi, gx, gz


def dist_to_wires(pts):
    dx = pts[0][:, None] - P.WIRE_X[None, :]
    dz = pts[1][:, None] - P.WIRE_Z
    return np.hypot(dx, dz).min(axis=1) - P.R_WIRE


def in_tube(pts):
    dx = pts[0][:, None] - P.TUBE_X[None, :]
    dz = pts[1][:, None] - P.TUBE_Z
    return np.hypot(dx, dz).min(axis=1) <= P.R_TUBE * 1.003


def outside(pts):
    return (pts[0] <= P.X_MIN + 1e-3) | (pts[0] >= P.X_MAX - 1e-3) | \
           (pts[1] <= P.Z_MIN + 1e-3) | (pts[1] >= P.Z_MAX - 1e-3)


def direction(basis, gx, gz, pts):
    """Unit vector along +grad(phi) (direction of negative-ion drift) and |grad phi|."""
    M = basis.probes(pts)
    vx, vz = M @ gx, M @ gz
    n = np.hypot(vx, vz) + 1e-300
    return np.vstack([vx / n, vz / n]), n


def trace(basis, gx, gz, start):
    """RK2 (midpoint) integration of all lines at once with a step that grows away from the wires."""
    pts = start.copy()
    paths = [[p.copy()] for p in start.T]
    length = np.zeros(start.shape[1])
    state = np.zeros(start.shape[1], int)       # 0 running, 1 reached tube, 2 left domain/too long
    for _ in range(MAX_STEPS):
        run = np.where(state == 0)[0]
        if run.size == 0:
            break
        p = pts[:, run]
        h = np.clip(0.15 * (dist_to_wires(p) + P.R_WIRE), 0.004, 0.25)
        d1, _ = direction(basis, gx, gz, p)
        pm = p + 0.5 * h * d1
        pm[0] = np.clip(pm[0], P.X_MIN + 1e-3, P.X_MAX - 1e-3)
        pm[1] = np.clip(pm[1], P.Z_MIN + 1e-3, P.Z_MAX - 1e-3)
        # A half-step can already land inside a tube (a hole in the mesh): stop those lines there.
        stop = in_tube(pm) | outside(pm) | (dist_to_wires(pm) < 0)
        pn = pm.copy()
        go = ~stop
        if go.any():
            d2, _ = direction(basis, gx, gz, pm[:, go])
            pn[:, go] = p[:, go] + h[go] * d2
        h = np.where(stop, 0.5 * h, h)
        bad = dist_to_wires(pn) < 0                  # never step into a wire hole
        state[run[bad]] = 2
        pts[:, run] = pn
        length[run] += h
        for k, i in enumerate(run):
            paths[i].append(pn[:, k].copy())
        state[run[in_tube(pn) & (state[run] == 0)]] = 1
        state[run[(outside(pn) | (length[run] > MAX_LEN)) & (state[run] == 0)]] = 2
    return np.array(paths, dtype=object), length, state


def main():
    mesh, basis, phi, gx, gz = load()
    wf = {int(r["wire"]): r for r in csv.DictReader(open(RES / "stage1_wire_fields.csv"))}

    # Seed points just outside each wire surface; also record the local surface field there.
    th = np.linspace(0, 2 * np.pi, N_SEED, endpoint=False)
    seeds, owner, w_field = [], [], []
    for i, xw in enumerate(P.WIRE_X):
        rs = P.R_WIRE * 1.05
        sx, sz = xw + rs * np.sin(th), P.WIRE_Z + rs * np.cos(th)   # angle 0 = +z (downstream)
        seeds.append(np.vstack([sx, sz]))
        owner += [i] * N_SEED
        _, gmag = direction(basis, gx, gz, np.vstack([sx, sz]))
        w_field.append(gmag * 1e3)                                  # V/m near the surface
    seeds = np.hstack(seeds)
    owner = np.array(owner)
    w_field = np.concatenate(w_field)

    paths, length, state = trace(basis, gx, gz, seeds)

    z_end = np.array([np.array(pth)[-1][1] for pth in paths])
    dz = z_end - P.WIRE_Z                       # net axial displacement of each ion path
    rows = []
    for i in range(P.N_WIRES):
        sel = (owner == i) & (state == 1)
        L, W, Z = length[sel], w_field[sel], dz[sel]
        front = sel & (np.cos(np.repeat(th[None, :], P.N_WIRES, 0).ravel()) > 0)
        rows.append(dict(wire=i, x_mm=round(float(P.WIRE_X[i]), 2), faces=wf[i]["faces"],
                         lines_reaching_tube=int(sel.sum()), of=N_SEED,
                         axial_disp_mm=round(float((Z * W).sum() / W.sum()), 2),
                         mean_len_mm=round(float(L.mean()), 2),
                         field_weighted_len_mm=round(float((L * W).sum() / W.sum()), 2),
                         front_half_len_mm=round(float(length[front].mean()), 2),
                         min_len_mm=round(float(L.min()), 2), max_len_mm=round(float(L.max()), 2)))
    reached = state == 1
    summary = dict(straight_gap_mm=P.GAP,
                   lines=int(seeds.shape[1]), reached_tube=int(reached.sum()),
                   axial_disp_field_weighted_mm=float((dz[reached] * w_field[reached]).sum() / w_field[reached].sum()),
                   axial_disp_over_gap=float((dz[reached] * w_field[reached]).sum() / w_field[reached].sum() / P.GAP),
                   mean_len_all_mm=float(length[reached].mean()),
                   field_weighted_len_all_mm=float((length[reached] * w_field[reached]).sum() / w_field[reached].sum()),
                   ratio_field_weighted_to_gap=float((length[reached] * w_field[reached]).sum()
                                                     / w_field[reached].sum() / P.GAP))
    with open(RES / "stage1_fieldline_lengths.csv", "w", newline="") as f:
        wr = csv.DictWriter(f, fieldnames=rows[0].keys())
        wr.writeheader()
        wr.writerows(rows)
    (RES / "stage1_fieldline_summary.json").write_text(json.dumps(summary, indent=2))
    for r in rows:
        print(r)
    print(json.dumps(summary, indent=2))

    plot(mesh, phi, paths, state, owner)


def plot(mesh, phi, paths, state, owner):
    from matplotlib.tri import Triangulation
    tri = Triangulation(mesh.p[0], mesh.p[1], mesh.t.T)
    pv = phi[: mesh.p.shape[1]] / 1e3                      # P2 vertex dofs come first (kV)

    for name, xl, zl, figsize in [("stage1_fieldlines_full.png", (P.X_MIN, P.X_MAX), (-20, 40), (11, 8)),
                                  ("stage1_fieldlines_zoom.png", (-14, 14), (-6, 18), (9, 8))]:
        fig, ax = plt.subplots(figsize=figsize)
        cs = ax.tricontourf(tri, pv, levels=np.linspace(P.V_WIRE / 1e3, 0, 23), cmap="viridis")
        fig.colorbar(cs, ax=ax, label="potential [kV]")
        for path, st in zip(paths, state):
            q = np.array(path)
            ax.plot(q[:, 0], q[:, 1], color="w" if st == 1 else "r", lw=0.5, alpha=0.8)
        for x in P.WIRE_X:
            ax.add_patch(Circle((x, P.WIRE_Z), max(P.R_WIRE, 0.15), color="orange", zorder=5))
        for x in P.TUBE_X:
            ax.add_patch(Circle((x, P.TUBE_Z), P.R_TUBE, color="#b87333", zorder=5))
        ax.set(xlim=xl, ylim=zl, xlabel="x [mm] (across the opening)", ylabel="z [mm] (flow direction)",
               aspect="equal",
               title="Stage 1: potential and ion paths (no space charge), wires at "
                     f"{P.V_WIRE/1e3:.0f} kV\nwhite = reaches a tube; wires drawn enlarged")
        fig.tight_layout()
        fig.savefig(FIG / name, dpi=160)
        plt.close(fig)

    # surface field per wire
    rows = list(csv.DictReader(open(RES / "stage1_wire_fields.csv")))
    x = [float(r["x_mm"]) for r in rows]
    fig, ax = plt.subplots(figsize=(8, 4.5))
    ax.bar(x, [float(r["E_max_MV_m"]) for r in rows], width=2.5, color=["#b87333" if r["faces"] == "tube"
                                                                         else "#5a7" for r in rows])
    ax.axhline(P.E_ONSET_SMOOTH / 1e6, color="k", ls="--", label="smooth-wire onset 14.6 MV/m")
    ax.set(xlabel="wire position x [mm]", ylabel="max surface field [MV/m]",
           title=f"Surface field per wire at {abs(P.V_WIRE)/1e3:.0f} kV (copper = faces a tube, green = faces a gap)")
    ax.legend()
    fig.tight_layout()
    fig.savefig(FIG / "stage1_surface_field_per_wire.png", dpi=160)
    plt.close(fig)


if __name__ == "__main__":
    main()
