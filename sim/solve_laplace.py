"""
Stage 1: electrostatic field of the real thruster geometry, no space charge.

11 wires at V_WIRE face 7 grounded tubes. Solves Laplace's equation, then reports for each
wire: the surface field (mean, max, where the max sits) and the charge per metre.

Why it matters: corona starts where the surface field reaches the onset field
(~14.6 MV/m for a smooth 36 AWG wire here). Comparing each wire's field with that value
shows which wires can emit and how far below onset a smooth wire sits at the supply voltage.
"""
import csv
import json
from pathlib import Path

import numpy as np

import params as P
from fem_tools import build_mesh, classify_boundaries, solve_laplace, surface_field

HERE = Path(__file__).parent
OUT = HERE / "results"
OUT.mkdir(exist_ok=True)


def thruster_circles():
    wires = [dict(name=f"wire{i}", x=float(x), z=P.WIRE_Z, r=P.R_WIRE, h=P.H_WIRE,
                  grow=P.WIRE_GROWTH_DIST) for i, x in enumerate(P.WIRE_X)]
    tubes = [dict(name=f"tube{j}", x=float(x), z=P.TUBE_Z, r=P.R_TUBE, h=P.H_TUBE,
                  grow=P.TUBE_GROWTH_DIST) for j, x in enumerate(P.TUBE_X)]
    return wires, tubes


def main():
    rect = (P.X_MIN, P.Z_MIN, P.X_MAX, P.Z_MAX)
    wires, tubes = thruster_circles()
    mesh = build_mesh(rect, wires + tubes, {"far": P.H_FAR})
    b = classify_boundaries(mesh, wires + tubes, rect)
    print(f"mesh: {mesh.p.shape[1]} nodes, {mesh.t.shape[1]} triangles")

    dirichlet = [(b[w["name"]], P.V_WIRE) for w in wires] + [(b[t["name"]], P.V_TUBE) for t in tubes]
    basis, phi = solve_laplace(mesh, dirichlet)

    # Analytic reference: an infinite row of wires (pitch s) in front of a flat plate at distance d.
    s, d, r = P.WIRE_PITCH, P.GAP, P.R_WIRE
    E_array = abs(P.V_WIRE) / (r * 1e-3 * (2 * np.pi * d / s + np.log(s / (2 * np.pi * r))))

    rows, q_wires = [], 0.0
    for i, w in enumerate(wires):
        sf = surface_field(mesh, phi, b[w["name"]], (w["x"], w["z"]))
        nearest = P.TUBE_X[np.argmin(np.abs(P.TUBE_X - w["x"]))]
        offset = w["x"] - nearest
        rows.append(dict(wire=i, x_mm=round(w["x"], 2), nearest_tube_x_mm=float(nearest),
                         lateral_offset_mm=round(offset, 2),
                         faces="tube" if abs(offset) < 1.5 else "gap",
                         E_mean_MV_m=round(sf["E_mean"] / 1e6, 3),
                         E_max_MV_m=round(sf["E_max"] / 1e6, 3),
                         angle_of_max_deg=round(sf["angle_max_deg"], 1),
                         E_max_over_onset=round(sf["E_max"] / P.E_ONSET_SMOOTH, 3),
                         charge_nC_per_m=round(sf["charge_per_m"] * 1e9, 3)))
        q_wires += sf["charge_per_m"]
    q_tubes = sum(surface_field(mesh, phi, b[t["name"]], (t["x"], t["z"]))["charge_per_m"] for t in tubes)

    with open(OUT / "stage1_wire_fields.csv", "w", newline="") as f:
        wr = csv.DictWriter(f, fieldnames=rows[0].keys())
        wr.writeheader()
        wr.writerows(rows)

    Emax = np.array([r_["E_max_MV_m"] for r_ in rows])
    summary = dict(V_wire=P.V_WIRE, nodes=int(mesh.p.shape[1]), triangles=int(mesh.t.shape[1]),
                   E_max_range_MV_m=[float(Emax.min()), float(Emax.max())],
                   E_mean_all_MV_m=float(np.mean([r_["E_mean_MV_m"] for r_ in rows])),
                   analytic_wire_row_vs_plate_MV_m=float(E_array / 1e6),
                   onset_smooth_MV_m=P.E_ONSET_SMOOTH / 1e6,
                   voltage_for_onset_on_hottest_wire_kV=float(abs(P.V_WIRE) * P.E_ONSET_SMOOTH / (Emax.max() * 1e6) / 1e3),
                   charge_wires_nC_m=q_wires * 1e9, charge_tubes_nC_m=q_tubes * 1e9)
    (OUT / "stage1_summary.json").write_text(json.dumps(summary, indent=2))

    print(f"{'wire':>4} {'x':>6} {'faces':>5} {'Emean':>7} {'Emax':>7} {'angle':>6} {'Emax/onset':>10} {'q nC/m':>8}")
    for r_ in rows:
        print(f"{r_['wire']:>4} {r_['x_mm']:6.1f} {r_['faces']:>5} {r_['E_mean_MV_m']:7.2f} {r_['E_max_MV_m']:7.2f} "
              f"{r_['angle_of_max_deg']:6.0f} {r_['E_max_over_onset']:10.2f} {r_['charge_nC_per_m']:8.2f}")
    print(json.dumps(summary, indent=2))

    np.savez(OUT / "stage1_solution.npz", p=mesh.p, t=mesh.t, phi=phi)


if __name__ == "__main__":
    main()
