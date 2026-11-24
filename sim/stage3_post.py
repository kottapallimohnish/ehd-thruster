"""
Stage 3 post-processing: time-averaged flow, velocity at the anemometer, momentum budget.

Key comparison. The stage 2 calibration assumed an IDEAL actuator disk: all the ion force
ends up as jet momentum, with no losses on the tubes. The CFD includes the tubes, frames
and jet mixing. So:
    r = (CFD velocity at the anemometer) / (ideal 2D actuator-disk velocity for the same force)
measures how much of the ideal jet survives to the anemometer. For the measured 1.21 m/s
to be reproduced, the ion force (and current) must be 1/r^2 times the ideal-disk value,
because jet speed scales with the square root of the force.
"""
import json
import re
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, Rectangle
from matplotlib.tri import Triangulation, LinearTriInterpolator

import params as P
import stage3_build_case as B

HERE = Path(__file__).parent
CASE = B.CASE
RES, FIG = HERE / "results", HERE / "figures"


def read_U(tdir):
    txt = (tdir / "U").read_text()
    i = txt.index("internalField")
    j = txt.index("(", txt.index("List<vector>", i))
    body = txt[j + 1: txt.index("\n)\n", j)]
    rows = re.findall(r"\(([^()]+)\)", body)
    return np.array([list(map(float, r.split())) for r in rows])


def main():
    C, V = B.cell_geometry()                      # mm, m^3
    area = V / (B.DEPTH * 1e-3)                   # m^2
    times = sorted([d for d in CASE.iterdir() if re.fullmatch(r"[0-9.e-]+", d.name) and float(d.name) > 0],
                   key=lambda d: float(d.name))
    use = [d for d in times if float(d.name) >= B.AVG_START - 1e-9]
    Us = np.array([read_U(d) for d in use])        # (nt, ncell, 3)
    Umean = Us.mean(axis=0)
    x, y = C[:, 0], C[:, 1]

    def band(y0, half=2.0, xhalf=B.VANE_HALF):
        return (np.abs(y - y0) <= half) & (np.abs(x) <= xhalf)

    sel = band(B.ANEMO_Y)
    w = area[sel]
    v_anemo_t = np.array([(U[sel, 1] * w).sum() / w.sum() for U in Us])
    v_anemo = float(v_anemo_t.mean())
    cl = band(B.ANEMO_Y, xhalf=1.0)
    v_centre = float((Umean[cl, 1] * area[cl]).sum() / area[cl].sum())

    # slot speed between tubes at the tube-centre plane
    slots = np.zeros_like(x, bool)
    for xm in (np.r_[P.TUBE_X[:-1]] + 4.5):
        slots |= (np.abs(x - xm) <= 1.0) & (np.abs(y - P.TUBE_Z) <= 0.5)
    v_slot = float((Umean[slots, 1] * area[slots]).sum() / area[slots].sum())

    # momentum flux of the mean jet across horizontal lines (per metre depth)
    T_force = json.loads((CASE / "force_mapping.json").read_text())["thrust_per_m_in_zones"]
    def mom_flux(y0, half=1.5):
        """rho * integral(Uy^2) dx across the line y = y0, from cells in a thin band (N per metre depth)."""
        s = (np.abs(y - y0) <= half) & (np.abs(x) <= 140)
        return float((P.RHO_AIR * Umean[s, 1] ** 2 * area[s]).sum() / (2 * half * 1e-3))
    flux = {f"{yy:.0f}mm": mom_flux(yy) for yy in (40.0, B.ANEMO_Y, 120.0)}

    v_ideal_2d = float(np.sqrt(2 * T_force / (P.RHO_AIR * 2 * B.FRAME_INNER * 1e-3)))
    r = v_anemo / v_ideal_2d
    I_ideal = json.loads((RES / "stage2_summary.json").read_text())["total_current_uA"] * 1e-6
    I_corr = I_ideal / r**2
    # time statistics of the (unsteady) anemometer signal
    tt = np.array([float(d.name) for d in use])
    slope = float(np.polyfit(tt, v_anemo_t, 1)[0])
    batches = np.array_split(v_anemo_t, 6)
    batch_se = float(np.std([b.mean() for b in batches], ddof=1) / np.sqrt(len(batches)))
    out = dict(
        averaging_window_s=[float(use[0].name), float(use[-1].name)], samples=len(use),
        force_per_m_N=T_force,
        ideal_2d_actuator_velocity_mps=v_ideal_2d,
        cfd_velocity_at_anemometer_mps=v_anemo,
        cfd_velocity_at_anemometer_std_over_time_mps=float(v_anemo_t.std()),
        cfd_centreline_velocity_at_anemometer_mps=v_centre,
        cfd_mean_slot_velocity_mps=v_slot,
        survival_ratio_r=r,
        jet_momentum_flux_N_per_m=flux,
        momentum_flux_fraction_of_force_at_anemometer=flux[f"{B.ANEMO_Y:.0f}mm"] / T_force,
        current_if_only_this_run_r_uA=I_corr * 1e6,
        efficiency_if_only_this_run_r_pct=100 * (0.5 * P.T_MEAS * P.V_JET) / (abs(P.V_WIRE) * I_corr),   # measured jet power / electrical power
        anemometer_trend_mps_per_s=slope,
        anemometer_mean_batch_se_mps=batch_se,
        anemometer_thirds_mps=[float(b.mean()) for b in np.array_split(v_anemo_t, 3)],
        note=("current_if_only_this_run_r = stage-2 current / r^2 using this single run's r; the headline "
              "current comes from sensitivity.py, which combines meshes, force effect, 2D->3D and input ranges. "
              "batch_se: standard error of the mean from 6 consecutive batch means (accounts for autocorrelation)."),
    )
    tag = f"_{B.VARIANT}" if B.VARIANT else ""
    (RES / f"stage3_summary{tag}.json").write_text(json.dumps(out, indent=2))
    print(json.dumps(out, indent=2))
    np.save(RES / f"stage3_anemometer_timeseries{tag}.npy", np.c_[[float(d.name) for d in use], v_anemo_t])
    np.savez_compressed(RES / f"stage3_mean_flow{tag}.npz", cell_centres_mm=C, cell_volume_m3=V, U_mean=Umean,
                        window_s=np.array([float(use[0].name), float(use[-1].name)]))
    if not B.VARIANT:
        plot(C, Umean, Us[-1])


def plot(C, Umean, Uinst):
    xy = C[:, :2]
    tri = Triangulation(xy[:, 0], xy[:, 1])
    # remove triangles spanning solids (tubes/frames): long edges
    pts = xy[tri.triangles]
    edge = np.max(np.linalg.norm(pts - np.roll(pts, 1, axis=1), axis=2), axis=1)
    cen = pts.mean(1)
    inside = np.zeros(len(cen), bool)
    for xt in P.TUBE_X:
        inside |= np.hypot(cen[:, 0] - xt, cen[:, 1] - P.TUBE_Z) < P.R_TUBE
    for ya, yb in (B.WIRE_FRAME_Y, B.TUBE_FRAME_Y):
        inside |= (np.abs(cen[:, 0]) > B.FRAME_INNER) & (np.abs(cen[:, 0]) < B.FRAME_OUTER) & (cen[:, 1] > ya) & (cen[:, 1] < yb)
    tri.set_mask(inside | (edge > 10))
    sp = np.hypot(Umean[:, 0], Umean[:, 1])

    fig, axs = plt.subplots(1, 2, figsize=(15, 7.5))
    for ax, field, title in ((axs[0], sp, "time-averaged speed [m/s]"),
                             (axs[1], np.hypot(Uinst[:, 0], Uinst[:, 1]), "instantaneous speed [m/s] (t = end)")):
        cs = ax.tricontourf(tri, field, levels=np.linspace(0, max(3.5, np.percentile(field, 99.8)), 36), cmap="turbo")
        fig.colorbar(cs, ax=ax, shrink=0.8, label="m/s")
        for xt in P.TUBE_X:
            ax.add_patch(Circle((xt, P.TUBE_Z), P.R_TUBE, color="#b87333"))
        for ya, yb in (B.WIRE_FRAME_Y, B.TUBE_FRAME_Y):
            for s in (-1, 1):
                xa = B.FRAME_INNER if s > 0 else -B.FRAME_OUTER
                ax.add_patch(Rectangle((xa, ya), B.FRAME_OUTER - B.FRAME_INNER, yb - ya, color="0.6"))
        ax.plot(P.WIRE_X, np.zeros_like(P.WIRE_X), ".", color="orange", ms=3)
        ax.add_patch(Rectangle((-B.VANE_HALF, B.ANEMO_Y - 2), 2 * B.VANE_HALF, 4, fill=False, ec="w", lw=1.5))
        ax.text(B.VANE_HALF + 2, B.ANEMO_Y, "anemometer", color="w", va="center", fontsize=9)
        ax.set(xlim=(-70, 70), ylim=(-30, 130), aspect="equal", xlabel="x [mm]", ylabel="flow direction [mm]", title=title)
    fig.suptitle("Stage 3: airflow driven by the ion force (2D, laminar, transient)")
    fig.tight_layout()
    fig.savefig(FIG / "stage3_flow.png", dpi=140)
    plt.close(fig)

    ts = np.load(RES / "stage3_anemometer_timeseries.npy")
    fig, ax = plt.subplots(figsize=(8, 3.5))
    ax.plot(ts[:, 0], ts[:, 1], "o-", label="CFD, vane-averaged")
    ax.axhline(P.V_JET, color="k", ls="--", label="measured 1.21 m/s")
    ax.axhline(ts[:, 1].mean(), color="C0", ls=":", label=f"CFD time mean {ts[:, 1].mean():.2f} m/s")
    ax.set(xlabel="time [s]", ylabel="velocity at anemometer [m/s]",
           title="Velocity 60 mm behind the tubes, at the ideal-jet (stage 2) force;\nthe shortfall from 1.21 m/s is the flow loss that raises the inferred current")
    ax.legend()
    fig.tight_layout()
    fig.savefig(FIG / "stage3_anemometer_timeseries.png", dpi=140)
    plt.close(fig)


if __name__ == "__main__":
    main()
