"""
EHD (ion-wind) thruster: analytic model vs measurements (V-I-velocity sweep).

NOTE: written before the data were known. Only --demo (synthetic data) runs; the real test
has one operating point and no current reading, so the analysis moved to infer_current.py
and the simulations in sim/. Kept as the template for a future test with a V-I sweep.

Usage:
    python ehd_model.py                 # uses data/geometry.json + data/measurements.csv
    python ehd_model.py --demo          # synthetic data, to check the pipeline runs

All equations are derived in theory.md (section numbers in comments).
SI units everywhere inside the code.
"""
import argparse
import json
import sys
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt

EPS0 = 8.854e-12   # F/m
E0_PEEK = 3.1e6    # V/m
K_PEEK = 0.0308    # m^0.5

ROOT = Path(__file__).parent
DATA = ROOT / "data"
FIGS = ROOT / "figures"


# ---------- physics (theory.md) ----------

def relative_air_density(p_kPa, T_C):
    return (p_kPa / 101.325) * (293.15 / (T_C + 273.15))


def air_density(p_kPa, T_C):
    return p_kPa * 1e3 / (287.05 * (T_C + 273.15))


def peek_onset_voltage(r, d, delta, m=1.0):
    """Sec. 1: corona onset voltage for a wire of radius r at distance d from a plane."""
    E_c = E0_PEEK * m * delta * (1 + K_PEEK / np.sqrt(delta * r))
    return E_c * r * np.log(2 * d / r), E_c


def fit_townsend(V, I):
    """Sec. 2: fit I = C V (V - V0) as a straight line I/V = C V - C V0."""
    slope, intercept = np.polyfit(V, I / V, 1)
    C = slope
    V0 = -intercept / slope
    return C, V0


def thrust_model(I, d, mu):
    """Sec. 3: T = I d / mu."""
    return I * d / mu


def jet_velocity_model(I, d, mu, rho, A):
    """Sec. 4: fully developed jet velocity v_w = sqrt(2 I d / (rho mu A))."""
    return np.sqrt(2 * I * d / (rho * mu * A))


def efficiency_model(v_w, V, d, mu):
    """Sec. 5: eta = d v_w / (2 mu V)."""
    return d * v_w / (2 * mu * V)


def mott_gurney_J(V, d, mu):
    """Sec. 6: collisional space-charge-limited current density."""
    return 9 / 8 * EPS0 * mu * V**2 / d**3


def child_langmuir_J(V, d, q_over_m):
    """Sec. 6: collisionless (vacuum) space-charge-limited current density."""
    return 4 / 9 * EPS0 * np.sqrt(2 * q_over_m) * V**1.5 / d**2


# ---------- data ----------

def load(demo):
    if demo:
        geo = {
            "_note": "SYNTHETIC DEMO GEOMETRY - not the real thruster",
            "wire_diameter_mm": 0.2, "gap_mm": 30.0, "span_mm": 200.0,
            "collector_width_mm": 30.0, "polarity": "positive",
            "pressure_kPa": 91.0, "temperature_C": 27.0,
            "velocity_location": "jet", "mobility_m2_per_Vs": 2.0e-4,
        }
        rng = np.random.default_rng(1)
        V = np.linspace(10e3, 25e3, 7)
        C_true, V0_true = 2.4e-13, 8e3
        I = C_true * V * (V - V0_true) * (1 + 0.05 * rng.standard_normal(V.size))
        rho = geo.get("air_density_moist_kg_m3") or air_density(geo["pressure_kPa"], geo["temperature_C"])
        A = geo["span_mm"] * geo["collector_width_mm"] * 1e-6
        v = jet_velocity_model(I, 0.03, 2.6e-4, rho, A) * (1 + 0.08 * rng.standard_normal(V.size))
        return geo, V, I, v, True

    geo = json.loads((DATA / "geometry.json").read_text())
    rows = np.genfromtxt(DATA / "measurements.csv", delimiter=",", names=True)
    V = np.atleast_1d(rows["voltage_kV"]) * 1e3
    I = np.atleast_1d(rows["current_uA"]) * 1e-6
    v = np.atleast_1d(rows["velocity_mps"])
    return geo, V, I, v, False


# ---------- analysis ----------

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--demo", action="store_true")
    args = ap.parse_args()

    geo, V, I, v_meas, synthetic = load(args.demo)
    tag = " [SYNTHETIC DEMO]" if synthetic else ""
    FIGS.mkdir(exist_ok=True)

    r = geo["wire_diameter_mm"] * 1e-3 / 2
    d = geo["gap_mm"] * 1e-3
    if "opening_diameter_mm" in geo:          # circular duct (bladeless-fan frame)
        A = np.pi * (geo["opening_diameter_mm"] * 1e-3) ** 2 / 4
    else:
        A = (geo["flow_area_mm2"] if "flow_area_mm2" in geo else geo["span_mm"] * geo["collector_width_mm"]) * 1e-6
    mu = geo["mobility_m2_per_Vs"]
    rho = geo.get("air_density_moist_kg_m3") or air_density(geo["pressure_kPa"], geo["temperature_C"])
    delta = relative_air_density(geo["pressure_kPa"], geo["temperature_C"])

    # Sec. 4: convert measured velocity to fully developed jet velocity
    v_w = 2 * v_meas if geo["velocity_location"] == "collector" else v_meas

    V0_peek, E_c = peek_onset_voltage(r, d, delta, geo.get("surface_factor_m", 1.0))
    lines = [f"EHD thruster analysis{tag}", "",
             f"air density rho = {rho:.3f} kg/m3, delta = {delta:.3f}",
             f"Peek critical field E_c = {E_c/1e6:.2f} MV/m, onset V0 = {V0_peek/1e3:.2f} kV"]

    if V.size >= 3:
        C, V0_fit = fit_townsend(V, I)
        lines.append(f"Townsend fit: C = {C:.3e} A/V^2, V0 = {V0_fit/1e3:.2f} kV "
                     f"(Peek predicts {V0_peek/1e3:.2f} kV)")
    else:
        C = V0_fit = None
        lines.append("Fewer than 3 voltage points: Townsend fit skipped.")

    T_model = thrust_model(I, d, mu)
    v_model = jet_velocity_model(I, d, mu, rho, A)
    T_from_v = 0.5 * rho * A * v_w**2                  # Sec. 4, jet form
    P_in = V * I
    eta_meas = 0.25 * rho * A * v_w**3 / P_in           # jet kinetic power / input
    eta_model = efficiency_model(v_model, V, d, mu)
    mu_eff = 2 * I * d / (rho * A * v_w**2)             # mobility that would fit each point

    lines += ["", f"{'V[kV]':>7} {'I[uA]':>8} {'v_w meas':>9} {'v_w model':>9} "
              f"{'T_v[mN]':>8} {'T_mod[mN]':>9} {'eta_meas%':>9} {'eta_mod%':>8} {'mu_eff':>9}"]
    for k in range(V.size):
        lines.append(f"{V[k]/1e3:7.2f} {I[k]*1e6:8.1f} {v_w[k]:9.2f} {v_model[k]:9.2f} "
                     f"{T_from_v[k]*1e3:8.2f} {T_model[k]*1e3:9.2f} "
                     f"{eta_meas[k]*100:9.2f} {eta_model[k]*100:8.2f} {mu_eff[k]:9.2e}")

    J_mg = mott_gurney_J(V, d, mu) * A
    lines += ["", "Mott-Gurney collisional limit current (whole collector area):",
              "  " + ", ".join(f"{x*1e6:.0f} uA" for x in J_mg)]
    summary = "\n".join(lines)
    print(summary)
    (FIGS / "summary.txt").write_text(summary + "\n")

    # ---- figures ----
    fig, ax = plt.subplots(1, 3, figsize=(15, 4.5))

    ax[0].plot(V / 1e3, I / V * 1e9, "o", label="measured")
    if C is not None:
        Vf = np.linspace(V.min(), V.max(), 100)
        ax[0].plot(Vf / 1e3, C * (Vf - V0_fit) * 1e9, "-", label=f"fit, V0={V0_fit/1e3:.1f} kV")
    ax[0].axvline(V0_peek / 1e3, ls="--", c="gray", label="Peek onset")
    ax[0].set(xlabel="V [kV]", ylabel="I/V [nA/V]", title="Townsend check (Sec. 2)")
    ax[0].legend()

    sI = np.sqrt(I)
    ax[1].plot(sI * 1e3, v_w, "o", label="measured")  # sqrt(A)*1e3 = sqrt(uA)
    sIf = np.linspace(0, sI.max() * 1.1, 100)
    for m_, ls in [(1.4e-4, ":"), (mu, "-"), (2.7e-4, ":")]:
        ax[1].plot(sIf * 1e3, jet_velocity_model(sIf**2, d, m_, rho, A), ls,
                   label=f"model, mu={m_:.1e}")
    ax[1].set(xlabel="sqrt(I) [sqrt(uA)]", ylabel="jet velocity [m/s]",
              title="v vs sqrt(I) (Sec. 4)")
    ax[1].legend()

    ax[2].plot(V / 1e3, eta_meas * 100, "o", label="measured")
    ax[2].plot(V / 1e3, eta_model * 100, "s", label="model")
    ax[2].set(xlabel="V [kV]", ylabel="efficiency [%]", title="Efficiency (Sec. 5)")
    ax[2].legend()

    fig.suptitle("EHD thruster: model vs measurement" + tag)
    fig.tight_layout()
    fig.savefig(FIGS / ("demo_analysis.png" if synthetic else "analysis.png"), dpi=150)


if __name__ == "__main__":
    sys.exit(main())
