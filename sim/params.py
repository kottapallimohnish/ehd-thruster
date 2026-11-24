"""
All simulation parameters in one place.

Length unit inside the mesh and solver: millimetres (mm).
Potential: volts. So grad(phi) comes out in V/mm; multiply by 1e3 for V/m.
"""
import numpy as np

# --- emitter: 11 bare 36 AWG copper wires ---------------------------------
R_WIRE = 0.0635                     # mm (0.127 mm diameter)
N_WIRES = 11
WIRE_PITCH = 4.3                    # mm, idealised uniform pitch
WIRE_X = -21.5 + WIRE_PITCH * np.arange(N_WIRES)   # mm
WIRE_Z = 0.0                        # wire plane

# --- collector: 7 copper tubes --------------------------------------------
R_TUBE = 3.0                        # mm (6 mm OD)
TUBE_X = np.array([-27.0, -18.0, -9.0, 0.0, 9.0, 18.0, 27.0])
TUBE_Z = 14.11                      # tube centres; front surface at 11.11 mm
GAP = TUBE_Z - R_TUBE               # 11.11 mm, wire plane -> tube front surface

# --- air domain (2D slice across the wires) -------------------------------
X_MIN, X_MAX = -40.0, 40.0
Z_MIN, Z_MAX = -20.0, 80.0

# --- electrical ----------------------------------------------------------
V_WIRE = -11_000.0                  # V, nominal peak supply (spark-gap estimate)
V_TUBE = 0.0
EPS0 = 8.854e-12                    # F/m

# --- reference values (for sanity checks, not inputs) ---------------------
E_ONSET_SMOOTH = 14.6e6             # V/m, Peek onset field for this wire (smooth, delta=0.944)

# --- mesh sizes (mm) ------------------------------------------------------
H_WIRE = 0.006                      # on wire surfaces (~66 segments per wire)
H_TUBE = 0.2                        # on tube surfaces
H_FAR = 2.0                         # far field
WIRE_GROWTH_DIST = 10.0             # distance over which size grows from H_WIRE to H_FAR
TUBE_GROWTH_DIST = 8.0

# --- stage 2: ions, air and the measurement -------------------------------
MU_ION = 2.0e-4                     # m^2/(V s), O2- in humid air (fit-range 1.8-2.2e-4)
WIRE_LEN_MM = np.array([48.728, 55.066, 59.082, 62.349, 64.363, 65.0,
                        64.363, 62.349, 59.082, 55.066, 48.728])   # active (nicked) length per wire
RHO_AIR = 1.124                     # kg/m^3, moist air, 15 Aug 2023 12-15 h
A_FLOW_MM2 = 3318.0                 # 65 mm opening
V_JET = 1.21                        # m/s, mean steady anemometer reading
T_MEAS = 0.5 * RHO_AIR * A_FLOW_MM2 * 1e-6 * V_JET**2   # N, actuator-disk thrust from the jet speed
