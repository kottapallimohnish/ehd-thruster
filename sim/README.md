# Stage 1 simulation: electrostatics (no space charge)

2D finite-element model of the real electrode layout: 11 wires (36 AWG, r = 0.0635 mm,
4.3 mm pitch) at −11 kV facing 7 grounded tubes (6 mm OD, 9 mm pitch, front surface 11.11 mm away).
gmsh mesh, scikit-fem solver with quadratic (P2) elements.

## Run
    pip install gmsh scikit-fem meshio matplotlib numpy     # Linux also needs libglu1-mesa
    python benchmark_wire_plane.py   # validation against the exact wire-above-plate solution
    python solve_laplace.py          # field on each wire  -> results/stage1_wire_fields.csv
    python fieldlines.py             # ion paths, lengths  -> results/, figures/

## Files
| file | what it does |
|---|---|
| `params.py` | every number in the model, in one place |
| `fem_tools.py` | mesh (gmsh), Laplace solve (P2 FEM), surface field and charge (Gauss's law) |
| `benchmark_wire_plane.py` | validation case with an exact answer |
| `solve_laplace.py` | thruster geometry, per-wire surface field and charge |
| `fieldlines.py` | traces negative-ion paths (along +grad φ) from each wire to the tubes |

## Validation
Wire above a flat plate, gap 11.11 mm, exact mean surface field V / (r·acosh(h/r)) = 29.57 MV/m:

| wire element size | elements | error in mean field | error in charge |
|---|---|---|---|
| 0.012 mm | 6 018 | −1.54 % | −1.68 % |
| **0.006 mm (used)** | 6 574 | **−0.67 %** | −0.71 % |
| 0.003 mm | 7 473 | −0.31 % | −0.32 % |

Error halves with element size, as expected. Main model at 0.003 mm changes wire fields by < 0.4 %.
Charge on wires (−431.6 nC/m) and tubes (+432.6 nC/m) balance to 0.2 % (Gauss's law check).

## Results at −11 kV
- **Surface field per wire:** end wires 15.9 MV/m, inner wires 10.3–12.2 MV/m.
  Smooth-wire onset is 14.6 MV/m, so **only the two outermost wires exceed onset**
  (the hottest wire reaches onset at 10.1 kV).
- **Neighbour shielding dominates, not tube alignment.** A wire facing a tube head-on and a wire
  facing a gap differ by < 1 % (10.33 vs 10.38 MV/m). The field is set by how many neighbours
  screen a wire: the end wires have neighbours on one side only.
- **A single isolated wire** at the same voltage and gap would have 29.6 MV/m (benchmark):
  the array cuts the surface field by about a factor of 3.
- **Ion paths:** net axial displacement (what sets thrust) = **12.1 mm**, 9 % more than the
  straight gap, because many ions land on the sides of the tubes. Mean arc length = 17.5 mm
  (ions from the back of the inner wires loop round their neighbours).

## What this means for the experiment
1. Consistent with the very faint hiss with smooth wires: at ~11 kV only the two end wires
   slightly exceed the smooth-wire onset. Wire ends and leads, not modelled, could also contribute.
2. The inner 7 wires need ≥ 14.5 kV at nominal inputs (≥ 11.8 kV with gap 9.86 mm and Peek m = 0.9),
   so smooth wires were at best marginal at the 9–14 kV supply; the nicks (local field enhancement,
   not modelled) are the likely reason flow appeared only after nicking.
3. Thrust model correction: use d_axial = 12.1 mm instead of 11.11 mm in T = I·d/μ.

## Limitations
- 2D: wires are treated as infinitely long; the real wires are chords of different lengths
  in a 65 mm circle, and their ends are not modelled.
- Printed frames, M8 rods and the module leads are not modelled. The air boundary at
  x = ±40, z = −20 mm has zero normal field. This mainly affects the back-emitted ion paths
  (one line from the centre wire leaves the domain), not the axial displacement.
- No space charge yet (stage 2): once ions fill the gap they shield the wires further and
  bend the field lines.
- Nicks are not resolved; their effect enters stage 2 through a calibrated onset field.

---

# Stage 2: space charge (`stage2_space_charge.py`)

## Model
- Poisson's equation −∇²φ = ρ/ε₀, coupled to ion charge conservation ∇·(ρμE) = 0
  (negative ions, μ = 2.0×10⁻⁴ m²/V·s; diffusion and air motion neglected, since the air at
  ~1 m/s is ≪ the ion drift at ~100–200 m/s).
- **Method of characteristics:** along an ion path 1/|ρ| = 1/|ρ₀| + τ/ε₀, with τ = ∫ds/|E|.
  Every mesh vertex is traced back along E to the wire its ions came from.
- **Emission model for nicked wires:** every wire emits the same current per metre of wire.
  The nicks, not the smooth-wire field, decide where ions are made, and they were spread
  along all wires. The total current is the **single calibration parameter**, set so the
  computed thrust ∫ρE_z dA (× each wire's real active length) equals the thrust inferred from
  the measured jet, T = ½ρ_air A v² = 2.730 mN (v = 1.21 m/s).
- Fixed-point iteration with 50 % under-relaxation, stopped after 30 iterations. The RMS change
  per iteration plateaus at a noise floor of about 3 V (< 0.03 % of 11 kV), set by re-tracing
  ~20 000 vertices each iteration. It never reaches the nominal 0.5 V tolerance, so the result is
  "converged to the noise floor", not to the tolerance. The integral results (thrust, current) agree
  to 0.02 % between runs, and a rerun is bit-identical.

## Results (calibrated)
| quantity | value |
|---|---|
| total corona current | **45.3 µA** (70 µA per metre of wire) |
| thrust | 2.73 mN (calibration target) |
| effective axial ion path, T·μ/I | **12.06 mm** (stage 1 Laplace estimate: 12.1 mm) |
| corona power V·I | 0.50 W |
| EHD efficiency ½Tv/(VI) | **0.33 %** |
| thrust per power | 5.5 mN/W |
| charge conservation (collected at tubes vs emitted) | 0.6 % |
| potential change from space charge, wire–tube gap | −270 V mean, −350 V max (≈ 2.5–3 % of 11 kV) |
| wire surface field (mean over the wire surface), Laplace → with space charge | end wires 14.95 → 13.07 MV/m; inner wires 9.6 → 8.2 MV/m (−13 to −15 %). Stage 1 and RESULTS.md quote the *maximum* surface field (15.9 / 10.3 MV/m), which is about 6 % higher |

## What this shows
1. **The simple thrust law is stable within the drift model.** With full space charge, T·μ/I = 12.06 mm,
   within 0.4 % of the field-line estimate (a model-internal result, not a test of the device). The ions barely change the field in the gap (~2.5 %), so
   T = I·d_axial/μ is accurate at this current.
2. **Space charge matters most at the wires:** with emission imposed, it lowers the computed
   surface field by 13–15 %. This is **not** extra evidence about onset. Below onset there are
   no ions, so onset is set by the Laplace field alone (stage 1). In real steady corona the
   surface field at an emitting point stays clamped near the onset value (Kaptzov's assumption).
   The low field here is therefore a property of the imposed-emission model. What it does imply
   is that, at the inferred current, emission must come from points where the local field is
   enhanced (the nicks): at the nominal current the model's inner-wire field (6.5 MV/m mean) is
   about 2.2× below the smooth-wire onset field.
3. The 45 µA here assumes an ideal jet. Stage 3 shows the collector and flow losses raise the
   inferred current to ~120 µA (see below).

## Checks
- **Box size:** rerun in a much larger air box (x ±70 mm, z −60…100 mm): thrust −0.2 %,
  effective path −0.2 %, wire surface fields within 2 % (`stage2_domain_check.py`).
  The largest potential change sits in the far field upstream and depends on the box
  (1.0 kV vs 3.1 kV), so it is not reported as physical; the gap values are box-independent.

## Limitations
- 2D, with wires treated as infinitely long; ions emitted towards the frame sides and ends are not modelled.
- Equal emission per metre of wire is an assumption about the nicks (they were not counted
  or mapped). Total current, thrust and efficiency depend on it only weakly, since every
  wire's ions travel about the same axial distance.
- Mobility and the actuator-disk thrust estimate carry the main uncertainties
  (±10 % and up to ×2 from the velocity accuracy); they scale the calibrated current directly.

---

# Stage 3: airflow (OpenFOAM v1912, `pimpleFoam`)

## Model
- 2D (one cell thick, 74 716 prism cells), laminar, transient. Air ν = 1.67×10⁻⁵ m²/s, ρ = 1.124 kg/m³.
- Geometry: 7 tubes (walls), both printed frames around the 65 mm opening (walls), and the open gap
  between the frames at the sides (it is open in the real build too). The 0.13 mm wires are not
  obstacles (3 % of the width). Open far-field boundaries: total pressure 0, so air can enter or leave.
- **Driving force:** the stage 2 ion force ρE at the calibrated 45.3 µA, mapped conservatively onto
  the CFD cells (total force preserved exactly: 0.04672 N per metre). The packaged OpenFOAM cannot
  compile coded sources, so the cell forces are grouped into 60 zones (k-means), each with its
  volume-weighted mean force. This keeps the total exactly and smooths the local distribution.
- Fine and coarse meshes run 0–3.0 s; fields written every 0.02 s; averages over **1.0–3.0 s**
  (101 snapshots). Before ~1 s the anemometer signal still trends upward: averaging from 0.8 s
  biased the fine-mesh mean low by about 2 %. Averages are computed in Python, because function
  objects crash in this packaged build.
- The anemometer is modelled as a 25 mm-wide vane (size not documented; assumption), 60 mm behind the tubes.
- The jet is unsteady (it flaps). The uncertainty of each time mean is a batch-means standard
  error (6 consecutive batches), which accounts for autocorrelation.

## Results (baseline force from stage 2, 45 µA)
| quantity | fine mesh (74 716 cells) | coarse mesh (37 892 cells) |
|---|---|---|
| ideal 2D actuator-disk jet speed (no losses) | 1.131 m/s | 1.131 m/s |
| speed at the anemometer, time mean 1.0–3.0 s | **0.654 ± 0.019 m/s** (time std 0.17; trend −0.01 m/s per s) | 0.690 ± 0.019 m/s (std 0.16; trend −0.02 m/s per s) |
| survival ratio r = anemometer / ideal | **0.578** | 0.610 |
| drag on tubes / frames, % of ion force | **62 % / 3 %** | 60 % / 2 % |
| momentum budget: (wall drag + mean jet flux) / ion force − 1, at 40 / 77 / 120 mm | +1 / +2 / +2 % | +3 / +4 / +5 % |

The budget does not close exactly because pressure and fluctuation (Reynolds-stress) terms are not
included in the jet flux.

## Re-run at a higher force (3.2×, coarse mesh, 0–2.0 s, averaged 1.0–2.0 s)
The current inferred from r is ~2.7× the stage 2 value, so the flow was re-run with 3.2× the ion force
to see how r depends on force:
- r = 0.639 vs 0.610 at the baseline force on the same mesh: a **+4.7 % force effect**. It was
  measured on the coarse mesh only and is applied to the fine-mesh r.
- Tube drag falls to **53 %** of the ion force at the higher Reynolds number.
- Speed at the anemometer: 1.29 m/s. Scaled to the nominal inferred force (2.7×, v ∝ √F) this is about
  1.19 m/s, against the measured 1.21 m/s. This is a **consistency check only, not an independent test**:
  the force was itself inferred from the measured speed. Also, the 2D force per unit width is 14.5 %
  lower than the 3D thrust per unit area, because the wires are chords of different lengths.

## What this gives (see `sensitivity.py`, `results/sensitivity.json`)
r = 0.605 nominal (fine mesh × force effect). The CFD range is 0.57–0.67: both meshes ±2 standard
errors, × force effect. This is then multiplied by 0.7–1.3, assumed for 2D → 3D. Using the measured jet:

| | nominal | 5–95 % range (Monte Carlo over all inputs) |
|---|---|---|
| corona current | **~120 µA** | 34–325 µA (anemometer accuracy alone: 35–262 µA) |
| electrical power (peak V × I) | ~1.35 W | 0.4–3.8 W |
| EHD efficiency, η_measured (theory.md §5) | **~0.12 %** | 0.05–0.25 % |
| jet thrust per electrical watt | ~2.0 mN/W | 1.0–3.7 mN/W (ion force per watt, ideal: 5.5 mN/W) |
| share of planar space-charge limit | ~21 % | 5–56 % |

Contribution of each input on its own (5–95 % factor on the current):

| input | factor |
|---|---|
| anemometer accuracy (±0.56 m/s) | ×0.34–2.0 |
| assumed 2D → 3D factor on r (±30 %) | ×0.62–1.9 |
| CFD r (meshes, time statistics) | ×0.82–1.10 |
| mobility (1.8–2.2×10⁻⁴) | ×0.91–1.09 |
| gap (9.9–12.4 mm) | ×0.91–1.11 |

The anemometer is the largest term, and the assumed 2D → 3D factor is nearly as large. Inviscid
actuator theory at the anemometer distance suggests the 3D jet keeps somewhat more of its speed
than the 2D one (r_3D / r_2D ≈ 1.1, which would lower the current by ~18 %). That lies inside the
assumed range but has not been computed.

`stage2_high_current_check.py` reruns the space-charge model at the nominal 122.5 µA:
- The effective axial path is 12.16 mm, vs 12.06 mm at 45 µA. T = I·d/μ is therefore stable within 1 %, *within the drift model*.
- The gap potential changes by ~−645 V (6 % of 11 kV).
- The inner-wire mean surface field falls to 6.5 MV/m. See stage 2, point 2, for what this does and does not mean.

## Robustness notes
- **Independent reproduction:** the packaged coarse case (`stage3_case_setup/ehd_coarse_case.zip`,
  0–2.0 s) was re-run on a different machine with a different OpenFOAM version (v2606 rc2 vs v1912).
  The anemometer velocity was 0.663 m/s vs 0.687 m/s from the same case and window here (0.8–2.0 s),
  a −3.5 % difference. That is within the time-statistics scatter of the unsteady jet.
- **Averaging window:** averaging from 0.4 s (as first reported) included the start-up transient.
  Averaging over 0.8–2.0 s was still slightly biased. The final window is 1.0–3.0 s, with no
  significant trend on either mesh.
- **Mesh:** the jet speed differs by ~5 % between meshes; the tube drag by 2 percentage points.
- **Reynolds number:** based on the 6 mm tube and the slot speed, Re ≈ 500 at the stage-2 force and
  ≈ 1000 at the higher force; the jet Re is ≈ 5000.

## Limitations
- 2D: the real jet comes from a 65 mm circle and spreads in 3D. A 2D jet flaps more and decays
  differently, and part of the 2D "loss" is geometric. The 2D → 3D factor on r is assumed (±30 %).
- Laminar 2D at these Reynolds numbers; the real 3D flow is partly turbulent.
- Two meshes only; no time-step study. The high-force case was run on the coarse mesh only, to 2.0 s.
- Force grouped into 60 zones (total exact, local distribution smoothed).
