# Final results: single source of numbers

Use only these numbers in the resume, SOP and write-up. Each row gives the value, its range,
and where it comes from. Updated 24 Nov 2026 after a second independent review, longer flow runs
(to 3.0 s, averaged over 1.0–3.0 s) and a rerun of every downstream step.

## Inputs (measured or reconstructed)
| quantity | value | range / uncertainty | source |
|---|---|---|---|
| emitter | 11 bare copper wires, 36 AWG (0.127 mm), knife-nicked | gauge from memory | user; photos |
| wire pitch (model) | 4.3 mm uniform | real 3.7–4.7 mm | CAD sketch |
| active emitter length | 644 mm | — | CAD chord lengths |
| collector | 7 tubes, 6 mm OD, 9 mm pitch, 33 % open area | — | CAD + measurement |
| gap, wire plane → tube surface | 11.11 mm | 9.9–12.4 mm | 13 M8 threads × 1.25 mm − 5.142 mm protrusion |
| air | 33 °C, 59 % RH, 99.9 kPa, ρ = 1.124 kg/m³ | ±0.3 % in ρ over 12–15 h | Open-Meteo archive; anemometer read 33.7 °C |
| supply peak voltage | ~11 kV, pulsed, waveform unknown | 9–14 kV | 7–8 mm sparkover, needle-gap tables |
| battery at test | 3.7 V open circuit (~half charged) | — | multimeter |
| jet velocity, ~60 mm downstream | 1.21 m/s (time std 0.06 over 39 s) | **instrument ±0.56 m/s (±46 %)** | video frames, Metravi ET-99A, hand-held |
| baseline, thruster off | 0.00 m/s (~48 s) | below vane start threshold | videos 1 and 3 |
| smooth vs nicked wires | faint hiss, no flow → flow | qualitative | user observation, incense smoke |
| ion mobility (O₂⁻, humid air) | 2.0×10⁻⁴ m²/V·s | 1.8–2.2×10⁻⁴ | literature range |

## Results
| result | value | range | how obtained |
|---|---|---|---|
| solver accuracy (wire above plate, exact solution) | −0.67 % | converges with mesh | `sim/benchmark_wire_plane.py` |
| maximum surface field at 11 kV, inner / end wires | 10.3–12.2 / 15.9 MV/m | ∝ V | stage 1 |
| single isolated wire, same V and gap (comparison) | 29.6 MV/m | — | benchmark |
| smooth-wire onset field (Peek, 63 µm wire: approximate) | 14.6 MV/m | 13.8 MV/m with Peek's parallel-wire constants | δ = 0.944, m = 1, coaxial constants |
| smooth-wire onset voltage, inner 7 wires | ≥ 14.5 kV | ≥ 11.8 kV (gap 9.9 mm, m = 0.9) | `sim/sensitivity.py` |
| wires above smooth onset (nominal inputs) | 0 of 11 at ≤ 10 kV; 2 at 11–13 kV; 4 at 14 kV | see above | `sim/sensitivity.py` |
| effective axial ion path (thrust) | 12.1 mm (gap 11.11 mm) | 12.06–12.16 mm (45–123 µA) | stages 1–2 |
| tube drag, % of ion force | **~55–60 %** | 53 % (3.2× force) – 62 % (baseline force) | stage 3, two meshes (frames add 1–3 %) |
| jet survival ratio r (anemometer / ideal 2D jet) | 0.60 | 0.57–0.67 (CFD: two meshes ±2 SE, force effect), × 0.7–1.3 assumed for 2D→3D | stage 3, averaged 1.0–3.0 s |
| independent reproduction of the coarse flow case | 0.663 m/s vs 0.687 m/s (−3.5 %) | — | user's laptop, OpenFOAM v2606 (rc2), 24 Nov 2026 |
| CFD jet speed at the nominal inferred force | ~1.19 m/s vs measured 1.21 | consistency check only (force was inferred from the measured speed) | coarse mesh 3.2× run, scaled to 2.7× |
| **corona current** | **~120 µA** | **34–325 µA** (5–95 %; anemometer alone 35–262 µA; 2D→3D factor alone ×0.62–1.9) | Monte Carlo, `sim/sensitivity.py` |
| electrical power (peak V × I) | ~1.35 W | 0.4–3.8 W | same |
| **EHD efficiency**, η_measured = ¼ρAv³ / (V·I) (theory.md §5) | **~0.12 %** | **0.05–0.25 %** | same |
| ideal efficiency if all ion force reached the jet, η_ideal = d·v / (2μV) | 0.33 % | — | theory.md §5 |
| thrust per electrical watt, delivered to the jet | ~2.0 mN/W | 1.0–3.7 mN/W | same (ion force per watt: 5.5 mN/W) |
| share of planar space-charge (Mott–Gurney) limit | ~21 % | 5–56 % | same; loose reference |
| smooth → nicked current increase | ~40–150× (order of magnitude) | smoke threshold 0.1–0.2 m/s; same flow losses assumed | (1.21 / v_smoke)² |

## Superseded estimates (do not quote)
| old value | replaced by | why |
|---|---|---|
| 45–49 µA, 0.31–0.33 % efficiency | ~120 µA, ~0.12 % | assumed an ideal jet; stage 3 showed collector losses |
| ~130 µA (40–370), ~0.11 % (0.04–0.22), ~1.5 W, r = 0.58 | ~120 µA (34–325), ~0.12 % (0.05–0.25), ~1.35 W, r = 0.60 | flow runs extended to 3 s; the 0.8–2.0 s window was still in the start-up trend; CFD range now from standard errors |
| "1.27 m/s at the inferred force, close to the measured 1.21" | consistency check only | the force was inferred from the measured speed, so agreement is built in |
| 100–150 µA, 0.08–0.18 % | 40–370 µA, 0.04–0.22 % | earlier ranges omitted anemometer accuracy, gap and 2D→3D |
| r = 0.59 (averaged from 0.4 s) | 0.58, then 0.60 | 0.4–0.8 s was start-up transient; force effect added; then window 1.0–3.0 s |
| "~60 % tube drag" | ~55–60 % | drag share falls at the higher (corrected) force |
| "≥ 40×" smooth → nicked | ~40–150× | follows from (1.21 / v_smoke)² |
| onset 20–29 kV (28 AWG) | 36 AWG array, see above | wire gauge corrected |
| "inner wires never reach onset across 9–14 kV" | only at nominal gap and m = 1 | pessimistic corner: 11.8 kV |

## Claims that are safe to make
1. In the 11-wire array, neighbouring wires shield each other (surface field cut to ~⅓ of a single
   wire). At nominal inputs, smooth wires would barely reach corona onset at the ~11 kV supply, which is
   consistent with the faint hiss and no flow observed before nicking. Surface sharpening was the likely
   reason flow appeared.
2. For this geometry, the drift model gives T = I·d/μ with d = 12.1 mm (net axial ion path), stable
   within 1 % when space charge is included.
3. In the 2D flow model, the round-tube collector took roughly half or more of the ion force as drag
   (~55–60 %): the largest loss in the model. Even without it, the physics caps the efficiency at
   η_ideal ≈ 0.33 %, because ions drift ~100× faster than the air moves.
4. The thruster ran in weak corona: of order 100 µA and of order 0.1 % efficiency, uncertain by about a
   factor of 3. The two largest terms are the anemometer accuracy and the assumed 2D → 3D flow
   factor.
5. Emission at that current required local field enhancement at the nicks. The model's wire field is
   about 2× below smooth-wire onset at the inferred current.

## Claims to avoid
- Any current or efficiency value without its range (no current was measured).
- "Validated" or "confirmed by simulation": the data support consistency with one operating point.
- "Only the nicks can explain it": wire ends and leads were not modelled.
- Exact polarity (the leads were unmarked; negative is the most likely).
- "1.21 ± 0.06 m/s": the 0.06 is time scatter; the instrument accuracy is ±0.56 m/s.
