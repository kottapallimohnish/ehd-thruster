# Physics model of a wire-to-collector EHD (ion-wind) thruster

Notation: V applied voltage, I corona current, d wire–collector gap, r wire radius,
L electrode span (length of wire), A = L·w flow area at the collector,
μ ion mobility in air, ρ air density, ε₀ vacuum permittivity.

The device has two regions:

1. **Ionisation region**: a thin sheath (≈ a few wire radii) around the wire where
   the field exceeds the breakdown field and electron avalanches produce ions.
2. **Drift region**: the rest of the gap. Ions of one sign (same polarity as the
   wire) drift to the collector, colliding with neutral molecules and transferring
   momentum to them. This is where thrust is produced.

---

## 1. Corona onset: Peek's law

The field at the surface of a thin wire parallel to a plane at distance d is

    E_surface ≈ V / ( r · ln(2d/r) )

Corona starts when E_surface reaches Peek's critical field for a cylinder:

    E_c = E₀ · m · δ · ( 1 + K / sqrt(δ · r) )
    E₀ ≈ 3.1 × 10⁶ V/m,  K ≈ 0.0308 m^½   (Peek's coaxial-cylinder constants)

Peek's parallel-wire form (E₀ = 3.0 × 10⁶ V/m, K = 0.0301 m^½) gives about 5 % lower fields.
This analysis uses the coaxial form. Neither form is exact for a wire array near a tube row,
so the onset field carries an uncertainty of at least ~5 % (see RESULTS.md).

δ = relative air density = (p/101.3 kPa)·(293 K/T); m ≈ 0.6–1 = wire surface
roughness factor (1 for a smooth, clean wire). So the onset voltage is

    V₀ = E_c · r · ln(2d/r)

Why the 1/sqrt(r) term: a thin wire's field falls off quickly, so an avalanche
needs a stronger surface field to grow over the short distance where the field
is high enough.

**Test conditions:** GITAM Visakhapatnam, 15 Aug 2023, 12:00–15:00 (weather archive): 99.9 kPa,
33 °C, 59 % RH, so δ = 0.944 and moist-air ρ = 1.124 kg/m³. For the 36 AWG wire (r = 0.0635 mm),
E_c ≈ 14.6 MV/m (smooth; 13.8 MV/m with the parallel-wire constants). In the 11-wire array the onset *voltage* is set by the field
simulation (sim/), not by the single-wire formula, because neighbouring wires shield each other.

## 2. Current–voltage relation (Townsend)

Above onset, the space charge of the drifting ions shields the wire. The current
is limited by how fast ions can leave, which gives the empirical/semi-analytic

    I = C · V · (V − V₀)

For the coaxial (wire-in-cylinder, radius R) case it can be derived in the near-onset approximation:
C = 8π ε₀ μ / (R² ln(R/r)) per unit length. For wire-to-plane, C is a geometric
constant that we **fit** from data. Linear check: plot I/V against V. It should be
a straight line with slope C and intercept −C·V₀.

## 3. Thrust: why T = I·d/μ

Force per unit volume on the air = force on the space charge, transferred by
collisions:

    f = ρ_e · E                          (ρ_e = ion charge density)

Ion current density with drift velocity μE (the air speed is ≪ μE, so it is ignored):

    j = ρ_e · μ · E   ⇒   ρ_e E = j / μ

So the force density is simply j/μ. Integrate over the gap volume. In a
quasi-1D picture (ions travel straight across distance d, current I crosses every
plane):

    T = ∫ (j/μ) dV = (1/μ) ∫₀ᵈ I dx = I · d / μ

**Mechanism:** mobility is in the denominator because a *less* mobile ion collides
more often per metre and hands its momentum to the air more effectively. Low
mobility = more thrust per ampere, but (see §5) also a higher voltage for the same
current. For curved ion paths, replace d with the **net axial displacement** of the ions,
not the arc length: sideways components of the force cancel, and an ion that first drifts
upstream pushes the air backwards for that part. (For this thruster: 12.1 mm, see sim/.)

## 4. Exit velocity (actuator-disk model)

Treat the electrode gap as an actuator disk of area A that adds a pressure jump
Δp = T/A to air starting from rest. Momentum and Bernoulli give

    velocity at the disk:      v_d = sqrt( T / (2ρA) )
    fully developed jet:       v_w = 2·v_d = sqrt( 2T / (ρA) )

Substitute T = I·d/μ:

    v_w = sqrt( 2·I·d / (ρ·μ·A) )        ⇒  v ∝ sqrt(I)

**Testable prediction:** measured velocity plotted against sqrt(I) should be a
straight line through the origin. The slope tells you the effective d/(μA).

**Measurement location matters by a factor of 2** between right at the collector
(v_d) and a few gap lengths downstream (v_w). It must go in the error budget.

Inverting to estimate thrust from a velocity measurement:

    from v at collector:  T ≈ 2 ρ A v_d²
    from v in the jet:    T ≈ ½ ρ A v_w²

## 5. Efficiency

Input power P = V·I. Useful power = kinetic energy flux of the jet
= ½ ṁ v_w² = ½ T v_w. With T = I d/μ:

    η_ideal = T·v_w / (2·V·I) = d·v_w / (2·μ·V) = v_w / (2·μ·Ē),   Ē = V/d

**Mechanism:** η_ideal is *half* the ratio of air speed to ion drift speed (μĒ). Ions drift at
~100–300 m/s and the air moves at a few m/s. The rest of the input power heats the
air through ion–neutral collisions (Joule heating of the drift region). For this thruster
(d = 12.1 mm, v_w = 1.21 m/s, V = 11 kV, μ = 2 × 10⁻⁴ m²/V·s) η_ideal ≈ 0.33 %. This is the
efficiency if **all** of the ion force ended up in the jet. It is not a construction flaw; it is
set by the physics.

The efficiency quoted in RESULTS.md (≈ 0.1 %) is a different quantity:

    η_measured = (¼ ρ A v_w³) / (V · I_inferred)

This is the measured jet power divided by the electrical power at the inferred current. It is
lower than η_ideal because part of the ion force is lost to drag on the collector tubes, so
more current (and power) is needed to produce the measured jet. Since I_inferred = I_ideal / r²,
η_measured = η_ideal × r², where r is the jet survival ratio from the flow simulation.

Thrust-to-power, independent of current:

    T / P = d / (μ·V) = 1 / (μ·Ē)

For this thruster 1/(μĒ) ≈ 5.5 mN/W. That is the ion force per watt. Net of collector drag,
the thrust delivered to the jet per watt is lower (see RESULTS.md). Thrust-to-power is the
standard figure of merit for EHD propulsion.

Lower field, higher T/P, but less thrust per unit area. That is the fundamental
trade-off of EHD propulsion.

## 6. Space-charge-limited current: collisional vs collisionless

The maximum current a gap can carry is set by the ions' own space charge, but the
scaling depends on how the ions move:

| | Ions move by | Limit | Scaling |
|---|---|---|---|
| EHD in air (collisional) | drift, v = μE | Mott–Gurney | J = (9/8) ε₀ μ V² / d³ |
| Gridded ion thruster (vacuum, collisionless) | free acceleration, v = sqrt(2qV/m) | Child–Langmuir | J = (4/9) ε₀ sqrt(2q/m) V^{3/2} / d² |

In vacuum, all of the electrical energy qV becomes ion kinetic energy and leaves as
exhaust momentum, so the efficiency can be high. In air, each ion reaches only its
drift speed and dumps its energy into many small collisions, which is why EHD
efficiency is set by v/(μE) as in §5. The same space-charge physics governs plasma
sheaths, which connects this device to plasma and fusion physics.

## Parameter ranges

| Quantity | Value | Note |
|---|---|---|
| μ, positive ions (N₂⁺/O₂⁺ clusters) | 1.4–2.2 × 10⁻⁴ m²/V·s | humidity lowers it |
| μ, negative ions (O₂⁻) | 1.8–2.7 × 10⁻⁴ m²/V·s | |
| ρ, air | 1.20 kg/m³ at 20 °C sea level; 1.124 kg/m³ at the test (33 °C, 59 % RH, 99.9 kPa) | |
| E₀ (Peek) | 3.1 × 10⁶ V/m | |
