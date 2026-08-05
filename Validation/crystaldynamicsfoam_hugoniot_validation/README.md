# crystaldynamicsFoam planar-impact Hugoniot validation (Mie–Grüneisen EOS mode)

## What this case tests, and why the elastic wave-speed case cannot

The existing two-wave benchmark is a **linear acoustics** test: infinitesimal
amplitude, so it probes exactly one number in the volumetric response — the
slope ∂p/∂ρ at the reference state. The legacy closure (p = Kμ, implied by
the trace of the hypo-elastic stress rate) and the Mie–Grüneisen EOS with
c₀ = √(K/ρ₀) share that slope **by construction**; they are built to agree
there. The acoustic limit is therefore degenerate for this purpose: every
equation of state with the same bulk modulus produces identical wave speeds,
and the case will pass with flying colours whether the shock physics is right
or absent.

The physics the EOS adds lives strictly at **finite amplitude**:

1. the curvature of the cold curve — the s in Uₛ = c₀ + s·uₚ (copper:
   1.489, versus the legacy closure's effective 0.5);
2. the thermal pressure Γ₀ρ₀e (release paths, reshock, thermal stress);
3. the cold/thermal energy partition (shock temperatures).

A planar impact measures (1) directly and cleanly: drive a sustained shock,
measure its speed against the particle velocity, and compare with the
material's measured Hugoniot. This is the canonical shock-physics validation
because Uₛ–uₚ data **is** how these EOS parameters were determined
experimentally in the first place.

The piston boundary condition used here (patch velocity = uₚ) is exactly
equivalent to a symmetric flyer impact at 2uₚ, without the flyer's
back-surface release wave contaminating the late-time trajectory.

## Expected results (copper, [100] loading)

| uₚ [m/s] | Uₛ Mie–Grüneisen [m/s] | p [GPa] | Uₛ legacy p=Kμ [m/s] | p [GPa] |
|---------:|-----------------------:|--------:|---------------------:|--------:|
|      200 |                 4237.8 |    7.57 |                 4019 |    7.18 |
|      400 |                 4535.6 |   16.20 |                 4123 |   14.73 |
|      600 |                 4833.4 |   25.90 |                 4229 |   22.66 |
|      800 |                 5131.2 |   36.66 |                 4338 |   30.99 |
|     1200 |                 5726.8 |   61.37 |                 4563 |   48.90 |
|     1600 |                 6322.4 |   90.32 |                 4799 |   68.58 |

Fits: **MieGruneisen** must recover c₀ ≈ 3940 m/s, s ≈ 1.489 (PASS tolerance
±80 m/s, ±0.12 — the strength offset and front-tracking discretisation
consume part of this). **linearElastic** lands near c₀ ≈ 3899 m/s,
effective s ≈ 0.56 over the six-point range (its exact locus is the curve
Uₛ = √(K/ρ₀)·√(1+μ), slightly convex): at 800 m/s the shock is 15% slow and the pressure 15%
low, growing without bound at stronger drive. That before/after pair is the
figure for the manuscript's verification section.

## Usage

```
./Allrun                    # both modes, 6 velocities each (12 serial runs
                            #   launched in parallel), then validation
./Allrun MieGruneisen       # one mode only
APP=/path/to/crystaldynamicsFoam ./Allrun both
./Allclean
```

Requires: the crystaldynamicsFoam build **with the EOS extension files**, blockMesh,
python3 (matplotlib optional — figures skipped gracefully without it).
Each run is 800 cells × ~2400 steps; minutes per run on one core.

Outputs: per-run logs in `runs/<mode>/up<value>/log.crystaldynamicsFoam`, fitted
Hugoniots printed to stdout, `hugoniot_UsUp.png` and per-mode profile
figures. The validator exits nonzero if the MieGruneisen fit misses the
input parameters.

## Design notes

**Strength is deliberately minimised** (CRSS 2 MPa, k = 5, no hardening →
flow stress ~50 MPa at shock strain rates) so the front kinematics isolate
the volumetric closure under test. A constant strength offset Δσ perturbs
Uₛ by Δσ/(ρ₀uₚ) — hyperbolic in uₚ — which with these parameters biases the
fitted s by under 0.05. Shock-regime strength physics (Steinberg–Guinan) is
a separate, later validation.

**uₚ = 200 m/s sits below the [100] overdrive threshold** (elastic precursor
at √(C₁₁/ρ₀) ≈ 4343 m/s; plastic shock overtakes it above uₚ ≈ 270 m/s).
The precursor amplitude at this CRSS is a few m/s, far under the half-max
tracking threshold, so the fit is unaffected — but you will see it in the
profiles, and it is physical.

**Legacy-mode profiles ring behind the front.** That is the known colocated
grad-then-div dispersion; the legacy path has no artificial viscosity by
design (it is byte-identical to the current solver). The EOS-mode profiles
are cleaned by the compression-only q, which also demonstrates that the
regularisation does not bias the measured Uₛ (the fit must still hit the
input line).

**Dictionary K vs ρ₀c₀².** C₁₁, C₁₂ give K = 137.1 GPa; the EOS carries
ρ₀c₀² = 138.6 GPa. The ~1% inconsistency is inherited from independent
experimental sources and is immaterial: the deviator uses C, the pressure
uses the EOS.

**Front tracking** uses the half-amplitude crossing of the uₓ profile with
the formation transient (first 35% of the record) excluded from the x–t fit.
Cell centres are reconstructed analytically from the known uniform mesh; if
you change the mesh, update `L_BAR`/`N_CELLS` at the top of
`validate_hugoniot.py`.

## What to run next (beyond this case)

The Hugoniot fit cannot see Γ₀ — on the Hugoniot, e = e_H identically. The
follow-up validations, in order: (i) a release test (drop the piston
velocity mid-run; the unloading wave rides the isentrope, which the legacy
closure gets wrong even at modest pressure); (ii) elastic-precursor +
plastic-wave structure at realistic copper strength against gas-gun VISAR
data, including the overdrive collapse of the two-wave structure.

## Experimental overlay (validation, not just verification)

`validate_hugoniot.py` reads `cu_hugoniot_experimental.dat` (or a file
given as its first argument): whitespace-separated `up Us source` rows,
m/s or km/s auto-detected, `#` comments. Populated points are drawn as
per-source scatter on the Uₛ–uₚ figure and the validator prints their rms
residual and bias against BOTH closures — data landing on the MG line and
off the legacy curve is the validation claim, made quantitative.

The file ships **empty of data by design**: transcribe real shots from
the sources listed in its header — Marsh's LASL compilation (1980),
Mitchell & Nellis, J. Appl. Phys. 52, 3363 (1981) (OFHC Cu, 30–430 GPa;
the uₚ = 1200 and 1600 m/s cases sit inside their range), or the IHED
RusBank database, which exports plain text in exactly this two-column
form. Do not re-fit; do not trust remembered numbers. Until then the
figure carries a ±0.9% band around the input line — the fractional
standard deviation Mitchell & Nellis report for gas-gun data about linear
Uₛ–uₚ fits — as the published-scatter reference.
