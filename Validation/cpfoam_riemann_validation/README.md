# CPFoam solid Riemann problem (Sod-type initial-value benchmark)

## What this uniquely certifies

The Hugoniot suite validates (c₀, s); the release benchmark validates Γ₀
and the state-dependent sound speed. Both leave one thing untested: the
**energy partition** — the split of e into stored cold-compression energy
and thermal energy that determines temperature. This benchmark's central
observable, the **entropy contact**, is made of exactly that.

Initial data (a pure initial-value problem — no piston, both halves at
rest, diaphragm at x₀ = 5 mm of a 10 mm bar):

- **left**: the principal-Hugoniot state of uₚ = 800 m/s, at rest —
  ρ_L = 10579.4 kg/m³, e_L = uₚ²/2 = 3.2×10⁵ J/kg (the Hugoniot energy
  identity), p_L = 36.66 GPa, *hot*;
- **right**: ambient copper (ρ₀, e = 0, p = 0).

At t = 0⁺ this breaks into a right-running **shock** into the ambient
material, a left-running **release fan** into the hot compressed material,
and between them a **contact**: pressure and velocity are continuous, but
density and temperature jump, because the two sides reach the same
(p*, u*) along different thermodynamic paths — single shock (right) versus
shock-then-release (left). The left side arrives hotter and lighter. A
solver only reproduces that jump if its energy bookkeeping — the e
equation, the p(ρ,e) coupling, and the e→T recovery — is right. As a
bonus, the case exercises the initial-condition pathway (the `0/e` read
and setFields regions) that the piston-driven benchmarks never touch.

## Exact solution and expected numbers

`validate_riemann.py` solves the star state by impedance matching: the
right-shock locus p_s(u) = ρ₀(c₀ + s·u)u against the left-release
isentrope p_r(u) integrated (RK4) through the left state with the
left-facing Riemann invariant du/dμ = −c/(1+μ). Temperature uses the
solver's own leading-order convention T = T_ref + (e − e_H(μ))/Cᵥ, so this
is verification against the model's exact solution — the honest framing.

|                              | MieGruneisen | linearElastic (legacy) |
|------------------------------|-------------:|-----------------------:|
| p_L implied by ρ_L [GPa]     |        36.66 |     25.32 (p = Kμ_L)   |
| u* [m/s]                     |        401.1 |                  332.1 |
| p* [GPa]                     |        16.25 |                  12.12 |
| Uₛ (right shock) [m/s]       |       4537.3 |                 4087.3 |
| contact Δρ [kg/m³]           |  **−34** (9762 \| 9796) |        **0** |
| contact ΔT [K]               | **+110** (410 \| 300)   |        **0** |

The legacy closure is athermal: same p* → same μ → same ρ on both sides —
it *cannot* produce the contact jump. The star state itself discriminates
by +34% in p* and +21% in u*. Note T_right = 300 K exactly is a property
of the leading-order energy split (all single-shock Hugoniot energy is
booked as stored); the *jump* ΔT is the meaningful observable, and both
solver and exact solution use the same convention.

PASS criteria (MieGruneisen): u* and p* within 3%/4%, Uₛ within 3%,
contact ΔT within 20% of +110 K, and contact continuity |Δu| < 1.5% of u*,
|Δp| < 2% of p* (a solver that generates spurious pressure blips at
contacts fails here — a classic multi-material pathology).

## Mode-specific initialisation (why two setFieldsDicts)

The Mie–Grüneisen state lives in (ρ, e): setFields writes the ρ and e
discontinuity and the EOS produces p_L = 36.66 GPa on the first corrector.
The legacy closure has no EOS — its pressure lives inside the *integrated*
stress — so initialising only ρ would leave it stress-free and nothing
would happen. Its setFieldsDict therefore writes σ.0 = −(Kμ_L)I = −25.32
GPa·I in the left region: each mode gets the initial condition that is
exact *within its own closure*, and each is validated against its own
exact Riemann solution. (The one-corrector momentum lag at start-up —
sigmaGlobal begins at zero — is a single-step transient absorbed by the
tolerances.)

## Usage

```
./Allrun                    # both modes (parallel) + validation + figures
./Allrun MieGruneisen
APP=/path/to/CPFoam ./Allrun both
./Allclean
```

Outputs: exact-solution summary and PASS/FAIL to stdout;
`riemann_profiles_<mode>.png` (four panels — u, p, ρ, T against x with the
exact solution overlaid; the T panel's step at the contact is the money
shot); `riemann_pu_diagram.png` (the classic Riemann diagram: shock locus,
release isentrope, their intersection, and the measured plateaus).

## Design notes

**Contact smearing.** TVD advection diffuses linearly degenerate fields,
so the ρ and T steps at the contact spread over ~10–20 cells by the final
time and keep spreading. This is a property of the advection scheme, not
the EOS — the measurement windows sit clear of the smear, and the plateau
values are the metrics. (Interface-sharpening or WENO reconstruction of
the advected fields is the known remedy if sharp contacts ever matter;
this is precisely the advected-history-field territory where WENO schemes
apply, unlike the stress-divergence terms.)

**Timing.** endTime 0.8 μs keeps everything inside the bar: right shock
reaches 8.63 mm, release head 0.74 mm, contact 5.32 mm. Neither boundary
is ever touched.

**Strength** is minimised as in the other benchmarks (~50 MPa flow stress)
so the comparison is against the hydrodynamic exact solution; the residual
deviatoric contribution is inside the tolerances.

**GrainProperties** carries the same fully annotated EOS parameter
documentation as the release benchmark.

## The verification trilogy

Hugoniot suite → (c₀, s). Release benchmark → Γ₀ and c(state). Riemann
problem → energy partition, contact physics, and the initial-value
pathway. Together they close verification of the full p(ρ, e) surface
and its coupling to the conservative transport — the complete EOS tier of
the CPC paper's verification section. The tier above (elastic precursor
structure and decay against single-crystal gas-gun data) is gated on
calibrated rate-dependent strength.
