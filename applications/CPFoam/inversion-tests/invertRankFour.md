# Rank-Four Stiffness Tensor Inversion

## Why We Need the Inverse

In the phase-field fracture model the mechanical contribution to the order-parameter
driving force is the elastic strain energy density evaluated at the current stress state:

    dfdp_mech = sigma : S : sigma

where `sigma` is the Cauchy stress tensor and `S = inv(C)` is the elastic compliance
tensor (the inverse of the stiffness `C`).  The compliance is needed because the
stress — not the strain — is the natural output of the constitutive update, so
computing the energy density requires contracting `sigma` through the inverse of the
stiffness.  This contraction is performed by `Material::doubleContractionTwice`.

---

## Why the Inverse Exists

The inverse exists if and only if `C` is non-singular as a linear map on the space of
symmetric second-order tensors.  For a physically admissible elastic material this is
guaranteed by two conditions:

1. **Positive definiteness** — the strain energy density `W =0.5 epsilon:C:epsilon > 0` for all
   non-zero symmetric strains `epsilon`.  This means all eigenvalues of the associated
   6×6 Mandel matrix are strictly positive, so the matrix is invertible.

2. **Full symmetry** — both minor and major symmetries must hold (see below).
   Without them the mapping from the 81-component tensor to a 6×6 matrix is
   ill-defined and inversion is meaningless.

---

## Assumptions

The following symmetries are*assumed to hold for any stiffness tensor passed into
the inversion pipeline:

| Symmetry | Condition | Physical origin |
|---|---|---|
| Minor (ij) | C_ijkl = C_jikl | Symmetry of the stress tensor |
| Minor (kl) | C_ijkl = C_ijlk | Symmetry of the strain tensor |
| Major | C_ijkl = C_klij | Existence of a strain energy density |

Together these reduce the 81 independent components to 21, which map one-to-one
onto the entries of a symmetric 6×6 Mandel matrix.

Positive definiteness is checked at the `mandelToFourth` stage by requiring all
diagonal entries of the inverted Mandel matrix to be positive.

---

## Symmetry Enforcement and Verification

### At initialisation (`createFields.H`)

1. The raw stiffness tensor is read component-by-component from `GrainProperties`.
   A `WarningInFunction` is issued for any group of four minor-symmetry-related
   components that are not all equal — this catches accidental user input errors.

2. `Material::castSymmetry` is called to project the raw tensor onto the symmetric
   subspace (see Method below).  The result is stored in `Cprime2[i]`.

3. `Material::assertSymmetry` performs strict `FatalErrorInFunction` checks
   (tolerance 1e-12) on `Cprime2[i]`.  If the tensor is genuinely non-symmetric
   after symmetrisation the simulation aborts immediately with a clear message.

### During the solve (`invertRankFour`)

Symmetry checks are **not repeated** every call.  After `castSymmetry` the tensor is
guaranteed symmetric by construction (see Floating-point note below), so re-checking
every cell every time step is unnecessary overhead.  A runtime safety check remains
inside `fourthToMandel` and will abort if it ever sees a violation, but under normal
operation it will never fire.

---

## Method

Inversion is performed by `Material::invertRankFour(C, verbosity)` in four steps:

### Step 1 — `castSymmetry(C)`

Projects the input onto the fully symmetric subspace by averaging all 8
symmetry-related components for each of the 21 unique Mandel pairs `(I ≤ J)`:

    avg = mean{ C(i,j,k,l), C(j,i,k,l), C(i,j,l,k), C(j,i,l,k),
                C(k,l,i,j), C(l,k,i,j), C(k,l,j,i), C(l,k,j,i) }

The same `double` value is then written to all 8 positions simultaneously.
This guarantees bitwise equality between all symmetry-equivalent entries — not
merely equality to within a tolerance — which is required for the checks in
`fourthToMandel` to pass.

### Step 2 — `fourthToMandel(C, verbosity)`

Converts the symmetric 3×3×3×3 tensor to a 6×6 Mandel matrix using the
orthonormal Mandel basis tensors `M_I`:

    CM(I,J) = Sum_{abcd} M_(I,a,b) · M_(J,c,d) · C(a,b,c,d)

If `verbosity == true`, the maximum symmetry residuals are printed before the
safety checks run:

    fourthToMandel symmetry diagnostics:
      max |C(i,j,k,l) - C(j,i,k,l)| = ...
      max |C(i,j,k,l) - C(i,j,l,k)| = ...
      max |C(i,j,k,l) - C(k,l,i,j)| = ...

### Step 3 — `Fastor::inverse(M)`

Inverts the 6×6 Mandel matrix using the Fastor library's dense matrix inverse.

### Step 4 — `mandelToFourth(M_inv, verbosity)`

Converts the inverted Mandel matrix back to a 3×3×3×3 tensor:

    C_inv(i,j,k,l) = Sum_{AB} M_(A,i,j) · M_(B,k,l) · M_inv(A,B)

Symmetry of the Mandel matrix and positive definiteness of its diagonal are
checked before conversion.

---

## Verbosity

The full pipeline — input tensor, symmetrised tensor, Mandel matrix, inverted Mandel
matrix, compliance tensor, and symmetry diagnostics — is printed by passing
`verbosity = true`:

```cpp
Fastor::Tensor<double,3,3,3,3> C_rot_inverse =
    Material::invertRankFour(C_rotated[active_phases[count_i]][cell], true);
```

---

## Floating-Point Noise After Crystal Rotation

The stiffness tensor is rotated each time step by a 4-index contraction with the
rotation matrix `R`:

    C_rot(i,j,k,l) = Sum_{mnpq} R(i,m) R(j,n) R(k,p) R(l,q) C(m,n,p,q)

This summation involves 81 terms.  Although `C_rot(i,j,k,l)` and `C_rot(j,i,k,l)`
are mathematically equal (rotation preserves symmetry), floating-point arithmetic
computes the 81-term sum in a different order for each index combination, producing
tiny differences.

**For stiffness components of order 1e+n, the resulting symmetry residuals are
approximately `epsilon × C_max ~ 2×10^(n-16)`.** For the Fracture test case
(`C_max ~ 1e10`, `n = 10`) this gives residuals of order 1e-6, confirmed by
running the diagnostic output (see `Eulerian_Mechanics/applications/CPFoam/inversion-tests/inversion-test-1/TomsError.log.CPFoam`).

The original (pre 11/03/2026) `fourthToMandel` symmetry check used a tolerance of 1e-12, which
is six orders of magnitude smaller than this noise - hence the spurious fatal error.

The fix in `castSymmetry` (Mandel pair enumeration) eliminates this entirely by
computing one average per equivalence class and writing the exact same bits to all
related positions, so the residuals entering `fourthToMandel` are identically zero.

---

## Why Mandel, Not Voigt

Both Mandel and Voigt representations reduce the symmetric rank-four tensor to a 6×6
matrix, but they differ in how shear components are scaled:

| | Normal components | Shear components |
|---|---|---|
| **Voigt** | 1 | 1 |
| **Mandel** | 1 | √2 |

The Mandel scaling makes the basis **orthonormal** under the standard Euclidean inner
product.  This has several important consequences for the inversion:

1. **Inner product preservation** — `A:B = A_M · B_M` exactly, with no correction
   factors.  With Voigt you must manually account for factors of 2 on shear terms.

2. **Correct double contraction** — `σ:S:σ = σ_M^T · S_M · σ_M` directly.  This is
   how `doubleContractionTwice` uses the inverted tensor.

3. **Inversion consistency** — the inverse of the Mandel stiffness matrix is directly
   the Mandel compliance matrix.  With Voigt, the raw matrix inverse requires
   additional scaling by factors of 2 and 4 on certain blocks to recover the true
   compliance.

4. **Spectral correctness** — the eigenvalues of the Mandel matrix equal those of the
   fourth-order stiffness operator acting on symmetric second-order tensors.
