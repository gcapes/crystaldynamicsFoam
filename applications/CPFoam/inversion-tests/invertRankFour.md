# Rank-Four Stiffness Tensor Inversion

The mechanical driving force for phase-field fracture is the elastic strain energy density
at the current stress state:

    dfdp_mech = sigma : S : sigma

where sigma is the Cauchy stress and S = C^{-1} is the compliance tensor. The compliance
is needed because the constitutive update outputs stress directly, so contracting back through
the stiffness to get strain is avoided by working with S instead. The contraction is done
by Material::doubleContractionTwice.

For a physically admissible elastic material the stiffness tensor is positive definite, so
the associated Mandel matrix has strictly positive eigenvalues and is invertible. Full minor
and major symmetry must also hold, otherwise the reduction from 81 components to a 6x6
matrix is not well-defined:

    minor symmetry:  C_ijkl = C_jikl = C_ijlk  (symmetry of stress and strain tensors)
    major symmetry:  C_ijkl = C_klij            (existence of a strain energy density)

These reduce the 81 components to 21 independent values.


## Symmetry enforcement

At initialisation (createFields.H) the raw tensor from GrainProperties is passed through
Material::castSymmetry, which computes one average per symmetry equivalence class and writes
the exact same double to all related positions. A WarningInFunction is issued for any group
of components that are not equal before averaging. Material::assertSymmetry then does a
strict fatal check (tolerance 1e-12) on the result.

Symmetry is not re-checked on every call to invertRankFour during the solve. After
castSymmetry the tensor is guaranteed symmetric by construction, so the per-cell check
would be wasted overhead. A safety check remains inside fourthToMandel as a backstop.


## Inversion pipeline

Material::invertRankFour(C, verbosity) runs four steps.

Step 1 - castSymmetry(C)
Projects onto the fully symmetric subspace. For each of the 21 Mandel pairs (I <= J):

    avg = mean{ C(i,j,k,l), C(j,i,k,l), C(i,j,l,k), C(j,i,l,k),
                C(k,l,i,j), C(l,k,i,j), C(k,l,j,i), C(l,k,j,i) }

All 8 positions are set to the same double, giving bitwise equality.

Step 2 - fourthToMandel(C)
Converts to a 6x6 Mandel matrix via the orthonormal Mandel basis tensors M_I:

    CM(I,J) = Sum_{abcd} M(I,a,b) * M(J,c,d) * C(a,b,c,d)

With verbosity = true the maximum symmetry residuals are printed before the safety checks.

Step 3 - Fastor::inverse(M)
Dense 6x6 matrix inversion.

Step 4 - mandelToFourth(M_inv)
Converts back to a 3x3x3x3 tensor:

    C_inv(i,j,k,l) = Sum_{AB} M(A,i,j) * M(B,k,l) * M_inv(A,B)


## Floating-point noise after crystal rotation

The stiffness is rotated each time step by:

    C_rot(i,j,k,l) = Sum_{mnpq} R(i,m) R(j,n) R(k,p) R(l,q) C(m,n,p,q)

This is an 81-term sum. Although C_rot(i,j,k,l) and C_rot(j,i,k,l) are mathematically
equal (rotation preserves symmetry), floating-point arithmetic computes the sum in a
different order for each index combination. For components of order 1e+n the residuals are
roughly eps * C_max ~ 2e(n-16). For the Fracture test case (C_max ~ 1e10) this is ~1e-6.

The original fourthToMandel check used a tolerance of 1e-12, six orders of magnitude below
this noise level, hence the spurious fatal error. The castSymmetry fix eliminates this: one
average per equivalence class, same bits written to all related positions, residuals entering
fourthToMandel are identically zero.

See inversion-test-1/TomsError.log.CPFoam for the diagnostic output confirming the noise level.


## Why Mandel, not Voigt

Both representations reduce the symmetric rank-four tensor to a 6x6 matrix. Voigt uses
scaling factors of 1 on all components; Mandel scales shear entries by sqrt(2), making the
basis orthonormal. This matters here for three reasons.

The double contraction sigma:S:sigma = sigma_M^T * S_M * sigma_M holds directly in Mandel
without correction factors on the shear blocks. In Voigt you need manual factors of 2 and 4.

The inverse of the Mandel stiffness matrix is directly the Mandel compliance matrix. The
Voigt inverse requires additional rescaling to recover the true compliance.

The eigenvalues of the Mandel matrix equal those of the fourth-order stiffness operator on
symmetric tensors, so spectral checks are meaningful.

doubleContractionTwice uses the Mandel inverse directly, so the representation must be
consistent end-to-end.
