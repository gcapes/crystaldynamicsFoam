/*---------------------------------------------------------------------------*\
  Polar Decomposition Utility for Fastor Tensors
  
  Computes F = R·U where:
  - R is the rotation tensor (proper orthogonal, det(R) = +1)
  - U is the right stretch tensor (symmetric positive definite)
  
  Method: Newton-Raphson iteration for U = sqrt(C) where C = F^T·F,
  then R = F·U^{-1}. The R matrix is orthogonalized via modified 
  Gram-Schmidt to handle numerical drift.
  
\*---------------------------------------------------------------------------*/

#ifndef polarDecomposition_H
#define polarDecomposition_H

#include "Fastor.h"
#include "fvCFD.H"

void polarDecomposition(
    const Fastor::Tensor<double, 3, 3>& F,
    Fastor::Tensor<double, 3, 3>& R,
    Fastor::Tensor<double, 3, 3>& U
)
{
    // Identity
    Fastor::Tensor<double, 3, 3> I_d(0.0);
    for (int i = 0; i < 3; i++) I_d(i, i) = 1.0;

    // Right Cauchy-Green tensor: C = F^T · F
    Fastor::Tensor<double, 3, 3> C = Fastor::matmul(Fastor::transpose(F), F);

    // Initial guess for U = sqrt(C)
    // For small deformation: U ≈ I + 0.5*(C - I)
    // For large deformation: start closer to C
    double traceC = C(0,0) + C(1,1) + C(2,2);
    if (traceC < 6.0)
    {
        U = I_d + 0.5 * (C - I_d);
    }
    else
    {
        U = I_d + 0.75 * (C - I_d);
    }

    // Newton-Raphson: U_{n+1} = 0.5 * (U_n + C · U_n^{-1})
    const int maxIter = 50;
    const double tol = 1.0e-11;

    for (int iter = 0; iter < maxIter; iter++)
    {
        Fastor::Tensor<double, 3, 3> U_inv = Fastor::inverse(U);
        Fastor::Tensor<double, 3, 3> U_new = 0.5 * (U + Fastor::matmul(C, U_inv));

        // Convergence check (Frobenius norm of change)
        double err = 0.0;
        for (int i = 0; i < 3; i++)
            for (int j = 0; j < 3; j++)
                err += Foam::sqr(U_new(i,j) - U(i,j));
        err = Foam::sqrt(err);

        U = U_new;

        if (err < tol) break;

        if (iter == maxIter - 1)
        {
            Info << "Warning: Polar decomposition NR did not converge! "
                 << "err = " << err << ", traceC = " << traceC << endl;
        }
    }

    // R = F · U^{-1}
    Fastor::Tensor<double, 3, 3> U_inv = Fastor::inverse(U);
    R = Fastor::matmul(F, U_inv);

    // Ensure det(R) = +1 (proper rotation)
    double detR = Fastor::determinant(R);
    if (detR < 0.0)
    {
        for (int i = 0; i < 3; i++)
            for (int j = 0; j < 3; j++)
                R(i,j) = -R(i,j);
    }

    // Orthogonalize R via modified Gram-Schmidt (columns)
    Fastor::Tensor<double, 3> r0(0.0), r1(0.0), r2(0.0);

    // Extract columns
    for (int i = 0; i < 3; i++)
    {
        r0(i) = R(i, 0);
        r1(i) = R(i, 1);
        r2(i) = R(i, 2);
    }

    // Normalize first column
    double norm0 = Fastor::norm(r0);
    if (norm0 > SMALL) r0 = r0 / norm0;

    // Orthogonalize second column
    double dot01 = Fastor::inner(r0, r1);
    r1 = r1 - dot01 * r0;
    double norm1 = Fastor::norm(r1);
    if (norm1 > SMALL) r1 = r1 / norm1;

    // Third column from cross product (ensures right-handed)
    r2(0) = r0(1)*r1(2) - r0(2)*r1(1);
    r2(1) = r0(2)*r1(0) - r0(0)*r1(2);
    r2(2) = r0(0)*r1(1) - r0(1)*r1(0);

    // Rebuild R from orthogonalized columns
    for (int i = 0; i < 3; i++)
    {
        R(i, 0) = r0(i);
        R(i, 1) = r1(i);
        R(i, 2) = r2(i);
    }
}

#endif

// ************************************************************************* //