/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "Material.H"
#include "Time.H"
#include "subCycle.H"
#include "MULES.H"
#include "surfaceInterpolate.H"
#include "fvcGrad.H"
#include "fvcSnGrad.H"
#include "fvcDiv.H"
#include "fvcFlux.H"


#include "fvCFD.H"

#include "fvmDdt.H"
#include "fvmLaplacian.H"
#include "fvmSup.H"


// Constructors are defaulted in Material.H (rule of zero): all members
// are value types with in-class initialisers, so the compiler-generated
// memberwise copy is correct.  The previous definitions here included an
// empty-bodied copy constructor that copied nothing.



// int Foam::Material::returnComponents()
//{
//     return N_slips_;
// }

void Foam::Material::Print(label length)
{
}

double Foam::Material::deltaG( scalar *mu_arr,scalar *xi_arr,const scalar& T,const label& length)
{
    // PhaseComposition[count_phase][count_i][cell];
    double G_Quad=0.0;
    return G_Quad;
}

tensor Foam::Material::RK4(scalar step_size)
{
    return tensor::zero;
}

const Fastor::Tensor<double,6,3,3> Material::M_ = []() 
{
    Fastor::Tensor<double,6,3,3> M; 
    M.zeros();

    const double s2 = std::sqrt(2.0);

    // I = 0 -> (11)
    M(0,0,0) = 1.0;
    // I = 1 -> (22)
    M(1,1,1) = 1.0;
    // I = 2 -> (33)
    M(2,2,2) = 1.0;

    // I = 3 -> (23) = (32)
    M(3,1,2) = s2/2.0;
    M(3,2,1) = s2/2.0;

    // I = 4 -> (13) = (31)
    M(4,0,2) = s2/2.0;
    M(4,2,0) = s2/2.0;

    // I = 5 -> (12) = (21)
    M(5,0,1) = s2/2.0;
    M(5,1,0) = s2/2.0;

    return M;
}();

/*
    Projects an arbitrary fourth‑order tensor onto the space of fully symmetric
    elasticity tensors.

    The projection enforces all required symmetries:
        - minor symmetry in the first index pair:  C_ijkl = C_jikl
        - minor symmetry in the second index pair: C_ijkl = C_ijlk
        - major symmetry between index pairs:      C_ijkl = C_klij

    After projection, the tensor satisfies:
        C_ijkl = C_jikl = C_ijlk = C_jilk = C_klij = C_lkij = C_klji = C_lkji

    This guarantees that the tensor has a unique and well‑defined 6×6
    Mandel/Voigt representation, and that the inverse mapping back to
    fourth‑order form is mathematically consistent.

    Mathematically, this operation corresponds to applying the fourth‑order
    symmetrisation operator S, 
        C_sym = S : C : S,
    which is the orthogonal projection onto the symmetry subspace.
*/
Fastor::Tensor<double,3,3,3,3> Foam::Material::castSymmetry(const Fastor::Tensor<double,3,3,3,3>& C_in)
{
    Fastor::Tensor<double,3,3,3,3> C_out(0.0);

    // Mandel index pairs (I -> (i,j)):
    //   0->(0,0)  1->(1,1)  2->(2,2)  3->(1,2)  4->(0,2)  5->(0,1)
    const int vi[6] = {0, 1, 2, 1, 0, 0};
    const int vj[6] = {0, 1, 2, 2, 2, 1};

    // Iterate over the 21 unique upper-triangular Mandel pairs (I <= J).
    // One average is computed per pair from a single fixed ordering of the 8
    // symmetry-related components, then the exact same double value is written
    // to all 8 tensor positions.  This guarantees bitwise equality between all
    // symmetry-equivalent entries so that fourthToMandel() never sees a false
    // symmetry violation.
    //
    // The previous 4-nested loop visited all 81 (i,j,k,l) independently.
    // For a pair such as (i,j,k,l) and (j,i,k,l) the same 8 values were
    // summed in a different order; IEEE 754 floating-point addition is not
    // associative so the averages differed by ~eps*C_max ~ 1e+(n-16) for
    // components of order 1e+n.  For a rotated cubic stiffness tensor with
    // components ~1e10 this produced residuals of ~1e-6, which exceeded the
    // 1e-12 tolerance in fourthToMandel() and caused a spurious fatal error.

    for (int I = 0; I < 6; ++I)
    {
        for (int J = I; J < 6; ++J)
        {
            const int i = vi[I], j = vj[I];
            const int k = vi[J], l = vj[J];

            // Collect all 8 symmetry-related components
            double vals[8] = {
                C_in(i,j,k,l),
                C_in(j,i,k,l),
                C_in(i,j,l,k),
                C_in(j,i,l,k),
                C_in(k,l,i,j),
                C_in(l,k,i,j),
                C_in(k,l,j,i),
                C_in(l,k,j,i)
            };

            // Average non-zero entries
            double sum = 0.0;
            int count = 0;

            for (int n=0; n<8; ++n)
            {
                if (std::abs(vals[n]) > 1e-14)
                {
                    sum += vals[n];
                    count++;
                }
            }

            const double avg = (count > 0 ? sum / count : 0.0);

            // Assign the SAME bitwise value to all 8 symmetry-related positions
            C_out(i,j,k,l) = C_out(j,i,k,l) = C_out(i,j,l,k) = C_out(j,i,l,k) = avg;
            C_out(k,l,i,j) = C_out(l,k,i,j) = C_out(k,l,j,i) = C_out(l,k,j,i) = avg;
        }
    }

    return C_out;
}

/*
    Strictly checks all minor and major symmetries of a rank-four tensor and
    aborts with a FatalError if any are violated beyond 1e-12.

    Call this ONCE at initialisation on the raw input stiffness tensor (before
    any rotation) to catch genuinely non-symmetric user input in GrainProperties.
    Do NOT call it inside the time loop — after castSymmetry() the tensor is
    guaranteed symmetric by construction, so re-checking every cell every step
    is unnecessary overhead.
*/
void Foam::Material::assertSymmetry(const Fastor::Tensor<double,3,3,3,3>& C)
{
    for (int i=0; i<3; ++i)
    {
        for (int j=0; j<3; ++j)
        {
            for (int k=0; k<3; ++k)
            {
                for (int l=0; l<3; ++l)
                {
                    const double Cijkl = C(i,j,k,l);

                    if (std::abs(Cijkl - C(j,i,k,l)) > 1e-12)
                    {
                        FatalErrorInFunction
                            << "Rank 4 tensor violates minor symmetry in (ij): "
                            << "C("<<i<<","<<j<<","<<k<<","<<l<<") != "
                            << "C("<<j<<","<<i<<","<<k<<","<<l<<")" << nl
                            << abort(FatalError);
                    }

                    if (std::abs(Cijkl - C(i,j,l,k)) > 1e-12)
                    {
                        FatalErrorInFunction
                            << "Rank 4 tensor violates minor symmetry in (kl): "
                            << "C("<<i<<","<<j<<","<<k<<","<<l<<") != "
                            << "C("<<i<<","<<j<<","<<l<<","<<k<<")" << nl
                            << abort(FatalError);
                    }

                    if (std::abs(Cijkl - C(k,l,i,j)) > 1e-12)
                    {
                        FatalErrorInFunction
                            << "Rank 4 tensor violates major symmetry: "
                            << "C("<<i<<","<<j<<","<<k<<","<<l<<") != "
                            << "C("<<k<<","<<l<<","<<i<<","<<j<<")" << nl
                            << abort(FatalError);
                    }
                }
            }
        }
    }
}

/*
    Casts a symmetric rank-four tensor to its 6x6 Mandel matrix representation.

    Mandel ordering:
        0 -> 11
        1 -> 22
        2 -> 33
        3 -> 23
        4 -> 13
        5 -> 12

    Shear components are scaled by sqrt(2) to make the basis orthonormal.
*/
Fastor::Tensor<double,6,6> Foam::Material::fourthToMandel(const Fastor::Tensor<double,3,3,3,3>& C, bool verbosity)
{
    if (verbosity)
    {
        double maxMinorIJ = 0.0, maxMinorKL = 0.0, maxMajor = 0.0;

        for (int i=0; i<3; ++i)
        {
            for (int j=0; j<3; ++j)
            {
                for (int k=0; k<3; ++k)
                {
                    for (int l=0; l<3; ++l)
                    {
                        const double Cijkl = C(i,j,k,l);
                        double d;
                        d = std::abs(Cijkl - C(j,i,k,l)); if (d > maxMinorIJ) maxMinorIJ = d;
                        d = std::abs(Cijkl - C(i,j,l,k)); if (d > maxMinorKL) maxMinorKL = d;
                        d = std::abs(Cijkl - C(k,l,i,j)); if (d > maxMajor)   maxMajor   = d;
                    }
                }
            }

        }

        Info << "fourthToMandel symmetry diagnostics:\n"
             << "  max |C(i,j,k,l) - C(j,i,k,l)| = " << maxMinorIJ << nl
             << "  max |C(i,j,k,l) - C(i,j,l,k)| = " << maxMinorKL << nl
             << "  max |C(i,j,k,l) - C(k,l,i,j)| = " << maxMajor   << nl;
    }

    // --- Symmetry safety checks ---
    for (int i=0; i<3; ++i)
    {    
        for (int j=0; j<3; ++j)
        {
            for (int k=0; k<3; ++k)
            {
                for (int l=0; l<3; ++l)
                {
                    const double Cijkl = C(i,j,k,l);

                    if (std::abs(Cijkl - C(j,i,k,l)) > 1e-12)
                    {
                        FatalErrorInFunction
                            << "Rank 4 tensor violates minor symmetry in (ij): "
                            << "C("<<i<<","<<j<<","<<k<<","<<l<<") != "
                            << "C("<<j<<","<<i<<","<<k<<","<<l<<")" << nl
                            << abort(FatalError);
                    }

                    if (std::abs(Cijkl - C(i,j,l,k)) > 1e-12)
                    {
                        FatalErrorInFunction
                            << "Rank 4 tensor violates minor symmetry in (kl): "
                            << "C("<<i<<","<<j<<","<<k<<","<<l<<") != "
                            << "C("<<i<<","<<j<<","<<l<<","<<k<<")" << nl
                            << abort(FatalError);
                    }

                    if (std::abs(Cijkl - C(k,l,i,j)) > 1e-12)
                    {
                        FatalErrorInFunction
                            << "Rank 4 tensor violates major symmetry: "
                            << "C("<<i<<","<<j<<","<<k<<","<<l<<") != "
                            << "C("<<k<<","<<l<<","<<i<<","<<j<<")" << nl
                            << abort(FatalError);
                    }
                }
            }
        }
  
    }

    // --- Mandel mapping using explicit loops ---
    Fastor::Tensor<double,6,6> CM(0.0);

    for (int I=0; I<6; ++I)
    {
        for (int J=0; J<6; ++J)
        {
            double sum = 0.0;

            for (int a=0; a<3; ++a)
            {
                for (int b=0; b<3; ++b)
                {
                    for (int c=0; c<3; ++c)
                    {
                        for (int d=0; d<3; ++d)
                        {
                            sum += M_(I,a,b) * M_(J,c,d) * C(a,b,c,d);
                        }
                    }
                }
            }
 
            CM(I,J) = sum;
        }
    }

    return CM;
}

/*
    Casts a 6x6 Mandel matrix back to a symmetric rank-four tensor.

    Mandel ordering:
        0 -> 11
        1 -> 22
        2 -> 33
        3 -> 23
        4 -> 13
        5 -> 12

    Shear components are scaled by 1/sqrt(2) to undo Mandel scaling.
*/
Fastor::Tensor<double,3,3,3,3> Foam::Material::mandelToFourth(const Fastor::Tensor<double,6,6>& M, bool verbosity)
{
    // --- Symmetry checks on Mandel matrix ---
    for (int I=0; I<6; ++I)
    {
        for (int J=0; J<6; ++J)
        {
            if (std::abs(M(I,J) - M(J,I)) > 1e-12)
            {
                FatalErrorInFunction
                    << "Input Mandel matrix is not symmetric: "
                    << "M("<<I<<","<<J<<") != M("<<J<<","<<I<<")" << nl
                    << abort(FatalError);
            }
        }
    }

    // --- Positive definiteness check ---
    for (int I=0; I<6; ++I)
    {
        if (M(I,I) <= 0.0)
        {
            FatalErrorInFunction
                << "Input Mandel matrix is not positive definite: "
                << "M("<<I<<","<<I<<") <= 0" << nl
                << abort(FatalError);
        }
    }

    // --- Inverse Mandel mapping ---
    Fastor::Tensor<double,3,3,3,3> C(0.0);

    for (int i=0; i<3; ++i)
    {
        for (int j=0; j<3; ++j)
        {
            for (int k=0; k<3; ++k)
            {
                for (int l=0; l<3; ++l)
                {
                    double sum = 0.0;

                    for (int A=0; A<6; ++A)
                    for (int B=0; B<6; ++B)
                    {
                        sum += M_(A,i,j) * M_(B,k,l) * M(A,B);
                    }

                    C(i,j,k,l) = sum;
                }
            }
        }
    }

    return C;
}

Fastor::Tensor<double,3,3,3,3> Foam::Material::invertRankFour(const Fastor::Tensor<double,3,3,3,3>& C, bool verbosity)
{
    if (verbosity)
    {
        Info << "\n==============================================\n";
        Info << "Foam::Material::invertRankFour() called\n";

        printRankFour(C, "C (input)");
        Info << nl;
    }

    // Project to symmetric subspace
    if (verbosity)
    {
        Info << "\tProjecting onto symmetric rank 4 space...\n\n";
    }

    auto Csym = castSymmetry(C);

    if (verbosity)
    {
        printRankFour(Csym, "C_sym");
        Info << nl;
    }

    // Convert to Mandel 6×6
    if (verbosity)
    {
        Info << "\tCasting to 6x6 Mandel representation...\n\n";
    }

    auto M = fourthToMandel(Csym, verbosity);

    if (verbosity)
    {
        printMandelMatrix(M, "M (6x6)");
        Info << nl;
    }

    // Invert the Mandel matrix
    if (verbosity)
    {
        Info << "\tInverting Mandel matrix...\n\n";
    }

    /* 
        Check for singularity of M before inverting.
        Two common causes in GrainProperties:
          (a) C11 == C12 → upper-left 3x3 Mandel block rank-1.
          (b) All shear components zero → lower-right 3x3 Mandel block rank-0.
        Fastor::inverse on a singular matrix produces a FPE or rubbish.
        We detect both cases and return zero compliance (physically correct
        for a stress-free void/crack phase).
    */
    {
        // --- Upper-left 3x3 (normal stiffness) ---
        const double m00 = M(0,0), m01 = M(0,1), m02 = M(0,2);
        const double m10 = M(1,0), m11 = M(1,1), m12 = M(1,2);
        const double m20 = M(2,0), m21 = M(2,1), m22 = M(2,2);
        const double det3 = m00*(m11*m22 - m12*m21)
                          - m01*(m10*m22 - m12*m20)
                          + m02*(m10*m21 - m11*m20);
        const double scale = std::abs(m00) * std::abs(m11) * std::abs(m22);
        const bool normalSingular = (scale > 0.0)
                                  ? (std::abs(det3 / scale) < 1e-6)
                                  : (std::abs(det3) < 1e-30);

        // --- Lower-right 3x3 (shear stiffness) ---
        const double m33 = M(3,3), m34 = M(3,4), m35 = M(3,5);
        const double m43 = M(4,3), m44 = M(4,4), m45 = M(4,5);
        const double m53 = M(5,3), m54 = M(5,4), m55 = M(5,5);
        const double det3s = m33*(m44*m55 - m45*m54)
                           - m34*(m43*m55 - m45*m53)
                           + m35*(m43*m54 - m44*m53);
        const double scales = std::abs(m33) * std::abs(m44) * std::abs(m55);
        const bool shearSingular = (scales > 0.0)
                                 ? (std::abs(det3s / scales) < 1e-6)
                                 : (std::abs(det3s) < 1e-30);

        if (normalSingular || shearSingular)
        {
            WarningInFunction
                << "Mandel matrix is singular ("
                << (normalSingular ? "normal block" : "shear block")
                << " det/scale = "
                << (normalSingular
                    ? (scale  > 0.0 ? det3 /scale  : det3 )
                    : (scales > 0.0 ? det3s/scales : det3s))
                << ").\n"
                << "Common causes: C11==C12 (normal block) or all shear\n"
                << "components zero (shear block) in GrainProperties.\n"
                << "Returning zero compliance tensor for this cell/phase.\n" << nl;
            return Fastor::Tensor<double,3,3,3,3>(0.0);
        }
    }

    auto Minv = Fastor::inverse(M);

    if (verbosity)
    {
        printMandelMatrix(Minv, "M_inv (6x6)");
        Info << nl;
    }

    // Convert back to rank‑4 tensor
    if (verbosity)
    {
        Info << "\tCasting inverse back to rank 4 tensor...\n\n";
    }

    auto Cinv = mandelToFourth(Minv, verbosity);

    if (verbosity)
    {
        printRankFour(Cinv, "C_inv");
        Info << nl;

        Info << "...complete\n";
        Info << "==============================================\n\n";
    }

    return Cinv;
}

double Foam::Material::doubleContractionTwice(const Fastor::Tensor<double, 3, 3>& second, const Fastor::Tensor<double,3,3,3,3>& fourth, bool verbosity)
{
    if (verbosity)
    {
        Info << "\nFoam::Material::doubleContractionTwice() called"
             << endl; 
    }

    double result = 0.0; 

    // compute first contraction - seconds order tensor from (second : fourth), B_{ij} = C_{ijkl} A_{kl}
    Fastor::Tensor<double, 3,3> secondPrime(0.0); 
    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            for (int k = 0; k < 3; k++)
            {
                for (int l = 0; l < 3; l++)
                {
                    secondPrime(i, j) += fourth(i, j, k, l) * second(k, l);
                }
            }
        }
    }

    if (verbosity)
    {
        // padding helper
        auto pad = [](const Foam::string& s, int w) -> Foam::string
        {
            int n = w - s.size();
            if (n > 0)
            {
                std::string spaces(n, ' ');
                return Foam::string(spaces + s);
            }
            return s;
        };

        Info << "\nfirst contraction result:\n\n";

        // Header row (no embedded quotes!)
        Info << pad("", 10);
        for (int j = 0; j < 3; ++j)
        {
            Foam::string label = "[" + Foam::name(j+1) + "]";
            Info << pad(label, 12);
        }
        Info << nl;

        // Body rows
        for (int i = 0; i < 3; ++i)
        {
            Foam::string rowLabel = "[" + Foam::name(i+1) + "]";
            Info << pad(rowLabel, 10);

            for (int j = 0; j < 3; ++j)
            {
                Foam::string val = Foam::name(secondPrime(i,j));
                Info << pad(val, 12);
            }
            Info << nl;
        }    
    }

    // Compute second contraction of secondPrime with second, B_{ij} = C_{ij} A_{ij}
    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            result += second(i, j) * secondPrime(i, j); 
        }
    }

    if (verbosity)
    {
        Info<< "\nsecond contraction result: " << nl
            << "   " << result << nl 
            << "   ... Done" << endl; 
    }

    return result; 
}

bool Foam::Material::isRankFourInverse(const Fastor::Tensor<double,3,3,3,3>& C, Fastor::Tensor<double,3,3,3,3>& S, bool verbosity)
{
    if (verbosity)
    {
        Info << "==============================================\n";
        Info << "Foam::Material::isRankFourInverse() called\n";
        printRankFour(C, "C");
        printRankFour(S, "S");
    }

    // Build 4th‑order symmetric identity
    Fastor::Tensor<double,3,3,3,3> I4(0.0);
    for (int i=0;i<3;i++)
    {
        for (int j=0;j<3;j++)
        {
            for (int m=0;m<3;m++)
            {
                for (int n=0;n<3;n++)
                {
                    I4(i,j,m,n) = 0.5*((i==m && j==n) + (i==n && j==m));
                }
            }
        }
    }

    if (verbosity)
    {
        printRankFour(I4, "I4");
    }

    // Compute R = C:S - I4 using explicit loops
    Fastor::Tensor<double,3,3,3,3> R(0.0);

    for (int i=0;i<3;i++)
    {
        for (int j=0;j<3;j++)
        {
            for (int m=0;m<3;m++)
            {
                for (int n=0;n<3;n++)
                {
                    double sum = 0.0;

                    for (int k=0;k<3;k++)
                    for (int l=0;l<3;l++)
                    {
                        sum += C(i,j,k,l) * S(k,l,m,n);
                    }

                    R(i,j,m,n) = sum - I4(i,j,m,n);
                }
            }
        }
    }

    if (verbosity)
    {
        printRankFour(R, "R");
    }

    // Frobenius norm of residual
    double err = 0.0;
    double normI = 0.0;

    for (int i=0;i<3;i++)
    {
        for (int j=0;j<3;j++)
        {
            for (int m=0;m<3;m++)
            {
                for (int n=0;n<3;n++)
                {
                    err   += R(i,j,m,n) * R(i,j,m,n);
                    normI += I4(i,j,m,n) * I4(i,j,m,n);
                }
            }

        }

    }

    err = std::sqrt(err);
    normI = std::sqrt(normI);
    double relErr = err / normI;

    // Find maximum component error
    double maxErr = 0.0;
    int mi=0, mj=0, mm=0, mn=0;

    for (int i=0;i<3;i++)
    {
        for (int j=0;j<3;j++)
        {
            for (int m=0;m<3;m++)
            {
                for (int n=0;n<3;n++)
                {
                    double a = std::abs(R(i,j,m,n));
                    if (a > maxErr)
                    {
                        maxErr = a;
                        mi=i; mj=j; mm=m; mn=n;
                    }
                }
            }
        }
    }

    if (verbosity)
    {
        Info << "\n\t----------------------------------------------\n";
        Info << "\tInverse Consistency Diagnostics\n";
        Info << "\t----------------------------------------------\n";
        Info << "\tFrobenius norm residual        : " << err << nl;
        Info << "\tRelative residual (||R||/||I||): " << relErr << nl;
        Info << "\tMaximum component error        : " << maxErr
             << "\t at (" << mi << "," << mj << "," << mm << "," << mn << ")\n";
        Info << "\t----------------------------------------------\n\n";
    }

    const double tol = 1e-8;
    bool ok = (relErr < tol);

    if (verbosity)
    {
        if (ok)
        {
            Info << "Result: PASS - C and S are inverses.\n";
        }
        else
        {
            Info << "Result: FAIL - C and S are NOT inverses.\n";
        }
        Info << "==============================================\n\n";
    }

    return ok;
}

bool Foam::Material::isMandelInverse(const Fastor::Tensor<double,3,3,3,3>& C, const Fastor::Tensor<double,3,3,3,3>& S, bool verbosity)
{
    auto M    = fourthToMandel(C);
    auto Minv = fourthToMandel(S);

    // Compute Minv * M and M * Minv, check against I6
    Fastor::Tensor<double,6,6> I(0.0);
    for (int a=0;a<6;++a) I(a,a) = 1.0;

    Fastor::Tensor<double,6,6> P1(0.0), P2(0.0);
    for (int a=0;a<6;++a)
    {
        for (int b=0;b<6;++b)
        {
            double s1 = 0.0, s2 = 0.0;
            for (int c=0;c<6;++c)
            {
                s1 += Minv(a,c)*M(c,b);
                s2 += M(a,c)*Minv(c,b);
            }
            P1(a,b) = s1;
            P2(a,b) = s2;
        }
    }

    double err2 = 0.0;
    for (int a=0;a<6;++a)
    {
        for (int b=0;b<6;++b)
        {
            double d1 = P1(a,b) - I(a,b);
            double d2 = P2(a,b) - I(a,b);
            err2 += d1*d1 + d2*d2;
        }
    }

    double frob = std::sqrt(err2);

    if (verbosity)
        Info << "Mandel inverse Frobenius error = " << frob << nl;

    return frob < 1e-10;
}

void Foam::Material::printRankTwo(const Fastor::Tensor<double,3,3>& A, const word& name) 
{
    Info << "\t" << name << " (3x3):" << nl;
    for (int i=0; i<3; ++i)
    {
        for (int j=0; j<3; ++j)
        {
            Info << "\t\t" << A(i,j) << " ";
        }
        Info << nl;
    }
    Info << nl;
}

void Foam::Material::printRankFour(const Fastor::Tensor<double,3,3,3,3>& A, const word& name)
{
    Info << "\t" << name << " (rank 4 3x3x3x3 tensor):" << nl;

    for (int i=0; i<3; ++i)
    {
        for (int j=0; j<3; ++j)
        {
            Info << "\t\tBlock (" << i << "," << j << "):" << nl;

            for (int k=0; k<3; ++k)
            {
                for (int l=0; l<3; ++l)
                {
                    Info << "\t\t\t" << A(i,j,k,l) << " ";
                }
                Info << nl;
            }

            Info << nl;
        }
    }
}

void Foam::Material::printMandelMatrix(const Fastor::Tensor<double,6,6>& M, const word& name)
{
    Info << "\t" << name << " (6x6 Mandel matrix):\n";

    for (int I=0; I<6; ++I)
    {
        for (int J=0; J<6; ++J)
        {
            Info << "\t\t" << M(I,J) << " ";
        }
        Info << "\n";
    }

    Info << "\n";
}