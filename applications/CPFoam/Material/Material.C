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


Foam::Material::Material (/*const int N_slips*/)
{
// N_slips_=N_slips;
}



Foam::Material::Material (const Material& M1)
{
// Bi = new Polynomial<8>[M1.N_comp_];
// *Bi=*M1.Bi; 
}

// Foam::Material::Material ()
// {
// }



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

    for (int i=0; i<3; ++i)
    for (int j=0; j<3; ++j)
    for (int k=0; k<3; ++k)
    for (int l=0; l<3; ++l)
    {
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

        double avg = (count > 0 ? sum / count : 0.0);

        // Assign the averaged value
        C_out(i,j,k,l) = avg;
    }

    return C_out;
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
Fastor::Tensor<double,6,6> Foam::Material::fourthToMandel(const Fastor::Tensor<double,3,3,3,3>& C)
{
    // --- Symmetry checks (unchanged) ---
    for (int i=0; i<3; ++i)
    for (int j=0; j<3; ++j)
    for (int k=0; k<3; ++k)
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

    // --- Mandel mapping using explicit loops ---
    Fastor::Tensor<double,6,6> CM(0.0);

    for (int I=0; I<6; ++I)
    for (int J=0; J<6; ++J)
    {
        double sum = 0.0;

        for (int a=0; a<3; ++a)
        for (int b=0; b<3; ++b)
        for (int c=0; c<3; ++c)
        for (int d=0; d<3; ++d)
        {
            sum += M_(I,a,b) * M_(J,c,d) * C(a,b,c,d);
        }

        CM(I,J) = sum;
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
Fastor::Tensor<double,3,3,3,3> Foam::Material::mandelToFourth(const Fastor::Tensor<double,6,6>& M)
{
    // --- Symmetry checks on Mandel matrix ---
    for (int I=0; I<6; ++I)
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

    // --- Positive definiteness check (simple diagonal test) ---
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
    for (int j=0; j<3; ++j)
    for (int k=0; k<3; ++k)
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

    // 2. Convert to Mandel 6×6
    if (verbosity)
    {
        Info << "\tCasting to 6x6 Mandel representation...\n\n";
    }

    auto M = fourthToMandel(Csym);

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

    auto Cinv = mandelToFourth(Minv);

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
    for (int j=0;j<3;j++)
    for (int m=0;m<3;m++)
    for (int n=0;n<3;n++)
    {
        I4(i,j,m,n) = 0.5*((i==m && j==n) + (i==n && j==m));
    }

    if (verbosity)
    {
        printRankFour(I4, "I4");
    }

    // Compute R = C:S - I4 using explicit loops
    Fastor::Tensor<double,3,3,3,3> R(0.0);

    for (int i=0;i<3;i++)
    for (int j=0;j<3;j++)
    for (int m=0;m<3;m++)
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

    if (verbosity)
    {
        printRankFour(R, "R");
    }

    // Frobenius norm of residual
    double err = 0.0;
    double normI = 0.0;

    for (int i=0;i<3;i++)
    for (int j=0;j<3;j++)
    for (int m=0;m<3;m++)
    for (int n=0;n<3;n++)
    {
        err   += R(i,j,m,n) * R(i,j,m,n);
        normI += I4(i,j,m,n) * I4(i,j,m,n);
    }

    err = std::sqrt(err);
    normI = std::sqrt(normI);
    double relErr = err / normI;

    // Find maximum component error
    double maxErr = 0.0;
    int mi=0, mj=0, mm=0, mn=0;

    for (int i=0;i<3;i++)
    for (int j=0;j<3;j++)
    for (int m=0;m<3;m++)
    for (int n=0;n<3;n++)
    {
        double a = std::abs(R(i,j,m,n));
        if (a > maxErr)
        {
            maxErr = a;
            mi=i; mj=j; mm=m; mn=n;
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

    double err2 = 0.0;
    for (int a=0;a<6;++a)
    for (int b=0;b<6;++b)
    {
        double d1 = P1(a,b) - I(a,b);
        double d2 = P2(a,b) - I(a,b);
        err2 += d1*d1 + d2*d2;
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
