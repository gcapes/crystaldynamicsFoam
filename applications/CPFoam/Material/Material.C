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



// int Foam::Material::returnComponents(){

//     return N_slips_;
// }



void Foam::Material::Print(label length){


}


double Foam::Material::deltaG( scalar *mu_arr,scalar *xi_arr,const scalar& T,const label& length){

// PhaseComposition[count_phase][count_i][cell];

double G_Quad=0.0;



    return G_Quad;
}

tensor Foam::Material::RK4(scalar step_size){

    return tensor::zero;
}


Fastor::Tensor<double,3,3,3,3> Foam::Material::invertRankFour(const Fastor::Tensor<double,3,3,3,3>& inp, bool verbosity) 
{
    if (verbosity)
    {
        Info << "Foam::Material::invertRankFour() called"
             << endl; 
    }

    // Map (i,j) to Voigt index alpha (1..6)
    static constexpr int pair_to_voigt[4][4] = {
    /* 0 */ {0, 0, 0, 0},
    /* 1 */ {0, 1, 6, 5}, // 11, 12, 13
    /* 2 */ {0, 6, 2, 4}, // 21, 22, 23
    /* 3 */ {0, 5, 4, 3}  // 31, 32, 33
    };

    // In pure Voigt 
    static constexpr double voigt_scale[7] = {
        0.0,
        1.0, // 1 -> 11
        1.0, // 2 -> 22
        1.0, // 3 -> 33
        1.0, // 4 -> 23
        1.0, // 5 -> 13
        1.0  // 6 -> 12
    };

    // Map input tensor to 6x6 voigt notation ... 
    Fastor::Tensor<double, 6, 6> C2(0.0);

    for (int i = 0; i < 3; ++i) 
    {
        for (int j = 0; j < 3; ++j) 
        {
            for (int k = 0; k < 3; ++k) 
            {
                for (int l = 0; l < 3; ++l) 
                {

                    int alpha = pair_to_voigt[i+1][j+1];
                    int beta  = pair_to_voigt[k+1][l+1];

                    if (alpha == 0 || beta == 0) continue;

                    double s_alpha = voigt_scale[alpha];
                    double s_beta  = voigt_scale[beta];

                    if (alpha != beta) 
                    {
                        C2(alpha-1, beta-1) += s_alpha * s_beta * inp(i,j,k,l);
                        C2(beta-1, alpha-1) += s_alpha * s_beta * inp(i,j,k,l);
                    }
                    else 
                    {
                        C2(alpha-1, beta-1) = s_alpha * s_beta * inp(i,j,k,l);
                    }
                }
            }
        }
    }

    // Invert Voigt stiffness -> S2
    Fastor::Tensor<double, 6, 6> S2;
    S2 = Fastor::inverse(C2);


    // Convert Voigt S2 back to 4th order

    Fastor::Tensor<double,3,3,3,3> out(0.0);

    for (int i = 0; i < 3; ++i) 
    {
        for (int j = 0; j < 3; ++j) 
        {
            for (int k = 0; k < 3; ++k) 
            {
                for (int l = 0; l < 3; ++l) 
                {

                    int alpha = pair_to_voigt[i+1][j+1];
                    int beta  = pair_to_voigt[k+1][l+1];

                    if (alpha == 0 || beta == 0) continue;

                    double s_alpha = voigt_scale[alpha];
                    double s_beta  = voigt_scale[beta];

                    double val = S2(alpha-1, beta-1) / (s_alpha * s_beta);

                    // Symmetry enforcement ... 
                    out(i,j,k,l) = val;
                    out(i,j,l,k) = val;
                    out(j,i,k,l) = val;
                    out(j,i,l,k) = val;
                }
            }
        }
    }

    if (verbosity)
    {
        Fastor::Tensor<double,6,6> R2 = 0.0;

        for (int a = 0; a < 6; ++a)
        for (int b = 0; b < 6; ++b)
        for (int c = 0; c < 6; ++c)
            R2(a,b) += C2(a,c) * S2(c,b);

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

        Info << "\nContraction R2 = C2 * S2:\n\n";

        // Header row (no embedded quotes!)
        Info << pad("", 10);
        for (int j = 0; j < 6; ++j)
        {
            Foam::string label = "[" + Foam::name(j+1) + "]";
            Info << pad(label, 12);
        }
        Info << nl;

        // Body rows
        for (int i = 0; i < 6; ++i)
        {
            Foam::string rowLabel = "[" + Foam::name(i+1) + "]";
            Info << pad(rowLabel, 10);

            for (int j = 0; j < 6; ++j)
            {
                Foam::string val = Foam::name(R2(i,j));
                Info << pad(val, 12);
            }
            Info << nl;
        }

        Info << "   ...Done" << endl;
    }


    return out; 
}




