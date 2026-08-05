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

Class
    Foam::multiphaseMixture

Description
    Dirty mess of functions and junk for solid dynamics and tensor maffs


SourceFiles
    multiphaseMixture.C

\*---------------------------------------------------------------------------*/

#ifndef material_H
#define material_H

#include "Interface.H"

#include <vector>
#include "Polynomial.H"
#include "IOdictionary.H"
#include "PtrDictionary.H"
#include "volFields.H"
#include "surfaceFields.H"

#include "Fastor.h"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
 class Material
 {
    private:
        // int N_slips_;
        // double test[6];
        // Polynomial<8>  * _arr;

        // Mandel tensor for converting between symmetric rank-four and Mandel representations
        static const Fastor::Tensor<double,6,3,3> M_;
    public:
        string Mat_type;
        scalar rho_0  = 0.0;
        scalar Cv     = 0.0;
        scalar kappa  = 0.0;
        int N_Slips   = 0;

        // Mie-Gruneisen EOS parameters 
        scalar c0        = 0.0;   // bulk sound speed at reference   [m/s]
        scalar s_UsUp    = 0.0;   // linear Us-up Hugoniot slope     [-]
        scalar Gamma0    = 0.0;   // Gruneisen parameter             [-]
        scalar cL0_bound = 0.0;   // derived CFL wave-speed bound    [m/s]

        std::vector<vector> SP_norm;
        std::vector<vector> SP_direction;

        double initial_CRSS = 0.0; // initial critically resolved shear stress
        double tau_S = 0.0;//saturation value for shear stress
        double h_0 = 0.0;//initial hardening rate
        double r = 0.0;//latent hardening parameter
        double   k_exponent = 0.0;//exponent on slip rate term
        double   ref_shear_rate = 0.0;//exponent on slip rate term

        // Polynomial<8> A;
        // std::vector<Polynomial<8>> Bi;
        // std::vector<Polynomial<8>> Cij;

        // std::vector<Polynomial<8>> eq_xi;

        // std::vector<Polynomial<8>> Qi;


        Material() = default;
        Material(const Material&) = default;
        Material& operator=(const Material&) = default;
        virtual ~Material(){}

        //MEMBER FUNCTIONS
        // int returnComponents();
        // Material Material::operator=

        void Print(label length);
        double deltaG(scalar *mu_arr,scalar *xi_arr,const scalar& T,const label& length);
        tensor RK4(scalar step_size);

        // Bespoke tensor functions
        static Fastor::Tensor<double,3,3,3,3> invertRankFour(const Fastor::Tensor<double,3,3,3,3>& C, bool verbosity = false);
        
        // Helper functions for inverting rank four ... 
        static Fastor::Tensor<double,3,3,3,3> castSymmetry(const Fastor::Tensor<double,3,3,3,3>& C_in);
        static void assertSymmetry(const Fastor::Tensor<double,3,3,3,3>& C);
        static Fastor::Tensor<double,6,6> fourthToMandel(const Fastor::Tensor<double,3,3,3,3>& C, bool verbosity = false);
        static Fastor::Tensor<double,3,3,3,3> mandelToFourth(const Fastor::Tensor<double, 6, 6>& M, bool verbosity = false);
        static bool isRankFourInverse(const Fastor::Tensor<double,3,3,3,3>& C, Fastor::Tensor<double,3,3,3,3>& S, bool verbosity = false);
        static bool isMandelInverse(const Fastor::Tensor<double,3,3,3,3>& C, const Fastor::Tensor<double,3,3,3,3>& S, bool verbosity = false);
        static double doubleContractionTwice(const Fastor::Tensor<double, 3, 3>& second, const Fastor::Tensor<double,3,3,3,3>& fourth, bool verbosity = false); 

        // Tensor utility functions 
        static void printRankTwo(const Fastor::Tensor<double,3,3>& A, const word& name="Rank 2 Tensor");
        static void printRankFour(const Fastor::Tensor<double,3,3,3,3>& A, const word& name="Rank 4 Tensor");
        static void printMandelMatrix(const Fastor::Tensor<double,6,6>& M, const word& name="Mandel Matrix");

        static Fastor::Tensor<double,3,3> symm(const Fastor::Tensor<double,3,3>& A);
};



























// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //