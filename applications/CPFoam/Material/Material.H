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

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{



 class Material
 {
    private:
    // int N_slips_;
    // double test[6];
    // Polynomial<8>  * _arr;
    
  public:
  string Mat_type;

  scalar rho_0;

  scalar Cv;

  scalar kappa;

  int N_Slips;

  std::vector<vector> SP_norm;
  std::vector<vector> SP_direction;

  double initial_CRSS; // initial critically resolved shear stress

  double tau_S;//saturation value for shear stress

  double h_0;//initial hardening rate

   double r;//latent hardening parameter

  double   k_exponent;//exponent on slip rate term

  double   ref_shear_rate;//exponent on slip rate term

//   Polynomial<8> A;
//   std::vector<Polynomial<8>> Bi;
//   std::vector<Polynomial<8>> Cij;

//   std::vector<Polynomial<8>> eq_xi;

//   std::vector<Polynomial<8>> Qi;

  Material (/*const int*/);
//   Material ();
  Material (const Material&);
  virtual ~Material(){}

  //MEMBER FUNCTIONS
//   int returnComponents();
//   Material Material::operator=
  void Print(label length);

  double deltaG(scalar *mu_arr,scalar *xi_arr,const scalar& T,const label& length);

  tensor RK4(scalar step_size);

};



























// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

