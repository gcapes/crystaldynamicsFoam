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

#include <vector>
#include "Polynomial.H"
#include "IOdictionary.H"
#include "PtrDictionary.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "Fastor/Fastor.h"

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


















class Interface
 :public IOdictionary
 {

    public:

    class interfacePair
        :
        public Pair<word>
    {
    public:

        class hash
            :
            public Hash<interfacePair>
        {
        public:

            hash()
            {}

            label operator()(const interfacePair& key) const
            {
                return word::hash()(key.first()) + word::hash()(key.second());
            }
        };


        // Constructors

        interfacePair()
        {}

        interfacePair(const word& alpha1Name, const word& alpha2Name)
            :
            Pair<word>(alpha1Name, alpha2Name)
        {}

        // interfacePair(const phase& alpha1, const phase& alpha2)
        // :
        //     Pair<word>(alpha1.name(), alpha2.name())
        // {}


        // Friend Operators

        friend bool operator==
        (
            const interfacePair& a,
            const interfacePair& b
        )
        {
            return
                (
                    ((a.first() == b.first()) && (a.second() == b.second()))
                );
        }

        friend bool operator!=
        (
            const interfacePair& a,
            const interfacePair& b
        )
        {
            return (!(a == b));
        }
    };
    
    private:
    const fvMesh& mesh_;
    const volScalarField& T_;
    word Interface_type;
  public:
    Interface (word,const volScalarField&);
    // Material (int,int);

    virtual ~Interface()
    {}

    //Member Functions
    void Assign_Energies(scalar **NRG_arr);
    void Assign_Mobilities(scalar **NRG_arr);
    void Assign_Anisotropy_Strength(scalar **Aniso_strength_arr);
    scalar aniso_energy(scalar *tensor_i, scalar *tensor_j, scalar *vec_i, scalar *vec_j, const scalar& phi_i, const scalar& phi_j,const scalar& int_width, const scalar& deltaN);
    // int area (void) {return (width*height);}
};











// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //