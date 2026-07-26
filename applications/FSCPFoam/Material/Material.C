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

#include "Fastor/Fastor.h"
#include "FastorOpenFOAMUtils.H"
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