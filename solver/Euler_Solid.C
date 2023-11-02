/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2021 OpenFOAM Foundation
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

Application
    Euler_Solid_MechFoam

Description
    
    
Authors
    Tom Flint
    Pratheek Shanthraj
    Philip Cardiff

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "EulerDdtScheme.H"
#include "localEulerDdtScheme.H"
#include "CrankNicolsonDdtScheme.H"
#include "subCycle.H"
//#include "pisoControl.H"
#include "pimpleControl.H"
#include "fvModels.H"
#include "fvConstraints.H"


#include <Fastor/Fastor.h> 
#include <random>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
#include "postProcess.H"

    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createMesh.H"

   #include "createDyMControls.H"

    
    #include "initContinuityErrs.H"
    #include "createFields.H"
    #include "CourantNo.H"

    #include "setInitialDeltaT.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    while (pimple.run(runTime))
    {
        #include "readDyMControls.H"


            #include "CourantNo.H"
         //   #include "alphaCourantNo.H"
       //  #include "setDeltaT.H"


        runTime++;

        Info<< "Time = " << runTime.userTimeName() << nl << endl;

        // --- Pressure-velocity PIMPLE corrector loop
        

        



sigma.correctBoundaryConditions();

        while (pimple.loop())
        {
         

            #include "UpdateProps.H"
            
   
            

            #include "UEqn.H"
            #include "rhoEqn.H"

            // #include "TEqn.H"

            

            

        //     while (pimple.correct())
        //     {

        //     // if (pimple.finalNonOrthogonalIter())
        // // {
        //     // #include "rhoEqn.H"
        // // }

        //     }

        }

        


        #include "TEqn.H"

// forAll(mesh.C(), cell)
// {
//         if(Hardening){

//   double log_term = 0.0;
//   if(epsilon_dot_P_VM[cell]>SMALL){
//     log_term=Foam::log(epsilon_dot_P_VM[cell]);
//   }

// sigma_Y[cell]=(sigma_Y_0.value()+(K_0*Foam::pow(Acc_Strain[cell],K_n)));//*(K_B*Foam::exp(-(K_Beta_0-K_Beta_1*log_term)*T[cell]));
// // Info<<sigma_Y[cell]<<endl;


// }
// else{
//   sigma_Y[cell]=sigma_Y_0.value();
// }
// }
        

        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }
    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
