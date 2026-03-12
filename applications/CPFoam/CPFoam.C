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
#include "dynamicFvMesh.H"
#include "EulerDdtScheme.H"
#include "localEulerDdtScheme.H"
#include "CrankNicolsonDdtScheme.H"
#include "subCycle.H"
//#include "pisoControl.H"
#include "pimpleControl.H"
#include "fvOptions.H"
#include "CorrectPhi.H"

#include "Material.H"


#include "Fastor.h" 
#include <random>
#include <math.h>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{

        argList::addNote
    (
        "Solver for"
    );

    
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

    while (runTime.run())
    {

        // int iCorr = 0;
        scalar initialResidual = 0;

        // int nCorr = 20;
        // scalar convergenceTolerance(1e-9);

        scalar maxTimeStep(readScalar(controlDict.lookup("maxDeltaT")));



        #include "readDyMControls.H"
        #include "CourantNo.H"
        #include "setDeltaT.H"



        ++runTime;

        Info<< "Time = " << runTime.timeName() << nl << endl;

        

        // volTensorField sigtemp(sigma);
          for (label i = 0; i < N_Ori; i++)
            {
                sigtemp[i]=sigma[i];
                sigtemp[i].correctBoundaryConditions();
            }
            sigmaGlobaltemp = sigmaGlobal;
            sigmaGlobaltemp.correctBoundaryConditions();
        label pimpleIter = 0;
        while (pimple.loop())//TOM - UEqn DEFINITELY goes inside the loop - need to take update rotation out of PIMPLE loop
        {
            pimpleIter++;

            volVectorField U_old = U;
            U_old.correctBoundaryConditions();

            
            #include "UEqn.H"
        // #include "setU.H"

            

            #include "UpdateProps.H"
            #include "StressUpdate.H"

        }

            #include "PF.H"
            #include "TEqn.H"
            #include "rhoEqn.H"
            

        #include "Advect.H"
        





        

        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }
    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
