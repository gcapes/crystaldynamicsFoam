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

    // Pre-populate C_rotated when stress solve is bypassed (stressUpdate=false).
    if (!stressUpdate)
    {
        Info<< "stressUpdate=false: initialising C_rotated\n" << endl;
        for (label h = 0; h < N_Ori; h++)
        {
            // Use cell 0 rotation angles (uniform per phase).
            const double p1 = Rot1[h][0];
            const double p2 = Rot2[h][0];
            const double p3 = Rot3[h][0];
            const double c1 = Foam::cos(p1), s1 = Foam::sin(p1);
            const double c2 = Foam::cos(p2), s2 = Foam::sin(p2);
            const double c3 = Foam::cos(p3), s3 = Foam::sin(p3);

            Fastor::Tensor<double,3,3> R(0.0);
            R(0,0) =  c2*c3;
            R(0,1) =  s1*s2*c3 - c1*s3;
            R(0,2) =  c1*s2*c3 + s1*s3;
            R(1,0) =  c2*s3;
            R(1,1) =  s1*s2*s3 + c1*c3;
            R(1,2) =  c1*s2*s3 - s1*c3;
            R(2,0) = -s2;
            R(2,1) =  s1*c2;
            R(2,2) =  c1*c2;

            Fastor::Tensor<double,3,3,3,3> C_rot(0.0);
            for (int i = 0; i < 3; i++)
                for (int j = 0; j < 3; j++)
                    for (int k = 0; k < 3; k++)
                        for (int l = 0; l < 3; l++)
                        {
                            double val = 0.0;
                            for (int m = 0; m < 3; m++)
                                for (int n = 0; n < 3; n++)
                                    for (int p = 0; p < 3; p++)
                                        for (int q = 0; q < 3; q++)
                                            val += R(i,m)*R(j,n)*R(k,p)*R(l,q)
                                                 * Cprime2[Phase2Material[h]](m,n,p,q);
                            C_rot(i,j,k,l) = val;
                        }

            for (label cell = 0; cell < mesh.nCells(); cell++)
                C_rotated[h][cell] = C_rot;
        }
    }

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

        

        for (label i = 0; i < N_Ori; i++)
        {
            sigtemp[i]=sigma[i];
            sigtemp[i].correctBoundaryConditions();
        }
        sigmaGlobaltemp = sigmaGlobal;
        sigmaGlobaltemp.correctBoundaryConditions();
        label pimpleIter = 0;
        while (pimple.loop())
        {
            pimpleIter++;

            volVectorField U_old = U;
            U_old.correctBoundaryConditions();
            
            #include "rhoEqn.H"

            #include "UEqn.H"

            #include "UpdateProps.H"
            #include "StressUpdate.H"
            #include "TEqn.H"

        }

            #include "PF.H"
            //#include "TEqn.H"
            //#include "rhoEqn.H"
            

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
