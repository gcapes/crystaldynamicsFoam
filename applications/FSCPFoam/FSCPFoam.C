/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           |
     \\/     M anipulation  |
-------------------------------------------------------------------------------
Application
    FSCPFoam

Description
    Finite Strain Eulerian Crystal Plasticity Solver.
    
    Implements multiplicative decomposition F = Fe·Fp with:
    - Global deformation gradient F advected through Eulerian grid
    - Per-phase plastic deformation gradient Fp advected with plastic source
    - 2nd Piola-Kirchhoff stress S = C:Ee in intermediate configuration
    - Push-forward to Cauchy stress σ = (1/J_e) Fe·S·Fe^T
    - Density from ρ = ρ_ref / det(F), where ρ_ref is phase-weighted
    
    Crystal orientation is encoded in Fp (initialized from Euler angles).
    No separate rotation fields are needed.
    
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
#include "pimpleControl.H"
#include "fvOptions.H"
#include "CorrectPhi.H"

#include "Material.H"
#include "FastorOpenFOAMUtils.H"

#include "Fastor.h"
#include <random>
#include <math.h>

#include "polarDecomposition.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Finite Strain Eulerian Crystal Plasticity Solver"
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

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info << "\nStarting time loop\n" << endl;

    while (runTime.run())
    {
        scalar initialResidual = 0;
        scalar maxTimeStep(readScalar(controlDict.lookup("maxDeltaT")));

        #include "readDyMControls.H"
        #include "CourantNo.H"
        #include "setDeltaT.H"

        ++runTime;

        Info << "Time = " << runTime.timeName() << nl << endl;

        // Store previous-iteration stress for PIMPLE relaxation
        for (label i = 0; i < N_Ori; i++)
        {
            sigtemp[i] = sigma[i];
            sigtemp[i].correctBoundaryConditions();
        }
        sigmaGlobaltemp = sigmaGlobal;
        sigmaGlobaltemp.correctBoundaryConditions();

        // ---- PIMPLE loop: momentum + constitutive update ----
        label pimpleIter = 0;
        while (pimple.loop())
        {
            pimpleIter++;

            volVectorField U_old = U;
            U_old.correctBoundaryConditions();

            #include "UEqn.H"

            // Velocity relaxation (non-final iterations)
            if (!pimple.finalIter())
            {
                U = pimple_relax_U * U + (1.0 - pimple_relax_U) * U_old;
                U.correctBoundaryConditions();
            }

            #include "UpdateProps.H"
            #include "StressUpdate.H"
        }

        // ---- Post-PIMPLE updates (once per time step) ----
        #include "PF.H"
        // #include "TEqn.H"
        #include "Advect.H"

        runTime.write();

        Info << "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
             << "  ClockTime = " << runTime.elapsedClockTime() << " s"
             << nl << endl;
    }

    Info << "End\n" << endl;
    return 0;
}

// ************************************************************************* //