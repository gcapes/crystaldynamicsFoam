/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | 
    \\  /    A nd           | 
     \\/     M anipulation  |
-------------------------------------------------------------------------------
Application
    diffuseHardBoundaries

Description
    Smooths the hard grain boundaries of the n.* fields by iterated face averaging
    and renormalises the phase sum.

Author
    Generated for CP phase field simulations

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Diffuse hard boundaries in phase fields using face averaging"
    );
    
    argList::addOption
    (
        "nIterations",
        "int",
        "number of smoothing iterations (default: 3)"
    );
    
    argList::addOption
    (
        "relaxation",
        "scalar",
        "relaxation factor for smoothing (default: 0.5)"
    );
    
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createMesh.H"
    
    // Read parameters
    
    label nIterations = 3;
    if (args.optionFound("nIterations"))
    {
        nIterations = args.optionRead<label>("nIterations");
    }
    
    scalar relaxation = 0.5;
    if (args.optionFound("relaxation"))
    {
        relaxation = args.optionRead<scalar>("relaxation");
    }
    
    Info<< "\nDiffusion parameters:" << nl
        << "  Number of iterations: " << nIterations << nl
        << "  Relaxation factor: " << relaxation << nl
        << endl;
    
    // Number of phases from GrainProperties
    
    IOdictionary grainProperties
    (
        IOobject
        (
            "GrainProperties",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );
    
    label N_Phases = readLabel(grainProperties.lookup("N_Phases"));
    Info<< "Number of phases: " << N_Phases << endl;
    
    // Read the phase fields
    
    Info<< "\nReading phase fields..." << endl;
    
    PtrList<volScalarField> phaseFields(N_Phases);
    
    forAll(phaseFields, i)
    {
        word fieldName = "n." + Foam::name(i);
        
        Info<< "  Reading " << fieldName << endl;
        
        phaseFields.set
        (
            i,
            new volScalarField
            (
                IOobject
                (
                    fieldName,
                    runTime.timeName(),
                    mesh,
                    IOobject::MUST_READ,
                    IOobject::AUTO_WRITE
                ),
                mesh
            )
        );
    }
    
    // Boundary diffusion
    
    Info<< "\nApplying boundary diffusion..." << endl;
    
    // Create temporary fields to store smoothed values
    PtrList<volScalarField> smoothedFields(N_Phases);
    
    forAll(smoothedFields, i)
    {
        smoothedFields.set
        (
            i,
            new volScalarField
            (
                IOobject
                (
                    "smoothed_" + Foam::name(i),
                    runTime.timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh,
                dimensionedScalar("zero", dimless, 0.0)
            )
        );
    }
    
    // Iterative smoothing
    for (label iter = 0; iter < nIterations; iter++)
    {
        Info<< "  Iteration " << iter + 1 << " / " << nIterations << endl;
        
        forAll(phaseFields, i)
        {
            volScalarField& phi = phaseFields[i];
            volScalarField& smoothed = smoothedFields[i];
            
            // face average of the interpolated field
            surfaceScalarField phiF = fvc::interpolate(phi);
            smoothed = fvc::average(phiF);
            
            // Apply relaxation to avoid over-smoothing
            phi = (1.0 - relaxation) * phi + relaxation * smoothed;
            
            // Ensure boundaries are treated correctly
            phi.correctBoundaryConditions();
        }
        
        // Renormalize to ensure sum of phases = 1.0 at each point
        volScalarField sumPhases
        (
            IOobject
            (
                "sumPhases",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedScalar("zero", dimless, 0.0)
        );
        
        forAll(phaseFields, i)
        {
            sumPhases += phaseFields[i];
        }
        
        // Normalize each field
        forAll(phaseFields, i)
        {
            forAll(mesh.C(), cellI)
            {
                if (sumPhases[cellI] > SMALL)
                {
                    phaseFields[i][cellI] /= sumPhases[cellI];
                }
            }
            
            // Also normalize boundary fields
            forAll(sumPhases.boundaryField(), patchI)
            {
                if (!sumPhases.boundaryField()[patchI].coupled()
                 && sumPhases.boundaryField()[patchI].size() > 0)
                {
                    forAll(sumPhases.boundaryField()[patchI], faceI)
                    {
                        if (sumPhases.boundaryField()[patchI][faceI] > SMALL)
                        {
                            phaseFields[i].boundaryFieldRef()[patchI][faceI] /= 
                                sumPhases.boundaryField()[patchI][faceI];
                        }
                    }
                }
            }
            
            phaseFields[i].correctBoundaryConditions();
        }
    }
    
    // Grain-boundary indicator
    
    Info<< "\nCalculating grain boundary indicator..." << endl;
    
    volScalarField GBREF
    (
        IOobject
        (
            "GBREF",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("zero", dimless, 0.0)
    );
    
    // GBREF = sum_{i<j} phi_i phi_j: large where phases coexist
    forAll(mesh.C(), cellI)
    {
        scalar gb = 0.0;
        for (label i = 0; i < N_Phases; i++)
        {
            for (label j = i+1; j < N_Phases; j++)
            {
                gb += phaseFields[i][cellI] * phaseFields[j][cellI];
            }
        }
        GBREF[cellI] = gb;
    }
    
    GBREF.correctBoundaryConditions();
    
    // Write
    
    Info<< "\nWriting diffused phase fields..." << endl;
    
    forAll(phaseFields, i)
    {
        phaseFields[i].write();
        Info<< "  Written " << phaseFields[i].name() << endl;
    }
    
    GBREF.write();
    Info<< "  Written GBREF (grain boundary indicator)" << endl;
    
    // Statistics
    
    Info<< "\nPhase field statistics:" << endl;
    
    volScalarField sumPhases
    (
        IOobject
        (
            "sumPhases",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("zero", dimless, 0.0)
    );
    
    forAll(phaseFields, i)
    {
        sumPhases += phaseFields[i];
        
        scalar minVal = gMin(phaseFields[i]);
        scalar maxVal = gMax(phaseFields[i]);
        scalar avgVal = gAverage(phaseFields[i]);
        
        Info<< "  " << phaseFields[i].name() << ": "
            << "min = " << minVal << ", "
            << "max = " << maxVal << ", "
            << "avg = " << avgVal << endl;
    }
    
    scalar minSum = gMin(sumPhases);
    scalar maxSum = gMax(sumPhases);
    
    Info<< "\nSum of phases: min = " << minSum 
        << ", max = " << maxSum << endl;
    
    if (mag(minSum - 1.0) > 1e-6 || mag(maxSum - 1.0) > 1e-6)
    {
        WarningInFunction
            << "Sum of phases deviates from 1.0" << nl
            << "  This may indicate numerical issues" << endl;
    }
    
    scalar maxGB = gMax(GBREF);
    scalar avgGB = gAverage(GBREF);
    
    Info<< "\nGrain boundary indicator (GBREF): "
        << "max = " << maxGB << ", avg = " << avgGB << endl;
    
    Info<< "\nDone!" << endl;
    
    return 0;
}

// ************************************************************************* //
