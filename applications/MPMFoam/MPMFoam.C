/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | 
    \\  /    A nd           | 
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

Application
    MPMFoam

Description
    Material Point Method (MPM) solver with crystal plasticity capability
    for large deformation solid mechanics problems.
    
    Features:
    - Large deformation kinematics with multiplicative F = Fe*Fp decomposition
    - Crystal plasticity with slip system evolution
    - Dynamic explicit time integration
    - GIMP/CPDI shape functions (optional)
    
Authors
    Tom Flint

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "pointMesh.H"
#include "pointFields.H"
#include "MaterialPointParcel.H"
#include "Material.H"

#include "IFstream.H"
#include "OFstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Material Point Method solver with crystal plasticity"
    );
    
    // #include "postProcess.H"
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createPointMesh.H"
    #include "initContinuityErrs.H"
    #include "readGravitationalAcceleration.H"
    
    // Create MPM fields
    #include "createMPMFields.H"
    
        // Read material properties
    #include "readMaterialProperties.H"

    #include "readGrainProperties.H"

    // Create material point cloud
    #include "createMPMCloud.H"
    

    
    // #include "CourantNo.H"
    #include "setInitialDeltaT.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting MPM time loop\n" << endl;

    while (runTime.run())
    {
        #include "readDyMControls.H"
        #include "MPMCourantNo.H"
        #include "setDeltaT.H"

        ++runTime;

        Info<< "Time = " << runTime.timeName() << nl << endl;

        scalar dt = runTime.deltaT().value();

        // =====================================================================
        // MPM Algorithm Steps
        // =====================================================================
        
        // Step 1: Reset nodal quantities
        Info<< "Step 1: Resetting nodal quantities" << endl;
        nodalMass = dimensionedScalar("zero", dimMass, 0.0);
        nodalMomentum = dimensionedVector("zero", dimMass*dimVelocity, Zero);
        nodalForce = dimensionedVector("zero", dimForce, Zero);
        nodalStress = dimensionedTensor("zero", dimPressure, Zero);
        
        // Step 2: Particle to Grid transfer (P2G)
        Info<< "Step 2: Particle to Grid transfer" << endl;
        forAllConstIter(materialPointCloud, particles, iter)
        {
            const MaterialPointParcel& p = iter();
            p.transferToGrid
            (
                mesh,
                nodalMass,
                nodalMomentum,
                nodalStress,
                nodalForce
            );
        }
        
        // Apply boundary conditions to nodal fields
        nodalMass.correctBoundaryConditions();
        nodalMomentum.correctBoundaryConditions();
        nodalForce.correctBoundaryConditions();
        



























// // CRITICAL DIAGNOSTIC - Add this to MPMFoam.C after Step 2 (P2G)

// Info<< "=== CRITICAL DIAGNOSTIC: Checking P2G Transfer ===" << endl;

// // 1. Check if particles exist
// label nParticles = particles.size();
// Info<< "  Number of particles: " << nParticles << endl;

// if (nParticles == 0)
// {
//     FatalErrorInFunction
//         << "No particles in cloud! Check particle initialization."
//         << exit(FatalError);
// }

// // 2. Check particle properties
// scalar totalParticleMass = 0;
// scalar maxParticleMass = 0;
// scalar minParticleMass = 1e10;

// forAllConstIter(materialPointCloud, particles, iter)
// {
//     scalar m = iter().mass();
//     totalParticleMass += m;
//     maxParticleMass = max(maxParticleMass, m);
//     minParticleMass = min(minParticleMass, m);
// }

// Info<< "  Total particle mass: " << totalParticleMass << " kg" << endl;
// Info<< "  Particle mass range: [" << minParticleMass << ", " 
//     << maxParticleMass << "] kg" << endl;

// if (totalParticleMass < SMALL)
// {
//     FatalErrorInFunction
//         << "Particles have zero mass! Check particle initialization."
//         << exit(FatalError);
// }

// // 3. Check nodal mass after P2G
// scalar totalNodalMass = 0;
// scalar maxNodalMass = 0;
// label nActiveNodes = 0;

// forAll(nodalMass, pointI)
// {
//     scalar m = nodalMass[pointI];
//     totalNodalMass += m;
//     if (m > SMALL)
//     {
//         nActiveNodes++;
//         maxNodalMass = max(maxNodalMass, m);
//     }
// }

// Info<< "  Total nodal mass (after P2G): " << totalNodalMass << " kg" << endl;
// Info<< "  Active nodes: " << nActiveNodes << " / " << nodalMass.size() << endl;
// Info<< "  Max nodal mass: " << maxNodalMass << " kg" << endl;

// // 4. CRITICAL CHECK
// if (totalNodalMass < SMALL)
// {
//     FatalErrorInFunction
//         << "\n\n"
//         << "========================================\n"
//         << "CRITICAL ERROR: Nodal mass is ZERO!\n"
//         << "========================================\n"
//         << "\n"
//         << "This means transferToGrid() is not working.\n"
//         << "\n"
//         << "Particles exist but aren't mapping to grid.\n"
//         << "\n"
//         << "Possible causes:\n"
//         << "1. transferToGrid() not implemented\n"
//         << "2. Shape functions returning zero\n"
//         << "3. Particles outside mesh bounds\n"
//         << "4. Node IDs incorrect\n"
//         << "\n"
//         << "Total particle mass: " << totalParticleMass << " kg\n"
//         << "Total nodal mass:    " << totalNodalMass << " kg\n"
//         << "\n"
//         << "ACTION: Implement transferToGrid() properly!\n"
//         << "========================================\n"
//         << exit(FatalError);
// }

// // 5. Check mass conservation
// scalar massError = mag(totalNodalMass - totalParticleMass) / totalParticleMass;
// Info<< "  Mass conservation error: " << massError * 100.0 << " %" << endl;

// if (massError > 0.01)
// {
//     WarningInFunction
//         << "Poor mass conservation in P2G transfer! Error = " 
//         << massError * 100.0 << " %" << endl;
// }

// // 6. Check momentum
// scalar totalNodalMomentum = 0;
// forAll(nodalMomentum, pointI)
// {
//     totalNodalMomentum += mag(nodalMomentum[pointI]);
// }
// Info<< "  Total nodal momentum magnitude: " << totalNodalMomentum << " kg⋅m/s" << endl;

// Info<< "===========================================" << endl;















































// Step 3: Compute total nodal forces
Info<< "Step 3: Computing nodal forces" << endl;

// Internal forces (from stress) already added in Step 2 P2G
// Add external body forces
forAll(nodalMass, pointI)
{
    if (nodalMass[pointI] > SMALL)
    {
        // F_external = m*g (gravity)
        nodalForce[pointI] += nodalMass[pointI] * g.value();
        
        // Could add other external forces here:
        // nodalForce[pointI] += tractionForce[pointI];
        // nodalForce[pointI] += dampingForce[pointI];
    }
}

// Step 4: Update nodal momentum (explicit forward Euler)
Info<< "Step 4: Updating nodal momentum (dynamic)" << endl;

forAll(nodalMass, pointI)
{
    if (nodalMass[pointI] > SMALL)
    {
        // CRITICAL: Update momentum first!
        // p^(n+1) = p^n + dt * F
        nodalMomentum[pointI] += dt * nodalForce[pointI];
        
        // Then compute velocity from updated momentum
        // v^(n+1) = p^(n+1) / m
        nodalVelocity[pointI] = nodalMomentum[pointI] / nodalMass[pointI];
        
        // Acceleration for G2P transfer
        // a = F/m
        nodalAcceleration[pointI] = nodalForce[pointI] / nodalMass[pointI];
    }
    else
    {
        nodalMomentum[pointI] = Zero;
        nodalVelocity[pointI] = Zero;
        nodalAcceleration[pointI] = Zero;
    }
}

// Apply boundary conditions AFTER momentum update
#include "applyMPMBoundaryConditions.H"




Info<< "Step 4b: Synchronizing 2D node pairs" << endl;

// Find Front and Back patches
label frontID = boundaryMesh.findPatchID("Front");
label backID = boundaryMesh.findPatchID("Back");

if (frontID != -1 && backID != -1)
{
    const labelList& frontPoints = boundaryMesh[frontID].meshPoints();
    const labelList& backPoints = boundaryMesh[backID].meshPoints();
    
    // For each front point, find corresponding back point
    forAll(frontPoints, i)
    {
        label frontPt = frontPoints[i];
        const vector& frontPos = mesh.points()[frontPt];
        
        // Find back point with same Y-Z
        forAll(backPoints, j)
        {
            label backPt = backPoints[j];
            const vector& backPos = mesh.points()[backPt];
            
            // Check if Y-Z match (X differs)
            if (mag(frontPos.y() - backPos.y()) < 1e-10 &&
                mag(frontPos.z() - backPos.z()) < 1e-10)
            {
                // Average properties between front-back pair
                scalar avgMass = 0.5 * (nodalMass[frontPt] + nodalMass[backPt]);
                vector avgVel = 0.5 * (nodalVelocity[frontPt] + nodalVelocity[backPt]);
                vector avgMom = 0.5 * (nodalMomentum[frontPt] + nodalMomentum[backPt]);
                vector avgAccel = 0.5 * (nodalAcceleration[frontPt] + nodalAcceleration[backPt]);
                
                // Synchronize
                nodalMass[frontPt] = avgMass;
                nodalMass[backPt] = avgMass;
                
                nodalVelocity[frontPt] = avgVel;
                nodalVelocity[backPt] = avgVel;
                
                nodalMomentum[frontPt] = avgMom;
                nodalMomentum[backPt] = avgMom;
                
                nodalAcceleration[frontPt] = avgAccel;
                nodalAcceleration[backPt] = avgAccel;
                
                break;
            }
        }
    }
}


        
        // Step 5: Update particle kinematics
        Info<< "Step 5: Updating particle kinematics" << endl;
        forAllIter(materialPointCloud, particles, iter)
        {
            MaterialPointParcel& p = iter();
            
            // Update deformation gradient
            p.updateDeformationGradient(mesh, nodalVelocity, dt);
            
            // Grid to Particle velocity transfer
            p.transferFromGrid(mesh, nodalVelocity, nodalAcceleration);
        }
        
        // Step 6: Advect particles
        Info<< "Step 6: Advecting particles" << endl;
        particle::trackingData td(particles);
        
        forAllIter(materialPointCloud, particles, iter)
        {
            MaterialPointParcel& p = iter();
            p.move(particles, td, dt);
        }
        
        // Step 7: Constitutive update (stress computation)
        Info<< "Step 7: Updating particle stresses" << endl;
        forAllIter(materialPointCloud, particles, iter)
        {
            MaterialPointParcel& p = iter();
            
            // Get material properties for this particle
            const dictionary& matProps = materialProperties[p.materialID()];
            
            // Perform constitutive update (includes crystal plasticity)
            p.constitutiveUpdate(dt, matProps);
            // p.computeElasticStress(matProps);//just elastic for debugging
        }
        
        // Step 8: Update history variables and check failure
        Info<< "Step 8: Checking material state" << endl;
        label nFailed = 0;
        scalar maxEqPlasticStrain = 0.0;
        
        forAllIter(materialPointCloud, particles, iter)
        {
            MaterialPointParcel& p = iter();
            
            // Check for excessive deformation
            scalar J = det(p.F());
            if (J < 0.1 || J > 10.0)
            {
                nFailed++;
                WarningInFunction
                    << "Particle " << iter().origId() 
                    << " has excessive deformation: J = " << J << endl;
            }
            
            // Track maximum plastic strain
            // (accessing private member would need getter method)
            // maxEqPlasticStrain = max(maxEqPlasticStrain, p.epsilonP_eq());
        }
        
        if (nFailed > 0)
        {
            Info<< "Warning: " << nFailed << " particles have excessive deformation" << endl;
        }
        
        // Write particle cloud
if (runTime.outputTime())
{
    Info<< "Writing particle cloud..." << endl;
    particles.write();
}

#include "writeMPMFields.H"
runTime.write();
        
        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }
    
    Info<< "End\n" << endl;

    return 0;
}

// ************************************************************************* //