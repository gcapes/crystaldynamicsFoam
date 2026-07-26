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
#include "pointBoundaryMesh.H"
#include "MaterialPointParcel.H"
#include "Material.H"

#include "IFstream.H"
#include "OFstream.H"

// Parallel communication headers for point fields
#include "processorPointPatch.H"
#include "OPstream.H"
#include "IPstream.H"

#include "syncTools.H"

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

    // Create material point cloud (reads grainProperties internally)
    #include "createMPMCloud.H"

    particles.storeGlobalPositions();
    
    // #include "diagnoseShapeFunctions.H"
    
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

            if (!p.active()) continue;

            
            p.transferToGrid
            (
                mesh,
                nodalMass,
                nodalMomentum,
                nodalStress,
                nodalForce
            );
        }
        

        #include "syncNodalFieldsAcrossProcessors.H"
        
        // Apply boundary conditions to nodal fields
        nodalMass.correctBoundaryConditions();
        nodalMomentum.correctBoundaryConditions();
        nodalForce.correctBoundaryConditions();
        










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
        // Update momentum: p^(n+1) = p^n + dt * F
        nodalMomentum[pointI] += dt * nodalForce[pointI];
        
        // Compute velocity: v^(n+1) = p^(n+1) / m
        nodalVelocity[pointI] = nodalMomentum[pointI] / nodalMass[pointI];
        
        // Acceleration: a = F/m
        nodalAcceleration[pointI] = nodalForce[pointI] / nodalMass[pointI];
    }
    else
    {
        nodalMomentum[pointI] = Zero;
        nodalVelocity[pointI] = Zero;
        nodalAcceleration[pointI] = Zero;
    }
}

// Apply boundary conditions
#include "applyMPMBoundaryConditions.H"


#include "syncNodalVelocitiesAcrossProcessors.H"



Info<< "Step 4c: Updating processor boundaries for velocity fields" << endl;

nodalVelocity.correctBoundaryConditions();
nodalAcceleration.correctBoundaryConditions();

// In parallel, ensure processor boundaries are properly communicated
if (Pstream::parRun())
{
    // This ensures processor boundary patches have correct values
    nodalVelocity.boundaryFieldRef().updateCoeffs();
    nodalAcceleration.boundaryFieldRef().updateCoeffs();
}

        
        // Step 5: Update particle kinematics
        Info<< "Step 5: Updating particle kinematics" << endl;
        forAllIter(materialPointCloud, particles, iter)
        {
            MaterialPointParcel& p = iter();

            if (!p.active()) continue;
            
            // Update deformation gradient
            p.updateDeformationGradient(mesh, nodalVelocity, dt);
            
            // Grid to Particle velocity transfer
            p.transferFromGrid(mesh, nodalVelocity, nodalAcceleration);
        }
        



        // // Step 6: Advect particles
        // Info<< "Step 6: Advecting particles" << endl;
        // particle::trackingData td(particles);
        
        // forAllIter(materialPointCloud, particles, iter)
        // {
        //     MaterialPointParcel& p = iter();
        //     p.move(particles, td, dt);
        // }
        



        Info<< "Step 6: Advecting particles" << endl;

particle::trackingData td(particles);

forAllIter(materialPointCloud, particles, iter)
{
    MaterialPointParcel& p = iter();
    if (!p.active()) continue;
    
    p.move(particles, td, dt);
}

// Parallel particle redistribution
// if (Pstream::parRun())
// {
//     // This transfers particles that crossed processor boundaries
//     particles.distribute();
// }

//         // Step 6: Advect particles
// 

// particle::trackingData td(particles);

// forAllIter(materialPointCloud, particles, iter)
// {
//     MaterialPointParcel& p = iter();
//     if (!p.active()) continue;
    
//     p.move(particles, td, dt);
// }

// // Handle parallel particle redistribution
// if (Pstream::parRun())
// {
//     // Redistribute particles that crossed processor boundaries
//     // This is critical for parallel runs
//     particles.autoMap(mesh);
// }













// Need this for isA<processorPolyPatch> check
#include "processorPolyPatch.H"

Info<< "Step 6b: Detecting outflow particles" << endl;

// -------------------------------------------------------------------------
// INITIALISE OUTFLOW DETECTION (ONCE, CONSISTENTLY ACROSS ALL PROCESSORS)
// -------------------------------------------------------------------------

static bool outflowInitialized = false;
static List<word> outflowPatchNames;
static List<word> outflowDirections;
static boundBox domainBounds;

if (!outflowInitialized)
{
    if (Pstream::master())
        Info<< "  Initializing outflow boundary detection..." << endl;

    // ----- Get global domain bounds -----
    domainBounds = mesh.bounds();

    if (Pstream::parRun())
    {
        vector globalMin = domainBounds.min();
        vector globalMax = domainBounds.max();
        reduce(globalMin, minOp<vector>());
        reduce(globalMax, maxOp<vector>());
        domainBounds = boundBox(globalMin, globalMax);
    }

    // ---------------------------------------------------------------------
    // PARALLEL-SAFE GLOBAL PATCH DISCOVERY AND GLOBAL PATCH-CENTER REDUCTION
    // ---------------------------------------------------------------------
    const polyBoundaryMesh& bMesh = mesh.boundaryMesh();
    const label nPatches = bMesh.size();

    // Local contributions for each physical patch
    List<vector> localCenterSum(nPatches, Zero);
    List<label> localFaceCount(nPatches, 0);

    // ---- Collect local sums ----
    for (label patchI = 0; patchI < nPatches; ++patchI)
    {
        const polyPatch& patch = bMesh[patchI];

        // Skip processor patches but still leave zero contributions
        if (isA<processorPolyPatch>(patch))
            continue;

        vector sum = Zero;
        label count = 0;

        forAll(patch, fI)
        {
            sum += patch.faceCentres()[fI];
            ++count;
        }

        localCenterSum[patchI] = sum;
        localFaceCount[patchI] = count;
    }

    // ---- Global reduction (ALL processors call ALL reduces → no deadlock) ----
//   const label nPatches = bMesh.size();
const label elemsPerPatch = 4; // px,py,pz, count
List<scalar> localBuf(nPatches * elemsPerPatch, scalar(0));

// fill local buffer
for (label patchI = 0; patchI < nPatches; ++patchI)
{
    const vector& s = localCenterSum[patchI];
    const label c = localFaceCount[patchI];

    const label base = patchI * elemsPerPatch;
    localBuf[base + 0] = s.x();
    localBuf[base + 1] = s.y();
    localBuf[base + 2] = s.z();
    localBuf[base + 3] = scalar(c); // store count as scalar to keep homogenous buffer
}

// perform a single global reduction of the whole buffer
List<scalar> globalBuf(localBuf); // will become the reduced result
if (Pstream::parRun())
{
    reduce(globalBuf, sumOp<List<scalar> >());
}

// unpack global results
List<vector> globalCenterSum(nPatches, Zero);
List<label> globalCount(nPatches, 0);

for (label patchI = 0; patchI < nPatches; ++patchI)
{
    const label base = patchI * elemsPerPatch;
    vector s(globalBuf[base + 0], globalBuf[base + 1], globalBuf[base + 2]);
    scalar cntScalar = globalBuf[base + 3];

    globalCenterSum[patchI] = s;
    // convert scalar back to integer safely (rounding)
    globalCount[patchI] = label(round(cntScalar));
}

    // ---------------------------------------------------------------------
    // BUILD GLOBAL LIST OF OUTFLOW PATCHES
    // ---------------------------------------------------------------------
    DynamicList<word> patchNames;
    DynamicList<word> patchDirs;

    for (label patchI = 0; patchI < nPatches; ++patchI)
    {
        const polyPatch& patch = bMesh[patchI];

        // Skip processor patches
        if (isA<processorPolyPatch>(patch))
            continue;

        const word& pName = patch.name();
        const word& pType = patch.type();

        // Identify outflow-like patches
        bool isOut =
            pType == "patch" ||
            pType == "outlet" ||
            pName.find("outlet") != string::npos ||
            pName.find("Outlet") != string::npos ||
            pName.find("outflow") != string::npos;

        if (!isOut)
            continue;

        // No faces globally → ignore
        if (globalCount[patchI] == 0)
            continue;

        // Compute global patch centre
        vector c = globalCenterSum[patchI] / scalar(globalCount[patchI]);

        // Determine direction
        scalar tol = 0.01 * mag(domainBounds.span());
        word direction = "unknown";

        if (mag(c.x() - domainBounds.min().x()) < tol) direction = "xMin";
        else if (mag(c.x() - domainBounds.max().x()) < tol) direction = "xMax";
        else if (mag(c.y() - domainBounds.min().y()) < tol) direction = "yMin";
        else if (mag(c.y() - domainBounds.max().y()) < tol) direction = "yMax";
        else if (mag(c.z() - domainBounds.min().z()) < tol) direction = "zMin";
        else if (mag(c.z() - domainBounds.max().z()) < tol) direction = "zMax";

        patchNames.append(pName);
        patchDirs.append(direction);
    }

    outflowPatchNames = patchNames;
    outflowDirections = patchDirs;

    if (Pstream::master())
    {
        Info<< "  Found " << outflowPatchNames.size() << " outflow boundaries:\n";
        forAll(outflowPatchNames, i)
            Info<< "    " << outflowPatchNames[i] << " (" << outflowDirections[i] << ")\n";
    }

    outflowInitialized = true;
}

// -------------------------------------------------------------------------
// DEACTIVATE PARTICLES OUTSIDE GLOBAL DOMAIN BOUNDS
// -------------------------------------------------------------------------

label nDeactivated = 0;

forAllIter(materialPointCloud, particles, iter)
{
    MaterialPointParcel& p = iter();

    if (!p.active()) continue;

    const vector pos = p.position();
    bool deactivate = false;

    forAll(outflowDirections, i)
    {
        const word& dir = outflowDirections[i];
        scalar eps = SMALL;

        if (dir == "xMin" && pos.x() < domainBounds.min().x() + eps) deactivate = true;
        else if (dir == "xMax" && pos.x() > domainBounds.max().x() - eps) deactivate = true;
        else if (dir == "yMin" && pos.y() < domainBounds.min().y() + eps) deactivate = true;
        else if (dir == "yMax" && pos.y() > domainBounds.max().y() - eps) deactivate = true;
        else if (dir == "zMin" && pos.z() < domainBounds.min().z() + eps) deactivate = true;
        else if (dir == "zMax" && pos.z() > domainBounds.max().z() - eps) deactivate = true;

        if (deactivate)
            break;
    }

    if (deactivate)
    {
        p.deactivate();
        ++nDeactivated;
    }
}

// Report
label nGlobal = returnReduce(nDeactivated, sumOp<label>());

if (Pstream::master() && nGlobal > 0)
    Info<< "  Deactivated " << nGlobal << " particles at outlets" << endl;


















        
        // // Step 6c: Add particles at inflow
        // #include "handleMPMInflow.H"
        
        // Step 6d: Delete marked particles
        #include "handleMPMOutflow.H"

























        // // Step 7: Constitutive update (stress computation)
        // Info<< "Step 7: Updating particle stresses" << endl;
        // forAllIter(materialPointCloud, particles, iter)
        // {
        //     MaterialPointParcel& p = iter();
            
        //     // Get material properties for this particle
        //     const dictionary& matProps = materialProperties[p.materialID()];
            
        //     // Perform constitutive update (includes crystal plasticity)
        //     p.constitutiveUpdate(dt, matProps);
        //     // p.computeElasticStress(matProps);//just elastic for debugging
        // }
        



        // Step 7: Constitutive update (stress computation)
// Step 7: Constitutive update (stress computation)
Info<< "Step 7: Updating particle stresses" << endl;

if (Pstream::parRun())
{
    // Pstream::scatter(label(0));  // Barrier before starting
        // Barrier - wait for all processors
    label dummy = 0;
    reduce(dummy, sumOp<label>());
}

label nProcessed = 0;

forAllIter(materialPointCloud, particles, iter)
{
    MaterialPointParcel& p = iter();
    if (!p.active()) continue;
    
    const Material& material = materials[p.materialID()];
    
    // Compute velocity gradient from shape functions
    scalarField N;
    List<vector> gradN;
    labelList nodeIDs;
    
    p.shapeFunctions(mesh, N, gradN, nodeIDs);
    

    if (N.size() == 0 || nodeIDs.size() == 0)
    {
        // Skip this particle - it's in an invalid state
        continue;
    }
    
    // Verify all node indices are in bounds
    bool allValid = true;
    forAll(nodeIDs, i)
    {
        if (nodeIDs[i] >= nodalVelocity.size())
        {
            allValid = false;
            break;
        }
    }
    
    if (!allValid)
    {
        WarningInFunction
            << "Particle has out-of-bounds node indices - skipping"
            << endl;
        continue;
    }
    
    tensor L = tensor::zero;
    forAll(nodeIDs, i)
    {
        label nodeI = nodeIDs[i];
        L += nodalVelocity[nodeI] * gradN[i];
    }
    
    p.constitutiveUpdate(material, L, dt);
    nProcessed++;
}

// Verify all processors processed some particles
if (Pstream::parRun())
{
    label nProcessedGlobal = returnReduce(nProcessed, sumOp<label>());
    
    if (Pstream::master())
    {
        Info<< "  Processed " << nProcessedGlobal 
            << " particles globally in constitutive update" << endl;
    }
}
























        // Step 8: Update history variables and check failure
        Info<< "Step 8: Checking material state" << endl;
        label nFailed = 0;
        scalar maxEqPlasticStrain = 0.0;
        
        forAllIter(materialPointCloud, particles, iter)
        {
            MaterialPointParcel& p = iter();

            if (!p.active()) continue;
            
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