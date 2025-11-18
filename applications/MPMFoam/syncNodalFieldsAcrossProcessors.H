/*---------------------------------------------------------------------------*\
  syncNodalFieldsAcrossProcessors.H - WORKING VERSION
  
  Save this file in your solver directory.
  This replaces the broken version with OPstream/IPstream.
  Uses OpenFOAM's syncTools to avoid deadlocks.
\*---------------------------------------------------------------------------*/

#include "syncTools.H"
#include "pointFields.H"

if (Pstream::parRun())
{
    Info<< "  Synchronizing nodal fields across " << Pstream::nProcs() 
        << " processors" << endl;
    
    // Extract internal field values as lists
    scalarList massInternal(nodalMass.primitiveField());
    List<vector> momentumInternal(nodalMomentum.primitiveField());
    List<vector> forceInternal(nodalForce.primitiveField());
    List<tensor> stressInternal(nodalStress.primitiveField());
    
    // CRITICAL: Use syncTools to synchronize across processor boundaries
    // This sums contributions from all processors at shared points
    
    syncTools::syncPointList
    (
        mesh,
        massInternal,
        plusEqOp<scalar>(),  // Sum operation
        0.0                   // Null value
    );
    
    syncTools::syncPointList
    (
        mesh,
        momentumInternal,
        plusEqOp<vector>(),
        vector::zero
    );
    
    syncTools::syncPointList
    (
        mesh,
        forceInternal,
        plusEqOp<vector>(),
        vector::zero
    );
    
    syncTools::syncPointList
    (
        mesh,
        stressInternal,
        plusEqOp<tensor>(),
        tensor::zero
    );
    
    // Update the point fields with synchronized values
    forAll(nodalMass, pointI)
    {
        nodalMass[pointI] = massInternal[pointI];
        nodalMomentum[pointI] = momentumInternal[pointI];
        nodalForce[pointI] = forceInternal[pointI];
        nodalStress[pointI] = stressInternal[pointI];
    }
    
    // Debug output to verify synchronization
    scalar totalMassLocal = gSum(massInternal);
    scalar minMass = gMin(massInternal);
    
    Info<< "    Total nodal mass after sync: " << totalMassLocal << " kg" << endl;
    
    if (minMass < -SMALL)
    {
        WarningInFunction
            << "Negative nodal mass detected: " << minMass 
            << " - check P2G transfer!" << endl;
    }
}

// ************************************************************************* //