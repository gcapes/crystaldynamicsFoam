/*---------------------------------------------------------------------------*\
  syncNodalVelocitiesAcrossProcessors.H - WORKING VERSION
  
  Save this file in your solver directory.
  Synchronizes velocities after momentum update (Step 4).
  Ensures consistent velocities at shared nodes.
\*---------------------------------------------------------------------------*/

#include "syncTools.H"

if (Pstream::parRun())
{
    // After momentum update, velocities at shared nodes need to be consistent
    // across processors
    
    // Extract velocity and acceleration as lists
    List<vector> velocityInternal(nodalVelocity.primitiveField());
    List<vector> accelerationInternal(nodalAcceleration.primitiveField());
    
    // For velocities, we want to average at shared nodes
    // since v = p/m has already been computed on each processor
    
    // First, sum the velocities
    syncTools::syncPointList
    (
        mesh,
        velocityInternal,
        plusEqOp<vector>(),
        vector::zero
    );
    
    syncTools::syncPointList
    (
        mesh,
        accelerationInternal,
        plusEqOp<vector>(),
        vector::zero
    );
    
    // Count how many processors share each point
    scalarList nProcsPerPoint(mesh.nPoints(), 1.0);
    
    syncTools::syncPointList
    (
        mesh,
        nProcsPerPoint,
        plusEqOp<scalar>(),
        0.0
    );
    
    // Average the velocities at shared points
    forAll(velocityInternal, pointI)
    {
        if (nProcsPerPoint[pointI] > 1.0)
        {
            velocityInternal[pointI] /= nProcsPerPoint[pointI];
            accelerationInternal[pointI] /= nProcsPerPoint[pointI];
        }
    }
    
    // Update point fields with averaged values
    forAll(nodalVelocity, pointI)
    {
        nodalVelocity[pointI] = velocityInternal[pointI];
        nodalAcceleration[pointI] = accelerationInternal[pointI];
    }
}

// ************************************************************************* //