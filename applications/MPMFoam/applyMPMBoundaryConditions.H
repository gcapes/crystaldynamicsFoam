/*---------------------------------------------------------------------------*\
  PROPER BOUNDARY CONDITIONS FOR MPM
  
  Replace applyMPMBoundaryConditions.H with this version
  
  Key insight: In explicit MPM, velocity BCs should be applied AFTER
  momentum update to properly couple the boundary to the domain.
\*---------------------------------------------------------------------------*/

// Apply boundary conditions to nodal fields

// Read boundary condition dictionary
IOdictionary mpmBoundaryDict
(
    IOobject
    (
        "mpmBoundaryConditions",
        runTime.constant(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    )
);

// Get boundary patches
const polyBoundaryMesh& boundaryMesh = mesh.boundaryMesh();

// Apply velocity boundary conditions if specified
if (mpmBoundaryDict.found("velocityBCs"))
{
    const dictionary& velocityBCs = mpmBoundaryDict.subDict("velocityBCs");
    
    forAllConstIter(dictionary, velocityBCs, iter)
    {
        const word& patchName = iter().keyword();
        const dictionary& patchDict = iter().dict();
        
        // Find patch
        label patchI = boundaryMesh.findPatchID(patchName);
        
        if (patchI != -1)
        {
            const polyPatch& patch = boundaryMesh[patchI];
            
            // Get boundary condition type
            word bcType(patchDict.lookup("type"));
            
            // if (bcType == "fixedValue")
            // {
            //     vector value(patchDict.lookup("value"));
                
            //     // Apply to nodes on this patch
            //     const labelList& patchPoints = patch.meshPoints();
                
            //     Info<< "  Applying fixedValue BC on " << patchName 
            //         << ": v = " << value << " to " << patchPoints.size() << " nodes" << endl;
                
            //     forAll(patchPoints, i)
            //     {
            //         label pointI = patchPoints[i];
                    
            //         // CRITICAL: Update both velocity AND momentum consistently
            //         nodalVelocity[pointI] = value;
                    
            //         if (nodalMass[pointI] > SMALL)
            //         {
            //             // Update momentum to be consistent: p = m*v
            //             nodalMomentum[pointI] = nodalMass[pointI] * value;
            //         }
                    
            //         // Zero out acceleration (no external forces on prescribed nodes)
            //         nodalAcceleration[pointI] = Zero;
            //     }
            // }
            if (bcType == "fixedValue")
{
    vector value(patchDict.lookup("value"));
    
    const labelList& patchPoints = patch.meshPoints();
    
    // Group nodes by Y-Z position (X is thin direction in 2D)
    Map<labelList> nodesAtPosition;
    
    forAll(patchPoints, i)
    {
        label pointI = patchPoints[i];
        const vector& pos = mesh.points()[pointI];
        
        // Create key from Y-Z position (ignore X)
        label key = label(1e6 * pos.y()) + label(1e3 * pos.z());
        
        if (!nodesAtPosition.found(key))
        {
            nodesAtPosition.insert(key, labelList());
        }
        
        nodesAtPosition[key].append(pointI);
    }
    
    // Apply BC to all nodes at each position
    forAllConstIter(Map<labelList>, nodesAtPosition, iter)
    {
        const labelList& nodes = iter();
        
        // Compute average properties for consistency
        scalar avgMass = 0;
        forAll(nodes, i)
        {
            avgMass += nodalMass[nodes[i]];
        }
        avgMass /= nodes.size();
        
        // Apply BC to all nodes
        forAll(nodes, i)
        {
            label pointI = nodes[i];
            
            nodalVelocity[pointI] = value;
            
            if (avgMass > SMALL)
            {
                nodalMomentum[pointI] = avgMass * value;
                nodalMass[pointI] = avgMass;  // Synchronize mass
            }
            
            nodalAcceleration[pointI] = Zero;
        }
    }
}
            else if (bcType == "symmetry" || bcType == "slip")
            {
                // Get patch normal
                const vectorField& patchNormals = patch.faceNormals();
                
                const labelList& patchPoints = patch.meshPoints();
                
                Info<< "  Applying symmetry BC on " << patchName 
                    << " to " << patchPoints.size() << " nodes" << endl;
                
                forAll(patchPoints, i)
                {
                    label pointI = patchPoints[i];
                    
                    // Average normal for this point
                    vector avgNormal = Zero;
                    label nFaces = 0;
                    
                    const labelList& pointFaces = patch.pointFaces()[i];
                    forAll(pointFaces, fi)
                    {
                        avgNormal += patchNormals[pointFaces[fi]];
                        nFaces++;
                    }
                    
                    if (nFaces > 0)
                    {
                        avgNormal /= mag(avgNormal) + SMALL;
                        
                        // Remove normal component of velocity
                        vector vn = (nodalVelocity[pointI] & avgNormal) * avgNormal;
                        nodalVelocity[pointI] -= vn;
                        
                        // Update momentum consistently
                        if (nodalMass[pointI] > SMALL)
                        {
                            nodalMomentum[pointI] = nodalMass[pointI] * nodalVelocity[pointI];
                        }
                        
                        // Remove normal component of acceleration
                        vector an = (nodalAcceleration[pointI] & avgNormal) * avgNormal;
                        nodalAcceleration[pointI] -= an;
                    }
                }
            }
            else if (bcType == "traction")
            {
                vector traction(patchDict.lookup("traction"));
                
                Info<< "  Applying traction BC on " << patchName 
                    << ": t = " << traction << endl;
                
                // Traction BCs should be applied in Step 3 (force calculation)
                // not here after momentum update
                WarningInFunction
                    << "Traction BC on " << patchName 
                    << " should be applied before momentum update!" << endl;
            }
            else if (bcType == "empty")
            {
                // Do nothing for empty patches (2D)
            }
            else
            {
                WarningInFunction
                    << "Unknown BC type: " << bcType << " on patch " << patchName << endl;
            }
        }
        else
        {
            WarningInFunction
                << "Patch " << patchName << " not found in mesh" << endl;
        }
    }
}

// ************************************************************************* //