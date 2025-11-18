/*---------------------------------------------------------------------------*\
  CORRECTED BOUNDARY CONDITIONS FOR MPM
  
  Key fixes:
  1. Remove node grouping/averaging - it interferes with proper mass mapping
  2. Apply velocity directly without mass manipulation
  3. Ensure consistency between velocity and momentum
  
  Replace applyMPMBoundaryConditions.H with this version.
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
            
            if (bcType == "fixedValue")
            {
                vector value(patchDict.lookup("value"));
                
                // Apply to nodes on this patch
                const labelList& patchPoints = patch.meshPoints();
                
                Info<< "  Applying fixedValue BC on " << patchName 
                    << ": v = " << value << " to " << patchPoints.size() << " nodes" << endl;
                
                // CRITICAL FIX: Apply BC to EACH node individually
                // Don't average or group nodes - this causes mass variations
                forAll(patchPoints, i)
                {
                    label pointI = patchPoints[i];
                    
                    // Only apply BC to nodes with mass (i.e., nodes that have particles)
                    if (nodalMass[pointI] > SMALL)
                    {
                        // Set velocity to prescribed value
                        nodalVelocity[pointI] = value;
                        
                        // Update momentum to be consistent: p = m*v
                        // KEEP the mass as-is from P2G transfer!
                        nodalMomentum[pointI] = nodalMass[pointI] * value;
                        
                        // Zero out acceleration (no additional forces on prescribed nodes)
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
                    
                    if (nodalMass[pointI] > SMALL)
                    {
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
                            nodalMomentum[pointI] = nodalMass[pointI] * nodalVelocity[pointI];
                            
                            // Remove normal component of acceleration
                            vector an = (nodalAcceleration[pointI] & avgNormal) * avgNormal;
                            nodalAcceleration[pointI] -= an;
                        }
                    }
                }
            }
            else if (bcType == "traction")
            {
                vector traction(patchDict.lookup("traction"));
                
                Info<< "  Applying traction BC on " << patchName 
                    << ": t = " << traction << endl;
                
                WarningInFunction
                    << "Traction BCs should be applied in force calculation (Step 3), not here!"
                    << endl;
            }
            else if (bcType == "empty")
            {
                // Do nothing for empty patches (2D simulations)
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

// =========================================================================
// OPTIONAL: 2D PLANE STRESS/STRAIN ENFORCEMENT
// =========================================================================
// For 2D simulations, you may want to ensure front/back faces stay synchronized
// But DO NOT average mass - this causes the velocity degradation issue!
// =========================================================================

// Check if this is a 2D case (thin in one direction)
vector meshExtent = mesh.bounds().max() - mesh.bounds().min();
scalar minExtent = min(meshExtent.x(), min(meshExtent.y(), meshExtent.z()));
scalar maxExtent = max(meshExtent.x(), max(meshExtent.y(), meshExtent.z()));

bool is2D = (minExtent / maxExtent < 0.1);  // One direction is <10% of others

if (is2D)
{
    Info<< "  2D case detected - synchronizing thin-direction node pairs" << endl;
    
    // Identify thin direction
    label thinDir = 0;
    if (meshExtent.x() < meshExtent.y() && meshExtent.x() < meshExtent.z()) thinDir = 0;
    else if (meshExtent.y() < meshExtent.x() && meshExtent.y() < meshExtent.z()) thinDir = 1;
    else thinDir = 2;
    
    // Find front and back patches
    label frontPatch = -1, backPatch = -1;
    
    if (thinDir == 0)
    {
        frontPatch = boundaryMesh.findPatchID("Front");
        backPatch = boundaryMesh.findPatchID("Back");
    }
    else if (thinDir == 1)
    {
        frontPatch = boundaryMesh.findPatchID("Top");
        backPatch = boundaryMesh.findPatchID("Bottom");
    }
    else
    {
        frontPatch = boundaryMesh.findPatchID("Front");
        backPatch = boundaryMesh.findPatchID("Back");
    }
    
    if (frontPatch != -1 && backPatch != -1)
    {
        const labelList& frontPoints = boundaryMesh[frontPatch].meshPoints();
        const labelList& backPoints = boundaryMesh[backPatch].meshPoints();
        
        // For each front point, find corresponding back point and synchronize VELOCITY only
        forAll(frontPoints, i)
        {
            label frontPt = frontPoints[i];
            const vector& frontPos = mesh.points()[frontPt];
            
            // Find back point with matching position in other 2 directions
            scalar bestDist = GREAT;
            label matchingBackPt = -1;
            
            forAll(backPoints, j)
            {
                label backPt = backPoints[j];
                const vector& backPos = mesh.points()[backPt];
                
                // Compute distance in the 2 non-thin directions
                vector diff = frontPos - backPos;
                if (thinDir == 0) diff.x() = 0;
                else if (thinDir == 1) diff.y() = 0;
                else diff.z() = 0;
                
                scalar dist = mag(diff);
                
                if (dist < bestDist && dist < 1e-6 * maxExtent)
                {
                    bestDist = dist;
                    matchingBackPt = backPt;
                }
            }
            
            if (matchingBackPt != -1)
            {
                // Synchronize velocities (average)
                // But DON'T touch mass - this is the key fix!
                vector avgVel = 0.5 * (nodalVelocity[frontPt] + nodalVelocity[matchingBackPt]);
                
                nodalVelocity[frontPt] = avgVel;
                nodalVelocity[matchingBackPt] = avgVel;
                
                // Update momentum using ORIGINAL mass from P2G
                nodalMomentum[frontPt] = nodalMass[frontPt] * avgVel;
                nodalMomentum[matchingBackPt] = nodalMass[matchingBackPt] * avgVel;
            }
        }
    }
}

// ************************************************************************* //