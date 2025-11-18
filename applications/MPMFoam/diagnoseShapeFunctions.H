/*---------------------------------------------------------------------------*\
  DIAGNOSTIC: Check Shape Function Node Mapping
  
  This diagnostic should be run ONCE at the start of your simulation to verify
  that the shape functions are correctly mapping to nodes.
  
  Add this to your MPMFoam.C right after creating the particle cloud.
\*---------------------------------------------------------------------------*/

Info<< "\n========================================" << endl;
Info<< "DIAGNOSTIC: Checking Shape Functions" << endl;
Info<< "========================================\n" << endl;

// Pick a test particle (first one)
if (particles.size() > 0)
{
    materialPointCloud::iterator iter = particles.begin();
    const MaterialPointParcel& testParticle = iter();
    
    // Get shape functions
    scalarField N;
    List<vector> gradN;
    labelList nodeIDs;
    
    testParticle.shapeFunctions(mesh, N, gradN, nodeIDs);
    
    Info<< "Test particle:" << endl;
    Info<< "  Position: " << testParticle.position() << endl;
    Info<< "  Cell: " << testParticle.cell() << endl;
    Info<< "  Number of nodes: " << nodeIDs.size() << endl;
    
    // Check partition of unity
    scalar sumN = 0.0;
    vector sumGrad = Zero;
    
    forAll(N, i)
    {
        sumN += N[i];
        sumGrad += gradN[i];
    }
    
    Info<< "\nPartition of unity checks:" << endl;
    Info<< "  Sum of N_i: " << sumN << " (should be 1.0)" << endl;
    Info<< "  Sum of grad(N_i): " << sumGrad << " (should be ~0)" << endl;
    
    if (mag(sumN - 1.0) > 1e-4)
    {
        FatalErrorInFunction
            << "Shape functions don't sum to 1!" << nl
            << "Sum = " << sumN << nl
            << "This indicates a serious error in shape function implementation!"
            << exit(FatalError);
    }
    
    if (mag(sumGrad) > 1e-4)
    {
        WarningInFunction
            << "Shape function gradients don't sum to zero!" << nl
            << "Sum = " << sumGrad << nl
            << "This may indicate incorrect gradient calculation."
            << endl;
    }
    
    // Check node positions relative to particle
    Info<< "\nNode information:" << endl;
    const pointField& points = mesh.points();
    
    forAll(nodeIDs, i)
    {
        label nodeI = nodeIDs[i];
        vector nodePos = points[nodeI];
        vector relPos = nodePos - testParticle.position();
        
        Info<< "  Node " << i << " (global " << nodeI << "):" << endl;
        Info<< "    Position: " << nodePos << endl;
        Info<< "    Relative: " << relPos << endl;
        Info<< "    N[" << i << "] = " << N[i] << endl;
        Info<< "    gradN[" << i << "] = " << gradN[i] << endl;
    }
    
    // Check if cell is hexahedral
    if (nodeIDs.size() == 8)
    {
        Info<< "\nHexahedral cell detected." << endl;
        
        // Get cell bounding box
        vector minPt = points[nodeIDs[0]];
        vector maxPt = points[nodeIDs[0]];
        
        forAll(nodeIDs, i)
        {
            minPt = min(minPt, points[nodeIDs[i]]);
            maxPt = max(maxPt, points[nodeIDs[i]]);
        }
        
        vector cellSize = maxPt - minPt;
        vector cellCenter = 0.5 * (minPt + maxPt);
        
        Info<< "  Cell bounds: " << minPt << " to " << maxPt << endl;
        Info<< "  Cell size: " << cellSize << endl;
        Info<< "  Cell center: " << cellCenter << endl;
        
        // Check if particle is actually inside cell
        vector pPos = testParticle.position();
        bool insideCell = 
            (pPos.x() >= minPt.x() - 1e-6) && (pPos.x() <= maxPt.x() + 1e-6) &&
            (pPos.y() >= minPt.y() - 1e-6) && (pPos.y() <= maxPt.y() + 1e-6) &&
            (pPos.z() >= minPt.z() - 1e-6) && (pPos.z() <= maxPt.z() + 1e-6);
        
        if (!insideCell)
        {
            WarningInFunction
                << "Test particle is OUTSIDE its cell bounding box!" << endl;
        }
        
        // Verify node ordering by checking corners
        Info<< "\nVerifying hex node ordering:" << endl;
        
        for (label i = 0; i < 8; i++)
        {
            const vector& pos = points[nodeIDs[i]];
            
            // Expected corner in standard ordering
            label expectedI = 0;
            if (mag(pos.x() - maxPt.x()) < 1e-6 * mag(cellSize)) expectedI += 1;
            if (mag(pos.y() - maxPt.y()) < 1e-6 * mag(cellSize)) expectedI += 2;
            if (mag(pos.z() - maxPt.z()) < 1e-6 * mag(cellSize)) expectedI += 4;
            
            Info<< "  Physical node " << i << " -> Expected local index " << expectedI;
            
            if (i != expectedI)
            {
                Info<< " [MISMATCH!]";
            }
            Info<< endl;
        }
    }
    else
    {
        Info<< "\nNon-hexahedral cell (n=" << nodeIDs.size() << " nodes)" << endl;
    }
}
else
{
    WarningInFunction
        << "No particles in cloud - cannot run diagnostics!"
        << endl;
}

Info<< "\n========================================" << endl;
Info<< "End of shape function diagnostics" << endl;
Info<< "========================================\n" << endl;

// ************************************************************************* //
