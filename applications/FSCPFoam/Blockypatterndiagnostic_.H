// ===================================================================
// FOCUSED BLOCKY PATTERN DIAGNOSTIC
// Add this RIGHT AFTER Advect.H in your main loop
// ===================================================================
{
    Info << "\n=== CHECKING FOR BLOCKY PATTERNS ===" << endl;
    
    // Check F field for identical adjacent cell values
    label nBlockyF = 0;
    label nCheckedF = 0;
    
    forAll(mesh.C(), cellI)
    {
        const labelList& neighbors = mesh.cellCells()[cellI];
        
        forAll(neighbors, neiI)
        {
            label neiCell = neighbors[neiI];
            
            // Check if F is EXACTLY identical (within machine precision)
            bool identical = true;
            scalar maxDiff = 0.0;
            for (int i = 0; i < 9; i++)
            {
                scalar diff = mag(F[cellI].component(i) - F[neiCell].component(i));
                maxDiff = max(maxDiff, diff);
                if (diff > 1e-15)
                {
                    identical = false;
                }
            }
            
            if (identical)
            {
                nBlockyF++;
                // Print first few occurrences
                if (nBlockyF <= 5)
                {
                    Info << "  Identical F found: cells " << cellI 
                         << " and " << neiCell 
                         << ", max diff = " << maxDiff << endl;
                }
            }
            
            nCheckedF++;
        }
    }
    
    scalar blockyPercentF = 100.0 * scalar(nBlockyF) / max(scalar(nCheckedF), 1.0);
    Info << "F: " << nBlockyF << " / " << nCheckedF 
         << " (" << blockyPercentF << "%) neighboring cell pairs have identical values" << endl;
    
    // Check sigma
    for (label h = 0; h < N_Ori; h++)
    {
        label nBlockySigma = 0;
        label nCheckedSigma = 0;
        
        forAll(mesh.C(), cellI)
        {
            const labelList& neighbors = mesh.cellCells()[cellI];
            
            forAll(neighbors, neiI)
            {
                label neiCell = neighbors[neiI];
                
                bool identical = true;
                for (int i = 0; i < 9; i++)
                {
                    if (mag(sigma[h][cellI].component(i) - sigma[h][neiCell].component(i)) > 1e-15)
                    {
                        identical = false;
                        break;
                    }
                }
                
                if (identical) nBlockySigma++;
                nCheckedSigma++;
            }
        }
        
        scalar blockyPercentSigma = 100.0 * scalar(nBlockySigma) / max(scalar(nCheckedSigma), 1.0);
        Info << "sigma[" << h << "]: " << blockyPercentSigma << "% blocky" << endl;
    }
    
    // Check rho
    label nBlockyRho = 0;
    label nCheckedRho = 0;
    
    forAll(mesh.C(), cellI)
    {
        const labelList& neighbors = mesh.cellCells()[cellI];
        
        forAll(neighbors, neiI)
        {
            label neiCell = neighbors[neiI];
            
            if (mag(rho[cellI] - rho[neiCell]) < 1e-15)
            {
                nBlockyRho++;
            }
            nCheckedRho++;
        }
    }
    
    scalar blockyPercentRho = 100.0 * scalar(nBlockyRho) / max(scalar(nCheckedRho), 1.0);
    Info << "rho: " << blockyPercentRho << "% blocky" << endl;
    
    // Check boundary vs interior consistency for F
    Info << "\n=== BOUNDARY vs INTERIOR for F ===" << endl;
    forAll(mesh.boundary(), patchI)
    {
        const fvPatch& patch = mesh.boundary()[patchI];
        
        if (!patch.coupled() && patch.size() > 0)
        {
            scalar avgDiff = 0.0;
            scalar maxDiff = 0.0;
            label nFaces = 0;
            
            forAll(patch, faceI)
            {
                label cellI = patch.faceCells()[faceI];
                
                for (int i = 0; i < 9; i++)
                {
                    scalar diff = mag(F[cellI].component(i) - 
                                     F.boundaryField()[patchI][faceI].component(i));
                    avgDiff += diff;
                    maxDiff = max(maxDiff, diff);
                }
                nFaces++;
            }
            
            avgDiff /= (nFaces * 9.0);
            
            Info << "Patch " << patch.name() 
                 << ": avg diff = " << avgDiff 
                 << ", max diff = " << maxDiff << endl;
            
            if (maxDiff > 0.01)
            {
                Warning << "  Large mismatch on patch " << patch.name() << "!" << endl;
            }
        }
    }
    
    Info << "========================================\n" << endl;
}
