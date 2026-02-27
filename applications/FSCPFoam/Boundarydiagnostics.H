// ===================================================================
// BOUNDARY UPDATE DIAGNOSTICS
// Add this after StressUpdate.H to verify boundaries are correct
// ===================================================================
{
    Info << "\n=== BOUNDARY UPDATE DIAGNOSTICS ===" << endl;
    
    // Check if interior and boundary values are consistent near boundaries
    forAll(mesh.boundary(), patchI)
    {
        const fvPatch& patch = mesh.boundary()[patchI];
        
        if (!patch.coupled() && patch.size() > 0)
        {
            Info << "\nPatch: " << patch.name() << endl;
            
            // Check F consistency
            scalar F_interior_avg = 0.0;
            scalar F_boundary_avg = 0.0;
            label nFaces = 0;
            
            forAll(patch, faceI)
            {
                label cellI = patch.faceCells()[faceI];
                
                // Compare interior cell value to boundary face value
                for (int i = 0; i < 9; i++)
                {
                    F_interior_avg += F[cellI].component(i);
                    F_boundary_avg += F.boundaryField()[patchI][faceI].component(i);
                }
                nFaces++;
            }
            
            F_interior_avg /= (nFaces * 9.0);
            F_boundary_avg /= (nFaces * 9.0);
            
            Info << "  F: interior avg = " << F_interior_avg 
                 << ", boundary avg = " << F_boundary_avg << endl;
            
            if (mag(F_interior_avg - F_boundary_avg) > 0.01)
            {
                Warning << "  Large F mismatch between interior and boundary!" << endl;
            }
            
            // Check rho consistency
            scalar rho_interior_avg = 0.0;
            scalar rho_boundary_avg = 0.0;
            
            forAll(patch, faceI)
            {
                label cellI = patch.faceCells()[faceI];
                rho_interior_avg += rho[cellI];
                rho_boundary_avg += rho.boundaryField()[patchI][faceI];
            }
            
            rho_interior_avg /= nFaces;
            rho_boundary_avg /= nFaces;
            
            Info << "  rho: interior avg = " << rho_interior_avg 
                 << ", boundary avg = " << rho_boundary_avg << endl;
            
            if (mag(rho_interior_avg - rho_boundary_avg) / rho_interior_avg > 0.01)
            {
                Warning << "  Large rho mismatch between interior and boundary!" << endl;
            }
            
            // Check phase fields
            for (label h = 0; h < N_Ori; h++)
            {
                scalar phase_interior_avg = 0.0;
                scalar phase_boundary_avg = 0.0;
                
                forAll(patch, faceI)
                {
                    label cellI = patch.faceCells()[faceI];
                    phase_interior_avg += PopBal[h][cellI];
                    phase_boundary_avg += PopBal[h].boundaryField()[patchI][faceI];
                }
                
                phase_interior_avg /= nFaces;
                phase_boundary_avg /= nFaces;
                
                Info << "  PopBal[" << h << "]: interior avg = " << phase_interior_avg 
                     << ", boundary avg = " << phase_boundary_avg << endl;
                
                if (mag(phase_interior_avg - phase_boundary_avg) > 0.01)
                {
                    Warning << "  Large PopBal[" << h << "] mismatch!" << endl;
                }
            }
            
            // Check sigma
            for (label h = 0; h < N_Ori; h++)
            {
                scalar sigma_interior_mag = 0.0;
                scalar sigma_boundary_mag = 0.0;
                
                forAll(patch, faceI)
                {
                    label cellI = patch.faceCells()[faceI];
                    
                    for (int i = 0; i < 9; i++)
                    {
                        sigma_interior_mag += mag(sigma[h][cellI].component(i));
                        sigma_boundary_mag += mag(sigma[h].boundaryField()[patchI][faceI].component(i));
                    }
                }
                
                sigma_interior_mag /= (nFaces * 9.0);
                sigma_boundary_mag /= (nFaces * 9.0);
                
                Info << "  |sigma[" << h << "]|: interior avg = " << sigma_interior_mag 
                     << ", boundary avg = " << sigma_boundary_mag << endl;
                
                if (mag(sigma_interior_mag - sigma_boundary_mag) / max(sigma_interior_mag, SMALL) > 0.1)
                {
                    Warning << "  Large sigma[" << h << "] mismatch!" << endl;
                }
            }
        }
    }
    
    // Check for "blocky" patterns - cells with identical values
    Info << "\n=== CHECKING FOR BLOCKY PATTERNS ===" << endl;
    
    // Sample some interior cells and check if adjacent cells have identical values
    label nBlocky = 0;
    label nChecked = 0;
    
    forAll(mesh.C(), cellI)
    {
        if (cellI % 100 == 0) // Sample every 100th cell
        {
            const labelList& neighbors = mesh.cellCells()[cellI];
            
            forAll(neighbors, neiI)
            {
                label neiCell = neighbors[neiI];
                
                // Check if F is EXACTLY identical
                bool identical = true;
                for (int i = 0; i < 9; i++)
                {
                    if (mag(F[cellI].component(i) - F[neiCell].component(i)) > SMALL)
                    {
                        identical = false;
                        break;
                    }
                }
                
                if (identical && mag(mesh.C()[cellI] - mesh.C()[neiCell]) > SMALL)
                {
                    nBlocky++;
                }
                
                nChecked++;
            }
        }
    }
    
    scalar blockyPercent = 100.0 * scalar(nBlocky) / max(scalar(nChecked), 1.0);
    
    Info << "Blocky pattern check: " << nBlocky << " / " << nChecked 
         << " (" << blockyPercent << "%) cell pairs have identical F" << endl;
    
    if (blockyPercent > 5.0)
    {
        Warning << "High percentage of blocky patterns detected!" << endl;
        Warning << "This suggests boundary conditions are not being applied correctly." << endl;
    }
    
    Info << "====================================\n" << endl;
}
