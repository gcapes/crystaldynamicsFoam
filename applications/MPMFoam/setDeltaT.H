/*---------------------------------------------------------------------------*\
  setDeltaT.H - Adaptive time stepping for MPM (Parallel Safe)
  
  Adjusts time step based on Courant number and wave speed
\*---------------------------------------------------------------------------*/

if (adjustTimeStep)
{
    // Calculate new time step based on Courant number
    scalar newDeltaT = GREAT;
    
    if (CoNum > SMALL)
    {
        newDeltaT = maxCo * minCellSize / (maxVelocityMag + SMALL);
    }
    
    // Estimate wave speed from material properties
    if (materials.size() > 0 && materials.set(0))
    {
        const Material& mat = materials[0];
        
        // Get elastic modulus (varies by material type)
        scalar E = 0.0;
        if (mat.elasticType == "isotropic")
        {
            E = mat.lambda + 2.0 * mat.mu;
        }
        else if (mat.elasticType == "cubic" || mat.elasticType == "anisotropic")
        {
            E = mat.C11;
        }
        else if (mat.elasticType == "hexagonal")
        {
            E = 0.5 * (mat.C11 + mat.C33);
        }
        else if (mat.elasticType == "orthorhombic")
        {
            E = (mat.C11 + mat.C22 + mat.C33) / 3.0;
        }
        
        scalar rho = mat.rho_0;
        
        if (E > SMALL && rho > SMALL)
        {
            scalar waveSpeed = Foam::sqrt(E/rho);
            
            // Time step based on wave speed (CFL condition)
            scalar dtWave = maxCo * minCellSize / waveSpeed;
            
            newDeltaT = min(newDeltaT, dtWave);
        }
    }
    
    // Limit time step change rate
    scalar deltaTFact = min(max(newDeltaT/runTime.deltaT().value(), 0.5), 1.2);
    
    newDeltaT = min
    (
        deltaTFact * runTime.deltaT().value(),
        maxDeltaT
    );
    
    runTime.setDeltaT(newDeltaT);
    
    Info<< "deltaT = " <<  runTime.deltaT().value() << endl;
}

// ************************************************************************* //

















// // Adjust time step based on Courant number

// if (adjustTimeStep && minCellSize > SMALL && maxVelocity > SMALL)
// {
//     // Consider both velocity and wave speed
//     scalar waveSpeed = 0.0;
    
//     forAll(materials, matI)
//     {
//         if (materials.set(matI))
//         {
//             const Material& mat = materials[matI];
            
//             // Get elastic modulus - depends on elastic type
//             scalar E = 0.0;
//             if (mat.elasticType == "isotropic")
//             {
//                 // E = lambda + 2*mu for isotropic materials (P-wave modulus)
//                 E = mat.lambda + 2.0 * mat.mu;
//             }
//             else if (mat.elasticType == "cubic" || mat.elasticType == "anisotropic")
//             {
//                 // Use C11 as representative modulus
//                 E = mat.C11;
//             }
//             else if (mat.elasticType == "hexagonal")
//             {
//                 // Average of C11 and C33
//                 E = 0.5 * (mat.C11 + mat.C33);
//             }
//             else if (mat.elasticType == "orthorhombic")
//             {
//                 // Average of all normal stiffnesses
//                 E = (mat.C11 + mat.C22 + mat.C33) / 3.0;
//             }
            
//             scalar rho = mat.rho_0;
//             if (E > SMALL && rho > SMALL)
//             {
//                 waveSpeed = max(waveSpeed, Foam::sqrt(E/rho));
//             }
//         }
//     }
    
//     // Use the more restrictive of velocity-based and wave-speed-based time step
//     scalar maxSpeed = max(maxVelocity, waveSpeed);
    
//     scalar deltaT = min
//     (
//         maxDeltaT,
//         maxCo * minCellSize / (maxSpeed + SMALL)
//     );
    
//     // Limit the rate of change of time step
//     scalar deltaTFact = min(max(deltaT/runTime.deltaT().value(), 0.5), 2.0);
    
//     runTime.setDeltaT(deltaTFact * runTime.deltaT().value());
    
//     Info<< "Adjusted time step to: " << runTime.deltaT().value() << endl;
// }

// // ************************************************************************* //