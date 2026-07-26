/*---------------------------------------------------------------------------*\
  MPMCourantNo.H - Calculate Courant number for MPM (Parallel Safe)
  
  Co = max(|v| * dt / dx)
  
  For explicit MPM, we want Co < 1 for stability
\*---------------------------------------------------------------------------*/

scalar maxVelocityMag = 0.0;
scalar minCellSize = GREAT;

// Find maximum particle velocity on this processor
forAllConstIter(materialPointCloud, particles, iter)
{
    scalar vmag = mag(iter().velocity());
    maxVelocityMag = max(maxVelocityMag, vmag);
}

// Find minimum cell size on this processor
forAll(mesh.V(), celli)
{
    scalar cellSize = Foam::pow(mesh.V()[celli], 1.0/3.0);
    minCellSize = min(minCellSize, cellSize);
}

// Reduce across processors to get global maximum velocity and minimum cell size
maxVelocityMag = returnReduce(maxVelocityMag, maxOp<scalar>());
minCellSize = returnReduce(minCellSize, minOp<scalar>());

scalar CoNum = 0.0;
scalar meanCoNum = 0.0;

if (minCellSize > SMALL && maxVelocityMag > SMALL)
{
    CoNum = maxVelocityMag * runTime.deltaT().value() / minCellSize;
    meanCoNum = CoNum;  // For MPM, max and mean are similar
}

Info<< "MPM Courant Number: mean = " << meanCoNum 
    << "  max = " << CoNum << endl;

// ************************************************************************* //










































// // Calculate Courant number for MPM based on particle velocities and wave speed

// scalar maxVelocity = 0.0;
// scalar minCellSize = GREAT;

// // Find maximum particle velocity
// forAllConstIter(materialPointCloud, particles, iter)
// {
//     const MaterialPointParcel& p = iter();
//     scalar vmag = mag(p.velocity());
//     maxVelocity = max(maxVelocity, vmag);
// }

// // Find minimum cell size
// forAll(mesh.C(), celli)
// {
//     scalar cellSize = Foam::pow(mesh.V()[celli], 1.0/3.0);
//     minCellSize = min(minCellSize, cellSize);
// }

// // Reduce across processors
// reduce(maxVelocity, maxOp<scalar>());
// reduce(minCellSize, minOp<scalar>());

// // Calculate Courant number
// scalar CoNum = 0.0;
// scalar meanCoNum = 0.0;

// if (minCellSize > SMALL && maxVelocity > SMALL)
// {
//     CoNum = maxVelocity * runTime.deltaT().value() / minCellSize;
    
//     // Also consider wave speed for elastic materials
//     // c = sqrt(E/rho) for elastic wave speed
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
//                 // E = lambda + 2*mu for isotropic materials
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
            
//             scalar rho = mat.rho_0;
//             if (E > SMALL && rho > SMALL)
//             {
//                 waveSpeed = max(waveSpeed, Foam::sqrt(E/rho));
//             }
//         }
//     }
    
//     if (waveSpeed > SMALL)
//     {
//         scalar waveCo = waveSpeed * runTime.deltaT().value() / minCellSize;
//         CoNum = max(CoNum, waveCo);
//     }
    
//     meanCoNum = CoNum; // Simplified - could calculate mean
// }

// Info<< "MPM Courant Number: mean = " << meanCoNum 
//     << " max = " << CoNum << endl;

// // ************************************************************************* //