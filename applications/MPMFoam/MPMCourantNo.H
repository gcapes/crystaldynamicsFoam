// Calculate Courant number for MPM based on particle velocities and wave speed

scalar maxVelocity = 0.0;
scalar minCellSize = GREAT;

// Find maximum particle velocity
forAllConstIter(materialPointCloud, particles, iter)
{
    const MaterialPointParcel& p = iter();
    scalar vmag = mag(p.velocity());
    maxVelocity = max(maxVelocity, vmag);
}

// Find minimum cell size
forAll(mesh.C(), celli)
{
    scalar cellSize = Foam::pow(mesh.V()[celli], 1.0/3.0);
    minCellSize = min(minCellSize, cellSize);
}

// Reduce across processors
reduce(maxVelocity, maxOp<scalar>());
reduce(minCellSize, minOp<scalar>());

// Calculate Courant number
scalar CoNum = 0.0;
scalar meanCoNum = 0.0;

if (minCellSize > SMALL && maxVelocity > SMALL)
{
    CoNum = maxVelocity * runTime.deltaT().value() / minCellSize;
    
    // Also consider wave speed for elastic materials
    // c = sqrt(E/rho) for elastic wave speed
    // This is a simplified estimate
    scalar waveSpeed = 0.0;
    
    forAll(materialProperties, matI)
    {
        if (materialProperties.set(matI))
        {
            scalar E = materialProperties[matI].lookupOrDefault<scalar>("E", 1e9);
            scalar rho = materialProperties[matI].lookupOrDefault<scalar>("rho", 1000);
            waveSpeed = max(waveSpeed, Foam::sqrt(E/rho));
        }
    }
    
    if (waveSpeed > SMALL)
    {
        scalar waveCo = waveSpeed * runTime.deltaT().value() / minCellSize;
        CoNum = max(CoNum, waveCo);
    }
    
    meanCoNum = CoNum; // Simplified - could calculate mean
}

Info<< "MPM Courant Number: mean = " << meanCoNum 
    << " max = " << CoNum << endl;
