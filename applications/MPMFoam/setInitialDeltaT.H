// Set initial time step

// Read time step control parameters if not already read
scalar maxCo = 0.5;
scalar maxDeltaT = GREAT;
bool adjustTimeStep = true;

if (mesh.solutionDict().found("PIMPLE"))
{
    const dictionary& pimple = mesh.solutionDict().subDict("PIMPLE");
    maxCo = pimple.lookupOrDefault<scalar>("maxCo", 0.5);
    maxDeltaT = pimple.lookupOrDefault<scalar>("maxDeltaT", GREAT);
    adjustTimeStep = pimple.lookupOrDefault<Switch>("adjustTimeStep", true);
}

if (adjustTimeStep)
{
    // Calculate initial time step based on Courant number
    // Use a conservative initial time step
    scalar initialDeltaT = maxDeltaT;
    
    // Estimate wave speed from materials if available
    if (materialProperties.size() > 0 && materialProperties.set(0))
    {
        scalar E = materialProperties[0].lookupOrDefault<scalar>("E", 1e9);
        scalar rho = materialProperties[0].lookupOrDefault<scalar>("rho", 1000);
        scalar waveSpeed = Foam::sqrt(E/rho);
        
        // Find minimum cell size
        scalar minCellSize = GREAT;
        forAll(mesh.V(), celli)
        {
            scalar cellSize = Foam::pow(mesh.V()[celli], 1.0/3.0);
            minCellSize = min(minCellSize, cellSize);
        }
        
        reduce(minCellSize, minOp<scalar>());
        
        // Calculate time step based on wave speed and Courant number
        if (minCellSize > SMALL && waveSpeed > SMALL)
        {
            initialDeltaT = min(maxDeltaT, maxCo * minCellSize / waveSpeed);
        }
    }
    
    runTime.setDeltaT(initialDeltaT);
}

Info<< "Initial time step: " << runTime.deltaT().value() << endl;
