/*---------------------------------------------------------------------------*\
  setInitialDeltaT.H - UPDATED FOR MATERIAL CLASS
\*---------------------------------------------------------------------------*/

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
    if (materials.size() > 0 && materials.set(0))
    {
        const Material& mat = materials[0];
        
        // Get elastic modulus
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
            
            // Find minimum cell size
            scalar minCellSize = GREAT;
            forAll(mesh.V(), celli)
            {
                scalar cellSize = Foam::pow(mesh.V()[celli], 1.0/3.0);
                minCellSize = min(minCellSize, cellSize);
            }
            
            reduce(minCellSize, minOp<scalar>());
            
            // Calculate time step based on wave speed and Courant number
            if (minCellSize > SMALL)
            {
                initialDeltaT = min(maxDeltaT, maxCo * minCellSize / waveSpeed);
            }
        }
    }
    
    runTime.setDeltaT(initialDeltaT);
}

Info<< "Initial time step: " << runTime.deltaT().value() << endl;

// ************************************************************************* //