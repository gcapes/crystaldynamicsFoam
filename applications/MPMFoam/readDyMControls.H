// Read dynamic mesh control parameters
// (Simplified for MPM - mostly for time step control)

const dictionary& pimple = mesh.solutionDict().subDict("PIMPLE");

// Maximum Courant number
scalar maxCo = pimple.lookupOrDefault<scalar>("maxCo", 0.5);

// Maximum time step
scalar maxDeltaT = pimple.lookupOrDefault<scalar>("maxDeltaT", GREAT);

// Adjustable time step
bool adjustTimeStep = pimple.lookupOrDefault<Switch>("adjustTimeStep", true);
