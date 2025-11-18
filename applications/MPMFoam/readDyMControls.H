

const dictionary& controlDict = runTime.controlDict();
// Maximum Courant number
scalar maxCo = controlDict.lookupOrDefault<scalar>("maxCo", 0.1);

// Maximum time step
scalar maxDeltaT = controlDict.lookupOrDefault<scalar>("maxDeltaT", GREAT);

// Adjustable time step
bool adjustTimeStep = controlDict.lookupOrDefault<Switch>("adjustTimeStep", true);
