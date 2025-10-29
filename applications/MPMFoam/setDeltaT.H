/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2022 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Global
    setDeltaT

Description
    Reset the timestep to maintain a constant maximum courant Number.
    Reduction of time-step is immediate, but increase is damped to avoid
    unstable oscillations.

\*---------------------------------------------------------------------------*/

// Adjust time step based on Courant number

if (adjustTimeStep && minCellSize > SMALL && maxVelocity > SMALL)
{
    // Consider both velocity and wave speed
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
    
    // Use the more restrictive of velocity-based and wave-speed-based time step
    scalar maxSpeed = max(maxVelocity, waveSpeed);
    
    scalar deltaT = min
    (
        maxDeltaT,
        maxCo * minCellSize / (maxSpeed + SMALL)
    );
    
    // Limit the rate of change of time step
    scalar deltaTFact = min(max(deltaT/runTime.deltaT().value(), 0.5), 2.0);
    
    runTime.setDeltaT(deltaTFact * runTime.deltaT().value());
    
    Info<< "Adjusted time step to: " << runTime.deltaT().value() << endl;
}

// ************************************************************************* //
