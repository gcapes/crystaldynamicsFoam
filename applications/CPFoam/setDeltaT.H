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
    MODIFIED FOR SOLID DYNAMICSS

\*---------------------------------------------------------------------------*/

if (adjustTimeStep)
{

    const scalar CoLimitingEOS =
        useEOS ? max(CoNum, acousticCoNum) : CoNum;

    if (useEOS)
    {
        Info<< "Time-step limiting Courant number (acoustic): "
            << CoLimitingEOS << endl;
    }

    scalar deltaT = maxCo*runTime.deltaTValue()/(CoLimitingEOS + SMALL);

    if (fracture)
    {
        static scalar deltaTPFdiff = -1.0;
        if (deltaTPFdiff < 0.0)
        {
            scalar MEmax = 0.0;
            for (label a = 0; a < N_Mats; a++)
            {
                for (label b = 0; b < N_Mats; b++)
                {
                    MEmax = max
                    (
                        MEmax,
                        PF_Mobility0_array[a][b]*PF_Energy0_array[a][b]
                    );
                }
            }
            const scalar dxMin = gMin(cellDx.primitiveField());
            deltaTPFdiff =
                0.5*sqr(dxMin)*sqr(pi.value())
               /(6.0*4.0*max(MEmax, SMALL)*interface_width.value());
            Info<< "PF stability: diffusive dt bound = "
                << deltaTPFdiff << " s" << endl;
        }

        const scalar maxDPhi
        (
            GrainProperties.lookupOrDefault<scalar>("PF_maxDeltaPhi", 0.05)//MaX CHANGE IN PHI in a  GIVEN TIMESTEP
        );
        scalar dpdtMax = 0.0;
        for (label h = 0; h < N_Ori; h++)
        {
            dpdtMax = max
            (
                dpdtMax,
                gMax(mag(dpdt[h])().primitiveField())
            );
        }
        deltaT = min
        (
            deltaT,
            min(deltaTPFdiff, maxDPhi/max(dpdtMax, SMALL))
        );
    }

    deltaT = min(deltaT, maxTimeStep/*fvModels.maxDeltaT()*/);
    deltaT = min
    (
        min(deltaT, runTime.deltaTValue() + 0.1*deltaT),
        1.2*runTime.deltaTValue()
    );
    runTime.setDeltaT(min(deltaT, maxDeltaT));

    Info<< "deltaT = " <<  runTime.deltaTValue() << endl;
}

// ************************************************************************* //