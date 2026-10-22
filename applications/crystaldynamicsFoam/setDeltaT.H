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

    // Acoustic Courant limit of the explicit stress coupling (both closures); off for
    // quasiStatic, for acousticCourantControl false and in the prescribed-stress mode.
    const bool acousticLimitActive =
        stressUpdate && (useEOS || acousticCourantControl) && !quasiStatic;

    const scalar CoLimitingEOS =
        acousticLimitActive ? max(CoNum, acousticCoNum) : CoNum;

    if (acousticLimitActive)
    {
        Info<< "Time-step limiting Courant number (acoustic): "
            << CoLimitingEOS << endl;
    }

    scalar deltaT = maxCo*runTime.deltaTValue()/(CoLimitingEOS + SMALL);

    if (PF_drivingForceActive)
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
            GrainProperties.lookupOrDefault<scalar>("PF_maxDeltaPhi", 0.05) // max change of any phi per step
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

    if (quasiStatic)
    {
        // Strain-increment bound from the interior strain rate and the boundary-face
        // velocity jump (active from the first step; coded BCs evaluated first).
        if (maxDeltaEps > 0.0)
        {
            U.correctBoundaryConditions();
            scalar epsDotMax = gMax(mag(epsilon_dot)().primitiveField());
            scalar epsDotBoundary = 0.0;
            forAll(U.boundaryField(), patchi)
            {
                const fvPatchVectorField& Up = U.boundaryField()[patchi];
                if (!Up.coupled() && Up.type() != "empty" && Up.size() > 0)
                {
                    epsDotBoundary = max(epsDotBoundary, max(mag(Up.snGrad())()));
                }
            }
            reduce(epsDotBoundary, maxOp<scalar>());
            epsDotMax = max(epsDotMax, epsDotBoundary);
            if (epsDotMax > SMALL)
            {
                deltaT = min(deltaT, maxDeltaEps/epsDotMax);
            }
            Info<< "quasiStatic: max strain rate " << epsDotMax
                << " 1/s (boundary faces " << epsDotBoundary
                << " 1/s), strain-increment dt bound "
                << (epsDotMax > SMALL ? maxDeltaEps/epsDotMax : GREAT) << " s" << endl;

            // A strain-increment bound below the acoustic step means the loading is not
            // quasi-static: warn.
            static scalar deltaTAcousticQS = -1.0;
            if (deltaTAcousticQS < 0.0)
            {
                scalar cLmaxQS = 0.0;
                for (label i = 0; i < N_Mats; i++)
                {
                    const scalar GammaQS =
                        GammaMaxMat[i]
                      + (
                            useEOS
                          ? max(0.0, MatList[i].rho_0*sqr(MatList[i].c0) - KCprimeMat[i])
                          : 0.0
                        );
                    cLmaxQS = max
                    (
                        cLmaxQS,
                        Foam::sqrt(GammaQS/max(MatList[i].rho_0, SMALL))
                    );
                }
                deltaTAcousticQS =
                    gMin(cellDx.primitiveField())/max(cLmaxQS, SMALL);
            }
            if (epsDotMax > SMALL && maxDeltaEps/epsDotMax < deltaTAcousticQS)
            {
                WarningInFunction
                    << "quasiStatic: the strain-increment bound "
                    << maxDeltaEps/epsDotMax << " s is below the acoustic time step "
                    << deltaTAcousticQS << " s. A strain rate of " << epsDotMax
                    << " 1/s is not quasi-static: check the velocity boundary"
                    << " conditions in 0/U, or run with momentumCoupling explicit." << endl;
            }
        }

        // Viscoplastic bound lambda*dt <= maxViscoplasticNumber (few sub-steps, well
        // conditioned outer iteration)
        if (maxViscoplasticNumber > 0.0)
        {
            scalar lamMax = 0.0;
            for (label h = 0; h < N_Ori; h++)
            {
                lamMax = max(lamMax, gMax(lamVp[h]));
            }
            if (lamMax > SMALL)
            {
                deltaT = min(deltaT, maxViscoplasticNumber/lamMax);
            }
            Info<< "quasiStatic: max viscoplastic rate " << lamMax
                << " 1/s, viscoplastic dt bound "
                << (lamMax > SMALL ? maxViscoplasticNumber/lamMax : GREAT) << " s" << endl;
        }
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