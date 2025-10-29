/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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

Class
    Foam::multiphaseMixture

Description


SourceFiles
    multiphaseMixture.C

\*---------------------------------------------------------------------------*/

// Map particle quantities to grid
void mapParticlesToGrid
(
    const fvMesh& mesh,
    const materialPointCloud& mpCloud,
    volTensorField& gridStress,
    volVectorField& gridMomentum,
    volScalarField& gridMass
)
{
    // Reset grid quantities
    gridStress = dimensionedTensor("zero", gridStress.dimensions(), tensor::zero);
    gridMomentum = dimensionedVector("zero", gridMomentum.dimensions(), vector::zero);
    gridMass = dimensionedScalar("zero", gridMass.dimensions(), 0.0);
    
    forAll(gridPhases, h)
    {
        gridPhases[h] = dimensionedScalar("zero", dimless, 0.0);
    }
    
    // Accumulate from particles
    forAllConstIter(materialPointCloud, mpCloud, iter)
    {
        const MaterialPointParcel& mp = iter();
        label celli = mp.cell();
        
        scalar weight = 1.0;
        
        gridStress[celli] += weight * mp.stress() * mp.volume();
        gridMomentum[celli] += weight * mp.mass() * mp.velocity();
        gridMass[celli] += weight * mp.mass();
    }
    
    // Normalize stress by volume
    forAll(mesh.C(), celli)
    {
        if(gridMass[celli] > SMALL)
        {
            gridStress[celli] /= mesh.cellVolumes()[celli];
        }
    }
}

void mapGridToParticles
(
    const fvMesh& mesh,
    materialPointCloud& mpCloud,
    const volVectorField& gridVelocity,
    const volTensorField& velocityGradient,
    const scalar deltaT
)
{
    forAllIter(materialPointCloud, mpCloud, iter)
    {
        MaterialPointParcel& mp = iter();
        label celli = mp.cell();
        
        // Update particle velocity
        mp.velocity() = gridVelocity[celli];
        
        // Update deformation gradient
        tensor L = velocityGradient[celli];
        mp.F() = mp.F() + deltaT * (L & mp.F());
    }
}
