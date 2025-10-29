// Map particle quantities to grid (P2G)

forAllConstIter(materialPointCloud, particles, iter)
{
    const MaterialPointParcel& p = iter();
    p.mapToGrid(gridStress, gridMomentum, gridMass);
}

// Normalize and compute velocity
forAll(mesh.C(), celli)
{
    if (gridMass[celli] > SMALL)
    {
        gridVelocity[celli] = gridMomentum[celli] / gridMass[celli];
        gridStress[celli] /= mesh.V()[celli]; // Normalize stress by cell volume
    }
    else
    {
        gridVelocity[celli] = vector::zero;
        gridStress[celli] = tensor::zero;
    }
}

gridVelocity.correctBoundaryConditions();
