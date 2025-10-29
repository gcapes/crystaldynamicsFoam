// Compute internal forces from stress divergence

gridForce = -fvc::div(gridStress);

// Add body forces (gravity)
forAll(mesh.C(), celli)
{
    if (gridMass[celli] > SMALL)
    {
        gridForce[celli] += (gridMass[celli] / mesh.V()[celli]) * g.value();
    }
}
