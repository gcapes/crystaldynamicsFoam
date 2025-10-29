// Solve momentum equation on grid (explicit for now)

scalar dt = runTime.deltaT().value();

forAll(mesh.C(), celli)
{
    if (gridMass[celli] > SMALL)
    {
        // Explicit update: v_new = v_old + dt * F / m
        gridVelocity[celli] += dt * gridForce[celli] * mesh.V()[celli] / gridMass[celli];
    }
}

// Apply boundary conditions
gridVelocity.correctBoundaryConditions();

// Compute velocity gradient for G2P
gridVelocityGradient = fvc::grad(gridVelocity);
