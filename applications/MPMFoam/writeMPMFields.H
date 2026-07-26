// Write MPM particle data to cell fields for visualization

// Create cell-centered fields from particle data
volScalarField cellMass
(
    IOobject
    (
        "cellMass",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("zero", dimMass, 0.0)
);

volVectorField cellVelocity
(
    IOobject
    (
        "cellVelocity",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedVector("zero", dimVelocity, Zero)
);

volTensorField cellStress
(
    IOobject
    (
        "cellStress",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedTensor("zero", dimPressure, Zero)
);

volTensorField cellF
(
    IOobject
    (
        "cellF",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedTensor("zero", dimless, tensor::I)
);

volScalarField cellJ
(
    IOobject
    (
        "cellJ",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("zero", dimless, 1.0)
);

// Accumulate particle data to cells
scalarField cellParticleCount(mesh.nCells(), 0.0);

forAllConstIter(materialPointCloud, particles, iter)
{
    const MaterialPointParcel& p = iter();
    label celli = p.cell();
    
    if (celli >= 0 && celli < mesh.nCells())
    {
        cellMass[celli] += p.mass();
        cellVelocity[celli] += p.mass() * p.velocity();
        cellStress[celli] += p.stress() * p.volume();
        cellF[celli] += p.F();
        cellParticleCount[celli] += 1.0;
    }
}

// Average the accumulated values
forAll(mesh.C(), celli)
{
    if (cellParticleCount[celli] > 0)
    {
        cellVelocity[celli] /= cellMass[celli] + SMALL;
        cellStress[celli] /= (cellParticleCount[celli] * mesh.V()[celli]);
        cellF[celli] /= cellParticleCount[celli];
        cellJ[celli] = det(cellF[celli]);
    }
}

cellMass.correctBoundaryConditions();
cellVelocity.correctBoundaryConditions();
cellStress.correctBoundaryConditions();
cellF.correctBoundaryConditions();
cellJ.correctBoundaryConditions();

// Calculate von Mises stress for visualization
volScalarField vonMisesStress
(
    IOobject
    (
        "vonMisesStress",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("zero", dimPressure, 0.0)
);

forAll(cellStress, celli)
{
    // Deviatoric stress
    tensor s = cellStress[celli] - (1.0/3.0) * tr(cellStress[celli]) * tensor::I;
    
    // Von Mises stress = sqrt(3/2 * s:s)
    vonMisesStress[celli] = Foam::sqrt(1.5 * (s && s));
}

vonMisesStress.correctBoundaryConditions();

Info<< "Wrote MPM fields" << endl;