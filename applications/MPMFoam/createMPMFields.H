// Create nodal (point) fields for MPM

Info<< "Creating nodal fields..." << endl;

// Nodal mass field
pointScalarField nodalMass
(
    IOobject
    (
        "nodalMass",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    pMesh,
    dimensionedScalar("zero", dimMass, 0.0)
);

// Nodal momentum field
pointVectorField nodalMomentum
(
    IOobject
    (
        "nodalMomentum",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    pMesh,
    dimensionedVector("zero", dimMass*dimVelocity, Zero)
);

// Nodal velocity field
pointVectorField nodalVelocity
(
    IOobject
    (
        "nodalVelocity",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    pMesh,
    dimensionedVector("zero", dimVelocity, Zero)
);

// Nodal acceleration field
pointVectorField nodalAcceleration
(
    IOobject
    (
        "nodalAcceleration",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    pMesh,
    dimensionedVector("zero", dimVelocity/dimTime, Zero)
);

// Nodal force field
pointVectorField nodalForce
(
    IOobject
    (
        "nodalForce",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    pMesh,
    dimensionedVector("zero", dimForce, Zero)
);

// Nodal stress field (for visualization/debugging)
pointTensorField nodalStress
(
    IOobject
    (
        "nodalStress",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    pMesh,
    dimensionedTensor("zero", dimPressure, Zero)
);

Info<< "Nodal fields created" << endl;
