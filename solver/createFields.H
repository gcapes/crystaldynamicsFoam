#include "createRDeltaT.H"

Info<< "Reading field p_rgh\n" << endl;
volScalarField p_rgh
(
    IOobject
    (
        "p_rgh",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field U\n" << endl;
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

volTensorField omega_dot
(
    IOobject
    (
        "omega_dot",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
        mesh,
dimensionedTensor("omega_dot",dimensionSet(0, 0, -1, 0, 0),tensor::zero)
);

volTensorField epsilon_dot
(
    IOobject
    (
        "epsilon_dot",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
        mesh,
dimensionedTensor("epsilon_dot",dimensionSet(0, 0, -1, 0, 0),tensor::zero)
);

volTensorField s
(
    IOobject
    (
        "s",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
        mesh,
dimensionedTensor("s",dimensionSet(1, -1, -2, 0, 0),tensor::zero)
);

volTensorField dsdt
(
    IOobject
    (
        "dsdt",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
        mesh,
dimensionedTensor("dsdt",dimensionSet(1, -1, -3, 0, 0),tensor::zero)
);



#include "createPhi.H"

Info<< "Constructing compressibleTwoPhaseMixture\n" << endl;
compressibleTwoPhaseMixture mixture(U, phi);

autoPtr<twoPhaseChangeModel> phaseChangePtr
(
    twoPhaseChangeModel::New(mixture)
);

volScalarField& alpha1(mixture.alpha1());
const volScalarField& rho = mixture.rho();

dimensionedScalar pMin
(
    "pMin",
    dimPressure,
    mixture
);

mesh.schemes().setFluxRequired(p_rgh.name());
mesh.schemes().setFluxRequired(alpha1.name());


#include "readGravitationalAcceleration.H"
#include "readhRef.H"
#include "gh.H"


pressureReference pressureReference(mixture.p(), p_rgh, pimple.dict(), false);


// Mass flux
// Initialisation does not matter because rhoPhi is reset after the
// alpha1 solution before it is used in the U equation.
surfaceScalarField rhoPhi
(
    IOobject
    (
        "rhoPhi",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    fvc::interpolate(rho)*phi
);

volScalarField::Internal dgdt
(
    IOobject
    (
        "dgdt",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    alpha1*fvc::div(phi)
);

#include "createAlphaFluxes.H"

// Construct compressible turbulence model
compressibleInterPhaseTransportModel turbulence
(
    rho,
    U,
    phi,
    rhoPhi,
    alphaPhi1,
    mixture
);

Info<< "Creating field kinetic energy K\n" << endl;
volScalarField K("K", 0.5*magSqr(U));

#include "createMRF.H"
#include "createFvModels.H"
#include "createFvConstraints.H"
