Info << "Reading physicalProperties\n"
     << endl;

IOdictionary physicalProperties(
    IOobject(
        "physicalProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE));

IOdictionary solutionDict(
    IOobject(
        "fvSolution",
        runTime.system(), // folder where the fvSolution file is located
        runTime,
        IOobject::MUST_READ));

IOdictionary controlDict(
    IOobject(
        "controlDict",
        runTime.system(), // folder where the fvSolution file is located
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE));

dictionary &PimpleDICT = solutionDict.subDict("PIMPLE");

int n_sub_steps(readLabel(PimpleDICT.lookup("nOuterCorrectors")));

// scalar maxTimeStep(readScalar(controlDict.lookup("maxDeltaT")));

// Info<<"TREST: "<<n_sub_steps<<endl;

// exit(0);

dimensionedScalar pi("pi", dimensionSet(0, 0, 0, 0, 0), constant::mathematical::pi);

dimensionedScalar Cv("Cv", dimensionSet(0, 2, -2, -1, 0), physicalProperties.lookup("Cv"));
dimensionedScalar kappa("kappa", dimensionSet(1, 1, -3, -1, 0), physicalProperties.lookup("kappa"));
dimensionedScalar rho_0("rho_0", dimensionSet(1, -3, 0, 0, 0), physicalProperties.lookup("rho_0"));

//  Info<<"TEST"<<(n_alpha[0] * m_alpha[0])<<endl;

// exit(0);



Info << "Reading field U\n"
     << endl;
volVectorField U(
    IOobject(
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE),
    mesh);



Info << "Reading field Temperature\n"
     << endl;
volScalarField T(
    IOobject(
        "T",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE),
    mesh);

volScalarField K(
    IOobject(
        "K",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE),
    mesh,
    dimensionedScalar("K", dimensionSet(0, 2, -2, 0, 0), 0.0));



volScalarField rho(
    IOobject(
        "rho",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE),
    mesh,
    dimensionedScalar("rho", dimensionSet(1, -3, 0, 0, 0), rho_0.value()),
    zeroGradientFvPatchScalarField::typeName);

volScalarField p(
    IOobject(
        "p",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE),
    mesh,
    dimensionedScalar("p", dimensionSet(1, -1, -2, 0, 0), 0.0),
    zeroGradientFvPatchScalarField::typeName);

surfaceScalarField rhoPhi(
    IOobject(
        "rhoPhi",
        runTime.timeName(),
        mesh),
    mesh,
    dimensionedScalar(dimMass / dimTime, 0));



volTensorField omega_dot(
    IOobject(
        "omega_dot",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE),
    mesh,
    dimensionedTensor("omega_dot", dimensionSet(0, 0, -1, 0, 0), tensor::zero)//,
    // zeroGradientFvPatchTensorField::typeName
    );

volTensorField epsilon_dot(
    IOobject(
        "epsilon_dot",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE),
    mesh,
    dimensionedTensor("epsilon_dot", dimensionSet(0, 0, -1, 0, 0), tensor::zero)//,
    // calculatedFvPatchTensorField::typeName
    );


volTensorField sigma(
    IOobject(
        "sigma",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE),
    mesh,
    dimensionedTensor("sigma", dimensionSet(1, -1, -2, 0, 0), tensor::zero)//,
    // zeroGradientFvPatchTensorField::typeName
    );


volVectorField divsigma(
    IOobject(
        "divsigma",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE),
    mesh,
    dimensionedVector("divsigma", dimensionSet(1, -2, -2, 0, 0), vector::zero) //,
                                                                               // zeroGradientFvPatchTensorField::typeName
);


volTensorField dsigmadt(
    IOobject(
        "dsigmadt",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE),
    mesh,
    dimensionedTensor("dsigmadt", dimensionSet(1, -1, -3, 0, 0), tensor::zero)//,
    // zeroGradientFvPatchTensorField::typeName
    );

// volTensorField dsigmadt_crystal
// (
//     IOobject
//     (
//         "dsigmadt_crystal",
//         runTime.timeName(),
//         mesh,
//         IOobject::READ_IF_PRESENT,
//         IOobject::AUTO_WRITE
//     ),
//         mesh,
// dimensionedTensor("dsigmadt_crystal",dimensionSet(1, -1, -3, 0, 0),tensor::zero),
//                 zeroGradientFvPatchTensorField::typeName
// );

#include "createPhi.H"


