
// ==================== DICTIONARIES ====================

IOdictionary solutionDict(
    IOobject(
        "fvSolution",
        runTime.system(),
        runTime,
        IOobject::MUST_READ));

IOdictionary controlDict(
    IOobject(
        "controlDict",
        runTime.system(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE));

dictionary &PimpleDICT = solutionDict.subDict("PIMPLE");

scalar pimple_relax_U = PimpleDICT.lookupOrDefault<scalar>("relaxU", 0.5);

dimensionedScalar pi("pi", dimensionSet(0, 0, 0, 0, 0), constant::mathematical::pi);


// ==================== PRIMARY FIELDS ====================

Info << "Reading field U\n" << endl;
volVectorField U(
    IOobject(
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE),
    mesh);

volVectorField D(
    IOobject(
        "D",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE),
    mesh,
    dimensionedVector("D", dimensionSet(0, 1, 0, 0, 0), vector::zero));

Info << "Reading field Temperature\n" << endl;
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

volScalarField Sum_Phases(
    IOobject(
        "Sum_Phases",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE),
    mesh,
    dimensionedScalar("Sum_Phases", dimensionSet(0, 0, 0, 0, 0), 0.0));


// ==================== GRAIN / MATERIAL PROPERTIES ====================

IOdictionary GrainProperties(
    IOobject(
        "GrainProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE));

label N_Ori(readLabel(GrainProperties.lookup("N_Phases")));
Info << "Number of Phase Field Order Parameters: " << N_Ori << endl;

label N_Mats(readLabel(GrainProperties.lookup("N_Materials")));
Info << "Number of Materials: " << N_Mats << endl;

label n_substeps(readLabel(GrainProperties.lookup("n_constitutive_substeps")));
Info << "Number of Constitutive Sub-cycling steps: " << n_substeps << endl;

labelList Phase2Material(N_Ori, 0);
GrainProperties.lookup("Phase2Material") >> Phase2Material;

for (label c = 0; c < N_Ori; c++)
{
    Info << "Phase[" << c << "] = Material[" << Phase2Material[c] << "]" << endl;
}

// Stiffness tensors in crystal coordinates (one per material)
Fastor::Tensor<double, 3, 3, 3, 3> Cprime2[N_Mats];

// Material objects
std::vector<Material> MatList(N_Mats, Material());

for (label i = 0; i < N_Mats; i++)
{
    word Mat_dict("Material." + name(i));
    dictionary Mat = GrainProperties.subDict(Mat_dict);

    // Read stiffness tensor
    Fastor::Tensor<double, 3, 3, 3, 3> Cprime_temp_Read(0.0);
    for (int m = 0; m < 3; m++)
        for (int n = 0; n < 3; n++)
            for (int p = 0; p < 3; p++)
                for (int q = 0; q < 3; q++)
                {
                    word Cijkl("C" + name(m) + name(n) + name(p) + name(q));
                    Cprime_temp_Read(m, n, p, q) = Mat.lookupOrDefault<scalar>(Cijkl, 0.0);
                }
    Cprime2[i] = Cprime_temp_Read;

    // Thermal properties
    MatList[i].rho_0  = readScalar(Mat.lookup("rho_0"));
    MatList[i].Cv     = readScalar(Mat.lookup("Cv"));
    MatList[i].kappa  = readScalar(Mat.lookup("kappa"));

    // Slip systems
    label N_slips_Temp(readLabel(Mat.lookup("N_slips")));
    MatList[i].N_Slips = N_slips_Temp;

    for (int j = 0; j < MatList[i].N_Slips; j++)
    {
        word Slip_direction_name = "m" + name(j);
        word Slip_normal_name = "n" + name(j);
        vector temp_m(Mat.lookup(Slip_direction_name));
        vector temp_n(Mat.lookup(Slip_normal_name));
        MatList[i].SP_norm.push_back(temp_n);
        MatList[i].SP_direction.push_back(temp_m);
    }

    // Hardening parameters
    MatList[i].initial_CRSS   = readScalar(Mat.lookup("initialCRSS"));
    MatList[i].tau_S           = readScalar(Mat.lookup("tau_S"));
    MatList[i].h_0             = readScalar(Mat.lookup("h_0"));
    MatList[i].r               = readScalar(Mat.lookup("r"));
    MatList[i].k_exponent      = readScalar(Mat.lookup("k_exponent"));
    MatList[i].ref_shear_rate  = readScalar(Mat.lookup("gammadot_0"));
}


scalar D_stab_value = 0.0;
for (label i = 0; i < N_Mats; i++)
{
    scalar C11 = Cprime2[i](0, 0, 0, 0);
    scalar C22 = Cprime2[i](1, 1, 1, 1);
    scalar C33 = Cprime2[i](2, 2, 2, 2);
    D_stab_value = max(D_stab_value, max(C11, max(C22, C33)));
}
dimensionedScalar D_stab("D_stab", dimensionSet(1, -1, -2, 0, 0), D_stab_value);
Info << "Stabilization diffusivity D_stab = " << D_stab_value << " Pa" << endl;

// Euler angles (used only for initializing Fp)
List<double> phi1 = GrainProperties.lookup("phi1");
List<double> phi2 = GrainProperties.lookup("phi2");
List<double> phi3 = GrainProperties.lookup("phi3");


// ==================== PHASE FIELDS ====================

PtrList<volScalarField> PopBal(N_Ori);
for (label count = 0; count < N_Ori; count++)
{
    PopBal.set(
        count,
        new volScalarField(
            IOobject(
                "n." + name(count),
                runTime.timeName(),
                mesh,
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE),
            mesh,
            dimensionedScalar("ni", dimensionSet(0, 0, 0, 0, 0), 0.0),
            zeroGradientFvPatchScalarField::typeName));
}

PtrList<volScalarField> PopBal_av(N_Ori);
for (label count = 0; count < N_Ori; count++)
{
    PopBal_av.set(
        count,
        new volScalarField(
            IOobject(
                "n_av" + name(count),
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE),
            mesh,
            dimensionedScalar("niav", dimensionSet(0, 0, 0, 0, 0), 0.0),
            zeroGradientFvPatchScalarField::typeName));
}

// Normalize phase fields onto Gibbs simplex
Sum_Phases *= 0.0;
for (label count = 0; count < N_Ori; count++)
{
    Sum_Phases += PopBal[count];
}
Sum_Phases.correctBoundaryConditions();

for (label count = 0; count < N_Ori; count++)
{
    volScalarField &ni = PopBal[count];
    ni = max(min(ni, 1.0), 0.0);
    volScalarField temp_holder = ni;

    forAll(mesh.C(), celli)
    {
        if (Sum_Phases[celli] > 0.0)
        {
            temp_holder[celli] *= (1.0 / Sum_Phases[celli]);
        }
    }

    forAll(Sum_Phases.boundaryField(), patchI)
    {
        if (!Sum_Phases.boundaryField()[patchI].coupled()
            && Sum_Phases.boundaryField()[patchI].size() > 0)
        {
            forAll(Sum_Phases.boundaryField()[patchI], faceI)
            {
                if (Sum_Phases.boundaryField()[patchI][faceI] > 0.0)
                {
                    temp_holder.boundaryFieldRef()[patchI][faceI] *= 
                        (1.0 / Sum_Phases.boundaryField()[patchI][faceI]);
                }
            }
        }
    }
    ni = temp_holder;
}


// ==================== THERMAL MIXTURE PROPERTIES ====================

volScalarField rho_ref(
    IOobject(
        "rho_ref",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE),
    mesh,
    dimensionedScalar("rho_ref", dimensionSet(1, -3, 0, 0, 0), 0.0),
    zeroGradientFvPatchScalarField::typeName);

volScalarField kappa(
    IOobject(
        "kappa",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE),
    mesh,
    dimensionedScalar("kappa", dimensionSet(1, 1, -3, -1, 0), 0.0));

volScalarField Cv(
    IOobject(
        "Cv",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE),
    mesh,
    dimensionedScalar("Cv", dimensionSet(0, 2, -2, -1, 0), 0.0));

// Initialize mixture properties from phase fields
for (label h = 0; h < N_Ori; h++)
{
    dimensionedScalar dimkappa("dimkappa", dimensionSet(1, 1, -3, -1, 0), 1.0);
    dimensionedScalar dimcv("dimcv", dimensionSet(0, 2, -2, -1, 0), 1.0);
    dimensionedScalar dimrho("dimrho", dimensionSet(1, -3, 0, 0, 0), 1.0);

    kappa   += PopBal[h] * MatList[Phase2Material[h]].kappa * dimkappa;
    Cv      += PopBal[h] * MatList[Phase2Material[h]].Cv * dimcv;
    rho_ref += PopBal[h] * MatList[Phase2Material[h]].rho_0 * dimrho;
}

Info << "Reading field rho\n" << endl;
volScalarField rho(
    IOobject(
        "rho",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE),
    rho_ref);




Info << "Reading field F (global deformation gradient)\n" << endl;
volTensorField F(
    IOobject(
        "F",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE),
    mesh);

// Jacobian J = det(F)
volScalarField J(
    IOobject(
        "J",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE),
    mesh,
    dimensionedScalar("J", dimensionSet(0, 0, 0, 0, 0), 1.0));




PtrList<volTensorField> F_plastic(N_Ori);
for (label count = 0; count < N_Ori; count++)
{
    F_plastic.set(
        count,
        new volTensorField(
            IOobject(
                "Fp." + Foam::name(count),
                runTime.timeName(),
                mesh,
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE),
            mesh,
            dimensionedTensor("Fp", dimensionSet(0, 0, 0, 0, 0), tensor::I),
            zeroGradientFvPatchTensorField::typeName));

    // If NOT reading from a restart, initialize Fp with crystal orientation
    if (!F_plastic[count].headerOk())
    {
        Info << "Initializing Fp." << count 
             << " with crystal orientation (phi1=" << phi1[count]
             << ", phi2=" << phi2[count] 
             << ", phi3=" << phi3[count] << ")" << endl;

        double p1 = phi1[count];
        double p2 = phi2[count];
        double p3 = phi3[count];

        // Build rotation matrix from Euler angles
        // Convention: extrinsic rotations about x, y, z axes
        tensor R0(
            Foam::cos(p2) * Foam::cos(p3),
            (Foam::sin(p1) * Foam::sin(p2) * Foam::cos(p3)) - (Foam::cos(p1) * Foam::sin(p3)),
            (Foam::cos(p1) * Foam::sin(p2) * Foam::cos(p3)) + (Foam::sin(p1) * Foam::sin(p3)),

            Foam::cos(p2) * Foam::sin(p3),
            (Foam::sin(p1) * Foam::sin(p2) * Foam::sin(p3)) + (Foam::cos(p1) * Foam::cos(p3)),
            (Foam::cos(p1) * Foam::sin(p2) * Foam::sin(p3)) - (Foam::sin(p1) * Foam::cos(p3)),

            -Foam::sin(p2),
            Foam::sin(p1) * Foam::cos(p2),
            Foam::cos(p1) * Foam::cos(p2)
        );

        // Set all cells
        forAll(mesh.C(), cell)
        {
            F_plastic[count][cell] = R0;
        }

        // Set all boundary faces
        forAll(F_plastic[count].boundaryField(), patchI)
        {
            if (!F_plastic[count].boundaryField()[patchI].coupled()
                && F_plastic[count].boundaryField()[patchI].size() > 0)
            {
                forAll(F_plastic[count].boundaryField()[patchI], faceI)
                {
                    F_plastic[count].boundaryFieldRef()[patchI][faceI] = R0;
                }
            }
        }
    }
}




PtrList<volTensorField> Lp(N_Ori);
for (label count = 0; count < N_Ori; count++)
{
    Lp.set(
        count,
        new volTensorField(
            IOobject(
                "Lp." + Foam::name(count),
                runTime.timeName(),
                mesh,
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE),
            mesh,
            dimensionedTensor("Lp", dimensionSet(0, 0, -1, 0, 0), tensor::zero),  // <-- ZERO!
            zeroGradientFvPatchTensorField::typeName));
}


// ==================== HARDENING FIELDS ====================

// g_alpha: critically resolved shear stress per slip system per phase
PtrList<PtrList<volScalarField>> g_alphas(N_Ori);
forAll(g_alphas, i)
{
    g_alphas.set(i, new PtrList<volScalarField>(MatList[Phase2Material[i]].N_Slips));
    forAll(g_alphas[i], j)
    {
        g_alphas[i].set(
            j,
            new volScalarField(
                IOobject(
                    "g_alpha." + Foam::name(i) + "." + Foam::name(j),
                    runTime.timeName(),
                    mesh,
                    IOobject::READ_IF_PRESENT,
                    IOobject::AUTO_WRITE),
                mesh,
                dimensionedScalar("g_alphas", dimless, MatList[Phase2Material[i]].initial_CRSS),
                zeroGradientFvPatchScalarField::typeName));
    }
}

// g_alpha_dots: hardening rates (RHS source terms)
PtrList<PtrList<volScalarField>> g_alpha_dots(N_Ori);
forAll(g_alpha_dots, i)
{
    g_alpha_dots.set(i, new PtrList<volScalarField>(MatList[Phase2Material[i]].N_Slips));
    forAll(g_alpha_dots[i], j)
    {
        g_alpha_dots[i].set(
            j,
            new volScalarField(
                IOobject(
                    "g_alpha_dots." + Foam::name(i) + "." + Foam::name(j),
                    runTime.timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE),
                mesh,
                dimensionedScalar("g_alpha_dots", dimensionSet(0, 0, -1, 0, 0), 0.0),
                zeroGradientFvPatchScalarField::typeName));
    }
}

// gamma_alpha_dots: slip rates per system
PtrList<PtrList<volScalarField>> gamma_alpha_dots(N_Ori);
forAll(gamma_alpha_dots, i)
{
    gamma_alpha_dots.set(i, new PtrList<volScalarField>(MatList[Phase2Material[i]].N_Slips));
    forAll(gamma_alpha_dots[i], j)
    {
        gamma_alpha_dots[i].set(
            j,
            new volScalarField(
                IOobject(
                    "gamma_alpha_dots." + Foam::name(i) + "." + Foam::name(j),
                    runTime.timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE),
                mesh,
                dimensionedScalar("gamma_alpha_dots", dimensionSet(0, 0, -1, 0, 0), 0.0),
                zeroGradientFvPatchScalarField::typeName));
    }
}

// gamma_alphas: accumulated slip per system (for hardening law)
PtrList<PtrList<volScalarField>> gamma_alphas(N_Ori);
forAll(gamma_alphas, i)
{
    gamma_alphas.set(i, new PtrList<volScalarField>(MatList[Phase2Material[i]].N_Slips));
    forAll(gamma_alphas[i], j)
    {
        gamma_alphas[i].set(
            j,
            new volScalarField(
                IOobject(
                    "gamma_alphas." + Foam::name(i) + "." + Foam::name(j),
                    runTime.timeName(),
                    mesh,
                    IOobject::READ_IF_PRESENT,
                    IOobject::AUTO_WRITE),
                mesh,
                dimensionedScalar("gamma_alphas", dimless, 0.0),
                zeroGradientFvPatchScalarField::typeName));
    }
}


// ==================== STRESS FIELDS ====================

// Per-phase Cauchy stress
PtrList<volTensorField> sigma(N_Ori);
for (label count = 0; count < N_Ori; count++)
{
    sigma.set(
        count,
        new volTensorField(
            IOobject(
                "sigma." + Foam::name(count),
                runTime.timeName(),
                mesh,
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE),
            mesh,
            dimensionedTensor("sigma", dimensionSet(1, -1, -2, 0, 0), tensor::zero),
            zeroGradientFvPatchTensorField::typeName));
}

// Temporary stress for sub-stepping
PtrList<volTensorField> sigtemp(N_Ori);
for (label count = 0; count < N_Ori; count++)
{
    sigtemp.set(
        count,
        new volTensorField(
            IOobject(
                "sigtemp." + Foam::name(count),
                runTime.timeName(),
                mesh,
                IOobject::READ_IF_PRESENT,
                IOobject::NO_WRITE),
            mesh,
            dimensionedTensor("sigtemp", dimensionSet(1, -1, -2, 0, 0), tensor::zero),
            zeroGradientFvPatchTensorField::typeName));
}

// Global (mixture) Cauchy stress
volTensorField sigmaGlobal(
    IOobject(
        "sigmaGlobal",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE),
    mesh,
    dimensionedTensor("sigmaGlobal", dimensionSet(1, -1, -2, 0, 0), tensor::zero),
    zeroGradientFvPatchTensorField::typeName);

volTensorField sigmaGlobaltemp(
    IOobject(
        "sigmaGlobaltemp",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE),
    mesh,
    dimensionedTensor("sigmaGlobaltemp", dimensionSet(1, -1, -2, 0, 0), tensor::zero),
    zeroGradientFvPatchTensorField::typeName);


// ==================== DIAGNOSTIC / AUXILIARY FIELDS ====================

volTensorField D_P(
    IOobject(
        "D_P",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE),
    mesh,
    dimensionedTensor("D_P", dimensionSet(0, 0, -1, 0, 0), tensor::zero));

volTensorField epsilon_dot(
    IOobject(
        "epsilon_dot",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE),
    mesh,
    dimensionedTensor("epsilon_dot", dimensionSet(0, 0, -1, 0, 0), tensor::zero));

volTensorField omega_dot(
    IOobject(
        "omega_dot",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE),
    mesh,
    dimensionedTensor("omega_dot", dimensionSet(0, 0, -1, 0, 0), tensor::zero));

volScalarField Acc_Strain(
    IOobject(
        "Acc_Strain",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE),
    mesh,
    dimensionedScalar("Acc_Strain", dimless, 0.0),
    zeroGradientFvPatchScalarField::typeName);

volVectorField divsigma(
    IOobject(
        "divsigma",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE),
    mesh,
    dimensionedVector("divsigma", dimensionSet(1, -2, -2, 0, 0), vector::zero));

volVectorField Euler_Angles(
    IOobject(
        "Euler_Angles",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE),
    mesh,
    dimensionedVector("Euler_Angles", dimless, vector::zero));

volScalarField GBREF(
    IOobject(
        "GBREF",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE),
    mesh,
    dimensionedScalar("GBREF", dimless, 0.0));

PtrList<volTensorField> Omega_P(N_Ori);
for (label count = 0; count < N_Ori; count++)
{
    Omega_P.set(
        count,
        new volTensorField(
            IOobject(
                "nOmega_P." + name(count),
                runTime.timeName(),
                mesh,
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE),
            mesh,
            dimensionedTensor("Omega_P", dimensionSet(0, 0, -1, 0, 0), tensor::zero)));
}

PtrList<volTensorField> Plastic_Strain(N_Ori);
for (label count = 0; count < N_Ori; count++)
{
    Plastic_Strain.set(
        count,
        new volTensorField(
            IOobject(
                "Plastic_Strain." + name(count),
                runTime.timeName(),
                mesh,
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE),
            mesh,
            dimensionedTensor("Plastic_Strain", dimensionSet(0, 0, -1, 0, 0), tensor::zero)));
}

PtrList<volTensorField> dsigmadt(N_Ori);
for (label count = 0; count < N_Ori; count++)
{
    dsigmadt.set(
        count,
        new volTensorField(
            IOobject(
                "dsigmadt." + Foam::name(count),
                runTime.timeName(),
                mesh,
                IOobject::READ_IF_PRESENT,
                IOobject::NO_WRITE),
            mesh,
            dimensionedTensor("dsigmadt", dimensionSet(1, -1, -3, 0, 0), tensor::zero)));
}

volTensorField testtraceeps(
    IOobject(
        "testtraceeps",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE),
    mesh,
    dimensionedTensor("testtraceeps", dimensionSet(1, -1, -3, 0, 0), tensor::zero));

volScalarField sumgamma(
    IOobject(
        "sumgamma",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE),
    mesh,
    dimensionedScalar("sumgamma", dimensionSet(0, 0, -1, 0, 0), 0.0),
    zeroGradientFvPatchScalarField::typeName);




#include "createPhi.H"