// Info << "Reading physicalProperties\n"
//      << endl;

// IOdictionary physicalProperties(
//     IOobject(
//         "physicalProperties",
//         runTime.constant(),
//         mesh,
//         IOobject::MUST_READ_IF_MODIFIED,
//         IOobject::NO_WRITE));

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

// int n_sub_steps(readLabel(PimpleDICT.lookup("nOuterCorrectors")));

scalar pimple_relax_U = PimpleDICT.lookupOrDefault<scalar>("relaxU", 0.5);

// scalar maxTimeStep(readScalar(controlDict.lookup("maxDeltaT")));

// Info<<"TREST: "<<n_sub_steps<<endl;

// exit(0);

dimensionedScalar pi("pi", dimensionSet(0, 0, 0, 0, 0), constant::mathematical::pi);



// dimensionedScalar Cv("Cv", dimensionSet(0, 2, -2, -1, 0), physicalProperties.lookup("Cv"));
// dimensionedScalar kappa("kappa", dimensionSet(1, 1, -3, -1, 0), physicalProperties.lookup("kappa"));
// dimensionedScalar rho_0("rho_0", dimensionSet(1, -3, 0, 0, 0), physicalProperties.lookup("rho_0"));

//  Info<<"TEST"<<(n_alpha[0] * m_alpha[0])<<endl;

// exit(0);

volTensorField D_P(
    IOobject(
        "D_P",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE),
    mesh,
    dimensionedTensor("D_P", dimensionSet(0, 0, -1, 0, 0), tensor::zero)//,
    // zeroGradientFvPatchTensorField::typeName
    );



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

volVectorField D
(
    IOobject
    (
        "D",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedVector("D",dimensionSet(0, 1, 0, 0, 0),vector::zero)
);





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

// volScalarField trsig
// (
//     IOobject
//     (
//         "trsig",
//         runTime.timeName(),
//         mesh,
//         IOobject::READ_IF_PRESENT,
//         IOobject::AUTO_WRITE
//     ),
//     mesh,
//     dimensionedScalar("trsig",dimensionSet(1, -1, -3, 0, 0),0.0)
// );

// volScalarField tredot
// (
//     IOobject
//     (
//         "tredot",
//         runTime.timeName(),
//         mesh,
//         IOobject::READ_IF_PRESENT,
//         IOobject::AUTO_WRITE
//     ),
//     mesh,
//     dimensionedScalar("tredot",dimensionSet(0, 0, -1, 0, 0),0.0)
// );

// volVectorField TESTADV
// (
//     IOobject
//     (
//         "TESTADV",
//         runTime.timeName(),
//         mesh,
//         IOobject::READ_IF_PRESENT,
//         IOobject::AUTO_WRITE
//     ),
//     mesh,
//     dimensionedVector("TESTADV",dimensionSet(1, -2, -2, 0, 0),vector::zero),
//                 zeroGradientFvPatchVectorField::typeName
// );

// volScalarField DP
// (
//     IOobject
//     (
//         "DP",
//         runTime.timeName(),
//         mesh,
//         IOobject::READ_IF_PRESENT,
//         IOobject::AUTO_WRITE
//     ),
//     mesh,
//     dimensionedScalar("DP",dimensionSet(1, -1, -3, 0, 0),0.0),
//                 zeroGradientFvPatchScalarField::typeName
// );

// volScalarField DP2
// (
//     IOobject
//     (
//         "DP2",
//         runTime.timeName(),
//         mesh,
//         IOobject::READ_IF_PRESENT,
//         IOobject::AUTO_WRITE
//     ),
//     mesh,
//     dimensionedScalar("DP2",dimensionSet(1, -1, -3, 0, 0),0.0),
//                 zeroGradientFvPatchScalarField::typeName
// );




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
    dimensionedScalar("rhoPhi", dimMass / dimTime, 0.0));

// volScalarField divU
// (
//     IOobject
//     (
//         "divU",
//         runTime.timeName(),
//         mesh,
//         IOobject::NO_READ,
//         IOobject::AUTO_WRITE
//     ),
//     mesh,
//     dimensionedScalar("divU",dimensionSet(0, 0, -1, 0, 0),0.0),
//                 zeroGradientFvPatchScalarField::typeName
// );

volScalarField Sum_Phases(
    IOobject(
        "Sum_Phases",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE),
    mesh,
    dimensionedScalar("Sum_Phases", dimensionSet(0, 0, 0, 0, 0), 0.0));



IOdictionary GrainProperties(
    IOobject(
        "GrainProperties",   // dictionary name
        runTime.constant(),  // dict is found in "constant"
        mesh,                // registry for the dict
        IOobject::MUST_READ, // must exist, otherwise failure
        IOobject::NO_WRITE   // dict is only read by the solver
        ));


dimensionedScalar interface_width("interface_width", dimensionSet(0, 1, 0, 0, 0), GrainProperties);

// const Switch PF_drivingForceActive(GrainProperties.lookup("PF_drivingForceActive"));

const Switch PF_drivingForceActive(GrainProperties.lookupOrDefault<Switch>("PF_drivingForceActive", false));

/*
lookupOrDefault<Switch>("Radial_Polar_HS", true)
*/

// const dimensionedScalar Mm
// (
//     "Mm",
//     dimensionSet(1, 0, 0, 0, -1),
//     transportProperties//.lookup("Mm")
// );

label N_Ori(readLabel(GrainProperties.lookup("N_Phases")));
Info << "Number of Phase Field Order Parameters: " << N_Ori << endl;

// Ready to read in multiple stiffness matrices

label N_Mats(readLabel(GrainProperties.lookup("N_Materials")));
Info << "Number of Materials: " << N_Mats << endl;

label n_substeps(readLabel(GrainProperties.lookup("n_constitutive_substeps")));
Info << "Number of Constitutive Sub-cycling steps: " << n_substeps << endl;


// Optional slip-rate bound |gammadot| <= maxSlipRate (1/s): the drag-limited regime the
// power law does not describe.  0 = none.  overstress_cap is obsolete and ignored.
const scalar maxSlipRate
(
    GrainProperties.lookupOrDefault<scalar>("maxSlipRate", 0.0)
);
if (maxSlipRate > 0.0)
{
    Info << "Slip-rate bound maxSlipRate: " << maxSlipRate << " 1/s" << endl;
}
else
{
    Info << "Slip-rate bound maxSlipRate: none" << endl;
}
if (GrainProperties.found("overstress_cap"))
{
    WarningInFunction
        << "overstress_cap is no longer used (it bounded |tau|/g, i.e. the slip rate at "
        << "gammadot_0*cap^k, which depends on k); use maxSlipRate (1/s) instead."
        << endl;
}



labelList Phase2Material(N_Ori, 0);

GrainProperties.lookup("Phase2Material") >> Phase2Material;

for (label c = 0; c < N_Ori; c++)
{
    Info << "Phase[" << c << "] = Material[" << Phase2Material[c] << "]" << endl;
}

Fastor::Tensor<double, 3, 3, 3, 3> Cprime2[N_Mats]; // array of tensor objects



std::vector<Material> MatList(N_Mats, Material());

for (label i = 0; i < N_Mats; i++)
{

    word Mat_dict("Material." + name(i));
    dictionary Mat = GrainProperties.subDict(Mat_dict);

    Fastor::Tensor<double, 3, 3, 3, 3> Cprime_temp_Read(0.0);

    for (int m = 0; m < 3; m++)
    {
        for (int n = 0; n < 3; n++)
        {
            for (int p = 0; p < 3; p++)
            {
                for (int q = 0; q < 3; q++)
                {
                    // Check for symmetry (minor and major) and warn the user ... 
                    word Cijkl("C" + name(m) + name(n) + name(p) + name(q));
                    word Cjikl("C" + name(n) + name(m) + name(p) + name(q));
                    word Cijlk("C" + name(m) + name(n) + name(q) + name(p));
                    word Cjilk("C" + name(n) + name(m) + name(q) + name(p));

                    double a = Mat.lookupOrDefault<scalar>(Cijkl, 0.0);
                    double b = Mat.lookupOrDefault<scalar>(Cjikl, 0.0);
                    double c = Mat.lookupOrDefault<scalar>(Cijlk, 0.0);
                    double d = Mat.lookupOrDefault<scalar>(Cjilk, 0.0);

                    double tol = 1e-12;
                    if (std::abs(a - b) > tol || std::abs(a - c) > tol || std::abs(a - d) > tol ||
                        std::abs(b - c) > tol || std::abs(b - d) > tol || std::abs(c - d) > tol)
                    {
                        WarningInFunction << "non-symmetric elastic stiffness detected:\n"
                                        << "      components:\n"
                                        << "         C_" << m << n << p << q << " = " << a << "\n"
                                        << "         C_" << n << m << p << q << " = " << b << "\n"
                                        << "         C_" << m << n << q << p << " = " << c << "\n"
                                        << "         C_" << n << m << q << p << " = " << d << "\n" << nl << endl;
                                        // << "      ...Taking the weighted average of the four..." << nl << endl;
                    }

                    Cprime_temp_Read(m, n, p, q) = a; 
                }
            }
        }
    }

    // Enforce minor and major symmetry (Mandel projection) before rotation and inversion
    Cprime2[i] = Foam::Material::castSymmetry(Cprime_temp_Read);

    // Strict symmetry check of the result
    Foam::Material::assertSymmetry(Cprime2[i]);

    // Inversion self-test, enabled by "verbose true;" in GrainProperties (Harry Smith)

    if (GrainProperties.lookupOrDefault<Switch>("verbose", false))//Harry Smith
    {
        Info << "\n=== invertRankFour self-test for Material." << i << " ===\n" << nl;
        auto S = Foam::Material::invertRankFour(Cprime2[i], true);
        bool isInverse = Foam::Material::isRankFourInverse(Cprime2[i], S, true);
        Info << "isRankFourInverse = " << isInverse << nl;
        Info << "...complete\n" << nl;
    }

    // // Quick check on the double contraction (twice) function 
    // Fastor::Tensor<double, 3,3,3,3> test_fourth_order(0.0); // Fourth order symmetry tensor 
    // Fastor::Tensor<double,3,3> test_second_order(0.0);      // Second order identity tensor 

    // for (int i = 0; i < 3; i++)
    // {
    //     for (int j = 0; j < 3; j++)
    //     {
    //         for (int k = 0; k < 3; k++)
    //         {
    //             for (int l = 0; l < 3; l++) 
    //             {
    //                 test_fourth_order(i,j,k,l) = 0.5 * ( (i==k && j==l) + (i==l && j==k) );
    //             }
    //         }
    //     }
    // }

    // for (int i = 0; i < 3; i++)
    // {
    //     test_second_order(i,i) = 1.0; 
    // }

    // double testContraction = Material::doubleContractionTwice(test_second_order, test_fourth_order, true); 
    
    // exit(0); 
// dimensionedScalar Cvtemp("Cvtemp", dimensionSet(0, 2, -2, -1, 0), Mat.lookup("Cv"));
// dimensionedScalar kappatemp("kappatemp", dimensionSet(1, 1, -3, -1, 0), Mat.lookup("kappa"));
// dimensionedScalar rho_0temp("rho_0temp", dimensionSet(1, -3, 0, 0, 0), Mat.lookup("rho_0"));

    scalar Cvtemp(readScalar(Mat.lookup("Cv")));
    scalar kappatemp(readScalar(Mat.lookup("kappa")));
    scalar rho_0temp(readScalar(Mat.lookup("rho_0")));

    // scalar rho
    MatList[i].rho_0 = rho_0temp;
    MatList[i].Cv = Cvtemp;
    MatList[i].kappa = kappatemp;


    label N_slips_Temp(readLabel(Mat.lookup("N_slips")));

    // Guard the fixed-size stack arrays of the constitutive kernel (MAX_SLIPS = 48)
    if (N_slips_Temp < 0 || N_slips_Temp > 48)
    {
        FatalErrorInFunction
            << Mat_dict << ": N_slips = " << N_slips_Temp
            << " outside the supported range 0..48 "
            << "(MAX_SLIPS stack arrays in StressUpdate.H)."
            << exit(FatalError);
    }

    MatList[i].N_Slips = N_slips_Temp;

    for (int j = 0; j < MatList[i].N_Slips; j++)
    {

        word Slip_direction_name = "m" + name(j);
        word Slip_normal_name = "n" + name(j);

        vector temp_m(Mat.lookup(Slip_direction_name));
        vector temp_n(Mat.lookup(Slip_normal_name));


        // Normalise the slip systems and check orthogonality
        const scalar magM = mag(temp_m);
        const scalar magN = mag(temp_n);
        if (magM < SMALL || magN < SMALL)
        {
            FatalErrorInFunction
                << Mat_dict << ": slip system " << j
                << " has a zero direction or normal vector."
                << exit(FatalError);
        }
        temp_m /= magM;
        temp_n /= magN;

        const scalar mDotN = temp_m & temp_n;
        if (mag(mDotN) > 1e-6)
        {
            FatalErrorInFunction
                << Mat_dict << ": slip system " << j
                << " is not orthogonal: normalised m.n = " << mDotN
                << nl
                << "For HCP pyramidal systems check that the c/a ratio "
                << "used to build m matches the one used for n."
                << exit(FatalError);
        }

        MatList[i].SP_norm.push_back(temp_n);
        MatList[i].SP_direction.push_back(temp_m);
    }

    scalar initialCRSS_temp(readScalar(Mat.lookup("initialCRSS")));
    MatList[i].initial_CRSS = initialCRSS_temp;

    scalar tau_S_temp(readScalar(Mat.lookup("tau_S")));
    MatList[i].tau_S = tau_S_temp;

    scalar h_0_temp(readScalar(Mat.lookup("h_0")));
    MatList[i].h_0 = h_0_temp;

    scalar r_temp(readScalar(Mat.lookup("r")));
    MatList[i].r = r_temp;

    scalar k_temp(readScalar(Mat.lookup("k_exponent")));
    MatList[i].k_exponent = k_temp;

    scalar gammadot_0_temp(readScalar(Mat.lookup("gammadot_0")));
    MatList[i].ref_shear_rate = gammadot_0_temp;
}

// Exit cleanly after the self-tests
if (GrainProperties.lookupOrDefault<Switch>("verbose", false))
{
    Info << "All material inversion self-tests complete. Exiting.\n" << endl;
    exit(0);
}













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
                dimensionedScalar("g_alphas", dimensionSet(0, 0, 0, 0, 0, 0, 0), MatList[Phase2Material[i]].initial_CRSS /*16.0e6*/), // should be initial for particular material
                zeroGradientFvPatchScalarField::typeName));
    }
}

// PtrList<volScalarField> g_alphas(N_Slips);

// for(label count = 0; count < N_Slips; count++)	//hardening field per slip system
// {

//     word name_g_alpha ("g_alpha." + name(count));

//     g_alphas.set
//     (
//         count,
//         new volScalarField
//         (
//             IOobject
//             (
//                 name_g_alpha,
//                 runTime.timeName(),
//                 mesh,
//                 IOobject::READ_IF_PRESENT,//IOobject::MUST_READ,
//                 IOobject::AUTO_WRITE
//             ),
//             mesh,
//             dimensionedScalar("g_alpha",dimensionSet(0, 0, 0, 0, 0),16.0e6),
//              zeroGradientFvPatchScalarField::typeName
//         )
//     );
// }

// PtrList<volScalarField> g_alpha_dots(N_Slips);//source terms on RHS of G_alpha_equation

// for(label count = 0; count < N_Slips; count++)	//hardening field per slip system
// {

//     word name_g_alpha_dot ("g_alpha_dot." + name(count));

//     g_alpha_dots.set
//     (
//         count,
//         new volScalarField
//         (
//             IOobject
//             (
//                 name_g_alpha_dot,
//                 runTime.timeName(),
//                 mesh,
//                 IOobject::READ_IF_PRESENT,//IOobject::MUST_READ,
//                 IOobject::NO_WRITE
//             ),
//             mesh,
//             dimensionedScalar("g_alpha_dot",dimensionSet(0, 0, -1, 0, 0),0.0),
//              zeroGradientFvPatchScalarField::typeName
//         )
//     );
// }

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
                dimensionedScalar("g_alpha_dots", dimensionSet(0, 0, -1, 0, 0, 0, 0), 0.0),
                zeroGradientFvPatchScalarField::typeName));
    }
}

// PtrList<volScalarField> gamma_alpha_dots(N_Slips);//source terms on RHS of G_alpha_equation

// for(label count = 0; count < N_Slips; count++)	//hardening field per slip system
// {

//     word name_gamma_alpha_dot ("gamma_alpha_dot." + name(count));

//     gamma_alpha_dots.set
//     (
//         count,
//         new volScalarField
//         (
//             IOobject
//             (
//                 name_gamma_alpha_dot,
//                 runTime.timeName(),
//                 mesh,
//                 IOobject::READ_IF_PRESENT,//IOobject::MUST_READ,
//                 IOobject::AUTO_WRITE
//             ),
//             mesh,
//             dimensionedScalar("gamma_alpha_dot",dimensionSet(0, 0, -1, 0, 0),0.0),
//              zeroGradientFvPatchScalarField::typeName
//         )
//     );
// }

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
                dimensionedScalar("gamma_alpha_dots", dimensionSet(0, 0, -1, 0, 0, 0, 0), 0.0),
                zeroGradientFvPatchScalarField::typeName));
    }
}

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
                dimensionedScalar("gamma_alphas", dimensionSet(0, 0, 0, 0, 0, 0, 0), 0.0),
                zeroGradientFvPatchScalarField::typeName));
    }
}

volScalarField sumgamma
(
    IOobject
    (
        "sumgamma",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("sumgamma",dimensionSet(0, 0, -1, 0, 0),0.0),
                zeroGradientFvPatchScalarField::typeName
);

PtrList<volScalarField> PopBal(N_Ori); // pointer to phase fields

for (label count = 0; count < N_Ori; count++) // how to assign the popbal array of fields
{

    word name_ni("n." + name(count));

    PopBal.set(
        count,
        new volScalarField(
            IOobject(
                name_ni,
                runTime.timeName(),
                mesh,
                IOobject::READ_IF_PRESENT, // IOobject::MUST_READ,
                IOobject::AUTO_WRITE),
            mesh,
            dimensionedScalar("ni", dimensionSet(0, 0, 0, 0, 0), 0.0),
            zeroGradientFvPatchScalarField::typeName));
}

PtrList<volScalarField> PopBal_av(N_Ori);

for(label count = 0; count < N_Ori; count++)
{

    word name_ni_av ("n_av" + name(count));

    PopBal_av.set
    (
        count,
        new volScalarField
        (
            IOobject
            (
                name_ni_av,
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedScalar("niav",dimensionSet(0, 0, 0, 0, 0),0.0),
            zeroGradientFvPatchScalarField::typeName
        )
    );
}











PtrList<volScalarField> Laplacian_n(N_Ori);

for(label count = 0; count < N_Ori; count++)
{

    word name_laplacian ("Laplacian" + name(count));

    Laplacian_n.set
    (
        count,
        new volScalarField
        (
            IOobject
            (
                name_laplacian,
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedScalar("Laplacian_n",dimensionSet(0, -2, 0, 0, 0),0.0),
            zeroGradientFvPatchScalarField::typeName
        )
    );
}



PtrList<volScalarField> dpdt(N_Ori);
for(label count = 0; count < N_Ori; count++)
{

    word name_dpdt ("dpdt" + name(count));

    dpdt.set
    (
        count,
        new volScalarField
        (
            IOobject
            (
                name_dpdt,
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedScalar("dpdt",dimensionSet(0, 0, -1, 0, 0),0.0),
            zeroGradientFvPatchScalarField::typeName
        )
    );
}

// Diagnostic fields: W_mech, W_interface, dpdt_crack
volScalarField W_mech
(
    IOobject
    (
        "W_mech",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("W_mech", dimensionSet(1, -1, -2, 0, 0), 0.0),
    zeroGradientFvPatchScalarField::typeName
);

volScalarField W_interface
(
    IOobject
    (
        "W_interface",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("W_interface", dimensionSet(1, -1, -2, 0, 0), 0.0),
    zeroGradientFvPatchScalarField::typeName
);

volScalarField dpdt_crack
(
    IOobject
    (
        "dpdt_crack",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("dpdt_crack", dimensionSet(0, 0, -1, 0, 0), 0.0),
    zeroGradientFvPatchScalarField::typeName
);









Sum_Phases *= 0.0;
for (label count = 0; count < N_Ori; count++)
{

    volScalarField &ni = PopBal[count];

    Sum_Phases += ni;
}
Sum_Phases.correctBoundaryConditions();

//  Project any out of bounds phases back onto the Gibbs Simplex
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
    if
    (
        !Sum_Phases.boundaryField()[patchI].coupled()
     && Sum_Phases.boundaryField()[patchI].size() > 0
    )
    {
        forAll(Sum_Phases.boundaryField()[patchI], faceI)
        {

            if (Sum_Phases.boundaryField()[patchI][faceI] > 0.0)
        {
            temp_holder.boundaryFieldRef()[patchI][faceI] *= (1.0 / Sum_Phases.boundaryField()[patchI][faceI]);
        }

        }

    }

}



    ni = temp_holder;
}


// momentumCoupling: explicit (default) | stabilised (Jasak-Weller implicit-stiffness
// split, UEqn.H) | quasiStatic (the split without the acoustic limit).  Legacy switches mapped.
word momentumCoupling
(
    GrainProperties.lookupOrDefault<word>("momentumCoupling", "explicit")
);
if (GrainProperties.lookupOrDefault<Switch>("quasiStatic", false))
{
    momentumCoupling = "quasiStatic";
}
else if
(
    GrainProperties.lookupOrDefault<Switch>("compactStabiliser", false)
 && momentumCoupling == "explicit"
)
{
    momentumCoupling = "stabilised";
}
if
(
    momentumCoupling != "explicit"
 && momentumCoupling != "stabilised"
 && momentumCoupling != "quasiStatic"
)
{
    FatalErrorInFunction
        << "Unknown momentumCoupling '" << momentumCoupling
        << "'. Allowed: explicit (default), stabilised, quasiStatic."
        << exit(FatalError);
}
Info<< "Momentum coupling: " << momentumCoupling << endl;

// quasiStatic: acoustic limit off; compactStabiliser: the split is assembled in UEqn.H
const Switch quasiStatic(momentumCoupling == "quasiStatic");
const Switch compactStabiliser(momentumCoupling != "explicit");

// Implicit stiffness = impKFactor * (phase-weighted longitudinal modulus bound);
// the outer iteration converges for impKFactor > 0.5 (0.6 suits cubic metals).
const scalar impKFactor
(
    GrainProperties.lookupOrDefault<scalar>("impKFactor", 0.6)
);

// Time-step controls applied when quasiStatic is true (a non-positive value disables one)
const scalar maxDeltaEps
(
    GrainProperties.lookupOrDefault<scalar>("maxDeltaEps", 1e-4)
);
const scalar maxViscoplasticNumber
(
    GrainProperties.lookupOrDefault<scalar>("maxViscoplasticNumber", 2.0)
);

// Adaptive sub-stepping: raise N per cell and phase so that sub_dt*lambda <= substepTarget
// (forward Euler needs < 2); n_constitutive_substeps is the minimum, maxSubsteps the cap.
const Switch adaptiveSubsteps
(
    GrainProperties.lookupOrDefault<Switch>("adaptiveSubsteps", quasiStatic)
);
const scalar substepTarget
(
    GrainProperties.lookupOrDefault<scalar>("substepTarget", 0.5)
);
const label maxSubsteps
(
    GrainProperties.lookupOrDefault<label>("maxSubsteps", 2000)
);

if (compactStabiliser && !quasiStatic)
{
    Info<< "momentumCoupling stabilised: implicit-stiffness split in the momentum "
        << "equation at the acoustic time step (impKFactor " << impKFactor << ")" << nl << endl;
}
if (compactStabiliser && impKFactor <= 0.5)
{
    FatalErrorInFunction
        << "impKFactor = " << impKFactor
        << " must exceed 0.5 for the implicit-stiffness outer iteration to converge."
        << exit(FatalError);
}
if (adaptiveSubsteps && (substepTarget <= 0.0 || maxSubsteps < n_substeps))
{
    FatalErrorInFunction
        << "adaptiveSubsteps requires substepTarget > 0 and "
        << "maxSubsteps >= n_constitutive_substeps."
        << exit(FatalError);
}

// Longitudinal modulus bound per material: largest acoustic-tensor eigenvalue over directions
scalarList GammaMaxMat(N_Mats, 0.0);

// Bulk modulus of C' per material (C'_aabb/9); in EOS mode UEqn.H raises the implicit
// stiffness by max(0, rho_0 c0^2 - KCprime).
scalarList KCprimeMat(N_Mats, 0.0);

// Schmid-projected moduli P_s : C' : P_s per material and slip system (crystal frame)
std::vector<std::vector<double>> PCPmat(N_Mats);

for (label i = 0; i < N_Mats; i++)
{
    const Fastor::Tensor<double,3,3,3,3>& Cm = Cprime2[i];

    const int nDirsQS = 800;
    scalar GamMaxQS = 0.0;
    for (int d = 0; d < nDirsQS; d++)
    {
        const scalar zf = 1.0 - 2.0*(d + 0.5)/nDirsQS;
        const scalar rf = Foam::sqrt(max(1.0 - zf*zf, 0.0));
        const scalar th = M_PI*(3.0 - Foam::sqrt(5.0))*d;
        const scalar nv[3] = {rf*Foam::cos(th), rf*Foam::sin(th), zf};

        double Gam[3][3] = {{0,0,0},{0,0,0},{0,0,0}};
        for (int a = 0; a < 3; a++)
            for (int b = 0; b < 3; b++)
                for (int jj = 0; jj < 3; jj++)
                    for (int ll = 0; ll < 3; ll++)
                        Gam[a][b] += Cm(a,jj,b,ll)*nv[jj]*nv[ll];

        double v[3] = {0.577, 0.577, 0.577};
        double lam = 0.0;
        for (int it = 0; it < 60; it++)
        {
            double w[3] =
            {
                Gam[0][0]*v[0] + Gam[0][1]*v[1] + Gam[0][2]*v[2],
                Gam[1][0]*v[0] + Gam[1][1]*v[1] + Gam[1][2]*v[2],
                Gam[2][0]*v[0] + Gam[2][1]*v[1] + Gam[2][2]*v[2]
            };
            lam = Foam::sqrt(w[0]*w[0] + w[1]*w[1] + w[2]*w[2]);
            if (lam < VSMALL) break;
            v[0] = w[0]/lam; v[1] = w[1]/lam; v[2] = w[2]/lam;
        }
        GamMaxQS = max(GamMaxQS, scalar(lam));
    }
    GammaMaxMat[i] = 1.02*GamMaxQS;

    for (int a = 0; a < 3; a++)
        for (int b = 0; b < 3; b++)
            KCprimeMat[i] += Cm(a,a,b,b)/9.0;

    PCPmat[i].assign(MatList[i].N_Slips, 0.0);
    double pcpMin = GREAT, pcpMax = 0.0;
    for (int s = 0; s < MatList[i].N_Slips; s++)
    {
        const vector& mS = MatList[i].SP_direction[s];
        const vector& nS = MatList[i].SP_norm[s];
        double P[3][3];
        for (int a = 0; a < 3; a++)
            for (int b = 0; b < 3; b++)
                P[a][b] = 0.5*(mS[a]*nS[b] + nS[a]*mS[b]);

        double pcp = 0.0;
        for (int a = 0; a < 3; a++)
            for (int b = 0; b < 3; b++)
                for (int c = 0; c < 3; c++)
                    for (int d = 0; d < 3; d++)
                        pcp += P[a][b]*Cm(a,b,c,d)*P[c][d];

        PCPmat[i][s] = pcp;
        pcpMin = min(pcpMin, pcp);
        pcpMax = max(pcpMax, pcp);
    }

    Info<< "Material." << i << ": longitudinal modulus bound "
        << GammaMaxMat[i]/1e9 << " GPa, bulk modulus of C' "
        << KCprimeMat[i]/1e9 << " GPa";
    if (MatList[i].N_Slips > 0)
    {
        Info<< ", Schmid-projected moduli " << pcpMin/1e9 << " - " << pcpMax/1e9 << " GPa";
    }
    Info<< endl;
}

// Implicit stiffness coefficient of the momentum equation (assembled in UEqn.H)
volScalarField impK
(
    IOobject
    (
        "impK",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("impK", dimensionSet(1, -1, -2, 0, 0), 0.0),
    zeroGradientFvPatchScalarField::typeName
);

// Viscoplastic relaxation rate per phase and cell (StressUpdate.H -> setDeltaT.H);
// allocated in quasi-static mode only.
List<scalarField> lamVp(N_Ori);
forAll(lamVp, h)
{
    lamVp[h].setSize(quasiStatic ? mesh.nCells() : 0, 0.0);
}

// Energy budget (energyBudget.H): per-phase plastic work rate, its mixture field and
// the running integrals of the domain energy balance.
List<scalarField> plasticWorkRatePhase(N_Ori);
forAll(plasticWorkRatePhase, h)
{
    plasticWorkRatePhase[h].setSize(mesh.nCells(), 0.0);
}

volScalarField plasticWorkRate
(
    IOobject
    (
        "plasticWorkRate",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("plasticWorkRate", dimensionSet(1, -1, -3, 0, 0), 0.0),
    zeroGradientFvPatchScalarField::typeName
);

const Switch energyBudgetActive
(
    GrainProperties.lookupOrDefault<Switch>("energyBudget", false)
);

// Boundary power terms from TEqn.H / EEqn.H on the final outer iteration
scalar ebBoundaryStressPower = 0.0;   // sum_b (sigma.U).S_f           [W]
scalar ebBoundaryConvective  = 0.0;   // sum_b rhoPhi_b (e_b + K_b)   [W]
scalar ebBoundaryHeat        = 0.0;   // sum_b kappa_b snGrad(T)|S_f| [W]
scalar ebBoundaryMassFlux    = 0.0;   // sum_b rhoPhi_b               [kg/s]

// Running integrals (J, kg) and reference values at the first step
scalar ebWork = 0.0, ebConv = 0.0, ebHeat = 0.0, ebDq = 0.0, ebDp = 0.0;
scalar ebMassIn = 0.0;
scalar ebE0 = 0.0, ebM0 = 0.0;
bool   ebInitialised = false;

if (quasiStatic)
{
    Info<< "momentumCoupling quasiStatic: implicit-stiffness coupling, acoustic limit off (impKFactor "
        << impKFactor << "); dt limited by maxDeltaEps " << maxDeltaEps
        << " and maxViscoplasticNumber " << maxViscoplasticNumber << nl << endl;
}
if (adaptiveSubsteps)
{
    Info<< "adaptiveSubsteps=true: constitutive sub-steps between " << n_substeps
        << " and " << maxSubsteps << " (target sub_dt*lambda = " << substepTarget << ")"
        << nl << endl;
}


volScalarField rhotemp(
    IOobject(
        "rhotemp",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE),
    mesh,
    dimensionedScalar("rhotemp", dimensionSet(1, -3, 0, 0, 0), 0.0),
    zeroGradientFvPatchScalarField::typeName
    );




volScalarField kappa
(
    IOobject
    (
        "kappa",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE//change to no write when working
    ),
    mesh,
    dimensionedScalar("kappa",dimensionSet(1, 1, -3, -1, 0),0.0)
);


volScalarField Cv
(
    IOobject
    (
        "Cv",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE//change to no write when working
    ),
    mesh,
    dimensionedScalar("Cv", dimensionSet(0, 2, -2, -1, 0), 0.0)
);



        

   for (label h = 0; h < N_Ori; h++)
    {

dimensionedScalar dimkappa("dimkappa", dimensionSet(1, 1, -3, -1, 0), 1.0);

dimensionedScalar dimcv("dimcv", dimensionSet(0, 2, -2, -1, 0), 1.0);

dimensionedScalar dimrho("dimrho", dimensionSet(1, -3, 0, 0, 0), 1.0);


kappa+=PopBal[h]*MatList[Phase2Material[h]].kappa*dimkappa;

Cv+=PopBal[h]*MatList[Phase2Material[h]].Cv*dimcv;

rhotemp+=PopBal[h]*MatList[Phase2Material[h]].rho_0*dimrho;
    }

Info << "Reading field rho\n"
     << endl;
volScalarField rho(
    IOobject(
        "rho",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE),
        rhotemp);
    
// volScalarField rho(
//     IOobject(
//         "rho",
//         runTime.timeName(),
//         mesh,
//         IOobject::READ_IF_PRESENT,
//         IOobject::AUTO_WRITE),
//     rhotemp
//     );

// rho.write();



volScalarField rhoCv
(
    IOobject("rhoCv", runTime.timeName(), mesh,
             IOobject::NO_READ, IOobject::NO_WRITE),
    rho*Cv
);

rhoCv.correctBoundaryConditions();


volScalarField GBREF
(
    IOobject
    (
        "GBREF",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("GBREF",dimensionSet(0,0,0,0,0),0.0)
);



PtrList<volTensorField> Omega_P(N_Ori); // pointer to phase fields

for (label count = 0; count < N_Ori; count++) // how to assign the popbal array of fields
{

    word name_Omega_P("nOmega_P." + name(count));

    Omega_P.set(
        count,
        new volTensorField(
            IOobject(
                name_Omega_P,
                runTime.timeName(),
                mesh,
                IOobject::READ_IF_PRESENT, // IOobject::MUST_READ,
                IOobject::AUTO_WRITE),
            mesh,
            dimensionedTensor("Omega_P", dimensionSet(0, 0, -1, 0, 0), tensor::zero)//,
            // zeroGradientFvPatchTensorField::typeName
            )
            );
}

PtrList<volTensorField> Plastic_Strain(N_Ori); // pointer to phase fields

for (label count = 0; count < N_Ori; count++) // how to assign the popbal array of fields
{

    word name_Plastic_Strain("Plastic_Strain." + name(count));

    Plastic_Strain.set(
        count,
        new volTensorField(
            IOobject(
                name_Plastic_Strain,
                runTime.timeName(),
                mesh,
                IOobject::READ_IF_PRESENT, // IOobject::MUST_READ,
                IOobject::AUTO_WRITE),
            mesh,
            dimensionedTensor("Plastic_Strain", dimensionSet(0, 0, -1, 0, 0), tensor::zero)//,
            // zeroGradientFvPatchTensorField::typeName
            )
            );
}

Random randObj2(clock::getTime());

// double phi1[N_Ori] = {0.0};

// double phi2[N_Ori] = {0.0};

// double phi3[N_Ori] = {0.0};

List<double> phi1 = GrainProperties.lookup("phi1");
List<double> phi2 = GrainProperties.lookup("phi2");
List<double> phi3 = GrainProperties.lookup("phi3");


PtrList<volScalarField> Rot1(N_Ori); // pointer to phase fields

for (label count = 0; count < N_Ori; count++) // how to assign the popbal array of fields
{

    word name_ROT("ROT.1." + name(count));

    Rot1.set(
        count,
        new volScalarField(
            IOobject(
                name_ROT,
                runTime.timeName(),
                mesh,
                IOobject::READ_IF_PRESENT, // IOobject::MUST_READ,
                IOobject::AUTO_WRITE),
            mesh,
            dimensionedScalar("ROT.1", dimensionSet(0, 0, 0, 0, 0), phi1[count]),
            zeroGradientFvPatchScalarField::typeName));
}

PtrList<volScalarField> Rot2(N_Ori); // pointer to phase fields

for (label count = 0; count < N_Ori; count++) // how to assign the popbal array of fields
{

    word name_ROT("ROT.2." + name(count));

    Rot2.set(
        count,
        new volScalarField(
            IOobject(
                name_ROT,
                runTime.timeName(),
                mesh,
                IOobject::READ_IF_PRESENT, // IOobject::MUST_READ,
                IOobject::AUTO_WRITE),
            mesh,
            dimensionedScalar("ROT.2", dimensionSet(0, 0, 0, 0, 0), phi2[count]),
            zeroGradientFvPatchScalarField::typeName));
}

PtrList<volScalarField> Rot3(N_Ori); // pointer to phase fields

for (label count = 0; count < N_Ori; count++) // how to assign the popbal array of fields
{

    word name_ROT("ROT.3." + name(count));

    Rot3.set(
        count,
        new volScalarField(
            IOobject(
                name_ROT,
                runTime.timeName(),
                mesh,
                IOobject::READ_IF_PRESENT, // IOobject::MUST_READ,
                IOobject::AUTO_WRITE),
            mesh,
            dimensionedScalar("ROT.3", dimensionSet(0, 0, 0, 0, 0), phi3[count]),
            zeroGradientFvPatchScalarField::typeName));
}

// Material TEST();

// Info<<MatList[0].RK4()<<endl;

// exit(0);


// PtrList<volScalarField> Total_slip_phi(N_Ori);

// for(label count = 0; count < N_Ori; count++)
// {

//     word name_Total_slip_phi ("Total_slip_phi." + name(count));

//     Total_slip_phi.set
//     (
//         count,
//         new volScalarField
//         (
//             IOobject
//             (
//                 name_Total_slip_phi,
//                 runTime.timeName(),
//                 mesh,
//                 IOobject::READ_IF_PRESENT,
//                 IOobject::AUTO_WRITE
//             ),
//             mesh,
//             dimensionedScalar("Total_slip_phi",dimensionSet(0, 0, 0, 0, 0),0.0)//,
//             // zeroGradientFvPatchScalarField::typeName
//         )
//     );
// }


volTensorField testtraceeps
(
    IOobject
    (
        "testtraceeps",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
        mesh,
dimensionedTensor("testtraceeps",dimensionSet(1, -1, -3, 0, 0),tensor::zero)//,
            //    zeroGradientFvPatchTensorField::typeName
);

volTensorField omega_dot(
    IOobject(
        "omega_dot",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE),
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
    // zeroGradientFvPatchTensorField::typeName
    );

// volTensorField omega_tot
// (
//     IOobject
//     (
//         "omega_tot",
//         runTime.timeName(),
//         mesh,
//         IOobject::READ_IF_PRESENT,
//         IOobject::AUTO_WRITE
//     ),
//         mesh,
// dimensionedTensor("omega_tot",dimensionSet(0, 0, 0, 0, 0),tensor::zero),
//                zeroGradientFvPatchTensorField::typeName
// );

volScalarField Acc_Strain(
    IOobject(
        "Acc_Strain",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE),
    mesh,
    dimensionedScalar("Acc_Strain", dimensionSet(0, 0, 0, 0, 0), 0.0),
    zeroGradientFvPatchScalarField::typeName
    );

volTensorField sigmaGlobal(
    IOobject(
        "sigmaGlobal",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE),
    mesh,
    dimensionedTensor("sigmaGlobal", dimensionSet(1, -1, -2, 0, 0), tensor::zero)//,
    // zeroGradientFvPatchTensorField::typeName
    );

    volTensorField sigmaGlobaltemp(
    IOobject(
        "sigmaGlobaltemp",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE),
    mesh,
    dimensionedTensor("sigmaGlobaltemp", dimensionSet(1, -1, -2, 0, 0), tensor::zero)//,
    // zeroGradientFvPatchTensorField::typeName
    );

// volTensorField sigma(
//     IOobject(
//         "sigma",
//         runTime.timeName(),
//         mesh,
//         IOobject::READ_IF_PRESENT,
//         IOobject::AUTO_WRITE),
//     mesh,
//     dimensionedTensor("sigma", dimensionSet(1, -1, -2, 0, 0), tensor::zero)//,
//     // zeroGradientFvPatchTensorField::typeName
//     );

// volTensorField sigmadiff(
//     IOobject(
//         "sigmadiff",
//         runTime.timeName(),
//         mesh,
//         IOobject::READ_IF_PRESENT,
//         IOobject::AUTO_WRITE),
//     mesh,
//     dimensionedTensor("sisigmadiffgmadiff", dimensionSet(1, -1, -2, 0, 0), tensor::zero)//,
//     // zeroGradientFvPatchTensorField::typeName
//     );

    // volTensorField sigmatemp(
    // IOobject(
    //     "sigmatemp",
    //     runTime.timeName(),
    //     mesh,
    //     IOobject::READ_IF_PRESENT,
    //     IOobject::AUTO_WRITE),
    // mesh,
    // dimensionedTensor("sigmatemp", dimensionSet(1, -1, -2, 0, 0), tensor::zero)//,
    // // zeroGradientFvPatchTensorField::typeName
    // );



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

volVectorField Euler_Angles
(
    IOobject
    (
        "Euler_Angles",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
        mesh,
dimensionedVector("Euler_Angles",dimensionSet(0, 0, 0, 0, 0),vector::zero)//,
                //zeroGradientFvPatchTensorField::typeName
);

// volTensorField testSIG
// (
//     IOobject
//     (
//         "testSIG",
//         runTime.timeName(),
//         mesh,
//         IOobject::READ_IF_PRESENT,
//         IOobject::AUTO_WRITE
//     ),
//         mesh,
// dimensionedTensor("testSIG",dimensionSet(1, -1, -2, 0, 0),tensor::zero)//,
//                 //zeroGradientFvPatchTensorField::typeName
// );

// volTensorField dsigmadt(
//     IOobject(
//         "dsigmadt",
//         runTime.timeName(),
//         mesh,
//         IOobject::NO_READ,
//         IOobject::AUTO_WRITE),
//     mesh,
//     dimensionedTensor("dsigmadt", dimensionSet(1, -1, -3, 0, 0), tensor::zero)//,
//     // zeroGradientFvPatchTensorField::typeName
//     );





PtrList<volTensorField> dsigmadt(N_Ori);

for (label count = 0; count < N_Ori; count++)
{
    word name_ni("dsigmadt." + Foam::name(count));
    
    dsigmadt.set(
        count,
        new volTensorField(
            IOobject(
                name_ni,
                runTime.timeName(),
                mesh,
                IOobject::READ_IF_PRESENT,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedTensor(
                "dsigmadt", 
                dimensionSet(1, -1, -3, 0, 0), 
                tensor::zero  // Note: tensor::zero, not Tensor::zero
            )//,
            // zeroGradientFvPatchTensorField::typeName
        )
    );
}


PtrList<volTensorField> sigtemp(N_Ori);

for (label count = 0; count < N_Ori; count++)
{
    word name_ni("sigtemp." + Foam::name(count));
    
    sigtemp.set(
        count,
        new volTensorField(
            IOobject(
                name_ni,
                runTime.timeName(),
                mesh,
                IOobject::READ_IF_PRESENT,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedTensor(
                "sigtemp", 
                dimensionSet(1, -1, -2, 0, 0), 
                tensor::zero  // Note: tensor::zero, not Tensor::zero
            )//,
            // zeroGradientFvPatchTensorField::typeName
        )
    );
}


PtrList<volTensorField> sigma(N_Ori);

for (label count = 0; count < N_Ori; count++)
{
    word name_ni("sigma." + Foam::name(count));
    
    sigma.set(
        count,
        new volTensorField(
            IOobject(
                name_ni,
                runTime.timeName(),
                mesh,
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            mesh,
            dimensionedTensor(
                "sigma", 
                dimensionSet(1, -1, -2, 0, 0), 
                tensor::zero  // Note: tensor::zero, not Tensor::zero
            ),
             zeroGradientFvPatchTensorField::typeName
        )
    );
}

#include "createPhi.H"






double **PF_Energy0_array;
PF_Energy0_array = new double *[N_Mats];
for(int i = 0; i <N_Mats; i++)
{
    PF_Energy0_array[i]=new double[N_Mats];
}

double **PF_Mobility0_array;
PF_Mobility0_array = new double *[N_Mats];
for(int i = 0; i <N_Mats; i++)
{
    PF_Mobility0_array[i]=new double[N_Mats];
}

Interface dummy("isotropic",T);
dummy.Assign_Energies(PF_Energy0_array);
dummy.Assign_Mobilities(PF_Mobility0_array);








// Rotated stiffness tensor per phase and cell
std::vector<std::vector<Fastor::Tensor<double,3,3,3,3>>> C_rotated(
    N_Ori,
    std::vector<Fastor::Tensor<double,3,3,3,3>>(mesh.nCells())
);



// PtrList<volScalarField> PopBal(N_Ori); // pointer to phase fields

// for (label count = 0; count < N_Ori; count++) // how to assign the popbal array of fields
// {

//     word name_ni("n." + name(count));

//     PopBal.set(
//         count,
//         new volScalarField(
//             IOobject(
//                 name_ni,
//                 runTime.timeName(),
//                 mesh,
//                 IOobject::READ_IF_PRESENT, // IOobject::MUST_READ,
//                 IOobject::AUTO_WRITE),
//             mesh,
//             dimensionedScalar("ni", dimensionSet(0, 0, 0, 0, 0), 0.0),
//             zeroGradientFvPatchScalarField::typeName));
// }







// volTensorField test_rot_diff
// (
//     IOobject
//     (
//         "test_rot_diff",
//         runTime.timeName(),
//         mesh,
//         IOobject::READ_IF_PRESENT,
//         IOobject::AUTO_WRITE
//     ),
//         mesh,
// dimensionedTensor("test_rot_diff",dimensionSet(0, 0, 0, 0, 0),(tensor::zero))//,
//                 //zeroGradientFvPatchTensorField::typeName
// );

// PtrList<volTensorField> test_rot_diff(N_Ori); // pointer to phase fields

// for (label count = 0; count < N_Ori; count++) // how to assign the popbal array of fields
// {

//     word name_ni("testRot." + name(count));

//     test_rot_diff.set(
//         count,
//         new volTensorField(
//             IOobject(
//                 name_ni,
//                 runTime.timeName(),
//                 mesh,
//                 IOobject::READ_IF_PRESENT, // IOobject::MUST_READ,
//                 IOobject::AUTO_WRITE),
//             mesh,
//             dimensionedTensor("ni", dimensionSet(0, 0, 0, 0, 0), (tensor::zero)),
//             zeroGradientFvPatchTensorField::typeName));
// }

const bool stressUpdate =
    GrainProperties.lookupOrDefault<Switch>("stressUpdate", true);

if (!stressUpdate)
{
    Info << "stressUpdate=false: prescribed-stress PF mode\n" << endl;
}

// Per-phase constant prescribed stress (with stressUpdate false only)
const bool perPhaseSigma =
    GrainProperties.lookupOrDefault<Switch>("PerPhaseSigma", false);

List<tensor> sigmaPrescribed;
if (perPhaseSigma)
{
    GrainProperties.lookup("sigmaPrescribed") >> sigmaPrescribed;

    if (sigmaPrescribed.size() != N_Ori)
    {
        FatalErrorInFunction
            << "sigmaPrescribed has " << sigmaPrescribed.size()
            << " entries but N_Phases = " << N_Ori
            << " (PerPhaseSigma requires one tensor per phase)."
            << exit(FatalError);
    }

    Info << "PerPhaseSigma=true: per-phase constant prescribed stress\n" << endl;
}

// Diagnostic: phase whose raw rate of change dpdt_crack (PF.H) reports; -1 = none
const label trackedPhase =
    GrainProperties.lookupOrDefault<label>("TrackedPhase", -1);

if (trackedPhase >= 0)
{
    if (trackedPhase >= N_Ori)
    {
        FatalErrorInFunction
            << "TrackedPhase = " << trackedPhase
            << " is out of range for N_Phases = " << N_Ori
            << exit(FatalError);
    }
    Info << "TrackedPhase=" << trackedPhase
         << ": dpdt_crack reports this phase's rate of change\n" << endl;
}

//   forAll(mesh.C(), cell)
//   {

//     for (label h = 0; h < N_Ori; h++)
//     {

//       double phi_1_total = phi1[h]; // ROT_1[cell];// phi1[h] ;//+ omega_tot[cell].component(5);

//       double phi_2_total = phi2[h]; // ROT_2[cell];// phi2[h] ;//+ omega_tot[cell].component(2);

//       double phi_3_total = phi3[h];

//       Fastor::Tensor<double, 3, 3> Rot_matrix(0.0);

//       Rot_matrix(0, 0) = Foam::cos(phi_2_total) * Foam::cos(phi_3_total);
//       Rot_matrix(0, 1) = (Foam::sin(phi_1_total) * Foam::sin(phi_2_total) * Foam::cos(phi_3_total)) - (Foam::cos(phi_1_total) * Foam::sin(phi_3_total));
//       Rot_matrix(0, 2) = (Foam::cos(phi_1_total) * Foam::sin(phi_2_total) * Foam::cos(phi_3_total)) + (Foam::sin(phi_1_total) * Foam::sin(phi_3_total));

//       Rot_matrix(1, 0) = Foam::cos(phi_2_total) * Foam::sin(phi_3_total);
//       Rot_matrix(1, 1) = (Foam::sin(phi_1_total) * Foam::sin(phi_2_total) * Foam::sin(phi_3_total)) + (Foam::cos(phi_1_total) * Foam::cos(phi_3_total));
//       Rot_matrix(1, 2) = (Foam::cos(phi_1_total) * Foam::sin(phi_2_total) * Foam::sin(phi_3_total)) - (Foam::sin(phi_1_total) * Foam::cos(phi_3_total));

//       Rot_matrix(2, 0) = -Foam::sin(phi_2_total);
//       Rot_matrix(2, 1) = Foam::sin(phi_1_total) * Foam::cos(phi_2_total);
//       Rot_matrix(2, 2) = Foam::cos(phi_1_total) * Foam::cos(phi_2_total);


//     // if(cell==240){
//     //     Info<<h<<endl;
//     // Fastor::print("Rot: ",Rot_matrix);
//     // }


//       test_rot_diff[h][cell].component(0)=Rot_matrix(0, 0);
//       test_rot_diff[h][cell].component(1)=Rot_matrix(0, 1);
//       test_rot_diff[h][cell].component(2)=Rot_matrix(0, 2);
//       test_rot_diff[h][cell].component(3)=Rot_matrix(1, 0);
//       test_rot_diff[h][cell].component(4)=Rot_matrix(1, 1);
//       test_rot_diff[h][cell].component(5)=Rot_matrix(1, 2);
//       test_rot_diff[h][cell].component(6)=Rot_matrix(2, 0);
//       test_rot_diff[h][cell].component(7)=Rot_matrix(2, 1);
//       test_rot_diff[h][cell].component(8)=Rot_matrix(2, 2);

//     }

//   }



// Fastor::Tensor<double, 3> TEST1(0.0);
// Fastor::Tensor<double, 3> TEST2(0.0);

// TEST1(0) = 2.0;
// TEST1(1) = 4.0;
// TEST1(2) = 5.0;

// TEST2(0) = 3.0;
// TEST2(1) = 7.0;
// TEST2(2) = 8.0;


// Fastor::print(TEST1);
// Fastor::print(TEST2);
// Fastor::print(outer(TEST1,TEST2));


// exit(0);


// forAll(mesh.C(), cell)
//   {
//     // test_rot_diff[cell].component(0)=1.0;
//     // test_rot_diff[cell].component(4)=2.0;
//     // test_rot_diff[cell].component(8)=3.0;

//     test_rot_diff[0][cell].component(0)=1.0/Foam::sqrt(2.0);
//     test_rot_diff[0][cell].component(1)=1.0/Foam::sqrt(6.0);
//     test_rot_diff[0][cell].component(2)=1.0/Foam::sqrt(3.0);
//     test_rot_diff[0][cell].component(3)=-1.0/Foam::sqrt(2.0);
//     test_rot_diff[0][cell].component(4)=1.0/Foam::sqrt(6.0);
//     test_rot_diff[0][cell].component(5)=1.0/Foam::sqrt(3.0);
//     test_rot_diff[0][cell].component(6)=0.0;
//     test_rot_diff[0][cell].component(7)=-2.0/Foam::sqrt(6.0);
//     test_rot_diff[0][cell].component(8)=1.0/Foam::sqrt(3.0);
//   }



// Fastor::Tensor<double, 3, 3> TEST_1(0.0);

// TEST_1(0,0) = -1.0;
// TEST_1(0,1) = 2.0;
// TEST_1(0,2) = 3.0;
// TEST_1(1,0) = 4.0;
// TEST_1(1,1) = -5.0;
// TEST_1(1,2) = 6.0;
// TEST_1(2,0) = 7.0;
// TEST_1(2,1) = 8.0;
// TEST_1(2,2) = 9.0;

// Fastor::Tensor<double, 3, 3> TEST_2(0.0);

// TEST_2(0,0) = 2.0;
// TEST_2(0,1) = -3.0;
// TEST_2(0,2) = -4.0;
// TEST_2(1,0) = 5.0;
// TEST_2(1,1) = -6.0;
// TEST_2(1,2) = 7.0;
// TEST_2(2,0) = 8.0;
// TEST_2(2,1) = -9.0;
// TEST_2(2,2) = 10.0;

// Fastor::print(Fastor::matmul(TEST_1,TEST_2));



// tensor T1(-1.0,2.0,3.0,4.0,-5.0,6.0,7.0,8.0,9.0);
// tensor T2(2.0,-3.0,-4.0,5.0,-6.0,7.0,8.0,-9.0,10.0);

// Info<<"Foam dot "<<(T1 & T2);

// exit(0);


// tensor T3(1.0,0.0,0.0,0.0,0.9397,-0.3420,0.0,0.3420,0.9397);
// tensor T4(1.0,0.0,0.0,0.0,2.0,0.0,0.0,0.0,3.0);

// Info<<"Foam dot "<<((T3 & T4) & T3.T());

// exit(0);


// //// TRYING TO MAKE SENSE OF ANGLES//
    // double phi_1_total1=phi1[1];//0.0314159;//0.5*M_PI;
    // double phi_2_total1=phi2[1];
    // double phi_3_total1=phi3[1];

    //   Fastor::Tensor<double, 3, 3> Rot_matrixtest1(0.0);

    //   Rot_matrixtest1(0, 0) = Foam::cos(phi_2_total1) * Foam::cos(phi_3_total1);
    //   Rot_matrixtest1(0, 1) = (Foam::sin(phi_1_total1) * Foam::sin(phi_2_total1) * Foam::cos(phi_3_total1)) - (Foam::cos(phi_1_total1) * Foam::sin(phi_3_total1));
    //   Rot_matrixtest1(0, 2) = (Foam::cos(phi_1_total1) * Foam::sin(phi_2_total1) * Foam::cos(phi_3_total1)) + (Foam::sin(phi_1_total1) * Foam::sin(phi_3_total1));

    //   Rot_matrixtest1(1, 0) = Foam::cos(phi_2_total1) * Foam::sin(phi_3_total1);
    //   Rot_matrixtest1(1, 1) = (Foam::sin(phi_1_total1) * Foam::sin(phi_2_total1) * Foam::sin(phi_3_total1)) + (Foam::cos(phi_1_total1) * Foam::cos(phi_3_total1));
    //   Rot_matrixtest1(1, 2) = (Foam::cos(phi_1_total1) * Foam::sin(phi_2_total1) * Foam::sin(phi_3_total1)) - (Foam::sin(phi_1_total1) * Foam::cos(phi_3_total1));

    //   Rot_matrixtest1(2, 0) = -Foam::sin(phi_2_total1);
    //   Rot_matrixtest1(2, 1) = Foam::sin(phi_1_total1) * Foam::cos(phi_2_total1);
    //   Rot_matrixtest1(2, 2) = Foam::cos(phi_1_total1) * Foam::cos(phi_2_total1);



    // double phi_1_total2=phi1[1]-0.7*M_PI;//0.5*M_PI;//0.0314159;//0.5*M_PI;//after 1e-7
    // double phi_2_total2=phi2[1];
    // double phi_3_total2=phi3[1];

    //   Fastor::Tensor<double, 3, 3> Rot_matrixtest(0.0);

    //   Rot_matrixtest(0, 0) = Foam::cos(phi_2_total2) * Foam::cos(phi_3_total2);
    //   Rot_matrixtest(0, 1) = (Foam::sin(phi_1_total2) * Foam::sin(phi_2_total2) * Foam::cos(phi_3_total2)) - (Foam::cos(phi_1_total2) * Foam::sin(phi_3_total2));
    //   Rot_matrixtest(0, 2) = (Foam::cos(phi_1_total2) * Foam::sin(phi_2_total2) * Foam::cos(phi_3_total2)) + (Foam::sin(phi_1_total2) * Foam::sin(phi_3_total2));

    //   Rot_matrixtest(1, 0) = Foam::cos(phi_2_total2) * Foam::sin(phi_3_total2);
    //   Rot_matrixtest(1, 1) = (Foam::sin(phi_1_total2) * Foam::sin(phi_2_total2) * Foam::sin(phi_3_total2)) + (Foam::cos(phi_1_total2) * Foam::cos(phi_3_total2));
    //   Rot_matrixtest(1, 2) = (Foam::cos(phi_1_total2) * Foam::sin(phi_2_total2) * Foam::sin(phi_3_total2)) - (Foam::sin(phi_1_total2) * Foam::cos(phi_3_total2));

    //   Rot_matrixtest(2, 0) = -Foam::sin(phi_2_total2);
    //   Rot_matrixtest(2, 1) = Foam::sin(phi_1_total2) * Foam::cos(phi_2_total2);
    //   Rot_matrixtest(2, 2) = Foam::cos(phi_1_total2) * Foam::cos(phi_2_total2);

    // Fastor::print("rot matrix initial: ",Rot_matrixtest1);
    // Fastor::print("rot matrix aftre - pi/2: ",Rot_matrixtest);

    //   exit(0);