//  Per-material entries required wih EOSMode MieGruneisen:
//
//      c0        3940;     // bulk sound speed at reference state   [m/s]
//      s_UsUp    1.489;    // linear Us-up Hugoniot slope           [-]
//      Gamma0    1.99;     // Gruneisen parameter at reference      [-]
//
//  Optional global entries (GrainProperties):
//
//      EOS_Tref     300.0;   // reference temperature: e = 0 at
//                            // (rho = rho0, T = EOS_Tref, p = 0)   [K]
//      EOS_qQuad    1.2;     // quadratic artificial-viscosity coeff [-]
//      EOS_qLinear  0.06;    // linear artificial-viscosity coeff    [-]
//      EOS_waveStabiliser 0.0;  // compact (Laplacian) momentum
//                            // stabiliser coefficient C_w [-]. 

const word EOSMode
(
    GrainProperties.lookupOrDefault<word>("EOSMode", "linearElastic")
);

const bool useEOS = (EOSMode == "MieGruneisen");

if (!useEOS && EOSMode != "linearElastic")
{
    FatalErrorInFunction
        << "Unknown EOSMode '" << EOSMode << "'. "
        << "Allowed: linearElastic (default), MieGruneisen."
        << exit(FatalError);
}

if (useEOS && !stressUpdate)
{
    FatalErrorInFunction
        << "EOSMode MieGruneisen requires stressUpdate = true "
        << "(the EOS closes the momentum/energy system; it is "
        << "meaningless in prescribed-stress PF mode)."
        << exit(FatalError);
}

Info<< "EOS mode: " << EOSMode << nl << endl;


if (useEOS)
{
    for (label i = 0; i < N_Mats; i++)
    {
        word Mat_dict_EOS("Material." + name(i));
        dictionary MatEOS = GrainProperties.subDict(Mat_dict_EOS);

        Material& matI = MatList[i];

        matI.c0     = readScalar(MatEOS.lookup("c0"));
        matI.s_UsUp = readScalar(MatEOS.lookup("s_UsUp"));
        matI.Gamma0 = readScalar(MatEOS.lookup("Gamma0"));

        if (matI.c0 <= 0.0 || matI.rho_0 <= 0.0)
        {
            FatalErrorInFunction
                << Mat_dict_EOS << ": c0 and rho_0 must be positive "
                << "when EOSMode is MieGruneisen."
                << exit(FatalError);
        }
        if (matI.s_UsUp < 1.0)
        {
            WarningInFunction
                << Mat_dict_EOS << ": s_UsUp = " << matI.s_UsUp
                << " < 1 is unusual for solids; proceeding." << endl;
        }

        // Reference longitudinal wave-speed bound for time-step control:

        {
            const Fastor::Tensor<double,3,3,3,3>& CmEOS = Cprime2[i];
            const int nDirsEOS = 800;
            scalar GamMaxEOS = 0.0;

            for (int d = 0; d < nDirsEOS; d++)
            {
                const scalar zf = 1.0 - 2.0*(d + 0.5)/nDirsEOS;
                const scalar rf = Foam::sqrt(max(1.0 - zf*zf, 0.0));
                const scalar th = M_PI*(3.0 - Foam::sqrt(5.0))*d;
                const scalar nv[3] =
                    {rf*Foam::cos(th), rf*Foam::sin(th), zf};

                double Gam[3][3] = {{0,0,0},{0,0,0},{0,0,0}};
                for (int a = 0; a < 3; a++)
                    for (int b = 0; b < 3; b++)
                        for (int jj = 0; jj < 3; jj++)
                            for (int ll = 0; ll < 3; ll++)
                                Gam[a][b] +=
                                    CmEOS(a,jj,b,ll)*nv[jj]*nv[ll];

                double v[3] = {0.577, 0.577, 0.577};
                double lam = 0.0;
                for (int it = 0; it < 60; it++)
                {
                    double w[3] =
                    {
                        Gam[0][0]*v[0]+Gam[0][1]*v[1]+Gam[0][2]*v[2],
                        Gam[1][0]*v[0]+Gam[1][1]*v[1]+Gam[1][2]*v[2],
                        Gam[2][0]*v[0]+Gam[2][1]*v[1]+Gam[2][2]*v[2]
                    };
                    lam = Foam::sqrt(w[0]*w[0]+w[1]*w[1]+w[2]*w[2]);
                    if (lam < VSMALL) break;
                    v[0] = w[0]/lam; v[1] = w[1]/lam; v[2] = w[2]/lam;
                }
                GamMaxEOS = max(GamMaxEOS, scalar(lam));
            }

            matI.cL0_bound = max
            (
                1.02*Foam::sqrt(max(GamMaxEOS, SMALL)/matI.rho_0),
                matI.c0
            );
        }

        Info<< "EOS " << Mat_dict_EOS
            << ": rho0 = "   << matI.rho_0
            << ", c0 = "     << matI.c0
            << ", s = "      << matI.s_UsUp
            << ", Gamma0 = " << matI.Gamma0
            << ", cL0(CFL bound) = " << matI.cL0_bound
            << ", K = rho0*c0^2 = " << matI.rho_0*matI.c0*matI.c0/1e9
            << " GPa" << endl;
    }
    Info<< endl;
}

// -------------------------------------------------------------------------
//  Numerical parameters
// -------------------------------------------------------------------------
const scalar Tref_EOS
(
    GrainProperties.lookupOrDefault<scalar>("EOS_Tref", 300.0)
);
const scalar qQuadCoeff
(
    GrainProperties.lookupOrDefault<scalar>("EOS_qQuad", 1.2)
);
const scalar qLinCoeff
(
    GrainProperties.lookupOrDefault<scalar>("EOS_qLinear", 0.06)
);

// Identity tensor helper for the volumetric-stress assembly
const dimensionedTensor I_EOS("I_EOS", dimless, tensor::I);

// Compact momentum-stabiliser coefficient (see header note above).
const scalar waveStabCoeff
(
    GrainProperties.lookupOrDefault<scalar>("EOS_waveStabiliser", 0.0)
);


volScalarField muWaveEOS
(
    IOobject
    (
        "muWaveEOS",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("muWaveEOS", dimensionSet(1, -1, -1, 0, 0), 0.0),
    zeroGradientFvPatchScalarField::typeName
);

// Acoustic Courant number, refreshed in EOS.H, consumed in setDeltaT.H
scalar acousticCoNum = 0.0;


volScalarField CvInitEOS
(
    IOobject
    (
        "CvInitEOS",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("CvInitEOS", dimensionSet(0, 2, -2, -1, 0), 0.0)
);
{
    const dimensionedScalar dimcvEOS
    (
        "dimcvEOS", dimensionSet(0, 2, -2, -1, 0), 1.0
    );
    for (label hEOS = 0; hEOS < N_Ori; hEOS++)
    {
        CvInitEOS +=
            PopBal[hEOS]*MatList[Phase2Material[hEOS]].Cv*dimcvEOS;
    }
}

// Specific internal energy, measured from  reference state
// (rho = rho0Mix, T = Tref_EOS, p = 0):  e(t=0) = Cv (T - Tref).
const bool eIntFileFound
(
    IOobject
    (
        "e",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT
    ).typeHeaderOk<volScalarField>(true)
);

volScalarField eInt
(
    IOobject
    (
        "e",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("e", dimensionSet(0, 2, -2, 0, 0), 0.0),
    zeroGradientFvPatchScalarField::typeName
);

if (!eIntFileFound)
{
    eInt = CvInitEOS
          *(T - dimensionedScalar("TrefEOS", dimTemperature, Tref_EOS));
    if (useEOS)
    {
        Info<< "No e field found: initialising e = Cv (T - Tref), "
            << "Tref = " << Tref_EOS << " K" << nl << endl;
    }
}
eInt.correctBoundaryConditions();

// Artificial (shock) viscosity, added to the pressure to stop unphysical oscillations.
volScalarField qAV
(
    IOobject
    (
        "qAV",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("qAV", dimensionSet(1, -1, -2, 0, 0), 0.0),
    zeroGradientFvPatchScalarField::typeName
);

// State-dependent bulk/longitudinal wave-speed estimate (CFL, q).
volScalarField cEOS
(
    IOobject
    (
        "cEOS",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("cEOS", dimVelocity, 0.0),
    zeroGradientFvPatchScalarField::typeName
);

// Phase-mixed EOS parameters (recomputed in EOS.H each call).
#define MAKE_EOS_MIX_FIELD(NAME, DIMS)                                      \
volScalarField NAME                                                         \
(                                                                           \
    IOobject                                                                \
    (                                                                       \
        #NAME,                                                              \
        runTime.timeName(),                                                 \
        mesh,                                                               \
        IOobject::NO_READ,                                                  \
        IOobject::NO_WRITE                                                  \
    ),                                                                      \
    mesh,                                                                   \
    dimensionedScalar(#NAME, DIMS, 0.0),                                    \
    zeroGradientFvPatchScalarField::typeName                                \
);

MAKE_EOS_MIX_FIELD(rho0Mix, dimDensity)
MAKE_EOS_MIX_FIELD(A0Mix,   dimensionSet(1, -1, -2, 0, 0))   // rho0*c0^2
MAKE_EOS_MIX_FIELD(sMix,    dimless)
MAKE_EOS_MIX_FIELD(G0Mix,   dimless)
MAKE_EOS_MIX_FIELD(cL0Mix,  dimVelocity)

#undef MAKE_EOS_MIX_FIELD

// Characteristic cell size for the artificial viscosity.

volScalarField cellDx
(
    IOobject
    (
        "cellDx",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("cellDx", dimLength, 0.0),
    zeroGradientFvPatchScalarField::typeName
);
{
    const cellList& cellFaceLists = mesh.cells();
    const vectorField& faceAreasEOS = mesh.faceAreas();
    scalarField& dxIF = cellDx.primitiveFieldRef();

    forAll(cellFaceLists, celli)
    {
        scalar Amax = SMALL;
        const cell& cFaces = cellFaceLists[celli];
        forAll(cFaces, fi)
        {
            Amax = max(Amax, mag(faceAreasEOS[cFaces[fi]]));
        }
        dxIF[celli] = mesh.V()[celli]/Amax;
    }
}
cellDx.correctBoundaryConditions();

// ************************************************************************* //
