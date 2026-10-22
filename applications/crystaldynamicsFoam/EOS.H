
// Mie-Gruneisen pressure p(rho, e) on the Hugoniot reference curve (compression) or its
// linear extrapolation (tension); sigma = -(p + q) I + sum_h phi_h s_h, q the shock viscosity.


if (useEOS)
{

    rho0Mix *= 0.0;
    A0Mix   *= 0.0;
    sMix    *= 0.0;
    G0Mix   *= 0.0;
    cL0Mix  *= 0.0;
    {
        const dimensionedScalar dimRho0EOS("dimRho0EOS", dimDensity, 1.0);
        const dimensionedScalar dimA0EOS
        (
            "dimA0EOS", dimensionSet(1, -1, -2, 0, 0), 1.0
        );
        const dimensionedScalar dimCEOS("dimCEOS", dimVelocity, 1.0);

        for (label hEOS = 0; hEOS < N_Ori; hEOS++)
        {
            const Material& matE = MatList[Phase2Material[hEOS]];

            rho0Mix += PopBal[hEOS]*matE.rho_0*dimRho0EOS;
            A0Mix   += PopBal[hEOS]
                      *(matE.rho_0*matE.c0*matE.c0)*dimA0EOS;
            sMix    += PopBal[hEOS]*matE.s_UsUp;
            G0Mix   += PopBal[hEOS]*matE.Gamma0;
            cL0Mix  += PopBal[hEOS]*matE.cL0_bound*dimCEOS;
        }
    }

    label nMuClamped = 0;
    {
        scalarField& pIF        = p.primitiveFieldRef();
        scalarField& cIF        = cEOS.primitiveFieldRef();
        scalarField& TIF        = T.primitiveFieldRef();
        const scalarField& rhoIF = rho.primitiveField();
        const scalarField& eIF   = eInt.primitiveField();
        const scalarField& r0IF  = rho0Mix.primitiveField();
        const scalarField& A0IF  = A0Mix.primitiveField();
        const scalarField& sIF   = sMix.primitiveField();
        const scalarField& G0IF  = G0Mix.primitiveField();
        const scalarField& cLIF  = cL0Mix.primitiveField();
        const scalarField& CvIF  = Cv.primitiveField();

        forAll(pIF, celli)
        {
            const scalar r0 = max(r0IF[celli], SMALL);
            const scalar A0 = max(A0IF[celli], SMALL);
            const scalar sU = sIF[celli];
            const scalar G0 = G0IF[celli];

            scalar mu = rhoIF[celli]/r0 - 1.0;

            if (sU > 1.0)
            {
                const scalar muMax = 0.95/(sU - 1.0);
                if (mu > muMax) { mu = muMax; nMuClamped++; }
            }

            scalar pH, eRef, dpHdmu;
            if (mu >= 0.0)
            {
                const scalar Dn = 1.0 - (sU - 1.0)*mu;
                pH     = A0*mu*(1.0 + mu)/(Dn*Dn);
                eRef   = 0.5*pH*mu/(r0*(1.0 + mu));
                dpHdmu = A0*((1.0 + 2.0*mu)*Dn
                       + 2.0*(sU - 1.0)*mu*(1.0 + mu))/(Dn*Dn*Dn);
            }
            else
            {
                // Tension
                pH     = A0*mu;
                eRef   = 0.5*(A0/r0)*mu*mu;
                dpHdmu = A0;
            }

            pIF[celli] = pH + G0*r0*(eIF[celli] - eRef);

            // Longitudinal wave-speed estimate, stiffened by the Hugoniot slope
            const scalar stiff = Foam::sqrt(max(dpHdmu/A0, 1.0));
            cIF[celli] = max(cLIF[celli]*stiff, Foam::sqrt(A0/r0));


            const scalar CvLoc = max(CvIF[celli], SMALL);
            TIF[celli] = max
            (
                Tref_EOS + (eIF[celli] - eRef)/CvLoc,
                1.0
            );
        }
    }

    if (returnReduce(nMuClamped, sumOp<label>()) > 0)
    {
        Info<< "EOS: mu clamped to 0.95/(s-1) in "
            << returnReduce(nMuClamped, sumOp<label>())
            << " cells (beyond Hugoniot-form validity)" << endl;
    }

    p.correctBoundaryConditions();
    cEOS.correctBoundaryConditions();
    T.correctBoundaryConditions();     // respects case-supplied T BCs


    // Artificial (shock) viscosity: compression only
    {
        const volScalarField trDEOS(tr(epsilon_dot));

        scalarField& qIF         = qAV.primitiveFieldRef();
        const scalarField& trDIF = trDEOS.primitiveField();
        const scalarField& dxIF  = cellDx.primitiveField();
        const scalarField& rhoIF = rho.primitiveField();
        const scalarField& cIF   = cEOS.primitiveField();

        forAll(qIF, celli)
        {
            const scalar comp = min(trDIF[celli], 0.0);
            qIF[celli] = rhoIF[celli]*dxIF[celli]*
            (
                qQuadCoeff*dxIF[celli]*comp*comp
              + qLinCoeff*cIF[celli]*mag(comp)
            );
        }
    }
    qAV.correctBoundaryConditions();

    if (waveStabCoeff > 0.0)
    {
        muWaveEOS = waveStabCoeff*rho*cEOS*cellDx;
        muWaveEOS.correctBoundaryConditions();
    }


}
else if (acousticCourantControl && stressUpdate && !quasiStatic)
{
    // Linear-elastic mode: c = sqrt(sum_h phi_h Gamma_max,h / rho) from the
    // Christoffel bound of each material (createFields.H).
    volScalarField GamMixEOS
    (
        IOobject("GamMixEOS", runTime.timeName(), mesh),
        mesh,
        dimensionedScalar("GamMixEOS", dimensionSet(1, -1, -2, 0, 0), 0.0)
    );
    for (label hEOS = 0; hEOS < N_Ori; hEOS++)
    {
        GamMixEOS += PopBal[hEOS]
                    *dimensionedScalar
                     (
                         "GamMax",
                         dimensionSet(1, -1, -2, 0, 0),
                         GammaMaxMat[Phase2Material[hEOS]]
                     );
    }
    cEOS = sqrt(GamMixEOS/max(rho, dimensionedScalar("rhoSmall", dimDensity, SMALL)));
    cEOS.correctBoundaryConditions();
}

// Acoustic Courant number for setDeltaT.H (linear-elastic mode: only with
// acousticCourantControl on and a non-quasi-static coupling).
if (useEOS || (acousticCourantControl && stressUpdate && !quasiStatic))
{
    const surfaceScalarField cfEOS(fvc::interpolate(cEOS));

    acousticCoNum =
        0.5*gMax
        (
            fvc::surfaceSum(mag(phi) + cfEOS*mesh.magSf())()
                .primitiveField()
           /mesh.V().field()
        )*runTime.deltaTValue();

    Info<< "Acoustic Courant Number max: " << acousticCoNum << endl;
}
