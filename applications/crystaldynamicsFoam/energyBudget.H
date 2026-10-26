// Domain energy, mass and dissipation budget on the end-of-step state, written to
// postProcessing/energyBudget.dat: E - E0 = Wb - Fb + Qb.  Switch: energyBudget.

if (stressUpdate && energyBudgetActive)
{
    const scalarField& Vc = mesh.V();
    const scalar dt = runTime.deltaTValue();

    // --- Kinetic, internal and total energy, mass ---
    const scalar KE = gSum(0.5*rho.primitiveField()
                          *magSqr(U.primitiveField())*Vc);

    scalar IE = 0.0;
    if (useEOS)
    {
        IE = gSum(rho.primitiveField()*eInt.primitiveField()*Vc);
    }
    else
    {
        IE = gSum(rho.primitiveField()*Cv.primitiveField()
                 *T.primitiveField()*Vc);
    }
    const scalar Etot = KE + IE;
    const scalar Mtot = gSum(rho.primitiveField()*Vc);

    // --- Irreversible work rates (W) ---
    scalar DqRate = 0.0;
    if (useEOS)
    {
        const scalarField trD(tr(epsilon_dot.primitiveField()));
        DqRate = gSum(qAV.primitiveField()*max(-trD, 0.0)*Vc);
    }

    plasticWorkRate.primitiveFieldRef() = 0.0;
    for (label h = 0; h < N_Ori; h++)
    {
        plasticWorkRate.primitiveFieldRef() +=
            PopBal[h].primitiveField()*plasticWorkRatePhase[h];
    }
    plasticWorkRate.correctBoundaryConditions();
    const scalar DpRate = gSum(plasticWorkRate.primitiveField()*Vc);

    // --- Time integrals ---
    if (!ebInitialised)
    {
        // First step: reconstruct the start-of-run energy and mass from this
        // step's increments.
        ebE0 = Etot
             - dt*(ebBoundaryStressPower - ebBoundaryConvective
                 + ebBoundaryHeat);
        ebM0 = Mtot + dt*ebBoundaryMassFlux;
        ebInitialised = true;

        if (Pstream::master())
        {
            mkDir(runTime.path()/"postProcessing");
            std::ofstream os
            (
                (runTime.path()/"postProcessing"/"energyBudget.dat").c_str()
            );
            os  << "# crystaldynamicsFoam energy budget (SI units; S_f outward)\n"
                << "# E = KE + IE; identity E - E0 = Wb - Fb + Qb\n"
                << "# 1:time 2:KE 3:IE 4:E 5:Wb(boundary stress work) "
                << "6:Fb(convective energy outflow) 7:Qb(conductive inflow) "
                << "8:Dq(artificial-viscosity work) 9:Dp(plastic work) "
                << "10:residual(normalised) 11:mass 12:massInflow "
                << "13:massResidual(normalised) 14:DqRate 15:DpRate\n";
        }
    }

    ebWork   += dt*ebBoundaryStressPower;
    ebConv   += dt*ebBoundaryConvective;
    ebHeat   += dt*ebBoundaryHeat;
    ebDq     += dt*DqRate;
    ebDp     += dt*DpRate;
    ebMassIn -= dt*ebBoundaryMassFlux;

    const scalar rhs = ebWork - ebConv + ebHeat;
    const scalar scaleE = max(max(mag(ebE0), mag(rhs)), mag(Etot - ebE0));
    const scalar residual = (Etot - ebE0 - rhs)/max(scaleE, SMALL);
    const scalar massResidual = (Mtot - ebM0 - ebMassIn)/max(mag(ebM0), SMALL);

    Info<< "Energy budget: KE = " << KE << " J, IE = " << IE
        << " J, E - E0 = " << Etot - ebE0
        << " J, Wb - Fb + Qb = " << rhs
        << " J, residual = " << residual
        << ", Dq = " << ebDq << " J, Dp = " << ebDp << " J" << endl;

    if (Pstream::master())
    {
        std::ofstream os
        (
            (runTime.path()/"postProcessing"/"energyBudget.dat").c_str(),
            std::ios::app
        );
        os.precision(12);
        os  << runTime.timeName() << ' '
            << KE << ' ' << IE << ' ' << Etot << ' '
            << ebWork << ' ' << ebConv << ' ' << ebHeat << ' '
            << ebDq << ' ' << ebDp << ' '
            << residual << ' '
            << Mtot << ' ' << ebMassIn << ' ' << massResidual << ' '
            << DqRate << ' ' << DpRate << '\n';
    }
}

// ************************************************************************* //
