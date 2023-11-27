{

    // Acc_Strain.storePrevIter();
    // TESD = fvc::div(phi, Acc_Strain);

    // fvScalarMatrix  AdvAccStrainEqn
    // (
    //     fvm::ddt(rho,Acc_Strain)
    //     + fvm::div(rhoPhi, Acc_Strain)
    //     ==
    //     (rho*sumgamma)//epsilon_dot_P_VM

    // );

    // AdvAccStrainEqn.relax();

    // AdvAccStrainEqn.solve();
    // Acc_Strain.correctBoundaryConditions();

    // Acc_Strain=max(Acc_Strain,0.0);

    // fvTensorMatrix  OmegaTotEqn
    // (
    //     fvm::ddt(omega_tot)
    //     + fvm::div(phi, omega_tot)
    //     ==
    //     Omega_E

    // );

    // OmegaTotEqn.relax();

    // OmegaTotEqn.solve();
    // // Acc_Strain.correctBoundaryConditions();

    for (label i = 0; i < N_Ori; i++)
    {
        for (label countS = 0; countS < MatList[Phase2Material[i]].N_Slips; countS++)
        {

            // volScalarField & ni = PopBal[count];
            volScalarField &g_a = g_alphas[i][countS];

            solve(
                fvm::ddt(rho, g_a) + fvm::div(rhoPhi, g_a) ==
                (rho * g_alpha_dots[i][countS])
                // RHS
            );
            // g_a.correctBoundaryConditions();
            // Sum_Phases += ni;

            // Info << "Max phase field(" << count << ") = " << gMax(ni) << endl;
            // Info << "Min phase field(" << count << ") = " << gMin(ni) << endl;
        }
    }

    for (label i = 0; i < N_Ori; i++)
    {
        for (label countS = 0; countS < MatList[Phase2Material[i]].N_Slips; countS++)
        {

            // volScalarField & ni = PopBal[count];
            volScalarField &gam_a = gamma_alphas[i][countS];

            solve(
                fvm::ddt(rho, gam_a) + fvm::div(rhoPhi, gam_a) ==
                (rho * gamma_alpha_dots[i][countS])
                // RHS
            );
            // g_a.correctBoundaryConditions();
            // Sum_Phases += ni;

            // Info << "Max phase field(" << count << ") = " << gMax(ni) << endl;
            // Info << "Min phase field(" << count << ") = " << gMin(ni) << endl;
        }
    }

    Sum_Phases *= 0.0;
    for (label count = 0; count < N_Ori; count++)
    {

        volScalarField &ni = PopBal[count];
        // volScalarField& Su=Sum_parts[count];
        // volScalarField & RHS = dpdt[count];

        // ni=ni+runTime.deltaT()*RHS;

        solve(
            fvm::ddt(rho, ni) + fvm::div(rhoPhi, ni)
            // - fvm::SuSp(fvc::ddt(rho) + fvc::div(rhoPhi), ni)
            // RHS
        );
        ni.correctBoundaryConditions();
        Sum_Phases += ni;

        Info << "Max phase field(" << count << ") = " << gMax(ni) << endl;
        Info << "Min phase field(" << count << ") = " << gMin(ni) << endl;
    }
    Sum_Phases.correctBoundaryConditions();

    //  Project any out of bounds phases back onto the Gibbs Simplex
    // for (label count = 0; count < N_Ori; count++)
    // {

    //     volScalarField & ni = PopBal[count];

    //     ni = max(min(ni, 1.0), 0.0);
    //     volScalarField temp_holder = ni;

    //     forAll(mesh.C(), celli)
    //     {
    //         if (Sum_Phases[celli] > 0.0)
    //         {
    //             temp_holder[celli] *= (1.0 / Sum_Phases[celli]);
    //         }
    //     }
    //     ni = temp_holder;
    // }

    for (label count = 0; count < N_Ori; count++)
    {

        volScalarField &Rot1_ = Rot1[count];
        volScalarField &Rot2_ = Rot2[count];
        volScalarField &Rot3_ = Rot3[count];

        solve(
            fvm::ddt(rho, Rot1_) + fvm::div(rhoPhi, Rot1_) 
            // ==
            // (rho * (omega_dot.component(7) - Omega_P[count].component(7))) // epsilon_dot_P_VM

        );
        Rot1_.correctBoundaryConditions();

        solve(
            fvm::ddt(rho, Rot2_) + fvm::div(rhoPhi, Rot2_) 
            // ==
            // (rho * (omega_dot.component(2) - Omega_P[count].component(2))) // epsilon_dot_P_VM

        );
        Rot2_.correctBoundaryConditions();

        solve(
            fvm::ddt(rho, Rot3_) + fvm::div(rhoPhi, Rot3_) 
            // ==
            // (rho * (omega_dot.component(3) - Omega_P[count].component(3))) // epsilon_dot_P_VM

        );
        Rot3_.correctBoundaryConditions();
    }

    // fvScalarMatrix  ROT1Eqn
    // (
    //     fvm::ddt(rho,ROT_1)
    //     + fvm::div(rhoPhi, ROT_1)
    //     ==
    //     (rho*(omega_dot.component(5) - Omega_P.component(5)))//epsilon_dot_P_VM

    // );

    // ROT1Eqn.relax();

    // ROT1Eqn.solve();

    // fvScalarMatrix  ROT2Eqn
    // (
    //     fvm::ddt(rho,ROT_2)
    //     + fvm::div(rhoPhi, ROT_2)
    //     ==
    //     (rho*(omega_dot.component(2) - Omega_P.component(2)))//epsilon_dot_P_VM

    // );

    // ROT2Eqn.relax();

    // ROT2Eqn.solve();

    //  fvScalarMatrix  ROT3Eqn
    // (
    //     fvm::ddt(rho,ROT_3)
    //     + fvm::div(rhoPhi, ROT_3)
    //     ==
    //     (rho*(omega_dot.component(1) - Omega_P.component(1)))//epsilon_dot_P_VM

    // );

    // ROT3Eqn.relax();

    // ROT3Eqn.solve();
}