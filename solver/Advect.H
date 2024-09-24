{



    fvScalarMatrix  AdvAccStrainEqn
    (
        fvm::ddt(rho,Acc_Strain)
        + fvm::div(rhoPhi, Acc_Strain)
        - fvm::SuSp(fvc::ddt(rho) + fvc::div(rhoPhi), Acc_Strain) 
        ==
        (rho*sumgamma)//epsilon_dot_P_VM

    );

    AdvAccStrainEqn.relax();

    AdvAccStrainEqn.solve();
    Acc_Strain.correctBoundaryConditions();


// //FOR HARDENING
    for (label i = 0; i < N_Ori; i++)//advect hardened fields
    {
        for (label countS = 0; countS < MatList[Phase2Material[i]].N_Slips; countS++)
        {

            // volScalarField & ni = PopBal[count];
            volScalarField &g_a = g_alphas[i][countS];

            solve(
                fvm::ddt(rho, g_a) + fvm::div(rhoPhi, g_a)
                - fvm::SuSp(fvc::ddt(rho) + fvc::div(rhoPhi), g_a) 
                 ==
                (rho * g_alpha_dots[i][countS])
                // RHS
            );
            // g_a.correctBoundaryConditions();
            // Sum_Phases += ni;

            // Info << "Max phase field(" << count << ") = " << gMax(ni) << endl;
            // Info << "Min phase field(" << count << ") = " << gMin(ni) << endl;
        }
    }
    // advect hardedned fields

    // for (label i = 0; i < N_Ori; i++)
    // {
    //     for (label countS = 0; countS < MatList[Phase2Material[i]].N_Slips; countS++)
    //     {

    //         // volScalarField & ni = PopBal[count];
    //         volScalarField &gam_a = gamma_alphas[i][countS];

    //         solve(
    //             fvm::ddt(rho, gam_a) + fvm::div(rhoPhi, gam_a) ==
    //             (rho * gamma_alpha_dots[i][countS])
    //             // RHS
    //         );
    //         // g_a.correctBoundaryConditions();
    //         // Sum_Phases += ni;

    //         // Info << "Max phase field(" << count << ") = " << gMax(ni) << endl;
    //         // Info << "Min phase field(" << count << ") = " << gMin(ni) << endl;
    //     }
    // }


    for (label count = 0; count < N_Ori; count++)
    {

        volScalarField &Rot1_ = Rot1[count];
        volScalarField &Rot2_ = Rot2[count];
        volScalarField &Rot3_ = Rot3[count];

        solve(
            fvm::ddt(rho, Rot1_) + fvm::div(rhoPhi, Rot1_)
            - fvm::SuSp(fvc::ddt(rho) + fvc::div(rhoPhi), Rot1_) 


        );
        Rot1_.correctBoundaryConditions();

        solve(
            fvm::ddt(rho, Rot2_) + fvm::div(rhoPhi, Rot2_) 
            - fvm::SuSp(fvc::ddt(rho) + fvc::div(rhoPhi), Rot2_)


        );
        Rot2_.correctBoundaryConditions();

        solve(
            fvm::ddt(rho, Rot3_) + fvm::div(rhoPhi, Rot3_)
            - fvm::SuSp(fvc::ddt(rho) + fvc::div(rhoPhi), Rot3_) 


        );
        Rot3_.correctBoundaryConditions();
    }

    
}