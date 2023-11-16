{




// Acc_Strain.storePrevIter();
// TESD = fvc::div(phi, Acc_Strain);



        fvScalarMatrix  AdvAccStrainEqn
        (
            fvm::ddt(rho,Acc_Strain)
            + fvm::div(rhoPhi, Acc_Strain)
            ==
            (rho*sumgamma)//epsilon_dot_P_VM

        );

        AdvAccStrainEqn.relax();

        AdvAccStrainEqn.solve();
        // Acc_Strain.correctBoundaryConditions();
        
// Acc_Strain=max(Acc_Strain,0.0);




        fvTensorMatrix  OmegaTotEqn
        (
            fvm::ddt(omega_tot)
            + fvm::div(phi, omega_tot)
            ==
            Omega_E

        );

        OmegaTotEqn.relax();

        OmegaTotEqn.solve();
        // Acc_Strain.correctBoundaryConditions();



for (label countS = 0; countS < N_Slips; countS++)
{

    // volScalarField & ni = PopBal[count];
    volScalarField & g_a =g_alphas[countS];



    solve
    (
        fvm::ddt(rho,g_a) 
        +fvm::div(rhoPhi, g_a)
        ==
        (rho*g_alpha_dots[countS])
        // RHS
    );
    // g_a.correctBoundaryConditions();
    // Sum_Phases += ni;

    // Info << "Max phase field(" << count << ") = " << gMax(ni) << endl;
    // Info << "Min phase field(" << count << ") = " << gMin(ni) << endl;

}


        


Sum_Phases*=0.0;
for (label count = 0; count < N_Ori; count++)
{

    volScalarField & ni = PopBal[count];
    // volScalarField& Su=Sum_parts[count];
    // volScalarField & RHS = dpdt[count];

    // ni=ni+runTime.deltaT()*RHS;

    solve
    (
        fvm::ddt(rho,ni) 
        +fvm::div(rhoPhi, ni)
        // RHS
    );
    ni.correctBoundaryConditions();
    Sum_Phases += ni;

    Info << "Max phase field(" << count << ") = " << gMax(ni) << endl;
    Info << "Min phase field(" << count << ") = " << gMin(ni) << endl;

}
Sum_Phases.correctBoundaryConditions();


//  Project any out of bounds phases back onto the Gibbs Simplex
for (label count = 0; count < N_Ori; count++)
{

    volScalarField & ni = PopBal[count];


    ni = max(min(ni, 1.0), 0.0);
    volScalarField temp_holder = ni;

    forAll(mesh.C(), celli)
    {
        if (Sum_Phases[celli] > 0.0)
        {
            temp_holder[celli] *= (1.0 / Sum_Phases[celli]);
        }
    }
    ni = temp_holder;
}



}