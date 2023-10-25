{

    // rho.storePrevIter();

        fvScalarMatrix  rhoEqn
        (
            fvm::ddt(rho)
            + fvm::div(phi, rho)

        );

        rhoEqn.relax();

        rhoEqn.solve();


// Acc_Strain.storePrevIter();
// TESD = fvc::div(phi, Acc_Strain);



        fvScalarMatrix  AdvAccStrainEqn
        (
            fvm::ddt(Acc_Strain)
            + fvm::div(phi, Acc_Strain)
            ==
            epsilon_dot_P_VM

        );

        AdvAccStrainEqn.relax();

        AdvAccStrainEqn.solve();
        // Acc_Strain.correctBoundaryConditions();


Acc_Strain=max(Acc_Strain,0.0);


Sum_Phases*=0.0;
for (label count = 0; count < N_Ori; count++)
{

    volScalarField & ni = PopBal[count];
    // volScalarField& Su=Sum_parts[count];
    // volScalarField & RHS = dpdt[count];

    // ni=ni+runTime.deltaT()*RHS;

    solve
    (
        fvm::ddt(ni) 
        +fvm::div(phi, ni)
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