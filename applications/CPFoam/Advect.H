// Advection of plastic slip/hardening fields: only meaningful when the CP solve is active.
if (stressUpdate)
{


        for (label h = 0; h < N_Ori; h++)
{
    sigma[h].oldTime() == sigma[h];
        Rot1[h].oldTime() == Rot1[h];
    Rot2[h].oldTime() == Rot2[h];
    Rot3[h].oldTime() == Rot3[h];
}




    // Advect hardening fields (CRSS per slip system)
    for (label i = 0; i < N_Ori; i++)
    {
        const Material& mat = MatList[Phase2Material[i]];
        const double tau_0   = mat.initial_CRSS;

        for (label countS = 0; countS < MatList[Phase2Material[i]].N_Slips; countS++)
        {
            volScalarField &g_a = g_alphas[i][countS];

            solve(
                fvm::ddt(rho, g_a) + fvm::div(rhoPhi, g_a)
                 ==
                (rho * g_alpha_dots[i][countS])
            );
            g_a.correctBoundaryConditions();

            
            g_a = max(g_a, dimensionedScalar("gMin", g_a.dimensions(), tau_0 * 0.001));
        }
    }

    // Advect accumulated slip per slip system
    for (label i = 0; i < N_Ori; i++)
    {
        for (label countS = 0; countS < MatList[Phase2Material[i]].N_Slips; countS++)
        {
            volScalarField &gam_a = gamma_alphas[i][countS];

            solve(
                fvm::ddt(rho, gam_a) + fvm::div(rhoPhi, gam_a) ==
                (rho * gamma_alpha_dots[i][countS])
            );
            gam_a.correctBoundaryConditions();  // <-- was commented out
        }
    }

    // Advect Euler-angle fields (pure transport, no source)
    for (label count = 0; count < N_Ori; count++)
    {
        volScalarField &Rot1_ = Rot1[count];
        volScalarField &Rot2_ = Rot2[count];
        volScalarField &Rot3_ = Rot3[count];

        solve(
            fvm::ddt(rho, Rot1_) + fvm::div(rhoPhi, Rot1_)
        );
        Rot1_.correctBoundaryConditions();

        solve(
            fvm::ddt(rho, Rot2_) + fvm::div(rhoPhi, Rot2_)
        );
        Rot2_.correctBoundaryConditions();

        solve(
            fvm::ddt(rho, Rot3_) + fvm::div(rhoPhi, Rot3_)
        );
        Rot3_.correctBoundaryConditions();
    }

    // Post-advection: compute mixture Euler angles and plastic strain
    Euler_Angles *= 0.0;
    D_P *= 0.0;

    forAll(mesh.C(), cell)
    {
        label largest_field = -1;
        double current_largest_value = -1.0;

        for (label h = 0; h < N_Ori; h++)
        {
            if (PopBal[h][cell] > current_largest_value)
            {
                current_largest_value = PopBal[h][cell];
                largest_field = h;
            }

            D_P[cell] += PopBal[h][cell] * Plastic_Strain[h][cell];
        }

        Euler_Angles[cell].component(0) = Rot1[largest_field][cell];
        Euler_Angles[cell].component(1) = Rot2[largest_field][cell];
        Euler_Angles[cell].component(2) = Rot3[largest_field][cell];
    }





    // Transport per-phase sigma (already Jaumann-updated in StressUpdate.H)
for (label h = 0; h < N_Ori; h++)
{
    solve(
        fvm::ddt(rho, sigma[h]) + fvm::div(rhoPhi, sigma[h])
    );
    sigma[h].correctBoundaryConditions();
}

// Re-form mixture stress from the transported per-phase sigmas
sigmaGlobal *= 0.0;
for (label h = 0; h < N_Ori; h++)
{
    sigmaGlobal += PopBal[h] * sigma[h];
}
sigmaGlobal.correctBoundaryConditions();





} // end if (stressUpdate)