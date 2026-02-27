
{
    // Continuity error: should be zero if rho satisfied discrete mass conservation
    volScalarField contErr
    (
        "contErr",
        fvc::ddt(rho) + fvc::div(rhoPhi)
    );


    // ==================== Advect GLOBAL deformation gradient F ====================
    {
        volTensorField L = fvc::grad(U);  // Velocity gradient (spatial)

        solve(
            fvm::ddt(rho, F) + fvm::div(rhoPhi, F)
            - fvm::SuSp(contErr, F)       // Continuity error correction
            ==
            rho * (L & F)                  // Source: L.F = (grad(v)).F
        );

        F.correctBoundaryConditions();

        // Info << "F: max(det(F)) = " << gMax(det(F))
        //      << ", min(det(F)) = " << gMin(det(F)) << endl;
    }


    // ==================== Advect Fp per phase ====================
    for (label h = 0; h < N_Ori; h++)
    {
        volTensorField& Fp = F_plastic[h];

        solve(
            fvm::ddt(rho, Fp) + fvm::div(rhoPhi, Fp)
            - fvm::SuSp(contErr, Fp)      // Continuity error correction
            ==
            rho * (Lp[h] & Fp)            // Source: Lp.Fp
        );

        Fp.correctBoundaryConditions();

        // Info << "Fp." << h << ": max(det(Fp)) = " << gMax(det(Fp))
        //      << ", min(det(Fp)) = " << gMin(det(Fp)) << endl;
    }
}