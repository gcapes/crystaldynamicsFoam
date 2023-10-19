{

        fvScalarMatrix  rhoEqn
        (
            fvm::ddt(rho)
            + fvm::div(phi, rho)

        );

        rhoEqn.relax();

        rhoEqn.solve();

}