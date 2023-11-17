{

    // rho.storePrevIter();

        fvScalarMatrix  rhoEqn
        (
            fvm::ddt(rho)
            + fvm::div(phi, rho)
            // - fvm::SuSp( fvc::div(phi), rho)

        );

        rhoEqn.relax();

        rhoEqn.solve();


}