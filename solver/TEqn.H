{

    
    const surfaceScalarField rhophicv(/*fvc::interpolate(Cv)*/ Cv*rhoPhi);

    const volScalarField rhoCv(rho*Cv);

    const volVectorField SigU(sigma & U);

    const volScalarField DP(sigma && fvc::grad(U));



        fvScalarMatrix  TEqn
        (
            fvm::ddt(rhoCv,T)
            +fvc::ddt(rho, K)
            + fvm::div(rhophicv, T)
            +fvc::div(rhoPhi,K)
            -fvc::div(SigU)
            // -DP
            - fvm::Sp(fvc::ddt(rhoCv) + fvc::div(rhophicv), T)
            - fvc::Sp(fvc::ddt(rho) + fvc::div(rhoPhi), K)
            - fvm::laplacian(kappa, T)

        );

        TEqn.relax();
        // fvConstraints.constrain(TEqn);

        TEqn.solve();

}





