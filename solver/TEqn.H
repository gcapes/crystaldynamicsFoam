{

    
    const surfaceScalarField rhophicv(/*fvc::interpolate(Cv)*/ Cv*rhoPhi);

    const volScalarField rhoCv(rho*Cv);

    const volVectorField SigU(sigmaGlobal & U);

    // const volScalarField 
    // DP = (sigma && fvc::grad(U));
    // DP.correctBoundaryConditions();

    // DP2 = fvc::div(SigU);



        fvScalarMatrix  TEqn
        (
            fvm::ddt(rhoCv,T)
            +fvc::ddt(rho, K)
            + fvm::div(rhophicv, T)
            +fvc::div(rhoPhi,K)
            -fvc::div(SigU,"div((sigma&U))")

            // - fvm::Sp(fvc::ddt(rhoCv) + fvc::div(rhophicv), T)
            // - fvc::Sp(fvc::ddt(rho) + fvc::div(rhoPhi), K)
            - fvm::laplacian(kappa, T)

        );

        TEqn.relax();
        // fvConstraints.constrain(TEqn);

        TEqn.solve();

}





