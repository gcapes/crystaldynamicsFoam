// Energy equation: only active when momentum/CP solve is running.
if (stressUpdate)
{


    const surfaceScalarField rhophicv(fvc::interpolate(Cv)*rhoPhi);

    // const volScalarField rhoCv(rho*Cv);// NEEDS TO NOT BE CONST OTHERWISE DROPS TIME DERVATIVE IN RHOCV!!!!!!!!!!
    rhoCv = rho*Cv;

    const volVectorField SigU(sigmaGlobal & U);

    // const volScalarField 
    // DP = (sigma && fvc::grad(U));
    // DP.correctBoundaryConditions();

    // DP2 = fvc::div(SigU);



        fvScalarMatrix  TEqn
        (
            fvm::ddt(rhoCv,T)
            +fvc::ddt(rho, K)
            + fvm::div(rhophicv, T,"div((Cv*rhoPhi),T)")
            +fvc::div(rhoPhi,K)
            -fvc::div(SigU,"div((sigma&U))")

            // - fvm::Sp(fvc::ddt(rhoCv) + fvc::div(rhophicv), T)
            // - fvc::Sp(fvc::ddt(rho) + fvc::div(rhoPhi), K)
            - fvm::laplacian(kappa, T)

        );

        TEqn.relax();
        // fvConstraints.constrain(TEqn);

        TEqn.solve();

} // end if (stressUpdate)





