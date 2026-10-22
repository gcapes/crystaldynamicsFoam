// Linear-elastic mode: conservative energy equation solved for T
if (stressUpdate && !useEOS)
{


    const surfaceScalarField rhophicv(fvc::interpolate(Cv)*rhoPhi);

    // const volScalarField rhoCv(rho*Cv); // must not be const: ddt(rhoCv) would be dropped
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

    // Boundary terms of the discrete energy balance from the same face fluxes
    // (final outer iteration; integrated in energyBudget.H).  S_f points outward.
    if (energyBudgetActive && pimple.finalIter())
    {
        ebBoundaryStressPower = 0.0;
        ebBoundaryConvective  = 0.0;
        ebBoundaryHeat        = 0.0;
        ebBoundaryMassFlux    = 0.0;
        forAll(mesh.boundary(), patchi)
        {
            const fvPatch& pp = mesh.boundary()[patchi];
            if (pp.coupled() || pp.size() == 0) continue;
            const vectorField& Sfp = mesh.Sf().boundaryField()[patchi];
            const scalarField& rhoPhip = rhoPhi.boundaryField()[patchi];
            ebBoundaryStressPower += sum(SigU.boundaryField()[patchi] & Sfp);
            ebBoundaryConvective  += sum
            (
                rhoPhip*(Cv.boundaryField()[patchi]*T.boundaryField()[patchi] + K.boundaryField()[patchi])
            );
            ebBoundaryMassFlux    += sum(rhoPhip);
            ebBoundaryHeat        += sum
            (
                kappa.boundaryField()[patchi]
               *T.boundaryField()[patchi].snGrad()*pp.magSf()
            );
        }
        reduce(ebBoundaryStressPower, sumOp<scalar>());
        reduce(ebBoundaryConvective,  sumOp<scalar>());
        reduce(ebBoundaryHeat,        sumOp<scalar>());
        reduce(ebBoundaryMassFlux,    sumOp<scalar>());
    }

} // end if (stressUpdate)





