// Mie-Gruneisen mode: conservative energy equation solved for e;
// T is recovered algebraically in EOS.H.
if (stressUpdate && useEOS)
{
    const volVectorField SigU(sigmaGlobaltemp & U);

    const volScalarField kByCv
    (
        "kByCv",
        kappa/max(Cv, dimensionedScalar("CvSmall", Cv.dimensions(), SMALL))
    );

    fvScalarMatrix EEqn
    (
        fvm::ddt(rho, eInt)
      + fvm::div(rhoPhi, eInt, "div(rhoPhi,e)")
      + fvc::ddt(rho, K)
      + fvc::div(rhoPhi, K)
      - fvc::div(SigU, "div((sigma&U))")
      - fvm::laplacian(kByCv, eInt)
      - (fvc::laplacian(kappa, T) - fvc::laplacian(kByCv, eInt))
    );

    EEqn.relax();

    EEqn.solve();

    eInt.correctBoundaryConditions();

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
                rhoPhip*(eInt.boundaryField()[patchi] + K.boundaryField()[patchi])
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

} // end if (stressUpdate && useEOS)

// ************************************************************************* //
