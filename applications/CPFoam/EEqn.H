//FOR STRONG SHOCK WE NEED TO EXPLICITALLY SOLVE FOR INTERNAL ENERGY

//  Temperature recoveed algebraically in EOS>H
//  (T = Tref + (e - e_ref(mu))/Cv)
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

} // end if (stressUpdate && useEOS)

// ************************************************************************* //
