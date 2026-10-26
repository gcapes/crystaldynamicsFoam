if (stressUpdate)
{

// const surfaceTensorField sigma_f(fvc::interpolate(sigma));

fvVectorMatrix UEqn(
    fvm::ddt(rho, U) 
    + fvm::div(rhoPhi, U) 
    // +fvc::div(rho * (U*U))
    // - fvc::div(sigma)
    - fvc::div(sigmaGlobaltemp,"div(sigma)")
    // ==

);

// Compact wave stabiliser (EOS mode, EOS_waveStabiliser > 0)
if (useEOS && waveStabCoeff > 0.0)
{
    UEqn -= fvm::laplacian(muWaveEOS, U);
}

// stabilised / quasiStatic: Jasak-Weller split, implicit -div(K dt grad U^(k+1)) plus the
// same operator on U^k explicitly; cancels at outer convergence, damps the checkerboard mode.
if (compactStabiliser)
{
    scalarField& impKIF = impK.primitiveFieldRef();
    impKIF = 0.0;
    for (label h = 0; h < N_Ori; h++)
    {
        const label m = Phase2Material[h];

        // EOS mode: the pressure replaces the volumetric part of C', so the bound rises
        // by max(0, rho_0 c0^2 - K_C').
        const scalar GammaQS =
            GammaMaxMat[m]
          + (useEOS ? max(0.0, MatList[m].rho_0*sqr(MatList[m].c0) - KCprimeMat[m]) : 0.0);

        impKIF += (impKFactor*GammaQS)*PopBal[h].primitiveField();
    }
    impK.correctBoundaryConditions();

    const volScalarField impKdt("impKdt", runTime.deltaT()*impK);

    UEqn -= fvm::laplacian(impKdt, U);
    UEqn += fvc::div(impKdt*fvc::grad(U), "div(sigma)");
}

UEqn.relax();



// volScalarField rAU(1.0 / UEqn.A());

// volVectorField HbyA("HbyA", U);
// HbyA = rAU * UEqn.H();
// surfaceScalarField phiHbyA(
//     "phiHbyA",
//     (fvc::interpolate(HbyA) & mesh.Sf())
//     // + fvc::ddtPhiCorr(rAU, U, phi)
// );

// phi = phiHbyA; // fvc::flux(U);

// UEqn.solve();
initialResidual = UEqn.solve().max().initialResidual();
U.correctBoundaryConditions();

phi = (fvc::interpolate(U) & mesh.Sf());
//rhoPhi = fvc::interpolate(rho) * phi;

rhoPhi = fvc::flux
(
    phi,
    rho,
    "div(phi,rho)"      // same scheme name as in fvSchemes
);


K = 0.5 * magSqr(U);

} // end if (stressUpdate)

// sigmatemp*=0.0;