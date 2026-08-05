
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