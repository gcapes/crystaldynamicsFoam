// Momentum solve is bypassed when stressUpdate=false (prescribed-stress mode).
// U, phi, and rhoPhi are held at their current values (zero for a fresh start),
// preventing the prescribed sigma from driving spurious acceleration.
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
