

fvVectorMatrix UEqn(
    fvm::ddt(rho, U) 
    + fvm::div(rhoPhi, U) 
    // +fvc::div(rho * (U*U))
    - fvc::div(sigma)
    // ==

);

UEqn.relax();



volScalarField rAU(1.0 / UEqn.A());

volVectorField HbyA("HbyA", U);
HbyA = rAU * UEqn.H();
surfaceScalarField phiHbyA(
    "phiHbyA",
    (fvc::interpolate(HbyA) & mesh.Sf())
    // + fvc::ddtPhiCorr(rAU, U, phi)
);

phi = phiHbyA; // fvc::flux(U);

// UEqn.solve();
initialResidual = UEqn.solve().max().initialResidual();
U.correctBoundaryConditions();


rhoPhi = fvc::interpolate(rho) * phi; // fvc::flux(U);
K = 0.5 * magSqr(U);
