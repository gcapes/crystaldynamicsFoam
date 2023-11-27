// MRF.correctBoundaryVelocity(U);

// fvc::interpolate(s);

// surfaceTensorField sf("sf", fvc::interpolate(s));
// // volVectorField YT =
// fvc::div(sf);

// Info<<"TEST"<<fvc::div(sf)<<endl;

// fvTensorMatrix  sigmaEqn
// (
//     fvm::ddt(sigma)
//     ==
//     dsigmadt

// );

// sigmaEqn.relax();

// sigmaEqn.solve();

divsigma = fvc::div(sigma);

// damping = ;

// TESTADV = fvc::div(rhoPhi, U);

fvVectorMatrix UEqn(
    fvm::ddt(rho, U) 
    + fvm::div(rhoPhi, U) 
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

UEqn.solve();
U.correctBoundaryConditions();


rhoPhi = fvc::interpolate(rho) * phi; // fvc::flux(U);
K = 0.5 * magSqr(U);
