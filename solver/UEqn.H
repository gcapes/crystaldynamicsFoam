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

// surfaceTensorField sigmaf("sigmaf", fvc::interpolate(sigma));

divsigma = fvc::div(sigma);

// divsigma = fvc::div(fvc::interpolate(sigma) & mesh.Sf());



// sigma_col1.component(0) = sigma.component(0);
// sigma_col1.component(1) = sigma.component(3);
// sigma_col1.component(2) = sigma.component(6);

// sigma_col2.component(0) = sigma.component(1);
// sigma_col2.component(1) = sigma.component(4);
// sigma_col2.component(2) = sigma.component(7);

// sigma_col3.component(0) = sigma.component(2);
// sigma_col3.component(1) = sigma.component(5);
// sigma_col3.component(2) = sigma.component(8);

// sigma_col1.correctBoundaryConditions();
// sigma_col2.correctBoundaryConditions();
// sigma_col3.correctBoundaryConditions();

// surfaceScalarField c1f = (fvc::flux(sigma_col1));
// surfaceScalarField c2f = (fvc::flux(sigma_col2));
// surfaceScalarField c3f = (fvc::flux(sigma_col3));

// divsigma.component(0) = fvc::div(c1f);
// divsigma.component(1) = fvc::div(c2f);
// divsigma.component(2) = fvc::div(c3f);



// damping = ;

// TESTADV = fvc::div(rho * (U*U));

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
