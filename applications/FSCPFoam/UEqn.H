

dimensionedScalar D_eff = D_stab * runTime.deltaT();

fvVectorMatrix UEqn
(
    fvm::ddt(rho, U) 
    + fvm::div(rhoPhi, U)
    - fvm::laplacian(D_eff, U)                     // Implicit stabilization
    ==
      fvc::div(sigmaGlobaltemp, "div(sigma)")       // Stress divergence 
    - fvc::laplacian(D_eff, U)                      // Cancels implicit term over many iterations
);

UEqn.relax();

initialResidual = UEqn.solve().max().initialResidual();
U.correctBoundaryConditions();


phi = (fvc::interpolate(U) & mesh.Sf());
rhoPhi = fvc::interpolate(rho) * phi;

K = 0.5 * magSqr(U);