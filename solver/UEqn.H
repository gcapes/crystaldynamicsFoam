   // MRF.correctBoundaryVelocity(U);

    // fvc::interpolate(s);

    // surfaceTensorField sf("sf", fvc::interpolate(s));
    // // volVectorField YT = 
    // fvc::div(sf);

    // Info<<"TEST"<<fvc::div(sf)<<endl;

// surfaceScalarField rhoUf("rhoUf", fvc::interpolate(rho*U));

    divsigma = fvc::div(sigma);

    // damping = ;


    fvVectorMatrix UEqn
    (
        fvm::ddt(rho, U) 
        + fvm::div(rhoPhi, U)
      -fvc::div(sigma)
    // ==

    );

    UEqn.relax();


    UEqn.solve();
    U.correctBoundaryConditions();

    K = 0.5*magSqr(U);

