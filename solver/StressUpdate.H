{


  sigma = sigma + (1.0 / n_sub_steps)*  runTime.deltaT() * (dsigmadt);
  sigma.correctBoundaryConditions();

D = D + (1.0 / n_sub_steps) * runTime.deltaT() * (U);
D.correctBoundaryConditions();

p=-(1.0/3.0)*tr(sigma);

testtraceeps = tr(epsilon_dot-D_P)*sigma;


/*

    fvTensorMatrix  sigmaEqn
    (
        fvm::ddt(sigma)
        ==
        (dsigmadt)

    );

    sigmaEqn.relax();

    sigmaEqn.solve();
    sigma.correctBoundaryConditions();

*/


}