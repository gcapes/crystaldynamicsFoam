{

// volTensorField sigmaold(sigtemp);

// dsigmadt +=(U & fvc::grad(sigma));
// dsigmadt -=fvc::div(phi,sigma);

  // sigma = sigma + (1.0 / n_sub_steps)*  runTime.deltaT() * (dsigmadt);
  // sigma.correctBoundaryConditions();

// fvc::div(phi,sigma);

  sigtemp = sigma +runTime.deltaT() * (dsigmadt);//+fvc::div(phi,sigma)
  sigtemp.correctBoundaryConditions();
// D = D + (1.0 / n_sub_steps) * runTime.deltaT() * (U);
// D.correctBoundaryConditions();

p=-(1.0/3.0)*tr(sigma);

testtraceeps = tr(epsilon_dot-D_P)*sigma;

// sigmadiff = sigmaold - sigtemp;

if (pimple.finalPimpleIter())
{
Info<<"Stress Field Updated"<<endl;
sigma=sigtemp;

D = D +  runTime.deltaT() * (U);
D.correctBoundaryConditions();


}
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