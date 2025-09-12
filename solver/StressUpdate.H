{



sigmaGlobaltemp*=0.0;

for (label h = 0; h < N_Ori; h++)
            {
  sigtemp[h] = sigma[h] +runTime.deltaT() * (dsigmadt[h]);//+fvc::div(phi,sigma)
  sigtemp[h].correctBoundaryConditions();

  sigmaGlobaltemp+=PopBal[h]*sigtemp[h];
            }
            sigmaGlobaltemp.correctBoundaryConditions();

// p=-(1.0/3.0)*tr(sigma);

// testtraceeps = tr(epsilon_dot-D_P)*sigma;

// sigmadiff = sigmaold - sigtemp;

if (pimple.finalPimpleIter())
{
Info<<"Stress Field Updated"<<endl;
for (label h = 0; h < N_Ori; h++)
            {
              sigma[h]=sigtemp[h];
            }
            sigmaGlobal = sigmaGlobaltemp;
            sigmaGlobal.correctBoundaryConditions();

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