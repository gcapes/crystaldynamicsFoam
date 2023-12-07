{

  // Info<<n_sub_steps<<endl;

  // divU=fvc::div(U);

  dimensionedScalar shear_mod ("shear_mod", dimensionSet(1,-1,-2,0,0,0,0),80e9);

  volTensorField TgradU(fvc::grad(U)().T());
  TgradU.correctBoundaryConditions(); // already called by line above

  omega_dot = 0.5 * (fvc::grad(U) - TgradU);   // this is the "symm" function, i.e. symm(gradU)
  omega_dot.correctBoundaryConditions();       // already called by line above
  epsilon_dot = 0.5 * (fvc::grad(U) + TgradU);//- (1.0/3.0)*(fvc::div(phi))*I; // this looks like dev(gradU)
  epsilon_dot.correctBoundaryConditions();     // already called by line above

  dsigmadt =(2.0*shear_mod*epsilon_dot) -(sigma & omega_dot) + (omega_dot & sigma);// - fvc::div(phi,sigma); // -(tr(epsilon_dot)*sigma);
  dsigmadt.correctBoundaryConditions();
  
  
  sigma = sigma + (1.0 / n_sub_steps) * runTime.deltaT() * (dsigmadt);
  sigma.correctBoundaryConditions();


}

