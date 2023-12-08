{


  volTensorField TgradU(fvc::grad(U)().T());
  // TgradU.correctBoundaryConditions(); // already called by line above

  omega_dot = 0.5 * (fvc::grad(U) - TgradU);   // this is the "symm" function, i.e. symm(gradU)
  // omega_dot.correctBoundaryConditions();       // already called by line above
  epsilon_dot = 0.5 * (fvc::grad(U) + TgradU);//- (1.0/3.0)*(fvc::div(phi))*I; // this looks like dev(gradU)
  // epsilon_dot.correctBoundaryConditions();     // already called by line above

  dsigmadt = -(sigma & omega_dot) + (omega_dot & sigma);// - fvc::div(phi,sigma); // -(tr(epsilon_dot)*sigma);



 dimensionedScalar shear_mod ("shear_mod", dimensionSet(1,-1,-2,0,0,0,0),80e9);


//  dsigmadt+=(2.0*shear_mod*epsilon_dot);//this gives what we want
// epsilon_dot.write();

// Info<<"ONE"<<epsilon_dot<<endl;
// epsilon_dot.correctBoundaryConditions();

// dsigmadt.write();



forAll(dsigmadt, cellI)
{
    // Info<<"One "<<epsilon_dot[cellI]<<endl;
    dsigmadt[cellI]+=2.0*shear_mod.value()*epsilon_dot[cellI];
    // epsilon_dot2[cellI]=epsilon_dot[cellI]
    // Info<<"Two "<<epsilon_dot[cellI]<<endl;
}

forAll(dsigmadt.boundaryField(), patchI)
{
    if
    (
        !dsigmadt.boundaryField()[patchI].coupled()
     && dsigmadt.boundaryField()[patchI].size() > 0
    )
    {
        forAll(dsigmadt.boundaryField()[patchI], faceI)
        {
            dsigmadt.boundaryFieldRef()[patchI][faceI] +=
                2.0*shear_mod.value()*epsilon_dot.boundaryField()[patchI][faceI];
        }
    }
}

dsigmadt.correctBoundaryConditions();






  sigma = sigma + (1.0 / n_sub_steps) * runTime.deltaT() * (dsigmadt);
  sigma.correctBoundaryConditions();



}
