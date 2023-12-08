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



forAll(mesh.C(), cell)
  {

    // Info<<"One "<<epsilon_dot[cell]<<endl;

    dsigmadt[cell]+=2.0*shear_mod.value()*epsilon_dot[cell];
    // epsilon_dot2[cell]=epsilon_dot[cell]
  // Info<<"Two "<<epsilon_dot[cell]<<endl;
  }
  
forAll(mesh.boundary(), patchI)
{
    if (!mesh.boundary()[patchI].coupled() && mesh.boundary()[patchI].size() > 0)
    {
        forAll(mesh.boundary()[patchI], faceI)
                 {
            dsigmadt.boundaryFieldRef[patchI][faceI]+=
                shear_mod.value()*epsilon_dot.boundaryFieldRef[patchI][faceI];
                 }
         }
}
dsigmadt.correctBoundaryConditions();






  sigma = sigma + (1.0 / n_sub_steps) * runTime.deltaT() * (dsigmadt);
  sigma.correctBoundaryConditions();



}