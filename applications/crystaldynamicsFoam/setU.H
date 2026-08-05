{

  forAll(mesh.C(), cell)
  {
    U[cell].component(1)=M_PI*10000.0*(mesh.C()[cell].component(2)-0.1e-3);
    U[cell].component(2)=M_PI*10000.0*(-mesh.C()[cell].component(1)+0.1e-3);
  }

U.correctBoundaryConditions();

phi = (fvc::interpolate(U) & mesh.Sf());

rhoPhi = fvc::interpolate(rho) * phi; // fvc::flux(U);

}
