{





  // gradU = fvc::grad(U);
  // gradU.correctBoundaryConditions();

  volTensorField TgradU(fvc::grad(U)().T());
  TgradU.correctBoundaryConditions(); // already called by line above

  omega_dot = 0.5 * (TgradU - fvc::grad(U));   
omega_dot.correctBoundaryConditions();
  //NOTE!!! Openfoam defines grad in the opposite manner to solid mechanics textbooks, it's actually grad transpose that we want so this has to be transposed again or just ahev a - sign

  epsilon_dot = 0.5 * (TgradU + fvc::grad(U));//- (1.0/3.0)*(fvc::div(phi))*I; // this looks like dev(gradU)
epsilon_dot.correctBoundaryConditions();


// dsigmadt*=0.0;



//
// gradU = TgradU;//fvc::grad(U);
// L = epsilon_dot + omega_dot;
//
kappa*=0.0;
Cv*=0.0;
GBREF*=0.0;
  for (label h = 0; h < N_Ori; h++)
  {
    dsigmadt[h]*=0.0;
    Omega_P[h] *= 0.0;
    Plastic_Strain[h]*=0.0;
    volScalarField & ni = PopBal[h];

    GBREF+=ni*ni;

    PopBal_av[h] = fvc::average(ni); //takes the average from all cells sharing a face - this way the information about phases in adjacent cells is shared with the local cell
    PopBal_av[h].correctBoundaryConditions();

    dpdt[h] *= 0.0;
    dpdt[h].correctBoundaryConditions();

    Laplacian_n[h] = fvc::laplacian(ni);
    Laplacian_n[h].correctBoundaryConditions();






  dimensionedScalar dimkappa("dimkappa", dimensionSet(1, 1, -3, -1, 0), 1.0);
  dimensionedScalar dimcv("dimcv", dimensionSet(0, 2, -2, -1, 0), 1.0);



kappa+=PopBal[h]*MatList[Phase2Material[h]].kappa*dimkappa;

Cv+=PopBal[h]*MatList[Phase2Material[h]].Cv*dimcv;


  }

  GBREF.correctBoundaryConditions();

  D_P *= 0.0;

  p *= 0.0;

  // for (label i = 0; i < N_Mats; i++)
  // {
  //   for (label countS = 0; countS < MatList[i].N_Slips; countS++)
  //   {
  //     g_alpha_dots[i][countS] *= 0.0;
  //     // gamma_alpha_dots[i][countS]*=0.0;
  //   }
  // }


    for (label i = 0; i < N_Ori; i++)
  {
    for (label countS = 0; countS < MatList[Phase2Material[i]].N_Slips; countS++)
    {
      g_alpha_dots[i][countS] *= 0.0;
      gamma_alpha_dots[i][countS]*= 0.0;
      // gamma_alpha_dots[i][countS].correctBoundaryConditions();
      // gamma_alphas[i][countS].correctBoundaryConditions();
      // g_alphas[i][countS].correctBoundaryConditions();
    }
  }






}
