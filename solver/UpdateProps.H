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
    volScalarField & ni = PopBal[h];

    GBREF+=ni*ni;

    PopBal_av[h] = fvc::average(ni); //takes the average from all cells sharing a face - this way the information about phases in adjacent cells is shared with the local cell
    PopBal_av[h].correctBoundaryConditions();

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
      ;
      // gamma_alpha_dots[i][countS].correctBoundaryConditions();
      // gamma_alphas[i][countS].correctBoundaryConditions();
      // g_alphas[i][countS].correctBoundaryConditions();
    }
  }



// // //   // // //update hardening rates
// // // // // // Info<<"HERE1"<<endl;

  forAll(mesh.C(), cell)
  {


        label active_phases[N_Ori] =
    {
        0
    };
    label n_active_phases = 0;



    for (label count = 0; count < N_Ori; count++)
    {
      
        if (PopBal_av[count][cell] > 0.0)
        {
            active_phases[n_active_phases++] = count;
        }

    }

    for (label h = 0; h < n_active_phases; h++)
    {

          double total_slip_on_phi = 0.0;
          for (int slip_it = 0; slip_it < MatList[Phase2Material[active_phases[h]]].N_Slips; slip_it++)
        {
          total_slip_on_phi+=mag(gamma_alphas[active_phases[h]][slip_it][cell]);
        }

        for (int S_alpha = 0; S_alpha < MatList[Phase2Material[active_phases[h]]].N_Slips; S_alpha++)
        {


          g_alpha_dots[active_phases[h]][S_alpha][cell] *= 0.0;

          for (int S_beta = 0; S_beta < MatList[Phase2Material[active_phases[h]]].N_Slips; S_beta++)
          {

            //   // double h_of_gamma = ho*Foam::pow(1.0/(Foam::cosh((ho*Acc_Strain[cell])/(ta_s-tau_0))),2.0);
            double h_alpha_beta = 0.0;
            // Pout<<Pstream::myProcNo()<<"test_strain: "<<Acc_Strain[cell]<<endl;

            double argument = (MatList[Phase2Material[active_phases[h]]].h_0 * total_slip_on_phi) / (MatList[Phase2Material[active_phases[h]]].tau_S - MatList[Phase2Material[active_phases[h]]].initial_CRSS);
            // argument = max(min(argument,20.0),-20.0);

            // double h_of_gamma = MatList[Phase2Material[h]].h_0 * Foam::pow(1.0 / (Foam::cosh((MatList[Phase2Material[h]].h_0 * gamma_alphas[h][S_alpha][cell]) / (MatList[Phase2Material[h]].tau_S - MatList[Phase2Material[h]].initial_CRSS))), 2.0);
              double h_of_gamma = MatList[Phase2Material[active_phases[h]]].h_0 * Foam::pow(1.0 / (Foam::cosh(argument)), 2.0);            

            if (S_beta == S_alpha)
            {
              h_alpha_beta = h_of_gamma * (MatList[Phase2Material[active_phases[h]]].r + (1.0 - MatList[Phase2Material[active_phases[h]]].r));
            }
            else
            {
              h_alpha_beta = h_of_gamma * (MatList[Phase2Material[active_phases[h]]].r);
            }

            g_alpha_dots[active_phases[h]][S_alpha][cell] += /*PopBal[h][cell]*/(h_alpha_beta)*gamma_alpha_dots[active_phases[h]][S_beta][cell];
          }
          // double deltaT_ =  runTime.deltaT().value(); 
          // g_alpha_dots[h][S_alpha][cell]=max(min(g_alpha_dots[h][S_alpha][cell],0.99/deltaT_),-0.99/deltaT_);
          //
        }

    // }

    }
  }


//Over Boundaries
forAll(T.boundaryField(), patchI)
{
    if
    (
        !T.boundaryField()[patchI].coupled()
     && T.boundaryField()[patchI].size() > 0
    )
    {
        forAll(T.boundaryField()[patchI], faceI)
        {






        label active_phases[N_Ori] =
    {
        0
    };
    label n_active_phases = 0;



    for (label count = 0; count < N_Ori; count++)
    {

        if (PopBal_av[count].boundaryField()[patchI][faceI] > 0.0)
        {
            active_phases[n_active_phases++] = count;
        }

    }





   for (label h = 0; h < n_active_phases; h++)
    {

          double total_slip_on_phi = 0.0;
          for (int slip_it = 0; slip_it < MatList[Phase2Material[active_phases[h]]].N_Slips; slip_it++)
        {
          total_slip_on_phi+=mag(gamma_alphas[active_phases[h]][slip_it].boundaryField()[patchI][faceI]);
        }

        for (int S_alpha = 0; S_alpha < MatList[Phase2Material[active_phases[h]]].N_Slips; S_alpha++)
        {


          g_alpha_dots[active_phases[h]][S_alpha].boundaryFieldRef()[patchI][faceI] *= 0.0;

          for (int S_beta = 0; S_beta < MatList[Phase2Material[active_phases[h]]].N_Slips; S_beta++)
          {


            double h_alpha_beta = 0.0;

            double argument = (MatList[Phase2Material[active_phases[h]]].h_0 * total_slip_on_phi) / (MatList[Phase2Material[active_phases[h]]].tau_S - MatList[Phase2Material[active_phases[h]]].initial_CRSS);
            // argument = max(min(argument,20.0),-20.0);

              double h_of_gamma = MatList[Phase2Material[active_phases[h]]].h_0 * Foam::pow(1.0 / (Foam::cosh(argument)), 2.0);            

            if (S_beta == S_alpha)
            {
              h_alpha_beta = h_of_gamma * (MatList[Phase2Material[active_phases[h]]].r + (1.0 - MatList[Phase2Material[active_phases[h]]].r));
            }
            else
            {
              h_alpha_beta = h_of_gamma * (MatList[Phase2Material[active_phases[h]]].r);
            }

            g_alpha_dots[active_phases[h]][S_alpha].boundaryFieldRef()[patchI][faceI] += /*PopBal[h].boundaryField()[patchI][faceI]*/(h_alpha_beta)*gamma_alpha_dots[active_phases[h]][S_beta].boundaryField()[patchI][faceI];
          }
          // double deltaT_ =  runTime.deltaT().value(); 
          // g_alpha_dots[h][S_alpha].boundaryFieldRef()[patchI][faceI]=max(min(g_alpha_dots[h][S_alpha].boundaryField()[patchI][faceI],0.99/deltaT_),-0.99/deltaT_);
          //
        }

    // }


    }



        }

    }

}


// //END update hardening rates
















  for (label i = 0; i < N_Ori; i++)
  {
    for (label countS = 0; countS < MatList[Phase2Material[i]].N_Slips; countS++)
    {
      // g_alpha_dots[i][countS].write();
      gamma_alpha_dots[i][countS]*=0.0;
    }
  }





// Euler_Angles*=0.0;

  forAll(mesh.C(), cell)
  {



        label active_phases[N_Ori] =
    {
        0
    };
    label n_active_phases = 0;



    for (label count = 0; count < N_Ori; count++)
    {
      
        if (PopBal_av[count][cell] > 0.0)
        {
            active_phases[n_active_phases++] = count;
        }

    }




    Fastor::Tensor<double, 3, 3> Fsigma(0.0); // CAST OPENFOAM VARIABLES INTO FASTOR VARIABLES
    Fastor::Tensor<double, 3, 3> Fepsilon_dot(0.0);
    // Fastor::Tensor<double, 3, 3> Fsigmaav(0.0);
    int count_index = 0;
    for (int i = 0; i < 3; i++)
    {
      for (int j = 0; j < 3; j++)
      {
        // Info<<i<<","<<j<<endl;
        Fepsilon_dot(i, j) = epsilon_dot[cell].component(count_index);
        // Fsigma(i, j) = sigtemp[cell].component(count_index);
        // Fsigmaav(i, j) = sigmaav[cell].component(count_index);

        count_index++;
      }
    }


    // Fastor::Tensor<double, 3, 3> sigAv(0.0);
    Fastor::Tensor<double, 3, 3> D_P_Average(0.0);

    label largest_field=-1;
    double current_largest_value = -1.0;

    // sumgamma[cell]*=0.0;
    for (label h = 0; h < n_active_phases; h++)
    {

      if(PopBal[active_phases[h]][cell]>current_largest_value){
        current_largest_value = PopBal[active_phases[h]][cell];
      
        largest_field = active_phases[h];
      }


          int count_index_siglocal = 0;
    for (int i = 0; i < 3; i++)
    {
      for (int j = 0; j < 3; j++)
      {
        // Info<<i<<","<<j<<endl;
        
        Fsigma(i, j) = sigtemp[active_phases[h]][cell].component(count_index_siglocal);
        

        count_index_siglocal++;
      }
    }



      Fastor::Tensor<double, 3, 3, 3, 3> Cprime(0.0); // global stiffness at this cell

      double phi_1_total = Rot1[active_phases[h]][cell]; // ROT_1[cell];// phi1[h] ;//+ omega_tot[cell].component(5);

      double phi_2_total = Rot2[active_phases[h]][cell]; // ROT_2[cell];// phi2[h] ;//+ omega_tot[cell].component(2);

      double phi_3_total = Rot3[active_phases[h]][cell]; // ROT_3[cell];// phi3[h] ;//+ omega_tot[cell].component(1);

    // if(PopBal[h][cell]>0.5){
      // Euler_Angles[cell].component(0) += PopBal[h][cell]*Rot1[h][cell];
      // Euler_Angles[cell].component(1) += PopBal[h][cell]*Rot2[h][cell];
      // Euler_Angles[cell].component(2) += PopBal[h][cell]*Rot3[h][cell];
    // }

      // vector Basis_1(1.0,0.0,0.0);
      // vector Basis_2(0.0,1.0,0.0);
      // vector Basis_3(0.0,0.0,1.0);


      Fastor::Tensor<double, 3, 3> Rot_matrix(0.0);


      Rot_matrix(0, 0) = Foam::cos(phi_2_total) * Foam::cos(phi_3_total);
      Rot_matrix(0, 1) = (Foam::sin(phi_1_total) * Foam::sin(phi_2_total) * Foam::cos(phi_3_total)) - (Foam::cos(phi_1_total) * Foam::sin(phi_3_total));
      Rot_matrix(0, 2) = (Foam::cos(phi_1_total) * Foam::sin(phi_2_total) * Foam::cos(phi_3_total)) + (Foam::sin(phi_1_total) * Foam::sin(phi_3_total));

      Rot_matrix(1, 0) = Foam::cos(phi_2_total) * Foam::sin(phi_3_total);
      Rot_matrix(1, 1) = (Foam::sin(phi_1_total) * Foam::sin(phi_2_total) * Foam::sin(phi_3_total)) + (Foam::cos(phi_1_total) * Foam::cos(phi_3_total));
      Rot_matrix(1, 2) = (Foam::cos(phi_1_total) * Foam::sin(phi_2_total) * Foam::sin(phi_3_total)) - (Foam::sin(phi_1_total) * Foam::cos(phi_3_total));

      Rot_matrix(2, 0) = -Foam::sin(phi_2_total);
      Rot_matrix(2, 1) = Foam::sin(phi_1_total) * Foam::cos(phi_2_total);
      Rot_matrix(2, 2) = Foam::cos(phi_1_total) * Foam::cos(phi_2_total);




        for (int i = 0; i < 3; i++)
        {
          for (int j = 0; j < 3; j++)
          {
            for (int k = 0; k < 3; k++)
            {
              for (int l = 0; l < 3; l++)
              {

                double value = 0.0;
                for (int m = 0; m < 3; m++)
                {
                  for (int n = 0; n < 3; n++)
                  {
                    for (int p = 0; p < 3; p++)
                    {
                      for (int q = 0; q < 3; q++)
                      {
                        // if (Cprime2[Phase2Material[h]](m, n, p, q) > SMALL)
                        // {
                          value += Rot_matrix(i, m) * Rot_matrix(j, n) * Rot_matrix(k, p) * Rot_matrix(l, q) * Cprime2[Phase2Material[active_phases[h]]](m, n, p, q); // C(m,n,p,q);
                        // }
                      }
                    }
                  }
                }

                // value*=PopBal[h][cell];

                Cprime(i, j, k, l) = value;
              }
            }
          }
        }

      // }//

        // If grain is there
        // Have the stiffness tensor now
        // need Dp and omegaP

        Fastor::Tensor<double, 3, 3> F_D_P(0.0);
        Fastor::Tensor<double, 3, 3> F_Omega_P(0.0);

        // double SUMGAMMA=0.0;
        for (int s = 0; s < MatList[Phase2Material[active_phases[h]]].N_Slips; s++)
        {

          Fastor::Tensor<double, 3> F_n_alpha(0.0);
          Fastor::Tensor<double, 3> F_m_alpha(0.0);
          Fastor::Tensor<double, 3, 3> P_alpha(0.0);
          Fastor::Tensor<double, 3, 3> omega_alpha(0.0);

          F_n_alpha(0) = MatList[Phase2Material[active_phases[h]]].SP_norm[s].component(0); // n_alpha[Phase2Material[h]][s].component(0);
          F_n_alpha(1) = MatList[Phase2Material[active_phases[h]]].SP_norm[s].component(1); // n_alpha[Phase2Material[h]][s].component(1);
          F_n_alpha(2) = MatList[Phase2Material[active_phases[h]]].SP_norm[s].component(2); // n_alpha[Phase2Material[h]][s].component(2);

          F_m_alpha(0) = MatList[Phase2Material[active_phases[h]]].SP_direction[s].component(0); // m_alpha[Phase2Material[h]][s].component(0);
          F_m_alpha(1) = MatList[Phase2Material[active_phases[h]]].SP_direction[s].component(1); // m_alpha[Phase2Material[h]][s].component(1);
          F_m_alpha(2) = MatList[Phase2Material[active_phases[h]]].SP_direction[s].component(2); // m_alpha[Phase2Material[h]][s].component(2);


          Fastor::Tensor<double, 3> F_n_alpha_Rot(0.0);
          Fastor::Tensor<double, 3> F_m_alpha_Rot(0.0);

          for (int i = 0; i < 3; i++)
          {
            for (int j = 0; j < 3; j++)
            { // Rotated slip information
              F_n_alpha_Rot(i) += Rot_matrix(i, j) * F_n_alpha(j);
              F_m_alpha_Rot(i) += Rot_matrix(i, j) * F_m_alpha(j);
            }
          }


          P_alpha = 0.5 * (outer(F_m_alpha_Rot, F_n_alpha_Rot) + outer(F_n_alpha_Rot, F_m_alpha_Rot));
          omega_alpha = 0.5 * (outer(F_m_alpha_Rot, F_n_alpha_Rot) - outer(F_n_alpha_Rot, F_m_alpha_Rot));

          Fastor::Tensor<double, 3, 3> Schmid_T(0.0);

          Schmid_T = outer(F_m_alpha_Rot, F_n_alpha_Rot);
          
          double tau_alpha = 0.0;
          for (int i = 0; i < 3; i++)
          {
            for (int j = 0; j < 3; j++)
            {
              // tau_alpha += Fsigma(i, j)* Schmid_T(i, j);//Fsigmaav(i,j) 
              tau_alpha += Fsigma(i, j)* P_alpha(i, j);
            }
          }

          double k_exponent = MatList[Phase2Material[active_phases[h]]].k_exponent;         // 10.0;//12.0
          double ref_shear_rate = MatList[Phase2Material[active_phases[h]]].ref_shear_rate; // 1e-3;
          

          gamma_alpha_dots[active_phases[h]][s][cell] = ref_shear_rate * Foam::pow((mag(tau_alpha / g_alphas[active_phases[h]][s][cell])), k_exponent) * sign(tau_alpha); // source term for RHS of g_alpha equations


          // sumgamma[cell]+= mag(gamma_alpha_dots[h][s][cell]);//for next time-step
          // SUMGAMMA+=mag(gamma_alpha_dots[h][s][cell]);//for next time-step;

          F_D_P += (P_alpha * gamma_alpha_dots[active_phases[h]][s][cell]);
          F_Omega_P += (omega_alpha * gamma_alpha_dots[active_phases[h]][s][cell]);
          // }

        } // end of loop over slip systems

        // sumgamma[cell]+=PopBal[h][cell]*SUMGAMMA;

        int OP_index = 0;
        for (int i = 0; i < 3; i++)
        {
          for (int j = 0; j < 3; j++)
          {
            Omega_P[active_phases[h]][cell].component(OP_index) = F_Omega_P(i, j);
            OP_index++;
          }
        }

        enum
        {
          I,J,K,L
        };
        Fastor::Tensor<double, 3, 3> sig = Fastor::einsum<Fastor::Index<I, J, K, L>, Fastor::Index<K, L> /*,Fastor::OIndex<I,J>*/>(Cprime, Fepsilon_dot - F_D_P); ////Double contraction
        // Fastor::Tensor<double, 3, 3> sig = 80e9*(Fepsilon_dot - F_D_P);
        // sigAv += (PopBal[h][cell]) * sig;///Sum_Phases[cell];
      // }
        // Omega_P_Average+=PopBal[h][cell]*F_Omega_P;
        D_P_Average += PopBal[active_phases[h]][cell] * F_D_P;

        // Fastor::Tensor<double, 3, 3> F_D_E(0.0);

        // F_D_E = Fepsilon_dot-F_D_P;

       label count_indexnew = 0;
      for (int i = 0; i < 3; i++)
      {
        for (int j = 0; j < 3; j++)
        {

          dsigmadt[active_phases[h]][cell].component(count_indexnew) = sig(i, j); // sig(i,j);
        count_indexnew++;
        }
      }

          dsigmadt[active_phases[h]][cell] -=((sigtemp[active_phases[h]][cell] & omega_dot[cell]) + (omega_dot[cell] & sigtemp[active_phases[h]][cell])) ;//-(tr(epsilon_dot)*sigma);;// 
          dsigmadt[active_phases[h]][cell] += ((sigtemp[active_phases[h]][cell] & Omega_P[active_phases[h]][cell]) - (Omega_P[active_phases[h]][cell] & sigtemp[active_phases[h]][cell]));
          // dsigmadt[cell] -= (PopBal[h][cell]) *(sigtemp[cell]*(F_D_E(0,0)+F_D_E(1,1)+F_D_E(2,2)));




    } // end loop over phase fields


      // Euler_Angles[cell].component(0) = Rot1[largest_field][cell];
      // Euler_Angles[cell].component(1) = Rot2[largest_field][cell];
      // Euler_Angles[cell].component(2) = Rot3[largest_field][cell];

    count_index = 0;
    for (int i = 0; i < 3; i++)
    {
      for (int j = 0; j < 3; j++)
      {

        // dsigmadt[cell].component(count_index) += sigAv(i, j); // sig(i,j);
        D_P[cell].component(count_index) = D_P_Average(i, j); // F_D_P(i,j);
        count_index++;
      }
    }
  }











//Over Boundaries
forAll(T.boundaryField(), patchI)
{
    if
    (
        !T.boundaryField()[patchI].coupled()
     && T.boundaryField()[patchI].size() > 0
    )
    {
        forAll(T.boundaryField()[patchI], faceI)
        {


        label active_phases[N_Ori] =
    {
        0
    };
    label n_active_phases = 0;



    for (label count = 0; count < N_Ori; count++)
    {

        if (PopBal_av[count].boundaryField()[patchI][faceI] > 0.0)
        {
            active_phases[n_active_phases++] = count;
        }

    }



    Fastor::Tensor<double, 3, 3> Fsigma(0.0); // CAST OPENFOAM VARIABLES INTO FASTOR VARIABLES
    Fastor::Tensor<double, 3, 3> Fepsilon_dot(0.0);
    // Fastor::Tensor<double, 3, 3> Fsigmaav(0.0);
    int count_index = 0;
    for (int i = 0; i < 3; i++)
    {
      for (int j = 0; j < 3; j++)
      {
        // Info<<i<<","<<j<<endl;
        Fepsilon_dot(i, j) = epsilon_dot.boundaryField()[patchI][faceI].component(count_index);
        // Fsigma(i, j) = sigtemp.boundaryField()[patchI][faceI].component(count_index);
        // Fsigmaav(i, j) = sigmaav[cell].component(count_index);

        count_index++;
      }
    }


    // Fastor::Tensor<double, 3, 3> sigAv(0.0);
    Fastor::Tensor<double, 3, 3> D_P_Average(0.0);


        label largest_field=-1;
    double current_largest_value = -1.0;

    // sumgamma.boundaryFieldRef()[patchI][faceI]*=0.0;
    for (label h = 0; h < n_active_phases; h++)
    {


        if(PopBal[active_phases[h]].boundaryField()[patchI][faceI]>current_largest_value){
        current_largest_value = PopBal[active_phases[h]].boundaryField()[patchI][faceI];
      
        largest_field = active_phases[h];
      }


int count_index = 0;
    for (int i = 0; i < 3; i++)
    {
      for (int j = 0; j < 3; j++)
      {
       
        Fsigma(i, j) = sigtemp[active_phases[h]].boundaryField()[patchI][faceI].component(count_index);
        // Fsigmaav(i, j) = sigmaav[cell].component(count_index);

        count_index++;
      }
    }


       Fastor::Tensor<double, 3, 3, 3, 3> Cprime(0.0); // global stiffness at this cell

      double phi_1_total = Rot1[active_phases[h]].boundaryField()[patchI][faceI]; // ROT_1[cell];// phi1[h] ;//+ omega_tot[cell].component(5);

      double phi_2_total = Rot2[active_phases[h]].boundaryField()[patchI][faceI]; // ROT_2[cell];// phi2[h] ;//+ omega_tot[cell].component(2);

      double phi_3_total = Rot3[active_phases[h]].boundaryField()[patchI][faceI]; // ROT_3[cell];// phi3[h] ;//+ omega_tot[cell].component(1);

          // if(PopBal[h].boundaryField()[patchI][faceI]>0.5){

    // }

    Fastor::Tensor<double, 3, 3> Rot_matrix(0.0);

      Rot_matrix(0, 0) = Foam::cos(phi_2_total) * Foam::cos(phi_3_total);
      Rot_matrix(0, 1) = (Foam::sin(phi_1_total) * Foam::sin(phi_2_total) * Foam::cos(phi_3_total)) - (Foam::cos(phi_1_total) * Foam::sin(phi_3_total));
      Rot_matrix(0, 2) = (Foam::cos(phi_1_total) * Foam::sin(phi_2_total) * Foam::cos(phi_3_total)) + (Foam::sin(phi_1_total) * Foam::sin(phi_3_total));

      Rot_matrix(1, 0) = Foam::cos(phi_2_total) * Foam::sin(phi_3_total);
      Rot_matrix(1, 1) = (Foam::sin(phi_1_total) * Foam::sin(phi_2_total) * Foam::sin(phi_3_total)) + (Foam::cos(phi_1_total) * Foam::cos(phi_3_total));
      Rot_matrix(1, 2) = (Foam::cos(phi_1_total) * Foam::sin(phi_2_total) * Foam::sin(phi_3_total)) - (Foam::sin(phi_1_total) * Foam::cos(phi_3_total));

      Rot_matrix(2, 0) = -Foam::sin(phi_2_total);
      Rot_matrix(2, 1) = Foam::sin(phi_1_total) * Foam::cos(phi_2_total);
      Rot_matrix(2, 2) = Foam::cos(phi_1_total) * Foam::cos(phi_2_total);
      // https://en.wikipedia.org/wiki/Rotation_matrix -
      //"represents an extrinsic rotation whose (improper) Euler angles are α, β, γ, about axes x, y, z" - phi1 ==alpha, phi2==beta,phi3==gamma




        for (int i = 0; i < 3; i++)
        {
          for (int j = 0; j < 3; j++)
          {
            for (int k = 0; k < 3; k++)
            {
              for (int l = 0; l < 3; l++)
              {

                double value = 0.0;
                for (int m = 0; m < 3; m++)
                {
                  for (int n = 0; n < 3; n++)
                  {
                    for (int p = 0; p < 3; p++)
                    {
                      for (int q = 0; q < 3; q++)
                      {
                        // if (Cprime2[Phase2Material[h]](m, n, p, q) > SMALL)
                        // {
                          value += Rot_matrix(i, m) * Rot_matrix(j, n) * Rot_matrix(k, p) * Rot_matrix(l, q) * Cprime2[Phase2Material[active_phases[h]]](m, n, p, q); // C(m,n,p,q);
                        // }
                      }
                    }
                  }
                }

                // value*=PopBal[h][cell];

                Cprime(i, j, k, l) = value;
              }
            }
          }
        }

      // }//


// If grain is there
        // Have the stiffness tensor now
        // need Dp and omegaP

        Fastor::Tensor<double, 3, 3> F_D_P(0.0);
        Fastor::Tensor<double, 3, 3> F_Omega_P(0.0);

      // double SUMGAMMA=0.0;
      for (int s = 0; s < MatList[Phase2Material[active_phases[h]]].N_Slips; s++)
        {

          Fastor::Tensor<double, 3> F_n_alpha(0.0);
          Fastor::Tensor<double, 3> F_m_alpha(0.0);
          Fastor::Tensor<double, 3, 3> P_alpha(0.0);
          Fastor::Tensor<double, 3, 3> omega_alpha(0.0);

          F_n_alpha(0) = MatList[Phase2Material[active_phases[h]]].SP_norm[s].component(0); // n_alpha[Phase2Material[h]][s].component(0);
          F_n_alpha(1) = MatList[Phase2Material[active_phases[h]]].SP_norm[s].component(1); // n_alpha[Phase2Material[h]][s].component(1);
          F_n_alpha(2) = MatList[Phase2Material[active_phases[h]]].SP_norm[s].component(2); // n_alpha[Phase2Material[h]][s].component(2);

          F_m_alpha(0) = MatList[Phase2Material[active_phases[h]]].SP_direction[s].component(0); // m_alpha[Phase2Material[h]][s].component(0);
          F_m_alpha(1) = MatList[Phase2Material[active_phases[h]]].SP_direction[s].component(1); // m_alpha[Phase2Material[h]][s].component(1);
          F_m_alpha(2) = MatList[Phase2Material[active_phases[h]]].SP_direction[s].component(2); // m_alpha[Phase2Material[h]][s].component(2);


          Fastor::Tensor<double, 3> F_n_alpha_Rot(0.0);
          Fastor::Tensor<double, 3> F_m_alpha_Rot(0.0);

          for (int i = 0; i < 3; i++)
          {
            for (int j = 0; j < 3; j++)
            { // Rotated slip information
              F_n_alpha_Rot(i) += Rot_matrix(i, j) * F_n_alpha(j);
              F_m_alpha_Rot(i) += Rot_matrix(i, j) * F_m_alpha(j);
            }
          }


          P_alpha = 0.5 * (outer(F_m_alpha_Rot, F_n_alpha_Rot) + outer(F_n_alpha_Rot, F_m_alpha_Rot));
          omega_alpha = 0.5 * (outer(F_m_alpha_Rot, F_n_alpha_Rot) - outer(F_n_alpha_Rot, F_m_alpha_Rot));

          Fastor::Tensor<double, 3, 3> Schmid_T(0.0);

          Schmid_T = outer(F_m_alpha_Rot, F_n_alpha_Rot);
          double tau_alpha = 0.0;
          for (int i = 0; i < 3; i++)
          {
            for (int j = 0; j < 3; j++)
            {
              // tau_alpha += Fsigma(i, j)* Schmid_T(i, j);//Fsigmaav(i,j) 
              tau_alpha += Fsigma(i, j)* P_alpha(i, j);//Fsigmaav(i,j) 
            }
          }

          double k_exponent = MatList[Phase2Material[active_phases[h]]].k_exponent;         // 10.0;//12.0
          double ref_shear_rate = MatList[Phase2Material[active_phases[h]]].ref_shear_rate; // 1e-3;

          gamma_alpha_dots[active_phases[h]][s].boundaryFieldRef()[patchI][faceI] = ref_shear_rate * Foam::pow((mag(tau_alpha / g_alphas[active_phases[h]][s].boundaryField()[patchI][faceI])), k_exponent) * sign(tau_alpha); // source term for RHS of g_alpha equations


          // sumgamma.boundaryFieldRef()[patchI][faceI]+= mag(gamma_alpha_dots[h][s].boundaryField()[patchI][faceI]);//for next time-step

          // SUMGAMMA+= mag(gamma_alpha_dots[h][s].boundaryField()[patchI][faceI]);//for next time-step

          F_D_P += (P_alpha * gamma_alpha_dots[active_phases[h]][s].boundaryField()[patchI][faceI]);
          F_Omega_P += (omega_alpha * gamma_alpha_dots[active_phases[h]][s].boundaryField()[patchI][faceI]);
          // }

        } // end of loop over slip systems
        // sumgamma.boundaryFieldRef()[patchI][faceI]+=PopBal[h].boundaryField()[patchI][faceI]*SUMGAMMA;

        int OP_index = 0;
        for (int i = 0; i < 3; i++)
        {
          for (int j = 0; j < 3; j++)
          {
            Omega_P[active_phases[h]].boundaryFieldRef()[patchI][faceI].component(OP_index) = F_Omega_P(i, j); //Philip, how to assign the boundary face values
            OP_index++;
          }
        }

        enum
        {
          I,J,K,L
        };
        Fastor::Tensor<double, 3, 3> sig = Fastor::einsum<Fastor::Index<I, J, K, L>, Fastor::Index<K, L> /*,Fastor::OIndex<I,J>*/>(Cprime, Fepsilon_dot - F_D_P); ////Double contraction
        // Fastor::Tensor<double, 3, 3> sig = 80e9*(Fepsilon_dot - F_D_P);

        // sigAv += (PopBal[h].boundaryField()[patchI][faceI]) * sig;///Sum_Phases[cell];
        // Omega_P_Average+=PopBal[h][cell]*F_Omega_P;
      // }
        D_P_Average += PopBal[active_phases[h]].boundaryField()[patchI][faceI]  * F_D_P;

        // Fastor::Tensor<double, 3, 3> F_D_E(0.0);

        // F_D_E = Fepsilon_dot-F_D_P;



          label  count_indexnew = 0;
    for (int i = 0; i < 3; i++)
    {
      for (int j = 0; j < 3; j++)
      {

        dsigmadt[active_phases[h]].boundaryFieldRef()[patchI][faceI].component(count_indexnew) = sig(i, j); // sig(i,j);
        // D_P.boundaryFieldRef()[patchI][faceI].component(count_indexnew) = D_P_Average(i, j); // F_D_P(i,j);
        count_indexnew++;
      }
    }

          dsigmadt[active_phases[h]].boundaryFieldRef()[patchI][faceI] -=((sigtemp[active_phases[h]].boundaryField()[patchI][faceI] & omega_dot.boundaryField()[patchI][faceI]) + (omega_dot.boundaryField()[patchI][faceI] & sigtemp[active_phases[h]].boundaryField()[patchI][faceI]));
          dsigmadt[active_phases[h]].boundaryFieldRef()[patchI][faceI] += ((sigtemp[active_phases[h]].boundaryField()[patchI][faceI] & Omega_P[active_phases[h]].boundaryField()[patchI][faceI]) - (Omega_P[active_phases[h]].boundaryField()[patchI][faceI] & sigtemp[active_phases[h]].boundaryField()[patchI][faceI]));//Philip, can we assign boundary fields with dot products in?
          // dsigmadt.boundaryFieldRef()[patchI][faceI] -= (PopBal[h].boundaryField()[patchI][faceI]) *(sigtemp.boundaryField()[patchI][faceI]*(F_D_E(0,0)+F_D_E(1,1)+F_D_E(2,2)));



    }// end loop over phase fields

      // Euler_Angles.boundaryFieldRef()[patchI][faceI].component(0) = Rot1[largest_field].boundaryField()[patchI][faceI];
      // Euler_Angles.boundaryFieldRef()[patchI][faceI].component(1) = Rot2[largest_field].boundaryField()[patchI][faceI];
      // Euler_Angles.boundaryFieldRef()[patchI][faceI].component(2) = Rot3[largest_field].boundaryField()[patchI][faceI];


    count_index = 0;
    for (int i = 0; i < 3; i++)
    {
      for (int j = 0; j < 3; j++)
      {

        // dsigmadt.boundaryFieldRef()[patchI][faceI].component(count_index) += sigAv(i, j); // sig(i,j);
        D_P.boundaryFieldRef()[patchI][faceI].component(count_index) = D_P_Average(i, j); // F_D_P(i,j);
        count_index++;
      }
    }



        
        
        }
    }
}

  




  for (label i = 0; i < N_Ori; i++)
  {
    for (label countS = 0; countS < MatList[Phase2Material[i]].N_Slips; countS++)
    {
      g_alpha_dots[i][countS].correctBoundaryConditions();
      gamma_alpha_dots[i][countS].correctBoundaryConditions();
      // gamma_alphas[i][countS].correctBoundaryConditions();
      g_alphas[i][countS].correctBoundaryConditions();
    }
  }


 for (label i = 0; i < N_Ori; i++)
  {
dsigmadt[i].correctBoundaryConditions();
  }





}
