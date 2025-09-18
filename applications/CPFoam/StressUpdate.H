{



sigmaGlobaltemp*=0.0;

for (label h = 0; h < N_Ori; h++)
            {

        if (pimpleIter == 1) {  // First PIMPLE iteration
        Info<<"initialising stress for sub-stepping"<<endl;
            sigtemp[h] = sigma[h];  // Start from committed stress
        }

              // int n_substeps = 10;
              

              Info<<"Number of Stress Integration Sub-Steps: "<<n_substeps<<endl;
              dimensionedScalar sub_dt = runTime.deltaT() / n_substeps;

volTensorField sigma_working = sigtemp[h];  // Start from current trial state
        

        // Sub-stepping loop
        for (int step = 0; step < n_substeps; step++) {
          // Info<<"Stress-Update Sub-Step: "<<step<<endl;

          sigma_working += sub_dt * dsigmadt[h];

          sigma_working -= sub_dt * ((sigma_working & omega_dot) + 
                                     (omega_dot & sigma_working));
            sigma_working += sub_dt * ((sigma_working & Omega_P[h]) - 
                                     (Omega_P[h] & sigma_working));

            sigma_working -= sub_dt *   sigma_working*tr(epsilon_dot-Plastic_Strain[h]);                      



        if (step == n_substeps - 1) {  // Last sub-step


    // for (label countS = 0; countS < MatList[Phase2Material[h]].N_Slips; countS++)//reset hardening rates
    // {
    //   g_alpha_dots[h][countS] *= 0.0;
    // }

    // // Recompute plastic rates with updated stress
    // forAll(mesh.C(), cell) {
    // //     // Recompute gamma_alpha_dots using updated sigma_working
    // //     // Copy the slip rate computation logic from UpdateProps.H here
    // //     // but use sigma_working instead of old sigtemp[h]




    //       double total_slip_on_phi = 0.0;
    //       for (int slip_it = 0; slip_it < MatList[Phase2Material[h]].N_Slips; slip_it++)
    //     {
    //       total_slip_on_phi+=mag(gamma_alphas[h][slip_it][cell]);
    //     }

    //     for (int S_alpha = 0; S_alpha < MatList[Phase2Material[h]].N_Slips; S_alpha++)
    //     {


    //       g_alpha_dots[h][S_alpha][cell] *= 0.0;

    //       for (int S_beta = 0; S_beta < MatList[Phase2Material[h]].N_Slips; S_beta++)
    //       {


    //         double h_alpha_beta = 0.0;
    //         double argument = (MatList[Phase2Material[h]].h_0 * total_slip_on_phi) / (MatList[Phase2Material[h]].tau_S - MatList[Phase2Material[h]].initial_CRSS);
    //         double h_of_gamma = MatList[Phase2Material[h]].h_0 * Foam::pow(1.0 / (Foam::cosh(argument)), 2.0);            

    //         if (S_beta == S_alpha)
    //         {
    //           h_alpha_beta = h_of_gamma * (MatList[Phase2Material[h]].r + (1.0 - MatList[Phase2Material[h]].r));
    //         }
    //         else
    //         {
    //           h_alpha_beta = h_of_gamma * (MatList[Phase2Material[h]].r);
    //         }

    //         g_alpha_dots[h][S_alpha][cell] += /*PopBal[h][cell]*/(h_alpha_beta)*mag(gamma_alpha_dots[h][S_beta][cell]);
    //       }

    //     }








    // Fastor::Tensor<double, 3, 3> Fsigma(0.0); // CAST OPENFOAM VARIABLES INTO FASTOR VARIABLES
    //           int count_index_siglocal = 0;
    // for (int i = 0; i < 3; i++)
    // {
    //   for (int j = 0; j < 3; j++)
    //   {
    //     // Info<<i<<","<<j<<endl;
        
    //     Fsigma(i, j) = sigma_working[cell].component(count_index_siglocal);
        

    //     count_index_siglocal++;
    //   }
    // }

    //      Fastor::Tensor<double, 3, 3, 3, 3> Cprime(0.0); // global stiffness at this cell

    //   double phi_1_total = Rot1[h][cell]; // ROT_1[cell];// phi1[h] ;//+ omega_tot[cell].component(5);

    //   double phi_2_total = Rot2[h][cell]; // ROT_2[cell];// phi2[h] ;//+ omega_tot[cell].component(2);

    //   double phi_3_total = Rot3[h][cell]; // ROT_3[cell];// phi3[h] ;//+ omega_tot[cell].component(1);


    //   Fastor::Tensor<double, 3, 3> Rot_matrix(0.0);


    //   Rot_matrix(0, 0) = Foam::cos(phi_2_total) * Foam::cos(phi_3_total);
    //   Rot_matrix(0, 1) = (Foam::sin(phi_1_total) * Foam::sin(phi_2_total) * Foam::cos(phi_3_total)) - (Foam::cos(phi_1_total) * Foam::sin(phi_3_total));
    //   Rot_matrix(0, 2) = (Foam::cos(phi_1_total) * Foam::sin(phi_2_total) * Foam::cos(phi_3_total)) + (Foam::sin(phi_1_total) * Foam::sin(phi_3_total));

    //   Rot_matrix(1, 0) = Foam::cos(phi_2_total) * Foam::sin(phi_3_total);
    //   Rot_matrix(1, 1) = (Foam::sin(phi_1_total) * Foam::sin(phi_2_total) * Foam::sin(phi_3_total)) + (Foam::cos(phi_1_total) * Foam::cos(phi_3_total));
    //   Rot_matrix(1, 2) = (Foam::cos(phi_1_total) * Foam::sin(phi_2_total) * Foam::sin(phi_3_total)) - (Foam::sin(phi_1_total) * Foam::cos(phi_3_total));

    //   Rot_matrix(2, 0) = -Foam::sin(phi_2_total);
    //   Rot_matrix(2, 1) = Foam::sin(phi_1_total) * Foam::cos(phi_2_total);
    //   Rot_matrix(2, 2) = Foam::cos(phi_1_total) * Foam::cos(phi_2_total);





    // for (int s = 0; s < MatList[Phase2Material[h]].N_Slips; s++)
    //     {

    //       Fastor::Tensor<double, 3> F_n_alpha(0.0);
    //       Fastor::Tensor<double, 3> F_m_alpha(0.0);
    //       Fastor::Tensor<double, 3, 3> P_alpha(0.0);
    //       Fastor::Tensor<double, 3, 3> omega_alpha(0.0);

    //       F_n_alpha(0) = MatList[Phase2Material[h]].SP_norm[s].component(0); // n_alpha[Phase2Material[h]][s].component(0);
    //       F_n_alpha(1) = MatList[Phase2Material[h]].SP_norm[s].component(1); // n_alpha[Phase2Material[h]][s].component(1);
    //       F_n_alpha(2) = MatList[Phase2Material[h]].SP_norm[s].component(2); // n_alpha[Phase2Material[h]][s].component(2);

    //       F_m_alpha(0) = MatList[Phase2Material[h]].SP_direction[s].component(0); // m_alpha[Phase2Material[h]][s].component(0);
    //       F_m_alpha(1) = MatList[Phase2Material[h]].SP_direction[s].component(1); // m_alpha[Phase2Material[h]][s].component(1);
    //       F_m_alpha(2) = MatList[Phase2Material[h]].SP_direction[s].component(2); // m_alpha[Phase2Material[h]][s].component(2);


    //       Fastor::Tensor<double, 3> F_n_alpha_Rot(0.0);
    //       Fastor::Tensor<double, 3> F_m_alpha_Rot(0.0);

    //       for (int i = 0; i < 3; i++)
    //       {
    //         for (int j = 0; j < 3; j++)
    //         { // Rotated slip information
    //           F_n_alpha_Rot(i) += Rot_matrix(i, j) * F_n_alpha(j);
    //           F_m_alpha_Rot(i) += Rot_matrix(i, j) * F_m_alpha(j);
    //         }
    //       }


    //       P_alpha = 0.5 * (outer(F_m_alpha_Rot, F_n_alpha_Rot) + outer(F_n_alpha_Rot, F_m_alpha_Rot));
    //       omega_alpha = 0.5 * (outer(F_m_alpha_Rot, F_n_alpha_Rot) - outer(F_n_alpha_Rot, F_m_alpha_Rot));

    //       Fastor::Tensor<double, 3, 3> Schmid_T(0.0);

    //       Schmid_T = outer(F_m_alpha_Rot, F_n_alpha_Rot);
          
    //       double tau_alpha = 0.0;
    //       for (int i = 0; i < 3; i++)
    //       {
    //         for (int j = 0; j < 3; j++)
    //         {
    //           // tau_alpha += Fsigma(i, j)* Schmid_T(i, j);//Fsigmaav(i,j) 
    //           tau_alpha += Fsigma(i, j)* P_alpha(i, j);
    //         }
    //       }

    //       double k_exponent = MatList[Phase2Material[h]].k_exponent;         // 10.0;//12.0
    //       double ref_shear_rate = MatList[Phase2Material[h]].ref_shear_rate; // 1e-3;
          

    //       gamma_alpha_dots[h][s][cell] = ref_shear_rate * Foam::pow((mag(tau_alpha / g_alphas[h][s][cell])), k_exponent) * sign(tau_alpha); // source term for RHS of g_alpha equations


    //     } // end of loop over slip systems






    // }
}

        }

        sigtemp[h] = sigma_working;
        sigtemp[h].correctBoundaryConditions();

  // sigtemp[h] = sigma[h] +runTime.deltaT() * (dsigmadt[h]);//+fvc::div(phi,sigma)
  // sigtemp[h].correctBoundaryConditions();

  sigmaGlobaltemp+=PopBal[h]*sigtemp[h];
            }
            sigmaGlobaltemp.correctBoundaryConditions();



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

p=-(1.0/3.0)*tr(sigmaGlobal);




// Info << "Recomputing slip rates with updated stress" << endl;
//     for (label i = 0; i < N_Ori; i++) {
//         for (label countS = 0; countS < MatList[Phase2Material[i]].N_Slips; countS++) {
//             gamma_alpha_dots[i][countS] *= 0.0;
//         }
//     }





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