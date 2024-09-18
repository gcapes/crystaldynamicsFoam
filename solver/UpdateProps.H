{


  // gradU = fvc::grad(U);

  volTensorField TgradU(fvc::grad(U)().T());
  // TgradU.correctBoundaryConditions(); // already called by line above

  omega_dot = 0.5 * (TgradU - fvc::grad(U));   

  //NOTE!!! Openfoam defines grad in the opposite manner to solid mechanics textbooks, it's actually grad transpose that we want so this has to be transposed again or just ahev a - sign

  epsilon_dot = 0.5 * (TgradU + fvc::grad(U));//- (1.0/3.0)*(fvc::div(phi))*I; // this looks like dev(gradU)

  // dsigmadt = -(sigma & omega_dot) + (omega_dot & sigma) ;//-(tr(epsilon_dot)*sigma);;// 

dsigmadt*=0.0;



//
// gradU = TgradU;//fvc::grad(U);
// L = epsilon_dot + omega_dot;
//


  for (label h = 0; h < N_Ori; h++)
  {
    Omega_P[h] *= 0.0;
    volScalarField & ni = PopBal[h];

    PopBal_av[h] = fvc::average(ni); //takes the average from all cells sharing a face - this way the information about phases in adjacent cells is shared with the local cell
    PopBal_av[h].correctBoundaryConditions();
  }

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



//   // // //update hardening rates
// // // // Info<<"HERE1"<<endl;

//   forAll(mesh.C(), cell)
//   {
//     //   label Mat_name =-1;
//     //   double max_phase = -1.0;
//     for (label h = 0; h < N_Ori; h++)
//     {
//       if (PopBal[h][cell] > SMALL)
//       {

//         for (int S_alpha = 0; S_alpha < MatList[Phase2Material[h]].N_Slips; S_alpha++)
//         {


//           g_alpha_dots[h][S_alpha][cell] *= 0.0;

//           for (int S_beta = 0; S_beta < MatList[Phase2Material[h]].N_Slips; S_beta++)
//           {

//             //   // double h_of_gamma = ho*Foam::pow(1.0/(Foam::cosh((ho*Acc_Strain[cell])/(ta_s-tau_0))),2.0);
//             double h_alpha_beta = 0.0;
//             // Pout<<Pstream::myProcNo()<<"test_strain: "<<Acc_Strain[cell]<<endl;

//             double argument = (MatList[Phase2Material[h]].h_0 * Acc_Strain[cell]/*gamma_alphas[h][S_alpha][cell]*/) / (MatList[Phase2Material[h]].tau_S - MatList[Phase2Material[h]].initial_CRSS);
//             argument = max(min(argument,20.0),-20.0);

//             // double h_of_gamma = MatList[Phase2Material[h]].h_0 * Foam::pow(1.0 / (Foam::cosh((MatList[Phase2Material[h]].h_0 * gamma_alphas[h][S_alpha][cell]) / (MatList[Phase2Material[h]].tau_S - MatList[Phase2Material[h]].initial_CRSS))), 2.0);
//               double h_of_gamma = MatList[Phase2Material[h]].h_0 * Foam::pow(1.0 / (Foam::cosh(argument)), 2.0);            

//             if (S_beta == S_alpha)
//             {
//               h_alpha_beta = h_of_gamma * (MatList[Phase2Material[h]].r + (1.0 - MatList[Phase2Material[h]].r));
//             }
//             else
//             {
//               h_alpha_beta = h_of_gamma * (MatList[Phase2Material[h]].r);
//             }

//             g_alpha_dots[h][S_alpha][cell] += /*PopBal[h][cell]*/(h_alpha_beta)*gamma_alpha_dots[h][S_beta][cell];
//           }
//           double deltaT_ = (1.0 / n_sub_steps) * runTime.deltaT().value(); 
//           g_alpha_dots[h][S_alpha][cell]=max(min(g_alpha_dots[h][S_alpha][cell],0.99/deltaT_),-0.99/deltaT_);
//           //
//         }

//     }

//     }
//   }


// //Over Boundaries
// forAll(dsigmadt.boundaryField(), patchI)
// {
//     if
//     (
//         !dsigmadt.boundaryField()[patchI].coupled()
//      && dsigmadt.boundaryField()[patchI].size() > 0
//     )
//     {
//         forAll(dsigmadt.boundaryField()[patchI], faceI)
//         {


//    for (label h = 0; h < N_Ori; h++)
//     {
//             if (PopBal[h].boundaryField()[patchI][faceI] > SMALL)
//       {


//         for (int S_alpha = 0; S_alpha < MatList[Phase2Material[h]].N_Slips; S_alpha++)
//         {


//           g_alpha_dots[h][S_alpha].boundaryFieldRef()[patchI][faceI] *= 0.0;

//           for (int S_beta = 0; S_beta < MatList[Phase2Material[h]].N_Slips; S_beta++)
//           {


//             double h_alpha_beta = 0.0;

//             double argument = (MatList[Phase2Material[h]].h_0 * Acc_Strain.boundaryField()[patchI][faceI]/*gamma_alphas[h][S_alpha].boundaryField()[patchI][faceI]*/) / (MatList[Phase2Material[h]].tau_S - MatList[Phase2Material[h]].initial_CRSS);
//             argument = max(min(argument,20.0),-20.0);

//               double h_of_gamma = MatList[Phase2Material[h]].h_0 * Foam::pow(1.0 / (Foam::cosh(argument)), 2.0);            

//             if (S_beta == S_alpha)
//             {
//               h_alpha_beta = h_of_gamma * (MatList[Phase2Material[h]].r + (1.0 - MatList[Phase2Material[h]].r));
//             }
//             else
//             {
//               h_alpha_beta = h_of_gamma * (MatList[Phase2Material[h]].r);
//             }

//             g_alpha_dots[h][S_alpha].boundaryFieldRef()[patchI][faceI] += /*PopBal[h][cell]*/(h_alpha_beta)*gamma_alpha_dots[h][S_beta].boundaryField()[patchI][faceI];
//           }
//           double deltaT_ = (1.0 / n_sub_steps) * runTime.deltaT().value(); 
//           g_alpha_dots[h][S_alpha].boundaryFieldRef()[patchI][faceI]=max(min(g_alpha_dots[h][S_alpha].boundaryField()[patchI][faceI],0.99/deltaT_),-0.99/deltaT_);
//           //
//         }

//     }


//     }



//         }

//     }

// }


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
        Fsigma(i, j) = sigma[cell].component(count_index);
        // Fsigmaav(i, j) = sigmaav[cell].component(count_index);

        count_index++;
      }
    }


    Fastor::Tensor<double, 3, 3> sigAv(0.0);
    Fastor::Tensor<double, 3, 3> D_P_Average(0.0);

    label largest_field=-1;
    double current_largest_value = -1.0;

    sumgamma[cell]*=0.0;
    for (label h = 0; h < N_Ori; h++)
    {

      if(PopBal[h][cell]>current_largest_value){
        current_largest_value = PopBal[h][cell];
      
        largest_field = h;
      }

      Fastor::Tensor<double, 3, 3, 3, 3> Cprime(0.0); // global stiffness at this cell

      double phi_1_total = Rot1[h][cell]; // ROT_1[cell];// phi1[h] ;//+ omega_tot[cell].component(5);

      double phi_2_total = Rot2[h][cell]; // ROT_2[cell];// phi2[h] ;//+ omega_tot[cell].component(2);

      double phi_3_total = Rot3[h][cell]; // ROT_3[cell];// phi3[h] ;//+ omega_tot[cell].component(1);

    // if(PopBal[h][cell]>0.5){
      // Euler_Angles[cell].component(0) += PopBal[h][cell]*Rot1[h][cell];
      // Euler_Angles[cell].component(1) += PopBal[h][cell]*Rot2[h][cell];
      // Euler_Angles[cell].component(2) += PopBal[h][cell]*Rot3[h][cell];
    // }

      // vector Basis_1(1.0,0.0,0.0);
      // vector Basis_2(0.0,1.0,0.0);
      // vector Basis_3(0.0,0.0,1.0);


      Fastor::Tensor<double, 3, 3> Rot_matrix(0.0);

      // Rot_matrix(0, 0) = test_rot_diff[h][cell].component(0);// Foam::cos(phi_2_total) * Foam::cos(phi_3_total);
      // Rot_matrix(0, 1) = test_rot_diff[h][cell].component(1);//(Foam::sin(phi_1_total) * Foam::sin(phi_2_total) * Foam::cos(phi_3_total)) - (Foam::cos(phi_1_total) * Foam::sin(phi_3_total));
      // Rot_matrix(0, 2) = test_rot_diff[h][cell].component(2);//(Foam::cos(phi_1_total) * Foam::sin(phi_2_total) * Foam::cos(phi_3_total)) + (Foam::sin(phi_1_total) * Foam::sin(phi_3_total));

      // Rot_matrix(1, 0) = test_rot_diff[h][cell].component(3);//Foam::cos(phi_2_total) * Foam::sin(phi_3_total);
      // Rot_matrix(1, 1) = test_rot_diff[h][cell].component(4);//(Foam::sin(phi_1_total) * Foam::sin(phi_2_total) * Foam::sin(phi_3_total)) + (Foam::cos(phi_1_total) * Foam::cos(phi_3_total));
      // Rot_matrix(1, 2) = test_rot_diff[h][cell].component(5);//(Foam::cos(phi_1_total) * Foam::sin(phi_2_total) * Foam::sin(phi_3_total)) - (Foam::sin(phi_1_total) * Foam::cos(phi_3_total));

      // Rot_matrix(2, 0) = test_rot_diff[h][cell].component(6);//-Foam::sin(phi_2_total);
      // Rot_matrix(2, 1) =test_rot_diff[h][cell].component(7);// Foam::sin(phi_1_total) * Foam::cos(phi_2_total);
      // Rot_matrix(2, 2) = test_rot_diff[h][cell].component(8);//Foam::cos(phi_1_total) * Foam::cos(phi_2_total);

      Rot_matrix(0, 0) = Foam::cos(phi_2_total) * Foam::cos(phi_3_total);
      Rot_matrix(0, 1) = (Foam::sin(phi_1_total) * Foam::sin(phi_2_total) * Foam::cos(phi_3_total)) - (Foam::cos(phi_1_total) * Foam::sin(phi_3_total));
      Rot_matrix(0, 2) = (Foam::cos(phi_1_total) * Foam::sin(phi_2_total) * Foam::cos(phi_3_total)) + (Foam::sin(phi_1_total) * Foam::sin(phi_3_total));

      Rot_matrix(1, 0) = Foam::cos(phi_2_total) * Foam::sin(phi_3_total);
      Rot_matrix(1, 1) = (Foam::sin(phi_1_total) * Foam::sin(phi_2_total) * Foam::sin(phi_3_total)) + (Foam::cos(phi_1_total) * Foam::cos(phi_3_total));
      Rot_matrix(1, 2) = (Foam::cos(phi_1_total) * Foam::sin(phi_2_total) * Foam::sin(phi_3_total)) - (Foam::sin(phi_1_total) * Foam::cos(phi_3_total));

      Rot_matrix(2, 0) = -Foam::sin(phi_2_total);
      Rot_matrix(2, 1) = Foam::sin(phi_1_total) * Foam::cos(phi_2_total);
      Rot_matrix(2, 2) = Foam::cos(phi_1_total) * Foam::cos(phi_2_total);


    //   // Fastor::print("Rot Matrix: ",Rot_matrix);

    //   // https://en.wikipedia.org/wiki/Rotation_matrix -
    //   //"represents an extrinsic rotation whose (improper) Euler angles are α, β, γ, about axes x, y, z" - phi1 ==alpha, phi2==beta,phi3==gamma

    //     double deltaT_ = /*(1.0 / n_sub_steps) */ runTime.deltaT().value(); 
    //     // Info<<deltaT_<<endl;
    //     Fastor::Tensor<double, 3, 3> Elastic_Spin_h_deltaT(0.0);
    //     Fastor::Tensor<double, 3, 3> Elastic_Spin_h_deltaT_sq(0.0);
    //     int count_index_E = 0;
    //     for (int ii = 0; ii < 3; ii++)
    //     {
    //       for (int jj = 0; jj < 3; jj++)
    //       {
    //     // Info<<i<<","<<j<<endl;
    //         Elastic_Spin_h_deltaT(ii, jj) = (omega_dot[cell].component(count_index_E) - Omega_P[h][cell].component(count_index_E));

    //         count_index_E++;
    //       }
    //     }
    //     // Info<<omega_dot[cell].component(7)<<endl;;
    //     // Fastor::print("TEST",Elastic_Spin_h_deltaT);
    //     Elastic_Spin_h_deltaT_sq = Fastor::matmul(Elastic_Spin_h_deltaT,Elastic_Spin_h_deltaT);  



    //  //Update rotation Matrix Based on Elastic spin

    //  Fastor::Tensor<double, 3, 3> Rot_matrix_N_plus_one(0.0);
    //  Fastor::Tensor<double, 3, 3> I_d(0.0);
    //  for(int set=0;set<3;set++){
    //   I_d(set,set)=1.0;
    //  }
      
    // //  Fastor::print("I",I_d);

    // double X_comps = Foam::sqrt(
    //                 // Foam::pow(Elastic_Spin_h_deltaT(0,0),2.0)
    //                 // +Foam::pow(Elastic_Spin_h_deltaT(0, 1),2.0)
    //                 // +Foam::pow(Elastic_Spin_h_deltaT(0, 2),2.0)
    //                 // +Foam::pow(Elastic_Spin_h_deltaT(1, 0),2.0)
    //                 // +Foam::pow(Elastic_Spin_h_deltaT(1, 1),2.0)
    //                 // +Foam::pow(Elastic_Spin_h_deltaT(1, 2),2.0)
    //                 // +Foam::pow(Elastic_Spin_h_deltaT(2, 0),2.0)
    //                 // +Foam::pow(Elastic_Spin_h_deltaT(2, 1),2.0)
    //                 // +Foam::pow(Elastic_Spin_h_deltaT(2,2),2.0)

    //                 Foam::pow(Elastic_Spin_h_deltaT(2,1),2.0)
    //                 +Foam::pow(Elastic_Spin_h_deltaT(0, 2),2.0)
    //                 +Foam::pow(Elastic_Spin_h_deltaT(1,0),2.0)
                    
    //                 );

    // //                 Fastor::print("Elastic_Spin_h_deltaT",Elastic_Spin_h_deltaT);

    // // Info<<"X_comps: "<<X_comps<<endl;
    // Fastor::Tensor<double, 3, 3> e_to_the_C(0.0);
    // e_to_the_C = I_d;
    // if(X_comps>SMALL){
    //   // Info<<"x_comps"<<X_comps<<endl;
    

    // e_to_the_C  +=  ((Foam::sin(X_comps*deltaT_)/X_comps)*Elastic_Spin_h_deltaT) + ((1.0-Foam::cos(X_comps*deltaT_))/Foam::pow(X_comps,2.0))* Elastic_Spin_h_deltaT_sq;


    
    // }

    // Rot_matrix_N_plus_one = Fastor::matmul(Rot_matrix,Fastor::transpose(e_to_the_C));//If we define euler angles relative to +ve rotation around axes then this update has to be consistent with this
    // //e_to_the_C is the rotation update, now need to rotate the rotation matrix by this

    //     // Rot_matrix_N_plus_one = Fastor::matmul(Fastor::transpose(e_to_the_C),Fastor::matmul(Rot_matrix,e_to_the_C));//incremental update of regular tensor

      

    // // test_rot_diff[h][cell].component(0) = Rot_matrix_N_plus_one(0,0);
    // // test_rot_diff[h][cell].component(1) = Rot_matrix_N_plus_one(0,1);
    // // test_rot_diff[h][cell].component(2) = Rot_matrix_N_plus_one(0,2);
    // // test_rot_diff[h][cell].component(3) = Rot_matrix_N_plus_one(1,0);
    // // test_rot_diff[h][cell].component(4) = Rot_matrix_N_plus_one(1,1);
    // // test_rot_diff[h][cell].component(5) = Rot_matrix_N_plus_one(1,2);
    // // test_rot_diff[h][cell].component(6) = Rot_matrix_N_plus_one(2,0);
    // // test_rot_diff[h][cell].component(7) = Rot_matrix_N_plus_one(2,1);
    // // test_rot_diff[h][cell].component(8) = Rot_matrix_N_plus_one(2,2);
    
    // // if(cell==240&&h==1){
    // //   Info<<h<<endl;
    // // Fastor::print("Rot: ",Rot_matrix);
    // // Fastor::print("e_to_the_C: ",e_to_the_C);

    // // Fastor::print("Rot_n+1: ",Rot_matrix_N_plus_one);
    // // // exit(0);
    // // }

    // // Rot1[h][cell] = ?
    // // Rot2[h][cell] = ?
    // // Rot3[h][cell] = ?

    // double theta_angle = 0.0;
    // double psi_angle = 0.0;
    // double phi_angle = 0.0;

    // if(mag(Rot_matrix_N_plus_one(2,0))!=1.0){
    //   // Info<<"HERE"<<endl;
    //   double theta_1 = -std::asin(Rot_matrix_N_plus_one(2,0));
    //   double theta_2 = M_PI - theta_1;
    //   double psi_1 = std::atan2((Rot_matrix_N_plus_one(2,1)/std::cos(theta_1)),(Rot_matrix_N_plus_one(2,2)/std::cos(theta_1)));
    //   double psi_2 = std::atan2((Rot_matrix_N_plus_one(2,1)/std::cos(theta_2)),(Rot_matrix_N_plus_one(2,2)/std::cos(theta_2)));
    //   double phi_1 = std::atan2((Rot_matrix_N_plus_one(1,0)/std::cos(theta_1)),(Rot_matrix_N_plus_one(0,0)/std::cos(theta_1)));
    //   double phi_2 = std::atan2((Rot_matrix_N_plus_one(1,0)/std::cos(theta_2)),(Rot_matrix_N_plus_one(0,0)/std::cos(theta_2)));

    //   // Info<<"angles"<<endl;
    //   // Info<<"theta1: "<<theta_1<<endl;
    //   // Info<<"theta_2: "<<theta_2<<endl;
    //   // Info<<"psi_1: "<<psi_1<<endl;
    //   // Info<<"psi_2: "<<psi_2<<endl;
    //   // Info<<"phi_1: "<<phi_1<<endl;
    //   // Info<<"phi_2: "<<phi_2<<endl;
    //   // Info<<""<<endl;

    //   theta_angle = theta_1;//max(theta_1,theta_2);
    //   psi_angle = psi_1;//max(psi_1,psi_2);
    //   phi_angle = phi_1;//max(phi_1,phi_2);

    //   // Info<<"TEST: "<<psi_1<<"\t"<<psi_2<<endl;

    // }
    // else{
    //   phi_angle = 0.0;
    //     if(Rot_matrix_N_plus_one(2,0)==-1.0){
    //       theta_angle = M_PI/2.0;
    //       psi_angle = phi_angle + std::atan2(Rot_matrix_N_plus_one(0,1),(Rot_matrix_N_plus_one(0,2)));
    //     }
    //     else{
    //       theta_angle = -M_PI/2.0;
    //       psi_angle = -phi_angle + std::atan2(-Rot_matrix_N_plus_one(0,1),-Rot_matrix_N_plus_one(0,2));
    //     }

    // }

    // Rot1[h][cell] = psi_angle;//update euler angles for advection
    // Rot2[h][cell] = theta_angle;
    // Rot3[h][cell] = phi_angle;


      if (PopBal_av[h][cell] > SMALL)
      { // if the phase is there, work out the stiffness contribution

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
                          value += Rot_matrix(i, m) * Rot_matrix(j, n) * Rot_matrix(k, p) * Rot_matrix(l, q) * Cprime2[Phase2Material[h]](m, n, p, q); // C(m,n,p,q);
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

      }//

        // If grain is there
        // Have the stiffness tensor now
        // need Dp and omegaP

        Fastor::Tensor<double, 3, 3> F_D_P(0.0);
        Fastor::Tensor<double, 3, 3> F_Omega_P(0.0);

        
        for (int s = 0; s < MatList[Phase2Material[h]].N_Slips; s++)
        {

          Fastor::Tensor<double, 3> F_n_alpha(0.0);
          Fastor::Tensor<double, 3> F_m_alpha(0.0);
          Fastor::Tensor<double, 3, 3> P_alpha(0.0);
          Fastor::Tensor<double, 3, 3> omega_alpha(0.0);

          F_n_alpha(0) = MatList[Phase2Material[h]].SP_norm[s].component(0); // n_alpha[Phase2Material[h]][s].component(0);
          F_n_alpha(1) = MatList[Phase2Material[h]].SP_norm[s].component(1); // n_alpha[Phase2Material[h]][s].component(1);
          F_n_alpha(2) = MatList[Phase2Material[h]].SP_norm[s].component(2); // n_alpha[Phase2Material[h]][s].component(2);

          F_m_alpha(0) = MatList[Phase2Material[h]].SP_direction[s].component(0); // m_alpha[Phase2Material[h]][s].component(0);
          F_m_alpha(1) = MatList[Phase2Material[h]].SP_direction[s].component(1); // m_alpha[Phase2Material[h]][s].component(1);
          F_m_alpha(2) = MatList[Phase2Material[h]].SP_direction[s].component(2); // m_alpha[Phase2Material[h]][s].component(2);


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

          double k_exponent = MatList[Phase2Material[h]].k_exponent;         // 10.0;//12.0
          double ref_shear_rate = MatList[Phase2Material[h]].ref_shear_rate; // 1e-3;

          gamma_alpha_dots[h][s][cell] = ref_shear_rate * Foam::pow((mag(tau_alpha / g_alphas[h][s][cell])), k_exponent) * sign(tau_alpha); // source term for RHS of g_alpha equations


          sumgamma[cell]+= mag(gamma_alpha_dots[h][s][cell]);//for next time-step

          F_D_P += (P_alpha * gamma_alpha_dots[h][s][cell]);
          F_Omega_P += (omega_alpha * gamma_alpha_dots[h][s][cell]);
          // }

        } // end of loop over slip systems

        int OP_index = 0;
        for (int i = 0; i < 3; i++)
        {
          for (int j = 0; j < 3; j++)
          {
            Omega_P[h][cell].component(OP_index) = F_Omega_P(i, j);
            OP_index++;
          }
        }

        enum
        {
          I,J,K,L
        };
        Fastor::Tensor<double, 3, 3> sig = Fastor::einsum<Fastor::Index<I, J, K, L>, Fastor::Index<K, L> /*,Fastor::OIndex<I,J>*/>(Cprime, Fepsilon_dot - F_D_P); ////Double contraction

        sigAv += (PopBal[h][cell]/Sum_Phases[cell]) * sig;///Sum_Phases[cell];
      // }
        // Omega_P_Average+=PopBal[h][cell]*F_Omega_P;
        D_P_Average += /*PopBal[h][cell] */ F_D_P;

        Fastor::Tensor<double, 3, 3> F_D_E(0.0);

        F_D_E = Fepsilon_dot-F_D_P;

          dsigmadt[cell] -=(PopBal[h][cell]/Sum_Phases[cell]) *((sigma[cell] & omega_dot[cell]) + (omega_dot[cell] & sigma[cell])) ;//-(tr(epsilon_dot)*sigma);;// 
          dsigmadt[cell] += (PopBal[h][cell]/Sum_Phases[cell]) *((sigma[cell] & Omega_P[h][cell]) - (Omega_P[h][cell] & sigma[cell]));
          dsigmadt[cell] -= (PopBal[h][cell]/Sum_Phases[cell]) *(sigma[cell]*(F_D_E(0,0)+F_D_E(1,1)+F_D_E(2,2)));




    } // end loop over phase fields


      // Euler_Angles[cell].component(0) = Rot1[largest_field][cell];
      // Euler_Angles[cell].component(1) = Rot2[largest_field][cell];
      // Euler_Angles[cell].component(2) = Rot3[largest_field][cell];

    count_index = 0;
    for (int i = 0; i < 3; i++)
    {
      for (int j = 0; j < 3; j++)
      {

        dsigmadt[cell].component(count_index) += sigAv(i, j); // sig(i,j);
        D_P[cell].component(count_index) = D_P_Average(i, j); // F_D_P(i,j);
        count_index++;
      }
    }
  }











//Over Boundaries
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
        Fsigma(i, j) = sigma.boundaryField()[patchI][faceI].component(count_index);
        // Fsigmaav(i, j) = sigmaav[cell].component(count_index);

        count_index++;
      }
    }


    Fastor::Tensor<double, 3, 3> sigAv(0.0);
    Fastor::Tensor<double, 3, 3> D_P_Average(0.0);


        label largest_field=-1;
    double current_largest_value = -1.0;

    sumgamma.boundaryFieldRef()[patchI][faceI]*=0.0;
    for (label h = 0; h < N_Ori; h++)
    {


        if(PopBal[h].boundaryField()[patchI][faceI]>current_largest_value){
        current_largest_value = PopBal[h].boundaryField()[patchI][faceI];
      
        largest_field = h;
      }



       Fastor::Tensor<double, 3, 3, 3, 3> Cprime(0.0); // global stiffness at this cell

      double phi_1_total = Rot1[h].boundaryField()[patchI][faceI]; // ROT_1[cell];// phi1[h] ;//+ omega_tot[cell].component(5);

      double phi_2_total = Rot2[h].boundaryField()[patchI][faceI]; // ROT_2[cell];// phi2[h] ;//+ omega_tot[cell].component(2);

      double phi_3_total = Rot3[h].boundaryField()[patchI][faceI]; // ROT_3[cell];// phi3[h] ;//+ omega_tot[cell].component(1);

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


    //     double deltaT_ = /*(1.0 / n_sub_steps) */ runTime.deltaT().value(); 

    //     Fastor::Tensor<double, 3, 3> Elastic_Spin_h_deltaT(0.0);
    //     Fastor::Tensor<double, 3, 3> Elastic_Spin_h_deltaT_sq(0.0);
    //     int count_index_E = 0;
    //     for (int ii = 0; ii < 3; ii++)
    //     {
    //       for (int jj = 0; jj < 3; jj++)
    //       {
    //     // Info<<i<<","<<j<<endl;
    //         Elastic_Spin_h_deltaT(ii, jj) = deltaT_*(omega_dot.boundaryField()[patchI][faceI].component(count_index_E) - Omega_P[h].boundaryField()[patchI][faceI].component(count_index_E));
    //         // Elastic_Spin_h_sq(ii, jj) = Foam::pow(omega_dot[cell].component(count_index_E) - Omega_P[h][cell].component(count_index_E),2.0);
    //         count_index_E++;
    //       }
    //     }

    //     Elastic_Spin_h_deltaT_sq = Fastor::matmul(Elastic_Spin_h_deltaT,Elastic_Spin_h_deltaT);  



    //  //Update rotation Matrix Based on Elastic spin

    //  Fastor::Tensor<double, 3, 3> Rot_matrix_N_plus_one(0.0);
    //  Fastor::Tensor<double, 3, 3> I_d(0.0);
    //  for(int set=0;set<3;set++){
    //   I_d(set,set)=1.0;
    //  }
      
    // //  Fastor::print("I",I_d);

    // double X_comps = Foam::sqrt(Foam::pow(Elastic_Spin_h_deltaT(2,1),2.0)
    //                 +Foam::pow(Elastic_Spin_h_deltaT(0, 2),2.0)
    //                 +Foam::pow(Elastic_Spin_h_deltaT(1,0),2.0));

    //                 // Fastor::print("Elastic_Spin_h_deltaT",Elastic_Spin_h_deltaT);

    // // Info<<"X_comps: "<<X_comps<<endl;
    // Fastor::Tensor<double, 3, 3> e_to_the_C(0.0);
    // e_to_the_C = I_d;
    // if(X_comps>SMALL){
    //   // Info<<"x_comps"<<X_comps<<endl;
    

    // e_to_the_C  += ((Foam::sin(X_comps)/X_comps)*Elastic_Spin_h_deltaT) + ((1.0-Foam::cos(X_comps))/Foam::pow(X_comps,2.0))* Elastic_Spin_h_deltaT_sq;


    
    // }
    // Rot_matrix_N_plus_one = Fastor::matmul(Rot_matrix,Fastor::transpose(e_to_the_C));//If we define euler angles relative to +ve rotation around axes then this update has to be consistent with this
    


    // // Rot1[h][cell] = ?
    // // Rot2[h][cell] = ?
    // // Rot3[h][cell] = ?

    // double theta_angle = 0.0;
    // double psi_angle = 0.0;
    // double phi_angle = 0.0;

    // if(mag(Rot_matrix_N_plus_one(2,0))!=1.0){
    //   double theta_1 = -std::asin(Rot_matrix_N_plus_one(2,0));
    //   double theta_2 = M_PI - theta_1;
    //   double psi_1 = std::atan2((Rot_matrix_N_plus_one(2,1)/std::cos(theta_1)),(Rot_matrix_N_plus_one(2,2)/std::cos(theta_1)));
    //   double psi_2 = std::atan2((Rot_matrix_N_plus_one(2,1)/std::cos(theta_2)),(Rot_matrix_N_plus_one(2,2)/std::cos(theta_2)));
    //   double phi_1 = std::atan2((Rot_matrix_N_plus_one(1,0)/std::cos(theta_1)),(Rot_matrix_N_plus_one(0,0)/std::cos(theta_1)));
    //   double phi_2 = std::atan2((Rot_matrix_N_plus_one(1,0)/std::cos(theta_2)),(Rot_matrix_N_plus_one(0,0)/std::cos(theta_2)));

    //   theta_angle = theta_1;//max(theta_1,theta_2);
    //   psi_angle = psi_1;//max(psi_1,psi_2);
    //   phi_angle = phi_1;//max(phi_1,phi_2);

    //   // Info<<"TEST: "<<psi_1<<"\t"<<psi_2<<endl;

    // }
    // else{
    //   phi_angle = 0.0;
    //     if(Rot_matrix_N_plus_one(2,0)==-1.0){
    //       theta_angle = M_PI/2.0;
    //       psi_angle = phi_angle + std::atan2(Rot_matrix_N_plus_one(0,1),(Rot_matrix_N_plus_one(0,2)));
    //     }
    //     else{
    //       theta_angle = -M_PI/2.0;
    //       psi_angle = -phi_angle + std::atan2(-Rot_matrix_N_plus_one(0,1),-Rot_matrix_N_plus_one(0,2));
    //     }

    // }

    // Rot1[h].boundaryFieldRef()[patchI][faceI] = psi_angle;//update euler angles for advection
    // Rot2[h].boundaryFieldRef()[patchI][faceI] = theta_angle;
    // Rot3[h].boundaryFieldRef()[patchI][faceI] = phi_angle;



      if (PopBal_av[h].boundaryField()[patchI][faceI] > SMALL)
      { // if the phase is there, work out the stiffness contribution

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
                          value += Rot_matrix(i, m) * Rot_matrix(j, n) * Rot_matrix(k, p) * Rot_matrix(l, q) * Cprime2[Phase2Material[h]](m, n, p, q); // C(m,n,p,q);
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

      }//


// If grain is there
        // Have the stiffness tensor now
        // need Dp and omegaP

        Fastor::Tensor<double, 3, 3> F_D_P(0.0);
        Fastor::Tensor<double, 3, 3> F_Omega_P(0.0);


      for (int s = 0; s < MatList[Phase2Material[h]].N_Slips; s++)
        {

          Fastor::Tensor<double, 3> F_n_alpha(0.0);
          Fastor::Tensor<double, 3> F_m_alpha(0.0);
          Fastor::Tensor<double, 3, 3> P_alpha(0.0);
          Fastor::Tensor<double, 3, 3> omega_alpha(0.0);

          F_n_alpha(0) = MatList[Phase2Material[h]].SP_norm[s].component(0); // n_alpha[Phase2Material[h]][s].component(0);
          F_n_alpha(1) = MatList[Phase2Material[h]].SP_norm[s].component(1); // n_alpha[Phase2Material[h]][s].component(1);
          F_n_alpha(2) = MatList[Phase2Material[h]].SP_norm[s].component(2); // n_alpha[Phase2Material[h]][s].component(2);

          F_m_alpha(0) = MatList[Phase2Material[h]].SP_direction[s].component(0); // m_alpha[Phase2Material[h]][s].component(0);
          F_m_alpha(1) = MatList[Phase2Material[h]].SP_direction[s].component(1); // m_alpha[Phase2Material[h]][s].component(1);
          F_m_alpha(2) = MatList[Phase2Material[h]].SP_direction[s].component(2); // m_alpha[Phase2Material[h]][s].component(2);


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

          double k_exponent = MatList[Phase2Material[h]].k_exponent;         // 10.0;//12.0
          double ref_shear_rate = MatList[Phase2Material[h]].ref_shear_rate; // 1e-3;

          gamma_alpha_dots[h][s].boundaryFieldRef()[patchI][faceI] = ref_shear_rate * Foam::pow((mag(tau_alpha / g_alphas[h][s].boundaryField()[patchI][faceI])), k_exponent) * sign(tau_alpha); // source term for RHS of g_alpha equations


          sumgamma.boundaryFieldRef()[patchI][faceI]+= mag(gamma_alpha_dots[h][s].boundaryField()[patchI][faceI]);//for next time-step

          F_D_P += (P_alpha * gamma_alpha_dots[h][s].boundaryField()[patchI][faceI]);
          F_Omega_P += (omega_alpha * gamma_alpha_dots[h][s].boundaryField()[patchI][faceI]);
          // }

        } // end of loop over slip systems

        int OP_index = 0;
        for (int i = 0; i < 3; i++)
        {
          for (int j = 0; j < 3; j++)
          {
            Omega_P[h].boundaryFieldRef()[patchI][faceI].component(OP_index) = F_Omega_P(i, j); //Philip, how to assign the boundary face values
            OP_index++;
          }
        }

        enum
        {
          I,J,K,L
        };
        Fastor::Tensor<double, 3, 3> sig = Fastor::einsum<Fastor::Index<I, J, K, L>, Fastor::Index<K, L> /*,Fastor::OIndex<I,J>*/>(Cprime, Fepsilon_dot - F_D_P); ////Double contraction


        sigAv += (PopBal[h].boundaryField()[patchI][faceI]/Sum_Phases.boundaryField()[patchI][faceI]) * sig;///Sum_Phases[cell];
        // Omega_P_Average+=PopBal[h][cell]*F_Omega_P;
      // }
        D_P_Average += /*PopBal[h][cell] */ F_D_P;

        Fastor::Tensor<double, 3, 3> F_D_E(0.0);

        F_D_E = Fepsilon_dot-F_D_P;



          dsigmadt.boundaryFieldRef()[patchI][faceI] -=(PopBal[h].boundaryField()[patchI][faceI]/Sum_Phases.boundaryField()[patchI][faceI]) *((sigma.boundaryField()[patchI][faceI] & omega_dot.boundaryField()[patchI][faceI]) + (omega_dot.boundaryField()[patchI][faceI] & sigma.boundaryField()[patchI][faceI]));
          dsigmadt.boundaryFieldRef()[patchI][faceI] += (PopBal[h].boundaryField()[patchI][faceI]/Sum_Phases.boundaryField()[patchI][faceI]) *((sigma.boundaryField()[patchI][faceI] & Omega_P[h].boundaryField()[patchI][faceI]) - (Omega_P[h].boundaryField()[patchI][faceI] & sigma.boundaryField()[patchI][faceI]));//Philip, can we assign boundary fields with dot products in?
          dsigmadt.boundaryFieldRef()[patchI][faceI] -= (PopBal[h].boundaryField()[patchI][faceI]/Sum_Phases.boundaryField()[patchI][faceI]) *(sigma.boundaryField()[patchI][faceI]*(F_D_E(0,0)+F_D_E(1,1)+F_D_E(2,2)));



    }// end loop over phase fields

      // Euler_Angles.boundaryFieldRef()[patchI][faceI].component(0) = Rot1[largest_field].boundaryField()[patchI][faceI];
      // Euler_Angles.boundaryFieldRef()[patchI][faceI].component(1) = Rot2[largest_field].boundaryField()[patchI][faceI];
      // Euler_Angles.boundaryFieldRef()[patchI][faceI].component(2) = Rot3[largest_field].boundaryField()[patchI][faceI];


    count_index = 0;
    for (int i = 0; i < 3; i++)
    {
      for (int j = 0; j < 3; j++)
      {

        dsigmadt.boundaryFieldRef()[patchI][faceI].component(count_index) += sigAv(i, j); // sig(i,j);
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



dsigmadt.correctBoundaryConditions();




// //RK4

// //RK4


//   sigma = sigma + (1.0 / n_sub_steps)*  runTime.deltaT() * (dsigmadt);
//   sigma.correctBoundaryConditions();

// D = D + (1.0 / n_sub_steps) * runTime.deltaT() * (U);
// D.correctBoundaryConditions();

// p=-(1.0/3.0)*tr(sigma);

// testtraceeps = tr(epsilon_dot-D_P)*sigma;


}
