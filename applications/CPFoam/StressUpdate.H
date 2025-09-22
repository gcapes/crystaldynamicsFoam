{
    sigmaGlobaltemp *= 0.0;

    for (label h = 0; h < N_Ori; h++)
    {
        // if (pimpleIter == 1) {
        //     Info << "Initialising stress for sub-stepping" << endl;
        //     sigtemp[h] = sigma[h];  // Start from committed stress
        // }

        Info << "Number of Stress Integration Sub-Steps: " << n_substeps << endl;
        dimensionedScalar sub_dt = runTime.deltaT() / n_substeps;
        
        volTensorField sigma_working = sigtemp[h];
        
        // Initialize plastic deformation for this phase
        volTensorField F_D_P_local(
            IOobject(
                "F_D_P_local",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedTensor("zero", dimensionSet(0, 0, -1, 0, 0), tensor::zero)
        );
        
        volTensorField F_Omega_P_local = F_D_P_local;  // Same dimensions



        volTensorField dsigmadt_local(
            IOobject(
                "dsigmadt_local",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedTensor("zero", dimensionSet(1, -1, -3, 0, 0), tensor::zero)
        );



        // Sub-stepping loop
        for (int step = 0; step < n_substeps; step++) 
        {
            // Step 1: Update CP constitutive response based on current stress
            forAll(mesh.C(), cell)
            {

                    // Extract current stress into Fastor tensor
                    Fastor::Tensor<double, 3, 3> Fsigma(0.0);
                    Fastor::Tensor<double, 3, 3> Fepsilon_dot(0.0);
                    int count_index = 0;
                    for (int i = 0; i < 3; i++)
                    {
                        for (int j = 0; j < 3; j++)
                        {
                            Fsigma(i, j) = sigma_working[cell].component(count_index);
                            Fepsilon_dot(i, j) = epsilon_dot[cell].component(count_index);
                            count_index++;
                        }
                    }
                    

                    Fastor::Tensor<double, 3, 3, 3, 3> Cprime(0.0); // global stiffness at this cell


                    // Get rotation matrix for this phase
                    double phi_1 = Rot1[h][cell];
                    double phi_2 = Rot2[h][cell];
                    double phi_3 = Rot3[h][cell];
                    
                    Fastor::Tensor<double, 3, 3> Rot_matrix(0.0);
                    Rot_matrix(0, 0) = Foam::cos(phi_2) * Foam::cos(phi_3);
                    Rot_matrix(0, 1) = (Foam::sin(phi_1) * Foam::sin(phi_2) * Foam::cos(phi_3)) - 
                                       (Foam::cos(phi_1) * Foam::sin(phi_3));
                    Rot_matrix(0, 2) = (Foam::cos(phi_1) * Foam::sin(phi_2) * Foam::cos(phi_3)) + 
                                       (Foam::sin(phi_1) * Foam::sin(phi_3));
                    Rot_matrix(1, 0) = Foam::cos(phi_2) * Foam::sin(phi_3);
                    Rot_matrix(1, 1) = (Foam::sin(phi_1) * Foam::sin(phi_2) * Foam::sin(phi_3)) + 
                                       (Foam::cos(phi_1) * Foam::cos(phi_3));
                    Rot_matrix(1, 2) = (Foam::cos(phi_1) * Foam::sin(phi_2) * Foam::sin(phi_3)) - 
                                       (Foam::sin(phi_1) * Foam::cos(phi_3));
                    Rot_matrix(2, 0) = -Foam::sin(phi_2);
                    Rot_matrix(2, 1) = Foam::sin(phi_1) * Foam::cos(phi_2);
                    Rot_matrix(2, 2) = Foam::cos(phi_1) * Foam::cos(phi_2);

                if (PopBal_av[h][cell] > 0.0)
                {

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
                          value += Rot_matrix(i, m) * Rot_matrix(j, n) * Rot_matrix(k, p) * Rot_matrix(l, q) * Cprime2[Phase2Material[h]](m, n, p, q); 
                        // }
                      }
                    }
                  }
                }

                // value*=PopBal[h][cell];
// Info<<value<<endl;
                Cprime(i, j, k, l) = value;
              }
            }
          }
        }




                    
                    // Calculate plastic deformation rate for this sub-step
                    Fastor::Tensor<double, 3, 3> F_D_P(0.0);
                    Fastor::Tensor<double, 3, 3> F_Omega_P(0.0);
                    
                    for (int s = 0; s < MatList[Phase2Material[h]].N_Slips; s++)
                    {
                        // Get slip system vectors
                        Fastor::Tensor<double, 3> F_n_alpha(0.0);
                        Fastor::Tensor<double, 3> F_m_alpha(0.0);
                        
                        F_n_alpha(0) = MatList[Phase2Material[h]].SP_norm[s].component(0);
                        F_n_alpha(1) = MatList[Phase2Material[h]].SP_norm[s].component(1);
                        F_n_alpha(2) = MatList[Phase2Material[h]].SP_norm[s].component(2);
                        
                        F_m_alpha(0) = MatList[Phase2Material[h]].SP_direction[s].component(0);
                        F_m_alpha(1) = MatList[Phase2Material[h]].SP_direction[s].component(1);
                        F_m_alpha(2) = MatList[Phase2Material[h]].SP_direction[s].component(2);
                        
                        // Rotate slip systems
                        Fastor::Tensor<double, 3> F_n_alpha_Rot(0.0);
                        Fastor::Tensor<double, 3> F_m_alpha_Rot(0.0);
                        
                        for (int i = 0; i < 3; i++)
                        {
                            for (int j = 0; j < 3; j++)
                            {
                                F_n_alpha_Rot(i) += Rot_matrix(i, j) * F_n_alpha(j);
                                F_m_alpha_Rot(i) += Rot_matrix(i, j) * F_m_alpha(j);
                            }
                        }
                        
                        // Calculate Schmid tensor and resolved shear stress
                        Fastor::Tensor<double, 3, 3> P_alpha = 0.5 * (outer(F_m_alpha_Rot, F_n_alpha_Rot) + 
                                                                      outer(F_n_alpha_Rot, F_m_alpha_Rot));
                        Fastor::Tensor<double, 3, 3> omega_alpha = 0.5 * (outer(F_m_alpha_Rot, F_n_alpha_Rot) - 
                                                                          outer(F_n_alpha_Rot, F_m_alpha_Rot));
                        
                        double tau_alpha = 0.0;
                        for (int i = 0; i < 3; i++)
                        {
                            for (int j = 0; j < 3; j++)
                            {
                                tau_alpha += Fsigma(i, j) * P_alpha(i, j);
                            }
                        }
                        
                        // Update slip rate with current stress and hardening
                        double k_exponent = MatList[Phase2Material[h]].k_exponent;
                        double ref_shear_rate = MatList[Phase2Material[h]].ref_shear_rate;

                                  double tau_ratio = mag(tau_alpha)/g_alphas[h][s][cell];



                        gamma_alpha_dots[h][s][cell] = ref_shear_rate * 
                            Foam::pow(mag(tau_alpha / g_alphas[h][s][cell]), k_exponent) * 
                            sign(tau_alpha);
                        
                        // // Update local plastic deformation rate
                        // gamma_alpha_dots[h][s][cell] = gamma_dot_alpha;
                        F_D_P += P_alpha * gamma_alpha_dots[h][s][cell];//gamma_dot_alpha;
                        F_Omega_P += omega_alpha * gamma_alpha_dots[h][s][cell];//gamma_dot_alpha;
                        
                        // Update hardening (simplified - you may want full hardening law)
                        // g_alphas[h][s][cell] += sub_dt.value() * g_alpha_dots[h][s][cell];
                    }







                    
                    // Store plastic deformation tensors
        


                    enum
        {
          I,J,K,L
        };
        Fastor::Tensor<double, 3, 3> dsigma_elastic = Fastor::einsum<Fastor::Index<I, J, K, L>, Fastor::Index<K, L> /*,Fastor::OIndex<I,J>*/>(Cprime, Fepsilon_dot - F_D_P); ////Double contraction
        
            count_index = 0;
            for (int i = 0; i < 3; i++)
            {
                for (int j = 0; j < 3; j++)
                {
                            dsigmadt_local[cell].component(count_index) = dsigma_elastic(i, j);
                            F_D_P_local[cell].component(count_index) = F_D_P(i, j);
                            F_Omega_P_local[cell].component(count_index) = F_Omega_P(i, j);
                            Plastic_Strain[h][cell].component(count_index) = F_D_P(i, j);
                            Omega_P[h][cell].component(count_index) = F_Omega_P(i, j);
                            count_index++;
                }
            }



                }





//////////////////////////////////////////////////////////////////////////////////////////////////////////////ROtation_update
        Fastor::Tensor<double, 3, 3> Elastic_Spin_h_deltaT(0.0);
        Fastor::Tensor<double, 3, 3> Elastic_Spin_h_deltaT_sq(0.0);
        int count_index_E = 0;
        for (int ii = 0; ii < 3; ii++)
        {
          for (int jj = 0; jj < 3; jj++)
          {

            Elastic_Spin_h_deltaT(ii, jj) = (omega_dot[cell].component(count_index_E) - F_Omega_P_local[cell].component(count_index_E));//shouldnt be the elastic spin, should be the total spin elastic + plastic

            count_index_E++;
          }
        }

        Elastic_Spin_h_deltaT_sq = Fastor::matmul(Elastic_Spin_h_deltaT,Elastic_Spin_h_deltaT);  



     //Update rotation Matrix Based on Elastic spin

     Fastor::Tensor<double, 3, 3> Rot_matrix_N_plus_one(0.0);
     Fastor::Tensor<double, 3, 3> I_d(0.0);
     for(int set=0;set<3;set++){
      I_d(set,set)=1.0;
     }
      
    //  Fastor::print("I",I_d);

    double X_comps = Foam::sqrt(

                    Foam::pow(Elastic_Spin_h_deltaT(2,1),2.0)
                    +Foam::pow(Elastic_Spin_h_deltaT(0, 2),2.0)
                    +Foam::pow(Elastic_Spin_h_deltaT(1,0),2.0)
                    
                    );

    //                 Fastor::print("Elastic_Spin_h_deltaT",Elastic_Spin_h_deltaT);

    // Info<<"X_comps: "<<X_comps<<endl;
    Fastor::Tensor<double, 3, 3> e_to_the_C(0.0);
    e_to_the_C = I_d;
    if(X_comps>SMALL){
      // Info<<"x_comps"<<X_comps<<endl;
    

    e_to_the_C  +=  ((Foam::sin(X_comps*sub_dt.value())/X_comps)*Elastic_Spin_h_deltaT) + ((1.0-Foam::cos(X_comps*sub_dt.value()))/Foam::pow(X_comps,2.0))* Elastic_Spin_h_deltaT_sq;

//sub_dt
    
    }

     Rot_matrix_N_plus_one = Fastor::matmul(e_to_the_C,Rot_matrix);


     double theta_angle = 0.0;
    double psi_angle = 0.0;
    double phi_angle = 0.0;

    if(mag(Rot_matrix_N_plus_one(2,0))!=1.0){
      // Info<<"HERE"<<endl;
      double theta_1 = -std::asin(Rot_matrix_N_plus_one(2,0));
      double theta_2 = M_PI - theta_1;
      double psi_1 = std::atan2((Rot_matrix_N_plus_one(2,1)/std::cos(theta_1)),(Rot_matrix_N_plus_one(2,2)/std::cos(theta_1)));
      double psi_2 = std::atan2((Rot_matrix_N_plus_one(2,1)/std::cos(theta_2)),(Rot_matrix_N_plus_one(2,2)/std::cos(theta_2)));
      double phi_1 = std::atan2((Rot_matrix_N_plus_one(1,0)/std::cos(theta_1)),(Rot_matrix_N_plus_one(0,0)/std::cos(theta_1)));
      double phi_2 = std::atan2((Rot_matrix_N_plus_one(1,0)/std::cos(theta_2)),(Rot_matrix_N_plus_one(0,0)/std::cos(theta_2)));

      theta_angle = theta_1;//max(theta_1,theta_2);
      psi_angle = psi_1;//max(psi_1,psi_2);
      phi_angle = phi_1;//max(phi_1,phi_2);

      // Info<<"TEST: "<<psi_1<<"\t"<<psi_2<<endl;

    }
    else{
      phi_angle = 0.0;
        if(Rot_matrix_N_plus_one(2,0)==-1.0){
          theta_angle = M_PI/2.0;
          psi_angle = phi_angle + std::atan2(Rot_matrix_N_plus_one(0,1),(Rot_matrix_N_plus_one(0,2)));
        }
        else{
          theta_angle = -M_PI/2.0;
          psi_angle = -phi_angle + std::atan2(-Rot_matrix_N_plus_one(0,1),-Rot_matrix_N_plus_one(0,2));
        }

    }

    Rot1[h][cell] = psi_angle;//update euler angles for advection
    Rot2[h][cell] = theta_angle;
    Rot3[h][cell] = phi_angle;

//////////////////////////////////////////////////////////////////////////////////////////////////////////////ROtation_update













            }//end over cell loop











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





                    // Extract current stress into Fastor tensor
                    Fastor::Tensor<double, 3, 3> Fsigma(0.0);
                    Fastor::Tensor<double, 3, 3> Fepsilon_dot(0.0);
                    int count_index = 0;
                    for (int i = 0; i < 3; i++)
                    {
                        for (int j = 0; j < 3; j++)
                        {
                            Fsigma(i, j) = sigma_working.boundaryField()[patchI][faceI].component(count_index);
                            Fepsilon_dot(i, j) = epsilon_dot.boundaryField()[patchI][faceI].component(count_index);
                            count_index++;
                        }
                    }
                    

                    Fastor::Tensor<double, 3, 3, 3, 3> Cprime(0.0); // global stiffness at this cell


                    // Get rotation matrix for this phase
                    double phi_1 = Rot1[h].boundaryField()[patchI][faceI];
                    double phi_2 = Rot2[h].boundaryField()[patchI][faceI];
                    double phi_3 = Rot3[h].boundaryField()[patchI][faceI];
                    
                    Fastor::Tensor<double, 3, 3> Rot_matrix(0.0);
                    Rot_matrix(0, 0) = Foam::cos(phi_2) * Foam::cos(phi_3);
                    Rot_matrix(0, 1) = (Foam::sin(phi_1) * Foam::sin(phi_2) * Foam::cos(phi_3)) - 
                                       (Foam::cos(phi_1) * Foam::sin(phi_3));
                    Rot_matrix(0, 2) = (Foam::cos(phi_1) * Foam::sin(phi_2) * Foam::cos(phi_3)) + 
                                       (Foam::sin(phi_1) * Foam::sin(phi_3));
                    Rot_matrix(1, 0) = Foam::cos(phi_2) * Foam::sin(phi_3);
                    Rot_matrix(1, 1) = (Foam::sin(phi_1) * Foam::sin(phi_2) * Foam::sin(phi_3)) + 
                                       (Foam::cos(phi_1) * Foam::cos(phi_3));
                    Rot_matrix(1, 2) = (Foam::cos(phi_1) * Foam::sin(phi_2) * Foam::sin(phi_3)) - 
                                       (Foam::sin(phi_1) * Foam::cos(phi_3));
                    Rot_matrix(2, 0) = -Foam::sin(phi_2);
                    Rot_matrix(2, 1) = Foam::sin(phi_1) * Foam::cos(phi_2);
                    Rot_matrix(2, 2) = Foam::cos(phi_1) * Foam::cos(phi_2);


                if (PopBal_av[h].boundaryField()[patchI][faceI] > 0.0)
                {




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



                Cprime(i, j, k, l) = value;
              }
            }
          }
        }




                    
                    // Calculate plastic deformation rate for this sub-step
                    Fastor::Tensor<double, 3, 3> F_D_P(0.0);
                    Fastor::Tensor<double, 3, 3> F_Omega_P(0.0);
                    
                    for (int s = 0; s < MatList[Phase2Material[h]].N_Slips; s++)
                    {
                        // Get slip system vectors
                        Fastor::Tensor<double, 3> F_n_alpha(0.0);
                        Fastor::Tensor<double, 3> F_m_alpha(0.0);
                        
                        F_n_alpha(0) = MatList[Phase2Material[h]].SP_norm[s].component(0);
                        F_n_alpha(1) = MatList[Phase2Material[h]].SP_norm[s].component(1);
                        F_n_alpha(2) = MatList[Phase2Material[h]].SP_norm[s].component(2);
                        
                        F_m_alpha(0) = MatList[Phase2Material[h]].SP_direction[s].component(0);
                        F_m_alpha(1) = MatList[Phase2Material[h]].SP_direction[s].component(1);
                        F_m_alpha(2) = MatList[Phase2Material[h]].SP_direction[s].component(2);
                        
                        // Rotate slip systems
                        Fastor::Tensor<double, 3> F_n_alpha_Rot(0.0);
                        Fastor::Tensor<double, 3> F_m_alpha_Rot(0.0);
                        
                        for (int i = 0; i < 3; i++)
                        {
                            for (int j = 0; j < 3; j++)
                            {
                                F_n_alpha_Rot(i) += Rot_matrix(i, j) * F_n_alpha(j);
                                F_m_alpha_Rot(i) += Rot_matrix(i, j) * F_m_alpha(j);
                            }
                        }
                        
                        // Calculate Schmid tensor and resolved shear stress
                        Fastor::Tensor<double, 3, 3> P_alpha = 0.5 * (outer(F_m_alpha_Rot, F_n_alpha_Rot) + 
                                                                      outer(F_n_alpha_Rot, F_m_alpha_Rot));
                        Fastor::Tensor<double, 3, 3> omega_alpha = 0.5 * (outer(F_m_alpha_Rot, F_n_alpha_Rot) - 
                                                                          outer(F_n_alpha_Rot, F_m_alpha_Rot));
                        
                        double tau_alpha = 0.0;
                        for (int i = 0; i < 3; i++)
                        {
                            for (int j = 0; j < 3; j++)
                            {
                                tau_alpha += Fsigma(i, j) * P_alpha(i, j);
                            }
                        }
                        
                        // Update slip rate with current stress and hardening
                        double k_exponent = MatList[Phase2Material[h]].k_exponent;
                        double ref_shear_rate = MatList[Phase2Material[h]].ref_shear_rate;

//                                   double tau_ratio = mag(tau_alpha)/g_alphas[h][s].boundaryField()[patchI][faceI];
//           double a0 = 1.0/k_exponent;  // This makes sinh argument ~1 at yield
//           double sinh_arg = tau_ratio/a0;
//           if (sinh_arg > 25.0) {  // sinh(20) ≈ 2.4e8
//     // Linear extrapolation beyond this point
//     gamma_alpha_dots[h][s].boundaryFieldRef()[patchI][faceI] = ref_shear_rate * Foam::sinh(25.0) * sign(tau_alpha);
// } else {
//     gamma_alpha_dots[h][s].boundaryFieldRef()[patchI][faceI] = ref_shear_rate * Foam::sinh(sinh_arg) * sign(tau_alpha);
// }


                        gamma_alpha_dots[h][s].boundaryFieldRef()[patchI][faceI] = ref_shear_rate * 
                            Foam::pow(mag(tau_alpha / g_alphas[h][s].boundaryField()[patchI][faceI]), k_exponent) * 
                            sign(tau_alpha);
                        
                        // // Update local plastic deformation rate
                        // gamma_alpha_dots[h][s][cell] = gamma_dot_alpha;
                        F_D_P += P_alpha * gamma_alpha_dots[h][s].boundaryField()[patchI][faceI];//gamma_dot_alpha;
                        F_Omega_P += omega_alpha * gamma_alpha_dots[h][s].boundaryField()[patchI][faceI];//gamma_dot_alpha;
                        
                        // Update hardening (simplified - you may want full hardening law)
                        // g_alphas[h][s][cell] += sub_dt.value() * g_alpha_dots[h][s][cell];
                    }
                    
                    // Store plastic deformation tensors
        
                    //THINK I SHOULD ADD MAIN CONTRIBUTION HERE TOO

                    enum
        {
          I,J,K,L
        };
        Fastor::Tensor<double, 3, 3> dsigma_elastic = Fastor::einsum<Fastor::Index<I, J, K, L>, Fastor::Index<K, L> /*,Fastor::OIndex<I,J>*/>(Cprime, Fepsilon_dot - F_D_P); ////Double contraction
        
            count_index = 0;
            for (int i = 0; i < 3; i++)
            {
                for (int j = 0; j < 3; j++)
                {
                            dsigmadt_local.boundaryFieldRef()[patchI][faceI].component(count_index) = dsigma_elastic(i, j);
                            F_D_P_local.boundaryFieldRef()[patchI][faceI].component(count_index) = F_D_P(i, j);
                            F_Omega_P_local.boundaryFieldRef()[patchI][faceI].component(count_index) = F_Omega_P(i, j);
                            Plastic_Strain[h].boundaryFieldRef()[patchI][faceI].component(count_index) = F_D_P(i, j);
                            Omega_P[h].boundaryFieldRef()[patchI][faceI].component(count_index) = F_Omega_P(i, j);
                            count_index++;
                }
            }



                }










////////////////////////////////////////////////////////////////////////////////////////////////////////////ROtation_update
        Fastor::Tensor<double, 3, 3> Elastic_Spin_h_deltaT(0.0);
        Fastor::Tensor<double, 3, 3> Elastic_Spin_h_deltaT_sq(0.0);
        int count_index_E = 0;
        for (int ii = 0; ii < 3; ii++)
        {
          for (int jj = 0; jj < 3; jj++)
          {

            Elastic_Spin_h_deltaT(ii, jj) = (omega_dot.boundaryField()[patchI][faceI].component(count_index_E) - F_Omega_P_local.boundaryField()[patchI][faceI].component(count_index_E));//

            count_index_E++;
          }
        }

        Elastic_Spin_h_deltaT_sq = Fastor::matmul(Elastic_Spin_h_deltaT,Elastic_Spin_h_deltaT);  



     //Update rotation Matrix Based on Elastic spin

     Fastor::Tensor<double, 3, 3> Rot_matrix_N_plus_one(0.0);
     Fastor::Tensor<double, 3, 3> I_d(0.0);
     for(int set=0;set<3;set++){
      I_d(set,set)=1.0;
     }
      
    //  Fastor::print("I",I_d);

    double X_comps = Foam::sqrt(

                    Foam::pow(Elastic_Spin_h_deltaT(2,1),2.0)
                    +Foam::pow(Elastic_Spin_h_deltaT(0, 2),2.0)
                    +Foam::pow(Elastic_Spin_h_deltaT(1,0),2.0)
                    
                    );

    //                 Fastor::print("Elastic_Spin_h_deltaT",Elastic_Spin_h_deltaT);

    // Info<<"X_comps: "<<X_comps<<endl;
    Fastor::Tensor<double, 3, 3> e_to_the_C(0.0);
    e_to_the_C = I_d;
    if(X_comps>SMALL){
      // Info<<"x_comps"<<X_comps<<endl;
    

    e_to_the_C  +=  ((Foam::sin(X_comps*sub_dt.value())/X_comps)*Elastic_Spin_h_deltaT) + ((1.0-Foam::cos(X_comps*sub_dt.value()))/Foam::pow(X_comps,2.0))* Elastic_Spin_h_deltaT_sq;

//sub_dt
    
    }

     Rot_matrix_N_plus_one = Fastor::matmul(e_to_the_C,Rot_matrix);


     double theta_angle = 0.0;
    double psi_angle = 0.0;
    double phi_angle = 0.0;

    if(mag(Rot_matrix_N_plus_one(2,0))!=1.0){
      // Info<<"HERE"<<endl;
      double theta_1 = -std::asin(Rot_matrix_N_plus_one(2,0));
      double theta_2 = M_PI - theta_1;
      double psi_1 = std::atan2((Rot_matrix_N_plus_one(2,1)/std::cos(theta_1)),(Rot_matrix_N_plus_one(2,2)/std::cos(theta_1)));
      double psi_2 = std::atan2((Rot_matrix_N_plus_one(2,1)/std::cos(theta_2)),(Rot_matrix_N_plus_one(2,2)/std::cos(theta_2)));
      double phi_1 = std::atan2((Rot_matrix_N_plus_one(1,0)/std::cos(theta_1)),(Rot_matrix_N_plus_one(0,0)/std::cos(theta_1)));
      double phi_2 = std::atan2((Rot_matrix_N_plus_one(1,0)/std::cos(theta_2)),(Rot_matrix_N_plus_one(0,0)/std::cos(theta_2)));

      theta_angle = theta_1;//max(theta_1,theta_2);
      psi_angle = psi_1;//max(psi_1,psi_2);
      phi_angle = phi_1;//max(phi_1,phi_2);

      // Info<<"TEST: "<<psi_1<<"\t"<<psi_2<<endl;

    }
    else{
      phi_angle = 0.0;
        if(Rot_matrix_N_plus_one(2,0)==-1.0){
          theta_angle = M_PI/2.0;
          psi_angle = phi_angle + std::atan2(Rot_matrix_N_plus_one(0,1),(Rot_matrix_N_plus_one(0,2)));
        }
        else{
          theta_angle = -M_PI/2.0;
          psi_angle = -phi_angle + std::atan2(-Rot_matrix_N_plus_one(0,1),-Rot_matrix_N_plus_one(0,2));
        }

    }

    Rot1[h].boundaryFieldRef()[patchI][faceI] = psi_angle;//update euler angles for advection
    Rot2[h].boundaryFieldRef()[patchI][faceI] = theta_angle;
    Rot3[h].boundaryFieldRef()[patchI][faceI] = phi_angle;

////////////////////////////////////////////////////////////////////////////////////////////////////////////ROtation_update







        }

    }

}






































            
            // Step 2: Calculate stress rate with updated plastic deformation
            // volTensorField dsigmadt_local = dsigmadt[h];
            
            // Recalculate elastic part of stress rate if needed
            // This is already in dsigmadt[h] from UpdateProps.H
            
            // Step 3: Update stress
            sigma_working += sub_dt * dsigmadt_local;
            
            // Add Jaumann correction terms
            sigma_working -= sub_dt * ((sigma_working & omega_dot) + 
                                       (omega_dot & sigma_working));
            sigma_working += sub_dt * ((sigma_working & F_Omega_P_local) - 
                                       (F_Omega_P_local & sigma_working));
            
            // Update hardening rates for next sub-step
            if (step < n_substeps - 1)  // Not the last step
            {
                // Recalculate hardening rates based on current slip rates
                forAll(mesh.C(), cell)
                {
                    if (PopBal_av[h][cell] > 0.0)
                    {
                        double total_slip = 0.0;
                        for (int s = 0; s < MatList[Phase2Material[h]].N_Slips; s++)
                        {
                            total_slip += mag(gamma_alphas[h][s][cell]);
                        }
                        
                        for (int S_alpha = 0; S_alpha < MatList[Phase2Material[h]].N_Slips; S_alpha++)
                        {
                            g_alpha_dots[h][S_alpha][cell] = 0.0;
                            
                            for (int S_beta = 0; S_beta < MatList[Phase2Material[h]].N_Slips; S_beta++)
                            {
                                double argument = (MatList[Phase2Material[h]].h_0 * total_slip) / 
                                                 (MatList[Phase2Material[h]].tau_S - 
                                                  MatList[Phase2Material[h]].initial_CRSS);
                                double h_of_gamma = MatList[Phase2Material[h]].h_0 * 
                                                   Foam::pow(1.0 / Foam::cosh(argument), 2.0);
                                
                                double h_alpha_beta = (S_beta == S_alpha) ? 
                                    h_of_gamma * (MatList[Phase2Material[h]].r + 
                                                 (1.0 - MatList[Phase2Material[h]].r)) :
                                    h_of_gamma * MatList[Phase2Material[h]].r;
                                
                                g_alpha_dots[h][S_alpha][cell] += h_alpha_beta * 
                                                                  mag(gamma_alpha_dots[h][S_beta][cell]);
                            }

                            double max_g_increment = 0.1 * g_alphas[h][S_alpha][cell];  // 10% change per step
double max_rate_from_increment = (n_substeps*max_g_increment) / runTime.deltaT().value();
g_alpha_dots[h][S_alpha][cell] = sign(g_alpha_dots[h][S_alpha][cell]) * 
    min(mag(g_alpha_dots[h][S_alpha][cell]), max_rate_from_increment);

                        }
                    }
                }
            }

        Rot1[h].correctBoundaryConditions();
        Rot2[h].correctBoundaryConditions();
        Rot3[h].correctBoundaryConditions();



        }//end over sub-stepping
        
        sigtemp[h] = sigma_working;
        sigtemp[h].correctBoundaryConditions();
        sigmaGlobaltemp += PopBal[h] * sigtemp[h];
    }
    
    sigmaGlobaltemp.correctBoundaryConditions();
    
    // Update committed state on final PIMPLE iteration
    if (pimple.finalPimpleIter())
    {
        Info << "Stress Field Updated" << endl;
        for (label h = 0; h < N_Ori; h++)
        {
            sigma[h] = sigtemp[h];
        }
        sigmaGlobal = sigmaGlobaltemp;
        sigmaGlobal.correctBoundaryConditions();
        
        D = D + runTime.deltaT() * U;
        D.correctBoundaryConditions();
        
        p = -(1.0/3.0) * tr(sigmaGlobal);
    }


        for (label h = 0; h < N_Ori; h++)
    {
        Rot1[h].correctBoundaryConditions();
        Rot2[h].correctBoundaryConditions();
        Rot3[h].correctBoundaryConditions();
    }



}