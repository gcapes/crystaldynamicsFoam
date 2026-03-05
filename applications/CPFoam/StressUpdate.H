{
    sigmaGlobaltemp *= 0.0;
    


    for (label h = 0; h < N_Ori; h++)
    {
        Info << "Phase " << h << " - Number of Stress Integration Sub-Steps: " << n_substeps << endl;
        dimensionedScalar sub_dt = runTime.deltaT() / n_substeps;
        
        // Working copy of stress for sub-stepping
        volTensorField sigma_working = sigtemp[h];
        
        // Create working copies of hardening variables for sub-stepping
        PtrList<volScalarField> g_alphas_working(MatList[Phase2Material[h]].N_Slips);
        for (int s = 0; s < MatList[Phase2Material[h]].N_Slips; s++)
        {
            g_alphas_working.set(
                s,
                new volScalarField(g_alphas[h][s])
            );
        }
        
        // Fields to accumulate plastic deformation over sub-steps
        volTensorField accumulated_plastic_strain(
            IOobject(
                "accumulated_plastic_strain",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedTensor("zero", dimensionSet(0, 0, 0, 0, 0), tensor::zero)
        );
        
        volTensorField accumulated_plastic_spin(
            IOobject(
                "accumulated_plastic_spin",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedTensor("zero", dimensionSet(0, 0, 0, 0, 0), tensor::zero)
        );
        
        // Sub-stepping loop for stress integration
        for (int step = 0; step < n_substeps; step++) 
        {
            // Local fields for this sub-step
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
            
            // Step 1: Calculate CP constitutive response for cells
            forAll(mesh.C(), cell)
            {
                if (PopBal_av[h][cell] > 0.0)
                {
                    // Extract current stress and strain rate
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
                    
                    // Use WORKING rotation for this phase
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
                    
                    // Rotate stiffness tensor
                    Fastor::Tensor<double, 3, 3, 3, 3> Cprime(0.0);
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
                                                    value += Rot_matrix(i, m) * Rot_matrix(j, n) * 
                                                            Rot_matrix(k, p) * Rot_matrix(l, q) * 
                                                            Cprime2[Phase2Material[h]](m, n, p, q);
                                                }
                                            }
                                        }
                                    }
                                    Cprime(i, j, k, l) = value;
                                }
                            }
                        }
                    }
                    
                    // Calculate plastic deformation rate
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
                        
                        // Calculate slip rate using current hardening
                        double k_exponent = MatList[Phase2Material[h]].k_exponent;
                        double ref_shear_rate = MatList[Phase2Material[h]].ref_shear_rate;
                        
                        gamma_alpha_dots[h][s][cell] = ref_shear_rate * 
                            Foam::pow(mag(tau_alpha / g_alphas_working[s][cell]), k_exponent) * 
                            sign(tau_alpha);
                        
                        // Accumulate plastic deformation
                        F_D_P += P_alpha * gamma_alpha_dots[h][s][cell];
                        F_Omega_P += omega_alpha * gamma_alpha_dots[h][s][cell];
                    }
                    
                    // Calculate elastic stress rate
                    enum { I, J, K, L };
                    Fastor::Tensor<double, 3, 3> dsigma_elastic = 
                        Fastor::einsum<Fastor::Index<I, J, K, L>, Fastor::Index<K, L>>(
                            Cprime, Fepsilon_dot - F_D_P);
                    
                    // Store stress rate and accumulate plastic deformation
                    count_index = 0;
                    for (int i = 0; i < 3; i++)
                    {
                        for (int j = 0; j < 3; j++)
                        {
                            dsigmadt_local[cell].component(count_index) = dsigma_elastic(i, j);
accumulated_plastic_strain[cell].component(count_index) += 
    (sub_dt * F_D_P(i, j)).value();
accumulated_plastic_spin[cell].component(count_index) += 
    (sub_dt * F_Omega_P(i, j)).value();
                            count_index++;
                        }
                    }
                    
                    // Update hardening evolution
                    double total_slip = 0.0;
                    for (int s = 0; s < MatList[Phase2Material[h]].N_Slips; s++)
                    {
                        total_slip += mag(gamma_alphas[h][s][cell]) + 
                                     sub_dt.value() * mag(gamma_alpha_dots[h][s][cell]);
                    }
                    
                    for (int S_alpha = 0; S_alpha < MatList[Phase2Material[h]].N_Slips; S_alpha++)
                    {
                        double g_dot = 0.0;
                        
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
                            
                            g_dot += h_alpha_beta * mag(gamma_alpha_dots[h][S_beta][cell]);
                        }
                        
                        g_alpha_dots[h][S_alpha][cell] = g_dot;
                        
                        // Limit hardening rate
                        double max_g_increment = 0.1 * g_alphas_working[S_alpha][cell];
                        double max_rate_from_increment = (n_substeps * max_g_increment) / 
                                                        runTime.deltaT().value();
                        g_alpha_dots[h][S_alpha][cell] = sign(g_alpha_dots[h][S_alpha][cell]) * 
                            min(mag(g_alpha_dots[h][S_alpha][cell]), max_rate_from_increment);
                    }
                }
            }
            
            // Boundary condition handling
            forAll(T.boundaryField(), patchI)
            {
                if (!T.boundaryField()[patchI].coupled() && 
                    T.boundaryField()[patchI].size() > 0)
                {
                    forAll(T.boundaryField()[patchI], faceI)
                    {
                        if (PopBal_av[h].boundaryField()[patchI][faceI] > 0.0)
                        {
                            // Extract current stress and strain rate at boundary
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
                            
                            // Use WORKING rotation for boundary
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
                            
                            // Rotate stiffness tensor for boundary
                            Fastor::Tensor<double, 3, 3, 3, 3> Cprime(0.0);
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
                                                            value += Rot_matrix(i, m) * Rot_matrix(j, n) * 
                                                                    Rot_matrix(k, p) * Rot_matrix(l, q) * 
                                                                    Cprime2[Phase2Material[h]](m, n, p, q);
                                                        }
                                                    }
                                                }
                                            }
                                            Cprime(i, j, k, l) = value;
                                        }
                                    }
                                }
                            }
                            
                            // Calculate plastic deformation rate for boundary
                            Fastor::Tensor<double, 3, 3> F_D_P(0.0);
                            Fastor::Tensor<double, 3, 3> F_Omega_P(0.0);
                            
                            for (int s = 0; s < MatList[Phase2Material[h]].N_Slips; s++)
                            {
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
                                
                                double k_exponent = MatList[Phase2Material[h]].k_exponent;
                                double ref_shear_rate = MatList[Phase2Material[h]].ref_shear_rate;
                                
                                gamma_alpha_dots[h][s].boundaryFieldRef()[patchI][faceI] = ref_shear_rate * 
                                    Foam::pow(mag(tau_alpha / g_alphas_working[s].boundaryField()[patchI][faceI]), k_exponent) * 
                                    sign(tau_alpha);
                                
                                F_D_P += P_alpha * gamma_alpha_dots[h][s].boundaryField()[patchI][faceI];
                                F_Omega_P += omega_alpha * gamma_alpha_dots[h][s].boundaryField()[patchI][faceI];
                            }
                            
                            // Calculate elastic stress rate for boundary
                            enum { I, J, K, L };
                            Fastor::Tensor<double, 3, 3> dsigma_elastic = 
                                Fastor::einsum<Fastor::Index<I, J, K, L>, Fastor::Index<K, L>>(
                                    Cprime, Fepsilon_dot - F_D_P);
                            
                            count_index = 0;
                            for (int i = 0; i < 3; i++)
                            {
                                for (int j = 0; j < 3; j++)
                                {
                                    dsigmadt_local.boundaryFieldRef()[patchI][faceI].component(count_index) = dsigma_elastic(i, j);
                                    accumulated_plastic_strain.boundaryFieldRef()[patchI][faceI].component(count_index) += 
    (sub_dt * F_D_P(i, j)).value();
accumulated_plastic_spin.boundaryFieldRef()[patchI][faceI].component(count_index) += 
    (sub_dt * F_Omega_P(i, j)).value();
                                    count_index++;
                                }
                            }
                        }
                    }
                }
            }
            
            // Step 2: Update stress for this sub-step
            // volTensorField omega_p_rate = accumulated_plastic_spin / (sub_dt.value() * (step + 1));
            volTensorField omega_p_rate = accumulated_plastic_spin / (sub_dt);
            omega_p_rate.correctBoundaryConditions();
            
            // sigma_working += sub_dt * dsigmadt_local;
            // sigma_working -= sub_dt * ((sigma_working & omega_dot) - (omega_dot & sigma_working));
            // sigma_working += sub_dt * ((sigma_working & omega_p_rate) - (omega_p_rate & sigma_working));

            volTensorField sigma_old = sigma_working;
            sigma_working = sigma_old
            + sub_dt * dsigmadt_local
            + sub_dt * ((omega_dot & sigma_old) - (sigma_old & omega_dot))
            + sub_dt * ((sigma_old & omega_p_rate) - (omega_p_rate & sigma_old));
            
            sigma_working.correctBoundaryConditions();

            // Step 3: Update hardening values for next sub-step
            for (int s = 0; s < MatList[Phase2Material[h]].N_Slips; s++)
            {
                g_alphas_working[s] += sub_dt * g_alpha_dots[h][s];
                g_alphas_working[s].correctBoundaryConditions();  // 
                gamma_alphas[h][s] += sub_dt * gamma_alpha_dots[h][s];
                gamma_alphas[h][s].correctBoundaryConditions();  //
            }
        } // End sub-stepping loop
        
        // Store converged stress
        sigtemp[h] = sigma_working;
        sigtemp[h].correctBoundaryConditions();
        
        // Update actual hardening fields with converged values
        for (label s = 0; s < MatList[Phase2Material[h]].N_Slips; s++)
        {
            g_alphas[h][s] = g_alphas_working[s];
            g_alphas[h][s].correctBoundaryConditions();  //
        }
        
        // Store plastic deformation rates (average over time step)
        Plastic_Strain[h] = accumulated_plastic_strain / runTime.deltaT();
        Omega_P[h] = accumulated_plastic_spin / runTime.deltaT();
        
        // ROTATION UPDATE - ONCE per time step, updates WORKING copies
        // Process cells
        forAll(mesh.C(), cell)
        {
            // if (PopBal_av[h][cell] > 0.0)
            // {
                // Calculate total elastic spin over full time step
                Fastor::Tensor<double, 3, 3> Elastic_Spin(0.0);
                int count_index_E = 0;
                for (int ii = 0; ii < 3; ii++)
                {
                    for (int jj = 0; jj < 3; jj++)
                    {
                        double total_plastic_spin = accumulated_plastic_spin[cell].component(count_index_E);
                        Elastic_Spin(ii, jj) = (omega_dot[cell].component(count_index_E)) - Omega_P[h][cell].component(count_index_E);
                    //   total_plastic_spin;
                        count_index_E++;
                    }
                }
                
                // Calculate rotation update using exponential map
                Fastor::Tensor<double, 3, 3> Elastic_Spin_sq = Fastor::matmul(Elastic_Spin, Elastic_Spin);
                Fastor::Tensor<double, 3, 3> I_d(0.0);
                for(int set = 0; set < 3; set++) {
                    I_d(set, set) = 1.0;
                }
                
                double X_comps = Foam::sqrt(
                    Foam::pow(Elastic_Spin(2,1), 2.0) +
                    Foam::pow(Elastic_Spin(0,2), 2.0) +
                    Foam::pow(Elastic_Spin(1,0), 2.0));
                //* runTime.deltaT().value()

                // Info<<X_comps<<endl;
                Fastor::Tensor<double, 3, 3> e_to_the_C = I_d;
                if (X_comps > SMALL) {
                    e_to_the_C += (Foam::sin(X_comps* runTime.deltaT().value()) / X_comps) * Elastic_Spin + 
                                 ((1.0 - Foam::cos(X_comps* runTime.deltaT().value())) / Foam::pow(X_comps, 2.0)) * Elastic_Spin_sq;
                }
                
                // Get current WORKING rotation matrix
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
                
                // Update rotation
                Fastor::Tensor<double, 3, 3> Rot_matrix_N_plus_one = Fastor::matmul(e_to_the_C, Rot_matrix);
                
                // Extract Euler angles from updated rotation matrix
                double theta_angle = 0.0;
                double psi_angle = 0.0;
                double phi_angle = 0.0;
                
                if (mag(Rot_matrix_N_plus_one(2,0)) != 1.0) {
                    double theta_1 = -std::asin(Rot_matrix_N_plus_one(2,0));
                    double theta_2 = M_PI - theta_1;
                    double psi_1 = std::atan2(Rot_matrix_N_plus_one(2,1)/std::cos(theta_1),
                                              Rot_matrix_N_plus_one(2,2)/std::cos(theta_1));
                    double psi_2 = std::atan2(Rot_matrix_N_plus_one(2,1)/std::cos(theta_2),
                                              Rot_matrix_N_plus_one(2,2)/std::cos(theta_2));
                    double phi_1 = std::atan2(Rot_matrix_N_plus_one(1,0)/std::cos(theta_1),
                                              Rot_matrix_N_plus_one(0,0)/std::cos(theta_1));
                    double phi_2 = std::atan2(Rot_matrix_N_plus_one(1,0)/std::cos(theta_2),
                                              Rot_matrix_N_plus_one(0,0)/std::cos(theta_2));
                    
                    theta_angle = theta_1;
                    psi_angle = psi_1;
                    phi_angle = phi_1;
                } else {
                    phi_angle = 0.0;
                    if (Rot_matrix_N_plus_one(2,0) == -1.0) {
                        theta_angle = M_PI/2.0;
                        psi_angle = phi_angle + std::atan2(Rot_matrix_N_plus_one(0,1),
                                                           Rot_matrix_N_plus_one(0,2));
                    } else {
                        theta_angle = -M_PI/2.0;
                        psi_angle = -phi_angle + std::atan2(-Rot_matrix_N_plus_one(0,1),
                                                            -Rot_matrix_N_plus_one(0,2));
                    }
                }
                
                // Update WORKING rotation fields
                Rot1[h][cell] = psi_angle;
                Rot2[h][cell] = theta_angle;
                Rot3[h][cell] = phi_angle;
            // }
        }
        
        // Process boundaries for rotation update
        forAll(T.boundaryField(), patchI)
        {
            if (!T.boundaryField()[patchI].coupled() && 
                T.boundaryField()[patchI].size() > 0)
            {
                forAll(T.boundaryField()[patchI], faceI)
                {
                    // if (PopBal_av[h].boundaryField()[patchI][faceI] > 0.0)
                    // {
                        // Calculate total elastic spin over full time step for boundary
                        Fastor::Tensor<double, 3, 3> Elastic_Spin(0.0);
                        int count_index_E = 0;
                        for (int ii = 0; ii < 3; ii++)
                        {
                            for (int jj = 0; jj < 3; jj++)
                            {
                                double total_plastic_spin = accumulated_plastic_spin.boundaryField()[patchI][faceI].component(count_index_E);
                                Elastic_Spin(ii, jj) = (omega_dot.boundaryField()[patchI][faceI].component(count_index_E) ) - Omega_P[h].boundaryField()[patchI][faceI].component(count_index_E);
                    //   total_plastic_spin;
                                count_index_E++;
                            }
                        }
                        
                        // Calculate rotation update using exponential map
                        Fastor::Tensor<double, 3, 3> Elastic_Spin_sq = Fastor::matmul(Elastic_Spin, Elastic_Spin);
                        Fastor::Tensor<double, 3, 3> I_d(0.0);
                        for(int set = 0; set < 3; set++) {
                            I_d(set, set) = 1.0;
                        }
                        
                        double X_comps = Foam::sqrt(
                            Foam::pow(Elastic_Spin(2,1), 2.0) +
                            Foam::pow(Elastic_Spin(0,2), 2.0) +
                            Foam::pow(Elastic_Spin(1,0), 2.0));
                        
                        Fastor::Tensor<double, 3, 3> e_to_the_C = I_d;
                        if (X_comps > SMALL) {
                            e_to_the_C += (Foam::sin(X_comps* runTime.deltaT().value()) / X_comps) * Elastic_Spin + 
                                         ((1.0 - Foam::cos(X_comps* runTime.deltaT().value())) / Foam::pow(X_comps, 2.0)) * Elastic_Spin_sq;
                        }
                        
                        // Get current WORKING rotation matrix for boundary
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
                        
                        // Update rotation
                        Fastor::Tensor<double, 3, 3> Rot_matrix_N_plus_one = Fastor::matmul(e_to_the_C, Rot_matrix);
                        
                        // Extract Euler angles from updated rotation matrix
                        double theta_angle = 0.0;
                        double psi_angle = 0.0;
                        double phi_angle = 0.0;
                        
                        if (mag(Rot_matrix_N_plus_one(2,0)) != 1.0) {
                            double theta_1 = -std::asin(Rot_matrix_N_plus_one(2,0));
                            double theta_2 = M_PI - theta_1;
                            double psi_1 = std::atan2(Rot_matrix_N_plus_one(2,1)/std::cos(theta_1),
                                                      Rot_matrix_N_plus_one(2,2)/std::cos(theta_1));
                            double psi_2 = std::atan2(Rot_matrix_N_plus_one(2,1)/std::cos(theta_2),
                                                      Rot_matrix_N_plus_one(2,2)/std::cos(theta_2));
                            double phi_1 = std::atan2(Rot_matrix_N_plus_one(1,0)/std::cos(theta_1),
                                                      Rot_matrix_N_plus_one(0,0)/std::cos(theta_1));
                            double phi_2 = std::atan2(Rot_matrix_N_plus_one(1,0)/std::cos(theta_2),
                                                      Rot_matrix_N_plus_one(0,0)/std::cos(theta_2));
                            
                            theta_angle = theta_1;
                            psi_angle = psi_1;
                            phi_angle = phi_1;
                        } else {
                            phi_angle = 0.0;
                            if (Rot_matrix_N_plus_one(2,0) == -1.0) {
                                theta_angle = M_PI/2.0;
                                psi_angle = phi_angle + std::atan2(Rot_matrix_N_plus_one(0,1),
                                                                   Rot_matrix_N_plus_one(0,2));
                            } else {
                                theta_angle = -M_PI/2.0;
                                psi_angle = -phi_angle + std::atan2(-Rot_matrix_N_plus_one(0,1),
                                                                    -Rot_matrix_N_plus_one(0,2));
                            }
                        }
                        
                        // Update WORKING rotation fields for boundary
                        Rot1[h].boundaryFieldRef()[patchI][faceI] = psi_angle;
                        Rot2[h].boundaryFieldRef()[patchI][faceI] = theta_angle;
                        Rot3[h].boundaryFieldRef()[patchI][faceI] = phi_angle;
                    // }
                }
            }
        }
        
        // Correct boundary conditions for working rotations
        Rot1[h].correctBoundaryConditions();
        Rot2[h].correctBoundaryConditions();
        Rot3[h].correctBoundaryConditions();
        
        // Accumulate global stress
        sigmaGlobaltemp += PopBal[h] * sigtemp[h];
        
    } // End loop over phases
    
    sigmaGlobaltemp.correctBoundaryConditions();
    
    // Update committed state ONLY on final PIMPLE iteration
    if (pimple.finalIter())
    {
        Info << "Updating committed fields on final PIMPLE iteration" << endl;
        for (label h = 0; h < N_Ori; h++)
        {
            // Commit stress with relaxation
            sigma[h] = sigtemp[h];



            sigma[h].correctBoundaryConditions();
            
            // Commit rotations
            // Rot1[h] = Rot1_working[h];
            // Rot2[h] = Rot2_working[h];
            // Rot3[h] = Rot3_working[h];

            Rot1[h].correctBoundaryConditions();
            Rot2[h].correctBoundaryConditions();
            Rot3[h].correctBoundaryConditions();

        }
        sigmaGlobal = sigmaGlobaltemp;
        sigmaGlobal.correctBoundaryConditions();
        
        D = D + runTime.deltaT() * U;
        D.correctBoundaryConditions();
        
        p = -(1.0/3.0) * tr(sigmaGlobal);
    }
}