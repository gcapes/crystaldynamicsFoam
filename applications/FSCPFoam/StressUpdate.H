/*---------------------------------------------------------------------------*\
  StressUpdate.H - Finite Strain Crystal Plasticity Constitutive Update
  
  ALGORITHM (per cell, per phase):
  ================================
  1. Extract F (global) and Fp (plastic, per-phase)
  2. Elastic part: Fe = F · inv(Fp)
  3. Polar decomposition: Fe = R · U  (R = rotation, U = stretch)
  4. Elastic Green strain: Ee = 0.5*(Fe^T · Fe - I)
  5. Rotate stiffness: C_rotated = R x R x R x R : C_crystal
  6. 2nd Piola-Kirchhoff stress: S = C_rotated : Ee
  7. For each slip system:
     a. Rotate slip system to intermediate config: m_rot = R.m, n_rot = R.n
     b. Schmid tensor: P = m_rot x n_rot 
     c. Resolved shear stress: tau = S : P
     d. Slip rate: gamma_dot = gamma_dot_0 * |tau/g|^k * sign(tau)
  8. Plastic velocity gradient: Lp = sum( gamma_dot * P )
  9. Sub-step Fp: Fp_new = Fp + dt_sub * Lp . Fp
  10. Cauchy stress: sigma = (1/Je) * Fe . S . Fe^T
  
\*---------------------------------------------------------------------------*/
{
    sigmaGlobaltemp *= 0.0;
    rho_ref *= 0.0;  // Will be rebuilt from phases

    for (label h = 0; h < N_Ori; h++)
    {
        Info << "Phase " << h << " - Finite strain constitutive update ("
             << n_substeps << " sub-steps)" << endl;

        dimensionedScalar sub_dt = runTime.deltaT() / n_substeps;

        // Working copy of Fp for sub-stepping (don't modify the advected field)
        volTensorField Fp_work = F_plastic[h];

        // Working copies of hardening variables
        PtrList<volScalarField> g_alphas_working(MatList[Phase2Material[h]].N_Slips);
        for (int s = 0; s < MatList[Phase2Material[h]].N_Slips; s++)
        {
            g_alphas_working.set(s, new volScalarField(g_alphas[h][s]));
        }

        // ============ Sub-stepping loop ============
        for (int step = 0; step < n_substeps; step++)
        {
            // ---------- Cell-center constitutive update ----------
            forAll(mesh.C(), cell)
            {
                // --- Extract F and Fp as Fastor tensors ---
                Fastor::Tensor<double, 3, 3> F_f = Foam::toFastor(F[cell]);
                Fastor::Tensor<double, 3, 3> Fp_f = Foam::toFastor(Fp_work[cell]);

                // --- Fe = F . inv(Fp) ---
                Fastor::Tensor<double, 3, 3> Fp_inv = Fastor::inverse(Fp_f);
                Fastor::Tensor<double, 3, 3> Fe_f = Fastor::matmul(F_f, Fp_inv);

                // --- Polar decomposition: Fe = R . U ---
                Fastor::Tensor<double, 3, 3> R(0.0), U_stretch(0.0);
                polarDecomposition(Fe_f, R, U_stretch);

                // --- Elastic Green strain: Ee = 0.5*(Fe^T.Fe - I) ---
                Fastor::Tensor<double, 3, 3> Ce = Fastor::matmul(
                    Fastor::transpose(Fe_f), Fe_f);
                Fastor::Tensor<double, 3, 3> Ee(0.0);
                for (int i = 0; i < 3; i++)
                    for (int j = 0; j < 3; j++)
                        Ee(i, j) = 0.5 * (Ce(i, j) - (i == j ? 1.0 : 0.0));

                // --- Rotate stiffness tensor ---
                Fastor::Tensor<double, 3, 3, 3, 3> Cprime(0.0);
                for (int i = 0; i < 3; i++)
                    for (int j = 0; j < 3; j++)
                        for (int k = 0; k < 3; k++)
                            for (int l = 0; l < 3; l++)
                            {
                                double value = 0.0;
                                for (int m = 0; m < 3; m++)
                                    for (int n = 0; n < 3; n++)
                                        for (int p = 0; p < 3; p++)
                                            for (int q = 0; q < 3; q++)
                                                value += R(i,m)*R(j,n)*R(k,p)*R(l,q)
                                                       * Cprime2[Phase2Material[h]](m,n,p,q);
                                Cprime(i, j, k, l) = value;
                            }

                // --- 2nd Piola-Kirchhoff stress: S = C : Ee ---
                enum { I, J, K, L };
                Fastor::Tensor<double, 3, 3> S =
                    Fastor::einsum<Fastor::Index<I,J,K,L>, Fastor::Index<K,L>>(Cprime, Ee);

                // --- Slip rates and plastic velocity gradient ---
                Fastor::Tensor<double, 3, 3> Lp_local(0.0);

                for (int s = 0; s < MatList[Phase2Material[h]].N_Slips; s++)
                {
                    // Slip system in reference config
                    Fastor::Tensor<double, 3> n_alpha = Foam::toFastor(
                        MatList[Phase2Material[h]].SP_norm[s]);
                    Fastor::Tensor<double, 3> m_alpha = Foam::toFastor(
                        MatList[Phase2Material[h]].SP_direction[s]);

                    // Rotate to intermediate config using R from polar decomp
                    Fastor::Tensor<double, 3> n_rot(0.0), m_rot(0.0);
                    for (int i = 0; i < 3; i++)
                        for (int j = 0; j < 3; j++)
                        {
                            n_rot(i) += R(i, j) * n_alpha(j);
                            m_rot(i) += R(i, j) * m_alpha(j);
                        }

                    // Schmid tensor: P = m x n (unsymmetrized for Lp)
                    Fastor::Tensor<double, 3, 3> P_alpha = Fastor::outer(m_rot, n_rot);

                    // Resolved shear stress: tau = S : P
                    double tau_alpha = 0.0;
                    for (int i = 0; i < 3; i++)
                        for (int j = 0; j < 3; j++)
                            tau_alpha += S(i, j) * P_alpha(i, j);

                    // Power-law slip rate
                    double k_exp = MatList[Phase2Material[h]].k_exponent;
                    double gamma_dot_0 = MatList[Phase2Material[h]].ref_shear_rate;
                    double g_a = g_alphas_working[s][cell];

                    gamma_alpha_dots[h][s][cell] = gamma_dot_0 
                        * Foam::pow(mag(tau_alpha / g_a), k_exp) 
                        * sign(tau_alpha);

                    // Accumulate Lp
                    Lp_local += P_alpha * gamma_alpha_dots[h][s][cell];
                }

                // --- Store Lp for advection equation ---
                Lp[h][cell] = Foam::toOpenFOAM(Lp_local);

                // --- Sub-step Fp: Fp_new = Fp + dt_sub * Lp . Fp ---
                Fastor::Tensor<double, 3, 3> Fp_dot = Fastor::matmul(Lp_local, Fp_f);
                Fp_f += Fp_dot * sub_dt.value();

                // --- Push forward to Cauchy: sigma = (1/Je) Fe.S.Fe^T ---
                double J_e = Fastor::determinant(Fe_f);
                Fastor::Tensor<double, 3, 3> sigma_f =
                    (1.0 / J_e) * Fastor::matmul(
                        Fastor::matmul(Fe_f, S), 
                        Fastor::transpose(Fe_f));

                // --- Store results ---
                Fp_work[cell] = Foam::toOpenFOAM(Fp_f);
                sigtemp[h][cell] = Foam::toOpenFOAM(sigma_f);

                // --- Update hardening ---
                double total_slip = 0.0;
                for (int s = 0; s < MatList[Phase2Material[h]].N_Slips; s++)
                {
                    total_slip += mag(gamma_alphas[h][s][cell])
                                + sub_dt.value() * mag(gamma_alpha_dots[h][s][cell]);
                }

                for (int sa = 0; sa < MatList[Phase2Material[h]].N_Slips; sa++)
                {
                    double g_dot = 0.0;
                    for (int sb = 0; sb < MatList[Phase2Material[h]].N_Slips; sb++)
                    {
                        double argument = (MatList[Phase2Material[h]].h_0 * total_slip)
                            / (MatList[Phase2Material[h]].tau_S 
                             - MatList[Phase2Material[h]].initial_CRSS);
                        double h_of_gamma = MatList[Phase2Material[h]].h_0 
                            * Foam::pow(1.0 / Foam::cosh(argument), 2.0);

                        double h_ab = (sb == sa) 
                            ? h_of_gamma * (MatList[Phase2Material[h]].r 
                                          + (1.0 - MatList[Phase2Material[h]].r))
                            : h_of_gamma * MatList[Phase2Material[h]].r;

                        g_dot += h_ab * mag(gamma_alpha_dots[h][sb][cell]);
                    }

                    g_alpha_dots[h][sa][cell] = g_dot;

                    // Limit hardening rate to avoid divergence
                    double max_g_inc = 0.1 * g_alphas_working[sa][cell];
                    double max_rate = (n_substeps * max_g_inc) / runTime.deltaT().value();
                    g_alpha_dots[h][sa][cell] = sign(g_alpha_dots[h][sa][cell])
                        * min(mag(g_alpha_dots[h][sa][cell]), max_rate);
                }

            } // End cell loop


            // ---------- Increment hardening for next sub-step (cells only) ----------
            // Field-level += operates on cells; boundaries set by correctBoundaryConditions
            for (int s = 0; s < MatList[Phase2Material[h]].N_Slips; s++)
            {
                g_alphas_working[s] += sub_dt * g_alpha_dots[h][s];
                gamma_alphas[h][s] += sub_dt * gamma_alpha_dots[h][s];
            }

            // ---------- Propagate cell values to boundaries ----------
            // zeroGradient BC: copies nearest cell center value to boundary face
            // This guarantees cell-boundary consistency for ALL fields
            Fp_work.correctBoundaryConditions();
            Lp[h].correctBoundaryConditions();
            sigtemp[h].correctBoundaryConditions();

            for (int s = 0; s < MatList[Phase2Material[h]].N_Slips; s++)
            {
                g_alphas_working[s].correctBoundaryConditions();
                gamma_alphas[h][s].correctBoundaryConditions();
                gamma_alpha_dots[h][s].correctBoundaryConditions();
                g_alpha_dots[h][s].correctBoundaryConditions();
            }

        } // End sub-stepping loop

        // ============ Commit converged values for this phase ============
        F_plastic[h] = Fp_work;
        F_plastic[h].correctBoundaryConditions();

        sigma[h] = sigtemp[h];
        sigma[h].correctBoundaryConditions();

        for (int s = 0; s < MatList[Phase2Material[h]].N_Slips; s++)
        {
            g_alphas[h][s] = g_alphas_working[s];
            g_alphas[h][s].correctBoundaryConditions();
        }

        // Accumulate phase-weighted reference density
        dimensionedScalar dimrho("dimrho", dimensionSet(1, -3, 0, 0, 0), 1.0);
        rho_ref += PopBal[h] * MatList[Phase2Material[h]].rho_0 * dimrho;

        // Accumulate global (mixture) stress
        sigmaGlobaltemp += PopBal[h] * sigtemp[h];

    } // End phase loop

    sigmaGlobaltemp.correctBoundaryConditions();
rho_ref.correctBoundaryConditions();

    // ============ Final updates (on last PIMPLE iteration) ============
    if (pimple.finalIter())
    {
        sigmaGlobal = sigmaGlobaltemp;
        sigmaGlobal.correctBoundaryConditions();

        // Jacobian from global deformation gradient
        J = det(F);
        J.correctBoundaryConditions();

        // Density
        rho = rho_ref / J;
        rho.correctBoundaryConditions();

        // Accumulated displacement
        D = D + runTime.deltaT() * U;
        D.correctBoundaryConditions();

        // Pressure (negative mean stress)
        p = -(1.0 / 3.0) * tr(sigmaGlobal);
        p.correctBoundaryConditions();
    }
}