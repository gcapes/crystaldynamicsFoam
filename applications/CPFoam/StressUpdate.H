{
    if (!stressUpdate)
    {
        // Prescribed-stress mode: bypass the full CP/momentum solve.
        // TWO sub-modes, which MUST NOT be mixed:
        //
        //  * PerPhaseSigma false (legacy): ONE prescribed field,
        //    sigmaGlobal, is copied into every solid phase each
        //    corrector; sigma[crack/void] = 0.  Driving-force contrast
        //    comes from per-phase orientation (C_rotated) under the
        //    shared stress.
        //
        //  * PerPhaseSigma true: sigma.<i> are CONSTANT per-phase
        //    labels, stamped ONCE by InitPrescribedStressMode.H and
        //    held for the whole run.  The copy below must therefore be
        //    SKIPPED -- running it wiped the labels on the first
        //    corrector (sigmaGlobal is zero in such cases) and killed
        //    the mechanical driving force entirely.
        if (!perPhaseSigma)
        {
        for (label h = 0; h < N_Ori; h++)
        {
            if (Phase2Material[h] == 0)
            {
                sigma[h] = sigmaGlobal;
                sigma[h].correctBoundaryConditions();
            }
            else
            {
                sigma[h] = dimensionedTensor(
                    "zero", sigma[h].dimensions(), tensor::zero);
                sigma[h].correctBoundaryConditions();
            }
        }
        }

        // Zero all flow fields so that rhoEqn, TEqn, and Advect.H see no
        // advection.  The tensile U boundary conditions (+/-1 m/s) remain active via correctBoundaryConditions
        phi    *= 0.0;
        rhoPhi *= 0.0;
    }
    else
    {

    sigmaGlobaltemp *= 0.0;

    // Fastor index helpers for einsum (C:epsilon contraction)
    enum FastorIdx { FI = 0, FJ = 1, FK = 2, FL = 3 };

    // Maximum slip systems (stack-allocated arrays)
    const int MAX_SLIPS = 48;

    for (label h = 0; h < N_Ori; h++)
    {
        // --- Cache material properties for this phase ---
        const Material& mat = MatList[Phase2Material[h]];
        const int N_Slips_h  = mat.N_Slips;
        const double k_exp   = mat.k_exponent;
        const double ref_shear = mat.ref_shear_rate;
        const double h_0     = mat.h_0;
        const double tau_S   = mat.tau_S;
        const double tau_0   = mat.initial_CRSS;
        const double r_param = mat.r;
        const double dt_val  = runTime.deltaT().value();
        const double sub_dt_val = dt_val / n_substeps;

        Info << "Phase " << h << " - Sub-Steps: " << n_substeps << endl;

        // =============================================================
        //  Constitutive integration kernel (shared by cells & boundary)
        //
        //  Inputs:  sigma_in[9], epsilon_in[9], omega_in[9],
        //           rot_angles[3], g_in[N_Slips], gamma_in[N_Slips],
        //           phase_active (bool)
        //
        //  Outputs: sigma_out[9], Dp_out[9], Wp_out[9],
        //           gdot_out[N_Slips], gamma_dot_out[N_Slips]
        // =============================================================
        auto constitutiveKernel = [&](
            const double* sigma_in,
            const double* epsilon_in,
            const double* omega_in,
            const double* rot_angles,
            const double* g_in,
            const double* gamma_in,
            bool phase_active,
            double* sigma_out,
            double* Dp_out,
            double* Wp_out,
            double* gdot_out,
            double* gamma_dot_out,
    Fastor::Tensor<double,3,3,3,3>* C_rot_out)
        {
            // Passthrough if phase not active
            if (!phase_active)
            {
                for (int c = 0; c < 9; c++)
                {
                    sigma_out[c] = sigma_in[c];
                    Dp_out[c] = 0.0;
                    Wp_out[c] = 0.0;
                }
                for (int s = 0; s < N_Slips_h; s++)
                {
                    gdot_out[s] = 0.0;
                    gamma_dot_out[s] = 0.0;
                }
                // Keep the stored rotated stiffness defined for inactive
                // cells so downstream consumers (the phase-field mechanical
                // driving force) never read stale or uninitialised memory.
                if (C_rot_out) (*C_rot_out).zeros();
                return;
            }

            // --- Extract strain rate & spin ---
            Fastor::Tensor<double,3,3> Fomega(0.0), Fepsilon(0.0);
            {
                int ci = 0;
                for (int i = 0; i < 3; i++)
                    for (int j = 0; j < 3; j++)
                    {
                        Fomega(i,j)   = omega_in[ci];
                        Fepsilon(i,j) = epsilon_in[ci];
                        ci++;
                    }
            }

            // --- Rotation matrix (ONCE) ---
            const double p1 = rot_angles[0], p2 = rot_angles[1], p3 = rot_angles[2];
            const double c1 = Foam::cos(p1), s1 = Foam::sin(p1);
            const double c2 = Foam::cos(p2), s2 = Foam::sin(p2);
            const double c3 = Foam::cos(p3), s3 = Foam::sin(p3);

            Fastor::Tensor<double,3,3> R(0.0);
            R(0,0) =  c2*c3;
            R(0,1) =  s1*s2*c3 - c1*s3;
            R(0,2) =  c1*s2*c3 + s1*s3;
            R(1,0) =  c2*s3;
            R(1,1) =  s1*s2*s3 + c1*c3;
            R(1,2) =  c1*s2*s3 - s1*c3;
            R(2,0) = -s2;
            R(2,1) =  s1*c2;
            R(2,2) =  c1*c2;

            // --- Rotate stiffness tensor (ONCE) ---
            Fastor::Tensor<double,3,3,3,3> C_rot(0.0);
            for (int i = 0; i < 3; i++)
                for (int j = 0; j < 3; j++)
                    for (int k = 0; k < 3; k++)
                        for (int l = 0; l < 3; l++)
                        {
                            double val = 0.0;
                            for (int m = 0; m < 3; m++)
                                for (int n = 0; n < 3; n++)
                                    for (int p = 0; p < 3; p++)
                                        for (int q = 0; q < 3; q++)
                                            val += R(i,m)*R(j,n)*R(k,p)*R(l,q)
                                                 * Cprime2[Phase2Material[h]](m,n,p,q);
                            C_rot(i,j,k,l) = val;
                        }

            // --- Pre-compute rotated Schmid & spin tensors (ONCE) ---

            // C_rotated[h][cellI] = C_rot;
            if (C_rot_out) *C_rot_out = C_rot;



            Fastor::Tensor<double,3,3> P_a[MAX_SLIPS];
            Fastor::Tensor<double,3,3> W_a[MAX_SLIPS];

            for (int s = 0; s < N_Slips_h; s++)
            {
                Fastor::Tensor<double,3> n_cr(0.0), m_cr(0.0);
                for (int d = 0; d < 3; d++)
                {
                    n_cr(d) = mat.SP_norm[s].component(d);
                    m_cr(d) = mat.SP_direction[s].component(d);
                }

                Fastor::Tensor<double,3> n_gl(0.0), m_gl(0.0);
                for (int i = 0; i < 3; i++)
                    for (int j = 0; j < 3; j++)
                    {
                        n_gl(i) += R(i,j) * n_cr(j);
                        m_gl(i) += R(i,j) * m_cr(j);
                    }

                auto mn = Fastor::outer(m_gl, n_gl);
                auto nm = Fastor::outer(n_gl, m_gl);
                P_a[s] = 0.5 * (mn + nm);
                W_a[s] = 0.5 * (mn - nm);
            }

            // --- Initialise working state ---
            Fastor::Tensor<double,3,3> Fsigma(0.0);
            {
                int ci = 0;
                for (int i = 0; i < 3; i++)
                    for (int j = 0; j < 3; j++)
                    {
                        Fsigma(i,j) = sigma_in[ci];
                        ci++;
                    }
            }

            // EOS mode: the kernel evolves the per-phase DEVIATOR only
            // (the pressure channel is owned by EOS.H).  Enforce the
            // invariant on entry: deviatorises any user-supplied initial
            // stress on the first visit, a no-op thereafter since every
            // increment below is then traceless.
            if (useEOS)
            {
                const double trsEOS =
                    (Fsigma(0,0) + Fsigma(1,1) + Fsigma(2,2))/3.0;
                Fsigma(0,0) -= trsEOS;
                Fsigma(1,1) -= trsEOS;
                Fsigma(2,2) -= trsEOS;
            }

            double local_g[MAX_SLIPS], local_gamma[MAX_SLIPS];
            for (int s = 0; s < N_Slips_h; s++)
            {
                local_g[s]     = g_in[s];
                local_gamma[s] = gamma_in[s];
            }

            double gamma_dot_sum[MAX_SLIPS] = {0.0};
            double g_dot_sum[MAX_SLIPS]     = {0.0};
            Fastor::Tensor<double,3,3> D_P_sum(0.0), Omega_P_sum(0.0);

            // --- Sub-stepping loop ---
            for (int step = 0; step < n_substeps; step++)
            {
                // Slip rates & plastic deformation
                Fastor::Tensor<double,3,3> F_D_P(0.0), F_Omega_P(0.0);
                double gdot_step[MAX_SLIPS] = {0.0};

                for (int s = 0; s < N_Slips_h; s++)
                {
                    double tau = 0.0;
                    for (int i = 0; i < 3; i++)
                        for (int j = 0; j < 3; j++)
                            tau += Fsigma(i,j) * P_a[s](i,j);

                    gdot_step[s] = ref_shear
                                 * Foam::pow(mag(tau / local_g[s]), k_exp)
                                 * Foam::sign(tau);

                    F_D_P     += P_a[s] * gdot_step[s];
                    F_Omega_P += W_a[s] * gdot_step[s];
                }

                for (int s = 0; s < N_Slips_h; s++)
                    gamma_dot_sum[s] += gdot_step[s];
                D_P_sum     += F_D_P;
                Omega_P_sum += F_Omega_P;

                // Elastic stress rate  C : (D - D_p)
                auto dsigma = Fastor::einsum<Fastor::Index<FI,FJ,FK,FL>,
                                             Fastor::Index<FK,FL>>(
                                  C_rot, Fepsilon - F_D_P);

                // EOS mode: deviatoric projection of the stiffness action.
                // The isotropic part (K tr(D-Dp) for cubic symmetry) is
                // removed here and supplied instead by the Mie-Gruneisen
                // p(rho,e) in EOS.H.  The subsequent -tr(De)*Fsigma term
                // is retained: it is exactly the deviatoric part of the
                // legacy volumetric-Jaumann equation, so at small
                // compression the deviator evolution is IDENTICAL to the
                // linearElastic path.
                if (useEOS)
                {
                    const double trdsEOS =
                        (dsigma(0,0) + dsigma(1,1) + dsigma(2,2))/3.0;
                    dsigma(0,0) -= trdsEOS;
                    dsigma(1,1) -= trdsEOS;
                    dsigma(2,2) -= trdsEOS;
                }

                // Jaumann update (elastic spin = W - W_p)
                auto W_e   = Fomega - F_Omega_P;
                auto WeSig = Fastor::matmul(W_e, Fsigma);
                auto SigWe = Fastor::matmul(Fsigma, W_e);

                // Fsigma = Fsigma + sub_dt_val * (dsigma + WeSig - SigWe);

                const double trDe = (Fepsilon(0,0) - F_D_P(0,0))
                                  + (Fepsilon(1,1) - F_D_P(1,1))
                                  + (Fepsilon(2,2) - F_D_P(2,2));

                Fsigma = Fsigma + sub_dt_val * (dsigma + WeSig - SigWe - trDe * Fsigma);


                // Hardening
                double total_slip = 0.0;
                for (int s = 0; s < N_Slips_h; s++)
                    total_slip += std::abs(local_gamma[s]);

                // Voce saturation guard
                double h_gamma = 0.0;
                if (h_0 > 0.0 && (tau_S - tau_0) > SMALL)
                {
                    double argument = std::min(
                        (h_0 * total_slip) / (tau_S - tau_0), 500.0);
                    h_gamma = h_0
                        * Foam::pow(1.0 / Foam::cosh(argument), 2.0);
                }

                for (int sa = 0; sa < N_Slips_h; sa++)
                {
                    double g_dot = 0.0;
                    for (int sb = 0; sb < N_Slips_h; sb++)
                    {
                        double h_ab = (sb == sa) ? h_gamma
                                                 : h_gamma * r_param;
                        g_dot += h_ab * std::abs(gdot_step[sb]);
                    }

                    g_dot = std::min(g_dot,
                                    0.1 * local_g[sa] / sub_dt_val);

                    g_dot_sum[sa] += g_dot;
                    local_g[sa]   += sub_dt_val * g_dot;
                    // local_g[sa] = std::max(local_g[sa] + sub_dt_val * g_dot, tau_0 * 0.001);
                    local_g[sa]    = std::max(local_g[sa], tau_0 * 0.001);
                }

                for (int s = 0; s < N_Slips_h; s++)
                    local_gamma[s] += sub_dt_val * gdot_step[s];

            } // end sub-stepping

            // --- Write outputs ---
            {
                int ci = 0;
                for (int i = 0; i < 3; i++)
                    for (int j = 0; j < 3; j++)
                    {
                        sigma_out[ci] = Fsigma(i,j);
                        Dp_out[ci]    = D_P_sum(i,j) / n_substeps;
                        Wp_out[ci]    = Omega_P_sum(i,j) / n_substeps;
                        ci++;
                    }
            }

            for (int s = 0; s < N_Slips_h; s++)
            {
                gamma_dot_out[s] = gamma_dot_sum[s] / n_substeps;
                gdot_out[s]      = g_dot_sum[s]     / n_substeps;
            }
        }; // end constitutiveKernel


        // =============================================================
        //  Rotation update kernel (shared by cells & boundary)
        // =============================================================
        auto rotationKernel = [&](
            const double* omega_in,
            const double* Wp_in,
            const double* rot_in,
            double* rot_out)
        {
            Fastor::Tensor<double,3,3> W_e(0.0);
            {
                int ci = 0;
                for (int ii = 0; ii < 3; ii++)
                    for (int jj = 0; jj < 3; jj++)
                    {
                        W_e(ii,jj) = omega_in[ci] - Wp_in[ci];
                        ci++;
                    }
            }

            auto W_e_sq = Fastor::matmul(W_e, W_e);

            double theta = Foam::sqrt(
                  Foam::pow(W_e(2,1), 2.0)
                + Foam::pow(W_e(0,2), 2.0)
                + Foam::pow(W_e(1,0), 2.0));

            Fastor::Tensor<double,3,3> expW(0.0);
            expW(0,0) = 1.0; expW(1,1) = 1.0; expW(2,2) = 1.0;

            if (theta > SMALL)
            {
                double tdt = theta * dt_val;
                expW += (Foam::sin(tdt) / theta) * W_e
                      + ((1.0 - Foam::cos(tdt)) / (theta*theta)) * W_e_sq;
            }

            const double rp1 = rot_in[0], rp2 = rot_in[1], rp3 = rot_in[2];
            const double rc1 = Foam::cos(rp1), rs1 = Foam::sin(rp1);
            const double rc2 = Foam::cos(rp2), rs2 = Foam::sin(rp2);
            const double rc3 = Foam::cos(rp3), rs3 = Foam::sin(rp3);

            Fastor::Tensor<double,3,3> Rold(0.0);
            Rold(0,0) =  rc2*rc3;
            Rold(0,1) =  rs1*rs2*rc3 - rc1*rs3;
            Rold(0,2) =  rc1*rs2*rc3 + rs1*rs3;
            Rold(1,0) =  rc2*rs3;
            Rold(1,1) =  rs1*rs2*rs3 + rc1*rc3;
            Rold(1,2) =  rc1*rs2*rs3 - rs1*rc3;
            Rold(2,0) = -rs2;
            Rold(2,1) =  rs1*rc2;
            Rold(2,2) =  rc1*rc2;

            auto Rnew = Fastor::matmul(expW, Rold);

            double psi_a = 0.0, theta_a = 0.0, phi_a = 0.0;

            if (mag(Rnew(2,0)) < 1.0 - 1e-12)
            {
                theta_a = -std::asin(Rnew(2,0));
                double ct = std::cos(theta_a);
                psi_a = std::atan2(Rnew(2,1)/ct, Rnew(2,2)/ct);
                phi_a = std::atan2(Rnew(1,0)/ct, Rnew(0,0)/ct);
            }
            else
            {
                phi_a = 0.0;
                if (Rnew(2,0) < 0.0)
                {
                    theta_a =  M_PI / 2.0;
                    psi_a   =  std::atan2(Rnew(0,1), Rnew(0,2));
                }
                else
                {
                    theta_a = -M_PI / 2.0;
                    psi_a   =  std::atan2(-Rnew(0,1), -Rnew(0,2));
                }
            }

            rot_out[0] = psi_a;
            rot_out[1] = theta_a;
            rot_out[2] = phi_a;
        }; // end rotationKernel


        // ================================================================
        //  CELL INTERIOR  –  constitutive integration
        // ================================================================
        forAll(mesh.C(), cell)
        {
            double sig_in[9], eps_in[9], omg_in[9], rot_ang[3];
            double g_in[MAX_SLIPS], gam_in[MAX_SLIPS];
            for (int ci = 0; ci < 9; ci++)
            {
                sig_in[ci] = sigtemp[h][cell].component(ci);
                eps_in[ci] = epsilon_dot[cell].component(ci);
                omg_in[ci] = omega_dot[cell].component(ci);
            }
            rot_ang[0] = Rot1[h][cell];
            rot_ang[1] = Rot2[h][cell];
            rot_ang[2] = Rot3[h][cell];
            for (int s = 0; s < N_Slips_h; s++)
            {
                g_in[s]   = g_alphas[h][s][cell];
                gam_in[s] = gamma_alphas[h][s][cell];
            }

            double sig_out[9], Dp_out[9], Wp_out[9];
            double gdot_out[MAX_SLIPS], gamma_dot_out[MAX_SLIPS];

            // Only compute on cells with phase fraction of n.N less than 0.95 
            constitutiveKernel(
                sig_in, eps_in, omg_in, rot_ang,
                g_in, gam_in, (PopBal_av[h][cell] > SMALL),
                sig_out, Dp_out, Wp_out, gdot_out, gamma_dot_out,&C_rotated[h][cell]);

            for (int ci = 0; ci < 9; ci++)
            {
                sigtemp[h][cell].component(ci)        = sig_out[ci];
                Plastic_Strain[h][cell].component(ci) = Dp_out[ci];
                Omega_P[h][cell].component(ci)        = Wp_out[ci];
            }
            for (int s = 0; s < N_Slips_h; s++)
            {
                gamma_alpha_dots[h][s][cell] = gamma_dot_out[s];
                g_alpha_dots[h][s][cell]     = gdot_out[s];
            }
        }


        // ================================================================
        //  BOUNDARY FACES  –  constitutive integration
        // ================================================================
        forAll(T.boundaryField(), patchI)
        {
            if (!T.boundaryField()[patchI].coupled()
             && T.boundaryField()[patchI].size() > 0)
            {
                forAll(T.boundaryField()[patchI], faceI)
                {
                    double sig_in[9], eps_in[9], omg_in[9], rot_ang[3];
                    double g_in[MAX_SLIPS], gam_in[MAX_SLIPS];
                    for (int ci = 0; ci < 9; ci++)
                    {
                        sig_in[ci] = sigtemp[h].boundaryField()[patchI][faceI].component(ci);
                        eps_in[ci] = epsilon_dot.boundaryField()[patchI][faceI].component(ci);
                        omg_in[ci] = omega_dot.boundaryField()[patchI][faceI].component(ci);
                    }
                    rot_ang[0] = Rot1[h].boundaryField()[patchI][faceI];
                    rot_ang[1] = Rot2[h].boundaryField()[patchI][faceI];
                    rot_ang[2] = Rot3[h].boundaryField()[patchI][faceI];
                    for (int s = 0; s < N_Slips_h; s++)
                    {
                        g_in[s]   = g_alphas[h][s].boundaryField()[patchI][faceI];
                        gam_in[s] = gamma_alphas[h][s].boundaryField()[patchI][faceI];
                    }

                    double sig_out[9], Dp_out[9], Wp_out[9];
                    double gdot_out[MAX_SLIPS], gamma_dot_out[MAX_SLIPS];

                    constitutiveKernel(
                        sig_in, eps_in, omg_in, rot_ang,
                        g_in, gam_in,
                        (PopBal_av[h].boundaryField()[patchI][faceI] > SMALL),
                        sig_out, Dp_out, Wp_out, gdot_out, gamma_dot_out,nullptr);//nullptr as only need cell values for mechanical driving force

                    for (int ci = 0; ci < 9; ci++)
                    {
                        sigtemp[h].boundaryFieldRef()[patchI][faceI].component(ci)        = sig_out[ci];
                        Plastic_Strain[h].boundaryFieldRef()[patchI][faceI].component(ci) = Dp_out[ci];
                        Omega_P[h].boundaryFieldRef()[patchI][faceI].component(ci)        = Wp_out[ci];
                    }
                    for (int s = 0; s < N_Slips_h; s++)
                    {
                        gamma_alpha_dots[h][s].boundaryFieldRef()[patchI][faceI] = gamma_dot_out[s];
                        g_alpha_dots[h][s].boundaryFieldRef()[patchI][faceI]     = gdot_out[s];
                    }
                }
            }
        }


        // Orientation update: apply the exponential map ONCE per time
        // step, on the final outer corrector 
        if (pimple.finalIter())
        {

        // ================================================================
        //  ROTATION UPDATE  –  cell interior
        // ================================================================
        forAll(mesh.C(), cell)
        {
            double omg_in[9], Wp_in[9], rot_in[3], rot_out[3];
            for (int ci = 0; ci < 9; ci++)
            {
                omg_in[ci] = omega_dot[cell].component(ci);
                Wp_in[ci]  = Omega_P[h][cell].component(ci);
            }
            rot_in[0] = Rot1[h][cell];
            rot_in[1] = Rot2[h][cell];
            rot_in[2] = Rot3[h][cell];

            rotationKernel(omg_in, Wp_in, rot_in, rot_out);

            Rot1[h][cell] = rot_out[0];
            Rot2[h][cell] = rot_out[1];
            Rot3[h][cell] = rot_out[2];
        }


        // ================================================================
        //  ROTATION UPDATE  –  boundary faces
        // ================================================================
        forAll(T.boundaryField(), patchI)
        {
            if (!T.boundaryField()[patchI].coupled()
             && T.boundaryField()[patchI].size() > 0)
            {
                forAll(T.boundaryField()[patchI], faceI)
                {
                    double omg_in[9], Wp_in[9], rot_in[3], rot_out[3];
                    for (int ci = 0; ci < 9; ci++)
                    {
                        omg_in[ci] = omega_dot.boundaryField()[patchI][faceI].component(ci);
                        Wp_in[ci]  = Omega_P[h].boundaryField()[patchI][faceI].component(ci);
                    }
                    rot_in[0] = Rot1[h].boundaryField()[patchI][faceI];
                    rot_in[1] = Rot2[h].boundaryField()[patchI][faceI];
                    rot_in[2] = Rot3[h].boundaryField()[patchI][faceI];

                    rotationKernel(omg_in, Wp_in, rot_in, rot_out);

                    Rot1[h].boundaryFieldRef()[patchI][faceI] = rot_out[0];
                    Rot2[h].boundaryFieldRef()[patchI][faceI] = rot_out[1];
                    Rot3[h].boundaryFieldRef()[patchI][faceI] = rot_out[2];
                }
            }
        }

        Rot1[h].correctBoundaryConditions();
        Rot2[h].correctBoundaryConditions();
        Rot3[h].correctBoundaryConditions();

        } // end if (pimple.finalIter()): orientation updated once per step

        // Accumulate global stress
        sigmaGlobaltemp += PopBal[h] * sigtemp[h];

    } // end loop over phases

    // EOS mode: sigtemp/sigma hold per-phase deviatoric part.  Add the
    // volumetric part later ....Mie-Gruneisen pressure plus
    // artificial shock viscosity
    if (useEOS)
    {
        sigmaGlobaltemp -= (p + qAV)*I_EOS;
    }

    sigmaGlobaltemp.correctBoundaryConditions();

    // ==============================================================
    //  Commit state on final PIMPLE iteration
    // ==============================================================
    if (pimple.finalIter())
    {
        Info << "Updating committed fields on final PIMPLE iteration"
             << endl;

        for (label h = 0; h < N_Ori; h++)
        {
            sigma[h] = sigtemp[h];
            sigma[h].correctBoundaryConditions();

            Rot1[h].correctBoundaryConditions();
            Rot2[h].correctBoundaryConditions();
            Rot3[h].correctBoundaryConditions();
        }

        sigmaGlobal = sigmaGlobaltemp;
        sigmaGlobal.correctBoundaryConditions();

        D = D + runTime.deltaT() * U;
        D.correctBoundaryConditions();

        // diagnostic only: in EOS mode p is the 
        // pressure (and tr(sigmaGlobal)/3 = -(p + q)
        if (!useEOS)
        {
            p = -(1.0/3.0) * tr(sigmaGlobal);
        }
    }

    } // end else (stressUpdate)
}
