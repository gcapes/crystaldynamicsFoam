{



    // Mechanical driving-force regularisation ..... both OFF by default


    // PF_drivingForceMax  =  cap on the
    //   mechanical driving force. 
    //
    // PF_mechPhaseMin =  trace-phase nucleation guard.  The
    //   mech term sqrt(phi_i*phi_j) amplifies spurious  third phases
    //   at a driven boundary......requiring BOTH cell values above this floor stops spurious
    //   decoration of moving boundaries while leaving genuinely
    //   present neighbours untouched.

    

    const scalar PF_drivingForceMax
    (
        GrainProperties.lookupOrDefault<scalar>("PF_drivingForceMax", -1.0)
    );
    const scalar PF_mechPhaseMin
    (
        GrainProperties.lookupOrDefault<scalar>("PF_mechPhaseMin", 0.0)
    );

    
        const Switch PF_amorSplit
    (
        GrainProperties.lookupOrDefault<Switch>("PF_amorSplit", false)
    );




    Sum_Phases.primitiveFieldRef() = 0.0;
GBREF.primitiveFieldRef() = 0.0;
W_mech.primitiveFieldRef()      = 0.0;
W_interface.primitiveFieldRef() = 0.0;
dpdt_crack.primitiveFieldRef()  = 0.0;

for (label count = 0; count < N_Ori; count++)
{
    volScalarField& ni = PopBal[count];

    Laplacian_n[count] = fvc::laplacian(ni);
    Laplacian_n[count].correctBoundaryConditions();

    PopBal_av[count] = fvc::average(ni);
    PopBal_av[count].correctBoundaryConditions();

    dpdt[count].primitiveFieldRef() = 0.0;  
    dpdt[count].correctBoundaryConditions();

    GBREF += PopBal_av[count] * PopBal_av[count];
}

GBREF.correctBoundaryConditions();







if (PF_drivingForceActive)
{
forAll(mesh.C(), cell)
{
    label active_phases[N_Ori] =
    {
        0
    };
    label n_active_phases = 0;



    for (label count = 0; count < N_Ori; count++)
    {
        dpdt[count][cell] = 0.0;
        if (PopBal_av[count][cell] > 0.0)
        {
            active_phases[n_active_phases++] = count;
        }

    }



    if (n_active_phases > 1)   //check if neighbour //add on neighbour active phases to active_phases list
    {
        // Info<<"HERE"<<endl;
        double dfdp[n_active_phases] =
        {
            0.0
        };
        double dfdp_mech[n_active_phases] =
        {
            0.0
        };
        double dfdp_aniso[n_active_phases] =
        {
            0.0
        };
        double dfdp_G[n_active_phases] =
        {
            0.0
        };
        for (label count_i = 0; count_i < n_active_phases; count_i++)
        {




            if (PopBal_av[active_phases[count_i]][cell]
                > (stressUpdate ? SMALL : SMALL))// should be small whatever but this is nasty, clean up
            {
            //
            // get stress in fastor format
            Fastor::Tensor<double,3,3> Fsigma_i(0.0);
            {
                int index_S = 0;
                for (int i = 0; i < 3; i++){
                    for (int j = 0; j < 3; j++)
                    {
                        Fsigma_i(i,j)   = sigma[active_phases[count_i]][cell].component(index_S);
                        
                        index_S++;
                    }

                }
            }

            // ---------------------------------------------------------
            // Zero-stiffness phases store no elastic
            // energy -    so inverting C_rotated would produce warnings.  The
            // crack grows by releasing the intact phase's energy; its
            // own dfdp_mech is exactly zero. . . . . .This also removes the invertrank4 warning
            // ---------------------------------------------------------
            const double CdiagPF = // CHECK IF THERE ARE ANY STIFFNESS , these compinents should always exist (at least one right?!?!)
                C_rotated[active_phases[count_i]][cell](0,0,0,0)
              + C_rotated[active_phases[count_i]][cell](1,1,1,1)
              + C_rotated[active_phases[count_i]][cell](2,2,2,2);

            if (CdiagPF > SMALL)
            {
            //need inverse of stiffness
            Fastor::Tensor<double,3,3,3,3> C_rot_inverse = Material::invertRankFour(C_rotated[active_phases[count_i]][cell]);
            // C_rot_inverse = Material::invertRankFour(C_rotated[active_phases[count_i]][cell]);
            //
            // !!! Negative sign: dfdp_mech = -(dF/dn)_mech (accumulated with -= not +=) !!!
            //
            // DEVIATORIC stored energy.  Under MieGruneisen the
            // per-phase sigma fields hold DEVIATORS because we use EOS for the hydrostatic buts there.  
            // For cubic C the compliance maps deviators to deviators and this contraction
            // is EXACT the deviatoric energy
            

            //TO DO THIS FULLY PROPERLY WE'LL NEED a finite strain formulation so its the deformation gradient were dealing with rether than stress - for future work - - 
            // for now this is only an issue with very anisotropic materials like Zn, i think this is the same approach other codes take, but we'll do better in future to get the pressure - shear coupling through F^e



           
            if (!useEOS && PF_amorSplit)
            {
                // Deviatoric energy always
                // drives; volumetric drives in tension only,
                double sigm = (Fsigma_i(0,0) + Fsigma_i(1,1)
                             + Fsigma_i(2,2))/3.0;
                Fastor::Tensor<double,3,3> s_dev = Fsigma_i;
                for (int d3 = 0; d3 < 3; d3++)
                {
                    s_dev(d3,d3) -= sigm;
                }

                dfdp_mech[count_i] =
                    -0.5*Material::doubleContractionTwice
                    (s_dev, C_rot_inverse);

                if (sigm > 0.0)   // tensile mean stress only
                {
                    double ISI = 0.0;
                    for (int a3 = 0; a3 < 3; a3++)
                    {
                        for (int b3 = 0; b3 < 3; b3++)
                        {
                            ISI += C_rot_inverse(a3,a3,b3,b3);
                        }
                    }
                    dfdp_mech[count_i] -= 0.5*sigm*sigm*ISI;
                }
            }
            else
            {
                //linearelastic but with compression driving too

            dfdp_mech [count_i] = -0.5*Material::doubleContractionTwice(Fsigma_i,C_rot_inverse);//SIGN IS RIGHT TOM YOU IDIOT.... SHOULD BE -ve, credit Harry Smith, anti-credit Tom Flint for the wrong sign
            }

            if (useEOS)
            {
                const Material& matPF =
                    MatList[Phase2Material[active_phases[count_i]]];
                const double A0PF = matPF.rho_0*matPF.c0*matPF.c0;
                const double pTens = max(-p[cell], 0.0);//ONLY TENSION
                if (A0PF > SMALL)
                {
                    dfdp_mech[count_i] -= 0.5*pTens*pTens/A0PF;
                }
            }

            } // end if (CdiagPF > SMALL): zero-stiffness phases keep
              // dfdp_mech = 0.0 from initialisation




        }//enf of PF_drivingForceActive conditional



        
            //  Info<<"HERE ... "<<active_phases[count_i]<<endl;
            for (label count_j = count_i + 1; count_j < n_active_phases; count_j++)
            {



                // Info<<"HERE2"<<endl;
                double int_energy = PF_Energy0_array[Phase2Material[active_phases[count_i]]][Phase2Material[active_phases[count_j]]];

                /*
                
                dfdp[count_i] accumulates terms involving phase j's field and Laplacian, and dfdp[count_j] accumulates terms involving phase i's.
                */

                dfdp[count_i] -= 4.0 * int_energy * PopBal[active_phases[count_j]][cell] / (interface_width.value()) +
                                 4.0 * int_energy * Laplacian_n[active_phases[count_j]][cell] * (interface_width.value()) / (pi.value() * pi.value());
                dfdp[count_j] -= 4.0 * int_energy * PopBal[active_phases[count_i]][cell] / (interface_width.value()) +
                                 4.0 * int_energy * Laplacian_n[active_phases[count_i]][cell] * (interface_width.value()) / (pi.value() * pi.value());

                
                
                                //  doubleContractionTwice(sigma_i. C_i);

                // dfdp_mech [count_i] -=0.0;

                // dfdp_mech [count_j] -=0.0;
    

                for (label count_k = count_j + 1; count_k < n_active_phases; count_k++)
                {
                    // double T_junction_energy = 0.0;
                    double int_energy_ij = PF_Energy0_array[Phase2Material[active_phases[count_i]]][Phase2Material[active_phases[count_j]]];
                    double int_energy_ik = PF_Energy0_array[Phase2Material[active_phases[count_i]]][Phase2Material[active_phases[count_k]]];
                    double int_energy_jk = PF_Energy0_array[Phase2Material[active_phases[count_j]]][Phase2Material[active_phases[count_k]]];

                    // T_junction_energy = T_junction_energy > int_energy_ij ? T_junction_energy : int_energy_ij;//condition ? result_if_true : result_if_false
                    // T_junction_energy = T_junction_energy > int_energy_ik ? T_junction_energy : int_energy_ik;
                    // T_junction_energy = T_junction_energy > int_energy_jk ? T_junction_energy : int_energy_jk;
                    // T_junction_energy = max(int_energy_ij, max(int_energy_ik, int_energy_jk));

                    // dfdp[count_i] -= 4.0 * T_junction_energy * (PopBal[active_phases[count_j]][cell] * PopBal[active_phases[count_k]][cell]) / interface_width.value();
                    // dfdp[count_j] -= 4.0 * T_junction_energy * (PopBal[active_phases[count_i]][cell] * PopBal[active_phases[count_k]][cell]) / interface_width.value();
                    // dfdp[count_k] -= 4.0 * T_junction_energy * (PopBal[active_phases[count_i]][cell] * PopBal[active_phases[count_j]][cell]) / interface_width.value();



                    dfdp[count_i] -= 4.0 * (int_energy_ij * int_energy_ik)
                                    / max(int_energy_ij + int_energy_ik, SMALL)
                                    * (PopBal[active_phases[count_j]][cell]
                                    * PopBal[active_phases[count_k]][cell])
                                    / interface_width.value();

                               
                    dfdp[count_j] -= 4.0 * (int_energy_ij * int_energy_jk)
                                    / max(int_energy_ij + int_energy_jk, SMALL)
                                    * (PopBal[active_phases[count_i]][cell]
                                    * PopBal[active_phases[count_k]][cell])
                                     / interface_width.value();


                    dfdp[count_k] -= 4.0 * (int_energy_ik * int_energy_jk)
                                    / max(int_energy_ik + int_energy_jk, SMALL)
                                    * (PopBal[active_phases[count_i]][cell]
                                    * PopBal[active_phases[count_j]][cell])
                                    / interface_width.value();



                }



            }


        }

        
        for (label count_i = 0; count_i < n_active_phases; count_i++)
        {
            //    Info<<"HERE"<<endl;
            for (label count_j = count_i + 1; count_j < n_active_phases; count_j++)
            {

                double int_mobility = PF_Mobility0_array[Phase2Material[active_phases[count_i]]][Phase2Material[active_phases[count_j]]];

                double driving_force_cut_off = 0.0;// PLACEHOLDER for chemistry from CIPHERFOMA

                // if(PopBal[active_phases[count_i]][cell]>0.1&&PopBal[active_phases[count_j]][cell]>0.1){
                    driving_force_cut_off= (dfdp_G[count_i] - dfdp_G[count_j]) * ((8.0) / pi.value()) * Foam::sqrt(PopBal[active_phases[count_i]][cell] * PopBal[active_phases[count_j]][cell]);//ORIGINAL

                
                    
                    // Pairwise mechanical driving force, regularised:
                    // (i) trace-phase guard, (ii) magnitude limiter.
                    double dMech =
                        dfdp_mech[count_i] - dfdp_mech[count_j];

                    if// NOT SURE IF NEEDED BECAUSE THE FIX WEAS DECREASING THE TIMESTEP FOR THE BLIPPING ON THE INTERFACE
                    (
                        PopBal[active_phases[count_i]][cell]
                            <= PF_mechPhaseMin
                     || PopBal[active_phases[count_j]][cell]
                            <= PF_mechPhaseMin
                    )
                    {
                        dMech = 0.0;
                    }

                    if (PF_drivingForceMax > 0.0)
                    {
                        dMech = std::max(std::min(dMech,
                                    PF_drivingForceMax),
                                   -PF_drivingForceMax);
                    }

                    double rhsval = (int_mobility) * 
                    (
                        (dfdp[count_i] - dfdp[count_j])
  + (dfdp_aniso[count_i] - dfdp_aniso[count_j])   // antisymmetric, like every other pair term (currently zero)
  + ((8.0)/pi.value())
      * Foam::sqrt(max(
            PopBal[active_phases[count_i]][cell]
          * PopBal[active_phases[count_j]][cell], 0.0))
      * dMech
  - driving_force_cut_off
);


                //   Info<<"RHS"<<rhsval<<endl;
                // Driving_Force_info[cell]+=driving_force_cut_off;

                dpdt[active_phases[count_i]][cell] += rhsval;
                dpdt[active_phases[count_j]][cell] -= rhsval;
            }
            dpdt[active_phases[count_i]][cell] /= n_active_phases;

            // Accumulate  diagnostics 
            // W_mech: phase-weighted elastic energy density 
            double ni_val = PopBal[active_phases[count_i]][cell];
            W_mech[cell]      += ni_val * (-dfdp_mech[count_i]);
            // W_interface: phase-weighted interfacial energy.
            // dfdp is negative for phases that want to shrink
            W_interface[cell] += ni_val * (-dfdp[count_i]);
        }

        // dpdt_crack: ROC of phase fields
        for (label count_i = 0; count_i < n_active_phases; count_i++)
        {
            if (Phase2Material[active_phases[count_i]] != 0)
            {
                dpdt_crack[cell] += dpdt[active_phases[count_i]][cell];
            }
        }




    }  // end if (n_active_phases > 1)

}  // end forAll(mesh.C(), cell)
} // end if (PF_drivingForceActive)






















    Sum_Phases *= 0.0;
    for (label count = 0; count < N_Ori; count++)
    {

        volScalarField &ni = PopBal[count];
        // volScalarField& Su=Sum_parts[count];
        volScalarField & RHS = dpdt[count];

        // ni=ni+runTime.deltaT()*RHS;

        solve(
            fvm::ddt(rho, ni) + fvm::div(rhoPhi, ni)
            ==
            rho*RHS

        );
        ni.correctBoundaryConditions();
        Sum_Phases += ni;

        Info << "Max phase field(" << count << ") = " << gMax(ni) << endl;
        Info << "Min phase field(" << count << ") = " << gMin(ni) << endl;
    }
    Sum_Phases.correctBoundaryConditions();


    // Remember each cell's dominant phase BEFORE clipping: if a cell's
    // solve undershoots EVERY phase to <= 0, clipping leaves Sum = 0 and
    // the normalisation below would skip it -- a permanent vacuum cell
    // nothing can regrow (nonzero phi is required for every pairwise
    // term).  Restore such cells to their pre-clip dominant phase.
    labelList domPhasePF(mesh.nCells(), label(0));
    {
        scalarField domValPF(mesh.nCells(), -GREAT);
        for (label count = 0; count < N_Ori; count++)
        {
            const scalarField& niIF = PopBal[count].primitiveField();
            forAll(niIF, celli)
            {
                if (niIF[celli] > domValPF[celli])
                {
                    domValPF[celli] = niIF[celli];
                    domPhasePF[celli] = count;
                }
            }
        }
    }

    for (label count = 0; count < N_Ori; count++)
    {
        PopBal[count] = max(min(PopBal[count], 1.0), 0.0);
    }

    Sum_Phases *= 0.0;
    for (label count = 0; count < N_Ori; count++)
    {
        Sum_Phases += PopBal[count];
    }
    Sum_Phases.correctBoundaryConditions();

    {
        label nVacuumPF = 0;
        forAll(mesh.C(), celli)
        {
            if (Sum_Phases[celli] < SMALL)
            {
                PopBal[domPhasePF[celli]][celli] = 1.0;
                Sum_Phases[celli] = 1.0;
                nVacuumPF++;
            }
        }
        if (returnReduce(nVacuumPF, sumOp<label>()) > 0)
        {
            Info<< "PF projection: restored "
                << returnReduce(nVacuumPF, sumOp<label>())
                << " vacuum cells to their dominant phase "
                << "(all-phase undershoot -- check dt bounds)" << endl;
        }
    }

    for (label count = 0; count < N_Ori; count++)
    {

        volScalarField & ni = PopBal[count];

        volScalarField temp_holder = ni;

        forAll(mesh.C(), celli)
        {
            if (Sum_Phases[celli] > 0.0)
            {
                temp_holder[celli] *= (1.0 / Sum_Phases[celli]);
            }
        }

            forAll(Sum_Phases.boundaryField(), patchI)
{
    if
    (
        !Sum_Phases.boundaryField()[patchI].coupled()
     && Sum_Phases.boundaryField()[patchI].size() > 0
    )
    {
        forAll(Sum_Phases.boundaryField()[patchI], faceI)
        {

            if (Sum_Phases.boundaryField()[patchI][faceI] > 0.0)
        {
            temp_holder.boundaryFieldRef()[patchI][faceI] *= (1.0 / Sum_Phases.boundaryField()[patchI][faceI]);
        }

        }

    }

}


        ni = temp_holder;
        ni.correctBoundaryConditions();
    }

W_mech.correctBoundaryConditions();
W_interface.correctBoundaryConditions();
dpdt_crack.correctBoundaryConditions();

}