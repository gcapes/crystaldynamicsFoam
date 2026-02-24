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
            // double mu_[N_Comps]={0.0};
            // for(label l = 0; l < N_Comps; l++){
            //     // PhaseComp_[l]=PhaseComposition[my_phase][l][cell];
            //     mu_[l]=Compbal[l][cell];
            // }

    ///////////////////////////////////////////////////////////Code to show how to output an example field read from material dictionaries
    // A_out[cell] *= 0.0;
    // for (label count = 0; count < N_Ori; count++) {
    // // Phase2Material[count]
    // A_out[cell]+=PopBal[count][cell]*test_A[Phase2Material[count]];

    // }
    ///////////////////////////////////////////////////////////Code to show how to output an example field read from material dictionaries



    if (n_active_phases > 1)   //check if neighbour //add on neighbour active phases to active_phases list
    {
        // Info<<"HERE"<<endl;
        double dfdp[n_active_phases] =
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
            //  Info<<"HERE ... "<<active_phases[count_i]<<endl;
            for (label count_j = count_i + 1; count_j < n_active_phases; count_j++)
            {

            //   double tens_i_[9]={0.0};
            //   double tens_j_[9]={0.0};
            // for(label l = 0; l < 9; l++){
            //     tens_i_[l]=Roti[active_phases[count_i]][cell].component(l);
            //     tens_j_[l]=-tens_i_[l];//Roti[active_phases[count_j]][cell].component(l);//-tens_i_[l];//
            // }   

            // double vec_i_[3]={0.0};
            // double vec_j_[3]={0.0};
            // for(label m = 0; m < 3; m++){
            //     vec_i_[m]=gradi[active_phases[count_i]][cell].component(m);
            //     vec_j_[m]=-vec_i_[m];//gradi[active_phases[count_j]][cell].component(m);//-vec_i_[m];//
            // } 


                // Info<<"HERE2"<<endl;
                double int_energy = PF_Energy0_array[Phase2Material[active_phases[count_i]]][Phase2Material[active_phases[count_j]]];
                // *(1.0+(0.05*((Foam::pow(vec_i_[0],4)+Foam::pow(vec_i_[1],4)+Foam::pow(vec_i_[2],4))-3.0)));

                double anisostrength=PF_Aniso_array[Phase2Material[active_phases[count_i]]][Phase2Material[active_phases[count_j]]];

                dfdp[count_i] -= 4.0 * int_energy * PopBal[active_phases[count_j]][cell] / (interface_width.value()) +
                                 4.0 * int_energy * Laplacian_n[active_phases[count_j]][cell] * (interface_width.value()) / (pi.value() * pi.value());
                dfdp[count_j] -= 4.0 * int_energy * PopBal[active_phases[count_i]][cell] / (interface_width.value()) +
                                 4.0 * int_energy * Laplacian_n[active_phases[count_i]][cell] * (interface_width.value()) / (pi.value() * pi.value());

                // if (n_active_phases < 3) {
                //     double temp_cap_aniso = /*BoundarySwitch[cell]*/anisostrength*int_energy*dummy.aniso_energy(tens_i_,tens_j_,vec_i_,vec_j_,PopBal[active_phases[count_i]][cell],PopBal[active_phases[count_j]][cell],interface_width.value(),/*deltaN3.value()*/1e-8/Foam::pow(mesh.V()[cell], 1.0/3.0));
                //             temp_cap_aniso=min(temp_cap_aniso,1e6);
                //             temp_cap_aniso=max(temp_cap_aniso,-1e6);
                // dfdp_aniso[count_i]-=/*2.0*Foam::sqrt(PopBal[active_phases[count_i]][cell]*PopBal[active_phases[count_j]][cell])*/temp_cap_aniso;///*BoundarySwitch[cell]*/anisostrength*int_energy*dummy.aniso_energy(tens_i_,tens_j_,vec_i_,vec_j_,PopBal[active_phases[count_i]][cell],PopBal[active_phases[count_j]][cell],interface_width.value(),/*deltaN3.value()*/1e-8/Foam::pow(mesh.V()[cell], 1.0/3.0));
                // // dfdp_aniso[count_j]-=anisostrength*int_energy*dummy.aniso_energy(tens_j_,tens_i_,vec_j_,vec_i_,PopBal[active_phases[count_j]][cell],PopBal[active_phases[count_i]][cell],interface_width.value(),deltaN3.value());
                // A_out[cell]+=    /*BoundarySwitch[cell]*/(dfdp_aniso[count_i]);
                // }
                // Info<<count_i<<": \t"<<dfdp_aniso[count_i]<<endl; 

                        

                for (label count_k = count_j + 1; count_k < n_active_phases; count_k++)
                {
                    double T_junction_energy = 0.0;
                    double int_energy_ij = PF_Energy0_array[Phase2Material[active_phases[count_i]]][Phase2Material[active_phases[count_j]]];
                    double int_energy_ik = PF_Energy0_array[Phase2Material[active_phases[count_i]]][Phase2Material[active_phases[count_k]]];
                    double int_energy_jk = PF_Energy0_array[Phase2Material[active_phases[count_j]]][Phase2Material[active_phases[count_k]]];

                    // T_junction_energy = T_junction_energy > int_energy_ij ? T_junction_energy : int_energy_ij;//condition ? result_if_true : result_if_false
                    // T_junction_energy = T_junction_energy > int_energy_ik ? T_junction_energy : int_energy_ik;
                    // T_junction_energy = T_junction_energy > int_energy_jk ? T_junction_energy : int_energy_jk;
                    T_junction_energy = max(int_energy_ij, max(int_energy_ik, int_energy_jk));

                    dfdp[count_i] -= 4.0 * T_junction_energy * (PopBal[active_phases[count_j]][cell] * PopBal[active_phases[count_k]][cell]) / interface_width.value();
                    dfdp[count_j] -= 4.0 * T_junction_energy * (PopBal[active_phases[count_i]][cell] * PopBal[active_phases[count_k]][cell]) / interface_width.value();
                    dfdp[count_k] -= 4.0 * T_junction_energy * (PopBal[active_phases[count_i]][cell] * PopBal[active_phases[count_j]][cell]) / interface_width.value();

                }



            }


            // A_out[cell]+=PopBal[active_phases[count_i]][cell]*test_A[Phase2Material[active_phases[count_i]]];//testing how to access phase properties based on mapping material array

            // if(active_phases[count_i]==N_Ori-1){
            //     dfdp_G[count_i] =1e4;;;//Fusion_PF_driving_force[cell];
            //     }


        }

        // for (label i = 0; i < n_active_phases; i++)
        // {

        //     // PhaseComposition[count_phase][count_i][cell];


        //     label my_phase = active_phases[i];
        //     label my_mat = Phase2Material[my_phase];

        //     // double PhaseComp_[N_Comps]={0.0};
        //     // // double mu_[N_Comps]={0.0};
        //     // for(label l = 0; l < N_Comps; l++){
        //     //         // loop over all lattices ..... populate phaseComp[l][k] to pass to function
        //     //     PhaseComp_[l]=PhaseComposition[my_phase][l][cell];
                
        //     // }
            


        //     //  // Info<<MatList[my_mat].deltaG(mu_,PhaseComp_,Temperature[cell],N_Comps)<<endl;;
        //     // if(MatList[my_mat].Mat_type=="Quadratic"){
        //     // dfdp_G[i] =(1.0/V_mol.value())*MatList[my_mat].deltaG(mu_,PhaseComp_,Temperature[cell],N_Comps);
        //     // }
        //     // if(MatList[my_mat].Mat_type=="TD_Model"){
        //     //     // Info<<"cell: "<<cell<<endl;
        //     // dfdp_G[i] = (1.0/V_mol.value())*MatList[my_mat].deltaG_TD(mu_,PhaseComp_,Temperature[cell],N_Comps);
        //     // }
        //     //  if(MatList[my_mat].Mat_type=="SL_Model"){
        //     //     label N_total = N_Comps + 1;

        //     //         // Get stored site fractions
        //     //     double sitefrac[2*N_total];
        //     //     for(label c = 0; c < N_total; c++) {
        //     //         sitefrac[c] = PhaseComposition2SL[my_phase][c][0][cell];
        //     //         sitefrac[N_total + c] = PhaseComposition2SL[my_phase][c][1][cell];
        //     //         }
    
        //     //     // Calculate driving force
        //     //     dfdp_G[i] = MatList[my_mat].delta_G_2SL(mu_, sitefrac, Temperature[cell], N_Comps);

        //     // }

        // }
        // A_out[cell] = dfdp_G[0] - dfdp_G[1];
        Driving_Force_info[cell]*=0.0;// dfdp_G[0] - dfdp_G[1];

        // Info<<dfdp[0]<<endl;
        for (label count_i = 0; count_i < n_active_phases; count_i++)
        {
            //    Info<<"HERE"<<endl;
            for (label count_j = count_i + 1; count_j < n_active_phases; count_j++)
            {

                double int_mobility = PF_Mobility0_array[Phase2Material[active_phases[count_i]]][Phase2Material[active_phases[count_j]]];

                double driving_force_cut_off = 0.0;

                // if(PopBal[active_phases[count_i]][cell]>0.1&&PopBal[active_phases[count_j]][cell]>0.1){
                    driving_force_cut_off= (dfdp_G[count_i] - dfdp_G[count_j]) * ((8.0) / pi.value()) * Foam::sqrt(PopBal[active_phases[count_i]][cell] * PopBal[active_phases[count_j]][cell]);//ORIGINAL

                double rhsval = (int_mobility) * ((dfdp[count_i] - dfdp[count_j])+(dfdp_aniso[count_i]/*-dfdp_aniso[count_j]*/) -
                                                  driving_force_cut_off
                                                 );
                //   Info<<"RHS"<<rhsval<<endl;
                Driving_Force_info[cell]+=driving_force_cut_off;

                dpdt[active_phases[count_i]][cell] += rhsval;
                dpdt[active_phases[count_j]][cell] -= rhsval;
            }
            dpdt[active_phases[count_i]][cell] /= n_active_phases;

        }


        //
        // for (label count_i = 0; count_i < n_active_phases; count_i++)
        // {
        //     PopBal[active_phases[count_i]][cell]+=runTime.deltaT().value()*dpdt[active_phases[count_i]][cell];
        // }

    }


}























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

    //  Project any out of bounds phases back onto the Gibbs Simplex
    for (label count = 0; count < N_Ori; count++)
    {

        volScalarField & ni = PopBal[count];

        ni = max(min(ni, 1.0), 0.0);
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
    }


}
