{


    


    Sum_Phases.primitiveFieldRef() = 0.0;  
GBREF.primitiveFieldRef() = 0.0;       

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


            if(fracture)
            {
            //
            // get stress in fastor format
            // --- Extract strain rate & spin ---
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

            //need inverse of stiffness
            Fastor::Tensor<double,3,3,3,3> C_rot_inverse = Material::invertRankFour(C_rotated[active_phases[count_i]][cell]);
            // C_rot_inverse = Material::invertRankFour(C_rotated[active_phases[count_i]][cell]);
            //
            dfdp_mech [count_i] = -0.5*Material::doubleContractionTwice(Fsigma_i,C_rot_inverse);




        }//enf of fracture conditional



        
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

                double driving_force_cut_off = 0.0;

                // if(PopBal[active_phases[count_i]][cell]>0.1&&PopBal[active_phases[count_j]][cell]>0.1){
                    driving_force_cut_off= (dfdp_G[count_i] - dfdp_G[count_j]) * ((8.0) / pi.value()) * Foam::sqrt(PopBal[active_phases[count_i]][cell] * PopBal[active_phases[count_j]][cell]);//ORIGINAL

                double rhsval = (int_mobility) * ((dfdp[count_i] - dfdp[count_j])+(dfdp_aniso[count_i]/*-dfdp_aniso[count_j]*/)+(dfdp_mech[count_i] -dfdp_mech[count_j]) -
                                                  driving_force_cut_off
                                                 );
                //   Info<<"RHS"<<rhsval<<endl;
                // Driving_Force_info[cell]+=driving_force_cut_off;

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
