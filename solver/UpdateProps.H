{

// Info<<n_sub_steps<<endl;



// divU=fvc::div(U);

volTensorField TgradU(fvc::grad(U)().T());
TgradU.correctBoundaryConditions(); // already called by line above


omega_dot = 0.5*(fvc::grad(U)-TgradU); // this is the "symm" function, i.e. symm(gradU)
omega_dot.correctBoundaryConditions(); // already called by line above
epsilon_dot = 0.5*(fvc::grad(U)+TgradU) ;//- (1.0/3.0)*(fvc::div(phi))*Ident; // this looks like dev(gradU)
epsilon_dot.correctBoundaryConditions(); // already called by line above


dsigmadt = -(sigma & omega_dot) + (omega_dot & sigma);// -(tr(epsilon_dot)*sigma);
dsigmadt.correctBoundaryConditions();
// dsigmadt_crystal*=0.0;

// epsilon_P_dot*=0.0;
Omega_P*=0.0;
D_P*=0.0;

p*=0.0;

// Euler_Angles*=0.0;
// epsilon_dot_P_VM*=0.0;
sumgamma*=0.0;


forAll(mesh.C(), cell)
{


  Fastor::Tensor<double,3,3> Fsigma(0.0); //CAST OPENFOAM VARIABLES INTO FASTOR VARIABLES
  Fastor::Tensor<double,3,3> Fepsilon_dot(0.0);
  int count_index=0;
  for(int i = 0 ; i < 3 ; i++ ){
    for(int j = 0 ; j < 3 ; j++ ){
        // Info<<i<<","<<j<<endl;
      Fepsilon_dot(i,j)=epsilon_dot[cell].component(count_index);
      Fsigma(i,j) = sigma[cell].component(count_index);

      count_index++;
    }
  }









  // if(PopBal[h][cell]>SMALL){// if the phase is there, work out the stiffness contribution
//  Info<<"Phase: "<<h<<"\t "<<PopBal[h][cell]<<endl;

Fastor::Tensor<double,3,3> sigAv(0.0);
Fastor::Tensor<double,3,3> Omega_P_Average(0.0);
Fastor::Tensor<double,3,3> D_P_Average(0.0);


for (label h = 0; h < N_Ori; h++)
{

  Fastor::Tensor<double,3,3,3,3> Cprime(0.0); //global stiffness at this cell



double phi_1_total =Rot1[h][cell];//ROT_1[cell];// phi1[h] ;//+ omega_tot[cell].component(5);

double phi_2_total =Rot2[h][cell];//ROT_2[cell];// phi2[h] ;//+ omega_tot[cell].component(2);

double phi_3_total =Rot3[h][cell];//ROT_3[cell];// phi3[h] ;//+ omega_tot[cell].component(1);






Fastor::Tensor<double,3,3> Rot_matrix(0.0);

  Rot_matrix(0,0)=Foam::cos(phi_2_total)*Foam::cos(phi_3_total);
  Rot_matrix(0,1)=(Foam::sin(phi_1_total)*Foam::sin(phi_2_total)*Foam::cos(phi_3_total))-(Foam::cos(phi_1_total)*Foam::sin(phi_3_total));
  Rot_matrix(0,2)=(Foam::cos(phi_1_total)*Foam::sin(phi_2_total)*Foam::cos(phi_3_total))+(Foam::sin(phi_1_total)*Foam::sin(phi_3_total));

  Rot_matrix(1,0)=Foam::cos(phi_2_total)*Foam::sin(phi_3_total);
  Rot_matrix(1,1)=(Foam::sin(phi_1_total)*Foam::sin(phi_2_total)*Foam::sin(phi_3_total))+(Foam::cos(phi_1_total)*Foam::cos(phi_3_total));
  Rot_matrix(1,2)=(Foam::cos(phi_1_total)*Foam::sin(phi_2_total)*Foam::sin(phi_3_total))-(Foam::sin(phi_1_total)*Foam::cos(phi_3_total));

  Rot_matrix(2,0)=-Foam::sin(phi_2_total);
  Rot_matrix(2,1)=Foam::sin(phi_1_total)*Foam::cos(phi_2_total);
  Rot_matrix(2,2)=Foam::cos(phi_1_total)*Foam::cos(phi_2_total);

  //https://en.wikipedia.org/wiki/Rotation_matrix - 
//"represents an extrinsic rotation whose (improper) Euler angles are α, β, γ, about axes x, y, z" - phi1 ==alpha, phi2==beta,phi3==gamma




  if(PopBal[h][cell]>SMALL){// if the phase is there, work out the stiffness contribution

  for(int i = 0 ; i < 3 ; i++ ){
  for(int j = 0 ; j < 3 ; j++ ){
    for(int k = 0 ; k < 3 ; k++ ){
      for(int l = 0 ;l  < 3 ; l++ ){
        

        double value =0.0;
        for(int m = 0 ; m < 3 ; m++ ){
          for(int n = 0 ; n < 3 ; n++ ){
            for(int p = 0 ; p < 3 ; p++ ){
              for(int q = 0 ;q  < 3 ; q++ ){
                    if(Cprime2[Phase2Material[h]](m,n,p,q)>SMALL){
                    value+=Rot_matrix(i,m)*Rot_matrix(j,n)*Rot_matrix(k,p)*Rot_matrix(l,q)*Cprime2[Phase2Material[h]](m,n,p,q);//C(m,n,p,q);
                        }
                    }
                  }
                }
              }

              // value*=PopBal[h][cell];

        Cprime(i,j,k,l)=value;
      
      }
    }
  }
  }

//If grain is there
//Have the stiffness tensor now
//need Dp and omegaP

  Fastor::Tensor<double,3,3> F_D_P(0.0);
  Fastor::Tensor<double,3,3> F_Omega_P(0.0);

for(int s=0;s<N_Slips;s++){

  Fastor::Tensor<double,3> F_n_alpha(0.0);
  Fastor::Tensor<double,3> F_m_alpha(0.0);
  Fastor::Tensor<double,3,3> P_alpha(0.0);
  Fastor::Tensor<double,3,3> omega_alpha(0.0);


F_n_alpha(0) = n_alpha[s].component(0);
F_n_alpha(1) = n_alpha[s].component(1);
F_n_alpha(2) = n_alpha[s].component(2);

F_m_alpha(0) = m_alpha[s].component(0);
F_m_alpha(1) = m_alpha[s].component(1);
F_m_alpha(2) = m_alpha[s].component(2);

Fastor::Tensor<double,3> F_n_alpha_Rot(0.0);
Fastor::Tensor<double,3> F_m_alpha_Rot(0.0);

for(int i=0;i<3;i++){
  for(int j=0;j<3;j++){//Rotated slip information
    F_n_alpha_Rot(i)+=Rot_matrix(i,j)*F_n_alpha(j);
    F_m_alpha_Rot(i)+=Rot_matrix(i,j)*F_m_alpha(j);
  }
}

P_alpha = 0.5*(outer(F_m_alpha_Rot,F_n_alpha_Rot) + outer(F_n_alpha_Rot,F_m_alpha_Rot));
omega_alpha = 0.5*(outer(F_m_alpha_Rot,F_n_alpha_Rot) - outer(F_n_alpha_Rot,F_m_alpha_Rot));

Fastor::Tensor<double,3,3> Schmid_T(0.0);

Schmid_T = outer(F_m_alpha_Rot,F_n_alpha_Rot);
double tau_alpha = 0.0;
for(int i = 0 ; i < 3 ; i++ ){
  for(int j = 0 ; j < 3 ; j++ ){
tau_alpha+=Fsigma(i,j) * Schmid_T(i,j);
  }
}


//Work Out g_alphas
double k_exponent = 5.0;//12.0
double ref_shear_rate = 1e-3;
double ta_s = 148e6;
double tau_0 = 16e6;
double r = 1.1;
double ho = 180e6;

g_alpha_dots[s][cell]*=0.0;

for(int beta=0;beta<N_Slips;beta++){

  double h_of_gamma = ho*Foam::pow(1.0/(Foam::cosh((ho*Acc_Strain[cell])/(ta_s-tau_0))),2.0);
  double h_alpha_beta = 0.0;
  if(beta==s){
    h_alpha_beta = h_of_gamma*(r+(1.0-r));
    }

    else{
      h_alpha_beta = h_of_gamma*(r);
    }

// h_alpha_beta = r*h_of_gamma;

  g_alpha_dots[s][cell]+=(h_alpha_beta)*gamma_alpha_dots[beta][cell];

}

// Info<<g_alphas[s][cell]<<endl;

gamma_alpha_dots[s][cell] = ref_shear_rate*Foam::pow((mag(tau_alpha/g_alphas[s][cell])),k_exponent)*sign(tau_alpha);//source term for RHS of g_alpha equations

sumgamma[cell]+= mag(gamma_alpha_dots[s][cell]);//for next time-step
//Work Out g_alphas


F_D_P+=(P_alpha*gamma_alpha_dots[s][cell]);
F_Omega_P+=(omega_alpha*gamma_alpha_dots[s][cell]);

}//end of loop over slip systems


enum {I,J,K,L};
Fastor::Tensor<double,3,3> sig= Fastor::einsum<Fastor::Index<I,J,K,L>,Fastor::Index<K,L>/*,Fastor::OIndex<I,J>*/>(Cprime,Fepsilon_dot-F_D_P);////Double contraction

sigAv+=PopBal[h][cell]*sig;
Omega_P_Average+=PopBal[h][cell]*F_Omega_P;
D_P_Average+=PopBal[h][cell]*F_D_P;

}//end loop over (if popbal[h]>SMALL)









}//end loop over phase fields











count_index=0;
for(int i = 0 ; i < 3 ; i++ ){
    for(int j = 0 ; j < 3 ; j++ ){

      dsigmadt/*_crystal*/[cell].component(count_index) +=sigAv(i,j);//sig(i,j);
      // epsilon_P_dot[cell].component(count_index) = Fepsilon_P_dot(i,j);

      D_P[cell].component(count_index) = D_P_Average(i,j);//F_D_P(i,j);
      Omega_P[cell].component(count_index) = Omega_P_Average(i,j);//F_Omega_P(i,j);
      // testSIG[cell].component(count_index) =sig(i,j);
  count_index++;
    }

}


}


/*




  double pressure=0.0;


  for(int i = 0 ; i < 3 ; i++ ){
   pressure-=Fsigma(i,i)/3.0; //trace(sigma)/3////

  }
   p[cell]=-pressure;
*/


sumgamma.correctBoundaryConditions();

for (label countS = 0; countS < N_Slips; countS++)
{
g_alpha_dots[countS].correctBoundaryConditions();
gamma_alpha_dots[countS].correctBoundaryConditions();
}


dsigmadt += (sigma & Omega_P) - (Omega_P & sigma);

dsigmadt.correctBoundaryConditions();
// epsilon_P_dot.correctBoundaryConditions();
// Info<<"HERE"<<endl;


// Omega_E = omega_dot - Omega_P;

// Omega_E.correctBoundaryConditions();



sigma=sigma + (1.0/n_sub_steps)*runTime.deltaT()*(dsigmadt);
sigma.correctBoundaryConditions();

// sigma = fvc::smoo(sigma);

// epsilon_dot_P_VM.correctBoundaryConditions();
// Acc_Strain = Acc_Strain + (1.0/n_sub_steps)*runTime.deltaT()*(epsilon_dot_P_VM);
Acc_Strain.correctBoundaryConditions();





// D = D+ (1.0/n_sub_steps)*runTime.deltaT()*(U);
// D.correctBoundaryConditions();
// Foam::label Foam::pimpleLoop::nCorrPimple() const;



























































//TESTING AREA

// Fastor::Tensor<double,3,3> teststrain={{1.0,2.0,3.0},
//                                  {4.0,5.0,6.0},
//                                  {7.0,8.0,9.0}};

// Fastor::Tensor<double,3,3> sigdiv={{2.0,4.0,6.0},
//                                  {8.0,10.0,12.0},
//                                  {14.0,16.0,18.0}};





// Fastor::Tensor<double,3,3,3,3> Stiffness_matrix(0.0);

//   // Stiffness_matrix(0,0,0,0)=2.0;
//   // Stiffness_matrix(0,0,1,1)=0.0;
//   // Stiffness_matrix(0,0,2,2)=0.0;
//   // Stiffness_matrix(1,1,0,0)=0.0;
//   // Stiffness_matrix(1,1,1,1)=2.0;
//   // Stiffness_matrix(1,1,2,2)=0.0;
//   // Stiffness_matrix(2,2,0,0)=0.0;
//   // Stiffness_matrix(2,2,1,1)=0.0;
//   // Stiffness_matrix(2,2,2,2)=2.0;

//   Stiffness_matrix(0,0,0,0)=2.0;
//   Stiffness_matrix(0,0,0,1)=3.0;
//   Stiffness_matrix(0,0,0,2)=4.0;
//   Stiffness_matrix(0,0,1,0)=5.0;
//   Stiffness_matrix(0,0,1,1)=6.0;
//   Stiffness_matrix(0,0,1,2)=7.0;
//   Stiffness_matrix(0,0,2,0)=8.0;
//   Stiffness_matrix(0,0,2,1)=9.0;
//   Stiffness_matrix(0,0,2,2)=10.0;

//   Stiffness_matrix(0,1,0,0)=12.0;
//   Stiffness_matrix(0,1,0,1)=13.0;
//   Stiffness_matrix(0,1,0,2)=14.0;
//   Stiffness_matrix(0,1,1,0)=15.0;
//   Stiffness_matrix(0,1,1,1)=16.0;
//   Stiffness_matrix(0,1,1,2)=17.0;
//   Stiffness_matrix(0,1,2,0)=18.0;
//   Stiffness_matrix(0,1,2,1)=19.0;
//   Stiffness_matrix(0,1,2,2)=110.0;

//   Stiffness_matrix(2,2,0,0)=12.0;
//   Stiffness_matrix(2,2,0,1)=13.0;
//   Stiffness_matrix(2,2,0,2)=14.0;
//   Stiffness_matrix(2,2,1,0)=15.0;
//   Stiffness_matrix(2,2,1,1)=16.0;
//   Stiffness_matrix(2,2,1,2)=17.0;
//   Stiffness_matrix(2,2,2,0)=18.0;
//   Stiffness_matrix(2,2,2,1)=19.0;
//   Stiffness_matrix(2,2,2,2)=110.0;


// enum {I,J,K,L};
// Fastor::Tensor<double,3,3> TEST_OUT= Fastor::einsum<Fastor::Index<I,J,K,L>,Fastor::Index<K,L>/*,Fastor::OIndex<I,J>*/>(Stiffness_matrix,teststrain);//

// Fastor::print("TEST_stress", TEST_OUT);

// exit(0);





// for(int i = 0 ; i < 3 ; i++ ){
//   for(int j = 0 ; j < 3 ; j++ ){
//     for(int k = 0 ; k < 3 ; k++ ){
//       for(int l = 0 ;l  < 3 ; l++ ){

//         double value =0.0;
//         for(int m = 0 ; m < 3 ; m++ ){
//           for(int n = 0 ; n < 3 ; n++ ){
//             for(int p = 0 ; p < 3 ; p++ ){
//               for(int q = 0 ;q  < 3 ; q++ ){

//                     value+=Rot_matrix(i,m)*Rot_matrix(j,n)*Rot_matrix(k,p)*Rot_matrix(l,q)*C(m,n,p,q);

//                     }
//                   }
//                 }
//               }

//         Cprime(i,j,k,l)=value;
      
//       }
//     }
//   }
// }

  // exit(0);

}

