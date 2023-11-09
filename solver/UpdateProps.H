{

// Info<<n_sub_steps<<endl;

phi=fvc::flux(U);

rhoPhi = fvc::interpolate(rho)*fvc::flux(U);

// divU=fvc::div(U);

volTensorField TgradU(fvc::grad(U)().T());
TgradU.correctBoundaryConditions(); // already called by line above


omega_dot = 0.5*(fvc::grad(U)-TgradU); // this is the "symm" function, i.e. symm(gradU)
omega_dot.correctBoundaryConditions(); // already called by line above
epsilon_dot = 0.5*(fvc::grad(U)+TgradU) ;//- (1.0/3.0)*(fvc::div(phi))*Ident; // this looks like dev(gradU)
epsilon_dot.correctBoundaryConditions(); // already called by line above


dsigmadt = - ((sigma & omega_dot) - (omega_dot & sigma));
dsigmadt.correctBoundaryConditions();
// dsigmadt_crystal*=0.0;

epsilon_P_dot*=0.0;

p*=0.0;

Euler_Angles*=0.0;
epsilon_dot_P_VM*=0.0;


forAll(mesh.C(), cell)
{



Fastor::Tensor<double,3,3,3,3> Cprime(0.0); //global stiffness at this cell



for (label h = 0; h < N_Ori; h++)
{
  if(PopBal[h][cell]>SMALL){// if the phase is there, work out the stiffness contribution
//  Info<<"Phase: "<<h<<"\t "<<PopBal[h][cell]<<endl;
double phi_1_total = phi1[h] + omega_tot[cell].component(5);

double phi_2_total = phi2[h] + omega_tot[cell].component(2);

double phi_3_total = phi3[h] + omega_tot[cell].component(1);

Euler_Angles[cell].component(0) += PopBal[h][cell]*phi_1_total;
Euler_Angles[cell].component(1) += PopBal[h][cell]*phi_2_total;
Euler_Angles[cell].component(2) += PopBal[h][cell]*phi_3_total;

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

        Cprime(i,j,k,l)+=PopBal[h][cell]*value;
      
      }
    }
  }
}



  }
}




// // Info<<"hang: "<<cell<<endl;
// for (label h = 0; h < N_Ori; h++)
// {

// Cprime += PopBal[h][cell]*Cprime2[h];

// }



  Fastor::Tensor<double,3,3> Fsigma(0.0);
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



Fastor::Tensor<double,3,3> Fepsilon_P_dot(0.0);

{


double m = 20.0;
double epsilon_P_dot_zero = 1e-3;
// double sigma_Y = 240e6;//1e9;





  Fastor::Tensor<double,3,3> sigma_D(0.0);

  

  double pressure=0.0;
  double sigma_VM=0.0;//von mises,J2
  double e_dot_DDP_s = 0.0; //double inner product of e_dot and deviatoric part of stress
  double s_DDP_s = 0.0; //double inner product of deviatoric part of stress and deviatoric part of stress
  // double epsilon_dot_P_VM=0.0;
  // epsilon_dot_P_VM[cell] = 0.0;

  for(int i = 0 ; i < 3 ; i++ ){
   pressure-=Fsigma(i,i)/3.0; //trace(sigma)/3////

  }
   p[cell]=-pressure;



  for(int i = 0 ; i < 3 ; i++ ){
    for(int j = 0 ; j < 3 ; j++ ){
      sigma_D(i,j) = Fsigma(i,j);
      if(i==j){
      sigma_D(i,j) = Fsigma(i,j)+pressure;
      }
    sigma_VM += sigma_D(i,j)*sigma_D(i,j);
    s_DDP_s += sigma_D(i,j)*sigma_D(i,j);
    e_dot_DDP_s += Fepsilon_dot(i,j)*sigma_D(i,j);


    }
  }
sigma_VM=Foam::sqrt(sigma_VM*(3.0/2.0));

// Info<<sigma_VM<<"\t"<<(Foam::pow(sigma_Y,2.0)/3.0)<<endl;




if(Hardening){

  double log_term = 0.0;
  if(epsilon_dot_P_VM[cell]>SMALL){
    log_term=Foam::log(epsilon_dot_P_VM[cell]);
  }

sigma_Y[cell]=(sigma_Y_0.value()+(K_0*Foam::pow(max(Acc_Strain[cell],0.0),K_n)))*(K_B*Foam::exp(-(K_Beta_0-K_Beta_1*log_term)*T[cell]));

}
else{
  sigma_Y[cell]=sigma_Y_0.value();
}




for(int i = 0 ; i < 3 ; i++ ){
    for(int j = 0 ; j < 3 ; j++ ){

// if(sigma_VM>=(sigma_Y[cell]) && (e_dot_DDP_s > SMALL)){
//   // Info<<"here"<<endl;
//         Fepsilon_P_dot(i,j)= (e_dot_DDP_s/s_DDP_s)*sigma_D(i,j); //epsilon_dot_P_VM[cell]*sigma_D(i,j)/sigma_VM;
//         }
//         else{
//           Fepsilon_P_dot(i,j)=0.0;
//         }

        // Pratheeks plasticity implementation
        if(sigma_VM>SMALL/*sigma_Y[cell]*/){
        Fepsilon_P_dot(i,j)=epsilon_P_dot_zero*Foam::pow((sigma_VM/sigma_Y[cell]),m)*sigma_D(i,j)/sigma_VM;
        }
        

    }
  }

double eps_P_dot_DDP_eps_P_dot = 0.0;
for(int i = 0 ; i < 3 ; i++ ){
    for(int j = 0 ; j < 3 ; j++ ){
      eps_P_dot_DDP_eps_P_dot+= Fepsilon_P_dot(i,j)*Fepsilon_P_dot(i,j);
    }
}

eps_P_dot_DDP_eps_P_dot*=(2.0/3.0);

epsilon_dot_P_VM[cell]=Foam::sqrt(eps_P_dot_DDP_eps_P_dot);


// Info<<"sig_y: "<<sigma_Y[cell]<<" Acc_strain: "<<Acc_Strain[cell]<<endl;

// epsilon_dot_P_VM[cell]=epsilon_P_dot_zero*Foam::pow((sigma_VM/sigma_Y),m);



//epsilon_dot_P_VM will need to be integrated dt, and the accumilated rate convected once we have hardening

// epsilon_P_dot_zero = constant, inpyut parameter.
// sigma_Y= constant, inpyut parameter. // will update later
// m= constant, inpyut parameter.


// for(int i = 0 ; i < 3 ; i++ ){
//     for(int j = 0 ; j < 3 ; j++ ){
// // Info<<"H1: "<<sigma_VM<<endl;
// if(sigma_VM>SMALL){
//         Fepsilon_P_dot(i,j)=epsilon_dot_P_VM[cell]*sigma_D(i,j)/sigma_VM;
//         }
// // Info<<"H1"<<endl;
//     }
//   }



}



enum {I,J,K,L};
Fastor::Tensor<double,3,3> sig= Fastor::einsum<Fastor::Index<I,J,K,L>,Fastor::Index<K,L>/*,Fastor::OIndex<I,J>*/>(Cprime,Fepsilon_dot-Fepsilon_P_dot);////Double contraction



count_index=0;
for(int i = 0 ; i < 3 ; i++ ){
    for(int j = 0 ; j < 3 ; j++ ){

      dsigmadt/*_crystal*/[cell].component(count_index) +=sig(i,j);
      epsilon_P_dot[cell].component(count_index) = Fepsilon_P_dot(i,j);
      // testSIG[cell].component(count_index) =sig(i,j);
  count_index++;
    }

}


}

dsigmadt.correctBoundaryConditions();
epsilon_P_dot.correctBoundaryConditions();
// Info<<"HERE"<<endl;




sigma=sigma + (1.0/n_sub_steps)*runTime.deltaT()*(dsigmadt);
sigma.correctBoundaryConditions();

// sigma = fvc::smoo(sigma);

epsilon_dot_P_VM.correctBoundaryConditions();
// Acc_Strain = Acc_Strain + (1.0/n_sub_steps)*runTime.deltaT()*(epsilon_dot_P_VM);
Acc_Strain.correctBoundaryConditions();





D = D+ (1.0/n_sub_steps)*runTime.deltaT()*(U);
D.correctBoundaryConditions();
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

