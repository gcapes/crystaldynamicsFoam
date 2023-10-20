{


phi=fvc::flux(U);

rhoPhi = fvc::interpolate(rho)*fvc::flux(U);

// divU=fvc::div(U);

volTensorField TgradU(fvc::grad(U)().T());
TgradU.correctBoundaryConditions(); // already called by line above



// volTensorField TgradD = (fvc::grad(D)().T());
// TgradD.correctBoundaryConditions();
// gradU = fvc::grad(fvc::interpolate(U));

// gradU=fvc::grad(U);

// tensor Ident(1,0,0,0,1,0,0,0,1);

omega_dot = 0.5*(fvc::grad(U)-TgradU); // this is the "symm" function, i.e. symm(gradU)
omega_dot.correctBoundaryConditions(); // already called by line above
epsilon_dot = 0.5*(fvc::grad(U)+TgradU) ;//- (1.0/3.0)*(fvc::div(phi))*Ident; // this looks like dev(gradU)
epsilon_dot.correctBoundaryConditions(); // already called by line above

// Info<<I<<endl;

// Info<<"IM"<<Ident<<endl;

// tredot = tr(epsilon_dot);


// epsilon = 0.5*(fvc::grad(D)-TgradD);
// epsilon.correctBoundaryConditions();
// dsdt = 2.0*shear_mod*(epsilon_dot/* - epsilon_P_dot*/);// - ((s & omega_dot) - (omega_dot & s));
// dsdt *=0.0;
dsigmadt = - ((sigma & omega_dot) - (omega_dot & sigma));

// dsigmadt_crystal*=0.0;

epsilon_P_dot*=0.0;

p*=0.0;

// p=tr(sigma);




// sigma && sigma;// double dot product

// Info<<"HERE2"<<endl;


forAll(mesh.C(), cell)
{



// Info<<"epsilon_dot: "<<epsilon_dot[cell]<<endl;

// Fastor::Tensor<double,3,3,3,3> C(0.0);

// // C(0,0,0,0) = 8.0e10;
// // C(0,1,0,1) = 8.0e10;
// // C(0,2,0,2) = 8.0e10;

// // C(1,0,1,0) = 8.0e10;
// // C(1,1,1,1) = 8.0e10;
// // C(1,2,1,2) = 8.0e10;

// // C(2,0,2,0) = 8.0e10;
// // C(2,1,2,1) = 8.0e10;
// // C(2,2,2,2) = 8.0e10;



//   double C11 = 261e9;
//   double C12 = 151e9;
//   double C44 = 132e9;




Fastor::Tensor<double,3,3,3,3> Cprime(0.0);

// // Info<<"hang: "<<cell<<endl;
for (label h = 0; h < N_Ori; h++)
{

Cprime += PopBal[h][cell]*Cprime2[h];

}

// Cprime = Cprime2[0];

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


double m = 1.0;
double epsilon_P_dot_zero = 1e-5;
double sigma_Y = 1e9;





  Fastor::Tensor<double,3,3> sigma_D(0.0);

  

  double pressure=0.0;
  double sigma_VM=0.0;
  // double epsilon_dot_P_VM=0.0;
  // epsilon_dot_P_VM[cell] = 0.0;

  for(int i = 0 ; i < 3 ; i++ ){
   pressure-=Fsigma(i,i)/3.0; //trace(sigma)/3////

  }
   p[cell]=-pressure;



  for(int i = 0 ; i < 3 ; i++ ){
    for(int j = 0 ; j < 3 ; j++ ){

      if(i==j){
      sigma_D(i,j) = Fsigma(i,j)+pressure;
      }
    sigma_VM += sigma_D(i,j)*sigma_D(i,j);


    }
  }
sigma_VM=Foam::sqrt(sigma_VM/2.0);


epsilon_dot_P_VM[cell]=epsilon_P_dot_zero*Foam::pow((sigma_VM/sigma_Y),m);

//epsilon_dot_P_VM will need to be integrated dt, and the accumilated rate convected once we have hardening

// epsilon_P_dot_zero = constant, inpyut parameter.
// sigma_Y= constant, inpyut parameter. // will update later
// m= constant, inpyut parameter.


for(int i = 0 ; i < 3 ; i++ ){
    for(int j = 0 ; j < 3 ; j++ ){
// Info<<"H1: "<<sigma_VM<<endl;
if(sigma_VM>SMALL){
        Fepsilon_P_dot(i,j)=epsilon_dot_P_VM[cell]*sigma_D(i,j)/sigma_VM;
        }
// Info<<"H1"<<endl;
    }
  }



}

// exit(0);

// Info<<"epsilon_dot: "<<epsilon_dot[cell]<<endl;

// Fastor::print("epsilon_dot:", Fepsilon_dot);


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
// Info<<"HERE"<<endl;



// dsigmadt_crystal.correctBoundaryConditions();
// dsigmadt+=dsigmadt_crystal;

// dsdt_crystal.write();
// exit(0);

// Info<<tr(dsigmadt)<<endl;
// trsig = tr(dsigmadt);


sigma=sigma + (1.0/n_sub_steps)*runTime.deltaT()*(dsigmadt);
sigma.correctBoundaryConditions();

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

