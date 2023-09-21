{

// rhoPhi = fvc::interpolate(rho)*fvc::flux(U);
// phi=fvc::flux(U);

dimensionedScalar shear_mod ("shear_mod", dimensionSet(1,-1,-2,0,0,0,0),80e9);

volTensorField TgradU = (fvc::grad(U)().T());

omega_dot = 0.5*(fvc::grad(U)-TgradU);
omega_dot.correctBoundaryConditions();
epsilon_dot = 0.5*(fvc::grad(U)-TgradU) - (1.0/3.0)*(fvc::div(phi))*I;
epsilon_dot.correctBoundaryConditions();
// dsdt = 2.0*shear_mod*(epsilon_dot/* - epsilon_P_dot*/) - ((s & omega_dot) - (omega_dot & s));
// dsdt *=0.0;
dsdt = - ((s & omega_dot) - (omega_dot & s));




dsdt_crystal*=0.0;
forAll(mesh.C(), cell)
{



Fastor::Tensor<double,3,3,3,3> C(0.0);

  double C_11 = 10.0e10;
  double C_12 = 9.0e10;
  double C_66 = 8.0e10;


  C(0,0,0,0)=C_11;
  C(1,1,1,1)=C_11;
  C(2,2,2,2)=C_11;

  C(0,0,1,1)=C_12;
  C(0,0,2,2)=C_12;
  C(1,1,2,2)=C_12;
  C(1,1,0,0)=C_12;
  C(2,2,0,0)=C_12;
  C(2,2,1,1)=C_12;

  C(1,2,1,2)=C_66;
  C(0,2,0,2)=C_66;
  C(0,1,0,1)=C_66;


  //apply the symmetries
  // double ct=0.0;
  // for(int i = 0 ; i < 3 ; i++ ){
  //   for(int j = 0 ; j < 3 ; j++ ){
  //     for(int k = 0 ; k < 3 ; k++ ){
  //       for(int l = 0 ;l  < 3 ; l++ ){

  //         C(i,j,k,l)=ct;

  // //         if(C(j,i,k,l)==0.0){
  // //                  C(j,i,k,l) =C(i,j,k,l);//ct;
  // //                }
  // //           if(C(i,j,l,k) ==0.0){
  // //                 C(i,j,l,k)  =C(i,j,k,l);//ct;
  // //                }

  // //           if(C(j,i,l,k)==0.0){
  // //                 C(j,i,l,k) =C(i,j,k,l);//ct;
  // //                 }
            
  // //           // C(i,j,l,k) =C(i,j,k,l);//ct;
  // //           //  C(j,i,l,k) =C(i,j,k,l);//ct;
  // //           // //C(k,l,i,j) =C(i,j,k,l);//ct;

  //          ct += 1.0;
          

  //        }
  //     }
  //   }
  // }
  //finish symmetry operations
//   Fastor::print("C: ",C);
// exit(0);
//   Info<<"epsilon_dot: "<<epsilon_dot[cell]<<endl;

//  Info<<"epsilon_dot 4: "<<epsilon_dot[cell].component(4)<<endl;
//  Info<<"epsilon_dot 5: "<<epsilon_dot[cell].component(5)<<endl;

  Fastor::Tensor<double,3,3> Fepsilon_dot(0.0);
  int count_index=0;
  for(int i = 0 ; i < 3 ; i++ ){
    for(int j = 0 ; j < 3 ; j++ ){

      Fepsilon_dot(i,j)=epsilon_dot[cell].component(count_index);

      count_index++;
    }
  }


  

// Fepsilon_dot(0) = epsilon_dot[cell].component(0);
// Fepsilon_dot(1) = epsilon_dot[cell].component(4);
// Fepsilon_dot(2) = epsilon_dot[cell].component(8);


// Fastor::print("Fepsilon_dot 1,1: ",Fepsilon_dot);
// Fastor::print("Fepsilon_dot 1,1: ",Fepsilon_dot(1,1));
// Fastor::print("Fepsilon_dot 1,2: ",Fepsilon_dot(1,2));

// Fastor::Tensor<double,9,9> reshapedC(0.0) ;
// reshapedC = Fastor::reshape<9,9>(C);

// Fastor::print("reshapedC: ",reshapedC);

// exit(0);
// Fastor::Tensor<double,9> flatedot(0.0) ;
// flatedot= Fastor::flatten(Fepsilon_dot);

// Fastor::Tensor<double,9,1> testmul(0.0);
// testmul= Fastor::matmul(reshapedC,flatedot);

enum {I,J,K,L};
Fastor::Tensor<double,3,3> sig= Fastor::einsum<Fastor::Index<I,J,K,L>,Fastor::Index<K,L>/*,Fastor::OIndex<I,J>*/>(C,Fepsilon_dot);//



count_index=0;
for(int i = 0 ; i < 3 ; i++ ){
    for(int j = 0 ; j < 3 ; j++ ){

      dsdt_crystal[cell].component(count_index) =sig(i,j);
  count_index++;
    }

}


}
dsdt_crystal.correctBoundaryConditions();
dsdt+=dsdt_crystal;




dsdt.correctBoundaryConditions();
s=s + runTime.deltaT()*(dsdt);

    // solve
    // (
    //     fvm::ddt(s) ==
    //     dsdt
    // );

s.correctBoundaryConditions();



// for (label count = 0; count < N_Ori; count++)
// {}




  





// Fastor::Tensor<double,3,3> Id3={{1.0,0.0,0.0},
//                                  {0.0,1.0,0.0},
//                                  {0.0,0.0,1.0}};



// // Id3(0,0)=1.0;
// // Id3(1,1)=1.0;
// // Id3(2,2)=1.0;

// print("Identity test: ",Id3);


// Fastor::Tensor<double,3,3,3,3> Stiffness_matrix(0.0);

//   Stiffness_matrix(0,0,0,0)=1.0;
//   Stiffness_matrix(0,0,1,1)=1.0;
//   Stiffness_matrix(0,0,2,2)=1.0;
//   Stiffness_matrix(1,1,0,0)=1.0;
//   Stiffness_matrix(1,1,1,1)=1.0;
//   Stiffness_matrix(1,1,2,2)=1.0;
//   Stiffness_matrix(2,2,0,0)=1.0;
//   Stiffness_matrix(2,2,1,1)=1.0;
//   Stiffness_matrix(2,2,2,2)=1.0;


//    Fastor::print("TEST: ", Stiffness_matrix);

//   Fastor::Tensor<double,9,9> reshaped = Fastor::reshape<9,9>(Stiffness_matrix);

// Fastor::print("TEST_reshaped: ", reshaped);

// Fastor::Tensor<double,9> flatI = Fastor::flatten(Id3);

// Fastor::print("TEST_flatI ", flatI);


// Fastor::Tensor<double,9,1> testmul= Fastor::matmul(reshaped,flatI);//double dot product of stiffness matrix with some rank 2 tensor

// Fastor::print("RESULT ", testmul);

// Fastor::Tensor<double,3,3> niceformat= Fastor::reshape<3,3>(testmul);

// Fastor::print("RESULT_nice_format ", niceformat);



// Fastor::Tensor<double,3,3> test1={{10.0,12.0,13.0},
//                                  {14.0,15.0,16.0},
//                                  {17.0,18.0,19.0}};

// Fastor::Tensor<double,3,3> test2={{1.0,2.0,3.0},
//                                  {4.0,5.0,6.0},
//                                  {7.0,8.0,9.0}};

// Fastor::print(Fastor::matmul(test1,test2));//dot product 
// /*

// Fastor::Tensor<double,3,3> Id3={{1.0,0.0,0.0},
//                                  {0.0,1.0,0.0},
//                                  {0.0,0.0,1.0}};




// Fastor::Tensor<double,3,3,3,3> Stiffness_matrix(0.0);

//   Stiffness_matrix(0,0,0,0)=1.0;
//   Stiffness_matrix(0,0,1,1)=1.0;
//   Stiffness_matrix(0,0,2,2)=1.0;
//   Stiffness_matrix(1,1,0,0)=1.0;
//   Stiffness_matrix(1,1,1,1)=1.0;
//   Stiffness_matrix(1,1,2,2)=1.0;
//   Stiffness_matrix(2,2,0,0)=1.0;
//   Stiffness_matrix(2,2,1,1)=1.0;
//   Stiffness_matrix(2,2,2,2)=1.0;


// enum {i,j,k,l,m,n};
// enum {I,J,K,L,M,N};
// Fastor::Tensor<double,3,3,3> A; 
// Fastor::Tensor<double,3,3,3,3> B;
// // fill A and B
// // A.random(); B.random();
 

// Fastor::Tensor<double,3,3> Cout= Fastor::einsum<Fastor::Index<0,1,2,3>,Fastor::Index<0,1>>(Stiffness_matrix,Id3);

// Fastor::print("TESTING DOT: ",Cout);
// exit(0);


// // or simply using constant integer numbers instead of an enum
// // auto C = einsum<Index<0,1,2>,Index<1,3,4,5>>(A,B);
// Fastor::print("TESTING DOT: ",C);
// */

// Fastor::Tensor<double,3,3> Rot_matrix={{1.0,0.0,0.0},
//                                  {0.0,1.0,0.0},
//                                  {0.0,0.0,1.0}};

// Fastor::Tensor<double,3,3,3,3> Cprime(0.0);

// double C_11 = 10.0e10;
// double C_12 = 9.0e10;
// double C_66 = 8.0e10;


// Fastor::Tensor<double,3,3,3,3> C(0.0);//stiffness tensor in crystal reference
//   C(0,0,0,0)=C_11;
//   C(1,1,1,1)=C_11;
//   C(2,2,2,2)=C_11;

//   C(0,0,1,1)=C_12;
//   C(0,0,2,2)=C_12;
//   C(1,1,2,2)=C_12;
//    C(1,1,0,0)=C_12;
//    C(2,2,0,0)=C_12;
//    C(2,2,1,1)=C_12;

//   C(1,2,1,2)=C_66;
//   C(0,2,0,2)=C_66;
//   C(0,1,0,1)=C_66;



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