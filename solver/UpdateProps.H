{

dimensionedScalar shear_mod ("shear_mod", dimensionSet(1,-1,-2,0,0,0,0),80e9);

volTensorField TgradU = (fvc::grad(U)().T());

omega_dot = 0.5*(fvc::grad(U)-TgradU);

epsilon_dot = 0.5*(fvc::grad(U)-TgradU) - (1.0/3.0)*(fvc::div(phi))*I;

dsdt = 2.0*shear_mod*(epsilon_dot/* - epsilon_P_dot*/) - ((s & omega_dot) - (omega_dot & s));

s=s + runTime.deltaT()*(dsdt);

    // solve
    // (
    //     fvm::ddt(s) ==
    //     dsdt
    // );




/*
tensor T1(1, 2, 3, 4, 5, 6, 7, 8, 9);

tensor T2(10, 11, 12, 13, 14, 15, 16, 17, 18);

tensor T3 = (T1 & T2);

tensor T4 = (T2 & T1);

Info<<"T3 ="<<T3<<endl;

Info<<"T4 ="<<T4<<endl;

exit(0);

*/


/*
  Fastor::Tensor<double> my_scalar = 2; // this is a scalar (double) with a value of 2
  // output this on the screen
  Fastor::print("my scalar",my_scalar); // print is a built-in Fastor function for printing

  Fastor::Tensor<double,3> my_vector = {1,2,3}; // this is a vector of 3 doubles with values {1,2,3}
  Fastor::print("my vector",my_vector);
  Fastor::Tensor<float,3,3> my_matrix = {{1,2,3},
                                 {4,5,6},
                                 {7,8,9}}; // this a 3x3 matrix of floats with values 1...9
  Fastor::print("my matrix",my_matrix);

  Fastor::Tensor<int,2,2,2,2,2> array_5d; // this a 5D array of ints with dimension 2x2x2x2x2
  array_5d.iota(1); // fill the 5d_array with sequentially ascending numbers from 1 to 2x2x2x2x2=32
  Fastor::print("my 5D array", array_5d);
*/


  





Fastor::Tensor<double,3,3> Id3={{1.0,0.0,0.0},
                                 {0.0,1.0,0.0},
                                 {0.0,0.0,1.0}};



// Id3(0,0)=1.0;
// Id3(1,1)=1.0;
// Id3(2,2)=1.0;

print("Identity test: ",Id3);


Fastor::Tensor<double,3,3,3,3> Stiffness_matrix(0.0);

  Stiffness_matrix(0,0,0,0)=1.0;
  Stiffness_matrix(0,0,1,1)=1.0;
  Stiffness_matrix(0,0,2,2)=1.0;
  Stiffness_matrix(1,1,0,0)=1.0;
  Stiffness_matrix(1,1,1,1)=1.0;
  Stiffness_matrix(1,1,2,2)=1.0;
  Stiffness_matrix(2,2,0,0)=1.0;
  Stiffness_matrix(2,2,1,1)=1.0;
  Stiffness_matrix(2,2,2,2)=1.0;


   Fastor::print("TEST: ", Stiffness_matrix);

  Fastor::Tensor<double,9,9> reshaped = Fastor::reshape<9,9>(Stiffness_matrix);

Fastor::print("TEST_reshaped: ", reshaped);

Fastor::Tensor<double,9> flatI = Fastor::flatten(Id3);

Fastor::print("TEST_flatI ", flatI);


Fastor::Tensor<double,9,1> testmul= Fastor::matmul(reshaped,flatI);//double dot product of stiffness matrix with some rank 2 tensor

Fastor::print("RESULT ", testmul);

Fastor::Tensor<double,3,3> niceformat= Fastor::reshape<3,3>(testmul);

Fastor::print("RESULT_nice_format ", niceformat);



Fastor::Tensor<double,3,3> test1={{10.0,12.0,13.0},
                                 {14.0,15.0,16.0},
                                 {17.0,18.0,19.0}};

Fastor::Tensor<double,3,3> test2={{1.0,2.0,3.0},
                                 {4.0,5.0,6.0},
                                 {7.0,8.0,9.0}};

Fastor::print(Fastor::matmul(test1,test2));//dot product 
/*
enum {i,j,k,l,m,n};
enum {I,J,K,L,M,N};
Fastor::Tensor<double,3,3,3> A; 
Fastor::Tensor<double,3,3,3,3> B;
// fill A and B
// A.random(); B.random();
 
Fastor::Tensor<double,3,3,3,3,3> C= Fastor::einsum<Fastor::Index<I,J,K>,Fastor::Index<J,L,M,N>>(A,B);


// or simply using constant integer numbers instead of an enum
// auto C = einsum<Index<0,1,2>,Index<1,3,4,5>>(A,B);
Fastor::print("TESTING DOT: ",C);
*/

Fastor::Tensor<double,3,3> Rot_matrix={{1.0,0.0,0.0},
                                 {0.0,1.0,0.0},
                                 {0.0,0.0,1.0}};

Fastor::Tensor<double,3,3,3,3> Cprime(0.0);

Fastor::Tensor<double,3,3,3,3> C(0.0);//stiffness tensor in crystal reference
  C(0,0,0,0)=1.0;
  C(1,1,1,1)=1.0;
  C(2,2,2,2)=1.0;


for(int i = 0 ; i < 3 ; i++ ){
  for(int j = 0 ; j < 3 ; j++ ){
    for(int k = 0 ; k < 3 ; k++ ){
      for(int l = 0 ;l  < 3 ; l++ ){

        double value =0.0;
        for(int m = 0 ; m < 3 ; m++ ){
          for(int n = 0 ; n < 3 ; n++ ){
            for(int p = 0 ; p < 3 ; p++ ){
              for(int q = 0 ;q  < 3 ; q++ ){

                    value+=Rot_matrix(i,m)*Rot_matrix(j,n)*Rot_matrix(k,p)*Rot_matrix(l,q)*C(m,n,p,q);

                    }
                  }
                }
              }

        Cprime(i,j,k,l)=value;
      
      }
    }
  }
}

  exit(0);

}