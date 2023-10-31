Info<< "Reading physicalProperties\n" << endl;

IOdictionary physicalProperties
(
    IOobject
    (
        "physicalProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);

IOdictionary solutionDict
(
    IOobject
    (
        "fvSolution",
        runTime.system(), //folder where the fvSolution file is located
        runTime,
        IOobject::MUST_READ
    )
);

dictionary& PimpleDICT = solutionDict.subDict("PIMPLE");

int n_sub_steps(readLabel(PimpleDICT.lookup("nOuterCorrectors")));

// Info<<"TREST: "<<n_sub_steps<<endl;

// exit(0);

dimensionedScalar pi("pi", dimensionSet(0, 0, 0, 0, 0), constant::mathematical::pi);



dimensionedScalar Cv("Cv",dimensionSet(0, 2, -2, -1, 0),physicalProperties.lookup("Cv"));
dimensionedScalar kappa("kappa",dimensionSet(1, 1, -3, -1, 0),physicalProperties.lookup("kappa"));
dimensionedScalar rho_0("rho_0",dimensionSet(1, -3, 0, 0, 0),physicalProperties.lookup("rho_0"));

const Switch Hardening(physicalProperties.lookupOrDefault<Switch>("Hardening", false));

// if(Hardening){

scalar K_0(physicalProperties.lookupOrDefault<scalar>("K_0", 0.0));   
scalar K_n(physicalProperties.lookupOrDefault<scalar>("K_n", 0.0));   
scalar K_B(physicalProperties.lookupOrDefault<scalar>("K_B", 0.0));   
scalar K_Beta_0(physicalProperties.lookupOrDefault<scalar>("K_Beta_0", 0.0));   
scalar K_Beta_1(physicalProperties.lookupOrDefault<scalar>("K_Beta_1", 0.0));    

// dimensionedScalar K_0("K_0",dimensionSet(1, -1, -2, 0, 0),physicalProperties.lookup("K_0"));
// dimensionedScalar K_n("K_n",dimensionSet(0, 0, 0, 0, 0),physicalProperties.lookup("K_n"));
// dimensionedScalar K_B("K_B",dimensionSet(0, 0, 0, 0, 0),physicalProperties.lookup("K_B"));
// dimensionedScalar K_Beta_0("K_Beta_0",dimensionSet(0, 0, 0, 0, 0),physicalProperties.lookup("K_Beta_0"));
// dimensionedScalar K_Beta_1("K_Beta_1",dimensionSet(0, 0, 0, 0, 0),physicalProperties.lookup("K_Beta_1"));

// }

dimensionedScalar sigma_Y_0("sigma_Y_0",dimensionSet(1, -1, -2, 0, 0),physicalProperties.lookup("sigma_Y_0"));


Info<< "Reading field U\n" << endl;
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

volVectorField D
(
    IOobject
    (
        "D",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedVector("D",dimensionSet(0, 1, 0, 0, 0),vector::zero)
);

Info<< "Reading field Temperature\n" << endl;
volScalarField T
(
    IOobject
    (
        "T",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

volScalarField K
(
    IOobject
    (
        "K",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("K",dimensionSet(0, 2, -2, 0, 0),0.0)
);


// volScalarField trsig
// (
//     IOobject
//     (
//         "trsig",
//         runTime.timeName(),
//         mesh,
//         IOobject::READ_IF_PRESENT,
//         IOobject::AUTO_WRITE
//     ),
//     mesh,
//     dimensionedScalar("trsig",dimensionSet(1, -1, -3, 0, 0),0.0)
// );

// volScalarField tredot
// (
//     IOobject
//     (
//         "tredot",
//         runTime.timeName(),
//         mesh,
//         IOobject::READ_IF_PRESENT,
//         IOobject::AUTO_WRITE
//     ),
//     mesh,
//     dimensionedScalar("tredot",dimensionSet(0, 0, -1, 0, 0),0.0)
// );

// volVectorField TESTADV
// (
//     IOobject
//     (
//         "TESTADV",
//         runTime.timeName(),
//         mesh,
//         IOobject::READ_IF_PRESENT,
//         IOobject::AUTO_WRITE
//     ),
//     mesh,
//     dimensionedVector("TESTADV",dimensionSet(1, -2, -2, 0, 0),vector::zero),
//                 zeroGradientFvPatchVectorField::typeName
// );

// volScalarField DP
// (
//     IOobject
//     (
//         "DP",
//         runTime.timeName(),
//         mesh,
//         IOobject::READ_IF_PRESENT,
//         IOobject::AUTO_WRITE
//     ),
//     mesh,
//     dimensionedScalar("DP",dimensionSet(1, -1, -3, 0, 0),0.0),
//                 zeroGradientFvPatchScalarField::typeName
// );

// volScalarField DP2
// (
//     IOobject
//     (
//         "DP2",
//         runTime.timeName(),
//         mesh,
//         IOobject::READ_IF_PRESENT,
//         IOobject::AUTO_WRITE
//     ),
//     mesh,
//     dimensionedScalar("DP2",dimensionSet(1, -1, -3, 0, 0),0.0),
//                 zeroGradientFvPatchScalarField::typeName
// );

volScalarField rho
(
    IOobject
    (
        "rho",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("rho",dimensionSet(1, -3, 0, 0, 0),rho_0.value()),
                zeroGradientFvPatchScalarField::typeName
);

volScalarField p
(
    IOobject
    (
        "p",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("p",dimensionSet(1, -1, -2, 0, 0),0.0),
                zeroGradientFvPatchScalarField::typeName
);


surfaceScalarField rhoPhi
            (
                IOobject
                (
                    "rhoPhi",
                    runTime.timeName(),
                    mesh
                ),
                mesh,
                dimensionedScalar(dimMass/dimTime, 0)
            );



// volScalarField divU
// (
//     IOobject
//     (
//         "divU",
//         runTime.timeName(),
//         mesh,
//         IOobject::NO_READ,
//         IOobject::AUTO_WRITE
//     ),
//     mesh,
//     dimensionedScalar("divU",dimensionSet(0, 0, -1, 0, 0),0.0),
//                 zeroGradientFvPatchScalarField::typeName
// );


volScalarField Sum_Phases
(
    IOobject
    (
        "Sum_Phases",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("Sum_Phases",dimensionSet(0, 0, 0, 0, 0),0.0)
);



// const label N_Ori=20;

IOdictionary GrainProperties
(
    IOobject
    (
        "GrainProperties",    // dictionary name
        runTime.constant(),     // dict is found in "constant"
        mesh,                   // registry for the dict
        IOobject::MUST_READ,    // must exist, otherwise failure
        IOobject::NO_WRITE      // dict is only read by the solver
    )
);


label N_Ori(readLabel(GrainProperties.lookup("N_Phases")));
Info<<"Number of Phase Field Order Parameters: "<<N_Ori<<endl;




  // Ready to read in multiple stiffness matrices

label N_Mats(readLabel(GrainProperties.lookup("N_Materials")));
Info<<"Number of Materials: "<<N_Mats<<endl;


labelList Phase2Material(N_Ori,0);

GrainProperties.lookup("Phase2Material") >> Phase2Material;

for(label c=0; c<N_Ori; c++)
{
    Info<<"Phase["<<c<<"] = Material["<<Phase2Material[c]<<"]"<<endl;
}


Fastor::Tensor<double,3,3,3,3> Cprime2[N_Mats];//array of tensor objects




for(label i = 0; i < N_Mats; i++)
{

    word Mat_dict ("Material." + name(i));
    dictionary Mat = GrainProperties.subDict(Mat_dict);

    Fastor::Tensor<double,3,3,3,3> Cprime_temp_Read(0.0);

    for(int m = 0 ; m < 3 ; m++ ){
          for(int n = 0 ; n < 3 ; n++ ){
            for(int p = 0 ; p < 3 ; p++ ){
              for(int q = 0 ;q  < 3 ; q++ ){

                word Cijkl ("C" + name(m)+ name(n)+ name(p)+ name(q));

                Cprime_temp_Read(m,n,p,q)= Mat.lookupOrDefault<scalar>(Cijkl, 0.0);


              }
            }
          }
    }

Cprime2[i] = Cprime_temp_Read;

}

// Fastor::print(Cprime2[0]);

// Fastor::print(Cprime2[1]);

// exit(0);







PtrList<volScalarField> PopBal(N_Ori);		//pointer to phase fields


for(label count = 0; count < N_Ori; count++)	// how to assign the popbal array of fields
{

    word name_ni ("n." + name(count));

    PopBal.set
    (
        count,
        new volScalarField
        (
            IOobject
            (
                name_ni,
                runTime.timeName(),
                mesh,
                IOobject::READ_IF_PRESENT,//IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            mesh,
            dimensionedScalar("ni",dimensionSet(0, 0, 0, 0, 0),0.0),
            zeroGradientFvPatchScalarField::typeName
        )
    );
}





double phi1[N_Ori]={0.0};

double phi2[N_Ori]={0.0};

double phi3[N_Ori]={0.0};


phi1[0]=0.0;
phi2[0]=0.0;
phi3[0]=0.0;

phi1[1]=2.0*pi.value();
phi2[1]=0.0;
phi3[1]=0.0;

phi1[2]=2.0*pi.value();
phi2[2]=0.0;
phi3[2]=0.0;

phi1[11]=0.25*pi.value();
phi2[11]=0.0;
phi3[11]=0.0;



// tensor rot[N_Ori]={tensor::zero};

// for (label i = 0; i < N_Ori; i++)
// {


//  tensor temprot(
//     Foam::cos(phi2[i])*Foam::cos(phi3[i]), 
//     (Foam::sin(phi1[i])*Foam::sin(phi2[i])*Foam::cos(phi3[i]))-(Foam::cos(phi1[i])*Foam::sin(phi3[i])),//-Foam::cos(phi3)*Foam::sin(phi2), 
//     (Foam::cos(phi1[i])*Foam::sin(phi2[i])*Foam::cos(phi3[i]))+(Foam::sin(phi1[i])*Foam::sin(phi3[i])),//Foam::sin(phi2)*Foam::sin(phi3),
   
//     Foam::cos(phi2[i])*Foam::sin(phi3[i]), 
//     (Foam::sin(phi1[i])*Foam::sin(phi2[i])*Foam::sin(phi3[i]))+(Foam::cos(phi1[i])*Foam::cos(phi3[i])), 
//     (Foam::cos(phi1[i])*Foam::sin(phi2[i])*Foam::sin(phi3[i]))-(Foam::sin(phi1[i])*Foam::cos(phi3[i])),
    
//     -Foam::sin(phi2[i]), 
//     Foam::sin(phi1[i])*Foam::cos(phi2[i]), 
//     Foam::cos(phi1[i])*Foam::cos(phi2[i])
// );//https://en.wikipedia.org/wiki/Rotation_matrix - 
// //"represents an extrinsic rotation whose (improper) Euler angles are α, β, γ, about axes x, y, z" - phi1 ==alpha, phi2==beta,phi3==gamma

// rot[i]=temprot;
// }


// Info<<"TEST: "<<rot[1]<<endl;
// Info<<"TEST: "<<rot[1].component(4)<<endl;
// exit(0);



Fastor::Tensor<double,3,3,3,3> C(0.0);//elasticity matrix in the grain reference




  double C11 = 107.3e9;
  double C12 = 60.9e9;
  double C44 = 128.3e9;

//Cubic

  C(0,0,0,0)=C11;
  C(1,1,1,1)=C11;
  C(2,2,2,2)=C11;

  C(1,0,1,0)=C44;
  C(2,0,2,0)=C44;
  C(0,0,1,1)=C12;
  C(1,0,0,1)=C44;
  C(0,0,2,2)=C12;
  C(2,0,0,2)=C44;
  C(0,1,1,0)=C44;
  C(1,1,0,0)=C12;
  C(0,1,0,1)=C44;
  C(2,1,2,1)=C44;
  C(1,1,2,2)=C12;
  C(2,1,1,2)=C44;
  C(0,2,2,0)=C44;
  C(2,2,0,0)=C12;
  C(1,2,2,1)=C44;
  C(2,2,1,1)=C12;
  C(0,2,0,2)=C44;
  C(1,2,1,2)=C44;


/*
//ISOTROPIC

C(0,0,0,0) = C11;
C(0,1,0,1) = C11;
C(0,2,0,2) = C11;

C(1,0,1,0) = C11;
C(1,1,1,1) = C11;
C(1,2,1,2) = C11;

C(2,0,2,0) = C11;
C(2,1,2,1) = C11;
C(2,2,2,2) = C11;
*/

// Fastor::Tensor<double,3,3,3,3> Cprime2[N_Ori];//array of tensor objects



// for (label h = 0; h < N_Ori; h++)
// {

// Fastor::Tensor<double,3,3,3,3> Cprime_temp(0.0);

//   Fastor::Tensor<double,3,3> Rot_matrix(0.0);

//   Rot_matrix(0,0)=rot[h].component(0);
//   Rot_matrix(0,1)=rot[h].component(1);
//   Rot_matrix(0,2)=rot[h].component(2);
//   Rot_matrix(1,0)=rot[h].component(3);
//   Rot_matrix(1,1)=rot[h].component(4);
//   Rot_matrix(1,2)=rot[h].component(5);
//   Rot_matrix(2,0)=rot[h].component(6);
//   Rot_matrix(2,1)=rot[h].component(7);
//   Rot_matrix(2,2)=rot[h].component(8);


//   // Rotate To Global needReference
//   for(int i = 0 ; i < 3 ; i++ ){
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

//               // value*=PopBal[h][cell];

//         Cprime_temp(i,j,k,l)=value;
      
//       }
//     }
//   }
// }

// // if(i==0){
// //     Cprime_temp(0,0,0,0)=2.35;
// // }
// // if(i==1){
// //     Cprime_temp(0,0,0,0)=247.65;
// // }


// Cprime2[h]=Cprime_temp;
// }
// Info<<"HERE5"<<endl;
//Cprime2[N_Ori] now has all the grain elasticity tensors


// Fastor::print("Test1",Cprime2[0](0,0,0,0));
// Fastor::print("Test1",Cprime2[1](0,0,0,0));
// Fastor::print("Test1",Cprime2[0](0,0,0,1));
// Fastor::print("Test1",Cprime2[2](0,0,0,0));

// exit(0);



// volTensorField Global_Orientation_Tensor
// (
//     IOobject
//     (
//         "Global_Orientation_Tensor",
//         runTime.timeName(),
//         mesh,
//         IOobject::READ_IF_PRESENT,
//         IOobject::AUTO_WRITE
//     ),
//         mesh,
// dimensionedTensor("Global_Orientation_Tensor",dimensionSet(0, 0, 0, 0, 0),tensor::zero),
//                zeroGradientFvPatchTensorField::typeName
// );





// exit(0);

// volScalarField detRM
// (
//     IOobject
//     (
//         "detRM",
//         runTime.timeName(),
//         mesh,
//         IOobject::READ_IF_PRESENT,
//         IOobject::AUTO_WRITE
//     ),
//         mesh,
// dimensionedScalar("detRM",dimensionSet(0, 0, 0, 0, 0),0.0),
//                zeroGradientFvPatchTensorField::typeName
// );


volTensorField omega_dot
(
    IOobject
    (
        "omega_dot",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
        mesh,
dimensionedTensor("omega_dot",dimensionSet(0, 0, -1, 0, 0),tensor::zero)//,
               // zeroGradientFvPatchTensorField::typeName
);

volTensorField epsilon_dot
(
    IOobject
    (
        "epsilon_dot",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
        mesh,
dimensionedTensor("epsilon_dot",dimensionSet(0, 0, -1, 0, 0),tensor::zero)//,
               // zeroGradientFvPatchTensorField::typeName
);

volTensorField omega_tot
(
    IOobject
    (
        "omega_tot",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
        mesh,
dimensionedTensor("omega_tot",dimensionSet(0, 0, 0, 0, 0),tensor::zero),
               zeroGradientFvPatchTensorField::typeName
);

volTensorField epsilon_P_dot
(
    IOobject
    (
        "epsilon_P_dot",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
        mesh,
dimensionedTensor("epsilon_P_dot",dimensionSet(0, 0, -1, 0, 0),tensor::zero)//,
               // zeroGradientFvPatchTensorField::typeName
);

volScalarField epsilon_dot_P_VM
(
    IOobject
    (
        "epsilon_dot_P_VM",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
        mesh,
dimensionedScalar("epsilon_dot_P_VM",dimensionSet(0, 0, -1, 0, 0),0.0),
               zeroGradientFvPatchTensorField::typeName
);

volScalarField Acc_Strain
(
    IOobject
    (
        "Acc_Strain",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
        mesh,
dimensionedScalar("Acc_Strain",dimensionSet(0, 0, 0, 0, 0),0.0),
                zeroGradientFvPatchScalarField::typeName
);


volScalarField TESD
(
    IOobject
    (
        "TESD",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
        mesh,
dimensionedScalar("TESD",dimensionSet(0, 0, -1, 0, 0),0.0),
                zeroGradientFvPatchScalarField::typeName
);

// volTensorField gradU
// (
//     IOobject
//     (
//         "gradU",
//         runTime.timeName(),
//         mesh,
//         IOobject::READ_IF_PRESENT,
//         IOobject::AUTO_WRITE
//     ),
//         mesh,
// dimensionedTensor("gradU",dimensionSet(0, 0, -1, 0, 0),tensor::zero)//,
//                // zeroGradientFvPatchTensorField::typeName
// );


// volTensorField epsilon
// (
//     IOobject
//     (
//         "epsilon",
//         runTime.timeName(),
//         mesh,
//         IOobject::READ_IF_PRESENT,
//         IOobject::AUTO_WRITE
//     ),
//         mesh,
// dimensionedTensor("epsilon",dimensionSet(0, 0, 0, 0, 0),tensor::zero)//,
//                // zeroGradientFvPatchTensorField::typeName
// );

volScalarField sigma_Y
(
    IOobject
    (
        "sigma_Y",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
        mesh,
dimensionedScalar("sigma_Y",dimensionSet(1, -1, -2, 0, 0),sigma_Y_0.value()),
                zeroGradientFvPatchScalarField::typeName
);

volTensorField sigma
(
    IOobject
    (
        "sigma",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
        mesh,
dimensionedTensor("sigma",dimensionSet(1, -1, -2, 0, 0),tensor::zero)//,
                //zeroGradientFvPatchTensorField::typeName
);

volVectorField divsigma
(
    IOobject
    (
        "divsigma",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
        mesh,
dimensionedVector("divsigma",dimensionSet(1, -2, -2, 0, 0),vector::zero)//,
                //zeroGradientFvPatchTensorField::typeName
);


volVectorField Euler_Angles
(
    IOobject
    (
        "Euler_Angles",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
        mesh,
dimensionedVector("Euler_Angles",dimensionSet(0, 0, 0, 0, 0),vector::zero)//,
                //zeroGradientFvPatchTensorField::typeName
);


// volTensorField testSIG
// (
//     IOobject
//     (
//         "testSIG",
//         runTime.timeName(),
//         mesh,
//         IOobject::READ_IF_PRESENT,
//         IOobject::AUTO_WRITE
//     ),
//         mesh,
// dimensionedTensor("testSIG",dimensionSet(1, -1, -2, 0, 0),tensor::zero)//,
//                 //zeroGradientFvPatchTensorField::typeName
// );


volTensorField dsigmadt
(
    IOobject
    (
        "dsigmadt",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
        mesh,
dimensionedTensor("dsigmadt",dimensionSet(1, -1, -3, 0, 0),tensor::zero),
                zeroGradientFvPatchTensorField::typeName
);

// volTensorField dsigmadt_crystal
// (
//     IOobject
//     (
//         "dsigmadt_crystal",
//         runTime.timeName(),
//         mesh,
//         IOobject::READ_IF_PRESENT,
//         IOobject::AUTO_WRITE
//     ),
//         mesh,
// dimensionedTensor("dsigmadt_crystal",dimensionSet(1, -1, -3, 0, 0),tensor::zero),
//                 zeroGradientFvPatchTensorField::typeName
// );

#include "createPhi.H"

/*
label pRefCell = 0;
scalar pRefValue = 0.0;
setRefCell(p, mesh.solution().dict().subDict("PISO"), pRefCell, pRefValue);
mesh.schemes().setFluxRequired(p.name());
*/
