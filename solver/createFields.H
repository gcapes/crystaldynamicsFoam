#include "createRDeltaT.H"

dimensionedScalar pi("pi", dimensionSet(0, 0, 0, 0, 0), constant::mathematical::pi);

Info<< "Reading field p_rgh\n" << endl;
volScalarField p_rgh
(
    IOobject
    (
        "p_rgh",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

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


// volScalarField damping
// (
//     IOobject
//     (
//         "damping",
//         runTime.timeName(),
//         mesh,
//         IOobject::READ_IF_PRESENT,
//         IOobject::AUTO_WRITE
//     ),
//     mesh,
//     dimensionedScalar("damping",dimensionSet(1,-3,-1,0,0),1e12)
// );




const label N_Ori=4;




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

phi1[1]=1.0*pi.value();
phi2[1]=0.0;
phi3[1]=0.0;

phi1[2]=2.0*pi.value();
phi2[2]=0.0;
phi3[2]=0.0;



tensor rot[N_Ori]={tensor::zero};

for (label i = 0; i < N_Ori; i++)
{


 tensor temprot(
    Foam::cos(phi2[i])*Foam::cos(phi3[i]), 
    (Foam::sin(phi1[i])*Foam::sin(phi2[i])*Foam::cos(phi3[i]))-(Foam::cos(phi1[i])*Foam::sin(phi3[i])),//-Foam::cos(phi3)*Foam::sin(phi2), 
    (Foam::cos(phi1[i])*Foam::sin(phi2[i])*Foam::cos(phi3[i]))+(Foam::sin(phi1[i])*Foam::sin(phi3[i])),//Foam::sin(phi2)*Foam::sin(phi3),
   
    Foam::cos(phi2[i])*Foam::sin(phi3[i]), 
    (Foam::sin(phi1[i])*Foam::sin(phi2[i])*Foam::sin(phi3[i]))+(Foam::cos(phi1[i])*Foam::cos(phi3[i])), 
    (Foam::cos(phi1[i])*Foam::sin(phi2[i])*Foam::sin(phi3[i]))-(Foam::sin(phi1[i])*Foam::cos(phi3[i])),
    
    -Foam::sin(phi2[i]), 
    Foam::sin(phi1[i])*Foam::cos(phi2[i]), 
    Foam::cos(phi1[i])*Foam::cos(phi2[i])
);//https://en.wikipedia.org/wiki/Rotation_matrix - 
//"represents an extrinsic rotation whose (improper) Euler angles are α, β, γ, about axes x, y, z" - phi1 ==alpha, phi2==beta,phi3==gamma

rot[i]=temprot;
}


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




Fastor::Tensor<double,3,3,3,3> Cprime2[4];//array of tensor objects



for (label h = 0; h < N_Ori; h++)
{

Fastor::Tensor<double,3,3,3,3> Cprime_temp(0.0);

  Fastor::Tensor<double,3,3> Rot_matrix(0.0);

  Rot_matrix(0,0)=rot[h].component(0);
  Rot_matrix(0,1)=rot[h].component(1);
  Rot_matrix(0,2)=rot[h].component(2);
  Rot_matrix(1,0)=rot[h].component(3);
  Rot_matrix(1,1)=rot[h].component(4);
  Rot_matrix(1,2)=rot[h].component(5);
  Rot_matrix(2,0)=rot[h].component(6);
  Rot_matrix(2,1)=rot[h].component(7);
  Rot_matrix(2,2)=rot[h].component(8);


  // Rotate To Global needReference
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

              // value*=PopBal[h][cell];

        Cprime_temp(i,j,k,l)=value;
      
      }
    }
  }
}

// if(i==0){
//     Cprime_temp(0,0,0,0)=2.35;
// }
// if(i==1){
//     Cprime_temp(0,0,0,0)=247.65;
// }


Cprime2[h]=Cprime_temp;
}

//Cprime2[N_Ori] now has all the grain elasticity tensors


// Fastor::print("Test1",Cprime2[0](0,0,0,0));
// Fastor::print("Test1",Cprime2[1](0,0,0,0));
// Fastor::print("Test1",Cprime2[0](0,0,0,1));
// Fastor::print("Test1",Cprime2[2](0,0,0,0));

// exit(0);

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

volTensorField s
(
    IOobject
    (
        "s",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
        mesh,
dimensionedTensor("s",dimensionSet(1, -1, -2, 0, 0),tensor::zero)//,
                //zeroGradientFvPatchTensorField::typeName
);

volVectorField divs
(
    IOobject
    (
        "divs",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
        mesh,
dimensionedVector("divs",dimensionSet(1, -2, -2, 0, 0),vector::zero)//,
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


volTensorField dsdt
(
    IOobject
    (
        "dsdt",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
        mesh,
dimensionedTensor("dsdt",dimensionSet(1, -1, -3, 0, 0),tensor::zero),
                zeroGradientFvPatchTensorField::typeName
);

volTensorField dsdt_crystal
(
    IOobject
    (
        "dsdt_crystal",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
        mesh,
dimensionedTensor("dsdt_crystal",dimensionSet(1, -1, -3, 0, 0),tensor::zero),
                zeroGradientFvPatchTensorField::typeName
);

#include "createPhi.H"


Info<< "Reading phaseProperties\n" << endl;
immiscibleIncompressibleTwoPhaseMixture mixture(U, phi);

autoPtr<twoPhaseChangeModel> phaseChangePtr
(
    twoPhaseChangeModel::New(mixture)
);


// Need to store rho for ddt(rho, U)
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
    mesh
);
// rho.oldTime();

// rho=mixture.rho();

// #include "readGravitationalAcceleration.H"
// #include "readhRef.H"
// #include "gh.H"


// volScalarField p
// (
//     IOobject
//     (
//         "p",
//         runTime.timeName(),
//         mesh,
//         IOobject::NO_READ,
//         IOobject::AUTO_WRITE
//     ),
//     p_rgh + rho*gh
// );

// pressureReference pressureReference(p, p_rgh, pimple.dict());

// if (p_rgh.needReference())
// {
//     p += dimensionedScalar
//     (
//         "p",
//         p.dimensions(),
//         pressureReference.refValue()
//       - getRefCellValue(p, pressureReference.refCell())
//     );
//     p_rgh = p - rho*gh;
// }

// mesh.schemes().setFluxRequired(p_rgh.name());

// volScalarField& alpha1(mixture.alpha1());
// mesh.schemes().setFluxRequired(alpha1.name());


// #include "createAlphaFluxes.H"

// // Construct incompressible turbulence model
// incompressibleInterPhaseTransportModel turbulence
// (
//     U,
//     phi,
//     alphaPhi1,
//     mixture
// );


#include "createMRF.H"
#include "createFvModels.H"
#include "createFvConstraints.H"
