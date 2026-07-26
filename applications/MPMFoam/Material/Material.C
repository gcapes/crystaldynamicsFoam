/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

\*---------------------------------------------------------------------------*/

#include "Material.H"
#include "Time.H"
#include "fvCFD.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::Material::Material()
:
    Mat_type(""),
    rho_0(0.0),
    Cv(0.0),
    kappa(0.0),
    elasticType("isotropic"),
    lambda(0.0),
    mu(0.0),
    E(0.0),
    nu(0.0),
    C11(0.0),
    C12(0.0),
    C44(0.0),
    C33(0.0),
    C13(0.0),
    C22(0.0),
    C23(0.0),
    C55(0.0),
    C66(0.0),
    C_tensor(),
    C_voigt(),
    fullTensorInitialized(false),
    crystalStructure("FCC"),
    N_Slips(0),
    initial_CRSS(0.0),
    tau_S(0.0),
    h_0(0.0),
    r(0.0),
    k_exponent(1.0),
    ref_shear_rate(1.0e-3)
{
    // Initialize with zeros
    C_tensor.zeros();
    C_voigt.zeros();
}

Foam::Material::Material(const Material& M1)
:
    Mat_type(M1.Mat_type),
    rho_0(M1.rho_0),
    Cv(M1.Cv),
    kappa(M1.kappa),
    elasticType(M1.elasticType),
    lambda(M1.lambda),
    mu(M1.mu),
    E(M1.E),
    nu(M1.nu),
    C11(M1.C11),
    C12(M1.C12),
    C44(M1.C44),
    C33(M1.C33),
    C13(M1.C13),
    C22(M1.C22),
    C23(M1.C23),
    C55(M1.C55),
    C66(M1.C66),
    C_tensor(M1.C_tensor),
    C_voigt(M1.C_voigt),
    fullTensorInitialized(M1.fullTensorInitialized),
    crystalStructure(M1.crystalStructure),
    N_Slips(M1.N_Slips),
    SP_norm(M1.SP_norm),
    SP_direction(M1.SP_direction),
    initial_CRSS(M1.initial_CRSS),
    tau_S(M1.tau_S),
    h_0(M1.h_0),
    r(M1.r),
    k_exponent(M1.k_exponent),
    ref_shear_rate(M1.ref_shear_rate),
    h_matrix(M1.h_matrix)
{
}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::Material::Print(label length)
{
    Info<< "Material Properties:" << nl
        << "  Type: " << Mat_type << nl
        << "  Density: " << rho_0 << " kg/m^3" << nl
        << "  Elastic type: " << elasticType << nl;
    
    if (elasticType == "isotropic")
    {
        Info<< "  Lambda: " << lambda << " Pa" << nl
            << "  Mu: " << mu << " Pa" << nl;
    }
    else if (elasticType == "cubic")
    {
        Info<< "  C11: " << C11 << " Pa" << nl
            << "  C12: " << C12 << " Pa" << nl
            << "  C44: " << C44 << " Pa" << nl;
    }
    else if (elasticType == "hexagonal")
    {
        Info<< "  C11: " << C11 << " Pa" << nl
            << "  C33: " << C33 << " Pa" << nl
            << "  C12: " << C12 << " Pa" << nl
            << "  C13: " << C13 << " Pa" << nl
            << "  C44: " << C44 << " Pa" << nl;
    }
    else if (elasticType == "orthorhombic")
    {
        Info<< "  C11: " << C11 << " Pa" << nl
            << "  C22: " << C22 << " Pa" << nl
            << "  C33: " << C33 << " Pa" << nl
            << "  C12: " << C12 << " Pa" << nl
            << "  C13: " << C13 << " Pa" << nl
            << "  C23: " << C23 << " Pa" << nl
            << "  C44: " << C44 << " Pa" << nl
            << "  C55: " << C55 << " Pa" << nl
            << "  C66: " << C66 << " Pa" << nl;
    }
    else if (elasticType == "anisotropic")
    {
        Info<< "  Using general anisotropic stiffness tensor" << nl;
    }
    
    if (N_Slips > 0)
    {
        Info<< "  Number of slip systems: " << N_Slips << nl
            << "  Initial CRSS: " << initial_CRSS << " Pa" << nl
            << "  Reference shear rate: " << ref_shear_rate << " 1/s" << nl
            << "  Rate sensitivity: " << k_exponent << nl;
    }
}

void Foam::Material::computeLameParameters()
{
    // Compute Lame parameters from Young's modulus and Poisson's ratio
    lambda = E * nu / ((1.0 + nu) * (1.0 - 2.0*nu));
    mu = E / (2.0 * (1.0 + nu));
    
    Info<< "Computed Lame parameters from E and nu:" << nl
        << "  lambda = " << lambda << " Pa" << nl
        << "  mu = " << mu << " Pa" << nl;
}

void Foam::Material::initializeElasticity(const dictionary& dict)
{
    // Read elastic type
    elasticType = dict.lookupOrDefault<word>("elasticType", "isotropic");
    
    if (elasticType == "isotropic")
    {
        // Try to read lambda and mu directly
        if (dict.found("lambda") && dict.found("mu"))
        {
            lambda = readScalar(dict.lookup("lambda"));
            mu = readScalar(dict.lookup("mu"));
        }
        // Or read E and nu and compute lambda and mu
        else if (dict.found("E") && dict.found("nu"))
        {
            E = readScalar(dict.lookup("E"));
            nu = readScalar(dict.lookup("nu"));
            computeLameParameters();
        }
        else
        {
            FatalErrorInFunction
                << "For isotropic elasticity, must provide either "
                << "(lambda, mu) or (E, nu)"
                << exit(FatalError);
        }
        
        buildStiffnessTensor();
    }
    else if (elasticType == "cubic")
    {
        // Read cubic crystal elastic constants
        C11 = readScalar(dict.lookup("C11"));
        C12 = readScalar(dict.lookup("C12"));
        C44 = readScalar(dict.lookup("C44"));
        
        Info<< "Cubic crystal elastic constants:" << nl
            << "  C11 = " << C11 << " Pa" << nl
            << "  C12 = " << C12 << " Pa" << nl
            << "  C44 = " << C44 << " Pa" << nl;
        
        buildStiffnessTensor();
    }
    else if (elasticType == "hexagonal")
    {
        // Read hexagonal crystal elastic constants (5 independent)
        C11 = readScalar(dict.lookup("C11"));
        C33 = readScalar(dict.lookup("C33"));
        C12 = readScalar(dict.lookup("C12"));
        C13 = readScalar(dict.lookup("C13"));
        C44 = readScalar(dict.lookup("C44"));
        
        // C66 is derived: C66 = (C11 - C12)/2
        C66 = (C11 - C12) / 2.0;
        
        Info<< "Hexagonal crystal elastic constants:" << nl
            << "  C11 = " << C11 << " Pa" << nl
            << "  C33 = " << C33 << " Pa" << nl
            << "  C12 = " << C12 << " Pa" << nl
            << "  C13 = " << C13 << " Pa" << nl
            << "  C44 = " << C44 << " Pa" << nl
            << "  C66 = " << C66 << " Pa (derived)" << nl;
        
        buildStiffnessTensor();
    }
    else if (elasticType == "orthorhombic")
    {
        // Read orthorhombic elastic constants (9 independent)
        C11 = readScalar(dict.lookup("C11"));
        C22 = readScalar(dict.lookup("C22"));
        C33 = readScalar(dict.lookup("C33"));
        C12 = readScalar(dict.lookup("C12"));
        C13 = readScalar(dict.lookup("C13"));
        C23 = readScalar(dict.lookup("C23"));
        C44 = readScalar(dict.lookup("C44"));
        C55 = readScalar(dict.lookup("C55"));
        C66 = readScalar(dict.lookup("C66"));
        
        Info<< "Orthorhombic crystal elastic constants:" << nl
            << "  C11 = " << C11 << " Pa" << nl
            << "  C22 = " << C22 << " Pa" << nl
            << "  C33 = " << C33 << " Pa" << nl
            << "  C12 = " << C12 << " Pa" << nl
            << "  C13 = " << C13 << " Pa" << nl
            << "  C23 = " << C23 << " Pa" << nl
            << "  C44 = " << C44 << " Pa" << nl
            << "  C55 = " << C55 << " Pa" << nl
            << "  C66 = " << C66 << " Pa" << nl;
        
        buildStiffnessTensor();
    }
    else if (elasticType == "anisotropic" || elasticType == "general")
    {
        // Read general anisotropic stiffness tensor
        readGeneralStiffness(dict);
        voigtToTensor();
        fullTensorInitialized = true;
    }
    else
    {
        FatalErrorInFunction
            << "Unknown elasticType: " << elasticType << nl
            << "Valid options: isotropic, cubic, hexagonal, orthorhombic, anisotropic"
            << exit(FatalError);
    }
}

void Foam::Material::readGeneralStiffness(const dictionary& dict)
{
    /*
    Read general stiffness components in the format:
    C0000 value;
    C1111 value;
    etc.
    
    Uses Voigt notation mapping:
    0 -> 11, 1 -> 22, 2 -> 33, 3 -> 23, 4 -> 13, 5 -> 12
    */
    
    // Initialize Voigt matrix
    C_voigt.zeros();
    
    // Define mapping from 4-index to Voigt notation
    // Voigt: 0=11, 1=22, 2=33, 3=23, 4=13, 5=12
    auto voigtIndex = [](label i, label j) -> label
    {
        if (i == j) return i;  // 00->0, 11->1, 22->2
        if ((i == 1 && j == 2) || (i == 2 && j == 1)) return 3;  // 12 or 21
        if ((i == 0 && j == 2) || (i == 2 && j == 0)) return 4;  // 02 or 20
        if ((i == 0 && j == 1) || (i == 1 && j == 0)) return 5;  // 01 or 10
        return 0;
    };
    
    // Read all possible components
    // The user can specify using notation like C0000, C1111, C0011, etc.
    
    Info<< "Reading general anisotropic stiffness tensor..." << nl;
    
    for (label i = 0; i < 3; i++)
    {
        for (label j = 0; j < 3; j++)
        {
            for (label k = 0; k < 3; k++)
            {
                for (label l = 0; l < 3; l++)
                {
                    word componentName = "C" 
                        + Foam::name(i) + Foam::name(j) 
                        + Foam::name(k) + Foam::name(l);
                    
                    if (dict.found(componentName))
                    {
                        scalar value = readScalar(dict.lookup(componentName));
                        
                        // Map to Voigt notation
                        label voigtI = voigtIndex(i, j);
                        label voigtJ = voigtIndex(k, l);
                        
                        C_voigt(voigtI, voigtJ) = value;
                        
                        Info<< "  " << componentName << " = " << value << " Pa" << nl;
                    }
                }
            }
        }
    }
    
    // Enforce symmetry
    for (label i = 0; i < 6; i++)
    {
        for (label j = i+1; j < 6; j++)
        {
            if (C_voigt(i,j) != 0.0 || C_voigt(j,i) != 0.0)
            {
                scalar avg = (C_voigt(i,j) + C_voigt(j,i)) / 2.0;
                C_voigt(i,j) = avg;
                C_voigt(j,i) = avg;
            }
        }
    }
    
    Info<< "Stiffness tensor in Voigt notation (6x6):" << nl;
    for (label i = 0; i < 6; i++)
    {
        Info<< "  ";
        for (label j = 0; j < 6; j++)
        {
            Info<< C_voigt(i,j) << " ";
        }
        Info<< nl;
    }
}

void Foam::Material::buildStiffnessTensor() const
{
    /*
    Build the 4th-order stiffness tensor C_ijkl based on the elastic type
    */
    
    C_tensor.zeros();
    C_voigt.zeros();
    
    if (elasticType == "isotropic")
    {
        // Isotropic: C_ijkl = lambda * delta_ij * delta_kl + mu * (delta_ik * delta_jl + delta_il * delta_jk)
        for (label i = 0; i < 3; i++)
        {
            for (label j = 0; j < 3; j++)
            {
                for (label k = 0; k < 3; k++)
                {
                    for (label l = 0; l < 3; l++)
                    {
                        scalar delta_ij = (i == j) ? 1.0 : 0.0;
                        scalar delta_kl = (k == l) ? 1.0 : 0.0;
                        scalar delta_ik = (i == k) ? 1.0 : 0.0;
                        scalar delta_jl = (j == l) ? 1.0 : 0.0;
                        scalar delta_il = (i == l) ? 1.0 : 0.0;
                        scalar delta_jk = (j == k) ? 1.0 : 0.0;
                        
                        C_tensor(i,j,k,l) = lambda * delta_ij * delta_kl 
                                          + mu * (delta_ik * delta_jl + delta_il * delta_jk);
                    }
                }
            }
        }
    }
    else if (elasticType == "cubic")
    {
        // Cubic symmetry: only C11, C12, C44 are independent
        // C_iiii = C11
        // C_iijj = C12 (i != j)
        // C_ijij = C44 (i != j)
        
        for (label i = 0; i < 3; i++)
        {
            C_tensor(i,i,i,i) = C11;
            
            for (label j = 0; j < 3; j++)
            {
                if (i != j)
                {
                    C_tensor(i,i,j,j) = C12;
                    C_tensor(i,j,i,j) = C44;
                    C_tensor(i,j,j,i) = C44;
                }
            }
        }
    }
    else if (elasticType == "hexagonal")
    {
        // Hexagonal symmetry (transverse isotropy in x-y plane, z is unique axis)
        // 5 independent constants: C11, C33, C12, C13, C44
        // C66 = (C11 - C12)/2
        
        // Diagonal terms
        C_tensor(0,0,0,0) = C11;
        C_tensor(1,1,1,1) = C11;
        C_tensor(2,2,2,2) = C33;
        
        // Off-diagonal coupling
        C_tensor(0,0,1,1) = C12;
        C_tensor(1,1,0,0) = C12;
        C_tensor(0,0,2,2) = C13;
        C_tensor(2,2,0,0) = C13;
        C_tensor(1,1,2,2) = C13;
        C_tensor(2,2,1,1) = C13;
        
        // Shear terms
        C_tensor(0,1,0,1) = C66;
        C_tensor(0,1,1,0) = C66;
        C_tensor(1,0,0,1) = C66;
        C_tensor(1,0,1,0) = C66;
        
        C_tensor(0,2,0,2) = C44;
        C_tensor(0,2,2,0) = C44;
        C_tensor(2,0,0,2) = C44;
        C_tensor(2,0,2,0) = C44;
        
        C_tensor(1,2,1,2) = C44;
        C_tensor(1,2,2,1) = C44;
        C_tensor(2,1,1,2) = C44;
        C_tensor(2,1,2,1) = C44;
    }
    else if (elasticType == "orthorhombic")
    {
        // Orthorhombic symmetry: 9 independent constants
        
        // Diagonal terms
        C_tensor(0,0,0,0) = C11;
        C_tensor(1,1,1,1) = C22;
        C_tensor(2,2,2,2) = C33;
        
        // Off-diagonal coupling
        C_tensor(0,0,1,1) = C12;
        C_tensor(1,1,0,0) = C12;
        C_tensor(0,0,2,2) = C13;
        C_tensor(2,2,0,0) = C13;
        C_tensor(1,1,2,2) = C23;
        C_tensor(2,2,1,1) = C23;
        
        // Shear terms
        C_tensor(1,2,1,2) = C44;
        C_tensor(1,2,2,1) = C44;
        C_tensor(2,1,1,2) = C44;
        C_tensor(2,1,2,1) = C44;
        
        C_tensor(0,2,0,2) = C55;
        C_tensor(0,2,2,0) = C55;
        C_tensor(2,0,0,2) = C55;
        C_tensor(2,0,2,0) = C55;
        
        C_tensor(0,1,0,1) = C66;
        C_tensor(0,1,1,0) = C66;
        C_tensor(1,0,0,1) = C66;
        C_tensor(1,0,1,0) = C66;
    }
    
    // Convert to Voigt notation for convenience
    tensorToVoigt();
    fullTensorInitialized = true;
}

void Foam::Material::voigtToTensor() const
{
    /*
    Convert from Voigt 6x6 matrix to full 4th-order tensor
    Voigt notation: 0=11, 1=22, 2=33, 3=23, 4=13, 5=12
    */
    
    C_tensor.zeros();
    
    // Mapping from Voigt indices to tensor indices
    label voigtToI[6] = {0, 1, 2, 1, 0, 0};
    label voigtToJ[6] = {0, 1, 2, 2, 2, 1};
    
    for (label v1 = 0; v1 < 6; v1++)
    {
        label i = voigtToI[v1];
        label j = voigtToJ[v1];
        
        for (label v2 = 0; v2 < 6; v2++)
        {
            label k = voigtToI[v2];
            label l = voigtToJ[v2];
            
            C_tensor(i,j,k,l) = C_voigt(v1, v2);
            
            // Add symmetries
            if (i != j)
            {
                C_tensor(j,i,k,l) = C_voigt(v1, v2);
            }
            if (k != l)
            {
                C_tensor(i,j,l,k) = C_voigt(v1, v2);
            }
            if (i != j && k != l)
            {
                C_tensor(j,i,l,k) = C_voigt(v1, v2);
            }
        }
    }
}

void Foam::Material::tensorToVoigt() const
{
    /*
    Convert from 4th-order tensor to Voigt 6x6 matrix
    */
    
    C_voigt.zeros();
    
    // Mapping
    label voigtToI[6] = {0, 1, 2, 1, 0, 0};
    label voigtToJ[6] = {0, 1, 2, 2, 2, 1};
    
    for (label v1 = 0; v1 < 6; v1++)
    {
        label i = voigtToI[v1];
        label j = voigtToJ[v1];
        
        for (label v2 = 0; v2 < 6; v2++)
        {
            label k = voigtToI[v2];
            label l = voigtToJ[v2];
            
            C_voigt(v1, v2) = C_tensor(i,j,k,l);
        }
    }
}

const Fastor::Tensor<double, 3, 3, 3, 3>& Foam::Material::getStiffnessTensor() const
{
    if (!fullTensorInitialized)
    {
        buildStiffnessTensor();
    }
    return C_tensor;
}

const Fastor::Tensor<double, 6, 6>& Foam::Material::getVoigtMatrix() const
{
    if (!fullTensorInitialized)
    {
        buildStiffnessTensor();
    }
    return C_voigt;
}

void Foam::Material::initializePlasticity(const dictionary& dict)
{
    // Read number of slip systems
    N_Slips = readInt(dict.lookup("N_Slips"));

    crystalStructure = dict.lookupOrDefault<word>("crystalStructure", "FCC");
    
    // Read flow rule parameters
    initial_CRSS = readScalar(dict.lookup("initial_CRSS"));
    tau_S = dict.lookupOrDefault<scalar>("tau_S", initial_CRSS * 2.0);
    h_0 = dict.lookupOrDefault<scalar>("h_0", initial_CRSS * 10.0);
    r = dict.lookupOrDefault<scalar>("r", 1.0);
    
    // Read slip rate parameters
    ref_shear_rate = readScalar(dict.lookup("ref_shear_rate"));
    k_exponent = readScalar(dict.lookup("k_exponent"));
    
    Info<< "Crystal plasticity parameters:" << nl
        << "  N_Slips = " << N_Slips << nl
        << "  initial_CRSS = " << initial_CRSS << " Pa" << nl
        << "  ref_shear_rate = " << ref_shear_rate << " 1/s" << nl
        << "  k_exponent = " << k_exponent << nl;
    
    // Initialize hardening matrix
    h_matrix.resize(N_Slips);
    for (int i = 0; i < N_Slips; i++)
    {
        h_matrix[i].resize(N_Slips);
        for (int j = 0; j < N_Slips; j++)
        {
            h_matrix[i][j] = (i == j) ? 1.0 : r;
        }
    }
}

// ************************************************************************* //