/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

Class
    Foam::Material

Description
    Material properties class with general anisotropic elasticity support

\*---------------------------------------------------------------------------*/

#ifndef material_H
#define material_H

#include <vector>
#include "Polynomial.H"
#include "IOdictionary.H"
#include "PtrDictionary.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "Fastor/Fastor.h"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class Material
{
private:
    
public:
    string Mat_type;
    
    // Basic thermophysical properties
    scalar rho_0;   // Density
    scalar Cv;      // Specific heat capacity
    scalar kappa;   // Thermal conductivity
    
    // ========================================================================
    // ELASTICITY PROPERTIES
    // ========================================================================
    
    // Elastic type: "isotropic", "cubic", "hexagonal", "orthorhombic", "anisotropic"
    word elasticType;
    
    // For isotropic elasticity
    scalar lambda;    // First Lame parameter
    scalar mu;        // Shear modulus (second Lame parameter)
    scalar E;         // Young's modulus (alternative to lambda/mu)
    scalar nu;        // Poisson's ratio (alternative to lambda/mu)
    
    // For cubic crystal elasticity (3 independent constants)
    scalar C11;       // C1111
    scalar C12;       // C1122
    scalar C44;       // C1212
    
    // For hexagonal crystal elasticity (5 independent constants)
    scalar C33;       // C3333
    scalar C13;       // C1133
    // C11, C12, C44 reused
    
    // For orthorhombic crystal elasticity (9 independent constants)
    scalar C22;       // C2222
    scalar C23;       // C2233
    scalar C55;       // C1313
    scalar C66;       // C2323
    // C11, C12, C13, C33 reused
    
    // Full 4th-order stiffness tensor (for arbitrary anisotropy)
    // This can store any elastic symmetry
    mutable Fastor::Tensor<double, 3, 3, 3, 3> C_tensor;
    
    // Voigt notation (6x6 matrix representation) - useful for some calculations
    mutable Fastor::Tensor<double, 6, 6> C_voigt;
    
    // Flag to indicate if full tensor has been initialized
    mutable bool fullTensorInitialized;
    
    // ========================================================================
    // CRYSTAL PLASTICITY PROPERTIES
    // ========================================================================
    
    // Crystal structure type for slip system initialization
    word crystalStructure;    // "FCC", "BCC", "HCP", etc.
    
    int N_Slips;
    
    std::vector<vector> SP_norm;       // Slip plane normals
    std::vector<vector> SP_direction;  // Slip directions
    
    double initial_CRSS;    // Initial critical resolved shear stress
    double tau_S;           // Saturation shear stress
    double h_0;             // Initial hardening modulus
    double r;               // Latent hardening parameter
    double k_exponent;      // Rate sensitivity exponent
    double ref_shear_rate;  // Reference shear rate
    
    // Hardening matrix h_matrix(alpha, beta): hardening of system alpha due to slip on beta
    std::vector<std::vector<scalar>> h_matrix;
    
    // Constructors
    Material();
    Material(const Material&);
    virtual ~Material() {}
    
    // Member Functions
    void Print(label length);
    
    // Initialize elastic properties from dictionary
    void initializeElasticity(const dictionary& dict);
    
    // Initialize crystal plasticity from dictionary
    void initializePlasticity(const dictionary& dict);
    
    // Compute lambda and mu from E and nu
    void computeLameParameters();
    
    // Build full 4th-order stiffness tensor based on elastic type
    void buildStiffnessTensor() const;
    
    // Read general stiffness components from dictionary
    void readGeneralStiffness(const dictionary& dict);
    
    // Convert Voigt notation to 4th-order tensor
    void voigtToTensor() const;
    
    // Convert 4th-order tensor to Voigt notation
    void tensorToVoigt() const;
    
    // Get stiffness tensor (builds it if not already done)
    const Fastor::Tensor<double, 3, 3, 3, 3>& getStiffnessTensor() const;
    
    // Get Voigt matrix (builds it if not already done)
    const Fastor::Tensor<double, 6, 6>& getVoigtMatrix() const;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //