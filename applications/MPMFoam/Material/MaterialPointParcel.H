/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | 
    \\  /    A nd           | 
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

Class
    Foam::MaterialPointParcel

Description
    Material point particle for MPM with crystal plasticity capability

\*---------------------------------------------------------------------------*/

#ifndef MaterialPointParcel_H
#define MaterialPointParcel_H

#include "particle.H"
#include "Cloud.H"
#include "IOstream.H"
#include "volFields.H"
#include "pointFields.H"
#include "fvMesh.H"
#include "polyMesh.H"
#include "labelList.H"
#include "cellList.H"
#include "faceList.H"
#include "HashSet.H"

namespace Foam
{

// Forward declaration
class MaterialPointParcel;

// Define cloud type
typedef Cloud<MaterialPointParcel> materialPointCloud;

class MaterialPointParcel
:
    public particle
{
public:
    // Crystal plasticity slip system data
    struct SlipSystem
    {
        vector m;  // Slip direction in reference config
        vector n;  // Slip plane normal in reference config
        tensor schmidTensor; // m ⊗ n
        
        SlipSystem() : m(Zero), n(Zero), schmidTensor(Zero) {}
    };

private:
    // Private member data
    
    // Basic MPM properties
    scalar mass_;
    scalar volume0_;  // Reference volume
    scalar volume_;   // Current volume
    vector velocity_;
    tensor stress_;   // Cauchy stress
    
    // Deformation kinematics
    tensor F_;        // Total deformation gradient
    tensor Fp_;       // Plastic deformation gradient
    tensor Fe_;       // Elastic deformation gradient: F = Fe*Fp
    tensor Cp_;       // Plastic right Cauchy-Green tensor: Fp^T*Fp
    
    // Crystal orientation
    label materialID_;
    label grainID_;
    scalar phi1_, phi2_, phi3_;  // Bunge Euler angles
    tensor Q0_;       // Initial rotation matrix from Euler angles
    
    // Crystal plasticity state variables
    List<scalar> gamma_;      // Slip on each system
    List<scalar> gamma_dot_;  // Slip rates
    List<scalar> tau_;        // Resolved shear stress on each system
    List<scalar> g_;          // Current hardness on each system
    List<SlipSystem> slipSystems_;  // Slip system geometry
    
    // Internal variables
    scalar epsilonP_eq_;      // Equivalent plastic strain
    scalar workPlastic_;      // Plastic work
    scalar temperature_;      // Temperature for thermal effects

public:
    
    // Constructors
    MaterialPointParcel
    (
        const polyMesh& mesh,
        const vector& position,
        const label celli
    );
    
    MaterialPointParcel
    (
        const polyMesh& mesh,
        Istream& is,
        bool readFields = true,
        bool newFormat = true
    );
    
    // Access functions
    inline scalar mass() const { return mass_; }
    inline scalar& mass() { return mass_; }
    
    inline scalar volume() const { return volume_; }
    inline scalar& volume() { return volume_; }

    inline scalar volume0() const { return volume0_; }
    inline scalar& volume0() { return volume0_; }
    
    inline const vector& velocity() const { return velocity_; }
    inline vector& velocity() { return velocity_; }
    
    inline const tensor& stress() const { return stress_; }
    inline tensor& stress() { return stress_; }
    
    inline const tensor& F() const { return F_; }
    inline tensor& F() { return F_; }
    
    inline const tensor& Fp() const { return Fp_; }
    inline const tensor& Fe() const { return Fe_; }
    
    inline label materialID() const { return materialID_; }
    inline label grainID() const { return grainID_; }
    
    // MPM specific functions
    
    // Shape function calculations
    void shapeFunctions
    (
        const fvMesh& mesh,
        scalarField& N,
        List<vector>& gradN,
        labelList& nodeIDs
    ) const;
    
    // GIMP shape functions (Generalized Interpolation Material Point)
    void shapeFunctionsGIMP
    (
        const fvMesh& mesh,
        const scalar& lp,  // Particle characteristic length
        scalarField& N,
        List<vector>& gradN,
        labelList& nodeIDs
    ) const;
    
    // Particle to grid transfer (P2G)
    void transferToGrid
    (
        const fvMesh& mesh,
        pointScalarField& nodalMass,
        pointVectorField& nodalMomentum,
        pointTensorField& nodalStress,
        pointVectorField& nodalForce
    ) const;
    
    // Grid to particle transfer (G2P)
    void transferFromGrid
    (
        const fvMesh& mesh,
        const pointVectorField& nodalVelocity,
        const pointVectorField& nodalAcceleration
    );
    
    // Update deformation gradient
    void updateDeformationGradient
    (
        const fvMesh& mesh,
        const pointVectorField& nodalVelocity,
        const scalar dt
    );
    
    // Crystal plasticity functions
    
    // Initialize slip systems for given crystal structure
    void initializeSlipSystems(const word& crystalType);
    
    // Initialize crystal orientation from Euler angles
    void setCrystalOrientation(scalar phi1, scalar phi2, scalar phi3);
    
    // Compute resolved shear stress on all slip systems
    void computeResolvedShearStress();
    
    // Update plastic deformation (return mapping algorithm)
    void updatePlasticDeformation
    (
        const scalar dt,
        const dictionary& materialProps
    );
    
    // Compute elastic stress (hyperelastic or hypoelastic)
    void computeElasticStress(const dictionary& materialProps);
    
    // Update hardness evolution
    void updateHardening
    (
        const scalar dt,
        const dictionary& materialProps
    );
    
    // Constitutive update (complete stress update)
    void constitutiveUpdate
    (
        const scalar dt,
        const dictionary& materialProps
    );
    
    // Tracking
    typedef particle::trackingData trackingData;
    
    bool move
    (
        materialPointCloud& cloud,
        trackingData& td,
        const scalar dt
    );
    
    // IO
    static void readFields(materialPointCloud& c);
    static void writeFields(const materialPointCloud& c);
    
    // IOstream operators
    friend Ostream& operator<<(Ostream&, const MaterialPointParcel&);
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

#endif

// ************************************************************************* //