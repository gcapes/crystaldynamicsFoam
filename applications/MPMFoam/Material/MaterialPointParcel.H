/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | 
    \\  /    A nd           | 
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

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
#include "Fastor/Fastor.h"
#include "Material.H"

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

    // Activity flag
    bool active_;  // True = active particle, False = inactive (don't process)
    
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

    // Private helper functions
    void polarDecomposition
    (
        const tensor& F,
        tensor& R,
        tensor& U
    ) const;
    
    void eigenDecomposition
    (
        const tensor& T,
        vector& eigenvalues,
        tensor& eigenvectors
    ) const;

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

    // ========================================================================
    // ACCESS FUNCTIONS
    // ========================================================================
    
    inline bool active() const { return active_; }
    inline bool& active() { return active_; }

    // Access to base class position and cell (for recycling)
inline void setPosition(const vector& pos)
{
    // Access the protected position_ member from particle base class
    static_cast<particle&>(*this).position() = pos;
}

inline void setCell(label cellI)
{
    // Access the protected cell_ member from particle base class
    static_cast<particle&>(*this).cell() = cellI;
}
    
    // Deactivate particle (e.g., when leaving domain)
    inline void deactivate() { active_ = false; }
    
    // Reactivate particle (e.g., when re-entering domain)
    inline void activate() { active_ = true; }

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
    inline tensor& Fp() { return Fp_; }
    
    inline const tensor& Fe() const { return Fe_; }
    inline tensor& Fe() { return Fe_; }
    
    inline const tensor& Cp() const { return Cp_; }
    inline tensor& Cp() { return Cp_; }
    
    inline label materialID() const { return materialID_; }
    inline label& materialID() { return materialID_; }
    
    inline label grainID() const { return grainID_; }
    inline label& grainID() { return grainID_; }

    inline scalar phi1() const { return phi1_; }
    inline scalar phi2() const { return phi2_; }
    inline scalar phi3() const { return phi3_; }
    inline scalar& phi1() { return phi1_; }
    inline scalar& phi2() { return phi2_; }
    inline scalar& phi3() { return phi3_; }

    inline const tensor& Q0() const { return Q0_; }
    inline tensor& Q0() { return Q0_; }

    // Crystal plasticity state accessors
    inline const List<scalar>& gamma() const { return gamma_; }
    inline List<scalar>& gamma() { return gamma_; }
    
    inline const List<scalar>& gamma_dot() const { return gamma_dot_; }
    inline List<scalar>& gamma_dot() { return gamma_dot_; }
    
    inline const List<scalar>& tau() const { return tau_; }
    inline List<scalar>& tau() { return tau_; }
    
    inline const List<scalar>& g() const { return g_; }
    inline List<scalar>& g() { return g_; }
    
    inline scalar epsilonP_eq() const { return epsilonP_eq_; }
    inline scalar& epsilonP_eq() { return epsilonP_eq_; }
    
    inline scalar workPlastic() const { return workPlastic_; }
    inline scalar& workPlastic() { return workPlastic_; }
    
    inline scalar temperature() const { return temperature_; }
    inline scalar& temperature() { return temperature_; }

    // ========================================================================
    // MPM SHAPE FUNCTIONS AND TRANSFER
    // ========================================================================
    
    // Standard shape function calculations
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
    
    // Update deformation gradient from grid velocity
    void updateDeformationGradient
    (
        const fvMesh& mesh,
        const pointVectorField& nodalVelocity,
        const scalar dt
    );


   bool hitPatch(materialPointCloud&, trackingData&);
   bool hitWallPatch(materialPointCloud&, trackingData&);

   void prepareForParallelTransfer();

   void correctAfterParallelTransfer
(
    const label patchi,
    trackingData& td
);

//    void transformProperties(const transformer&);

    // ========================================================================
    // CRYSTAL PLASTICITY INITIALIZATION
    // ========================================================================
    
    // Initialize slip systems for given crystal structure
    void initializeSlipSystems(const word& crystalType);
    
    // Initialize crystal orientation from Euler angles
    void setCrystalOrientation(scalar phi1, scalar phi2, scalar phi3);
    
    // Build rotation matrix from Bunge Euler angles
    Fastor::Tensor<double, 3, 3> buildRotationMatrix
    (
        const scalar phi1,
        const scalar phi2,
        const scalar phi3
    ) const;

    // ========================================================================
    // STIFFNESS TENSOR CREATION
    // ========================================================================
    
    // Create 4th-order stiffness tensor for cubic crystals
    Fastor::Tensor<double, 3, 3, 3, 3> createCrystalStiffness
    (
        const Material& material
    ) const;
    
    // Create 4th-order isotropic stiffness tensor
    Fastor::Tensor<double, 3, 3, 3, 3> createIsotropicStiffness
    (
        const scalar lambda,
        const scalar mu
    ) const;
    
    // Rotate stiffness tensor using current crystal orientation
    Fastor::Tensor<double, 3, 3, 3, 3> rotateStiffnessTensor
    (
        const Material& material
    ) const;

    // ========================================================================
    // CONSTITUTIVE UPDATE (UNIFIED FUNCTION)
    // ========================================================================
    
    // Complete constitutive update: compute elastic stress and update plasticity
    void constitutiveUpdate
    (
        const Material& material,
        const tensor& L,        // Velocity gradient
        const scalar& dt
    );

    // ========================================================================
    // ELASTIC STRESS COMPUTATION
    // ========================================================================
    
    // Compute elastic stress from Fe using Fastor
    void computeElasticStress
    (
        const Material& material
    );

    // ========================================================================
    // PLASTIC DEFORMATION UPDATE
    // ========================================================================
    
    // Update plastic deformation using return mapping
    void updatePlasticDeformation
    (
        const Material& material,
        const scalar& dt
    );
    
    // Compute resolved shear stress on all slip systems
    void computeResolvedShearStress();
    
    // Update hardness evolution
    void updateHardening
    (
        const Material& material,
        const scalar& dt
    );
    
    // Update crystal rotation based on lattice spin
    // Update crystal lattice orientation (elastic deformation only)
    void updateCrystalRotation(const scalar dt);
    
    // Update crystal orientation incrementally using elastic spin (for plasticity)
    void updateCrystalRotationIncremental
    (
        const Material& material,
        const tensor& L,
        const scalar& dt
    );

    // ========================================================================
    // TRACKING AND IO
    // ========================================================================
    
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