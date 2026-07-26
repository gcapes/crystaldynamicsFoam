/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | 
    \\  /    A nd           | 
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

Description
    Material Point Parcel with Crystal Plasticity Implementation

\*---------------------------------------------------------------------------*/
#include "Fastor/Fastor.h"
#include "FastorOpenFOAMUtils.H"
#include "Material.H"

#include "MaterialPointParcel.H"
#include "meshSearch.H"
#include "mathematicalConstants.H"
#include "polyMesh.H"
#include "IOField.H"




using namespace Foam::constant::mathematical;
using Fastor::Index;


// Declare einsum indices for Fastor tensor contractions
constexpr size_t i = 0, j = 1, k = 2, l = 3;
constexpr size_t m = 4, n = 5, p = 6, q = 7;

// inline scalar sign(scalar x)
// {
//     return (x > 0.0) ? 1.0 : ((x < 0.0) ? -1.0 : 0.0);
// }



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::MaterialPointParcel::MaterialPointParcel
(
    const polyMesh& mesh,
    const vector& position,
    const label celli
)
:
    particle(mesh, position, celli),
    active_(true),
    mass_(0.0),
    volume0_(0.0),
    volume_(0.0),
    velocity_(Zero),
    stress_(Zero),
    F_(tensor::I),
    Fp_(tensor::I),
    Fe_(tensor::I),
    Cp_(tensor::I),
    materialID_(0),
    grainID_(0),
    phi1_(0.0),
    phi2_(0.0),
    phi3_(0.0),
    Q0_(tensor::I),
    epsilonP_eq_(0.0),
    workPlastic_(0.0),
    temperature_(300.0)
{}

Foam::MaterialPointParcel::MaterialPointParcel
(
    const polyMesh& mesh,
    Istream& is,
    bool readFields,
    bool newFormat
)
:
    particle(mesh, is, readFields, newFormat)
{
    if (readFields)
    {
        // Read all material point data
        is  >> active_ >> mass_ >> volume0_ >> volume_ >> velocity_ 
            >> stress_ >> F_ >> Fp_ >> Fe_ >> Cp_
            >> materialID_ >> grainID_
            >> phi1_ >> phi2_ >> phi3_ >> Q0_
            >> epsilonP_eq_ >> workPlastic_ >> temperature_;
        
        // Read slip system data if present
        label nSlips;
        is >> nSlips;
        
        if (nSlips > 0)
        {
            gamma_.setSize(nSlips);
            gamma_dot_.setSize(nSlips);
            tau_.setSize(nSlips);
            g_.setSize(nSlips);
            
            is >> gamma_ >> gamma_dot_ >> tau_ >> g_;
            
            // Note: slipSystems_ geometry would need to be reconstructed
            // based on crystal type - not transferred
        }
    }
    
    // Check for valid data
    if (mass_ <= 0) mass_ = SMALL;
    if (volume_ <= 0) volume_ = SMALL;
    if (det(F_) <= 0) F_ = tensor::I;
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::MaterialPointParcel::initializeSlipSystems(const word& crystalType)
{
    /*
    Initialize slip systems based on crystal structure
    
    FCC: 12 slip systems {111}<110>
    BCC: 12 slip systems {110}<111> (can extend to 24 or 48)
    HCP: Multiple families (basal, prismatic, pyramidal)
    */
    
    slipSystems_.clear();
    
    if (crystalType == "FCC")
    {
        // FCC: {111}<110> - 12 slip systems
        slipSystems_.setSize(12);
        
        // Four {111} planes
        List<vector> planes(4);
        planes[0] = vector(1, 1, 1);
        planes[1] = vector(-1, 1, 1);
        planes[2] = vector(1, -1, 1);
        planes[3] = vector(1, 1, -1);
        
        // Three <110> directions per plane
        List<List<vector>> directions(4);
        directions[0].setSize(3);
        directions[0][0] = vector(1, -1, 0);
        directions[0][1] = vector(-1, 0, 1);
        directions[0][2] = vector(0, 1, -1);
        
        directions[1].setSize(3);
        directions[1][0] = vector(1, 1, 0);
        directions[1][1] = vector(1, 0, 1);
        directions[1][2] = vector(0, 1, 1);
        
        directions[2].setSize(3);
        directions[2][0] = vector(1, 1, 0);
        directions[2][1] = vector(1, 0, -1);
        directions[2][2] = vector(0, 1, 1);
        
        directions[3].setSize(3);
        directions[3][0] = vector(1, -1, 0);
        directions[3][1] = vector(1, 0, 1);
        directions[3][2] = vector(0, 1, -1);
        
        label sysIndex = 0;
        forAll(planes, i)
        {
            vector n = planes[i] / mag(planes[i]);
            
            forAll(directions[i], j)
            {
                vector m = directions[i][j];
                m = m / mag(m);
                
                // Ensure orthogonality (m · n should be ~0)
                m = m - (m & n) * n;
                m = m / mag(m);
                
                slipSystems_[sysIndex].n = n;
                slipSystems_[sysIndex].m = m;
                slipSystems_[sysIndex].schmidTensor = m * n;
                
                sysIndex++;
            }
        }
    }
    else if (crystalType == "BCC")
    {
        // BCC: {110}<111> - 12 slip systems
        slipSystems_.setSize(12);
        
        List<vector> planes(6);
        planes[0] = vector(1, 1, 0);
        planes[1] = vector(1, -1, 0);
        planes[2] = vector(0, 1, 1);
        planes[3] = vector(0, 1, -1);
        planes[4] = vector(1, 0, 1);
        planes[5] = vector(1, 0, -1);
        
        List<vector> directions(4);
        directions[0] = vector(1, 1, 1);
        directions[1] = vector(-1, 1, 1);
        directions[2] = vector(1, -1, 1);
        directions[3] = vector(1, 1, -1);
        
        label sysIndex = 0;
        forAll(planes, i)
        {
            vector n = planes[i] / mag(planes[i]);
            
            forAll(directions, j)
            {
                vector m = directions[j] / mag(directions[j]);
                
                // Only use directions that lie in the plane
                if (mag(m & n) < 0.1)  // Nearly orthogonal
                {
                    m = m - (m & n) * n;
                    m = m / mag(m);
                    
                    slipSystems_[sysIndex].n = n;
                    slipSystems_[sysIndex].m = m;
                    slipSystems_[sysIndex].schmidTensor = m * n;
                    
                    sysIndex++;
                    
                    if (sysIndex >= 12) break;
                }
            }
            if (sysIndex >= 12) break;
        }
    }
    else if (crystalType == "HCP")
    {
        // HCP: Basal (3) + Prismatic (3) + Pyramidal (6) = 12 systems
        slipSystems_.setSize(12);
        
        scalar c_a = 1.633;  // c/a ratio for ideal HCP
        
        // Basal slip: {0001}<11-20> - 3 systems
        vector basalNormal(0, 0, 1);
        
        slipSystems_[0].n = basalNormal;
        slipSystems_[0].m = vector(1, 0, 0) / mag(vector(1, 0, 0));
        slipSystems_[0].schmidTensor = slipSystems_[0].m * slipSystems_[0].n;
        
        slipSystems_[1].n = basalNormal;
        slipSystems_[1].m = vector(-0.5, 0.866, 0) / mag(vector(-0.5, 0.866, 0));
        slipSystems_[1].schmidTensor = slipSystems_[1].m * slipSystems_[1].n;
        
        slipSystems_[2].n = basalNormal;
        slipSystems_[2].m = vector(-0.5, -0.866, 0) / mag(vector(-0.5, -0.866, 0));
        slipSystems_[2].schmidTensor = slipSystems_[2].m * slipSystems_[2].n;
        
        // Prismatic slip: {10-10}<11-20> - 3 systems
        slipSystems_[3].n = vector(1, 0, 0);
        slipSystems_[3].m = vector(0, 1, 0);
        slipSystems_[3].schmidTensor = slipSystems_[3].m * slipSystems_[3].n;
        
        slipSystems_[4].n = vector(-0.5, 0.866, 0);
        slipSystems_[4].m = vector(-0.866, -0.5, 0);
        slipSystems_[4].schmidTensor = slipSystems_[4].m * slipSystems_[4].n;
        
        slipSystems_[5].n = vector(-0.5, -0.866, 0);
        slipSystems_[5].m = vector(0.866, -0.5, 0);
        slipSystems_[5].schmidTensor = slipSystems_[5].m * slipSystems_[5].n;
        
        // Pyramidal <c+a> slip: {10-11}<11-23> - 6 systems (simplified)
        for (label i = 6; i < 12; i++)
        {
            scalar angle = (i - 6) * pi / 3.0;
            vector n = vector(cos(angle), sin(angle), 1.0/c_a);
            n = n / mag(n);
            
            vector m = vector(-sin(angle), cos(angle), c_a);
            m = m - (m & n) * n;
            m = m / mag(m);
            
            slipSystems_[i].n = n;
            slipSystems_[i].m = m;
            slipSystems_[i].schmidTensor = m * n;
        }
    }
    else
    {
        FatalErrorInFunction
            << "Unknown crystal type: " << crystalType << nl
            << "Valid types: FCC, BCC, HCP"
            << exit(FatalError);
    }
    
    // Initialize slip system state variables
    label nSlip = slipSystems_.size();
    gamma_.setSize(nSlip, 0.0);
    gamma_dot_.setSize(nSlip, 0.0);
    tau_.setSize(nSlip, 0.0);
    g_.setSize(nSlip, 0.0);
    
    Info<< "Initialized " << nSlip << " slip systems for " 
        << crystalType << " crystal" << endl;
}

void Foam::MaterialPointParcel::setCrystalOrientation
(
    scalar phi1, 
    scalar phi2, 
    scalar phi3
)
{
    /*
    Set crystal orientation from Bunge Euler angles (ZXZ convention)
    
    phi1: Rotation about Z
    Phi:  Rotation about X'
    phi2: Rotation about Z''
    */
    
    phi1_ = phi1;
    phi2_ = phi2;
    phi3_ = phi3;
    
    // // Convert to radians
    // scalar phi1Rad = phi1;// * pi / 180.0;
    // scalar phi2Rad = phi2;// * pi / 180.0;
    // scalar phi3Rad = phi3;// * pi / 180.0;
    
    // // Compute rotation matrix using Bunge Euler angles (ZXZ)
    // scalar c1 = cos(phi1Rad);
    // scalar s1 = sin(phi1Rad);
    // scalar c2 = cos(phi2Rad);
    // scalar s2 = sin(phi2Rad);
    // scalar c3 = cos(phi3Rad);
    // scalar s3 = sin(phi3Rad);
    
    auto R = buildRotationMatrix(phi1, phi2, phi3);

    Q0_ = toOpenFOAM(R);


}

void Foam::MaterialPointParcel::computeResolvedShearStress()
{
    /*
    Compute resolved shear stress on each slip system
    
    tau^alpha = (Fe * S0^alpha * Fe^T) : sigma
    
    where S0^alpha = m^alpha ⊗ n^alpha (Schmid tensor in reference config)
          sigma = Cauchy stress
          Fe = elastic deformation gradient
    */
    
    // Current rotation tensor from polar decomposition of Fe
    tensor R, U;
    polarDecomposition(Fe_, R, U);
    
    forAll(slipSystems_, alpha)
    {
        // Rotate slip system to current configuration
        vector m_current = R & (Q0_ & slipSystems_[alpha].m);
        vector n_current = R & (Q0_ & slipSystems_[alpha].n);
        
        // Schmid tensor in current configuration
        tensor S_current = m_current * n_current;
        
        // Resolved shear stress: tau = S : sigma
        tau_[alpha] = S_current && stress_;
    }
}

void Foam::MaterialPointParcel::updatePlasticDeformation
(
    const Material& material,
    const scalar& dt
)
{
    /*
    Update plastic deformation using rate-dependent crystal plasticity
    
    Power law flow rule:
    gamma_dot^alpha = gamma_dot_0 * |tau^alpha / g^alpha|^n * sign(tau^alpha)
    
    SPECIAL CASE: If ref_shear_rate = 0, plasticity is DISABLED
    */
    
    // =======================================================================
    // CHECK IF PLASTICITY IS DISABLED
    // =======================================================================
    
    // Method 1: Explicitly disabled by setting ref_shear_rate to zero
    if (material.ref_shear_rate <= 0.0)
    {
        // Plasticity is disabled - set all slip rates to zero
        forAll(gamma_dot_, alpha)
        {
            gamma_dot_[alpha] = 0.0;
        }
        
        // Early return - no plastic deformation
        return;
    }
    
    // Method 2: No slip systems defined
    if (slipSystems_.size() == 0)
    {
        // No slip systems - no plasticity possible
        return;
    }
    
    
    bool needsHardenessInit = false;
    forAll(g_, alpha)
    {
        if (g_[alpha] < material.initial_CRSS * 0.1)
        {
            needsHardenessInit = true;
            break;
        }
    }
    
    if (needsHardenessInit)
    {
        forAll(g_, alpha)
        {
            g_[alpha] = material.initial_CRSS;
        }
    }
    
    // =======================================================================
    // MATERIAL PARAMETERS (with validation)
    // =======================================================================
    
    scalar gamma_dot_0 = material.ref_shear_rate;
    scalar n_exp = material.k_exponent;
    
    // Validate parameters (but don't replace zero - it means disabled!)
    if (gamma_dot_0 < VSMALL && gamma_dot_0 > 0.0)
    {
        // Non-zero but very small - probably a mistake
        WarningInFunction
            << "ref_shear_rate is very small (but non-zero): " << gamma_dot_0 << nl
            << "This will make plasticity extremely slow. Did you mean 0.0 (disabled)?" << endl;
    }
    
    if (n_exp < 0.1)
    {
        WarningInFunction
            << "k_exponent is very small: " << n_exp << nl
            << "Clamping to 0.1 to avoid numerical issues" << endl;
        n_exp = 0.1;
    }
    else if (n_exp > 100.0)
    {
        WarningInFunction
            << "k_exponent is very large: " << n_exp << nl
            << "Clamping to 100 to avoid overflow" << endl;
        n_exp = 100.0;
    }
    
    // =======================================================================
    // COMPUTE SLIP RATES ON ALL SYSTEMS
    // =======================================================================
    
    scalar maxGammaDot = 0.0;
    
    forAll(slipSystems_, alpha)
    {
        gamma_dot_[alpha] = 0.0;  // Default to zero
        
        // Check hardness is valid
        if (g_[alpha] < material.initial_CRSS * 0.001)
        {
            continue;  // Skip this system
        }
        
        // Compute normalized resolved shear stress
        scalar tau_g = tau_[alpha] / g_[alpha];
        scalar tau_g_mag = mag(tau_g);
        
        // Check if slip system is active (τ > 0.5*g)
        if (tau_g_mag < 0.5)
        {
            gamma_dot_[alpha] = 0.0;
            continue;
        }
        
        // =======================================================================
        // ROBUST POWER OPERATION
        // =======================================================================
        
        scalar gamma_dot_magnitude = 0.0;
        
        // Clamp tau_g_mag to prevent overflow
        scalar tau_g_clamped = min(tau_g_mag, 100.0);
        
        // For large n_exp, use logarithms to avoid overflow
        if (n_exp > 10.0 && tau_g_clamped > 2.0)
        {
            // log(a^n) = n * log(a)
            scalar log_result = n_exp * Foam::log(tau_g_clamped);
            
            // Clamp to prevent overflow when exponentiating
            if (log_result > 100.0)
            {
                gamma_dot_magnitude = gamma_dot_0 * Foam::exp(100.0);
            }
            else
            {
                gamma_dot_magnitude = gamma_dot_0 * Foam::exp(log_result);
            }
        }
        else
        {
            // Safe to use pow directly
            gamma_dot_magnitude = gamma_dot_0 * Foam::pow(tau_g_clamped, n_exp);
        }
        
        // Apply sign of shear stress
        gamma_dot_[alpha] = gamma_dot_magnitude * sign(tau_g);
        
        maxGammaDot = max(maxGammaDot, mag(gamma_dot_[alpha]));
    }
    
    // =======================================================================
    // SAFETY CHECK ON SLIP RATES
    // =======================================================================
    
    if (maxGammaDot > 1.0e10)
    {
        WarningInFunction
            << "Very large slip rate detected: " << maxGammaDot << " 1/s" << nl
            << "Clamping to prevent instability" << endl;
        
        forAll(gamma_dot_, alpha)
        {
            gamma_dot_[alpha] = sign(gamma_dot_[alpha]) * 
                               min(mag(gamma_dot_[alpha]), 1.0e6);
        }
    }
    
    // =======================================================================
    // UPDATE ACCUMULATED SLIP
    // =======================================================================
    
    forAll(gamma_, alpha)
    {
        scalar dgamma = gamma_dot_[alpha] * dt;
        
        if (mag(dgamma) > 10.0)
        {
            WarningInFunction
                << "Excessive slip increment: " << dgamma << nl
                << "Time step may be too large!" << endl;
        }
        
        gamma_[alpha] += dgamma;
    }
    
    // =======================================================================
    // UPDATE PLASTIC DEFORMATION GRADIENT
    // =======================================================================
    
    tensor Lp = tensor::zero;
    
    forAll(slipSystems_, alpha)
    {
        if (mag(gamma_dot_[alpha]) < VSMALL) continue;
        
        // Rotate Schmid tensor to current configuration
        tensor R, U;
        polarDecomposition(Fe_, R, U);
        
        vector m_current = R & (Q0_ & slipSystems_[alpha].m);
        vector n_current = R & (Q0_ & slipSystems_[alpha].n);
        
        // Normalize
        m_current /= mag(m_current) + VSMALL;
        n_current /= mag(n_current) + VSMALL;
        
        tensor S_current = m_current * n_current;
        
        Lp += gamma_dot_[alpha] * S_current;
    }
    
    // Update Fp using forward Euler
    tensor dFp = dt * (Lp & Fp_);
    Fp_ = Fp_ + dFp;
    
    // Safety check on determinant
    scalar J_p = det(Fp_);
    if (J_p < 0.1 || J_p > 10.0)
    {
        WarningInFunction
            << "Plastic deformation gradient has extreme determinant: " << J_p << endl;
    }
    
    // Update elastic deformation gradient
    Fe_ = F_ & inv(Fp_);
    
    // Update plastic right Cauchy-Green tensor
    Cp_ = Fp_.T() & Fp_;
    
    // Update equivalent plastic strain
    scalar dgamma_eq = 0.0;
    forAll(gamma_dot_, alpha)
    {
        dgamma_eq += mag(gamma_dot_[alpha]);
    }
    epsilonP_eq_ += dgamma_eq * dt;
}

void Foam::MaterialPointParcel::updateCrystalRotation(const scalar dt)
{
    tensor Re, Ue;
    polarDecomposition(Fe_, Re, Ue);
    
    tensor Q_current = Re & Q0_;
    
    // FIX: Use ZYZ convention (NOT ZXZ!)
    // R = Rz(α) * Ry(β) * Rz(γ)
    
    scalar beta = Foam::acos(Foam::max(-1.0, Foam::min(1.0, Q_current.zz())));
    scalar alpha, gamma;
    
    if (Foam::mag(Foam::sin(beta)) > 1e-6)
    {
        // General case
        alpha = Foam::atan2(Q_current.yz(), Q_current.xz());   // FIX
        gamma = Foam::atan2(Q_current.zy(), -Q_current.zx());  // FIX
    }
    else
    {
        // Gimbal lock
        if (beta < constant::mathematical::pi / 2.0)
        {
            gamma = 0.0;
            alpha = Foam::atan2(-Q_current.xy(), Q_current.xx());
        }
        else
        {
            gamma = 0.0;
            alpha = Foam::atan2(Q_current.xy(), Q_current.xx());
        }
    }
    
    phi1_ = alpha;
    phi2_ = beta;
    phi3_ = gamma;
    
    Q0_ = Q_current;
}



void Foam::MaterialPointParcel::updateCrystalRotationIncremental
(
    const Material& material,
    const tensor& L,
    const scalar& dt
)
{
    /*
    Update crystal orientation using exponential map (Rodrigues formula)
    
    This is MORE ACCURATE than first-order Euler integration!
    
    Key improvements:
    1. Uses exp(We*dt) instead of I + We*dt
    2. Stores angles in RADIANS (no wrapping at 180°)
    3. Uses ZYZ Euler angle convention (more robust)
    */
    
    // Skip if no plastic deformation
    if (material.N_Slips == 0 || slipSystems_.size() == 0)
    {
        return;
    }
    
    // Compute plastic velocity gradient Lp from slip systems
    tensor Lp = tensor::zero;
    
    tensor R, U;
    polarDecomposition(Fe_, R, U);
    
    forAll(slipSystems_, alpha)
    {
        // Rotate slip system vectors to current configuration
        vector m_current = R & (Q0_ & slipSystems_[alpha].m);
        vector n_current = R & (Q0_ & slipSystems_[alpha].n);
        
        // Schmid tensor in current configuration
        tensor S_current = m_current * n_current;
        
        // Accumulate plastic velocity gradient
        Lp += gamma_dot_[alpha] * S_current;
    }
    
    // Extract plastic spin (antisymmetric part)
    tensor Wp = 0.5 * (Lp - Lp.T());
    
    // Total spin from velocity gradient
    tensor W = 0.5 * (L - L.T());
    
    // Elastic (lattice) spin - this rotates the crystal
    tensor We = W - Wp;
    
    // ========================================================================
    // EXPONENTIAL MAP FOR ROTATION UPDATE (Rodrigues Formula)
    // ========================================================================
    
    // Compute ||We|| = magnitude of rotation
    scalar X_comps = Foam::sqrt(
        Foam::sqr(We.zy()) +  // We(2,1)
        Foam::sqr(We.xz()) +  // We(0,2)
        Foam::sqr(We.yx())    // We(1,0)
    );
    
    // Compute We² for second-order term
    tensor We_sq = We & We;
    
    // Exponential map: exp(We*dt) = I + (sin(θ)/θ)*We*dt + ((1-cos(θ))/θ²)*(We*dt)²
    // where θ = ||We||*dt
    tensor exp_We_dt = tensor::I;
    
    if (X_comps > SMALL)
    {
        scalar theta = X_comps * dt;  // Total rotation angle
        
        // First-order term: (sin(θ)/θ) * We*dt
        scalar c1 = Foam::sin(theta) / X_comps;
        
        // Second-order term: ((1-cos(θ))/θ²) * (We*dt)²
        scalar c2 = (1.0 - Foam::cos(theta)) / Foam::sqr(X_comps);

        exp_We_dt += (c1 * We)  + (c2 * We_sq);
    }
    // If X_comps ≈ 0, exp(We*dt) ≈ I (no rotation)
    
    // Update rotation: Q_new = exp(We*dt) * Q_old
    tensor Q_current = exp_We_dt & Q0_;
    // toFastor(exp_We_dt);

    // Fastor::Tensor<double, 3, 3> Rot_matrix_N_plus_one = Fastor::matmul(toFastor(exp_We_dt), toFastor(Q0_));
    
    // tensor Q_current = toOpenFOAM(Rot_matrix_N_plus_one);
    // ========================================================================
    // ORTHONORMALIZATION (To maintain rotation matrix properties)
    // ========================================================================
    
    // Extract column vectors
    vector v1(Q_current.xx(), Q_current.yx(), Q_current.zx());
    vector v2(Q_current.xy(), Q_current.yy(), Q_current.zy());
    vector v3(Q_current.xz(), Q_current.yz(), Q_current.zz());
    
    // Gram-Schmidt orthonormalization
    v1 = v1 / (mag(v1) + VSMALL);
    v2 = v2 - (v2 & v1) * v1;
    v2 = v2 / (mag(v2) + VSMALL);
    v3 = v1 ^ v2;  // Cross product ensures right-handedness
    
    // Reconstruct rotation matrix
    Q_current = tensor
    (
        v1.x(), v2.x(), v3.x(),
        v1.y(), v2.y(), v3.y(),
        v1.z(), v2.z(), v3.z()
    );
    
    // Update reference orientation
    Q0_ = Q_current;
    
    // ========================================================================
// EXTRACT EULER ANGLES (ZYZ Convention) - STORE IN RADIANS!
// ========================================================================

scalar theta_angle, psi_angle, phi_angle;

// For ZYZ with R = Rz(α) * Ry(β) * Rz(γ):
// R(2,2) = cos(β)
// R(0,2) = cos(α)sin(β)
// R(1,2) = sin(α)sin(β)
// R(2,0) = -sin(β)cos(γ)
// R(2,1) = sin(β)sin(γ)

// Extract β (theta) from R(2,2)
scalar beta = Foam::acos(Foam::max(-1.0, Foam::min(1.0, Q_current.zz())));

// Check for gimbal lock: sin(β) ≈ 0
if (Foam::mag(Foam::sin(beta)) > SMALL)
{
    // General case: sin(β) ≠ 0
    
    // Extract α (psi) from R(0,2) and R(1,2)
    // α = atan2(R(1,2), R(0,2)) = atan2(sin(α)sin(β), cos(α)sin(β))
    scalar alpha = Foam::atan2(Q_current.yz(), Q_current.xz());
    
    // Extract γ (phi) from R(2,0) and R(2,1)
    // γ = atan2(R(2,1), -R(2,0)) = atan2(sin(β)sin(γ), sin(β)cos(γ))
    scalar gamma = Foam::atan2(Q_current.zy(), -Q_current.zx());
    
    psi_angle = alpha;
    theta_angle = beta;
    phi_angle = gamma;
}
else
{
    // Gimbal lock: sin(β) ≈ 0, so β ≈ 0 or β ≈ π
    
    if (beta < constant::mathematical::pi / 2.0)
    {
        // β ≈ 0: R becomes Rz(α+γ)
        // R(0,0) = cos(α+γ), R(0,1) = -sin(α+γ)
        // Set γ = 0 by convention
        phi_angle = 0.0;
        psi_angle = Foam::atan2(-Q_current.xy(), Q_current.xx());
        theta_angle = beta;  // Small positive value
    }
    else
    {
        // β ≈ π: R becomes Rz(α-γ)
        // R(0,0) = cos(α-γ), R(0,1) = -sin(α-γ)
        // Set γ = 0 by convention
        phi_angle = 0.0;
        psi_angle = Foam::atan2(-Q_current.xy(), Q_current.xx());
        theta_angle = beta;  // Close to π
    }
}

phi1_ = psi_angle;    // α (first Z rotation) in RADIANS
phi2_ = theta_angle;  // β (Y rotation) in RADIANS
phi3_ = phi_angle;    // γ (second Z rotation) in RADIANS


}


void Foam::MaterialPointParcel::updateHardening
(
    const Material& material,
    const scalar& dt
)
{
    /*
    Update hardness evolution on slip systems
    
    Voce hardening law:
    dg^alpha/dt = sum_beta( h_alpha_beta * |gamma_dot^beta| )
    
    h_alpha_beta = q * h_0 * (1 - g^beta/g_sat)  for alpha != beta (latent)
    h_alpha_alpha = h_0 * (1 - g^alpha/g_sat)     (self)
    
    where:
    - h_0: initial hardening modulus
    - g_sat: saturation hardness
    - q: latent hardening ratio
    */
    
    // Material parameters from Material object
    scalar g_0 = material.initial_CRSS;
    scalar g_sat = material.tau_S;
    scalar h_0 = material.h_0;
    scalar q = material.r;
    
    // CRITICAL SAFETY CHECK: Prevent division by zero
    if (g_sat < SMALL)
    {
        WarningInFunction
            << "tau_S (saturation strength) is too small or zero: " << g_sat << nl
            << "This will cause division by zero in hardening calculation." << nl
            << "Check your materialProperties dictionary!" << endl;
        return;  // Skip hardening update
    }
    
    if (g_0 < SMALL)
    {
        WarningInFunction
            << "initial_CRSS is too small or zero: " << g_0 << nl
            << "Slip system hardening will not work properly." << endl;
    }
    
    // Initialize hardness if needed
    bool needsInit = false;
    forAll(g_, alpha)
    {
        if (g_[alpha] < SMALL)
        {
            needsInit = true;
            break;
        }
    }
    
    if (needsInit)
    {
        forAll(g_, alpha)
        {
            g_[alpha] = g_0;
        }
    }
    
    // Compute hardening rates
    List<scalar> dg(slipSystems_.size(), 0.0);
    
    forAll(slipSystems_, alpha)
    {
        forAll(slipSystems_, beta)
        {
            scalar h_alpha_beta;
            
            if (alpha == beta)
            {
                // Self hardening
                h_alpha_beta = h_0 * (1.0 - g_[beta] / g_sat);
            }
            else
            {
                // Latent hardening
                h_alpha_beta = q * h_0 * (1.0 - g_[beta] / g_sat);
            }
            
            // Ensure hardening rate doesn't go negative
            h_alpha_beta = max(h_alpha_beta, 0.0);
            
            dg[alpha] += h_alpha_beta * mag(gamma_dot_[beta]);
        }
    }
    
    // Update hardness
    forAll(g_, alpha)
    {
        g_[alpha] += dg[alpha] * dt;
        
        // Ensure hardness stays within bounds
        g_[alpha] = max(g_[alpha], g_0);
        g_[alpha] = min(g_[alpha], g_sat);
    }
}


Fastor::Tensor<double, 3, 3> 
Foam::MaterialPointParcel::buildRotationMatrix
(
    const scalar phi1,  // α - first rotation about Z (RADIANS)
    const scalar phi2,  // β - rotation about Y (RADIANS)
    const scalar phi3   // γ - second rotation about Z (RADIANS)
) const
{
    /*
    Build rotation matrix from ZYZ Euler angles (in RADIANS)
    
    Convention: R = Rz(α) * Ry(β) * Rz(γ)
                  = Rz(phi1) * Ry(phi2) * Rz(phi3)
    
    This is the PROPER ZYZ Euler angle convention.
    
    NOTE: Angles are in RADIANS!
    */
    
    Fastor::Tensor<double, 3, 3> R(0.0);
    
    // Compute trig functions
    const scalar c1 = Foam::cos(phi1);  // cos(α)
    const scalar s1 = Foam::sin(phi1);  // sin(α)
    const scalar c2 = Foam::cos(phi2);  // cos(β)
    const scalar s2 = Foam::sin(phi2);  // sin(β)
    const scalar c3 = Foam::cos(phi3);  // cos(γ)
    const scalar s3 = Foam::sin(phi3);  // sin(γ)
    
    // ========================================================================
    // CORRECTED FORMULAS FOR ZYZ: R = Rz(α) * Ry(β) * Rz(γ)
    // ========================================================================
    
    // First row
    R(0, 0) = c1*c2*c3 - s1*s3;         // ← FIXED: was c2*c3
    R(0, 1) = -c1*c2*s3 - s1*c3;        // ← FIXED: was s1*s2*c3 - c1*s3
    R(0, 2) = c1*s2;                    // ← FIXED: was c1*s2*c3 + s1*s3
    
    // Second row
    R(1, 0) = s1*c2*c3 + c1*s3;         // ← FIXED: was c2*s3
    R(1, 1) = -s1*c2*s3 + c1*c3;        // ← FIXED: was s1*s2*s3 + c1*c3
    R(1, 2) = s1*s2;                    // ← FIXED: was c1*s2*s3 - s1*c3
    
    // Third row (these were mostly correct)
    R(2, 0) = -s2*c3;                   // ← FIXED: was -s2
    R(2, 1) = s2*s3;                    // ← FIXED: was s1*c2
    R(2, 2) = c2;                       // ← FIXED: was c1*c2
    
    return R;
}





// Fastor::Tensor<double, 3, 3> 
// Foam::MaterialPointParcel::buildRotationMatrix//THIS IS BUNGE CONVENTION
// (
//     const scalar phi1,
//     const scalar phi2,
//     const scalar phi3
// ) const
// {
//     Fastor::Tensor<double, 3, 3> R(0.0);
    
//     const scalar c1 = Foam::cos(phi1);
//     const scalar s1 = Foam::sin(phi1);
//     const scalar c2 = Foam::cos(phi2);
//     const scalar s2 = Foam::sin(phi2);
//     const scalar c3 = Foam::cos(phi3);
//     const scalar s3 = Foam::sin(phi3);
    
//     R(0, 0) = c2*c3;
//     R(0, 1) = s1*s2*c3 - c1*s3;
//     R(0, 2) = c1*s2*c3 + s1*s3;
//     R(1, 0) = c2*s3;
//     R(1, 1) = s1*s2*s3 + c1*c3;
//     R(1, 2) = c1*s2*s3 - s1*c3;
//     R(2, 0) = -s2;
//     R(2, 1) = s1*c2;
//     R(2, 2) = c1*c2;
    
//     return R;
// }

Fastor::Tensor<double, 3, 3, 3, 3> 
Foam::MaterialPointParcel::createCrystalStiffness
(
    const Material& material
) const
{
    Fastor::Tensor<double, 3, 3, 3, 3> C(0.0);
    
    const scalar C11 = material.C11;
    const scalar C12 = material.C12;
    const scalar C44 = material.C44;
    
    // Diagonal
    C(0,0,0,0) = C11; C(1,1,1,1) = C11; C(2,2,2,2) = C11;
    
    // Off-diagonal
    C(0,0,1,1) = C12; C(1,1,0,0) = C12;
    C(0,0,2,2) = C12; C(2,2,0,0) = C12;
    C(1,1,2,2) = C12; C(2,2,1,1) = C12;
    
    // Shear
    C(0,1,0,1) = C44; C(0,1,1,0) = C44; C(1,0,0,1) = C44; C(1,0,1,0) = C44;
    C(0,2,0,2) = C44; C(0,2,2,0) = C44; C(2,0,0,2) = C44; C(2,0,2,0) = C44;
    C(1,2,1,2) = C44; C(1,2,2,1) = C44; C(2,1,1,2) = C44; C(2,1,2,1) = C44;
    
    return C;
}

Fastor::Tensor<double, 3, 3, 3, 3> 
Foam::MaterialPointParcel::createIsotropicStiffness
(
    const scalar lambda,
    const scalar mu
) const
{
    Fastor::Tensor<double, 3, 3, 3, 3> C(0.0);
    
    for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
    for (int k = 0; k < 3; k++)
    for (int l = 0; l < 3; l++)
    {
        scalar dij = (i==j) ? 1.0 : 0.0;
        scalar dkl = (k==l) ? 1.0 : 0.0;
        scalar dik = (i==k) ? 1.0 : 0.0;
        scalar djl = (j==l) ? 1.0 : 0.0;
        scalar dil = (i==l) ? 1.0 : 0.0;
        scalar djk = (j==k) ? 1.0 : 0.0;
        
        C(i,j,k,l) = lambda*dij*dkl + mu*(dik*djl + dil*djk);
    }
    
    return C;
}

Fastor::Tensor<double, 3, 3, 3, 3> 
Foam::MaterialPointParcel::rotateStiffnessTensor
(
    const Material& material
) const
{
    /*
    Rotate the stiffness tensor from the crystal/material reference frame
    to the spatial/sample frame using the current crystal orientation
    
    C'_ijkl = R_im R_jn R_kp R_lq C_mnpq
    
    where R is the rotation matrix from Bunge Euler angles
    */
    
    // Get stiffness in material frame
    Fastor::Tensor<double, 3, 3, 3, 3> C_material = material.getStiffnessTensor();
    
    // For isotropic materials, rotation has no effect
    if (material.elasticType == "isotropic")
    {
        return C_material;
    }
    
    // Build rotation matrix from current orientation
    Fastor::Tensor<double, 3, 3> R = buildRotationMatrix(phi1_, phi2_, phi3_);
    
    // Rotate stiffness tensor: C'_ijkl = R_im R_jn R_kp R_lq C_mnpq
    Fastor::Tensor<double, 3, 3, 3, 3> C_spatial = 
        Fastor::einsum<
            Index<i,m>, Index<j,n>, Index<k,p>, Index<l,q>, Index<m,n,p,q>
        >(R, R, R, R, C_material);
    
    return C_spatial;
}

void Foam::MaterialPointParcel::computeElasticStress
(
    const Material& material
)
{
    // Convert Fe to Fastor
    Fastor::Tensor<double, 3, 3> FFe = toFastor(Fe_);
    
    // Compute elastic strain: Ee = 0.5*(Fe^T*Fe - I)
    Fastor::Tensor<double, 3, 3> I(0.0);
    I(0, 0) = 1.0;
    I(1, 1) = 1.0;
    I(2, 2) = 1.0;

    Fastor::Tensor<double, 3, 3> FEe = 0.5 * (
        Fastor::transpose(FFe) % FFe - I
    );
    
    // Get stiffness tensor (this handles all symmetries: isotropic, cubic, hexagonal, orthorhombic, anisotropic)
    Fastor::Tensor<double, 3, 3, 3, 3> FStiffness = material.getStiffnessTensor();
    
    // If anisotropic or crystal elasticity, rotate to current orientation
    if (material.elasticType != "isotropic")
    {
        // Rotation matrix from current orientation
        auto R = buildRotationMatrix(phi1_, phi2_, phi3_);
        
        // Rotate stiffness to sample frame: C' = R_im R_jn R_kp R_lq C_mnpq
        FStiffness = Fastor::einsum<
            Index<i,m>, Index<j,n>, Index<k,p>, Index<l,q>, Index<m,n,p,q>
        >(R, R, R, R, FStiffness);
    }
    
    // // Compute stress: sigma = C:Ee
    // Fastor::Tensor<double, 3, 3> FSigma = 
    //     Fastor::einsum<Index<i,j,k,l>, Index<k,l>>(FStiffness, FEe);
    
    // Compute 2nd PK stress: S = C:Ee
    Fastor::Tensor<double, 3, 3> S = 
        Fastor::einsum<Index<i,j,k,l>, Index<k,l>>(FStiffness, FEe);
    
        // PUSH-FORWARD to Cauchy stress: sigma = (1/J) * Fe * S * Fe^T
    scalar J = Fastor::determinant(FFe);

    Fastor::Tensor<double, 3, 3> FSigma = 
        (1.0/J) * Fastor::einsum<Index<i,k>, Index<k,l>, Index<j,l>>(
            FFe, S, FFe
        );
        
    // Convert back to OpenFOAM
    stress_ = toOpenFOAM(FSigma);
}





// void Foam::MaterialPointParcel::computeElasticStress
// (
//     const dictionary& materialProps
// )
// {
//     /*
//     Compute Cauchy stress from elastic deformation
    
//     For anisotropic elasticity:
//     S = C : Ee  (2nd Piola-Kirchhoff stress)
//     sigma = (1/J) * Fe * S * Fe^T  (Cauchy stress)
    
//     where Ee = 0.5*(Ce - I), Ce = Fe^T * Fe
//     */
    
//     // Check if anisotropic or isotropic
//     bool anisotropic = materialProps.lookupOrDefault<bool>("anisotropic", false);
    
//     if (anisotropic)
//     {
//         // Anisotropic elasticity (cubic symmetry)
//         scalar C11 = materialProps.lookupOrDefault<scalar>("C11", 168e9);
//         scalar C12 = materialProps.lookupOrDefault<scalar>("C12", 121e9);
//         scalar C44 = materialProps.lookupOrDefault<scalar>("C44", 75e9);
        
//         // Elastic right Cauchy-Green tensor
//         tensor Ce = Fe_.T() & Fe_;
        
//         // Green-Lagrange strain
//         tensor Ee = 0.5 * (Ce - tensor::I);
        
//         // Rotate strain to crystal frame
//         tensor Ee_crystal = Q0_.T() & Ee & Q0_;
        
//         // Compute 2nd PK stress in crystal frame using cubic elasticity
//         tensor S_crystal = tensor::zero;
        
//         // Voigt notation components
//         scalar E11 = Ee_crystal.xx();
//         scalar E22 = Ee_crystal.yy();
//         scalar E33 = Ee_crystal.zz();
//         scalar E23 = Ee_crystal.yz();
//         scalar E13 = Ee_crystal.xz();
//         scalar E12 = Ee_crystal.xy();
        
//         // Stress components
//         S_crystal.xx() = C11*E11 + C12*E22 + C12*E33;
//         S_crystal.yy() = C12*E11 + C11*E22 + C12*E33;
//         S_crystal.zz() = C12*E11 + C12*E22 + C11*E33;
//         S_crystal.yz() = 2.0*C44*E23;
//         S_crystal.xz() = 2.0*C44*E13;
//         S_crystal.xy() = 2.0*C44*E12;
        
//         // Make symmetric
//         S_crystal.zy() = S_crystal.yz();
//         S_crystal.zx() = S_crystal.xz();
//         S_crystal.yx() = S_crystal.xy();
        
//         // Rotate back to global frame
//         tensor S = Q0_ & S_crystal & Q0_.T();
        
//         // Push forward to Cauchy stress
//         scalar J = det(Fe_);
//         stress_ = (1.0 / J) * (Fe_ & S & Fe_.T());
//     }
//     else
//     {
//         // Isotropic neo-Hookean or St. Venant-Kirchhoff
//         scalar E = materialProps.lookupOrDefault<scalar>("E", 200e9);
//         scalar nu = materialProps.lookupOrDefault<scalar>("nu", 0.3);
        
//         scalar lambda = E * nu / ((1.0 + nu) * (1.0 - 2.0*nu));
//         scalar mu = E / (2.0 * (1.0 + nu));
        
//         // Elastic right Cauchy-Green tensor
//         tensor Ce = Fe_.T() & Fe_;
        
//         // Green-Lagrange strain
//         tensor Ee = 0.5 * (Ce - tensor::I);
        
//         // 2nd PK stress (St. Venant-Kirchhoff)
//         tensor S = lambda * tr(Ee) * tensor::I + 2.0 * mu * Ee;
        
//         // Push forward to Cauchy stress
//         scalar J = det(Fe_);
//         stress_ = (1.0 / J) * (Fe_ & S & Fe_.T());
//     }
    
//     // Update volume based on total deformation
//     scalar J_total = det(F_);
//     volume_ = volume0_ * J_total;
// }

void Foam::MaterialPointParcel::constitutiveUpdate
(
    const Material& material,
    const tensor& L,
    const scalar& dt
)
{
    /*
    Complete constitutive update for crystal plasticity
    
    FIXED: Ensures crystal rotation update happens for all materials
    
    Steps:
    1. Compute elastic stress
    2. Update crystal rotation for elastic materials
    3. Compute resolved shear stresses (if plastic)
    4. Update plastic deformation if active
    5. Recompute elastic stress if plastic
    6. Update hardening
    7. Update crystal rotation for plastic materials
    */
    
    // Step 1: Compute elastic stress
    computeElasticStress(material);
    

    if (material.N_Slips == 0 || slipSystems_.size() == 0)
    {
        // Elastic only - extract rotation from Fe
        updateCrystalRotation(dt);
        
        // Elastic-only material - we're done
        return;
    }
    
    // ========== PLASTIC MATERIALS ONLY BELOW THIS LINE ==========
    
    // Step 3: Compute RSS on all slip systems
    computeResolvedShearStress();
    
    // Step 4: Check if any system is active
    bool plasticActive = false;
    forAll(tau_, alpha)
    {
        if (g_[alpha] > SMALL && mag(tau_[alpha]) > 0.9 * g_[alpha])
        {
            plasticActive = true;
            break;
        }
    }
    
    // Step 5: If plastic, update Fp and recompute stress
    if (plasticActive)
    {
        // Sub-stepping for stability
        label nSubsteps = 10;
        scalar subDt = dt / nSubsteps;
        
        for (label sub = 0; sub < nSubsteps; sub++)
        {
            computeResolvedShearStress();
            updatePlasticDeformation(material, subDt);
            computeElasticStress(material);
        }
    }
    
    // Step 6: Update hardening
    updateHardening(material, dt);
    
    // Step 7: Update crystal rotation for plastic materials
    // Uses incremental method with elastic spin We = W - Wp
    updateCrystalRotationIncremental(material, L, dt);
    
    // Step 8: Compute plastic work
    scalar dWp = 0.0;
    forAll(tau_, alpha)
    {
        dWp += tau_[alpha] * gamma_dot_[alpha];
    }
    workPlastic_ += dWp * volume_ * dt;
}

// Helper function for polar decomposition
void Foam::MaterialPointParcel::polarDecomposition
(
    const tensor& F,
    tensor& R,
    tensor& U
) const
{
    /*
    Perform polar decomposition: F = R * U
    where R is rotation, U is right stretch
    */
    
    // Compute right Cauchy-Green tensor: C = F^T * F
    tensor C = F.T() & F;
    
    // Compute eigenvalues and eigenvectors of C
    vector eigenvalues;
    tensor eigenvectors;
    eigenDecomposition(C, eigenvalues, eigenvectors);
    
    // U = sqrt(C) = Q * sqrt(Lambda) * Q^T
    tensor Lambda_sqrt = tensor
    (
        Foam::sqrt(max(eigenvalues.x(), 0.0)), 0, 0,
        0, Foam::sqrt(max(eigenvalues.y(), 0.0)), 0,
        0, 0, Foam::sqrt(max(eigenvalues.z(), 0.0))
    );
    
    U = eigenvectors & Lambda_sqrt & eigenvectors.T();
    
    // R = F * U^-1
    R = F & inv(U);
}

void Foam::MaterialPointParcel::eigenDecomposition
(
    const tensor& T,
    vector& eigenvalues,
    tensor& eigenvectors
) const
{
    /*
    Compute eigenvalues and eigenvectors of symmetric 3x3 tensor
    Using Jacobi iteration method
    
    Simplified version using direct tensor component access
    */
    
    // Start with identity for eigenvectors
    eigenvectors = tensor::I;
    tensor A = T;
    
    const label maxIter = 50;
    const scalar tol = 1e-12;
    
    for (label iter = 0; iter < maxIter; iter++)
    {
        // Find largest off-diagonal element
        scalar maxVal = 0.0;
        label p = 0, q = 1;
        
        // Check xy
        if (mag(A.xy()) > maxVal)
        {
            maxVal = mag(A.xy());
            p = 0; q = 1;
        }
        // Check xz
        if (mag(A.xz()) > maxVal)
        {
            maxVal = mag(A.xz());
            p = 0; q = 2;
        }
        // Check yz
        if (mag(A.yz()) > maxVal)
        {
            maxVal = mag(A.yz());
            p = 1; q = 2;
        }
        
        // Check convergence
        if (maxVal < tol) break;
        
        // Get diagonal and off-diagonal elements
        scalar app, aqq, apq;
        
        if (p == 0 && q == 1)
        {
            app = A.xx(); aqq = A.yy(); apq = A.xy();
        }
        else if (p == 0 && q == 2)
        {
            app = A.xx(); aqq = A.zz(); apq = A.xz();
        }
        else // p == 1 && q == 2
        {
            app = A.yy(); aqq = A.zz(); apq = A.yz();
        }
        
        // Compute rotation angle
        scalar theta = 0.5 * atan2(2.0*apq, aqq - app);
        scalar c = cos(theta);
        scalar s = sin(theta);
        
        // Build rotation matrix G
        tensor G = tensor::I;
        
        if (p == 0 && q == 1)
        {
            G.xx() = c; G.xy() = s;
            G.yx() = -s; G.yy() = c;
        }
        else if (p == 0 && q == 2)
        {
            G.xx() = c; G.xz() = s;
            G.zx() = -s; G.zz() = c;
        }
        else // p == 1 && q == 2
        {
            G.yy() = c; G.yz() = s;
            G.zy() = -s; G.zz() = c;
        }
        
        // Apply rotation
        A = G.T() & A & G;
        eigenvectors = eigenvectors & G;
    }
    
    // Extract eigenvalues from diagonal
    eigenvalues.x() = A.xx();
    eigenvalues.y() = A.yy();
    eigenvalues.z() = A.zz();
}

// * * * * * * * * * * * * * * Shape Functions * * * * * * * * * * * * * * * //

/*---------------------------------------------------------------------------*\
  CORRECTED SHAPE FUNCTIONS - Use OpenFOAM's native hex node ordering
  
  Replace the shapeFunctions method in MaterialPointParcel.C with this version.
  
  Key fix: Trust OpenFOAM's node ordering instead of trying to remap based on
  positions. OpenFOAM hex cells have a consistent ordering we can rely on.
\*---------------------------------------------------------------------------*/

void Foam::MaterialPointParcel::shapeFunctions
(
    const fvMesh& mesh,
    scalarField& N,
    List<vector>& gradN,
    labelList& nodeIDs
) const
{


        // Skip if inactive
    if (!active_)
    {
        N.clear();
        gradN.clear();
        nodeIDs.clear();
        return;
    }


    
    label celli = cell();

    // Validate cell before accessing mesh
    if (celli < 0 || celli >= mesh.nCells())
    {
        // Particle is outside mesh or has invalid cell index
        N.setSize(0);
        gradN.setSize(0);
        nodeIDs.setSize(0);
        
        WarningInFunction
            << "Particle at position " << position() 
            << " has invalid cell index " << celli 
            << " (mesh has " << mesh.nCells() << " cells)" << endl;
        
        return;
    }


    const labelList& cellPoints = mesh.cellPoints()[celli];
    const pointField& points = mesh.points();
    
    label nNodes = cellPoints.size();
    nodeIDs = cellPoints;
    
    vector xp = position();
    
    // =========================================================================
    // HEXAHEDRAL CELL (8 nodes)
    // =========================================================================
    
    if (nNodes == 8)
{
    N.setSize(8);
    gradN.setSize(8);
    nodeIDs = cellPoints; // keep original ordering

    // compute axis-aligned bounding box (component-wise)
    vector xMin = points[cellPoints[0]];
    vector xMax = points[cellPoints[0]];
    forAll(cellPoints, i)
    {
        const vector& pt = points[cellPoints[i]];
        xMin = min(xMin, pt);
        xMax = max(xMax, pt);
    }
    vector dx = xMax - xMin;

    // degenerate check
    if (dx.x() < VSMALL || dx.y() < VSMALL || dx.z() < VSMALL)
    {
        WarningInFunction << "Degenerate cell detected! Cell " << celli
                          << " has dimensions " << dx << endl;
        N = scalar(1.0 / scalar(nNodes));
        gradN = Zero;
        return;
    }

    // particle natural coordinates in [0,1]
    vector xp = position();
    scalar xi   = (xp.x() - xMin.x()) / dx.x();
    scalar eta  = (xp.y() - xMin.y()) / dx.y();
    scalar zeta = (xp.z() - xMin.z()) / dx.z();

    xi   = max(0.0, min(1.0, xi));
    eta  = max(0.0, min(1.0, eta));
    zeta = max(0.0, min(1.0, zeta));

    // Prepare shape values for the 8 local *corner patterns* using bit encoding
    // but we will assign them to OpenFOAM node indices explicitly below.
    scalar N_local[8];
    vector gradN_local[8];

    for (label iLocal = 0; iLocal < 8; ++iLocal)
    {
        label i_xi   = (iLocal & 1);        // bit 0
        label i_eta  = (iLocal & 2) >> 1;   // bit 1
        label i_zeta = (iLocal & 4) >> 2;   // bit 2

        scalar f_xi   = (i_xi == 0) ? (1.0 - xi) : xi;
        scalar f_eta  = (i_eta == 0) ? (1.0 - eta) : eta;
        scalar f_zeta = (i_zeta == 0) ? (1.0 - zeta) : zeta;

        N_local[iLocal] = f_xi * f_eta * f_zeta;

        scalar dN_dxi   = ((i_xi == 0) ? -1.0 : 1.0) * f_eta * f_zeta;
        scalar dN_deta  = ((i_eta == 0) ? -1.0 : 1.0) * f_xi * f_zeta;
        scalar dN_dzeta = ((i_zeta == 0) ? -1.0 : 1.0) * f_xi * f_eta;

        gradN_local[iLocal] = vector(dN_dxi / dx.x(),
                                     dN_deta / dx.y(),
                                     dN_dzeta / dx.z());
    }

    // Build the canonical corner coordinates in OpenFOAM hex ordering:
    // (0) (xMin,yMin,zMin)
    // (1) (xMax,yMin,zMin)
    // (2) (xMax,yMax,zMin)
    // (3) (xMin,yMax,zMin)
    // (4) (xMin,yMin,zMax)
    // (5) (xMax,yMin,zMax)
    // (6) (xMax,yMax,zMax)
    // (7) (xMin,yMax,zMax)
    point corners[8];
    corners[0] = point(xMin.x(), xMin.y(), xMin.z());
    corners[1] = point(xMax.x(), xMin.y(), xMin.z());
    corners[2] = point(xMax.x(), xMax.y(), xMin.z());
    corners[3] = point(xMin.x(), xMax.y(), xMin.z());
    corners[4] = point(xMin.x(), xMin.y(), xMax.z());
    corners[5] = point(xMax.x(), xMin.y(), xMax.z());
    corners[6] = point(xMax.x(), xMax.y(), xMax.z());
    corners[7] = point(xMin.x(), xMax.y(), xMax.z());

    // Map the 8 canonical corners to the bit-encoded local indices used above.
    // The bit encoding we used (iLocal) corresponds to:
    // iLocalBitPattern -> (xBit, yBit, zBit) where bit==1 -> max
    // Create mapping from canonical FOAM corner index -> bit-encoded index:
    label foamToBitIndex[8];
    // For each foam corner index determine its bit index:
    for (label ifoam = 0; ifoam < 8; ++ifoam)
    {
        // determine bit pattern for that corner (compare to xMin/xMax)
        bool isXmax = ( fabs(corners[ifoam].x() - xMax.x()) < 1e-9*mag(dx) );
        bool isYmax = ( fabs(corners[ifoam].y() - xMax.y()) < 1e-9*mag(dx) );
        bool isZmax = ( fabs(corners[ifoam].z() - xMax.z()) < 1e-9*mag(dx) );

        label bitIndex = 0;
        if (isXmax) bitIndex += 1;
        if (isYmax) bitIndex += 2;
        if (isZmax) bitIndex += 4;
        foamToBitIndex[ifoam] = bitIndex;
    }

    // Now for each actual physical node in cellPoints find which canonical corner it matches,
    // and pick the corresponding shape function from N_local[bitIndex].
    scalar tol = 1e-6 * mag(dx);
    forAll(cellPoints, iPhys)
    {
        const point& pos = points[cellPoints[iPhys]];
        bool matched = false;

        for (label ifoam = 0; ifoam < 8; ++ifoam)
        {
            if ( mag(vector(pos - corners[ifoam])) < tol )
            {
                label bitIdx = foamToBitIndex[ifoam];
                N[iPhys] = N_local[bitIdx];
                gradN[iPhys] = gradN_local[bitIdx];
                matched = true;
                break;
            }
        }
        if (!matched)
        {
            // fallback: assign small weight to avoid leaving uninitialised
            N[iPhys] = 0.0;
            gradN[iPhys] = Zero;
            WarningInFunction << "Could not match hexahedron vertex " << iPhys
                              << " at " << pos << " to expected corner positions."
                              << " Using zero weight for this node." << endl;
        }
    }
} // end hex
    
    // =========================================================================
    // TETRAHEDRAL CELL (4 nodes)
    // =========================================================================
    
    else if (nNodes == 4)
    {
        N.setSize(4);
        gradN.setSize(4);
        
        vector v0 = points[cellPoints[0]];
        vector v1 = points[cellPoints[1]];
        vector v2 = points[cellPoints[2]];
        vector v3 = points[cellPoints[3]];
        
        vector e1 = v1 - v0;
        vector e2 = v2 - v0;
        vector e3 = v3 - v0;
        
        scalar V6 = mag(e1 & (e2 ^ e3));
        scalar V = V6 / 6.0;
        
        if (V < VSMALL)
        {
            N = 0.25;
            gradN = Zero;
            return;
        }
        
        scalar V0 = mag((xp - v1) & ((v2 - v1) ^ (v3 - v1))) / V6;
        scalar V1 = mag((xp - v0) & ((v2 - v0) ^ (v3 - v0))) / V6;
        scalar V2 = mag((xp - v0) & ((v1 - v0) ^ (v3 - v0))) / V6;
        scalar V3 = mag((xp - v0) & ((v1 - v0) ^ (v2 - v0))) / V6;
        
        N[0] = V0;
        N[1] = V1;
        N[2] = V2;
        N[3] = V3;
        
        scalar sumN = N[0] + N[1] + N[2] + N[3];
        if (sumN > VSMALL)
        {
            N /= sumN;
        }
        
        gradN[0] = ((v2 - v1) ^ (v3 - v1)) / V6;
        gradN[1] = ((v3 - v0) ^ (v2 - v0)) / V6;
        gradN[2] = ((v1 - v0) ^ (v3 - v0)) / V6;
        gradN[3] = ((v2 - v0) ^ (v1 - v0)) / V6;
    }
    
    // =========================================================================
    // OTHER CELL TYPES
    // =========================================================================
    
    else
    {
        WarningInFunction
            << "Unsupported cell type with " << nNodes << " nodes. "
            << "Using inverse distance weighting." << endl;
        
        N.setSize(nNodes);
        gradN.setSize(nNodes);
        
        scalarField weights(nNodes);
        scalar sumWeights = 0.0;
        
        forAll(cellPoints, i)
        {
            vector dx = xp - points[cellPoints[i]];
            scalar dist = mag(dx) + VSMALL;
            weights[i] = 1.0 / dist;
            sumWeights += weights[i];
        }
        
        N = weights / sumWeights;
        
        scalar h = Foam::pow(mesh.V()[celli], 1.0/3.0);
        forAll(cellPoints, i)
        {
            vector dx = xp - points[cellPoints[i]];
            gradN[i] = -N[i] * dx / (h * h + VSMALL);
        }
    }
}

// ************************************************************************* //

// * * * * * * * * * * * * * * P2G Transfer * * * * * * * * * * * * * * * * //

void Foam::MaterialPointParcel::transferToGrid
(
    const fvMesh& mesh,
    pointScalarField& nodalMass,
    pointVectorField& nodalMomentum,
    pointTensorField& nodalStress,
    pointVectorField& nodalForce
) const
{
    /*
    Transfer particle data to grid nodes (P2G)
    
    Mass: m_i = sum_p( N_ip * m_p )
    Momentum: (mv)_i = sum_p( N_ip * m_p * v_p )
    Force: f_i = sum_p( -V_p * sigma_p * grad(N_ip) )
    */

        if (!active_)
    {
        return;
    }


    
    // Get shape functions
    scalarField N;
    List<vector> gradN;
    labelList nodeIDs;
    
    shapeFunctions(mesh, N, gradN, nodeIDs);
    
    // Transfer to each node
    forAll(nodeIDs, i)
    {
        label nodeI = nodeIDs[i];
        
        // Mass transfer
        nodalMass[nodeI] += N[i] * mass_;
        
        // Momentum transfer
        nodalMomentum[nodeI] += N[i] * mass_ * velocity_;
        
        // Internal force from stress divergence
        // f_int = -V * sigma * grad(N)
        nodalForce[nodeI] -= volume_ * (stress_ & gradN[i]);
        
        // Stress transfer (for visualization)
        nodalStress[nodeI] += N[i] * stress_ * volume_;
    }
}

// * * * * * * * * * * * * * * G2P Transfer * * * * * * * * * * * * * * * * //

void Foam::MaterialPointParcel::transferFromGrid
(
    const fvMesh& mesh,
    const pointVectorField& nodalVelocity,
    const pointVectorField& nodalAcceleration
)
{
    /*
    Transfer grid data back to particles (G2P)
    
    v_p = sum_i( N_ip * v_i )
    a_p = sum_i( N_ip * a_i )
    */
    
    // Get shape functions
    scalarField N;
    List<vector> gradN;
    labelList nodeIDs;
    
    shapeFunctions(mesh, N, gradN, nodeIDs);
    
    // Interpolate velocity and acceleration
    velocity_ = Zero;
    vector acceleration = Zero;
    
    forAll(nodeIDs, i)
    {
        label nodeI = nodeIDs[i];
        
        velocity_ += N[i] * nodalVelocity[nodeI];
        acceleration += N[i] * nodalAcceleration[nodeI];
    }
}

// * * * * * * * * * * * * * * Deformation Gradient * * * * * * * * * * * * //

void Foam::MaterialPointParcel::updateDeformationGradient
(
    const fvMesh& mesh,
    const pointVectorField& nodalVelocity,
    const scalar dt
)
{
    /*
    Update total deformation gradient
    
    F_new = (I + dt * grad(v)) * F_old
    
    where grad(v) = sum_i( v_i ⊗ grad(N_i) )
    */
    
    // Get shape functions
    scalarField N;
    List<vector> gradN;
    labelList nodeIDs;
    
    shapeFunctions(mesh, N, gradN, nodeIDs);
    
    // Compute velocity gradient
    tensor L = tensor::zero;
    
    forAll(nodeIDs, i)
    {
        label nodeI = nodeIDs[i];
        L += nodalVelocity[nodeI] * gradN[i];
    }
    
    // Update F using forward Euler
    tensor dF = dt * (L & F_);
    F_ = F_ + dF;
    
    // Update Fe (Fp is updated in constitutive update)
    Fe_ = F_ & inv(Fp_);
}




// bool Foam::MaterialPointParcel::hitPatch
// (
//     materialPointCloud& cloud,
//     trackingData& td
// )
// {
//     /*
//     Called when particle hits ANY boundary patch.
    
//     Return false: particle leaves domain (will be deleted)
//     Return true:  particle stays in domain
//     */
    
//     const polyMesh& mesh = cloud.pMesh();
//     const polyBoundaryMesh& boundaryMesh = mesh.boundaryMesh();
    
//     // Get patch information
//     label patchI = patch();
    
//     if (patchI < 0 || patchI >= boundaryMesh.size())
//     {
//         // Invalid patch - shouldn't happen, but be safe
//         WarningInFunction
//             << "Invalid patch index: " << patchI << endl;
//         td.keepParticle = true;
//         return true;
//     }
    
//     const polyPatch& pp = boundaryMesh[patchI];
//     word patchName = pp.name();
//     word patchType = pp.type();
    
//     // Try to read particle boundary configuration
//     IOdictionary boundaryDict
//     (
//         IOobject
//         (
//             "particleBoundaryConditions",
//             mesh.time().constant(),
//             mesh,
//             IOobject::READ_IF_PRESENT,
//             IOobject::NO_WRITE
//         )
//     );
    
//     // Determine behavior type
//     word behaviorType = patchType;  // Default to mesh patch type
    
//     if (boundaryDict.headerOk() && boundaryDict.found("particleBoundaryTypes"))
//     {
//         const dictionary& boundaryTypes = 
//             boundaryDict.subDict("particleBoundaryTypes");
        
//         if (boundaryTypes.found(patchName))
//         {
//             const dictionary& patchDict = boundaryTypes.subDict(patchName);
//             behaviorType = word(patchDict.lookup("type"));
//         }
//         else if (boundaryDict.found("defaultBehavior"))
//         {
//             const dictionary& defaultDict = boundaryDict.subDict("defaultBehavior");
//             behaviorType = word(defaultDict.lookup("type"));
//         }
//     }
    
//     // =========================================================================
//     // HANDLE DIFFERENT BOUNDARY TYPES
//     // =========================================================================
    
//     if (behaviorType == "outflow" || behaviorType == "outlet")
//     {
//         // ---------------------------------------------------------------------
//         // OUTFLOW: Particle exits domain
//         // ---------------------------------------------------------------------
        
//         Info<< "    Particle leaving domain through " << patchName << endl;
        
//         // Mark particle as outside domain
//         cell() = -1;
        
//         // Tell tracking system to delete particle
//         td.keepParticle = false;
//         td.switchProcessor = false;
        
//         return false;  // Don't keep particle
//     }
//     else if (behaviorType == "wall")
//     {
//         // ---------------------------------------------------------------------
//         // WALL: Delegate to wall handler
//         // ---------------------------------------------------------------------
        
//         return hitWallPatch(cloud, td);
//     }
//     else if (behaviorType == "symmetry" || behaviorType == "symmetryPlane")
//     {
//         // ---------------------------------------------------------------------
//         // SYMMETRY: Mirror particle velocity
//         // ---------------------------------------------------------------------
        
//         // Get face normal
//         label faceI = face();
//         label patchFaceI = faceI - pp.start();
        
//         if (patchFaceI < 0 || patchFaceI >= pp.size())
//         {
//             WarningInFunction
//                 << "Invalid face index on patch " << patchName << endl;
//             td.keepParticle = true;
//             return true;
//         }
        
//         vector n = pp.faceNormals()[patchFaceI];
//         n /= mag(n) + VSMALL;
        
//         // Reflect velocity: v' = v - 2*(v·n)*n
//         scalar vn = velocity_ & n;
//         if (vn < 0)  // Moving into boundary
//         {
//             velocity_ -= 2.0 * vn * n;
//         }
        
//         td.keepParticle = true;
//         return true;
//     }
//     else if (behaviorType == "empty" || behaviorType == "wedge")
//     {
//         // ---------------------------------------------------------------------
//         // EMPTY/WEDGE: For 2D/axisymmetric cases
//         // ---------------------------------------------------------------------
        
//         // No action needed - particle stays in domain
//         td.keepParticle = true;
//         return true;
//     }
//     else if (behaviorType == "processor" || behaviorType == "processorCyclic")
//     {
//         // ---------------------------------------------------------------------
//         // PROCESSOR: Parallel domain decomposition
//         // ---------------------------------------------------------------------
        
//         // Let OpenFOAM handle processor boundary transfer
//         td.keepParticle = true;
//         td.switchProcessor = true;
//         return true;
//     }
//     else if (behaviorType == "cyclic" || behaviorType == "cyclicAMI")
//     {
//         // ---------------------------------------------------------------------
//         // CYCLIC: Periodic boundaries
//         // ---------------------------------------------------------------------
        
//         // Let OpenFOAM handle cyclic transfer
//         td.keepParticle = true;
//         return true;
//     }
//     else
//     {
//         // ---------------------------------------------------------------------
//         // UNKNOWN: Default to wall behavior for safety
//         // ---------------------------------------------------------------------
        
//         WarningInFunction
//             << "Unknown boundary type '" << behaviorType 
//             << "' for patch " << patchName << nl
//             << "Treating as wall (particle will stick)" << endl;
        
//         return hitWallPatch(cloud, td);
//     }
// }






// bool Foam::MaterialPointParcel::hitWallPatch
// (
//     materialPointCloud& cloud,
//     trackingData& td
// )
// {
//     /*
//     Called when particle hits a WALL patch.
    
//     Behavior options:
//     - stick:   particle stops at wall (v = 0)
//     - reflect: particle bounces off wall
//     - slip:    particle slides along wall (vn = 0, vt maintained)
//     */
    
//     const polyMesh& mesh = cloud.pMesh();
//     const polyBoundaryMesh& boundaryMesh = mesh.boundaryMesh();
    
//     label patchI = patch();
//     const polyPatch& pp = boundaryMesh[patchI];
//     word patchName = pp.name();
    
//     // Read wall behavior from configuration
//     IOdictionary boundaryDict
//     (
//         IOobject
//         (
//             "particleBoundaryConditions",
//             mesh.time().constant(),
//             mesh,
//             IOobject::READ_IF_PRESENT,
//             IOobject::NO_WRITE
//         )
//     );
    
//     // Default wall behavior
//     word wallBehavior = "stick";
//     scalar restitution = 0.0;  // Coefficient of restitution (0=inelastic, 1=elastic)
//     scalar friction = 0.0;     // Tangential friction coefficient
    
//     if (boundaryDict.headerOk() && boundaryDict.found("particleBoundaryTypes"))
//     {
//         const dictionary& boundaryTypes = 
//             boundaryDict.subDict("particleBoundaryTypes");
        
//         if (boundaryTypes.found(patchName))
//         {
//             const dictionary& patchDict = boundaryTypes.subDict(patchName);
            
//             wallBehavior = patchDict.lookupOrDefault<word>("behavior", "stick");
//             restitution = patchDict.lookupOrDefault<scalar>("restitution", 0.0);
//             friction = patchDict.lookupOrDefault<scalar>("friction", 0.0);
//         }
//     }
    
//     // Get face normal
//     label faceI = face();
//     label patchFaceI = faceI - pp.start();
    
//     if (patchFaceI < 0 || patchFaceI >= pp.size())
//     {
//         WarningInFunction
//             << "Invalid face index on patch " << patchName << endl;
//         velocity_ = Zero;
//         td.keepParticle = true;
//         return true;
//     }
    
//     vector n = pp.faceNormals()[patchFaceI];
//     n /= mag(n) + VSMALL;
    
//     // Decompose velocity into normal and tangential components
//     scalar vn = velocity_ & n;              // Normal component
//     vector vt = velocity_ - vn * n;         // Tangential component
    
//     // Apply wall behavior
//     if (wallBehavior == "stick")
//     {
//         // =====================================================================
//         // STICK: Particle adheres to wall (zero velocity)
//         // =====================================================================
        
//         velocity_ = Zero;
        
//         td.keepParticle = true;
//         return true;
//     }
//     else if (wallBehavior == "reflect")
//     {
//         // =====================================================================
//         // REFLECT: Elastic or inelastic collision
//         // =====================================================================
        
//         if (vn < 0)  // Moving INTO wall
//         {
//             // Reflect normal component with restitution
//             scalar vn_new = -restitution * vn;
            
//             // Apply tangential friction if specified
//             vector vt_new = vt;
//             if (friction > SMALL && mag(vt) > SMALL)
//             {
//                 // Simple Coulomb friction model
//                 scalar vt_mag = mag(vt);
//                 scalar vt_reduced = max(0.0, vt_mag - friction * mag(vn));
//                 vt_new = (vt / vt_mag) * vt_reduced;
//             }
            
//             // Reconstruct velocity
//             velocity_ = vn_new * n + vt_new;
//         }
        
//         td.keepParticle = true;
//         return true;
//     }
//     else if (wallBehavior == "slip")
//     {
//         // =====================================================================
//         // SLIP: Free slip (remove normal component, keep tangential)
//         // =====================================================================
        
//         if (vn < 0)  // Moving INTO wall
//         {
//             // Zero normal component, maintain tangential
//             velocity_ = vt;
//         }
        
//         td.keepParticle = true;
//         return true;
//     }
//     else
//     {
//         // =====================================================================
//         // UNKNOWN: Default to stick for safety
//         // =====================================================================
        
//         WarningInFunction
//             << "Unknown wall behavior '" << wallBehavior 
//             << "' for patch " << patchName << nl
//             << "Using 'stick' as default" << endl;
        
//         velocity_ = Zero;
//         td.keepParticle = true;
//         return true;
//     }
// }










bool Foam::MaterialPointParcel::hitPatch
(
    materialPointCloud&,
    trackingData& td
)
{
    // Get the patch that was hit
    const polyPatch& patch = mesh().boundaryMesh()[this->patch()];
    
    // Check patch type
    if (isA<processorPolyPatch>(patch))
    {
        // Processor boundary - flag for transfer
        td.switchProcessor = true;
        td.keepParticle = true;
        return true;
    }
    else
    {
        // Physical boundary - deactivate particle
        active_ = false;
        td.keepParticle = true;
        td.switchProcessor = false;
        return true;
    }
}

bool Foam::MaterialPointParcel::hitWallPatch
(
    materialPointCloud&,
    trackingData& td
)
{
    // For MPM, particles hitting walls are deactivated
    active_ = false;
    td.keepParticle = true;
    td.switchProcessor = false;
    return true;
}

void Foam::MaterialPointParcel::prepareForParallelTransfer()
{
    // This is called BEFORE particle data is sent via MPI
    // Nothing special needed for MPM - all data is already in member variables
}

void Foam::MaterialPointParcel::correctAfterParallelTransfer
(
    const label patchi,
    trackingData& td
)
{
    // This is called AFTER particle data is received via MPI
    // Update particle's patch index if needed
    particle::correctAfterParallelTransfer(patchi, td);
}










// void Foam::MaterialPointParcel::transformProperties(const transformer& transform)
// {
//     /*
//     Transform particle properties when crossing processor/periodic boundaries.
    
//     This is called automatically by OpenFOAM's particle tracking when a particle
//     moves between processors or through periodic boundaries.
//     */
    
//     // Transform base particle properties
//     particle::transformProperties(transform);
    
//     // Transform vectors
//     velocity_ = transform.transform(velocity_);
    
//     // Transform second-order tensors (stress, deformation gradients)
//     stress_ = transform.transform(stress_);
//     F_ = transform.transform(F_);
//     Fp_ = transform.transform(Fp_);
//     Fe_ = transform.transform(Fe_);
//     Cp_ = transform.transform(Cp_);
//     Q0_ = transform.transform(Q0_);
    
//     // Transform slip system geometry (if initialized)
//     forAll(slipSystems_, i)
//     {
//         slipSystems_[i].m = transform.transform(slipSystems_[i].m);
//         slipSystems_[i].n = transform.transform(slipSystems_[i].n);
//         slipSystems_[i].schmidTensor = transform.transform(slipSystems_[i].schmidTensor);
//     }
    
//     // Note: Scalar quantities (mass, volume, temperature, etc.) don't need transformation
// }




// * * * * * * * * * * * * * * Particle Movement * * * * * * * * * * * * * * //

bool Foam::MaterialPointParcel::move
(
    materialPointCloud& cloud,
    trackingData& td,
    const scalar dt
)
{
    // Skip inactive particles
    if (!active_)
    {
        td.keepParticle = true;
        td.switchProcessor = false;
        return true;
    }
    
    td.keepParticle = true;
    td.switchProcessor = false;
    
    // Calculate displacement
    vector displacement = velocity_ * dt;
    
    // Skip if displacement is negligible
    if (mag(displacement) < VSMALL)
    {
        return true;
    }
    
    // Use OpenFOAM's built-in tracking - this handles cell crossing properly
    scalar trackFraction = 1.0;
    
    while (trackFraction > SMALL)
    {
        // Calculate remaining displacement
        vector remainingDisp = displacement * trackFraction;
        
        // Track particle using base class function
        // This automatically updates position and cell
        trackFraction = track(remainingDisp, 0);
        
        // Check if we hit a boundary
        if (trackFraction > SMALL)
        {
            // We didn't complete the move - hit something
            
            if (onBoundaryFace())
            {
                // Get boundary face
                label facei = face();
                label patchi = patch();
                
                if (patchi < 0 || patchi >= mesh().boundaryMesh().size())
                {
                    // Invalid patch
                    active_ = false;
                    td.keepParticle = true;
                    return true;
                }
                
                const polyPatch& pp = mesh().boundaryMesh()[patchi];
                
                if (isA<processorPolyPatch>(pp))
                {
                    // Hit processor boundary
                    td.switchProcessor = true;
                    td.keepParticle = true;
                    return true;
                }
                else
                {
                    // Hit physical boundary - deactivate
                    active_ = false;
                    td.keepParticle = true;
                    return true;
                }
            }
            else
            {
                // Not on boundary but didn't complete move - stuck
                // This shouldn't happen, but handle it
                active_ = false;
                td.keepParticle = true;
                return true;
            }
        }
    }
    
    // Successfully completed move
    return true;
}



// // Alternative simpler version if the above has compilation issues:
// bool Foam::MaterialPointParcel::move
// (
//     materialPointCloud& cloud,
//     trackingData& td,
//     const scalar trackTime
// )
// {
//     td.switchProcessor = false;
//     td.keepParticle = true;

//     const scalar tEnd = (1.0 - stepFraction())*trackTime;
//     scalar dtMax = tEnd;

//     while (td.keepParticle && dtMax > ROOTVSMALL)
//     {
//         // Set the Lagrangian time-step
//         scalar dt = min(dtMax, tEnd);

//         // Remember which cell the particle is in
//         // since this will change if a boundary is hit
//         const label celli = cell();
        
//         const scalar magU = mag(velocity_);
//         if (magU < ROOTVSMALL)
//         {
//             // Particle is stationary
//             break;
//         }

//         // Update particle position
//         scalar f = trackToAndHitFace(dt * velocity_, 0.0, cloud, td);

//         // Check if particle is on a boundary
//         if (onBoundaryFace())
//         {
//             // Check processor patch
//             if (isA<processorPolyPatch>(mesh().boundaryMesh()[patch()]))
//             {
//                 td.switchProcessor = true;
//             }
            
//             if (!td.switchProcessor)
//             {
//                 // Handle other boundaries
//                 hitPatch(cloud, td);
//             }
//         }

//         dtMax = (1.0 - f)*dt;
//     }

//     return td.keepParticle;
// }


// * * * * * * * * * * * * * * * IO Functions * * * * * * * * * * * * * * * //

void Foam::MaterialPointParcel::readFields(materialPointCloud& c)
{
    if (!c.size())
    {
        return;
    }
    
    particle::readFields(c);
    
    // Read custom fields
    IOField<bool> active(c.fieldIOobject("active", IOobject::MUST_READ));
    IOField<scalar> mass(c.fieldIOobject("mass", IOobject::MUST_READ));
    IOField<scalar> volume(c.fieldIOobject("volume", IOobject::MUST_READ));
    IOField<scalar> volume0(c.fieldIOobject("volume0", IOobject::MUST_READ));
    IOField<vector> velocity(c.fieldIOobject("velocity", IOobject::MUST_READ));
    IOField<tensor> stress(c.fieldIOobject("stress", IOobject::MUST_READ));
    IOField<tensor> F(c.fieldIOobject("F", IOobject::MUST_READ));
    IOField<tensor> Fp(c.fieldIOobject("Fp", IOobject::MUST_READ));
    IOField<tensor> Fe(c.fieldIOobject("Fe", IOobject::MUST_READ));
    IOField<label> materialID(c.fieldIOobject("materialID", IOobject::MUST_READ));
    IOField<label> grainID(c.fieldIOobject("grainID", IOobject::MUST_READ));
    IOField<scalar> phi1(c.fieldIOobject("phi1", IOobject::MUST_READ));
    IOField<scalar> phi2(c.fieldIOobject("phi2", IOobject::MUST_READ));
    IOField<scalar> phi3(c.fieldIOobject("phi3", IOobject::MUST_READ));
    IOField<tensor> Q0(c.fieldIOobject("Q0", IOobject::MUST_READ));
    IOField<scalar> epsilonP_eq(c.fieldIOobject("epsilonP_eq", IOobject::MUST_READ));
    
    
    // Assign to particles
    label i = 0;
    forAllIter(materialPointCloud, c, iter)
    {
        MaterialPointParcel& p = iter();
        
        p.active_ = active[i];
        p.mass_ = mass[i];
        p.volume_ = volume[i];
        p.volume0_ = volume0[i];
        p.velocity_ = velocity[i];
        p.stress_ = stress[i];
        p.F_ = F[i];
        p.Fp_ = Fp[i];
        p.Fe_ = Fe[i];
        p.materialID_ = materialID[i];
        p.grainID_ = grainID[i];
        p.phi1_ = phi1[i];
        p.phi2_ = phi2[i];
        p.phi3_ = phi3[i];
        p.Q0_ = Q0[i];
        p.epsilonP_eq_ = epsilonP_eq[i];
        
        
        i++;
    }
}

void Foam::MaterialPointParcel::writeFields(const materialPointCloud& c)
{
    particle::writeFields(c);
    
    const label np = c.size();
    
    // Create IOFields for all custom data
    IOField<bool> active(c.fieldIOobject("active", IOobject::NO_READ), np);
    IOField<scalar> mass(c.fieldIOobject("mass", IOobject::NO_READ), np);
    IOField<scalar> volume(c.fieldIOobject("volume", IOobject::NO_READ), np);
    IOField<scalar> volume0(c.fieldIOobject("volume0", IOobject::NO_READ), np);
    IOField<vector> velocity(c.fieldIOobject("velocity", IOobject::NO_READ), np);
    IOField<tensor> stress(c.fieldIOobject("stress", IOobject::NO_READ), np);
    IOField<tensor> F(c.fieldIOobject("F", IOobject::NO_READ), np);
    IOField<tensor> Fp(c.fieldIOobject("Fp", IOobject::NO_READ), np);
    IOField<tensor> Fe(c.fieldIOobject("Fe", IOobject::NO_READ), np);
    IOField<tensor> Cp(c.fieldIOobject("Cp", IOobject::NO_READ), np);
    IOField<label> materialID(c.fieldIOobject("materialID", IOobject::NO_READ), np);
    IOField<label> grainID(c.fieldIOobject("grainID", IOobject::NO_READ), np);
    IOField<scalar> phi1(c.fieldIOobject("phi1", IOobject::NO_READ), np);
    IOField<scalar> phi2(c.fieldIOobject("phi2", IOobject::NO_READ), np);
    IOField<scalar> phi3(c.fieldIOobject("phi3", IOobject::NO_READ), np);
    IOField<tensor> Q0(c.fieldIOobject("Q0", IOobject::NO_READ), np);
    IOField<scalar> epsilonP_eq(c.fieldIOobject("epsilonP_eq", IOobject::NO_READ), np);
    IOField<scalar> workPlastic(c.fieldIOobject("workPlastic", IOobject::NO_READ), np);
    
    
    // Write slip system data if particles have slip systems initialized
    bool hasSlipSystems = false;
    label nSlipSystems = 0;
    
    forAllConstIter(materialPointCloud, c, iter)
    {
        if (iter().gamma_.size() > 0)
        {
            hasSlipSystems = true;
            nSlipSystems = iter().gamma_.size();
            break;
        }
    }
    
    // Collect data from particles
    label i = 0;
    forAllConstIter(materialPointCloud, c, iter)
    {
        const MaterialPointParcel& p = iter();
        
        active[i] = p.active_;
        mass[i] = p.mass_;
        volume[i] = p.volume_;
        volume0[i] = p.volume0_;
        velocity[i] = p.velocity_;
        stress[i] = p.stress_;
        F[i] = p.F_;
        Fp[i] = p.Fp_;
        Fe[i] = p.Fe_;
        Cp[i] = p.Cp_;
        materialID[i] = p.materialID_;
        grainID[i] = p.grainID_;
        phi1[i] = p.phi1_;
        phi2[i] = p.phi2_;
        phi3[i] = p.phi3_;
        Q0[i] = p.Q0_;
        epsilonP_eq[i] = p.epsilonP_eq_;
        workPlastic[i] = p.workPlastic_;
        
        
        i++;
    }
    
    // Write all fields
    active.write();
    mass.write();
    volume.write();
    volume0.write();
    velocity.write();
    stress.write();
    F.write();
    Fp.write();
    Fe.write();
    Cp.write();
    materialID.write();
    grainID.write();
    phi1.write();
    phi2.write();
    phi3.write();
    Q0.write();
    epsilonP_eq.write();
    workPlastic.write();
    
    
    // Write slip system data if available
    if (hasSlipSystems && nSlipSystems > 0)
    {
        // Write accumulated slip for each system
        for (label alpha = 0; alpha < nSlipSystems; alpha++)
        {
            word fieldName = "gamma" + Foam::name(alpha);
            IOField<scalar> gamma(c.fieldIOobject(fieldName, IOobject::NO_READ), np);
            
            i = 0;
            forAllConstIter(materialPointCloud, c, iter)
            {
                if (iter().gamma_.size() > alpha)
                {
                    gamma[i] = iter().gamma_[alpha];
                }
                else
                {
                    gamma[i] = 0.0;
                }
                i++;
            }
            
            gamma.write();
        }
        
        // Write slip rates
        for (label alpha = 0; alpha < nSlipSystems; alpha++)
        {
            word fieldName = "gamma_dot" + Foam::name(alpha);
            IOField<scalar> gamma_dot(c.fieldIOobject(fieldName, IOobject::NO_READ), np);
            
            i = 0;
            forAllConstIter(materialPointCloud, c, iter)
            {
                if (iter().gamma_dot_.size() > alpha)
                {
                    gamma_dot[i] = iter().gamma_dot_[alpha];
                }
                else
                {
                    gamma_dot[i] = 0.0;
                }
                i++;
            }
            
            gamma_dot.write();
        }
        
        // Write resolved shear stresses
        for (label alpha = 0; alpha < nSlipSystems; alpha++)
        {
            word fieldName = "tau" + Foam::name(alpha);
            IOField<scalar> tau(c.fieldIOobject(fieldName, IOobject::NO_READ), np);
            
            i = 0;
            forAllConstIter(materialPointCloud, c, iter)
            {
                if (iter().tau_.size() > alpha)
                {
                    tau[i] = iter().tau_[alpha];
                }
                else
                {
                    tau[i] = 0.0;
                }
                i++;
            }
            
            tau.write();
        }
        
        // Write hardness
        for (label alpha = 0; alpha < nSlipSystems; alpha++)
        {
            word fieldName = "g" + Foam::name(alpha);
            IOField<scalar> g(c.fieldIOobject(fieldName, IOobject::NO_READ), np);
            
            i = 0;
            forAllConstIter(materialPointCloud, c, iter)
            {
                if (iter().g_.size() > alpha)
                {
                    g[i] = iter().g_[alpha];
                }
                else
                {
                    g[i] = 0.0;
                }
                i++;
            }
            
            g.write();
        }
        
        Info<< "    Wrote " << nSlipSystems << " slip system fields per particle" << endl;
    }
    
    Info<< "    Wrote " << np << " particle fields" << endl;
}

// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

Foam::Ostream& Foam::operator<<(Ostream& os, const MaterialPointParcel& p)
{
    // Write base particle data
    os  << static_cast<const particle&>(p);
    
    // Write material point data
    os  
        << token::SPACE << p.active_
        << token::SPACE << p.mass_
        << token::SPACE << p.volume0_
        << token::SPACE << p.volume_
        << token::SPACE << p.velocity_
        << token::SPACE << p.stress_
        << token::SPACE << p.F_
        << token::SPACE << p.Fp_
        << token::SPACE << p.Fe_
        << token::SPACE << p.Cp_
        << token::SPACE << p.materialID_
        << token::SPACE << p.grainID_
        << token::SPACE << p.phi1_
        << token::SPACE << p.phi2_
        << token::SPACE << p.phi3_
        << token::SPACE << p.Q0_
        << token::SPACE << p.epsilonP_eq_
        << token::SPACE << p.workPlastic_
        << token::SPACE << p.temperature_;;
    
    // Write slip system data
    os  << token::SPACE << p.gamma_.size();
    
    if (p.gamma_.size() > 0)
    {
        os  << token::SPACE << p.gamma_
            << token::SPACE << p.gamma_dot_
            << token::SPACE << p.tau_
            << token::SPACE << p.g_;
    }
    
    return os;
}

// ************************************************************************* //