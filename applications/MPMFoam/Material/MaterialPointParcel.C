/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | 
    \\  /    A nd           | 
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

\*---------------------------------------------------------------------------*/

#include "MaterialPointParcel.H"
#include "IOstreams.H"

// * * * * * * * * * * * * * * Constructors * * * * * * * * * * * * * * //

Foam::MaterialPointParcel::MaterialPointParcel
(
    const polyMesh& mesh,
    const vector& position,
    const label celli
)
:
    particle(mesh, position, celli),
    mass_(0.0),
    volume0_(0.0),
    volume_(0.0),
    velocity_(Zero),
    stress_(Zero),
    F_(I),
    Fp_(I),
    Fe_(I),
    Cp_(I),
    materialID_(0),
    grainID_(0),
    phi1_(0.0),
    phi2_(0.0),
    phi3_(0.0),
    Q0_(I),
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
    // particle(mesh, is, readFields, newFormat)
    particle(mesh, is, readFields)
{
    if (readFields)
    {
        is >> mass_ >> volume0_ >> volume_ >> velocity_ >> stress_ 
           >> F_ >> Fp_ >> Fe_ >> Cp_
           >> materialID_ >> grainID_
           >> phi1_ >> phi2_ >> phi3_ >> Q0_
           >> gamma_ >> gamma_dot_ >> tau_ >> g_
           >> epsilonP_eq_ >> workPlastic_ >> temperature_;
    }
}

// * * * * * * * * * * * Shape Function Calculations * * * * * * * * * * //



void Foam::MaterialPointParcel::shapeFunctions
(
    const fvMesh& mesh,
    scalarField& N,
    List<vector>& gradN,
    labelList& nodeIDs
) const
{
    label cellI = cell();
    // const cell& c = mesh.cells()[cellI];
    
    // Get cell vertices (nodes)
    const labelList& cellPoints = mesh.cellPoints()[cellI];
    label nNodes = cellPoints.size();
    
    nodeIDs = cellPoints;
    N.setSize(nNodes);
    gradN.setSize(nNodes);
    
    const vector& particlePos = position();
    
    // Get cell bounds
    vector cellMin(GREAT, GREAT, GREAT);
    vector cellMax(-GREAT, -GREAT, -GREAT);
    
    forAll(cellPoints, i)
    {
        const vector& pt = mesh.points()[cellPoints[i]];
        cellMin = min(cellMin, pt);
        cellMax = max(cellMax, pt);
    }
    
    vector cellSize = cellMax - cellMin;
    
    // Normalized particle position within cell [-1,1]
    // Component-wise division
    vector xi;
    xi.x() = 2.0 * (particlePos.x() - cellMin.x()) / (cellSize.x() + SMALL) - 1.0;
    xi.y() = 2.0 * (particlePos.y() - cellMin.y()) / (cellSize.y() + SMALL) - 1.0;
    xi.z() = 2.0 * (particlePos.z() - cellMin.z()) / (cellSize.z() + SMALL) - 1.0;
    
    // Clamp to ensure particle is within cell
    xi.x() = max(-1.0, min(1.0, xi.x()));
    xi.y() = max(-1.0, min(1.0, xi.y()));
    xi.z() = max(-1.0, min(1.0, xi.z()));
    
    // Trilinear shape functions for hexahedron
    // Node numbering: 0-3 bottom, 4-7 top (standard OpenFOAM hex)
    if (nNodes == 8)
    {
        // Shape functions
        N[0] = 0.125 * (1-xi.x()) * (1-xi.y()) * (1-xi.z());
        N[1] = 0.125 * (1+xi.x()) * (1-xi.y()) * (1-xi.z());
        N[2] = 0.125 * (1+xi.x()) * (1+xi.y()) * (1-xi.z());
        N[3] = 0.125 * (1-xi.x()) * (1+xi.y()) * (1-xi.z());
        N[4] = 0.125 * (1-xi.x()) * (1-xi.y()) * (1+xi.z());
        N[5] = 0.125 * (1+xi.x()) * (1-xi.y()) * (1+xi.z());
        N[6] = 0.125 * (1+xi.x()) * (1+xi.y()) * (1+xi.z());
        N[7] = 0.125 * (1-xi.x()) * (1+xi.y()) * (1+xi.z());
        
        // Shape function derivatives in parametric space
        // dN/dxi, dN/deta, dN/dzeta
        List<vector> dN_dxi(8);
        
        dN_dxi[0] = vector(-0.125*(1-xi.y())*(1-xi.z()), -0.125*(1-xi.x())*(1-xi.z()), -0.125*(1-xi.x())*(1-xi.y()));
        dN_dxi[1] = vector( 0.125*(1-xi.y())*(1-xi.z()), -0.125*(1+xi.x())*(1-xi.z()), -0.125*(1+xi.x())*(1-xi.y()));
        dN_dxi[2] = vector( 0.125*(1+xi.y())*(1-xi.z()),  0.125*(1+xi.x())*(1-xi.z()), -0.125*(1+xi.x())*(1+xi.y()));
        dN_dxi[3] = vector(-0.125*(1+xi.y())*(1-xi.z()),  0.125*(1-xi.x())*(1-xi.z()), -0.125*(1-xi.x())*(1+xi.y()));
        dN_dxi[4] = vector(-0.125*(1-xi.y())*(1+xi.z()), -0.125*(1-xi.x())*(1+xi.z()),  0.125*(1-xi.x())*(1-xi.y()));
        dN_dxi[5] = vector( 0.125*(1-xi.y())*(1+xi.z()), -0.125*(1+xi.x())*(1+xi.z()),  0.125*(1+xi.x())*(1-xi.y()));
        dN_dxi[6] = vector( 0.125*(1+xi.y())*(1+xi.z()),  0.125*(1+xi.x())*(1+xi.z()),  0.125*(1+xi.x())*(1+xi.y()));
        dN_dxi[7] = vector(-0.125*(1+xi.y())*(1+xi.z()),  0.125*(1-xi.x())*(1+xi.z()),  0.125*(1-xi.x())*(1+xi.y()));
        
        // Transform derivatives to physical space: ∇N = J^(-T) * ∇_ξ N
        // where J is the Jacobian matrix dx/dxi
        
        // Simple approximation: assume uniform cell
        // J^(-T) ≈ diag(2/Lx, 2/Ly, 2/Lz)
        forAll(gradN, nodeI)
        {
            gradN[nodeI].x() = dN_dxi[nodeI].x() * 2.0 / (cellSize.x() + SMALL);
            gradN[nodeI].y() = dN_dxi[nodeI].y() * 2.0 / (cellSize.y() + SMALL);
            gradN[nodeI].z() = dN_dxi[nodeI].z() * 2.0 / (cellSize.z() + SMALL);
        }
    }
    else
    {
        // For other cell shapes, use simple inverse distance weighting
        forAll(N, nodeI)
        {
            vector toNode = mesh.points()[cellPoints[nodeI]] - particlePos;
            scalar dist = mag(toNode) + SMALL;
            
            // Inverse distance weighting
            N[nodeI] = 1.0 / (dist + SMALL);
            gradN[nodeI] = -toNode / (dist * dist * dist + SMALL);
        }
        
        // Normalize
        scalar sumN = sum(N) + SMALL;
        N /= sumN;
    }
}

// void Foam::MaterialPointParcel::shapeFunctions
// (
//     const fvMesh& mesh,
//     scalarField& N,
//     List<vector>& gradN,
//     labelList& nodeIDs
// ) const
// {
//     // Get particle position and cell
//     const vector& xp = position();
//     label celli = cell();
    
//     // Validate cell index
//     if (celli < 0 || celli >= mesh.nCells())
//     {
//         FatalErrorInFunction
//             << "Particle has invalid cell index: " << celli << nl
//             << "Particle position: " << xp << nl
//             << exit(FatalError);
//     }
    
//     // Get nodes of this cell
//     const labelList& cellPoints = mesh.cellPoints()[celli];
//     const label nNodes = cellPoints.size();
    
//     // Resize output arrays
//     N.setSize(nNodes);
//     gradN.setSize(nNodes);
//     nodeIDs.setSize(nNodes);
    
//     // =================================================================
//     // INVERSE DISTANCE WEIGHTING
//     // N_i = (1/r_i) / sum(1/r_j)
//     // where r_i = distance from particle to node i
//     // =================================================================
    
//     scalar sumWeights = 0.0;
//     scalar minDist = GREAT;
//     label closestNode = -1;
    
//     // First pass: compute weights
//     forAll(cellPoints, i)
//     {
//         label pointI = cellPoints[i];
//         nodeIDs[i] = pointI;
        
//         // Get node position
//         const vector& xn = mesh.points()[pointI];
        
//         // Distance from particle to node
//         scalar dist = mag(xp - xn);
        
//         // Track closest node
//         if (dist < minDist)
//         {
//             minDist = dist;
//             closestNode = i;
//         }
        
//         // Inverse distance weight (with small regularization)
//         scalar weight = 1.0 / (dist + 1e-20);
//         N[i] = weight;
//         sumWeights += weight;
//     }
    
//     // Special case: particle very close to a node
//     if (minDist < 1e-10)
//     {
//         // Use delta function (all weight on closest node)
//         N = 0.0;
//         N[closestNode] = 1.0;
//         gradN = vector::zero;
//         return;
//     }
    
//     // Normalize shape functions (now they sum to 1)
//     N /= sumWeights;
    
//     // =================================================================
//     // COMPUTE GRADIENTS
//     // grad(N_i) = -N_i * (xp - xi) / |xp - xi|^2
//     // =================================================================
    
//     forAll(cellPoints, i)
//     {
//         label pointI = cellPoints[i];
//         const vector& xn = mesh.points()[pointI];
        
//         // Vector from node to particle
//         vector dx = xp - xn;
//         scalar r2 = magSqr(dx) + 1e-20;
        
//         // Gradient of inverse distance weight
//         gradN[i] = -N[i] * dx / r2;
//     }
    
//     // Optional: Verify shape functions (can comment out after debugging)
//     #ifdef DEBUG_SHAPE_FUNCTIONS
//     scalar sumN = sum(N);
//     if (mag(sumN - 1.0) > 1e-6)
//     {
//         WarningInFunction
//             << "Shape functions don't sum to 1: sum = " << sumN << nl
//             << "Particle: " << xp << ", Cell: " << celli << endl;
//     }
//     #endif
// }


/*---------------------------------------------------------------------------*\
  ALTERNATIVE: Equal weights (simplest possible, for debugging)
  
  Use this version if you just want to test that P2G mapping works.
  Less accurate but absolutely guaranteed to work.
  
  To use: Replace shapeFunctions() with this shapeFunctionsEqual()
          and rename it to shapeFunctions()
\*---------------------------------------------------------------------------*/
/*
void Foam::MaterialPointParcel::shapeFunctionsEqual
(
    const fvMesh& mesh,
    scalarField& N,
    List<vector>& gradN,
    labelList& nodeIDs
) const
{
    const vector& xp = position();
    label celli = cell();
    
    const labelList& cellPoints = mesh.cellPoints()[celli];
    const label nNodes = cellPoints.size();
    
    N.setSize(nNodes);
    gradN.setSize(nNodes);
    nodeIDs.setSize(nNodes);
    
    // Equal weight for all nodes
    const scalar weight = 1.0 / scalar(nNodes);
    
    forAll(cellPoints, i)
    {
        label pointI = cellPoints[i];
        nodeIDs[i] = pointI;
        
        // Shape function: equal for all
        N[i] = weight;
        
        // Gradient: simple approximation
        const vector& xn = mesh.points()[pointI];
        vector dx = xp - xn;
        scalar r2 = magSqr(dx) + 1e-20;
        gradN[i] = -weight * dx / r2;
    }
}
*/



/*
void Foam::MaterialPointParcel::shapeFunctions
(
    const fvMesh& mesh,
    scalarField& N,
    List<vector>& gradN,
    labelList& nodeIDs
) const
{
    // Get particle position and cell
    const vector& xp = position();
    label celli = cell();
    
    if (celli < 0)
    {
        FatalErrorInFunction
            << "Particle has invalid cell index: " << celli
            << exit(FatalError);
    }
    
    // Get cell points (nodes of the cell containing particle)
    const labelList& cellPoints = mesh.cellPoints()[celli];
    label nNodes = cellPoints.size();
    
    // Resize arrays
    N.setSize(nNodes);
    gradN.setSize(nNodes);
    nodeIDs.setSize(nNodes);
    
    // Get cell center and characteristic size
    const vector& xc = mesh.C()[celli];
    const scalar h = pow(mesh.V()[celli], 1.0/3.0);  // Characteristic cell size
    
    // For each node in the cell
    forAll(cellPoints, i)
    {
        label pointI = cellPoints[i];
        nodeIDs[i] = pointI;
        
        // Get node position
        const vector& xn = mesh.points()[pointI];
        
        // Relative position
        vector xi = (xp - xn) / h;  // Normalized by cell size
        
        // Linear shape function (trilinear interpolation)
        // N(ξ) = (1 - |ξx|)(1 - |ξy|)(1 - |ξz|)  for |ξi| < 1
        
        scalar Nx = max(0.0, 1.0 - mag(xi.x()));
        scalar Ny = max(0.0, 1.0 - mag(xi.y()));
        scalar Nz = max(0.0, 1.0 - mag(xi.z()));
        
        N[i] = Nx * Ny * Nz;
        
        // Gradient of shape function
        // ∂N/∂x = ∂N/∂ξ * ∂ξ/∂x = ∂N/∂ξ / h
        
        scalar dNdx = -sign(xi.x()) * Ny * Nz / h;
        scalar dNdy = -sign(xi.y()) * Nx * Nz / h;
        scalar dNdz = -sign(xi.z()) * Nx * Ny / h;
        
        gradN[i] = vector(dNdx, dNdy, dNdz);
    }
    
    // Normalize shape functions (should sum to 1)
    scalar sumN = sum(N);
    if (sumN > SMALL)
    {
        N /= sumN;
        // Note: gradN normalization is more complex, skipped for now
    }
    else
    {
        WarningInFunction
            << "Shape functions sum to zero for particle at " << xp << nl
            << "Cell: " << celli << ", sumN = " << sumN << endl;
    }
}
*/

/*
void Foam::MaterialPointParcel::shapeFunctions
(
    const fvMesh& mesh,
    scalarField& N,
    List<vector>& gradN,
    labelList& nodeIDs
) const
{
    const label celli = cell();
    
    if (celli < 0)
    {
        N.clear();
        gradN.clear();
        nodeIDs.clear();
        return;
    }
    
    // Get mesh connectivity
    const pointField& points = mesh.points();
    
    // Collect unique vertices for this cell using pointCells
    // This is a simpler approach that avoids the cell type issue
    labelHashSet uniquePoints;
    
    // Get all points in the cell - simpler approach
    const labelList& cellPoints = mesh.cellPoints()[celli];
    
    forAll(cellPoints, pointi)
    {
        uniquePoints.insert(cellPoints[pointi]);
    }
    
    nodeIDs = uniquePoints.toc();
    const label nNodes = nodeIDs.size();
    
    N.setSize(nNodes);
    gradN.setSize(nNodes);
    
    // For hex cells (8 nodes), use trilinear shape functions
    if (nNodes == 8)
    {
        // Find min/max coordinates to establish reference element
        vector minPt(GREAT, GREAT, GREAT);
        vector maxPt(-GREAT, -GREAT, -GREAT);
        
        forAll(nodeIDs, i)
        {
            const vector& pt = points[nodeIDs[i]];
            minPt = min(minPt, pt);
            maxPt = max(maxPt, pt);
        }
        
        // Map particle position to reference coordinates
        const vector& xp = position();
        vector xi = cmptDivide(xp - minPt, maxPt - minPt);
        
        // Clamp to [0,1]
        xi.x() = max(0.0, min(1.0, xi.x()));
        xi.y() = max(0.0, min(1.0, xi.y()));
        xi.z() = max(0.0, min(1.0, xi.z()));
        
        // Order nodes according to hex convention (simplified)
        // This assumes nodes are already reasonably ordered
        
        // Standard trilinear shape functions
        N[0] = (1.0 - xi.x()) * (1.0 - xi.y()) * (1.0 - xi.z());
        N[1] = xi.x() * (1.0 - xi.y()) * (1.0 - xi.z());
        N[2] = xi.x() * xi.y() * (1.0 - xi.z());
        N[3] = (1.0 - xi.x()) * xi.y() * (1.0 - xi.z());
        N[4] = (1.0 - xi.x()) * (1.0 - xi.y()) * xi.z();
        N[5] = xi.x() * (1.0 - xi.y()) * xi.z();
        N[6] = xi.x() * xi.y() * xi.z();
        N[7] = (1.0 - xi.x()) * xi.y() * xi.z();
        
        // Shape function gradients
        vector h = maxPt - minPt;
        
        // Protect against division by zero
        if (mag(h.x()) < SMALL) h.x() = SMALL;
        if (mag(h.y()) < SMALL) h.y() = SMALL;
        if (mag(h.z()) < SMALL) h.z() = SMALL;
        
        gradN[0] = vector(-(1.0 - xi.y()) * (1.0 - xi.z()) / h.x(),
                         -(1.0 - xi.x()) * (1.0 - xi.z()) / h.y(),
                         -(1.0 - xi.x()) * (1.0 - xi.y()) / h.z());
        
        gradN[1] = vector((1.0 - xi.y()) * (1.0 - xi.z()) / h.x(),
                         -xi.x() * (1.0 - xi.z()) / h.y(),
                         -xi.x() * (1.0 - xi.y()) / h.z());
        
        gradN[2] = vector(xi.y() * (1.0 - xi.z()) / h.x(),
                         xi.x() * (1.0 - xi.z()) / h.y(),
                         -xi.x() * xi.y() / h.z());
        
        gradN[3] = vector(-xi.y() * (1.0 - xi.z()) / h.x(),
                         (1.0 - xi.x()) * (1.0 - xi.z()) / h.y(),
                         -(1.0 - xi.x()) * xi.y() / h.z());
        
        gradN[4] = vector(-(1.0 - xi.y()) * xi.z() / h.x(),
                         -(1.0 - xi.x()) * xi.z() / h.y(),
                         (1.0 - xi.x()) * (1.0 - xi.y()) / h.z());
        
        gradN[5] = vector((1.0 - xi.y()) * xi.z() / h.x(),
                         -xi.x() * xi.z() / h.y(),
                         xi.x() * (1.0 - xi.y()) / h.z());
        
        gradN[6] = vector(xi.y() * xi.z() / h.x(),
                         xi.x() * xi.z() / h.y(),
                         xi.x() * xi.y() / h.z());
        
        gradN[7] = vector(-xi.y() * xi.z() / h.x(),
                         (1.0 - xi.x()) * xi.z() / h.y(),
                         (1.0 - xi.x()) * xi.y() / h.z());
    }
    else
    {
        // For non-hex cells, use inverse distance weighting (simplified)
        scalar sumW = 0.0;
        forAll(nodeIDs, i)
        {
            scalar dist = mag(position() - points[nodeIDs[i]]);
            if (dist < SMALL)
            {
                N = scalarField(nNodes, 0.0);
                N[i] = 1.0;
                gradN = List<vector>(nNodes, Zero);
                return;
            }
            N[i] = 1.0 / dist;
            sumW += N[i];
        }
        
        // Normalize
        if (sumW > SMALL)
        {
            N /= sumW;
        }
        
        // Simplified gradient calculation
        forAll(nodeIDs, i)
        {
            vector diff = position() - points[nodeIDs[i]];
            scalar dist = mag(diff);
            if (dist > SMALL)
            {
                gradN[i] = -diff / (dist * dist * dist * sumW);
            }
            else
            {
                gradN[i] = Zero;
            }
        }
    }
}
*/




void Foam::MaterialPointParcel::setCrystalOrientation
(
    scalar phi1_deg,
    scalar Phi_deg,
    scalar phi2_deg
)
{
    // Store Euler angles (in degrees for reference)
    phi1_ = phi1_deg;
    phi2_ = Phi_deg;   // Note: phi2_ stores Phi
    phi3_ = phi2_deg;  // Note: phi3_ stores phi2
    
    // Convert to radians for computation
    const scalar pi = constant::mathematical::pi;
    scalar phi1 = phi1_deg * pi / 180.0;
    scalar Phi  = Phi_deg  * pi / 180.0;
    scalar phi2 = phi2_deg * pi / 180.0;
    
    // Compute rotation matrix using Bunge convention
    // R = Rz(phi2) * Rx(Phi) * Rz(phi1)
    //
    // This transforms from crystal coordinates to sample coordinates
    
    scalar c1 = Foam::cos(phi1);
    scalar s1 = Foam::sin(phi1);
    scalar cP = Foam::cos(Phi);
    scalar sP = Foam::sin(Phi);
    scalar c2 = Foam::cos(phi2);
    scalar s2 = Foam::sin(phi2);
    
    // First rotation: Rz(phi1)
    tensor Rz1
    (
        c1, -s1,  0,
        s1,  c1,  0,
         0,   0,  1
    );
    
    // Second rotation: Rx(Phi)
    tensor Rx
    (
        1,   0,   0,
        0,  cP, -sP,
        0,  sP,  cP
    );
    
    // Third rotation: Rz(phi2)
    tensor Rz2
    (
        c2, -s2,  0,
        s2,  c2,  0,
         0,   0,  1
    );
    
    // Combined rotation matrix
    Q0_ = Rz2 & Rx & Rz1;
    
    // Ensure it's a proper rotation (det = 1, orthogonal)
    scalar detQ = det(Q0_);
    if (mag(detQ - 1.0) > 1e-6)
    {
        WarningInFunction
            << "Rotation matrix determinant = " << detQ 
            << " (should be 1.0)" << endl;
    }
}

// Alternative implementation using compact form:
/*
void Foam::MaterialPointParcel::setCrystalOrientation
(
    scalar phi1_deg,
    scalar Phi_deg,
    scalar phi2_deg
)
{
    // Store angles
    phi1_ = phi1_deg;
    phi2_ = Phi_deg;
    phi3_ = phi2_deg;
    
    // Convert to radians
    const scalar pi = constant::mathematical::pi;
    scalar p1 = phi1_deg * pi / 180.0;
    scalar P  = Phi_deg  * pi / 180.0;
    scalar p2 = phi2_deg * pi / 180.0;
    
    // Bunge Euler angles to rotation matrix (compact form)
    scalar c1 = cos(p1), s1 = sin(p1);
    scalar cP = cos(P),  sP = sin(P);
    scalar c2 = cos(p2), s2 = sin(p2);
    
    Q0_ = tensor
    (
        c1*c2 - s1*s2*cP,  -c1*s2 - s1*c2*cP,   s1*sP,
        s1*c2 + c1*s2*cP,  -s1*s2 + c1*c2*cP,  -c1*sP,
        s2*sP,              c2*sP,              cP
    );
}
*/


// Simplified GIMP shape functions stub
void Foam::MaterialPointParcel::shapeFunctionsGIMP
(
    const fvMesh& mesh,
    const scalar& lp,
    scalarField& N,
    List<vector>& gradN,
    labelList& nodeIDs
) const
{
    // For now, just use standard shape functions
    shapeFunctions(mesh, N, gradN, nodeIDs);
    // Full GIMP implementation would go here
}

// * * * * * * * * * Particle-Grid Transfer Operations * * * * * * * * * //

void Foam::MaterialPointParcel::transferToGrid
(
    const fvMesh& mesh,
    pointScalarField& nodalMass,
    pointVectorField& nodalMomentum,
    pointTensorField& nodalStress,
    pointVectorField& nodalForce
) const
{
    // Get particle properties
    const vector& pos = position();
    label cellI = cell();
    
    if (cellI < 0) return; // Particle outside mesh
    
    // Calculate shape functions and gradients
    scalarField N;
    List<vector> gradN;
    labelList nodeIDs;
    
    shapeFunctions(mesh, N, gradN, nodeIDs);
    
    // Transfer mass and momentum
    forAll(nodeIDs, i)
    {
        label nodeI = nodeIDs[i];
        
        // Mass: m_i += N_i * m_p
        nodalMass[nodeI] += N[i] * mass_;
        
        // Momentum: p_i += N_i * m_p * v_p
        nodalMomentum[nodeI] += N[i] * mass_ * velocity_;
        
        // CRITICAL: Internal force from stress
        // F_internal,i = -V_p * σ_p · ∇N_i
        // This is how stress creates forces on the grid!
        vector forceContrib = -volume_ * (stress_ & gradN[i]);
        nodalForce[nodeI] += forceContrib;
        
        // Optional: Transfer stress for visualization
        nodalStress[nodeI] += N[i] * stress_;
    }
}

void Foam::MaterialPointParcel::transferFromGrid
(
    const fvMesh& mesh,
    const pointVectorField& nodalVelocity,
    const pointVectorField& nodalAcceleration
)
{
    // Get shape functions
    scalarField N;
    List<vector> gradN;
    labelList nodeIDs;
    
    shapeFunctions(mesh, N, gradN, nodeIDs);
    
    if (nodeIDs.empty()) return;
    
    // Update particle velocity
    velocity_ = Zero;
    
    forAll(nodeIDs, i)
    {
        label nodeI = nodeIDs[i];
        
        // Check bounds
        if (nodeI < 0 || nodeI >= nodalVelocity.size()) continue;
        
        velocity_ += N[i] * nodalVelocity[nodeI];
    }
}


void Foam::MaterialPointParcel::updateDeformationGradient
(
    const fvMesh& mesh,
    const pointVectorField& nodalVelocity,
    const scalar dt
)
{
    // Get shape functions
    scalarField N;
    List<vector> gradN;
    labelList nodeIDs;
    
    shapeFunctions(mesh, N, gradN, nodeIDs);
    
    if (nodeIDs.empty()) return;
    
    // Compute velocity gradient L = ∂v/∂x = Σ v_i ⊗ ∇N_i
    tensor L = Zero;
    
    forAll(nodeIDs, i)
    {
        label nodeI = nodeIDs[i];
        
        // Check bounds
        if (nodeI < 0 || nodeI >= nodalVelocity.size()) continue;
        
        // Outer product: v ⊗ ∇N
        L += nodalVelocity[nodeI] * gradN[i];
    }
    
    // SAFETY CHECK 1: Detect extremely large velocity gradients
    scalar L_mag = mag(L);
    if (L_mag * dt > 0.5)
    {
        WarningInFunction
            << "Large velocity gradient detected for particle " << origId() << nl
            << "  |L|*dt = " << L_mag * dt << " (should be << 1)" << nl
            << "  Time step may be too large!" << endl;
    }
    
    // Store old F for potential rollback
    tensor F_old = F_;
    
    // Update deformation gradient: F^(n+1) = (I + dt*L) * F^n
    // This is the correct multiplicative update for large deformation
    F_ = (tensor::I + dt*L) & F_;
    
    // SAFETY CHECK 2: Check for inversion or extreme deformation
    scalar J = det(F_);
    
    if (J <= 0)
    {
        // CRITICAL: Particle inverted!
        WarningInFunction
            << "Particle " << origId() << " INVERTED: J = " << J << nl
            << "  Position: " << position() << nl
            << "  |L|*dt = " << L_mag * dt << nl
            << "  Reverting to previous F" << endl;
        
        // Emergency: revert to old F
        F_ = F_old;
        J = det(F_);
    }
    else if (J < 0.1 || J > 10.0)
    {
        // Warning: approaching trouble
        WarningInFunction
            << "Particle " << origId() << " has extreme J = " << J << nl
            << "  Position: " << position() << endl;
    }
    
    // Update volume from Jacobian
    // J = det(F) = V / V0, so V = V0 * J
    if (J > 0.01 && J < 100.0)
    {
        volume_ = volume0_ * J;
    }
    else
    {
        // Something is very wrong - keep old volume
        WarningInFunction
            << "Invalid J = " << J << " for particle " << origId()
            << ", keeping old volume" << endl;
    }
}

/*
void Foam::MaterialPointParcel::updateDeformationGradient
(
    const fvMesh& mesh,
    const pointVectorField& nodalVelocity,
    const scalar dt
)
{
    // Get shape functions
    scalarField N;
    List<vector> gradN;
    labelList nodeIDs;
    
    shapeFunctions(mesh, N, gradN, nodeIDs);
    
    if (nodeIDs.empty()) return;
    
    // Compute velocity gradient L = ∂v/∂x
    tensor L = Zero;
    
    forAll(nodeIDs, i)
    {
        label nodeI = nodeIDs[i];
        
        // Check bounds
        if (nodeI < 0 || nodeI >= nodalVelocity.size()) continue;
        
        L += nodalVelocity[nodeI] * gradN[i];
    }
    
    // Update deformation gradient: F_new = (I + dt*L) * F_old
    F_ = (I + dt*L) & F_;
    
    // Update Jacobian and volume
    scalar J = det(F_);
    
    // Limit Jacobian to prevent numerical issues
    if (J > 0.01 && J < 100.0)
    {
        volume_ = volume0_ * J;
    }
}

*/

// * * * * * * * Crystal Plasticity Functions * * * * * * * * * * * * //

void Foam::MaterialPointParcel::initializeSlipSystems(const word& crystalType)
{
    if (crystalType == "FCC")
    {
        // FCC has 12 slip systems: {111}<110>
        slipSystems_.setSize(12);
        
        // Slip plane normals {111}
        vector n1(1, 1, 1);   n1 /= mag(n1);
        vector n2(-1, 1, 1);  n2 /= mag(n2);
        vector n3(1, -1, 1);  n3 /= mag(n3);
        vector n4(1, 1, -1);  n4 /= mag(n4);
        
        // Build slip systems
        label sysId = 0;
        
        // Plane (111) - 3 systems
        slipSystems_[sysId].n = n1;
        slipSystems_[sysId].m = vector(1, -1, 0) / sqrt(2.0);
        slipSystems_[sysId++].schmidTensor = slipSystems_[sysId-1].m * slipSystems_[sysId-1].n;
        
        slipSystems_[sysId].n = n1;
        slipSystems_[sysId].m = vector(-1, 0, 1) / sqrt(2.0);
        slipSystems_[sysId++].schmidTensor = slipSystems_[sysId-1].m * slipSystems_[sysId-1].n;
        
        slipSystems_[sysId].n = n1;
        slipSystems_[sysId].m = vector(0, 1, -1) / sqrt(2.0);
        slipSystems_[sysId++].schmidTensor = slipSystems_[sysId-1].m * slipSystems_[sysId-1].n;
        
        // Plane (-111) - 3 systems
        slipSystems_[sysId].n = n2;
        slipSystems_[sysId].m = vector(1, 1, 0) / sqrt(2.0);
        slipSystems_[sysId++].schmidTensor = slipSystems_[sysId-1].m * slipSystems_[sysId-1].n;
        
        slipSystems_[sysId].n = n2;
        slipSystems_[sysId].m = vector(1, 0, 1) / sqrt(2.0);
        slipSystems_[sysId++].schmidTensor = slipSystems_[sysId-1].m * slipSystems_[sysId-1].n;
        
        slipSystems_[sysId].n = n2;
        slipSystems_[sysId].m = vector(0, 1, 1) / sqrt(2.0);
        slipSystems_[sysId++].schmidTensor = slipSystems_[sysId-1].m * slipSystems_[sysId-1].n;
        
        // Plane (1-11) - 3 systems
        slipSystems_[sysId].n = n3;
        slipSystems_[sysId].m = vector(1, 1, 0) / sqrt(2.0);
        slipSystems_[sysId++].schmidTensor = slipSystems_[sysId-1].m * slipSystems_[sysId-1].n;
        
        slipSystems_[sysId].n = n3;
        slipSystems_[sysId].m = vector(-1, 0, 1) / sqrt(2.0);
        slipSystems_[sysId++].schmidTensor = slipSystems_[sysId-1].m * slipSystems_[sysId-1].n;
        
        slipSystems_[sysId].n = n3;
        slipSystems_[sysId].m = vector(0, 1, 1) / sqrt(2.0);
        slipSystems_[sysId++].schmidTensor = slipSystems_[sysId-1].m * slipSystems_[sysId-1].n;
        
        // Plane (11-1) - 3 systems
        slipSystems_[sysId].n = n4;
        slipSystems_[sysId].m = vector(1, -1, 0) / sqrt(2.0);
        slipSystems_[sysId++].schmidTensor = slipSystems_[sysId-1].m * slipSystems_[sysId-1].n;
        
        slipSystems_[sysId].n = n4;
        slipSystems_[sysId].m = vector(1, 0, 1) / sqrt(2.0);
        slipSystems_[sysId++].schmidTensor = slipSystems_[sysId-1].m * slipSystems_[sysId-1].n;
        
        slipSystems_[sysId].n = n4;
        slipSystems_[sysId].m = vector(0, 1, 1) / sqrt(2.0);
        slipSystems_[sysId++].schmidTensor = slipSystems_[sysId-1].m * slipSystems_[sysId-1].n;
        
        // Initialize state variables
        gamma_.setSize(12, 0.0);
        gamma_dot_.setSize(12, 0.0);
        tau_.setSize(12, 0.0);
        g_.setSize(12, 50e6);  // Initial CRSS in Pa
    }
    else
    {
        FatalErrorIn("MaterialPointParcel::initializeSlipSystems")
            << "Crystal type " << crystalType << " not yet implemented" 
            << exit(FatalError);
    }
}

// void Foam::MaterialPointParcel::setCrystalOrientation
// (
//     scalar phi1, 
//     scalar phi2, 
//     scalar phi3
// )
// {
//     phi1_ = phi1;
//     phi2_ = phi2;  
//     phi3_ = phi3;
    
//     // Convert Bunge Euler angles to rotation matrix
//     scalar c1 = Foam::cos(phi1);
//     scalar s1 = Foam::sin(phi1);
//     scalar c2 = Foam::cos(phi2);
//     scalar s2 = Foam::sin(phi2);
//     scalar c3 = Foam::cos(phi3);
//     scalar s3 = Foam::sin(phi3);
    
//     Q0_ = tensor
//     (
//         c1*c3 - s1*s3*c2,  -c1*s3 - s1*c3*c2,  s1*s2,
//         s1*c3 + c1*s3*c2,  -s1*s3 + c1*c3*c2,  -c1*s2,
//         s3*s2,             c3*s2,              c2
//     );
    
//     // Rotate slip systems to crystal orientation
//     forAll(slipSystems_, i)
//     {
//         // Note: Keep reference slip systems and rotate during stress calculation
//         // This allows for proper rotation with deformation
//     }
// }

void Foam::MaterialPointParcel::computeResolvedShearStress()
{
    // Compute Mandel stress in intermediate configuration
    scalar J = det(F_);
    
    // Protect against singular F
    if (J < 0.01 || J > 100.0) return;
    
    tensor M = J * (Fe_.T() & stress_ & Fe_);
    
    // Compute resolved shear stress on each slip system
    forAll(tau_, i)
    {
        // Rotate slip system to current configuration
        tensor rotatedSchmid = Q0_ & slipSystems_[i].schmidTensor & Q0_.T();
        tau_[i] = M && rotatedSchmid;
    }
}

void Foam::MaterialPointParcel::updatePlasticDeformation
(
    const scalar dt,
    const dictionary& materialProps
)
{
    // Power law flow rule parameters
    scalar gamma0_dot = materialProps.lookupOrDefault<scalar>("referenceShearRate", 0.001);
    scalar m = materialProps.lookupOrDefault<scalar>("rateExponent", 20.0);
    
    // Update slip rates
    forAll(gamma_dot_, i)
    {
        if (mag(tau_[i]) > SMALL && g_[i] > SMALL)
        {
            scalar ratio = mag(tau_[i])/g_[i];
            if (ratio > 0.001)  // Threshold to avoid numerical issues
            {
                gamma_dot_[i] = gamma0_dot * sign(tau_[i]) 
                              * Foam::pow(ratio, 1.0/m);
            }
            else
            {
                gamma_dot_[i] = 0.0;
            }
        }
        else
        {
            gamma_dot_[i] = 0.0;
        }
    }
    
    // Compute plastic velocity gradient Lp
    tensor Lp = Zero;
    forAll(gamma_dot_, i)
    {
        // Rotate slip system to current configuration
        tensor rotatedSchmid = Q0_ & slipSystems_[i].schmidTensor & Q0_.T();
        Lp += gamma_dot_[i] * rotatedSchmid;
    }
    
    // Update plastic deformation gradient
    Fp_ = (I + dt*Lp) & Fp_;
    
    // Update elastic deformation gradient
    scalar detFp = det(Fp_);
    if (detFp > 0.01 && detFp < 100.0)
    {
        Fe_ = F_ & inv(Fp_);
    }
    
    // Update plastic strains
    forAll(gamma_, i)
    {
        gamma_[i] += gamma_dot_[i] * dt;
    }
    
    // Update equivalent plastic strain
    scalar gamma_dot_eq = 0.0;
    forAll(gamma_dot_, i)
    {
        gamma_dot_eq += mag(gamma_dot_[i]);
    }
    epsilonP_eq_ += gamma_dot_eq * dt;
}

void Foam::MaterialPointParcel::updateHardening
(
    const scalar dt,
    const dictionary& materialProps
)
{
    // Voce hardening law parameters
    scalar g0 = materialProps.lookupOrDefault<scalar>("initialCRSS", 50e6);
    scalar gs = materialProps.lookupOrDefault<scalar>("saturationStress", 500e6);
    scalar h0 = materialProps.lookupOrDefault<scalar>("hardeningRate", 1e9);
    scalar q = materialProps.lookupOrDefault<scalar>("latentHardeningRatio", 1.4);
    
    // Update hardness for each slip system
    forAll(g_, alpha)
    {
        scalar g_dot = 0.0;
        
        // Self hardening
        scalar h_self = h0 * Foam::pow(max(1.0 - g_[alpha]/gs, 0.0), 2.0);
        g_dot += h_self * mag(gamma_dot_[alpha]);
        
        // Latent hardening from other systems
        forAll(g_, beta)
        {
            if (beta != alpha)
            {
                scalar h_latent = q * h0 * Foam::pow(max(1.0 - g_[beta]/gs, 0.0), 2.0);
                g_dot += h_latent * mag(gamma_dot_[beta]);
            }
        }
        
        g_[alpha] += g_dot * dt;
        g_[alpha] = min(g_[alpha], gs);
        g_[alpha] = max(g_[alpha], g0);
    }
}

void Foam::MaterialPointParcel::computeElasticStress
(
    const dictionary& materialProps
)
{
    scalar E = materialProps.lookupOrDefault<scalar>("youngModulus", 200e9);
    scalar nu = materialProps.lookupOrDefault<scalar>("poissonRatio", 0.3);
    
    // Lame parameters
    scalar mu = E / (2.0 * (1.0 + nu));
    scalar lambda = E * nu / ((1.0 + nu) * (1.0 - 2.0*nu));
    
    // Check for valid Fe
    scalar Je = det(Fe_);
    if (Je < 0.01 || Je > 100.0) return;
    
    // Neo-Hookean model
    tensor be = Fe_ & Fe_.T();
    
    // Kirchhoff stress
    tensor tau_k = mu * dev(be) + lambda * Foam::log(Je) * I;
    
    // Cauchy stress
    stress_ = (1.0/Je) * tau_k;
}

void Foam::MaterialPointParcel::constitutiveUpdate
(
    const scalar dt,
    const dictionary& materialProps
)
{
    // Complete constitutive update
    
    // 1. Compute resolved shear stress
    computeResolvedShearStress();
    
    // 2. Update plastic deformation
    updatePlasticDeformation(dt, materialProps);
    
    // 3. Update hardening
    updateHardening(dt, materialProps);
    
    // 4. Compute elastic stress
    computeElasticStress(materialProps);
}

// * * * * * * * * * * Tracking Functions * * * * * * * * * * * * * * //

bool Foam::MaterialPointParcel::move
(
    materialPointCloud& cloud,
    trackingData& td,
    const scalar dt
)
{
    // Simple particle advection through mesh
    vector displacement = dt * velocity_;
    
    // Track particle through mesh cells
    scalar f = 1.0;
    label nTrackingSteps = 0;
    const label maxSteps = 100;
    
    // Track displacement through mesh
    while (f > SMALL && nTrackingSteps < maxSteps)
    {
        // trackToFace returns fraction of displacement completed
        // and updates particle position and cell
        f *= trackToFace(f*displacement, f);
        nTrackingSteps++;
    }
    
    // Return true to keep particle, false to delete
    // For MPM, we keep all particles in the domain
    return true;
}

// * * * * * * * * * * * IO Functions * * * * * * * * * * * * * * * * //

void Foam::MaterialPointParcel::readFields(materialPointCloud& c)
{
    if (!c.size())
    {
        return;
    }
    
    particle::readFields(c);
    
    // Note: Actual field reading would be implemented here
    // This is a simplified version
}

void Foam::MaterialPointParcel::writeFields(const materialPointCloud& c)
{
    if (!c.size())
    {
        return;
    }
    
    // Write base particle fields
    particle::writeFields(c);
    
    // Write MPM-specific fields using only public accessors
    IOField<scalar> mass
    (
        c.fieldIOobject("mass", IOobject::NO_READ),
        c.size()
    );
    
    IOField<scalar> volume
    (
        c.fieldIOobject("volume", IOobject::NO_READ),
        c.size()
    );
    
    IOField<vector> velocity
    (
        c.fieldIOobject("velocity", IOobject::NO_READ),
        c.size()
    );
    
    IOField<tensor> stress
    (
        c.fieldIOobject("stress", IOobject::NO_READ),
        c.size()
    );
    
    IOField<tensor> F
    (
        c.fieldIOobject("F", IOobject::NO_READ),
        c.size()
    );
    
    IOField<tensor> Fp
    (
        c.fieldIOobject("Fp", IOobject::NO_READ),
        c.size()
    );
    
    IOField<tensor> Fe
    (
        c.fieldIOobject("Fe", IOobject::NO_READ),
        c.size()
    );
    
    IOField<label> materialID
    (
        c.fieldIOobject("materialID", IOobject::NO_READ),
        c.size()
    );
    
    IOField<label> grainID
    (
        c.fieldIOobject("grainID", IOobject::NO_READ),
        c.size()
    );

        // Write Euler angles
    IOField<scalar> phi1
    (
        c.fieldIOobject("phi1", IOobject::NO_READ),
        c.size()
    );
    
    IOField<scalar> phi2
    (
        c.fieldIOobject("phi2", IOobject::NO_READ),
        c.size()
    );
    
    IOField<scalar> phi3
    (
        c.fieldIOobject("phi3", IOobject::NO_READ),
        c.size()
    );
    
    IOField<tensor> Q0
    (
        c.fieldIOobject("Q0", IOobject::NO_READ),
        c.size()
    );
    
    // Fill the fields using ONLY public accessors
    label np = 0;
    forAllConstIter(materialPointCloud, c, iter)
    {
        const MaterialPointParcel& p = iter();
        
        mass[np] = p.mass();
        volume[np] = p.volume();
        velocity[np] = p.velocity();
        stress[np] = p.stress();
        F[np] = p.F();
        Fp[np] = p.Fp();
        Fe[np] = p.Fe();
        materialID[np] = p.materialID();
        grainID[np] = p.grainID();
        phi1[np] = p.phi1();      
        phi2[np] = p.phi2();      
        phi3[np] = p.phi3();      
        Q0[np] = p.Q0();          
        
        np++;
    }
    
    // Write all fields to disk
    mass.write();
    volume.write();
    velocity.write();
    stress.write();
    F.write();
    Fp.write();
    Fe.write();
    materialID.write();
    grainID.write();
    phi1.write();         
    phi2.write();        
    phi3.write();         
    Q0.write();           
    
    Info<< "    Wrote " << np << " particles to lagrangian cloud" << endl;
}

Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const MaterialPointParcel& p
)
{
    os << static_cast<const particle&>(p)
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
       << token::SPACE << p.gamma_
       << token::SPACE << p.gamma_dot_
       << token::SPACE << p.tau_
       << token::SPACE << p.g_
       << token::SPACE << p.epsilonP_eq_
       << token::SPACE << p.workPlastic_
       << token::SPACE << p.temperature_;
       
    return os;
}

// ************************************************************************* //