// Assign grain structure to particles

Info<< "Assigning grain structure to particles..." << endl;

// Read grain properties
label nGrains = readLabel(particleDict.lookup("numberOfGrains"));

Info<< "Creating " << nGrains << " grains" << endl;

// Generate random grain centers
Random rndGen(12345); // Fixed seed for reproducibility

List<vector> grainCenters(nGrains);
List<tensor> grainRotations(nGrains);

boundBox bounds(mesh.points());

forAll(grainCenters, graini)
{
    // Random position in domain
    grainCenters[graini] = vector
    (
        rndGen.position<scalar>(bounds.min().x(), bounds.max().x()),
        rndGen.position<scalar>(bounds.min().y(), bounds.max().y()),
        rndGen.position<scalar>(bounds.min().z(), bounds.max().z())
    );
    
    // Random orientation (Euler angles)
    scalar phi1 = rndGen.position<scalar>(0.0, 2.0 * constant::mathematical::pi);
    scalar cosTheta = rndGen.position<scalar>(-1.0, 1.0);
    scalar Phi = Foam::acos(cosTheta);
    scalar phi2 = rndGen.position<scalar>(0.0, 2.0 * constant::mathematical::pi);
    
    // Build rotation matrix from Euler angles (Bunge convention)
    scalar c1 = Foam::cos(phi1);
    scalar s1 = Foam::sin(phi1);
    scalar c  = Foam::cos(Phi);
    scalar s  = Foam::sin(Phi);
    scalar c2 = Foam::cos(phi2);
    scalar s2 = Foam::sin(phi2);
    
    grainRotations[graini] = tensor
    (
        c1*c2 - s1*s2*c,  s1*c2 + c1*s2*c, s2*s,
       -c1*s2 - s1*c2*c, -s1*s2 + c1*c2*c, c2*s,
        s1*s,            -c1*s,             c
    );
}

// Assign each particle to nearest grain
forAllIter(materialPointCloud, particles, iter)
{
    MaterialPointParcel& p = iter();
    
    // Find nearest grain center
    scalar minDist = GREAT;
    label nearestGrain = 0;
    
    forAll(grainCenters, graini)
    {
        scalar dist = mag(p.position() - grainCenters[graini]);
        if (dist < minDist)
        {
            minDist = dist;
            nearestGrain = graini;
        }
    }
    
    // Assign grain properties
    p.materialID_ = nearestGrain;
    
    // Extract Euler angles from rotation matrix
    const tensor& R = grainRotations[nearestGrain];
    
    p.phi1_ = Foam::atan2(R.zy(), -R.zx());
    p.phi2_ = Foam::acos(R.zz());
    p.phi3_ = Foam::atan2(R.yz(), R.xz());
}

Info<< "Grain structure assigned" << endl;