/*---------------------------------------------------------------------------*\
  createMPMCloud.H - Initialize material point cloud with grain assignments
\*---------------------------------------------------------------------------*/

Info<< "Creating material point cloud" << endl;

// Read particle positions using IOList (no IFstream needed)
List<vector> particlePositions
(
    IOList<vector>
    (
        IOobject
        (
            "particlePositions",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    )
);

label nParticles = particlePositions.size();
Info<< "  Number of particles to create: " << nParticles << endl;

// Create the cloud
materialPointCloud particles(mesh, "materialPointCloud");

// Read positions and create particles
label particlesCreated = 0;
label particlesOutsideMesh = 0;

forAll(particlePositions, i)
{
    const vector& position = particlePositions[i];
    
    // Find which cell this position is in
    label cellI = mesh.findCell(position);
    
    if (cellI >= 0)
    {
        // Create particle
        particles.addParticle
        (
            new MaterialPointParcel(mesh, position, cellI)
        );
        particlesCreated++;
    }
    else
    {
        particlesOutsideMesh++;
    }
}

Info<< "  Particles created: " << particlesCreated << endl;
if (particlesOutsideMesh > 0)
{
    WarningIn("createMPMCloud")
        << particlesOutsideMesh << " particles were outside the mesh and not created"
        << endl;
}

// Now assign grain IDs and properties based on spatial location
Info<< "  Assigning grain IDs to particles..." << endl;

label particlesAssigned = 0;
label particlesUnassigned = 0;

forAllIter(materialPointCloud, particles, iter)
{
    MaterialPointParcel& p = iter();
    vector pos = p.position();
    
    // Find which grain this particle belongs to based on spatial region
    bool assigned = false;
    
    forAll(grainDefinitions, grainI)
    {
        const GrainInfo& grain = grainDefinitions[grainI];
        
        // Check if particle is within grain's bounding box
        if (pos.x() >= grain.minPoint.x() && pos.x() <= grain.maxPoint.x() &&
            pos.y() >= grain.minPoint.y() && pos.y() <= grain.maxPoint.y() &&
            pos.z() >= grain.minPoint.z() && pos.z() <= grain.maxPoint.z())
        {
            // Assign grain properties to particle
            p.grainID() = grain.grainID;
            p.materialID() = grain.materialID;
            
            // Set crystal orientation
            p.phi1() = grain.phi1;
            p.phi2() = grain.Phi;  // Note: readGrainProperties uses Phi
            p.phi3() = grain.phi2;
            
            // Initialize rotation matrix from Euler angles
            p.setCrystalOrientation(grain.phi1, grain.Phi, grain.phi2);
            
            assigned = true;
            particlesAssigned++;
            break;
        }
    }
    
    if (!assigned)
    {
        // Particle not in any grain region - assign to default grain 0
        p.grainID() = 0;
        p.materialID() = 0;
        p.phi1() = 0.0;
        p.phi2() = 0.0;
        p.phi3() = 0.0;
        p.setCrystalOrientation(0.0, 0.0, 0.0);
        particlesUnassigned++;
    }
}

Info<< "  Particles assigned to grains: " << particlesAssigned << endl;
if (particlesUnassigned > 0)
{
    WarningIn("createMPMCloud")
        << particlesUnassigned << " particles not in any grain region (assigned to default grain 0)"
        << endl;
}

// Initialize particle volumes and masses
Info<< "  Initializing particle properties..." << endl;

scalar domainVolume = gSum(mesh.V());
scalar particleVolume = domainVolume / particles.size();

Info<< "    Domain volume: " << domainVolume << " m^3" << endl;
Info<< "    Particle volume: " << particleVolume << " m^3" << endl;

forAllIter(materialPointCloud, particles, iter)
{
    MaterialPointParcel& p = iter();
    
    // Set volume
    p.volume() = particleVolume;
    p.volume0() = particleVolume;
    
    // Get material properties for this particle's material
    const dictionary& matProps = materialProperties[p.materialID()];
    scalar rho = readScalar(matProps.lookup("rho"));
    
    // Set mass
    p.mass() = rho * particleVolume;
    
    // Initialize slip systems based on crystal structure
    word crystalType = word(matProps.lookup("crystalStructure"));
    p.initializeSlipSystems(crystalType);
}

// Print grain statistics
Info<< "\nGrain distribution:" << endl;
labelList particlesPerGrain(grainDefinitions.size(), 0);

forAllConstIter(materialPointCloud, particles, iter)
{
    label gid = iter().grainID();
    if (gid >= 0 && gid < grainDefinitions.size())
    {
        particlesPerGrain[gid]++;
    }
}

forAll(grainDefinitions, grainI)
{
    const GrainInfo& grain = grainDefinitions[grainI];
    Info<< "  Grain " << grain.grainID 
        << " (Material " << grain.materialID << "): "
        << particlesPerGrain[grainI] << " particles"
        << ", Orientation(" << grain.phi1 << ", " << grain.Phi << ", " << grain.phi2 << ")" 
        << endl;
}

Info<< "\nMaterial point cloud created successfully" << nl << endl;

// ************************************************************************* //