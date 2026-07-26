Info<< "Creating material point cloud" << endl;

// Create empty cloud
materialPointCloud particles(mesh, "materialPointCloud", false);

// =========================================================================
// CHECK IF RESTARTING (particle data already exists)
// =========================================================================

bool restarting = false;

// Check if particle cloud files exist in current time directory
word cloudName = "materialPointCloud";
fileName cloudPath = runTime.timeName()/cloud::prefix/cloudName;

if (isDir(cloudPath))
{
    fileName positionsFile = cloudPath/"positions";
    if (isFile(positionsFile))
    {
        restarting = true;
        if (Pstream::master())
        {
            Info<< "  Detected existing particle data - RESTARTING" << endl;
            Info<< "  Reading from: " << cloudPath << endl;
        }
    }
}

if (restarting)
{
    // =====================================================================
    // RESTART: Read existing particles (decomposed per processor)
    // =====================================================================
    
    Info<< "  Processor " << Pstream::myProcNo() 
        << ": Reading existing particles..." << endl;
    
    // Each processor reads its own positions
    IOList<vector> positions
    (
        IOobject
        (
            "positions",
            runTime.timeName(),
            cloud::prefix/cloudName,
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );
    
    Info<< "  Processor " << Pstream::myProcNo() 
        << ": Found " << positions.size() << " particles" << endl;
    
    // Create particles from positions
    forAll(positions, i)
    {
        label cellI = mesh.findCell(positions[i]);
        
        if (cellI >= 0)
        {
            particles.addParticle
            (
                new MaterialPointParcel(mesh, positions[i], cellI)
            );
        }
    }
    
    Info<< "  Processor " << Pstream::myProcNo() 
        << ": Created " << particles.size() << " particles" << endl;
    
    // Read all particle fields
    Info<< "  Processor " << Pstream::myProcNo() 
        << ": Reading particle fields..." << endl;
    MaterialPointParcel::readFields(particles);
    
    
    // Get global statistics
    label nParticlesGlobal = returnReduce(particles.size(), sumOp<label>());
    
    if (Pstream::master())
    {
        Info<< "  Total particles loaded globally: " << nParticlesGlobal << endl;
        Info<< "  Particle cloud successfully restarted" << nl << endl;
    }
}
else
{
    // =====================================================================
    // FRESH START: Initialize particles from scratch
    // =====================================================================
    
    if (Pstream::master())
    {
        Info<< "  No existing particle data - FRESH START" << endl;
    }
    
    // =====================================================================
    // Read grain properties (ALL PROCESSORS - same data)
    // =====================================================================
    
    // Read grain properties - OpenFOAM will automatically find this in
    // the decomposed processor*/constant/ directories after decomposePar
    IOdictionary grainPropertiesDict
    (
        IOobject
        (
            "grainProperties",
            runTime.constant(),
            mesh,  // Use mesh as registry - handles parallel correctly
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );
    
    label nGrains = readLabel(grainPropertiesDict.lookup("nGrains"));
    
    if (Pstream::master())
    {
        Info<< "  Reading " << nGrains << " grain definitions..." << endl;
    }
    
    // Read grain data
    List<dictionary> grainDicts(grainPropertiesDict.lookup("grains"));
    
    struct GrainData
    {
        label grainID;
        label materialID;
        vector seedPosition;
        scalar phi1, Phi, phi2;
    };
    
    List<GrainData> grains(nGrains);
    
    forAll(grainDicts, grainI)
    {
        const dictionary& grainDict = grainDicts[grainI];
        
        grains[grainI].grainID = readLabel(grainDict.lookup("grainID"));
        grains[grainI].materialID = readLabel(grainDict.lookup("materialID"));
        grains[grainI].seedPosition = vector(grainDict.lookup("seedPosition"));
        grains[grainI].phi1 = readScalar(grainDict.lookup("phi1"));
        grains[grainI].Phi = readScalar(grainDict.lookup("Phi"));
        grains[grainI].phi2 = readScalar(grainDict.lookup("phi2"));
        
        if (Pstream::master())
        {
            Info<< "    Grain " << grains[grainI].grainID
                << ": Material " << grains[grainI].materialID
                << ", Seed at " << grains[grainI].seedPosition
                << ", Orientation (" 
                << grains[grainI].phi1 * 180.0 / constant::mathematical::pi << ", "
                << grains[grainI].Phi * 180.0 / constant::mathematical::pi << ", "
                << grains[grainI].phi2 * 180.0 / constant::mathematical::pi << ") deg"
                << endl;
        }
    }
    
    // =====================================================================
    // PARALLEL PARTICLE INITIALIZATION - EFFICIENT MPI APPROACH
    // =====================================================================
    
    List<vector> particlePositions;
    
    if (Pstream::master())
    {
        Info<< "  Master processor reading particle positions from case root..." << endl;
        
        // Read using IOList which properly handles OpenFOAM file format
        // Use explicit path to case root to avoid reading from processor directories
        fileName constantDir = runTime.globalPath()/"constant";
        
        IOList<vector> positions
        (
            IOobject
            (
                "particlePositions",
                constantDir,           // Read from case root constant/
                runTime,              // Use runTime as registry
                IOobject::MUST_READ,
                IOobject::NO_WRITE,
                false                 // Don't register - avoid conflicts
            )
        );
        
        particlePositions = positions;
        
        Info<< "  Master read " << particlePositions.size() 
            << " particle positions" << endl;
    }
    
    // Broadcast number of particles to all processors
    label nParticlesTotal = particlePositions.size();
    Pstream::broadcast(nParticlesTotal);
    
    if (!Pstream::master())
    {
        particlePositions.setSize(nParticlesTotal);
    }
    
    // Broadcast positions to all processors
    Pstream::broadcast(particlePositions);
    
    Info<< "  Processor " << Pstream::myProcNo() 
        << " received " << nParticlesTotal << " particle positions" << endl;
    
    // =====================================================================
    // Create particles on each processor (MPI-efficient)
    // =====================================================================
    
    label particlesCreated = 0;
    
    forAll(particlePositions, i)
    {
        const vector& position = particlePositions[i];
        
        // Find which cell this position is in ON THIS PROCESSOR
        label cellI = mesh.findCell(position);
        
        if (cellI >= 0)
        {
            // This particle belongs to this processor
            particles.addParticle
            (
                new MaterialPointParcel(mesh, position, cellI)
            );
            particlesCreated++;
        }
    }
    
    Info<< "  Processor " << Pstream::myProcNo() 
        << ": Created " << particlesCreated << " particles" << endl;
    
    // MPI reduction to get global statistics
    label nParticlesGlobal = returnReduce(particlesCreated, sumOp<label>());
    
    if (Pstream::master())
    {
        Info<< "  Total particles created globally: " << nParticlesGlobal << endl;
        
        if (nParticlesGlobal < nParticlesTotal)
        {
            WarningInFunction
                << (nParticlesTotal - nParticlesGlobal) 
                << " particles were not assigned to any processor" << endl;
        }
    }
    
    // =====================================================================
    // Assign grain IDs using Voronoi tessellation
    // =====================================================================
    
    Info<< "  Processor " << Pstream::myProcNo() 
        << ": Assigning grain IDs..." << endl;
    
    forAllIter(materialPointCloud, particles, iter)
    {
        MaterialPointParcel& p = iter();
        vector pos = p.position();
        
        // Find nearest grain seed (Voronoi assignment)
        scalar minDist = GREAT;
        label nearestGrain = 0;
        
        forAll(grains, grainI)
        {
            scalar dist = mag(pos - grains[grainI].seedPosition);
            if (dist < minDist)
            {
                minDist = dist;
                nearestGrain = grainI;
            }
        }
        
        // Assign grain properties
        const GrainData& grain = grains[nearestGrain];
        
        p.grainID() = grain.grainID;
        p.materialID() = grain.materialID;
        p.phi1() = grain.phi1;
        p.phi2() = grain.Phi;
        p.phi3() = grain.phi2;
        
        p.setCrystalOrientation(grain.phi1, grain.Phi, grain.phi2);
    }
    
    // =====================================================================
    // Initialize particle properties
    // =====================================================================
    
    if (Pstream::master())
    {
        Info<< "  Initializing particle properties..." << endl;
    }
    
    // Calculate particle volume based on GLOBAL domain volume
    scalar domainVolumeLocal = gSum(mesh.V());
    scalar domainVolumeGlobal = returnReduce(domainVolumeLocal, sumOp<scalar>());
    scalar particleVolume = domainVolumeGlobal / nParticlesGlobal;
    
    if (Pstream::master())
    {
        Info<< "    Domain volume (global): " << domainVolumeGlobal << " m^3" << endl;
        Info<< "    Particle volume: " << particleVolume << " m^3" << endl;
    }
    
    forAllIter(materialPointCloud, particles, iter)
    {
        MaterialPointParcel& p = iter();
        
        p.active() = true;

        p.volume() = particleVolume;
        p.volume0() = particleVolume;
        
        const Material& material = materials[p.materialID()];
        
        scalar rho = material.rho_0;
        p.mass() = rho * particleVolume;
        
        // Initialize slip systems if material has plasticity
        if (material.N_Slips > 0)
        {
            word crystalType = material.crystalStructure;
            p.initializeSlipSystems(crystalType);
            
            p.g().setSize(material.N_Slips, material.initial_CRSS);
            p.gamma().setSize(material.N_Slips, 0.0);
            p.gamma_dot().setSize(material.N_Slips, 0.0);
            p.tau().setSize(material.N_Slips, 0.0);
        }
    }
    
    // =====================================================================
    // Verification
    // =====================================================================
    
    label nWithSlipSystems = 0;
    forAllConstIter(materialPointCloud, particles, iter)
    {
        if (iter().gamma().size() > 0)
        {
            nWithSlipSystems++;
        }
    }
    
    label nWithSlipSystemsGlobal = returnReduce(nWithSlipSystems, sumOp<label>());
    
    if (Pstream::master())
    {
        Info<< "  Particles with slip systems (global): " 
            << nWithSlipSystemsGlobal << " / " << nParticlesGlobal << endl;
        
        if (nWithSlipSystemsGlobal == 0 && materials[0].N_Slips > 0)
        {
            WarningInFunction
                << "No particles have slip systems initialized!" << nl
                << "Crystal plasticity will not work." << endl;
        }
    }
}

if (Pstream::master())
{
    Info<< "\nMaterial point cloud initialization complete" << nl << endl;
}

// ************************************************************************* //