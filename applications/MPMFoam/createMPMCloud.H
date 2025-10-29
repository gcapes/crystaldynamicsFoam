/*---------------------------------------------------------------------------*\
  createMPMCloud.H - Particle Cloud Initialization
  
  This file creates and initializes the MPM particle cloud with proper
  mass and volume values based on material properties.
\*---------------------------------------------------------------------------*/

Info<< "Creating MPM particle cloud..." << endl;

// Create the cloud
materialPointCloud particles
(
    mesh,
    "mpmCloud",
    IDLList<MaterialPointParcel>()
);

// Read cloud properties
IOdictionary mpmCloudProperties
(
    IOobject
    (
        "mpmCloudProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);

// Get initialization parameters
word initType(mpmCloudProperties.lookup("initializationType"));
label nParticlesPerCell = readLabel(mpmCloudProperties.lookup("nParticlesPerCell"));

Info<< "  Initialization type: " << initType << endl;
Info<< "  Particles per cell: " << nParticlesPerCell << endl;

// Material region for particle generation
const dictionary& regionDict = mpmCloudProperties.subDict("materialRegion");
word regionType(regionDict.lookup("type"));

vector minPoint(regionDict.lookup("min"));
vector maxPoint(regionDict.lookup("max"));

Info<< "  Material region: " << regionType << endl;
Info<< "    Min: " << minPoint << endl;
Info<< "    Max: " << maxPoint << endl;

// Get material properties
label materialID = readLabel(mpmCloudProperties.lookup("materialID"));
const dictionary& matProps = materialProperties[materialID];
scalar rho = readScalar(matProps.lookup("rho"));

Info<< "  Material ID: " << materialID << endl;
Info<< "  Density: " << rho << " kg/m³" << endl;

// Generate particles
label nParticlesGenerated = 0;
scalar totalMassGenerated = 0.0;
scalar totalVolumeGenerated = 0.0;

if (initType == "generate")
{
    Info<< "\nGenerating particles..." << endl;
    
    // Loop over all cells
    forAll(mesh.C(), celli)
    {
        const vector& cellCenter = mesh.C()[celli];
        scalar cellVolume = mesh.V()[celli];
        
        // Check if cell center is in material region
        bool inRegion = false;
        
        if (regionType == "box")
        {
            inRegion = 
                cellCenter.x() >= minPoint.x() && cellCenter.x() <= maxPoint.x() &&
                cellCenter.y() >= minPoint.y() && cellCenter.y() <= maxPoint.y() &&
                cellCenter.z() >= minPoint.z() && cellCenter.z() <= maxPoint.z();
        }
        
        if (inRegion)
        {
            // Compute particle volume and mass
            scalar particleVolume = cellVolume / scalar(nParticlesPerCell);
            scalar particleMass = rho * particleVolume;
            
            // Generate particles in this cell
            for (label pi = 0; pi < nParticlesPerCell; pi++)
            {
                // Particle position (simple grid for multiple particles per cell)
                vector particlePos = cellCenter;
                
                if (nParticlesPerCell > 1)
                {
                    // Offset from cell center for multiple particles
                    scalar dx = cellVolume / nParticlesPerCell;
                    scalar offset = (pi - 0.5*nParticlesPerCell) * 0.3 * dx;
                    particlePos.x() += offset;
                }
                
                // Create particle
                MaterialPointParcel* pPtr = new MaterialPointParcel
                (
                    mesh,
                    particlePos,
                    celli
                );
                
                // CRITICAL: Set mass and volume!
                pPtr->mass() = particleMass;
                pPtr->volume() = particleVolume;
                pPtr->velocity() = Zero;
                pPtr->stress() = Zero;
                pPtr->F() = tensor::I;
                
                // Add to cloud
                particles.addParticle(pPtr);
                
                nParticlesGenerated++;
                totalMassGenerated += particleMass;
                totalVolumeGenerated += particleVolume;
            }
        }
    }
    
    Info<< "  Generated " << nParticlesGenerated << " particles" << endl;
    Info<< "  Total mass: " << totalMassGenerated << " kg" << endl;
    Info<< "  Total volume: " << totalVolumeGenerated << " m³" << endl;
    Info<< "  Average particle mass: " 
        << (nParticlesGenerated > 0 ? totalMassGenerated/nParticlesGenerated : 0) 
        << " kg" << endl;
}
else if (initType == "file")
{
    FatalErrorInFunction
        << "Initialization from file not yet implemented" << nl
        << exit(FatalError);
}
else
{
    FatalErrorInFunction
        << "Unknown initialization type: " << initType << nl
        << "Valid types: generate, file" << nl
        << exit(FatalError);
}

if (nParticlesGenerated == 0)
{
    FatalErrorInFunction
        << "No particles were generated!" << nl
        << "Check material region bounds in mpmCloudProperties" << nl
        << "Region: [" << minPoint << "] to [" << maxPoint << "]" << nl
        << exit(FatalError);
}

Info<< "Particle cloud created successfully\n" << endl;