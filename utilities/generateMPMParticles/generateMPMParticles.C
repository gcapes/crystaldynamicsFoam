/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | 
    \\  /    A nd           | 
     \\/     M anipulation  |
-------------------------------------------------------------------------------
Application
    generateMPMParticles

Description
    Generates material point particles with Voronoi tessellation for grain
    assignment. Reads domain bounds from blockMeshDict and generates:
    
    1. N particles per cell (user specified)
    2. M Voronoi grains (user specified)
    3. Assigns each particle to nearest Voronoi seed
    4. Assigns random crystal orientations (Euler angles) per grain
    5. Assigns material types per grain
    
    Output files:
    - constant/particlePositions (List<vector>)
    - constant/grainProperties (IOdictionary)

Author
    Tom Flint

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "IFstream.H"
#include "OFstream.H"
#include "Random.H"
#include "boundBox.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

struct VoronoiSeed
{
    label grainID;
    label materialID;
    vector position;
    scalar phi1;  // Euler angle (radians)
    scalar Phi;   // Euler angle (radians)  
    scalar phi2;  // Euler angle (radians)
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Generate MPM particles with Voronoi grain assignment"
    );
    
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createMesh.H"
    
    // =========================================================================
    // Read generation parameters
    // =========================================================================
    
    IOdictionary particleGenDict
    (
        IOobject
        (
            "particleGenerationProperties",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );
    
    // Number of particles per cell
    label particlesPerCell = readLabel(particleGenDict.lookup("particlesPerCell"));
    
    // Number of Voronoi grains
    label nGrains = readLabel(particleGenDict.lookup("nGrains"));
    
    // Material assignment for each grain
    // Can be:
    // 1. Single material for all grains: materialIDs (0);
    // 2. Random materials: materialIDs uniform;
    // 3. Specific list: materialIDs (0 1 0 1 ...);
    word materialAssignment = particleGenDict.lookupOrDefault<word>("materialAssignment", "single");
    label defaultMaterialID = particleGenDict.lookupOrDefault<label>("defaultMaterialID", 0);
    label nMaterialTypes = particleGenDict.lookupOrDefault<label>("nMaterialTypes", 1);
    
    // Orientation assignment
    word orientationMode = particleGenDict.lookupOrDefault<word>("orientationMode", "random");
    
    // Optional: specify material IDs explicitly
    labelList materialIDs;
    if (particleGenDict.found("materialIDs"))
    {
        materialIDs = labelList(particleGenDict.lookup("materialIDs"));
        if (materialIDs.size() != nGrains)
        {
            FatalErrorInFunction
                << "materialIDs list size (" << materialIDs.size() 
                << ") does not match nGrains (" << nGrains << ")"
                << exit(FatalError);
        }
    }
    
    // Optional: Domain bounds override (otherwise read from blockMeshDict)
    boundBox domain;
    bool useMeshBounds = particleGenDict.lookupOrDefault<bool>("useMeshBounds", true);
    
    if (useMeshBounds)
    {
        // Get domain bounds from mesh
        domain = boundBox(mesh.points(), false);
        Info<< "Using mesh bounds: " << domain << endl;
    }
    else
    {
        // Read explicit bounds
        vector minPoint(particleGenDict.lookup("minPoint"));
        vector maxPoint(particleGenDict.lookup("maxPoint"));
        domain = boundBox(minPoint, maxPoint);
        Info<< "Using specified bounds: " << domain << endl;
    }
    
    Info<< "\nGeneration parameters:" << nl
        << "  Particles per cell: " << particlesPerCell << nl
        << "  Number of grains: " << nGrains << nl
        << "  Material assignment: " << materialAssignment << nl
        << "  Orientation mode: " << orientationMode << nl
        << "  Domain: " << domain << nl
        << endl;
    
    // =========================================================================
    // Generate Voronoi seeds (grain centers)
    // =========================================================================
    
    Info<< "Generating " << nGrains << " Voronoi grain seeds..." << endl;
    
    Random rndGen(123456);  // Fixed seed for reproducibility
    
    List<VoronoiSeed> grainSeeds(nGrains);
    
    forAll(grainSeeds, grainI)
    {
        VoronoiSeed& seed = grainSeeds[grainI];
        
        seed.grainID = grainI;
        
        // Random position within domain
        seed.position.x() = domain.min().x() + rndGen.sample01<scalar>() * domain.span().x();
        seed.position.y() = domain.min().y() + rndGen.sample01<scalar>() * domain.span().y();
        seed.position.z() = domain.min().z() + rndGen.sample01<scalar>() * domain.span().z();
        
        // Assign material ID
        if (materialAssignment == "single")
        {
            seed.materialID = defaultMaterialID;
        }
        else if (materialAssignment == "random")
        {
            seed.materialID = label(rndGen.sample01<scalar>() * nMaterialTypes);
            if (seed.materialID >= nMaterialTypes) seed.materialID = nMaterialTypes - 1;
        }
        else if (materialAssignment == "specified" && materialIDs.size() > 0)
        {
            seed.materialID = materialIDs[grainI];
        }
        else
        {
            seed.materialID = grainI % nMaterialTypes;
        }
        
        // Assign crystal orientation (Euler angles in radians)
        if (orientationMode == "random")
        {
            // Random orientation using ZYZ Euler angles
            seed.phi1 = rndGen.sample01<scalar>() * 2.0 * constant::mathematical::pi;  // [0, 2π]
            seed.Phi = Foam::acos(2.0 * rndGen.sample01<scalar>() - 1.0);  // [0, π] uniform on sphere
            seed.phi2 = rndGen.sample01<scalar>() * 2.0 * constant::mathematical::pi;  // [0, 2π]
        }
        else if (orientationMode == "randomCubic")
        {
            // Random orientation restricted to cubic fundamental zone
            seed.phi1 = rndGen.sample01<scalar>() * constant::mathematical::pi / 2.0;  // [0, π/2]
            seed.Phi = Foam::acos(2.0 * rndGen.sample01<scalar>() - 1.0);
            seed.phi2 = rndGen.sample01<scalar>() * constant::mathematical::pi / 2.0;  // [0, π/2]
        }
        else if (orientationMode == "aligned")
        {
            // All grains have same orientation (useful for testing)
            seed.phi1 = 0.0;
            seed.Phi = 0.0;
            seed.phi2 = 0.0;
        }
        else
        {
            FatalErrorInFunction
                << "Unknown orientationMode: " << orientationMode << nl
                << "Valid options: random, randomCubic, aligned"
                << exit(FatalError);
        }
        
        Info<< "  Grain " << grainI 
            << ": pos = " << seed.position
            << ", material = " << seed.materialID
            << ", orientation (deg) = (" 
            << seed.phi1 * 180.0 / constant::mathematical::pi << ", "
            << seed.Phi * 180.0 / constant::mathematical::pi << ", "
            << seed.phi2 * 180.0 / constant::mathematical::pi << ")"
            << endl;
    }
    
    // =========================================================================
    // Generate particles in each cell
    // =========================================================================
    
    Info<< "\nGenerating particles in mesh cells..." << endl;
    
    label nCells = mesh.nCells();
    label totalParticles = nCells * particlesPerCell;
    
    Info<< "  Total particles to generate: " << totalParticles << endl;
    
    List<vector> particlePositions(totalParticles);
    labelList particleGrainIDs(totalParticles);
    
    label pIndex = 0;
    
    forAll(mesh.C(), cellI)
    {
        vector cellCenter = mesh.C()[cellI];
        scalar cellSize = Foam::pow(mesh.V()[cellI], 1.0/3.0);
        
        // Generate particlesPerCell particles in this cell
        for (label p = 0; p < particlesPerCell; p++)
        {
            // Random offset within cell (using simple cubic distribution)
            scalar dx = (rndGen.sample01<scalar>() - 0.5) * cellSize;
            scalar dy = (rndGen.sample01<scalar>() - 0.5) * cellSize;
            scalar dz = (rndGen.sample01<scalar>() - 0.5) * cellSize;
            
            vector particlePos = cellCenter + vector(dx, dy, dz);
            
            // Find nearest Voronoi seed
            scalar minDist = GREAT;
            label nearestGrain = 0;
            
            forAll(grainSeeds, grainI)
            {
                scalar dist = mag(particlePos - grainSeeds[grainI].position);
                if (dist < minDist)
                {
                    minDist = dist;
                    nearestGrain = grainI;
                }
            }
            
            particlePositions[pIndex] = particlePos;
            particleGrainIDs[pIndex] = nearestGrain;
            
            pIndex++;
        }
        
        if (cellI % 1000 == 0)
        {
            Info<< "  Processed " << cellI << " / " << nCells << " cells" << endl;
        }
    }
    
    Info<< "  Generated " << pIndex << " particles" << endl;
    
    // =========================================================================
    // Print grain statistics
    // =========================================================================
    
    Info<< "\nGrain statistics:" << endl;
    
    labelList particlesPerGrain(nGrains, 0);
    
    forAll(particleGrainIDs, i)
    {
        particlesPerGrain[particleGrainIDs[i]]++;
    }
    
    forAll(particlesPerGrain, grainI)
    {
        scalar percentage = 100.0 * scalar(particlesPerGrain[grainI]) / scalar(totalParticles);
        Info<< "  Grain " << grainI 
            << " (Material " << grainSeeds[grainI].materialID << "): "
            << particlesPerGrain[grainI] << " particles ("
            << percentage << " %)" << endl;
    }
    
    // =========================================================================
    // Write particle positions
    // =========================================================================
    
    Info<< "\nWriting particle positions to constant/particlePositions..." << endl;
    
    IOList<vector> particlePositionsIO
    (
        IOobject
        (
            "particlePositions",
            runTime.constant(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        particlePositions
    );
    
    particlePositionsIO.write();
    
    // =========================================================================
    // Write grain properties
    // =========================================================================
    
    Info<< "Writing grain properties to constant/grainProperties..." << endl;
    
    OFstream grainPropsFile
    (
        runTime.constant()/"grainProperties"
    );
    
    grainPropsFile
        << "/*--------------------------------*- C++ -*----------------------------------*\\" << nl
        << "| =========                 |                                                 |" << nl
        << "| \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |" << nl
        << "|  \\\\    /   O peration     | Version:  v2012                                 |" << nl
        << "|   \\\\  /    A nd           | Website:  www.openfoam.com                      |" << nl
        << "|    \\\\/     M anipulation  |                                                 |" << nl
        << "\\*---------------------------------------------------------------------------*/" << nl
        << "FoamFile" << nl
        << "{" << nl
        << "    version     2.0;" << nl
        << "    format      ascii;" << nl
        << "    class       dictionary;" << nl
        << "    object      grainProperties;" << nl
        << "}" << nl
        << "// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //" << nl
        << nl;
    
    grainPropsFile << "nGrains " << nGrains << ";" << nl << nl;
    
    grainPropsFile << "grains" << nl << "(" << nl;
    
    forAll(grainSeeds, grainI)
    {
        const VoronoiSeed& seed = grainSeeds[grainI];
        
        grainPropsFile
            << "    {" << nl
            << "        grainID      " << seed.grainID << ";" << nl
            << "        materialID   " << seed.materialID << ";" << nl
            << "        seedPosition " << seed.position << ";" << nl
            << "        phi1         " << seed.phi1 << ";  // radians" << nl
            << "        Phi          " << seed.Phi << ";  // radians" << nl
            << "        phi2         " << seed.phi2 << ";  // radians" << nl
            << "        minPoint     " << domain.min() << ";" << nl
            << "        maxPoint     " << domain.max() << ";" << nl
            << "    }" << nl;
    }
    
    grainPropsFile << ");" << nl << nl;
    grainPropsFile << "// ************************************************************************* //" << nl;
    
    Info<< "\nParticle generation complete!" << nl
        << "  Particle positions: constant/particlePositions" << nl
        << "  Grain properties: constant/grainProperties" << nl
        << endl;
    
    // =========================================================================
    // Optional: Write VTK visualization of grain seeds
    // =========================================================================
    
    Info<< "Writing grain seed visualization..." << endl;
    
    fileName vtkFile = runTime.constant()/"grainSeeds.vtk";
    OFstream vtk(vtkFile);
    
    vtk << "# vtk DataFile Version 2.0" << nl
        << "Voronoi Grain Seeds" << nl
        << "ASCII" << nl
        << "DATASET POLYDATA" << nl
        << "POINTS " << nGrains << " float" << nl;
    
    forAll(grainSeeds, grainI)
    {
        vtk << grainSeeds[grainI].position.x() << " "
            << grainSeeds[grainI].position.y() << " "
            << grainSeeds[grainI].position.z() << nl;
    }
    
    vtk << nl << "POINT_DATA " << nGrains << nl
        << "SCALARS materialID int 1" << nl
        << "LOOKUP_TABLE default" << nl;
    
    forAll(grainSeeds, grainI)
    {
        vtk << grainSeeds[grainI].materialID << nl;
    }
    
    Info<< "  Grain seeds visualization: constant/grainSeeds.vtk" << nl << endl;
    
    Info<< "Done!" << endl;
    
    return 0;
}

// ************************************************************************* //