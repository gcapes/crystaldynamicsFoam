// Initialize material point particles

IOdictionary particleDict
(
    IOobject
    (
        "particleProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
);

// Check if we should read from file or initialize
word initMethod(particleDict.lookup("initializationMethod"));

if (initMethod == "fromFile")
{
    // Particles will be read automatically by Cloud constructor
    Info<< "Reading particles from " << runTime.timeName() << endl;
}
else if (initMethod == "uniform")
{
    // Uniform distribution in domain
    scalar particlesPerCell = readScalar(particleDict.lookup("particlesPerCell"));
    
    Info<< "Initializing " << particlesPerCell << " particles per cell" << endl;
    
    // Loop over cells
    forAll(mesh.C(), celli)
    {
        const vector& cellCenter = mesh.C()[celli];
        scalar cellVol = mesh.V()[celli];
        scalar dx = Foam::pow(cellVol / particlesPerCell, 1.0/3.0);
        
        // Create grid of particles within cell
        label nPerSide = label(Foam::pow(particlesPerCell, 1.0/3.0) + 0.5);
        
        for (label i = 0; i < nPerSide; i++)
        {
            for (label j = 0; j < nPerSide; j++)
            {
                for (label k = 0; k < nPerSide; k++)
                {
                    vector offset
                    (
                        (i - nPerSide/2.0) * dx,
                        (j - nPerSide/2.0) * dx,
                        (k - nPerSide/2.0) * dx
                    );
                    
                    vector particlePos = cellCenter + offset;
                    
                    // Create particle
                    MaterialPointParcel* pPtr = new MaterialPointParcel
                    (
                        mesh,
                        particlePos,
                        celli
                    );
                    
                    // Initialize particle properties
                    pPtr->mass_ = cellVol * 7850.0 / particlesPerCell; // Density kg/m^3
                    pPtr->volume_ = cellVol / particlesPerCell;
                    pPtr->velocity_ = vector::zero;
                    pPtr->stress_ = tensor::zero;
                    pPtr->plasticStrain_ = tensor::zero;
                    pPtr->F_ = tensor::I;
                    pPtr->materialID_ = 0; // Will be set by grain structure
                    pPtr->phaseID_ = 0;
                    
                    // Initialize crystal plasticity variables
                    label nSlips = 12; // FCC slip systems
                    pPtr->g_alphas_.setSize(nSlips, 50e6); // Initial CRSS (Pa)
                    pPtr->gamma_alphas_.setSize(nSlips, 0.0); // No initial slip
                    
                    // Euler angles will be set by grain initialization
                    pPtr->phi1_ = 0.0;
                    pPtr->phi2_ = 0.0;
                    pPtr->phi3_ = 0.0;
                    
                    particles.addParticle(pPtr);
                }
            }
        }
    }
    
    // Assign grain structure
    #include "assignGrainStructure.H"
}

Info<< "Initialized " << particles.size() << " particles" << endl;
