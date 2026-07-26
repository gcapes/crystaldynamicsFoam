/*---------------------------------------------------------------------------*\
  readGrainProperties.H - Read grain definitions (CORRECTED VERSION)
\*---------------------------------------------------------------------------*/

Info<< "Reading grain properties..." << endl;

IOdictionary grainPropertiesDict
(
    IOobject
    (
        "grainProperties",
        runTime.constant(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    )
);

// Storage for grain information
struct GrainInfo
{
    label grainID;
    label materialID;
    vector minPoint;
    vector maxPoint;
    scalar phi1;
    scalar Phi;
    scalar phi2;
    
    GrainInfo()
    :
        grainID(0), materialID(0),
        minPoint(vector::zero), maxPoint(vector::zero),
        phi1(0.0), Phi(0.0), phi2(0.0)
    {}
};

PtrList<GrainInfo> grainDefinitions;

if (grainPropertiesDict.headerOk())
{
    Info<< "  Found grainProperties file" << endl;
    
    word assignmentMethod(grainPropertiesDict.lookup("assignmentMethod"));
    Info<< "  Assignment method: " << assignmentMethod << endl;
    
    if (assignmentMethod == "regions")
    {
        // Read grains subdictionary (not a list!)
        const dictionary& grainsDict = grainPropertiesDict.subDict("grains");
        
        // Count number of grains
        label nGrains = grainsDict.size();
        grainDefinitions.setSize(nGrains);
        
        Info<< "  Number of grains: " << nGrains << endl;
        
        // Iterate over all entries in grains dictionary
        label grainIndex = 0;
        forAllConstIter(dictionary, grainsDict, iter)
        {
            const dictionary& grainDict = iter().dict();
            
            GrainInfo* grain = new GrainInfo();
            
            // Read required entries
            grain->grainID = readLabel(grainDict.lookup("grainID"));
            grain->materialID = readLabel(grainDict.lookup("materialID"));
            
            // Read region subdictionary
            const dictionary& regionDict = grainDict.subDict("region");
            grain->minPoint = vector(regionDict.lookup("min"));
            grain->maxPoint = vector(regionDict.lookup("max"));
            
            // Read orientation subdictionary
            const dictionary& orientDict = grainDict.subDict("orientation");
            grain->phi1 = readScalar(orientDict.lookup("phi1"));
            grain->Phi = readScalar(orientDict.lookup("Phi"));
            grain->phi2 = readScalar(orientDict.lookup("phi2"));
            
            grainDefinitions.set(grainIndex, grain);
            
            Info<< "  Grain " << grain->grainID 
                << ": Material " << grain->materialID
                << ", Orient(" << grain->phi1 << ", " 
                << grain->Phi << ", " << grain->phi2 << ")" << endl;
            
            grainIndex++;
        }
        
        Info<< "  Total grains loaded: " << grainDefinitions.size() << endl;
    }
}
else
{
    Info<< "  No grainProperties file - using default single grain" << endl;
    
    grainDefinitions.setSize(1);
    GrainInfo* defaultGrain = new GrainInfo();
    defaultGrain->grainID = 0;
    defaultGrain->materialID = 0;
    defaultGrain->minPoint = vector(-GREAT, -GREAT, -GREAT);
    defaultGrain->maxPoint = vector(GREAT, GREAT, GREAT);
    defaultGrain->phi1 = 0.0;
    defaultGrain->Phi = 0.0;
    defaultGrain->phi2 = 0.0;
    grainDefinitions.set(0, defaultGrain);
}

Info<< endl;

// ************************************************************************* //