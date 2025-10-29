// Read material properties for crystal plasticity

Info<< "Reading material properties..." << endl;

IOdictionary materialPropertiesDict
(
    IOobject
    (
        "materialProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
);

// Store material properties for each material ID
PtrList<dictionary> materialProperties;

// Get number of materials
label nMaterials = materialPropertiesDict.lookupOrDefault<label>("nMaterials", 1);
materialProperties.setSize(nMaterials);

// Read each material's properties
for (label matI = 0; matI < nMaterials; matI++)
{
    word materialName = "material" + Foam::name(matI);
    
    if (materialPropertiesDict.found(materialName))
    {
        materialProperties.set
        (
            matI,
            new dictionary(materialPropertiesDict.subDict(materialName))
        );
        
        Info<< "Material " << matI << ": " << materialName << endl;
        
        // Read and display density
        scalar rho = materialProperties[matI].lookupOrDefault<scalar>("rho", 1000.0);
        Info<< "  Density: " << rho << " kg/m^3" << endl;
    }
    else
    {
        FatalErrorInFunction
            << "Material " << materialName << " not found in materialProperties"
            << exit(FatalError);
    }
}

Info<< "Material properties read for " << nMaterials << " materials" << endl;