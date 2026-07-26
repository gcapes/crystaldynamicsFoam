/*---------------------------------------------------------------------------*\
  readMaterialProperties.H - UPDATED FOR IMPROVED MATERIAL CLASS
\*---------------------------------------------------------------------------*/

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

// Store Material objects (not dictionaries!) for each material ID
PtrList<Material> materials;

// Get number of materials
label nMaterials = materialPropertiesDict.lookupOrDefault<label>("nMaterials", 1);
materials.setSize(nMaterials);

Info<< "  Number of materials: " << nMaterials << nl << endl;

// Read and initialize each material
for (label matI = 0; matI < nMaterials; matI++)
{
    word materialName = "material" + Foam::name(matI);
    
    if (materialPropertiesDict.found(materialName))
    {
        const dictionary& matDict = materialPropertiesDict.subDict(materialName);
        
        // Create Material object
        Material* material = new Material();
        
        // Read basic properties
        material->Mat_type = matDict.lookupOrDefault<word>("Mat_type", "generic");
        material->rho_0 = matDict.lookupOrDefault<scalar>("rho", 1000.0);
        material->Cv = matDict.lookupOrDefault<scalar>("Cv", 1000.0);
        material->kappa = matDict.lookupOrDefault<scalar>("kappa", 1.0);
        
        Info<< "Material " << matI << ": " << materialName << endl;
        Info<< "  Type: " << material->Mat_type << endl;
        Info<< "  Density: " << material->rho_0 << " kg/m^3" << endl;
        
        // Initialize elasticity (handles all symmetries: isotropic, cubic, hexagonal, etc.)
        material->initializeElasticity(matDict);
        
        // Read crystal structure for slip system initialization
        if (matDict.found("crystalStructure"))
        {
            material->crystalStructure = word(matDict.lookup("crystalStructure"));
            Info<< "  Crystal structure: " << material->crystalStructure << endl;
        }
        else
        {
            material->crystalStructure = "FCC";  // Default
        }
        
        // Initialize plasticity if specified
        if (matDict.found("N_Slips"))
        {
            material->initializePlasticity(matDict);
        }
        
        // Print material properties
        material->Print(0);
        
        materials.set(matI, material);
        
        Info<< endl;
    }
    else
    {
        FatalErrorInFunction
            << "Material " << materialName << " not found in materialProperties"
            << exit(FatalError);
    }
}

Info<< "Material properties initialized for " << nMaterials << " materials" << nl << endl;

// ************************************************************************* //