Info<< "Reading material properties\n" << endl;

IOdictionary materialsDict
(
    IOobject
    (
        "materialProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);

// Get list of materials
wordList materialNames = materialsDict.toc();

Info<< "Found " << materialNames.size() << " materials:" << endl;

// Store material properties for each material ID
PtrList<dictionary> materialProperties(materialNames.size());

forAll(materialNames, matI)
{
    Info<< "  Material " << matI << ": " << materialNames[matI] << endl;
    
    materialProperties.set
    (
        matI,
        new dictionary(materialsDict.subDict(materialNames[matI]))
    );
    
    // Display key properties
    const dictionary& matDict = materialProperties[matI];
    
    Info<< "    Crystal type: " 
        << matDict.lookupOrDefault<word>("crystalType", "FCC") << endl;
    
    Info<< "    Young's modulus: " 
        << matDict.lookupOrDefault<scalar>("youngModulus", 200e9)/1e9 << " GPa" << endl;
        
    Info<< "    Initial CRSS: "
        << matDict.lookupOrDefault<scalar>("initialCRSS", 50e6)/1e6 << " MPa" << endl;
        
    Info<< "    Saturation stress: "
        << matDict.lookupOrDefault<scalar>("saturationStress", 500e6)/1e6 << " MPa" << endl;
}

Info<< endl;
