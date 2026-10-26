/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | 
    \\  /    A nd           | 
     \\/     M anipulation  |
-------------------------------------------------------------------------------
Application
    generateFromDREAM3D

Description
    Reads a DREAM3D (HDF5) microstructure from constant/ and writes a matching system/blockMeshDict,
    constant/GrainProperties (orientations, materials) and, once the mesh exists, the n.* fields.

Author
    Tom Flint

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "IFstream.H"
#include "OFstream.H"
#include "boundBox.H"
#include <hdf5.h>
#include <vector>
#include <map>
#include <set>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

struct GrainData
{
    label grainID;
    label phaseID;
    scalar phi1;  // Euler angle (radians)
    scalar Phi;   // Euler angle (radians)  
    scalar phi2;  // Euler angle (radians)
    label cellCount;  // Number of cells in this grain
};

// Read an HDF5 dataset
template<typename T>
std::vector<T> readHDF5Dataset(hid_t file_id, const std::string& dataset_path, std::vector<hsize_t>& dims)
{
    hid_t dataset = H5Dopen(file_id, dataset_path.c_str(), H5P_DEFAULT);
    if (dataset < 0)
    {
        FatalErrorInFunction
            << "Cannot open dataset: " << dataset_path
            << exit(FatalError);
    }
    
    hid_t dataspace = H5Dget_space(dataset);
    int rank = H5Sget_simple_extent_ndims(dataspace);
    dims.resize(rank);
    H5Sget_simple_extent_dims(dataspace, dims.data(), NULL);
    
    hsize_t total_size = 1;
    for (int i = 0; i < rank; i++)
    {
        total_size *= dims[i];
    }
    
    std::vector<T> data(total_size);
    
    hid_t memtype;
    if (std::is_same<T, int>::value)
        memtype = H5T_NATIVE_INT;
    else if (std::is_same<T, float>::value)
        memtype = H5T_NATIVE_FLOAT;
    else if (std::is_same<T, double>::value)
        memtype = H5T_NATIVE_DOUBLE;
    
    H5Dread(dataset, memtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, data.data());
    
    H5Sclose(dataspace);
    H5Dclose(dataset);
    
    return data;
}

// Read the image geometry
void readGeometry(hid_t file_id, const std::string& container_name, 
                  vector& origin, vector& spacing, labelList& dimensions)
{
    std::vector<hsize_t> dims;
    
    // Read dimensions
    std::string dims_path = "/DataContainers/" + container_name + "/_SIMPL_GEOMETRY/DIMENSIONS";
    std::vector<int> dims_data = readHDF5Dataset<int>(file_id, dims_path, dims);
    dimensions.setSize(3);
    dimensions[0] = dims_data[0];
    dimensions[1] = dims_data[1];
    dimensions[2] = dims_data[2];
    
    // Read origin
    std::string origin_path = "/DataContainers/" + container_name + "/_SIMPL_GEOMETRY/ORIGIN";
    std::vector<float> origin_data = readHDF5Dataset<float>(file_id, origin_path, dims);
    origin = vector(origin_data[0], origin_data[1], origin_data[2]);
    
    // Read spacing/resolution
    std::string spacing_path = "/DataContainers/" + container_name + "/_SIMPL_GEOMETRY/SPACING";
    std::vector<float> spacing_data = readHDF5Dataset<float>(file_id, spacing_path, dims);
    spacing = vector(spacing_data[0], spacing_data[1], spacing_data[2]);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Generate Eulerian CP fields from DREAM3D microstructure"
    );
    
    argList::addOption
    (
        "file",
        "filename",
        "DREAM3D file to read (default: microstructure.dream3d)"
    );
    
    argList::addOption
    (
        "container",
        "name",
        "Data container name (default: SyntheticVolumeDataContainer)"
    );
    
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createMesh.H"
    
    // Read parameters
    
    word dream3dFile = args.optionLookupOrDefault<word>("file", "microstructure.dream3d");
    word containerName = args.optionLookupOrDefault<word>("container", "SyntheticVolumeDataContainer");
    
    fileName dream3dPath = runTime.constant()/dream3dFile;
    
    Info<< "\nReading DREAM3D file: " << dream3dPath << endl;
    Info<< "Data container: " << containerName << endl;
    
    if (!isFile(dream3dPath))
    {
        FatalErrorInFunction
            << "DREAM3D file not found: " << dream3dPath << nl
            << "Place your .dream3d file in the constant/ directory"
            << exit(FatalError);
    }
    
    // Open the HDF5 file
    
    hid_t file_id = H5Fopen(dream3dPath.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    if (file_id < 0)
    {
        FatalErrorInFunction
            << "Cannot open HDF5 file: " << dream3dPath
            << exit(FatalError);
    }
    
    Info<< "Successfully opened HDF5 file" << endl;
    
    // Read geometry
    
    Info<< "\nReading geometry..." << endl;
    
    vector origin(0, 0, 0);
    vector spacing(1, 1, 1);
    labelList dimensions(3);
    
    readGeometry(file_id, containerName.c_str(), origin, spacing, dimensions);
    
    Info<< "  Grid dimensions: " << dimensions[0] << " x " 
        << dimensions[1] << " x " << dimensions[2] << endl;
    Info<< "  Origin: " << origin << endl;
    Info<< "  Spacing: " << spacing << endl;
    
    label nx = dimensions[0];
    label ny = dimensions[1];
    label nz = dimensions[2];
    label totalCells = nx * ny * nz;
    
    // Calculate domain bounds
    vector minPoint = origin;
    vector maxPoint = origin + vector(
        nx * spacing.x(),
        ny * spacing.y(),
        nz * spacing.z()
    );
    
    Info<< "  Domain: " << minPoint << " to " << maxPoint << endl;
    
    // Read cell data
    
    Info<< "\nReading cell data..." << endl;
    
    std::vector<hsize_t> dims;
    
    // Read FeatureIds (grain IDs)
    std::string featureIds_path = "/DataContainers/" + word(containerName) + "/CellData/FeatureIds";
    std::vector<int> featureIds = readHDF5Dataset<int>(file_id, featureIds_path, dims);
    
    Info<< "  FeatureIds: " << featureIds.size() << " cells" << endl;
    
    // Read EulerAngles
    std::string eulerAngles_path = "/DataContainers/" + word(containerName) + "/CellData/EulerAngles";
    std::vector<float> eulerAngles = readHDF5Dataset<float>(file_id, eulerAngles_path, dims);
    
    Info<< "  EulerAngles: " << eulerAngles.size() / 3 << " cells x 3 components" << endl;
    
    // Read Phases
    std::string phases_path = "/DataContainers/" + word(containerName) + "/CellData/Phases";
    std::vector<int> phases = readHDF5Dataset<int>(file_id, phases_path, dims);
    
    Info<< "  Phases: " << phases.size() << " cells" << endl;
    
    H5Fclose(file_id);
    
    Info<< "Successfully read DREAM3D data" << endl;
    
    // Grain structure
    
    Info<< "\nAnalyzing grain structure..." << endl;
    
    // Build map of grains
    std::map<label, GrainData> grainMap;
    
    for (label i = 0; i < totalCells; i++)
    {
        label grainID = featureIds[i];
        
        if (grainID > 0)  // Skip grain ID 0 (typically void/background)
        {
            if (grainMap.find(grainID) == grainMap.end())
            {
                // New grain - initialize with first occurrence
                GrainData grain;
                grain.grainID = grainID;
                grain.phaseID = phases[i];
                grain.phi1 = eulerAngles[3*i + 0];
                grain.Phi = eulerAngles[3*i + 1];
                grain.phi2 = eulerAngles[3*i + 2];
                grain.cellCount = 1;
                grainMap[grainID] = grain;
            }
            else
            {
                grainMap[grainID].cellCount++;
            }
        }
    }
    
    label nGrains = grainMap.size();
    Info<< "  Found " << nGrains << " grains" << endl;
    
    // Get unique phases
    std::set<label> uniquePhases;
    for (const auto& pair : grainMap)
    {
        uniquePhases.insert(pair.second.phaseID);
    }
    
    label nPhases = uniquePhases.size();
    Info<< "  Found " << nPhases << " unique phases" << endl;
    
    // Print grain statistics
    for (const auto& pair : grainMap)
    {
        const GrainData& grain = pair.second;
        scalar percentage = 100.0 * scalar(grain.cellCount) / scalar(totalCells);
        Info<< "    Grain " << grain.grainID 
            << " (Phase " << grain.phaseID << "): "
            << grain.cellCount << " cells ("
            << percentage << " %), "
            << "Euler (deg) = (" 
            << grain.phi1 * 180.0 / constant::mathematical::pi << ", "
            << grain.Phi * 180.0 / constant::mathematical::pi << ", "
            << grain.phi2 * 180.0 / constant::mathematical::pi << ")"
            << endl;
    }
    
    // Update blockMeshDict
    
    Info<< "\nUpdating system/blockMeshDict..." << endl;
    
    fileName blockMeshPath = runTime.system()/"blockMeshDict";
    
    // Read existing blockMeshDict
    IOdictionary blockMeshDict
    (
        IOobject
        (
            "blockMeshDict",
            runTime.system(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );
    
    // Update vertices
    List<point> vertices(8);
    vertices[0] = point(minPoint.x(), minPoint.y(), minPoint.z());
    vertices[1] = point(maxPoint.x(), minPoint.y(), minPoint.z());
    vertices[2] = point(maxPoint.x(), maxPoint.y(), minPoint.z());
    vertices[3] = point(minPoint.x(), maxPoint.y(), minPoint.z());
    vertices[4] = point(minPoint.x(), minPoint.y(), maxPoint.z());
    vertices[5] = point(maxPoint.x(), minPoint.y(), maxPoint.z());
    vertices[6] = point(maxPoint.x(), maxPoint.y(), maxPoint.z());
    vertices[7] = point(minPoint.x(), maxPoint.y(), maxPoint.z());
    
    // Write updated blockMeshDict
    OFstream blockMeshFile(blockMeshPath);
    
    blockMeshFile
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
        << "    object      blockMeshDict;" << nl
        << "}" << nl
        << "// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //" << nl
        << nl
        << "// Generated from DREAM3D file: " << dream3dFile << nl
        << "// Grid: " << nx << " x " << ny << " x " << nz << nl
        << nl
        << "scale   1;" << nl
        << nl
        << "vertices" << nl
        << "(" << nl;
    
    forAll(vertices, i)
    {
        blockMeshFile << "    (" << vertices[i].x() << " " 
                      << vertices[i].y() << " " << vertices[i].z() << ")" << nl;
    }
    
    blockMeshFile
        << ");" << nl
        << nl
        << "blocks" << nl
        << "(" << nl
        << "    hex (0 1 2 3 4 5 6 7) (" << nx << " " << ny << " " << nz << ") simpleGrading (1 1 1)" << nl
        << ");" << nl
        << nl
        << "edges" << nl
        << "(" << nl
        << ");" << nl
        << nl
        << "boundary" << nl
        << "(" << nl
        << "    xMin" << nl
        << "    {" << nl
        << "        type patch;" << nl
        << "        faces" << nl
        << "        (" << nl
        << "            (0 4 7 3)" << nl
        << "        );" << nl
        << "    }" << nl
        << "    xMax" << nl
        << "    {" << nl
        << "        type patch;" << nl
        << "        faces" << nl
        << "        (" << nl
        << "            (1 2 6 5)" << nl
        << "        );" << nl
        << "    }" << nl
        << "    yMin" << nl
        << "    {" << nl
        << "        type patch;" << nl
        << "        faces" << nl
        << "        (" << nl
        << "            (0 1 5 4)" << nl
        << "        );" << nl
        << "    }" << nl
        << "    yMax" << nl
        << "    {" << nl
        << "        type patch;" << nl
        << "        faces" << nl
        << "        (" << nl
        << "            (3 7 6 2)" << nl
        << "        );" << nl
        << "    }" << nl
        << "    zMin" << nl
        << "    {" << nl
        << "        type patch;" << nl
        << "        faces" << nl
        << "        (" << nl
        << "            (0 3 2 1)" << nl
        << "        );" << nl
        << "    }" << nl
        << "    zMax" << nl
        << "    {" << nl
        << "        type patch;" << nl
        << "        faces" << nl
        << "        (" << nl
        << "            (4 5 6 7)" << nl
        << "        );" << nl
        << "    }" << nl
        << ");" << nl
        << nl
        << "mergePatchPairs" << nl
        << "(" << nl
        << ");" << nl
        << nl
        << "// ************************************************************************* //" << nl;
    
    Info<< "  Updated blockMeshDict with domain: " << minPoint << " to " << maxPoint << endl;
    Info<< "  Grid cells: " << nx << " x " << ny << " x " << nz << endl;
    
    // Write GrainProperties
    
    Info<< "\nWriting constant/GrainProperties..." << endl;
    
    OFstream grainPropsFile(runTime.constant()/"GrainProperties");
    
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
        << "    location    \"constant\";" << nl
        << "    object      GrainProperties;" << nl
        << "}" << nl
        << "// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //" << nl
        << nl
        << "// Generated from DREAM3D file: " << dream3dFile << nl
        << nl
        << "N_Phases\t" << nGrains << ";" << nl
        << nl
        << "N_Materials\t" << nPhases << ";" << nl
        << nl;
    
    // Phase2Material: DREAM3D grain IDs need not be sequential, so map the sorted IDs to indices
    std::vector<label> grainIndices;
    for (const auto& pair : grainMap)
    {
        grainIndices.push_back(pair.first);
    }
    std::sort(grainIndices.begin(), grainIndices.end());
    
    grainPropsFile << "Phase2Material (";
    for (size_t i = 0; i < grainIndices.size(); i++)
    {
        label grainID = grainIndices[i];
        label phaseID = grainMap[grainID].phaseID;
        // Assume phase IDs start from 1, map to 0-based material index
        label materialID = max(0, phaseID - 1);
        grainPropsFile << materialID;
        if (i < grainIndices.size() - 1)
            grainPropsFile << " ";
    }
    grainPropsFile << ");" << nl << nl;
    
    // Write Euler angles
    grainPropsFile << "phi1 (";
    for (size_t i = 0; i < grainIndices.size(); i++)
    {
        label grainID = grainIndices[i];
        grainPropsFile << grainMap[grainID].phi1;
        if (i < grainIndices.size() - 1)
            grainPropsFile << " ";
    }
    grainPropsFile << ");" << nl << nl;
    
    grainPropsFile << "phi2 (";
    for (size_t i = 0; i < grainIndices.size(); i++)
    {
        label grainID = grainIndices[i];
        grainPropsFile << grainMap[grainID].Phi;
        if (i < grainIndices.size() - 1)
            grainPropsFile << " ";
    }
    grainPropsFile << ");" << nl << nl;
    
    grainPropsFile << "phi3 (";
    for (size_t i = 0; i < grainIndices.size(); i++)
    {
        label grainID = grainIndices[i];
        grainPropsFile << grainMap[grainID].phi2;
        if (i < grainIndices.size() - 1)
            grainPropsFile << " ";
    }
    grainPropsFile << ");" << nl << nl;
    
    grainPropsFile << "// Add material properties below as needed" << nl;
    grainPropsFile << "// See Material.0, Material.1, etc." << nl;
    grainPropsFile << nl
        << "// ************************************************************************* //" << nl;
    
    Info<< "  Wrote " << nGrains << " grain orientations" << endl;
    
    // Initialise the n.* fields
    
    Info<< "\nInitializing PopBal fields..." << endl;
    Info<< "  NOTE: You must run 'blockMesh' first to create the mesh!" << endl;
    Info<< "  After running blockMesh, run this utility again to initialize fields." << endl;
    
    // Check if mesh exists and matches DREAM3D dimensions
    if (mesh.nCells() != totalCells)
    {
        WarningInFunction
            << "Mesh cell count (" << mesh.nCells() 
            << ") does not match DREAM3D (" << totalCells << ")" << nl
            << "Run 'blockMesh' with the updated blockMeshDict, then run this utility again"
            << endl;
    }
    else
    {
        Info<< "  Mesh matches DREAM3D dimensions, initializing fields..." << endl;
        
        // Create mapping from DREAM3D grain IDs to field indices (0 to nGrains-1)
        std::map<label, label> grainIDToIndex;
        for (size_t i = 0; i < grainIndices.size(); i++)
        {
            grainIDToIndex[grainIndices[i]] = i;
        }
        
        // Initialize PopBal fields
        for (label h = 0; h < nGrains; h++)
        {
            word fieldName = "n." + name(h);
            
            volScalarField PopBal
            (
                IOobject
                (
                    fieldName,
                    runTime.timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                mesh,
                dimensionedScalar("ni", dimensionSet(0, 0, 0, 0, 0), 0.0),
                zeroGradientFvPatchScalarField::typeName
            );
            
            // Initialize all cells to 0
            PopBal = dimensionedScalar("zero", dimensionSet(0, 0, 0, 0, 0), 0.0);
            
            // Set cells belonging to this grain to 1
            label grainID = grainIndices[h];
            label cellsSet = 0;
            
            forAll(mesh.C(), cellI)
            {
                // Calculate DREAM3D cell index from OpenFOAM cell center
                vector pos = mesh.C()[cellI];
                
                // Convert position to DREAM3D indices
                label ix = label((pos.x() - minPoint.x()) / spacing.x());
                label iy = label((pos.y() - minPoint.y()) / spacing.y());
                label iz = label((pos.z() - minPoint.z()) / spacing.z());
                
                // Clamp to valid range
                ix = max(0, min(ix, nx - 1));
                iy = max(0, min(iy, ny - 1));
                iz = max(0, min(iz, nz - 1));
                
                // DREAM3D stores in Z-Y-X order (Fortran-style)
                label dream3dIdx = iz * (ny * nx) + iy * nx + ix;
                
                if (dream3dIdx < totalCells && featureIds[dream3dIdx] == grainID)
                {
                    PopBal[cellI] = 1.0;
                    cellsSet++;
                }
            }
            
            PopBal.correctBoundaryConditions();
            PopBal.write();
            
            Info<< "    Field " << fieldName << ": " << cellsSet 
                << " cells set to 1.0 (Grain " << grainID << ")" << endl;
        }
        
        Info<< "  Successfully initialized " << nGrains << " PopBal fields" << endl;
    }
    
    // Summary
    
    Info<< "\n========================================" << endl;
    Info<< "DREAM3D import complete!" << endl;
    Info<< "========================================" << endl;
    Info<< "\nNext steps:" << endl;
    Info<< "  1. Run 'blockMesh' to create the mesh" << endl;
    Info<< "  2. Run this utility again to initialize PopBal fields" << endl;
    Info<< "  3. Add material properties to constant/GrainProperties" << endl;
    Info<< "  4. Initialize other required fields (U, T, etc.)" << endl;
    Info<< "  5. Run your Eulerian CP solver" << endl;
    Info<< "\nGenerated files:" << endl;
    Info<< "  - system/blockMeshDict (domain: " << nx << "x" << ny << "x" << nz << ")" << endl;
    Info<< "  - constant/GrainProperties (" << nGrains << " grains, " << nPhases << " phases)" << endl;
    if (mesh.nCells() == totalCells)
    {
        Info<< "  - 0/n.* (PopBal fields initialized)" << endl;
    }
    Info<< endl;
    
    return 0;
}

// ************************************************************************* //
