/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           |
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    ebsdToFoam

Description
    Maps an EBSD export (X Y Euler1 Euler2 Euler3, Bunge) onto the mesh and writes n.i,
    ROT.k.i and constant/<grainPropertiesInclude>; controlled by system/ebsdToFoamDict.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "Tuple2.H"
#include "OFstream.H"
#include "IStringStream.H"

#include <iomanip>
#include <sstream>

#include "orientation.H"
#include "ipfColour.H"
#include "ebsdMap.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace
{

ebsd::Mat3 toMat3(const tensor& t)
{
    ebsd::Mat3 m;
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
            m(i, j) = t(i, j);
    return m;
}

ebsd::Vec3 toVec3(const vector& v)
{
    return ebsd::Vec3(v.x(), v.y(), v.z());
}

//- Rotation taking Euler-frame components to map-frame components (MTEX "EulerCorrection")
ebsd::Mat3 vendorCorrection
(
    const word& name,
    const vector& customAxis,
    const scalar customAngleDeg
)
{
    if (name == "none" || name == "identity")
    {
        return ebsd::Mat3::identity();
    }
    else if (name == "oxford" || name == "ctf" || name == "channel5" || name == "aztec")
    {
        return ebsd::Mat3::axisAngle(ebsd::Vec3(1, 0, 0), ebsd::PI);
    }
    else if (name == "edax1" || name == "tsl1")
    {
        return ebsd::Mat3::axisAngle(ebsd::Vec3(1, 1, 0), ebsd::PI);
    }
    else if (name == "edax2" || name == "tsl2" || name == "edax" || name == "tsl")
    {
        return ebsd::Mat3::axisAngle(ebsd::Vec3(1, -1, 0), ebsd::PI);
    }
    else if (name == "edax3" || name == "tsl3")
    {
        return ebsd::Mat3::axisAngle(ebsd::Vec3(1, 0, 0), ebsd::PI);
    }
    else if (name == "edax4" || name == "tsl4")
    {
        return ebsd::Mat3::axisAngle(ebsd::Vec3(0, 1, 0), ebsd::PI);
    }
    else if (name == "custom")
    {
        return ebsd::Mat3::axisAngle(toVec3(customAxis), ebsd::deg2rad(customAngleDeg));
    }

    FatalErrorInFunction
        << "Unknown vendorCorrection '" << name << "'. Valid: none, oxford, "
        << "edax1, edax2, edax3, edax4, custom" << exit(FatalError);
    return ebsd::Mat3::identity();
}

//- Steinbach profile 0.5(1 + sin(pi d/eta)), d the signed distance into the grain; 1 beyond eta/2
scalar interfaceProfile(const scalar d, const scalar eta)
{
    if (eta <= 0) return 1.0;
    if (d >= 0.5*eta) return 1.0;
    if (d <= -0.5*eta) return 0.0;
    return 0.5*(1.0 + Foam::sin(constant::mathematical::pi*d/eta));
}

//- Snap round-off (|x| < 1e-12) to zero for printing
double clean(double x) { return (std::fabs(x) < 1e-12) ? 0.0 : x; }

std::string fmtAngles(double a, double b, double c)
{
    std::ostringstream os;
    os << std::fixed << std::setprecision(2)
       << "(" << std::setw(7) << ebsd::rad2deg(a) << " " << std::setw(6) << ebsd::rad2deg(b)
       << " " << std::setw(7) << ebsd::rad2deg(c) << ")";
    return os.str();
}

} // anonymous namespace


int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Map an EBSD orientation map (X Y Euler1 Euler2 Euler3 text export) onto the mesh\n"
        "and write the crystaldynamicsFoam phase-field order parameters n.i and\n"
        "lattice orientation fields ROT.1.i, ROT.2.i, ROT.3.i."
    );
    argList::addOption("dict", "file", "alternative dictionary (default system/ebsdToFoamDict)");
    argList::addBoolOption("dryRun", "analyse the EBSD map and report; write nothing");

    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    const bool dryRun = args.optionFound("dryRun");

    // --- dictionary
    autoPtr<IOdictionary> dictPtr;
    if (args.optionFound("dict"))
    {
        fileName dictPath(args["dict"]);
        dictPath.expand();
        if (!dictPath.isAbsolute()) dictPath = runTime.path()/dictPath;
        dictPtr.reset
        (
            new IOdictionary
            (
                IOobject(dictPath, mesh, IOobject::MUST_READ, IOobject::NO_WRITE)
            )
        );
    }
    else
    {
        dictPtr.reset
        (
            new IOdictionary
            (
                IOobject
                (
                    "ebsdToFoamDict",
                    runTime.system(),
                    mesh,
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE
                )
            )
        );
    }
    const IOdictionary& dict = dictPtr();

    // --- input
    fileName ebsdFile(dict.lookup("file"));
    ebsdFile.expand();
    if (!ebsdFile.isAbsolute()) ebsdFile = runTime.path()/ebsdFile;

    const scalar lengthUnit = dict.lookupOrDefault<scalar>("lengthUnit", 1e-6);
    const word angleUnit = dict.lookupOrDefault<word>("angleUnit", "degrees");
    const bool zeroEulerIsUnindexed = dict.lookupOrDefault<Switch>("zeroEulerIsUnindexed", true);
    const word symmetryName = dict.lookupOrDefault<word>("symmetry", "cubic");

    // --- reference frames
    const word correctionName = dict.lookupOrDefault<word>("vendorCorrection", "oxford");
    const vector correctionAxis = dict.lookupOrDefault<vector>("correctionAxis", vector(1, 0, 0));
    const scalar correctionAngle = dict.lookupOrDefault<scalar>("correctionAngle", 180.0);
    const word meshFrame = dict.lookupOrDefault<word>("meshFrame", "euler");
    const vector frameRotationAxis = dict.lookupOrDefault<vector>("frameRotationAxis", vector(0, 0, 1));
    const scalar frameRotationAngle = dict.lookupOrDefault<scalar>("frameRotationAngle", 0.0);
    const vector origin = dict.lookupOrDefault<vector>("origin", vector::zero);

    // --- grain reconstruction
    const scalar tolDeg = dict.lookupOrDefault<scalar>("misorientationTolerance", 5.0);
    const label minGrainPixels = dict.lookupOrDefault<label>("minGrainPixels", 10);
    const word orientationMode = dict.lookupOrDefault<word>("orientationMode", "grainMean");
    const scalar outlierAngle = dict.lookupOrDefault<scalar>("outlierAngle", 3.0*tolDeg);

    // --- phase-field order parameters
    const word phaseAssignment = dict.lookupOrDefault<word>("phaseAssignment", "onePerGrain");
    const scalar minPhaseSeparation = dict.lookupOrDefault<scalar>("minPhaseSeparation", 4.0);
    const scalar interfaceWidth = dict.lookupOrDefault<scalar>("interfaceWidth", 0.0);
    const label materialIndex = dict.lookupOrDefault<label>("material", 0);

    // --- buffer: the mesh outside the map becomes one extra order parameter of another material
    const dictionary bufferDict(dict.subOrEmptyDict("buffer"));
    const bool useBuffer = bufferDict.lookupOrDefault<Switch>("enabled", false);
    const label bufferMaterial = bufferDict.lookupOrDefault<label>("material", 1);
    const scalar bufferMargin = bufferDict.lookupOrDefault<scalar>("margin", 0.0);
    const vector bufferColour = bufferDict.lookupOrDefault<vector>("ipfColour", vector(0.6, 0.6, 0.6));
    if (useBuffer && bufferMaterial == materialIndex)
    {
        FatalErrorInFunction
            << "buffer/material must differ from the material of the grains (" << materialIndex << ")"
            << exit(FatalError);
    }

    // --- output
    const bool writeIPF = dict.lookupOrDefault<Switch>("writeIPF", true);
    const bool writeGrainID = dict.lookupOrDefault<Switch>("writeGrainID", true);
    const word ipfKeyName = dict.lookupOrDefault<word>("ipfKey", "hsv");
    const scalar ipfSatPower = dict.lookupOrDefault<scalar>("ipfSaturationPower", 1.0);
    const word includeName = dict.lookupOrDefault<word>("grainPropertiesInclude", "ebsdGrains");

    List<Tuple2<word, vector>> ipfDirections;
    if (dict.found("ipfDirections"))
    {
        dict.lookup("ipfDirections") >> ipfDirections;
    }
    else
    {
        ipfDirections.setSize(3);
        ipfDirections[0] = Tuple2<word, vector>("X", vector(1, 0, 0));
        ipfDirections[1] = Tuple2<word, vector>("Y", vector(0, 1, 0));
        ipfDirections[2] = Tuple2<word, vector>("Z", vector(0, 0, 1));
    }

    if (orientationMode != "grainMean" && orientationMode != "perPixel")
    {
        FatalErrorInFunction << "orientationMode must be grainMean or perPixel" << exit(FatalError);
    }
    if (phaseAssignment != "onePerGrain" && phaseAssignment != "graphColouring")
    {
        FatalErrorInFunction << "phaseAssignment must be onePerGrain or graphColouring" << exit(FatalError);
    }
    if (meshFrame != "euler" && meshFrame != "map")
    {
        FatalErrorInFunction << "meshFrame must be euler or map" << exit(FatalError);
    }
    if (mesh.nSolutionD() != 2)
    {
        WarningInFunction
            << "The mesh has " << mesh.nSolutionD()
            << " solution directions; ebsdToFoam expects a one-cell-thick mesh in the plane of the map "
            << "(the x-y plane unless frameRotation moves it) and ignores the thickness coordinate."
            << nl << endl;
    }

    const ebsd::Symmetry sym = ebsd::symmetryFromName(symmetryName);
    const ebsd::IpfKey ipfKey = ebsd::ipfKeyFromName(ipfKeyName);

    // --- read EBSD map
    ebsd::ReadOptions ropt;
    ropt.degrees = (angleUnit == "degrees" || angleUnit == "deg");
    ropt.zeroEulerIsUnindexed = zeroEulerIsUnindexed;
    if (dict.found("columns"))
    {
        const dictionary& cd = dict.subDict("columns");
        ropt.xName  = word(cd.lookupOrDefault<word>("x", "X"));
        ropt.yName  = word(cd.lookupOrDefault<word>("y", "Y"));
        ropt.e1Name = word(cd.lookupOrDefault<word>("euler1", "Euler1"));
        ropt.e2Name = word(cd.lookupOrDefault<word>("euler2", "Euler2"));
        ropt.e3Name = word(cd.lookupOrDefault<word>("euler3", "Euler3"));
        ropt.xCol  = cd.lookupOrDefault<label>("xColumn", -1);
        ropt.yCol  = cd.lookupOrDefault<label>("yColumn", -1);
        ropt.e1Col = cd.lookupOrDefault<label>("euler1Column", -1);
        ropt.e2Col = cd.lookupOrDefault<label>("euler2Column", -1);
        ropt.e3Col = cd.lookupOrDefault<label>("euler3Column", -1);
    }

    Info<< "Reading EBSD map " << ebsdFile << nl << endl;
    ebsd::EbsdMap map;
    try
    {
        map = ebsd::readEbsdText(ebsdFile, ropt, std::cout);
        ebsd::fillUnindexed(map, std::cout);
    }
    catch (const std::exception& e)
    {
        FatalErrorInFunction << e.what() << exit(FatalError);
    }

    // --- grains
    Info<< nl << "Reconstructing grains (tolerance " << tolDeg << " deg, symmetry "
        << ebsd::symmetryName(sym).c_str() << ")" << endl;
    ebsd::segmentGrains(map, ebsd::deg2rad(tolDeg), sym);
    Info<< "  " << map.nGrains << " grains before merging" << endl;
    ebsd::mergeSmallGrains(map, minGrainPixels, std::cout);
    Info<< "  " << map.nGrains << " grains" << endl;

    std::vector<ebsd::Quat> grainMean = ebsd::grainMeanOrientations(map, sym);
    std::vector<int> grainPixels = ebsd::grainSizes(map);

    if (orientationMode == "perPixel")
    {
        const std::size_t nOut = ebsd::cleanOutliers(map, grainMean, ebsd::deg2rad(outlierAngle), sym);
        Info<< "  perPixel orientations: " << label(nOut) << " outlier pixels (> " << outlierAngle
            << " deg from the grain mean) replaced by the grain mean" << endl;
        grainMean = ebsd::grainMeanOrientations(map, sym);
    }

    // intragranular spread (diagnostic)
    std::vector<double> maxSpread(map.nGrains, 0.0);
    for (std::size_t p = 0; p < map.size(); p++)
    {
        const double a = ebsd::misorientationAngle(map.q[p], grainMean[map.grain[p]], sym);
        if (a > maxSpread[map.grain[p]]) maxSpread[map.grain[p]] = a;
    }

    // --- phases
    std::vector<int> phaseOfGrain(map.nGrains, 0);
    label nPhases = 0;
    if (phaseAssignment == "onePerGrain")
    {
        for (int g = 0; g < map.nGrains; g++) phaseOfGrain[g] = g;
        nPhases = map.nGrains;
    }
    else
    {
        int nColours = 0;
        phaseOfGrain = ebsd::graphColouring(map, minPhaseSeparation, nColours);
        nPhases = nColours;
        Info<< "  graph colouring with a minimum separation of " << minPhaseSeparation
            << " pixels: " << nPhases << " order parameters for " << map.nGrains << " grains" << endl;
    }
    // representative grain per phase (largest)
    std::vector<int> representative(nPhases, -1);
    for (int g = 0; g < map.nGrains; g++)
    {
        const int p = phaseOfGrain[g];
        if (representative[p] < 0 || grainPixels[g] > grainPixels[representative[p]]) representative[p] = g;
    }

    // --- frames
    // Q: Euler-frame components -> map-frame components
    const ebsd::Mat3 Q = vendorCorrection(correctionName, correctionAxis, correctionAngle);
    // P: optional rigid rotation of the whole configuration (positions and orientations together)
    const ebsd::Mat3 P =
        (mag(frameRotationAngle) > 0)
      ? ebsd::Mat3::axisAngle(toVec3(frameRotationAxis), ebsd::deg2rad(frameRotationAngle))
      : ebsd::Mat3::identity();
    // Q_FM: map -> mesh,  Q_FE: Euler -> mesh
    const ebsd::Mat3 Q_FM = P*((meshFrame == "euler") ? Q.T() : ebsd::Mat3::identity());
    const ebsd::Mat3 Q_FE = P*((meshFrame == "euler") ? ebsd::Mat3::identity() : Q);

    // physical extent of the map in the mesh frame (pixel rectangle expanded by half a pixel)
    const double L = lengthUnit;
    const double xLo = (map.x0 - 0.5*map.dx)*L, xHi = (map.x0 + (map.nx - 0.5)*map.dx)*L;
    const double yLo = (map.y0 - 0.5*map.dy)*L, yHi = (map.y0 + (map.ny - 0.5)*map.dy)*L;
    ebsd::Vec3 lo(1e300, 1e300, 1e300), hi(-1e300, -1e300, -1e300);
    {
        const double xs[2] = {xLo, xHi}, ys[2] = {yLo, yHi};
        for (int a = 0; a < 2; a++)
            for (int b = 0; b < 2; b++)
            {
                const ebsd::Vec3 r = Q_FM*ebsd::Vec3(xs[a], ys[b], 0.0);
                for (int k = 0; k < 3; k++)
                {
                    lo[k] = std::min(lo[k], r[k]);
                    hi[k] = std::max(hi[k], r[k]);
                }
            }
    }
    for (int k = 0; k < 3; k++) { lo[k] = clean(lo[k]); hi[k] = clean(hi[k]); }
    // translation placing the lower corner of the map at 'origin'
    const ebsd::Vec3 t = toVec3(origin) - lo;
    const ebsd::Mat3 Q_MF = Q_FM.T();

    // mesh-frame directions of the map axes and of the map normal (the thickness direction)
    ebsd::Vec3 mapXdir = Q_FM*ebsd::Vec3(1, 0, 0);
    ebsd::Vec3 mapYdir = Q_FM*ebsd::Vec3(0, 1, 0);
    ebsd::Vec3 mapNdir = Q_FM*ebsd::Vec3(0, 0, 1);
    for (int k = 0; k < 3; k++)
    {
        mapXdir[k] = clean(mapXdir[k]); mapYdir[k] = clean(mapYdir[k]); mapNdir[k] = clean(mapNdir[k]);
    }

    Info<< nl << "Reference frames" << nl
        << "  vendorCorrection " << correctionName << " (Euler frame -> map frame)" << nl
        << "  meshFrame        " << meshFrame << nl
        << "  frameRotation    " << frameRotationAngle << " deg about " << frameRotationAxis << nl
        << "  map x axis (image right) along mesh direction (" << mapXdir.x << " " << mapXdir.y << " " << mapXdir.z << ")" << nl
        << "  map y axis (image down)  along mesh direction (" << mapYdir.x << " " << mapYdir.y << " " << mapYdir.z << ")" << nl
        << "  map z axis (into image)  along mesh direction (" << mapNdir.x << " " << mapNdir.y << " " << mapNdir.z << ")" << nl
        << "  map extent in the mesh frame: x = [" << lo.x + t.x << ", " << hi.x + t.x
        << "]  y = [" << lo.y + t.y << ", " << hi.y + t.y
        << "]  z = [" << lo.z + t.z << ", " << hi.z + t.z << "]  (m)" << nl
        << "  pixel size (" << map.dx*L << ", " << map.dy*L << ") m" << nl << endl;
    if (useBuffer)
    {
        Info<< "Buffer: cells outside the map extent";
        if (bufferMargin != 0) Info<< " inset by " << bufferMargin << " m";
        Info<< " form order parameter " << nPhases << " of Material." << bufferMaterial
            << " (" << (interfaceWidth > 0 ? "diffuse" : "sharp") << " interface to the grains)"
            << nl << endl;
    }

    // mesh-frame orientation of a pixel/grain quaternion
    auto meshR = [&](const ebsd::Quat& q) -> ebsd::Mat3 { return Q_FE*q.toMatrix(); };

    // --- report
    {
        const double totalPixels = static_cast<double>(map.size());
        Info<< "Grains (orientations: Bunge = as exported, in the Euler frame; code = ROT.1 ROT.2 ROT.3 in the mesh frame)" << nl
            << "  grain  phase   pixels   area%   spread   Bunge(phi1 Phi phi2) deg      code(p1 p2 p3) deg"
            << "       crystal dir || mesh x        crystal dir || mesh y        crystal dir || mesh z" << endl;
        for (int g = 0; g < map.nGrains; g++)
        {
            double b1, b2, b3, c1, c2, c3;
            ebsd::matrixToBunge(grainMean[g].toMatrix().T(), b1, b2, b3);
            const ebsd::Mat3 R = meshR(grainMean[g]);
            ebsd::matrixToCodeAngles(R, c1, c2, c3);
            std::ostringstream line;
            line << std::fixed << std::setprecision(2);
            line << "  " << std::setw(5) << g << "  " << std::setw(5) << phaseOfGrain[g]
                 << "  " << std::setw(7) << grainPixels[g]
                 << "  " << std::setw(6) << 100.0*grainPixels[g]/totalPixels
                 << "  " << std::setw(6) << ebsd::rad2deg(maxSpread[g])
                 << "   " << fmtAngles(b1, b2, b3)
                 << "   " << fmtAngles(c1, c2, c3);
            const ebsd::Vec3 ax[3] = {ebsd::Vec3(1,0,0), ebsd::Vec3(0,1,0), ebsd::Vec3(0,0,1)};
            for (int k = 0; k < 3; k++)
            {
                const ebsd::Vec3 d = ebsd::toFundamentalSector(R.T()*ax[k], sym);
                line << "   [" << std::setw(5) << d.x << " " << std::setw(5) << d.y << " " << std::setw(5) << d.z << "]";
            }
            Info<< line.str().c_str() << endl;
        }
        Info<< endl;
    }

    if (dryRun)
    {
        Info<< "dryRun: nothing written" << nl << "End" << endl;
        return 0;
    }

    // --- cells -> pixels
    const label nCells = mesh.nCells();
    std::vector<std::size_t> cellPixel(nCells);
    std::vector<double> cellU(nCells), cellV(nCells);
    // crystalline fraction per cell: 1 without a buffer, else the Steinbach profile of the
    // signed distance to the (inset) map rectangle
    std::vector<double> cellF(nCells, 1.0);
    label nOutside = 0, nBuffer = 0;
    const double cxMap = 0.5*(xLo + xHi), cyMap = 0.5*(yLo + yHi);
    const double hxMap = 0.5*(xHi - xLo), hyMap = 0.5*(yHi - yLo);
    forAll(mesh.C(), celli)
    {
        const vector& C = mesh.C()[celli];
        const ebsd::Vec3 rM = Q_MF*(toVec3(C) - t);
        const double u = (rM.x/L - map.x0)/map.dx;
        const double v = (rM.y/L - map.y0)/map.dy;
        int i = static_cast<int>(std::floor(u + 0.5));
        int j = static_cast<int>(std::floor(v + 0.5));
        if (i < 0 || i >= map.nx || j < 0 || j >= map.ny) nOutside++;
        i = std::max(0, std::min(map.nx - 1, i));
        j = std::max(0, std::min(map.ny - 1, j));
        cellPixel[celli] = map.index(i, j);
        cellU[celli] = u; cellV[celli] = v;

        if (useBuffer)
        {
            // signed distance to the map rectangle (positive outside), then to the inset one
            const double qx = std::fabs(rM.x - cxMap) - hxMap;
            const double qy = std::fabs(rM.y - cyMap) - hyMap;
            const double dOut =
                std::sqrt(sqr(std::max(qx, 0.0)) + sqr(std::max(qy, 0.0)))
              + std::min(std::max(qx, qy), 0.0);
            const double dIn = -dOut - bufferMargin;   // positive inside the crystalline region
            cellF[celli] = (interfaceWidth > 0) ? interfaceProfile(dIn, interfaceWidth)
                                                : (dIn >= 0.0 ? 1.0 : 0.0);
            if (cellF[celli] < 0.5) nBuffer++;
        }
    }
    if (useBuffer)
    {
        Info<< "  " << nBuffer << " of " << nCells << " cells are buffer (crystalline fraction < 0.5), "
            << nOutside << " cell centres lie outside the map" << nl << endl;
        if (nBuffer == 0)
        {
            WarningInFunction
                << "buffer enabled but no cell lies outside the (inset) map: "
                << "enlarge the mesh or set origin/margin in " << dict.name() << nl << endl;
        }
    }
    else if (nOutside)
    {
        WarningInFunction
            << nOutside << " of " << nCells << " cell centres lie outside the EBSD map; "
            << "they take the values of the nearest pixel." << nl << endl;
    }

    // own grain, own orientation, interface profile
    std::vector<int> cellGrain(nCells), cellOther(nCells, -1);
    std::vector<double> cellNOwn(nCells, 1.0);
    ebsd::BoundaryDistance bd;
    if (interfaceWidth > 0) bd.build(map);

    forAll(mesh.C(), celli)
    {
        const std::size_t p = cellPixel[celli];
        const int g = map.grain[p];
        cellGrain[celli] = g;
        if (interfaceWidth > 0)
        {
            int other = -1;
            const double d = bd.query(map, cellU[celli], cellV[celli], g, other)*L;
            cellNOwn[celli] = interfaceProfile(d, interfaceWidth);
            cellOther[celli] = (other == g) ? -1 : other;
        }
    }

    const bool perPixel = (orientationMode == "perPixel");
    auto pixelOrientation = [&](std::size_t p) -> const ebsd::Quat&
    {
        return perPixel ? map.q[p] : grainMean[map.grain[p]];
    };

    // --- write fields
    const word patchType = zeroGradientFvPatchScalarField::typeName;

    Info<< "Writing fields for " << nPhases << " grain order parameters"
        << (useBuffer ? " and one buffer order parameter" : "") << " to " << runTime.timeName() << endl;

    // remove n.i / ROT.k.i files of an earlier run with more phases, so the time directory
    // holds exactly the N_Phases fields declared in constant/<grainPropertiesInclude>
    {
        const label nPhasesTotalWrite = nPhases + (useBuffer ? 1 : 0);
        label nStale = 0;
        const fileNameList files(readDir(runTime.timePath(), fileName::FILE));
        forAll(files, fi)
        {
            const std::string s = files[fi];
            label idx = -1;
            if (s.compare(0, 2, "n.") == 0 && s.size() > 2 && s.find_first_not_of("0123456789", 2) == std::string::npos)
            {
                idx = atoi(s.c_str() + 2);
            }
            else if (s.compare(0, 4, "ROT.") == 0 && s.size() > 6 && s[5] == '.'
                  && s.find_first_not_of("0123456789", 6) == std::string::npos)
            {
                idx = atoi(s.c_str() + 6);
            }
            if (idx >= nPhasesTotalWrite)
            {
                rm(runTime.timePath()/files[fi]);
                nStale++;
            }
        }
        if (nStale)
        {
            Info<< "  removed " << nStale << " stale n.i / ROT.k.i files of order parameters >= "
                << nPhasesTotalWrite << endl;
        }
    }

    std::vector<ebsd::Vec3> phaseRepAngles(nPhases);

    for (label ph = 0; ph < nPhases; ph++)
    {
        // nearest pixel belonging to this phase (orientation extension into the far field)
        std::vector<char> feature(map.size(), 0);
        for (std::size_t p = 0; p < map.size(); p++) feature[p] = (phaseOfGrain[map.grain[p]] == ph);
        std::vector<int> nearest; std::vector<double> dist;
        ebsd::edt2D(map.nx, map.ny, map.dx, map.dy, feature, nearest, dist);

        volScalarField n
        (
            IOobject(word("n." + name(ph)), runTime.timeName(), mesh, IOobject::NO_READ, IOobject::NO_WRITE),
            mesh, dimensionedScalar("n", dimless, 0.0), patchType
        );
        volScalarField rot1
        (
            IOobject(word("ROT.1." + name(ph)), runTime.timeName(), mesh, IOobject::NO_READ, IOobject::NO_WRITE),
            mesh, dimensionedScalar("ROT.1", dimless, 0.0), patchType
        );
        volScalarField rot2
        (
            IOobject(word("ROT.2." + name(ph)), runTime.timeName(), mesh, IOobject::NO_READ, IOobject::NO_WRITE),
            mesh, dimensionedScalar("ROT.2", dimless, 0.0), patchType
        );
        volScalarField rot3
        (
            IOobject(word("ROT.3." + name(ph)), runTime.timeName(), mesh, IOobject::NO_READ, IOobject::NO_WRITE),
            mesh, dimensionedScalar("ROT.3", dimless, 0.0), patchType
        );

        forAll(mesh.C(), celli)
        {
            const std::size_t p = cellPixel[celli];
            const int g = cellGrain[celli];

            // orientation of this phase at this cell
            std::size_t src = p;
            if (phaseOfGrain[g] != ph)
            {
                src = (nearest[p] >= 0) ? static_cast<std::size_t>(nearest[p]) : p;
            }
            double p1, p2, p3;
            ebsd::matrixToCodeAngles(meshR(pixelOrientation(src)), p1, p2, p3);
            rot1[celli] = p1; rot2[celli] = p2; rot3[celli] = p3;

            // order parameter (scaled by the crystalline fraction when there is a buffer)
            scalar val = 0.0;
            if (phaseOfGrain[g] == ph) val += cellNOwn[celli];
            const int other = cellOther[celli];
            if (other >= 0 && phaseOfGrain[other] == ph) val += 1.0 - cellNOwn[celli];
            n[celli] = val*cellF[celli];
        }

        n.correctBoundaryConditions();
        rot1.correctBoundaryConditions();
        rot2.correctBoundaryConditions();
        rot3.correctBoundaryConditions();
        n.write(); rot1.write(); rot2.write(); rot3.write();

        double p1, p2, p3;
        ebsd::matrixToCodeAngles(meshR(grainMean[representative[ph]]), p1, p2, p3);
        phaseRepAngles[ph] = ebsd::Vec3(p1, p2, p3);
    }

    const label nPhasesTotal = nPhases + (useBuffer ? 1 : 0);

    if (useBuffer)
    {
        // buffer order parameter (index nPhases) with the identity orientation
        volScalarField n
        (
            IOobject(word("n." + name(nPhases)), runTime.timeName(), mesh, IOobject::NO_READ, IOobject::NO_WRITE),
            mesh, dimensionedScalar("n", dimless, 0.0), patchType
        );
        forAll(mesh.C(), celli) n[celli] = 1.0 - cellF[celli];
        n.correctBoundaryConditions();
        n.write();
        for (int k = 1; k <= 3; k++)
        {
            volScalarField rot
            (
                IOobject(word("ROT." + name(k) + "." + name(nPhases)), runTime.timeName(), mesh, IOobject::NO_READ, IOobject::NO_WRITE),
                mesh, dimensionedScalar("ROT", dimless, 0.0), patchType
            );
            rot.write();
        }

        volScalarField crystalFraction
        (
            IOobject("crystalFraction", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::NO_WRITE),
            mesh, dimensionedScalar("crystalFraction", dimless, 1.0), patchType
        );
        forAll(mesh.C(), celli) crystalFraction[celli] = cellF[celli];
        crystalFraction.correctBoundaryConditions();
        crystalFraction.write();
        Info<< "  wrote n." << nPhases << ", ROT.k." << nPhases << " (buffer, Material." << bufferMaterial
            << ") and crystalFraction" << endl;
    }

    if (writeGrainID)
    {
        volScalarField grainID
        (
            IOobject("grainID", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::NO_WRITE),
            mesh, dimensionedScalar("grainID", dimless, 0.0), patchType
        );
        volScalarField phaseID
        (
            IOobject("phaseID", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::NO_WRITE),
            mesh, dimensionedScalar("phaseID", dimless, 0.0), patchType
        );
        forAll(mesh.C(), celli)
        {
            const bool inBuffer = useBuffer && cellF[celli] < 0.5;
            grainID[celli] = inBuffer ? -1 : cellGrain[celli];
            phaseID[celli] = inBuffer ? nPhases : phaseOfGrain[cellGrain[celli]];
        }
        grainID.correctBoundaryConditions();
        phaseID.correctBoundaryConditions();
        grainID.write();
        phaseID.write();
    }

    if (writeIPF)
    {
        forAll(ipfDirections, di)
        {
            const word fieldName("IPF_" + ipfDirections[di].first());
            const ebsd::Vec3 dS = toVec3(ipfDirections[di].second()).normalised();
            volVectorField ipf
            (
                IOobject(fieldName, runTime.timeName(), mesh, IOobject::NO_READ, IOobject::NO_WRITE),
                mesh, dimensionedVector("IPF", dimless, vector::zero), patchType
            );
            forAll(mesh.C(), celli)
            {
                if (useBuffer && cellF[celli] < 0.5)
                {
                    ipf[celli] = bufferColour;
                    continue;
                }
                const ebsd::RGB c = ebsd::ipfColour
                (
                    meshR(pixelOrientation(cellPixel[celli])), dS, sym, ipfKey, ipfSatPower
                );
                ipf[celli] = vector(c.r, c.g, c.b);
            }
            ipf.correctBoundaryConditions();
            ipf.write();
            Info<< "  wrote " << fieldName << " (sample direction " << ipfDirections[di].second()
                << ", key " << ipfKeyName << ")" << endl;
        }
    }

    // --- GrainProperties fragment
    {
        const fileName incPath = runTime.path()/runTime.constant()/includeName;
        OFstream os(incPath);
        os.precision(12);
        os  << "// Generated by ebsdToFoam from " << ebsdFile << nl
            << "// #include this file from constant/GrainProperties and remove any existing" << nl
            << "// N_Phases, Phase2Material, phi1, phi2, phi3 entries there." << nl
            << "//" << nl
            << "// Orientation convention: phi1/phi2/phi3 are the crystaldynamicsFoam angles" << nl
            << "// (ROT.1, ROT.2, ROT.3; radians; R = Rz(phi3) Ry(phi2) Rx(phi1), v_sample = R v_crystal)" << nl
            << "// of the largest grain of each order parameter, expressed in the mesh frame." << nl
            << "// The per-cell fields " << runTime.timeName() << "/ROT.k.i written alongside take precedence." << nl
            << "//" << nl
            << "// vendorCorrection " << correctionName << ", meshFrame " << meshFrame
            << ", frameRotation " << frameRotationAngle << " deg about " << frameRotationAxis
            << ", symmetry " << ebsd::symmetryName(sym).c_str() << ", orientationMode " << orientationMode << nl
            << "// phaseAssignment " << phaseAssignment << ", interfaceWidth " << interfaceWidth << nl;
        if (useBuffer)
        {
            os  << "// buffer: order parameter " << nPhases << " is the mesh outside the map (Material."
                << bufferMaterial << ", identity orientation)" << nl;
        }
        os  << "//" << nl
            << "// phase  grain(largest)  pixels   Bunge(phi1 Phi phi2) deg   code(ROT.1 ROT.2 ROT.3) deg" << nl;
        for (label ph = 0; ph < nPhases; ph++)
        {
            const int g = representative[ph];
            double b1, b2, b3;
            ebsd::matrixToBunge(grainMean[g].toMatrix().T(), b1, b2, b3);
            os  << "//  " << ph << "  " << g << "  " << grainPixels[g] << "   "
                << fmtAngles(b1, b2, b3).c_str() << "   "
                << fmtAngles(phaseRepAngles[ph].x, phaseRepAngles[ph].y, phaseRepAngles[ph].z).c_str() << nl;
        }
        if (useBuffer)
        {
            os  << "//  " << nPhases << "  buffer   " << nBuffer << " cells   (0 0 0)   (0 0 0)" << nl;
        }
        os  << nl << "N_Phases " << nPhasesTotal << ";" << nl << nl;
        os  << "Phase2Material" << nl << nPhasesTotal << nl << "(" << nl;
        for (label ph = 0; ph < nPhases; ph++) os << "    " << materialIndex << nl;
        if (useBuffer) os << "    " << bufferMaterial << nl;
        os  << ");" << nl << nl;
        const char* names[3] = {"phi1", "phi2", "phi3"};
        for (int k = 0; k < 3; k++)
        {
            os  << names[k] << nl << nPhasesTotal << nl << "(" << nl;
            for (label ph = 0; ph < nPhases; ph++) os << "    " << phaseRepAngles[ph][k] << nl;
            if (useBuffer) os << "    0" << nl;
            os  << ");" << nl << nl;
        }
        Info<< nl << "Wrote " << incPath << nl
            << "  -> add '#include \"" << includeName << "\"' to constant/GrainProperties "
            << "(and remove its own N_Phases, Phase2Material, phi1, phi2, phi3)." << endl;
    }

    Info<< nl << "End" << nl << endl;
    return 0;
}


// ************************************************************************* //