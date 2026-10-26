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
    foamToIPF

Description
    Inverse-pole-figure colouring of the lattice orientation from n.i and ROT.k.i: writes
    IPF_<name> (r g b in [0,1]) per sample direction, dominantPhase and crystalFraction.

    Usage
        foamToIPF [-latestTime | -time <ranges>] [-directions "((RD (1 0 0)) (ND (0 1 0)))"]
                  [-key hsv|tsl] [-symmetry cubic|hexagonal] [-nPhases N] [-blend]
                  [-bufferMaterial M] [-bufferColour "(0.6 0.6 0.6)"]

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "timeSelector.H"
#include "Tuple2.H"
#include "IStringStream.H"

#include "orientation.H"
#include "ipfColour.H"

using namespace Foam;

int main(int argc, char *argv[])
{
    timeSelector::addOptions();
    argList::addNote
    (
        "Write IPF colour fields IPF_<name> (r g b) from the crystaldynamicsFoam\n"
        "orientation fields ROT.k.i and order parameters n.i."
    );
    argList::addOption
    (
        "directions", "list",
        "named sample directions, e.g. '((RD (1 0 0)) (ND (0 1 0)) (TD (0 0 1)))' (default X Y Z)"
    );
    argList::addOption("key", "word", "colour key: hsv (white centre, default) | tsl");
    argList::addOption("symmetry", "word", "cubic (default) | hexagonal");
    argList::addOption("nPhases", "label", "number of order parameters (default: N_Phases from constant/GrainProperties)");
    argList::addOption("satPower", "scalar", "saturation exponent of the hsv key (default 1)");
    argList::addBoolOption("blend", "blend the colours of all phases by phase fraction instead of taking the dominant phase");
    argList::addOption("bufferMaterial", "label", "material index of the non-crystalline buffer order parameter(s) (Phase2Material)");
    argList::addOption("bufferColour", "vector", "r g b in [0,1] painted on buffer cells (default (0.6 0.6 0.6))");

    #include "setRootCase.H"
    #include "createTime.H"
    instantList timeDirs = timeSelector::select0(runTime, args);
    #include "createMesh.H"

    List<Tuple2<word, vector>> directions;
    if (args.optionFound("directions"))
    {
        IStringStream is(args["directions"]);
        is >> directions;
    }
    else
    {
        directions.setSize(3);
        directions[0] = Tuple2<word, vector>("X", vector(1, 0, 0));
        directions[1] = Tuple2<word, vector>("Y", vector(0, 1, 0));
        directions[2] = Tuple2<word, vector>("Z", vector(0, 0, 1));
    }

    const word keyName = args.optionFound("key") ? word(args["key"]) : word("hsv");
    const word symName = args.optionFound("symmetry") ? word(args["symmetry"]) : word("cubic");
    const scalar satPower = args.optionFound("satPower") ? readScalar(IStringStream(args["satPower"])()) : 1.0;
    const bool blend = args.optionFound("blend");

    const ebsd::Symmetry sym = ebsd::symmetryFromName(symName);
    const ebsd::IpfKey key = ebsd::ipfKeyFromName(keyName);

    const label bufferMaterial =
        args.optionFound("bufferMaterial") ? readLabel(IStringStream(args["bufferMaterial"])()) : -1;
    vector bufferColour(0.6, 0.6, 0.6);
    if (args.optionFound("bufferColour"))
    {
        IStringStream is(args["bufferColour"]);
        is >> bufferColour;
    }

    label nPhases = -1;
    labelList phaseMaterial;
    if (args.optionFound("nPhases") && bufferMaterial < 0)
    {
        nPhases = readLabel(IStringStream(args["nPhases"])());
    }
    else
    {
        IOdictionary GrainProperties
        (
            IOobject
            (
                "GrainProperties",
                runTime.constant(),
                mesh,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            )
        );
        nPhases = readLabel(GrainProperties.lookup("N_Phases"));
        if (bufferMaterial >= 0)
        {
            GrainProperties.lookup("Phase2Material") >> phaseMaterial;
            if (phaseMaterial.size() != nPhases)
            {
                FatalErrorInFunction
                    << "Phase2Material has " << phaseMaterial.size() << " entries, N_Phases is "
                    << nPhases << exit(FatalError);
            }
        }
    }
    // crystalline (IPF-coloured) order parameters
    boolList crystalline(nPhases, true);
    label nBufferPhases = 0;
    if (bufferMaterial >= 0)
    {
        forAll(phaseMaterial, ph)
        {
            if (phaseMaterial[ph] == bufferMaterial) { crystalline[ph] = false; nBufferPhases++; }
        }
    }
    Info<< "Number of order parameters: " << nPhases;
    if (nBufferPhases)
    {
        Info<< " (" << nBufferPhases << " buffer, Material." << bufferMaterial << ")";
    }
    Info<< nl << "Colour key: " << keyName << ", symmetry: " << symName
        << (blend ? ", phase-fraction blended" : ", dominant phase") << nl << endl;

    const word patchType = zeroGradientFvPatchScalarField::typeName;

    forAll(timeDirs, timei)
    {
        runTime.setTime(timeDirs[timei], timei);
        Info<< "Time = " << runTime.timeName() << endl;
        mesh.readUpdate();

        PtrList<volScalarField> n(nPhases), rot1(nPhases), rot2(nPhases), rot3(nPhases);
        bool ok = true;
        for (label ph = 0; ph < nPhases; ph++)
        {
            const word nName("n." + name(ph));
            const word r1Name("ROT.1." + name(ph));
            if (!isFile(runTime.timePath()/nName) || !isFile(runTime.timePath()/r1Name))
            {
                Info<< "    " << nName << " or " << r1Name << " not found at this time, skipping" << endl;
                ok = false;
                break;
            }
            n.set(ph, new volScalarField(IOobject(nName, runTime.timeName(), mesh, IOobject::MUST_READ, IOobject::NO_WRITE), mesh));
            rot1.set(ph, new volScalarField(IOobject(r1Name, runTime.timeName(), mesh, IOobject::MUST_READ, IOobject::NO_WRITE), mesh));
            rot2.set(ph, new volScalarField(IOobject(word("ROT.2." + name(ph)), runTime.timeName(), mesh, IOobject::MUST_READ, IOobject::NO_WRITE), mesh));
            rot3.set(ph, new volScalarField(IOobject(word("ROT.3." + name(ph)), runTime.timeName(), mesh, IOobject::MUST_READ, IOobject::NO_WRITE), mesh));
        }
        if (!ok) continue;

        // dominant phase per cell
        volScalarField dominantPhase
        (
            IOobject("dominantPhase", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::NO_WRITE),
            mesh, dimensionedScalar("dominantPhase", dimless, 0.0), patchType
        );
        forAll(mesh.C(), celli)
        {
            label best = 0;
            scalar bestVal = -1.0;
            for (label ph = 0; ph < nPhases; ph++)
            {
                if (n[ph][celli] > bestVal) { bestVal = n[ph][celli]; best = ph; }
            }
            dominantPhase[celli] = best;
        }
        dominantPhase.correctBoundaryConditions();
        dominantPhase.write();

        if (nBufferPhases)
        {
            volScalarField crystalFraction
            (
                IOobject("crystalFraction", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::NO_WRITE),
                mesh, dimensionedScalar("crystalFraction", dimless, 0.0), patchType
            );
            forAll(mesh.C(), celli)
            {
                scalar f = 0.0;
                for (label ph = 0; ph < nPhases; ph++)
                {
                    if (crystalline[ph]) f += max(n[ph][celli], 0.0);
                }
                crystalFraction[celli] = min(f, 1.0);
            }
            crystalFraction.correctBoundaryConditions();
            crystalFraction.write();
        }

        forAll(directions, di)
        {
            const word fieldName("IPF_" + directions[di].first());
            const ebsd::Vec3 dS = ebsd::Vec3
            (
                directions[di].second().x(), directions[di].second().y(), directions[di].second().z()
            ).normalised();

            volVectorField ipf
            (
                IOobject(fieldName, runTime.timeName(), mesh, IOobject::NO_READ, IOobject::NO_WRITE),
                mesh, dimensionedVector("IPF", dimless, vector::zero), patchType
            );

            forAll(mesh.C(), celli)
            {
                if (blend)
                {
                    double r = 0, g = 0, b = 0, w = 0;
                    for (label ph = 0; ph < nPhases; ph++)
                    {
                        const scalar f = max(n[ph][celli], 0.0);
                        if (f < SMALL) continue;
                        if (!crystalline[ph])
                        {
                            r += f*bufferColour.x(); g += f*bufferColour.y(); b += f*bufferColour.z(); w += f;
                            continue;
                        }
                        const ebsd::Mat3 R = ebsd::codeAnglesToMatrix(rot1[ph][celli], rot2[ph][celli], rot3[ph][celli]);
                        const ebsd::RGB c = ebsd::ipfColour(R, dS, sym, key, satPower);
                        r += f*c.r; g += f*c.g; b += f*c.b; w += f;
                    }
                    if (w > SMALL) { r /= w; g /= w; b /= w; }
                    ipf[celli] = vector(r, g, b);
                }
                else
                {
                    const label ph = label(dominantPhase[celli]);
                    if (!crystalline[ph])
                    {
                        ipf[celli] = bufferColour;
                        continue;
                    }
                    const ebsd::Mat3 R = ebsd::codeAnglesToMatrix(rot1[ph][celli], rot2[ph][celli], rot3[ph][celli]);
                    const ebsd::RGB c = ebsd::ipfColour(R, dS, sym, key, satPower);
                    ipf[celli] = vector(c.r, c.g, c.b);
                }
            }
            ipf.correctBoundaryConditions();
            ipf.write();
            Info<< "    wrote " << fieldName << " (sample direction " << directions[di].second() << ")" << endl;
        }
    }

    Info<< nl << "End" << nl << endl;
    return 0;
}


// ************************************************************************* //
