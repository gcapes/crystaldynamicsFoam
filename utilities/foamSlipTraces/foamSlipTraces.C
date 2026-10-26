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
    foamSlipTraces

Description
    Slip activity from n.i, ROT.k.i and gamma_alphas.i.s in the form of HR-DIC slip-band maps:
    slip totals, plastic strain, dominant system and plane, slip-trace angles, per-grain table.

    Usage
        foamSlipTraces [-latestTime | -time <ranges>] [-normal "(0 0 1)"]
                       [-loading "(0 1 0)"] [-noTable] [-parallel]

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "timeSelector.H"
#include "IStringStream.H"
#include "OFstream.H"
#include "symmTensor.H"

#include <vector>
#include <cmath>

using namespace Foam;

namespace
{

//- crystaldynamicsFoam orientation matrix R = Rz(p3) Ry(p2) Rx(p1), v_sample = R v_crystal
tensor codeR(const scalar p1, const scalar p2, const scalar p3)
{
    const scalar c1 = Foam::cos(p1), s1 = Foam::sin(p1);
    const scalar c2 = Foam::cos(p2), s2 = Foam::sin(p2);
    const scalar c3 = Foam::cos(p3), s3 = Foam::sin(p3);
    return tensor
    (
        c2*c3,  s1*s2*c3 - c1*s3,  c1*s2*c3 + s1*s3,
        c2*s3,  s1*s2*s3 + c1*c3,  c1*s2*s3 - s1*c3,
        -s2,    s1*c2,             c1*c2
    );
}

struct SlipMaterial
{
    label nSlips = 0;
    List<vector> n;      // plane normals (crystal frame, unit)
    List<vector> m;      // slip directions (crystal frame, unit)
    labelList plane;     // plane index of every system (systems grouped by +/- normal)
    label nPlanes = 0;
};

//- angle of an in-plane direction t (deg, in [0, 180)) in the basis (e1, e2)
scalar traceAngle(const vector& t, const vector& e1, const vector& e2)
{
    scalar a = Foam::atan2(t & e2, t & e1)*180.0/constant::mathematical::pi;
    while (a < 0) a += 180.0;
    while (a >= 180.0) a -= 180.0;
    return a;
}

} // anonymous namespace


int main(int argc, char *argv[])
{
    timeSelector::addOptions();
    argList::addNote
    (
        "Slip-trace and slip-activity fields from the crystaldynamicsFoam order parameters,\n"
        "orientations and accumulated slips, for comparison with HR-DIC slip-band maps."
    );
    argList::addOption("normal", "vector", "normal of the observation plane (default (0 0 1))");
    argList::addOption("loading", "vector", "uniaxial loading direction for the Schmid factors in the table (default (0 1 0))");
    argList::addOption("minFraction", "scalar", "ignore phases below this fraction in a cell (default 1e-3)");
    argList::addBoolOption("noTable", "do not write the per-grain table");

    #include "setRootCase.H"
    #include "createTime.H"
    instantList timeDirs = timeSelector::select0(runTime, args);
    #include "createMesh.H"

    vector N(0, 0, 1);
    if (args.optionFound("normal")) { IStringStream is(args["normal"]); is >> N; }
    N /= mag(N);
    vector L(0, 1, 0);
    if (args.optionFound("loading")) { IStringStream is(args["loading"]); is >> L; }
    L /= mag(L);
    const scalar minFraction =
        args.optionFound("minFraction") ? readScalar(IStringStream(args["minFraction"])()) : 1e-3;
    const bool writeTable = !args.optionFound("noTable");

    // in-plane basis: e1 = the mesh axis most nearly perpendicular to N, e2 = N x e1
    vector e1(1, 0, 0);
    {
        const vector ax[3] = {vector(1,0,0), vector(0,1,0), vector(0,0,1)};
        scalar best = GREAT;
        for (int k = 0; k < 3; k++)
        {
            if (mag(ax[k] & N) < best) { best = mag(ax[k] & N); e1 = ax[k]; }
        }
        e1 -= (e1 & N)*N; e1 /= mag(e1);
    }
    const vector e2 = N ^ e1;

    // --- GrainProperties
    IOdictionary GrainProperties
    (
        IOobject("GrainProperties", runTime.constant(), mesh, IOobject::MUST_READ, IOobject::NO_WRITE)
    );
    const label nPhases = readLabel(GrainProperties.lookup("N_Phases"));
    labelList phaseMaterial;
    GrainProperties.lookup("Phase2Material") >> phaseMaterial;
    const label nMats = readLabel(GrainProperties.lookup("N_Materials"));
    if (phaseMaterial.size() != nPhases)
    {
        FatalErrorInFunction << "Phase2Material has " << phaseMaterial.size()
            << " entries, N_Phases is " << nPhases << exit(FatalError);
    }

    List<SlipMaterial> mats(nMats);
    label maxSlips = 0;
    for (label mi = 0; mi < nMats; mi++)
    {
        const dictionary& md = GrainProperties.subDict("Material." + name(mi));
        SlipMaterial& sm = mats[mi];
        sm.nSlips = readLabel(md.lookup("N_slips"));
        sm.n.setSize(sm.nSlips); sm.m.setSize(sm.nSlips); sm.plane.setSize(sm.nSlips, -1);
        for (label s = 0; s < sm.nSlips; s++)
        {
            vector nv(md.lookup("n" + name(s)));
            vector mv(md.lookup("m" + name(s)));
            sm.n[s] = nv/mag(nv); sm.m[s] = mv/mag(mv);
        }
        // group the systems by plane (normals equal up to sign)
        for (label s = 0; s < sm.nSlips; s++)
        {
            if (sm.plane[s] >= 0) continue;
            sm.plane[s] = sm.nPlanes;
            for (label r = s + 1; r < sm.nSlips; r++)
            {
                if (sm.plane[r] < 0 && mag(mag(sm.n[s] & sm.n[r]) - 1.0) < 1e-6) sm.plane[r] = sm.nPlanes;
            }
            sm.nPlanes++;
        }
        maxSlips = max(maxSlips, sm.nSlips);
        Info<< "Material." << mi << ": " << sm.nSlips << " slip systems on " << sm.nPlanes << " planes" << endl;
    }
    Info<< "Observation normal " << N << ", in-plane axes e1 " << e1 << " e2 " << e2
        << "; loading direction " << L << nl << endl;

    const word patchType = zeroGradientFvPatchScalarField::typeName;

    forAll(timeDirs, timei)
    {
        runTime.setTime(timeDirs[timei], timei);
        Info<< "Time = " << runTime.timeName() << endl;
        mesh.readUpdate();

        // --- order parameters and the dominant phase per cell
        PtrList<volScalarField> n(nPhases);
        bool ok = true;
        for (label h = 0; h < nPhases; h++)
        {
            const word nName("n." + name(h));
            if (!isFile(runTime.timePath()/nName)) { Info<< "    " << nName << " not found, skipping" << endl; ok = false; break; }
            n.set(h, new volScalarField(IOobject(nName, runTime.timeName(), mesh, IOobject::MUST_READ, IOobject::NO_WRITE), mesh));
        }
        if (!ok) continue;

        labelList dominant(mesh.nCells(), 0);
        forAll(dominant, celli)
        {
            scalar best = -1;
            for (label h = 0; h < nPhases; h++)
            {
                if (n[h][celli] > best) { best = n[h][celli]; dominant[celli] = h; }
            }
        }

        // --- output fields
        volScalarField slipTotal(IOobject("slipTotal", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::NO_WRITE), mesh, dimensionedScalar("0", dimless, 0.0), patchType);
        volSymmTensorField plasticStrain(IOobject("plasticStrain", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::NO_WRITE), mesh, dimensionedSymmTensor("0", dimless, symmTensor::zero), patchType);
        volScalarField shearStrainInPlane(IOobject("shearStrainInPlane", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::NO_WRITE), mesh, dimensionedScalar("0", dimless, 0.0), patchType);
        volScalarField slipDominantSystem(IOobject("slipDominantSystem", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::NO_WRITE), mesh, dimensionedScalar("0", dimless, -1.0), patchType);
        volScalarField slipDominantPlane(IOobject("slipDominantPlane", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::NO_WRITE), mesh, dimensionedScalar("0", dimless, -1.0), patchType);
        volScalarField slipDominantFraction(IOobject("slipDominantFraction", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::NO_WRITE), mesh, dimensionedScalar("0", dimless, 0.0), patchType);
        volScalarField slipTraceAngle(IOobject("slipTraceAngle", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::NO_WRITE), mesh, dimensionedScalar("0", dimless, -1.0), patchType);
        volScalarField slipTraceAngle2(IOobject("slipTraceAngle2", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::NO_WRITE), mesh, dimensionedScalar("0", dimless, -1.0), patchType);
        volVectorField slipTraceDir(IOobject("slipTraceDir", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::NO_WRITE), mesh, dimensionedVector("0", dimless, vector::zero), patchType);

        // --- per-grain accumulators (cells where the phase is dominant)
        List<label> gCount(nPhases, 0);
        List<scalarList> gSumA(nPhases), gSumSchmid(nPhases);
        List<scalar> gCos2(nPhases, 0.0), gSin2(nPhases, 0.0), gSumEff(nPhases, 0.0);
        List<symmTensor> gSumEps(nPhases, symmTensor::zero);
        List<scalarList> gPlaneA(nPhases);
        for (label h = 0; h < nPhases; h++)
        {
            const SlipMaterial& sm = mats[phaseMaterial[h]];
            gSumA[h].setSize(sm.nSlips, 0.0); gSumSchmid[h].setSize(sm.nSlips, 0.0);
            gPlaneA[h].setSize(max(sm.nPlanes, 1), 0.0);
        }

        // --- phase by phase
        for (label h = 0; h < nPhases; h++)
        {
            const SlipMaterial& sm = mats[phaseMaterial[h]];
            if (sm.nSlips == 0) continue;

            const word r1("ROT.1." + name(h)), r2("ROT.2." + name(h)), r3("ROT.3." + name(h));
            if (!isFile(runTime.timePath()/r1)) { Info<< "    " << r1 << " not found, phase skipped" << endl; continue; }
            volScalarField rot1(IOobject(r1, runTime.timeName(), mesh, IOobject::MUST_READ, IOobject::NO_WRITE), mesh);
            volScalarField rot2(IOobject(r2, runTime.timeName(), mesh, IOobject::MUST_READ, IOobject::NO_WRITE), mesh);
            volScalarField rot3(IOobject(r3, runTime.timeName(), mesh, IOobject::MUST_READ, IOobject::NO_WRITE), mesh);
            PtrList<volScalarField> gam(sm.nSlips);
            bool okG = true;
            for (label s = 0; s < sm.nSlips; s++)
            {
                const word gName("gamma_alphas." + name(h) + "." + name(s));
                if (!isFile(runTime.timePath()/gName)) { Info<< "    " << gName << " not found, phase skipped" << endl; okG = false; break; }
                gam.set(s, new volScalarField(IOobject(gName, runTime.timeName(), mesh, IOobject::MUST_READ, IOobject::NO_WRITE), mesh));
            }
            if (!okG) continue;

            scalarList a(sm.nSlips);
            List<vector> nGl(sm.nSlips), mGl(sm.nSlips);

            forAll(dominant, celli)
            {
                const scalar f = n[h][celli];
                const bool isDom = (dominant[celli] == h);
                if (f < minFraction && !isDom) continue;

                const tensor R = codeR(rot1[celli], rot2[celli], rot3[celli]);
                scalar total = 0.0;
                symmTensor eps(symmTensor::zero);
                for (label s = 0; s < sm.nSlips; s++)
                {
                    nGl[s] = R & sm.n[s];
                    mGl[s] = R & sm.m[s];
                    const scalar g = gam[s][celli];
                    a[s] = mag(g);
                    total += a[s];
                    eps += g*symm(mGl[s]*nGl[s]);
                }
                const scalar w = max(f, 0.0);
                slipTotal[celli] += w*total;
                plasticStrain[celli] += w*eps;

                if (!isDom) continue;

                // dominant system of the dominant phase, and the runner-up on another plane
                label sDom = -1, sSec = -1;
                scalar aDom = 0.0, aSec = 0.0;
                for (label s = 0; s < sm.nSlips; s++)
                {
                    if (a[s] > aDom) { aDom = a[s]; sDom = s; }
                }
                if (sDom >= 0)
                {
                    for (label s = 0; s < sm.nSlips; s++)
                    {
                        if (sm.plane[s] != sm.plane[sDom] && a[s] > aSec) { aSec = a[s]; sSec = s; }
                    }
                }
                if (sDom >= 0 && aDom > 0.0)
                {
                    slipDominantSystem[celli] = sDom;
                    slipDominantPlane[celli] = sm.plane[sDom];
                    slipDominantFraction[celli] = aDom/max(total, VSMALL);
                    vector t = N ^ nGl[sDom];
                    if (mag(t) > SMALL)
                    {
                        t /= mag(t);
                        const scalar ang = traceAngle(t, e1, e2);
                        slipTraceAngle[celli] = ang;
                        slipTraceDir[celli] = aDom*t;
                        const scalar a2 = 2.0*ang*constant::mathematical::pi/180.0;
                        gCos2[h] += Foam::cos(a2); gSin2[h] += Foam::sin(a2);
                    }
                    if (sSec >= 0 && aSec > 0.0)
                    {
                        vector t2 = N ^ nGl[sSec];
                        if (mag(t2) > SMALL) slipTraceAngle2[celli] = traceAngle(t2/mag(t2), e1, e2);
                    }
                }

                // table accumulators
                gCount[h]++;
                for (label s = 0; s < sm.nSlips; s++)
                {
                    gSumA[h][s] += a[s];
                    gPlaneA[h][sm.plane[s]] += a[s];
                    gSumSchmid[h][s] += (mGl[s] & L)*(nGl[s] & L);
                }
                gSumEps[h] += eps;
                const scalar e11 = e1 & (eps & e1), e22 = e2 & (eps & e2), e12 = e1 & (eps & e2);
                gSumEff[h] += Foam::sqrt(sqr(0.5*(e11 - e22)) + sqr(e12));
            }
        }

        forAll(dominant, celli)
        {
            const symmTensor& eps = plasticStrain[celli];
            const scalar e11 = e1 & (eps & e1), e22 = e2 & (eps & e2), e12 = e1 & (eps & e2);
            shearStrainInPlane[celli] = Foam::sqrt(sqr(0.5*(e11 - e22)) + sqr(e12));
        }

        slipTotal.correctBoundaryConditions(); slipTotal.write();
        plasticStrain.correctBoundaryConditions(); plasticStrain.write();
        shearStrainInPlane.correctBoundaryConditions(); shearStrainInPlane.write();
        slipDominantSystem.correctBoundaryConditions(); slipDominantSystem.write();
        slipDominantPlane.correctBoundaryConditions(); slipDominantPlane.write();
        slipDominantFraction.correctBoundaryConditions(); slipDominantFraction.write();
        slipTraceAngle.correctBoundaryConditions(); slipTraceAngle.write();
        slipTraceAngle2.correctBoundaryConditions(); slipTraceAngle2.write();
        slipTraceDir.correctBoundaryConditions(); slipTraceDir.write();
        Info<< "    wrote slipTotal, plasticStrain, shearStrainInPlane, slipDominantSystem, slipDominantPlane,"
            << " slipDominantFraction, slipTraceAngle, slipTraceAngle2, slipTraceDir" << endl;

        if (!writeTable) continue;

        // --- reduce the per-grain sums over processors
        for (label h = 0; h < nPhases; h++)
        {
            reduce(gCount[h], sumOp<label>());
            reduce(gCos2[h], sumOp<scalar>()); reduce(gSin2[h], sumOp<scalar>()); reduce(gSumEff[h], sumOp<scalar>());
            reduce(gSumEps[h], sumOp<symmTensor>());
            forAll(gSumA[h], s) { reduce(gSumA[h][s], sumOp<scalar>()); reduce(gSumSchmid[h][s], sumOp<scalar>()); }
            forAll(gPlaneA[h], pl) reduce(gPlaneA[h][pl], sumOp<scalar>());
        }

        if (Pstream::master())
        {
            const fileName dir = runTime.path()/"postProcessing"/"foamSlipTraces"/runTime.timeName();
            mkDir(dir);
            OFstream os(dir/"grains.csv");
            os  << "# foamSlipTraces, time " << runTime.timeName() << ", observation normal " << N
                << ", in-plane axes e1 " << e1 << " e2 " << e2 << ", loading " << L << nl
                << "# angles in degrees, counter-clockwise from e1 seen from the normal (mesh frame);"
                << " imageFrame = (180 - meshFrame) mod 180 for a picture with y down" << nl
                << "# gammaMean_s: mean |gamma| of system s over the cells where the phase is dominant;"
                << " planeShare_p: share of plane p in the total slip; schmid_s for the loading direction" << nl;
            os  << "phase,material,cells,slipTotalMean,shearStrainInPlaneMean,"
                << "eps11,eps22,eps33,eps12,eps13,eps23,"
                << "dominantPlane,traceAngleMesh,traceAngleImage,traceAngleSpread";
            for (label s = 0; s < maxSlips; s++) os << ",gammaMean_" << s;
            for (label s = 0; s < maxSlips; s++) os << ",schmid_" << s;
            for (label pl = 0; pl < 8; pl++) os << ",planeShare_" << pl;
            os << nl;
            for (label h = 0; h < nPhases; h++)
            {
                const SlipMaterial& sm = mats[phaseMaterial[h]];
                const label c = gCount[h];
                const scalar inv = (c > 0) ? 1.0/c : 0.0;
                scalar total = 0.0;
                forAll(gSumA[h], s) total += gSumA[h][s];
                label plDom = -1; scalar plBest = -1;
                forAll(gPlaneA[h], pl) { if (gPlaneA[h][pl] > plBest) { plBest = gPlaneA[h][pl]; plDom = pl; } }
                const scalar rl = Foam::sqrt(sqr(gCos2[h]) + sqr(gSin2[h]));
                scalar angMesh = -1, angImg = -1, spread = -1;
                if (rl > 0)
                {
                    angMesh = 0.5*Foam::atan2(gSin2[h], gCos2[h])*180.0/constant::mathematical::pi;
                    while (angMesh < 0) angMesh += 180.0;
                    angImg = 180.0 - angMesh; if (angImg >= 180.0) angImg -= 180.0;
                    spread = 0.5*Foam::sqrt(-2.0*Foam::log(max(rl*inv, 1e-12)))*180.0/constant::mathematical::pi;
                }
                const symmTensor epsM = gSumEps[h]*inv;
                os  << h << "," << phaseMaterial[h] << "," << c << "," << total*inv << "," << gSumEff[h]*inv << ","
                    << epsM.xx() << "," << epsM.yy() << "," << epsM.zz() << "," << epsM.xy() << "," << epsM.xz() << "," << epsM.yz() << ","
                    << (sm.nSlips ? plDom : -1) << "," << angMesh << "," << angImg << "," << spread;
                for (label s = 0; s < maxSlips; s++) os << "," << (s < sm.nSlips ? gSumA[h][s]*inv : 0.0);
                for (label s = 0; s < maxSlips; s++) os << "," << (s < sm.nSlips ? gSumSchmid[h][s]*inv : 0.0);
                for (label pl = 0; pl < 8; pl++) os << "," << (pl < sm.nPlanes && total > 0 ? gPlaneA[h][pl]/total : 0.0);
                os << nl;
            }
            Info<< "    wrote " << dir/"grains.csv" << endl;
        }
    }

    Info<< nl << "End" << nl << endl;
    return 0;
}


// ************************************************************************* //
