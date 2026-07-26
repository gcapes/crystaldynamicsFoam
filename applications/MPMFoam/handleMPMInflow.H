/*---------------------------------------------------------------------------*\
  handleMPMInflow.H - TWO-PASS VERSION (PARALLEL-SAFE)
  
  Pass 1: Identify particles to reactivate (read-only)
  Pass 2: Modify particles (after iteration)
\*---------------------------------------------------------------------------*/

IOdictionary inflowDict
(
    IOobject
    (
        "inflowBoundaries",
        runTime.constant(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    )
);

if (inflowDict.headerOk() && inflowDict.found("inflows"))
{
    const dictionary& inflowsDict = inflowDict.subDict("inflows");
    
    // Storage for reactivation instructions
    struct ReactivationData
    {
        MaterialPointParcel* particle;
        vector position;
        label cellI;
        vector velocity;
        label materialID;
        label grainID;
        scalar phi1, Phi, phi2;
        scalar volume;
        scalar mass;
    };
    
    DynamicList<ReactivationData> reactivationList;
    
    // Process each inflow boundary
    forAllConstIter(dictionary, inflowsDict, iter)
    {
        const word& inflowName = iter().keyword();
        const dictionary& inflowParams = iter().dict();
        
        word patchName(inflowParams.lookup("patch"));
        vector inflowVelocity(inflowParams.lookup("velocity"));
        label materialID = readLabel(inflowParams.lookup("materialID"));
        scalar targetParticlesPerCell = readScalar(inflowParams.lookup("particlesPerCell"));
        scalar minDensityFraction = inflowParams.lookupOrDefault<scalar>("minDensityFraction", 0.5);
        
        scalar phi1 = inflowParams.lookupOrDefault<scalar>("phi1", 0.0);
        scalar Phi = inflowParams.lookupOrDefault<scalar>("Phi", 0.0);
        scalar phi2 = inflowParams.lookupOrDefault<scalar>("phi2", 0.0);
        label grainID = inflowParams.lookupOrDefault<label>("grainID", 0);
        
        label patchI = mesh.boundaryMesh().findPatchID(patchName);
        
        if (patchI == -1 || isA<processorPolyPatch>(mesh.boundaryMesh()[patchI]))
        {
            continue;
        }
        
        const polyPatch& patch = mesh.boundaryMesh()[patchI];
        
        // =====================================================================
        // PASS 1: READ-ONLY - Identify what needs to be done
        // =====================================================================
        
        // Get boundary cells
        labelHashSet boundaryCells;
        forAll(patch, faceI)
        {
            boundaryCells.insert(patch.faceCells()[faceI]);
        }
        
        // Count active particles (read-only)
        labelList particlesInCell(mesh.nCells(), 0);
        forAllConstIter(materialPointCloud, particles, pIter)
        {
            if (pIter().active() && pIter().cell() >= 0)
            {
                particlesInCell[pIter().cell()]++;
            }
        }
        
        // Find cells needing particles
        DynamicList<label> cellsNeedingParticles;
        DynamicList<label> nNeededPerCell;
        
        forAllConstIter(labelHashSet, boundaryCells, cellIter)
        {
            label cellI = cellIter.key();
            scalar currentCount = particlesInCell[cellI];
            
            if (currentCount < minDensityFraction * targetParticlesPerCell)
            {
                label nNeeded = label(targetParticlesPerCell - currentCount);
                cellsNeedingParticles.append(cellI);
                nNeededPerCell.append(nNeeded);
            }
        }
        
        // Find inactive particles (read-only)
        DynamicList<MaterialPointParcel*> inactiveParticles;
        forAllIter(materialPointCloud, particles, pIter)
        {
            if (!pIter().active())
            {
                inactiveParticles.append(&pIter());
            }
        }
        
        // Build reactivation list
        const Material& material = materials[materialID];
        label particleIndex = 0;
        
        forAll(cellsNeedingParticles, i)
        {
            label cellI = cellsNeedingParticles[i];
            label nNeeded = nNeededPerCell[i];
            
            vector cellCenter = mesh.C()[cellI];
            scalar cellVolume = mesh.V()[cellI];
            scalar particleVolume = cellVolume / targetParticlesPerCell;
            
            for (label j = 0; j < nNeeded && particleIndex < inactiveParticles.size(); j++)
            {
                ReactivationData data;
                data.particle = inactiveParticles[particleIndex];
                data.position = cellCenter;
                data.cellI = cellI;
                data.velocity = inflowVelocity;
                data.materialID = materialID;
                data.grainID = grainID;
                data.phi1 = phi1;
                data.Phi = Phi;
                data.phi2 = phi2;
                data.volume = particleVolume;
                data.mass = material.rho_0 * particleVolume;
                
                reactivationList.append(data);
                particleIndex++;
            }
        }
    }
    
    // =========================================================================
    // PASS 2: MODIFY - Now apply all changes (outside iteration)
    // =========================================================================
    
    forAll(reactivationList, i)
    {
        const ReactivationData& data = reactivationList[i];
        MaterialPointParcel* p = data.particle;
        
        // Access base class members
p->setPosition(data.position);
p->setCell(data.cellI);
        
        // Set properties
        p->active() = true;
        p->volume() = data.volume;
        p->volume0() = data.volume;
        p->mass() = data.mass;
        p->velocity() = data.velocity;
        p->stress() = tensor::zero;
        
        p->materialID() = data.materialID;
        p->grainID() = data.grainID;
        
        p->phi1() = data.phi1;
        p->phi2() = data.Phi;
        p->phi3() = data.phi2;
        p->setCrystalOrientation(data.phi1, data.Phi, data.phi2);
        
        p->F() = tensor::I;
        p->Fp() = tensor::I;
        p->Fe() = tensor::I;
        p->Cp() = tensor::I;
        
        const Material& material = materials[data.materialID];
        
        if (material.N_Slips > 0)
        {
            if (p->g().size() != material.N_Slips)
            {
                p->g().setSize(material.N_Slips);
                p->gamma().setSize(material.N_Slips);
                p->gamma_dot().setSize(material.N_Slips);
                p->tau().setSize(material.N_Slips);
            }
            
            forAll(p->g(), sysI)
            {
                p->g()[sysI] = material.initial_CRSS;
                p->gamma()[sysI] = 0.0;
                p->gamma_dot()[sysI] = 0.0;
                p->tau()[sysI] = 0.0;
            }
        }
    }
    
    // Report statistics
    label nReactivated = reactivationList.size();
    
    if (Pstream::parRun())
    {
        reduce(nReactivated, sumOp<label>());
    }
    
    if (Pstream::master() && nReactivated > 0)
    {
        Info<< "  Reactivated " << nReactivated << " particles at inflow" << endl;
    }
}

// ************************************************************************* //