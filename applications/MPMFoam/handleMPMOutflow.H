/*---------------------------------------------------------------------------*\
  handleMPMOutflow.H - STATISTICS ONLY (No deletion needed!)
  ===========================================================
  
  Since we're using active flags instead of deleting particles,
  this file just reports statistics.
  
\*---------------------------------------------------------------------------*/

// Count active and inactive particles
label nActive = 0;
label nInactive = 0;

forAllConstIter(materialPointCloud, particles, iter)
{
    if (iter().active())
        nActive++;
    else
        nInactive++;
}

// Report statistics
if (Pstream::parRun())
{
    label nActiveGlobal = returnReduce(nActive, sumOp<label>());
    label nInactiveGlobal = returnReduce(nInactive, sumOp<label>());
    
    if (Pstream::master())
    {
        Info<< "Particle status:" << endl;
        Info<< "  Active:   " << nActiveGlobal << endl;
        Info<< "  Inactive: " << nInactiveGlobal << endl;
    }
}
else
{
    Info<< "Particle status:" << endl;
    Info<< "  Active:   " << nActive << endl;
    Info<< "  Inactive: " << nInactive << endl;
}

// ************************************************************************* //