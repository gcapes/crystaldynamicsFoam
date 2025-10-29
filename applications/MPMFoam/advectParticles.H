// Advect particles through mesh

// scalar dt = runTime.deltaT().value();

// Create tracking data
particle::trackingData td(particles);

forAllIter(materialPointCloud, particles, iter)
{
    MaterialPointParcel& p = iter();
    p.move(particles, td, dt);
}