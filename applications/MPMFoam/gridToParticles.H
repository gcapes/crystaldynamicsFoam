// Map grid velocity to particles (G2P)

forAllIter(materialPointCloud, particles, iter)
{
    MaterialPointParcel& p = iter();
    p.mapFromGrid(gridVelocity, gridVelocityGradient);
}
