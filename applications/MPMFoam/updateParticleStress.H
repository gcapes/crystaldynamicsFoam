// Update particle stress using elastic constitutive model

// scalar dt = runTime.deltaT().value();

forAllIter(materialPointCloud, particles, iter)
{
    MaterialPointParcel& p = iter();
    
    label celli = p.cell();
    if (celli < 0) continue;
    
    // Get velocity gradient at particle location
    tensor L = gridVelocityGradient[celli];
    
    // Update deformation gradient
    // F_new = (I + dt*L) * F_old
    tensor Fold = p.F_;
    p.F_ = (tensor::I + dt*L) & Fold;
    
    // Elastic constitutive model
    // For now: Neo-Hookean or linear elastic
    
    // Material properties (should come from Material class)
    scalar E = 200e9;  // Young's modulus (Pa)
    scalar nu = 0.3;   // Poisson's ratio
    
    // Lame parameters
    scalar lambda = E*nu / ((1.0 + nu)*(1.0 - 2.0*nu));
    scalar mu = E / (2.0*(1.0 + nu));
    
    // Green strain: E = 0.5*(F^T*F - I)
    tensor C = p.F_.T() & p.F_;
    tensor greenStrain = 0.5 * (C - tensor::I);
    
    // Second Piola-Kirchhoff stress
    scalar trE = tr(greenStrain);
    tensor S = lambda*trE*tensor::I + 2.0*mu*greenStrain;
    
    // Cauchy stress: sigma = (1/J) * F * S * F^T
    scalar J = det(p.F_);
    p.stress_ = (1.0/J) * (p.F_ & S & p.F_.T());
}
