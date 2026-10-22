// Pre-populate C_rotated when the stress solve is bypassed
for (label h = 0; h < N_Ori; h++)
{
    // Use cell 0 rotation angles (uniform per phase).
    const double p1 = Rot1[h][0];
    const double p2 = Rot2[h][0];
    const double p3 = Rot3[h][0];
    const double c1 = Foam::cos(p1), s1 = Foam::sin(p1);
    const double c2 = Foam::cos(p2), s2 = Foam::sin(p2);
    const double c3 = Foam::cos(p3), s3 = Foam::sin(p3);

    Fastor::Tensor<double,3,3> R(0.0);
    R(0,0) =  c2*c3;
    R(0,1) =  s1*s2*c3 - c1*s3;
    R(0,2) =  c1*s2*c3 + s1*s3;
    R(1,0) =  c2*s3;
    R(1,1) =  s1*s2*s3 + c1*c3;
    R(1,2) =  c1*s2*s3 - s1*c3;
    R(2,0) = -s2;
    R(2,1) =  s1*c2;
    R(2,2) =  c1*c2;

    Fastor::Tensor<double,3,3,3,3> C_rot(0.0);
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
            for (int k = 0; k < 3; k++)
                for (int l = 0; l < 3; l++)
                {
                    double val = 0.0;
                    for (int m = 0; m < 3; m++)
                        for (int n = 0; n < 3; n++)
                            for (int p = 0; p < 3; p++)
                                for (int q = 0; q < 3; q++)
                                    val += R(i,m)*R(j,n)*R(k,p)*R(l,q)
                                            * Cprime2[Phase2Material[h]](m,n,p,q);
                    C_rot(i,j,k,l) = val;
                }

    for (label cell = 0; cell < mesh.nCells(); cell++)
        C_rotated[h][cell] = C_rot;
}

// PerPhaseSigma: set sigma[h] once from sigmaPrescribed; StressUpdate.H never
// touches it again.
if (perPhaseSigma)
{
    Info<< "PerPhaseSigma=true: initialising sigma.<i> from sigmaPrescribed\n" << endl;
    for (label h = 0; h < N_Ori; h++)
    {
        sigma[h] == dimensionedTensor(
            "sigmaPrescribed",
            sigma[h].dimensions(),
            sigmaPrescribed[h]
        );
        sigma[h].correctBoundaryConditions();
    }
}
