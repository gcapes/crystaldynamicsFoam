
{
    // Velocity gradient decomposition (still useful for diagnostics)
    volTensorField TgradU(fvc::grad(U)().T());
    TgradU.correctBoundaryConditions();

    omega_dot = 0.5 * (TgradU - fvc::grad(U));
    omega_dot.correctBoundaryConditions();

    epsilon_dot = 0.5 * (TgradU + fvc::grad(U));
    epsilon_dot.correctBoundaryConditions();

    // Reset mixture properties for recomputation
    kappa *= 0.0;
    Cv *= 0.0;
    GBREF *= 0.0;

    for (label h = 0; h < N_Ori; h++)
    {
        dsigmadt[h] *= 0.0;
        Omega_P[h] *= 0.0;
        Plastic_Strain[h] *= 0.0;

        volScalarField& ni = PopBal[h];
        GBREF += ni * ni;

        // Phase-averaged field (shares info with adjacent cells)
        PopBal_av[h] = fvc::average(ni);
        PopBal_av[h].correctBoundaryConditions();

        // Mixture thermal properties
        dimensionedScalar dimkappa("dimkappa", dimensionSet(1, 1, -3, -1, 0), 1.0);
        dimensionedScalar dimcv("dimcv", dimensionSet(0, 2, -2, -1, 0), 1.0);

        kappa += PopBal[h] * MatList[Phase2Material[h]].kappa * dimkappa;
        Cv    += PopBal[h] * MatList[Phase2Material[h]].Cv * dimcv;
    }
    GBREF.correctBoundaryConditions();

    D_P *= 0.0;
    p *= 0.0;

    // Zero out rate fields before recomputation in StressUpdate
    for (label i = 0; i < N_Ori; i++)
    {
        for (label s = 0; s < MatList[Phase2Material[i]].N_Slips; s++)
        {
            g_alpha_dots[i][s] *= 0.0;
            gamma_alpha_dots[i][s] *= 0.0;
        }
    }
}