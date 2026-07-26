
{
    // Compute the discrete continuity error
    // This should be zero if rho satisfies mass conservation exactly.
    // Since rho = rho_ref/J is algebraic, it won't be zero discretely.
    volScalarField contErr
    (
        "contErr",
        fvc::ddt(rho) + fvc::div(rhoPhi)
    );

    Sum_Phases *= 0.0;
    for (label count = 0; count < N_Ori; count++)
    {
        volScalarField &ni = PopBal[count];

        solve(
            fvm::ddt(rho, ni) + fvm::div(rhoPhi, ni)
            - fvm::SuSp(contErr, ni)  // <-- Continuity error correction
        );
        ni.correctBoundaryConditions();
        Sum_Phases += ni;

        Info << "Max phase field(" << count << ") = " << gMax(ni) << endl;
        Info << "Min phase field(" << count << ") = " << gMin(ni) << endl;
    }
    Sum_Phases.correctBoundaryConditions();

    // Project any out-of-bounds phases back onto the Gibbs simplex
    for (label count = 0; count < N_Ori; count++)
    {
        volScalarField &ni = PopBal[count];

        ni = max(min(ni, 1.0), 0.0);
        volScalarField temp_holder = ni;

        forAll(mesh.C(), celli)
        {
            if (Sum_Phases[celli] > 0.0)
            {
                temp_holder[celli] *= (1.0 / Sum_Phases[celli]);
            }
        }

        forAll(Sum_Phases.boundaryField(), patchI)
        {
            if (!Sum_Phases.boundaryField()[patchI].coupled()
                && Sum_Phases.boundaryField()[patchI].size() > 0)
            {
                forAll(Sum_Phases.boundaryField()[patchI], faceI)
                {
                    if (Sum_Phases.boundaryField()[patchI][faceI] > 0.0)
                    {
                        temp_holder.boundaryFieldRef()[patchI][faceI] *= 
                            (1.0 / Sum_Phases.boundaryField()[patchI][faceI]);
                    }
                }
            }
        }

        ni = temp_holder;
    }
}