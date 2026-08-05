/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | CIPHERFoam - Interface Class Implementation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
\*---------------------------------------------------------------------------*/

#include "Interface.H"
#include "mathematicalConstants.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * //

Interface::Interface(word interfaceType, const volScalarField& T)
:
    IOdictionary
    (
        IOobject
        (
            "GrainProperties",
            T.time().constant(),
            T.db(),
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    ),
    mesh_(T.mesh()),
    T_(T),
    Interface_type_(interfaceType)
{
    Info << "Interface keyword is: " << interfaceType << endl;
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Interface::Assign_Energies(scalar** NRG_arr)
{
    IOdictionary PhaseFieldProperties
    (
        IOobject
        (
            "GrainProperties",
            mesh_.time().constant(),   // FIXED: was time().constant()
            mesh_,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    label N_Mats(readLabel(PhaseFieldProperties.lookup("N_Materials")));

    typedef HashTable<scalar, Interface::interfacePair, Interface::interfacePair::hash>
        PF_Energytable;

    PF_Energytable PF_Energytable_;

    PhaseFieldProperties.lookup("Energy") >> PF_Energytable_;

    for (label count_i = 0; count_i < N_Mats; count_i++)
    {
        for (label count_j = 0; count_j < N_Mats; count_j++)
        {
            word base1 = "Material.";
            word base2 = "Material.";
            word iter_i = std::to_string(min(count_i, count_j));
            word iter_j = std::to_string(max(count_i, count_j));
            word phase1 = base1.append(iter_i);
            word phase2 = base2.append(iter_j);

            PF_Energytable::const_iterator energy_iterator =
                PF_Energytable_.find(Interface::interfacePair(phase1, phase2));

            if (energy_iterator == PF_Energytable_.end())
            {
                FatalErrorInFunction
                    << "Cannot find interface " << Interface::interfacePair(phase1, phase2)
                    << " in list of Energy values"
                    << exit(FatalError);
            }

            NRG_arr[count_i][count_j] = energy_iterator();
        }
    }
}


// void Interface::Assign_Anisotropy_Strength(scalar** Aniso_strength_arr)
// {
//     IOdictionary PhaseFieldProperties
//     (
//         IOobject
//         (
//             "GrainProperties",
//             mesh_.time().constant(),   // FIXED: was time().constant()
//             mesh_,
//             IOobject::MUST_READ,
//             IOobject::NO_WRITE
//         )
//     );

//     label N_Mats(readLabel(PhaseFieldProperties.lookup("N_Materials")));

//     typedef HashTable<scalar, Interface::interfacePair, Interface::interfacePair::hash>
//         PF_Aniso_table;

//     PF_Aniso_table PF_Aniso_table_;

//     PhaseFieldProperties.lookup("Anisotropy") >> PF_Aniso_table_;

//     for (label count_i = 0; count_i < N_Mats; count_i++)
//     {
//         for (label count_j = 0; count_j < N_Mats; count_j++)
//         {
//             word base1 = "Material.";
//             word base2 = "Material.";
//             word iter_i = std::to_string(min(count_i, count_j));
//             word iter_j = std::to_string(max(count_i, count_j));
//             word phase1 = base1.append(iter_i);
//             word phase2 = base2.append(iter_j);

//             PF_Aniso_table::const_iterator energy_iterator =
//                 PF_Aniso_table_.find(Interface::interfacePair(phase1, phase2));

//             if (energy_iterator == PF_Aniso_table_.end())
//             {
//                 FatalErrorInFunction
//                     << "Cannot find interface " << Interface::interfacePair(phase1, phase2)
//                     << " in list of Anisotropy values"
//                     << exit(FatalError);
//             }

//             Aniso_strength_arr[count_i][count_j] = energy_iterator();
//         }
//     }
// }


void Interface::Assign_Mobilities(scalar** Mob_arr)
{
    IOdictionary PhaseFieldProperties
    (
        IOobject
        (
            "GrainProperties",
            mesh_.time().constant(),   // FIXED: was time().constant()
            mesh_,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    label N_Mats(readLabel(PhaseFieldProperties.lookup("N_Materials")));

    typedef HashTable<scalar, Interface::interfacePair, Interface::interfacePair::hash>
        PF_Mobtable;

    PF_Mobtable PF_Mobtable_;

    PhaseFieldProperties.lookup("Mobility") >> PF_Mobtable_;

    for (label count_i = 0; count_i < N_Mats; count_i++)
    {
        for (label count_j = 0; count_j < N_Mats; count_j++)
        {
            word base1 = "Material.";
            word base2 = "Material.";
            word iter_i = std::to_string(min(count_i, count_j));
            word iter_j = std::to_string(max(count_i, count_j));
            word phase1 = base1.append(iter_i);
            word phase2 = base2.append(iter_j);

            PF_Mobtable::const_iterator mob_iterator =
                PF_Mobtable_.find(Interface::interfacePair(phase1, phase2));

            if (mob_iterator == PF_Mobtable_.end())
            {
                FatalErrorInFunction
                    << "Cannot find interface " << Interface::interfacePair(phase1, phase2)
                    << " in list of Mobility values"
                    << exit(FatalError);
            }

            Mob_arr[count_i][count_j] = mob_iterator();
        }
    }
}


scalar Interface::aniso_energy(
    scalar* tensor_i,
    scalar* tensor_j,
    scalar* vec_i,
    scalar* vec_j,
    const scalar& phi_i,
    const scalar& phi_j,
    const scalar& int_width,
    const scalar& deltaN
)
{
    scalar B = (4.0 / int_width);
    scalar A = -(4.0 * int_width / (constant::mathematical::pi * constant::mathematical::pi));

    double ni = phi_i;

    double ni_x = vec_i[0];
    double ni_y = vec_i[1];
    double ni_z = vec_i[2];

    double ni_xx = tensor_i[0];
    double ni_xy = tensor_i[1];
    double ni_xz = tensor_i[2];
    double ni_yy = tensor_i[4];
    double ni_yz = tensor_i[5];
    double ni_zz = tensor_i[8];

    double nj = phi_j;

    double nj_x = vec_j[0];
    double nj_y = vec_j[1];
    double nj_z = vec_j[2];

    double nj_xx = tensor_j[0];
    double nj_xy = tensor_j[1];
    double nj_xz = tensor_j[2];
    double nj_yy = tensor_j[4];
    double nj_yz = tensor_j[5];
    double nj_zz = tensor_j[8];

    scalar Aniso_term = 0.0;

    // df/dpi - df/dpj
    Aniso_term = (
        -B*ni*(-3 + (4*(Foam::pow(ni_z,4.0) + Foam::pow(ni_y,4.0) + Foam::pow(ni_x,4.0)))
            / ((Foam::pow((Foam::pow(ni_z,2.0) + Foam::pow(ni_y,2.0) + Foam::pow(ni_x,2.0)),2.0)) + deltaN))
        + A*ni_zz*(-3 + (4*(Foam::pow(ni_z,4.0) + Foam::pow(ni_y,4.0) + Foam::pow(ni_x,4.0)))
            / ((Foam::pow((Foam::pow(ni_z,2.0) + Foam::pow(ni_y,2.0) + Foam::pow(ni_x,2.0)),2.0)) + deltaN))
        + A*ni_yy*(-3 + (4*(Foam::pow(ni_z,4.0) + Foam::pow(ni_y,4.0) + Foam::pow(ni_x,4.0)))
            / ((Foam::pow((Foam::pow(ni_z,2.0) + Foam::pow(ni_y,2.0) + Foam::pow(ni_x,2.0)),2.0)) + deltaN))
        + (4*A*ni_z*(-4*(Foam::pow(ni_z,4.0) + Foam::pow(ni_y,4.0) + Foam::pow(ni_x,4.0))
            * (ni_z*ni_zz + ni_y*ni_yz + ni_x*ni_xz)
            + 4*(Foam::pow(ni_z,2.0) + Foam::pow(ni_y,2.0) + Foam::pow(ni_x,2.0))
            * (Foam::pow(ni_z,3.0)*ni_zz + Foam::pow(ni_y,3.0)*ni_yz + Foam::pow(ni_x,3.0)*ni_xz)))
            / (Foam::pow((Foam::pow(ni_z,2.0) + Foam::pow(ni_y,2.0) + Foam::pow(ni_x,2.0)),3.0) + deltaN)
        + (4*A*ni_y*(-4*(Foam::pow(ni_z,4.0) + Foam::pow(ni_y,4.0) + Foam::pow(ni_x,4.0))
            * (ni_z*ni_yz + ni_y*ni_yy + ni_x*ni_xy)
            + 4*(Foam::pow(ni_z,2.0) + Foam::pow(ni_y,2.0) + Foam::pow(ni_x,2.0))
            * (Foam::pow(ni_z,3.0)*ni_yz + Foam::pow(ni_y,3.0)*ni_yy + Foam::pow(ni_x,3.0)*ni_xy)))
            / (Foam::pow((Foam::pow(ni_z,2.0) + Foam::pow(ni_y,2.0) + Foam::pow(ni_x,2.0)),3.0) + deltaN)
        + A*(-3 + (4*(Foam::pow(ni_z,4.0) + Foam::pow(ni_y,4.0) + Foam::pow(ni_x,4.0)))
            / ((Foam::pow((Foam::pow(ni_z,2.0) + Foam::pow(ni_y,2.0) + Foam::pow(ni_x,2.0)),2.0)) + deltaN))*ni_xx
        + (4*A*ni_x*(-4*(Foam::pow(ni_z,4.0) + Foam::pow(ni_y,4.0) + Foam::pow(ni_x,4.0))
            * (ni_z*ni_xz + ni_y*ni_xy + ni_x*ni_xx)
            + 4*(Foam::pow(ni_z,2.0) + Foam::pow(ni_y,2.0) + Foam::pow(ni_x,2.0))
            * (Foam::pow(ni_z,3.0)*ni_xz + Foam::pow(ni_y,3.0)*ni_xy + Foam::pow(ni_x,3.0)*ni_xx)))
            / (Foam::pow((Foam::pow(ni_z,2.0) + Foam::pow(ni_y,2.0) + Foam::pow(ni_x,2.0)),3.0) + deltaN)
        + (-B*nj*Foam::pow((Foam::pow(ni_z,2.0) + Foam::pow(ni_y,2.0) + Foam::pow(ni_x,2.0)),2.0)
            * (3*Foam::pow((Foam::pow(ni_z,2.0) + Foam::pow(ni_y,2.0) + Foam::pow(ni_x,2.0)),2.0)
            - 4*(Foam::pow(ni_z,4.0) + Foam::pow(ni_y,4.0) + Foam::pow(ni_x,4.0)))
            + A*nj_zz*Foam::pow((Foam::pow(ni_z,2.0) + Foam::pow(ni_y,2.0) + Foam::pow(ni_x,2.0)),2.0)
            * (3*Foam::pow((Foam::pow(ni_z,2.0) + Foam::pow(ni_y,2.0) + Foam::pow(ni_x,2.0)),2.0)
            - 4*(Foam::pow(ni_z,4.0) + Foam::pow(ni_y,4.0) + Foam::pow(ni_x,4.0)))
            + A*nj_yy*Foam::pow((Foam::pow(ni_z,2.0) + Foam::pow(ni_y,2.0) + Foam::pow(ni_x,2.0)),2.0)
            * (3*Foam::pow((Foam::pow(ni_z,2.0) + Foam::pow(ni_y,2.0) + Foam::pow(ni_x,2.0)),2.0)
            - 4*(Foam::pow(ni_z,4.0) + Foam::pow(ni_y,4.0) + Foam::pow(ni_x,4.0)))
            + 4*A*nj_z*(Foam::pow(ni_z,2.0) + Foam::pow(ni_y,2.0) + Foam::pow(ni_x,2.0))
            * (4*(Foam::pow(ni_z,4.0) + Foam::pow(ni_y,4.0) + Foam::pow(ni_x,4.0))
            * (ni_z*ni_zz + ni_y*ni_yz + ni_x*ni_xz)
            - 4*(Foam::pow(ni_z,2.0) + Foam::pow(ni_y,2.0) + Foam::pow(ni_x,2.0))
            * (Foam::pow(ni_z,3.0)*ni_zz + Foam::pow(ni_y,3.0)*ni_yz + Foam::pow(ni_x,3.0)*ni_xz))
            - 16*(B*ni*nj + A*(ni_z*nj_z + ni_y*nj_y + ni_x*nj_x))
            * (3*Foam::pow(ni_z,2.0)*ni_zz*Foam::pow((Foam::pow(ni_z,2.0) + Foam::pow(ni_y,2.0) + Foam::pow(ni_x,2.0)),2.0)
            - ni_zz*(Foam::pow(ni_z,2.0) + Foam::pow(ni_y,2.0) + Foam::pow(ni_x,2.0))
            * (Foam::pow(ni_z,4.0) + Foam::pow(ni_y,4.0) + Foam::pow(ni_x,4.0))
            - 4*Foam::pow(ni_z,3.0)*(Foam::pow(ni_z,2.0) + Foam::pow(ni_y,2.0) + Foam::pow(ni_x,2.0))
            * (ni_z*ni_zz + ni_y*ni_yz + ni_x*ni_xz)
            + 6*ni_z*(Foam::pow(ni_z,4.0) + Foam::pow(ni_y,4.0) + Foam::pow(ni_x,4.0))
            * (ni_z*ni_zz + ni_y*ni_yz + ni_x*ni_xz)
            - 4*ni_z*(Foam::pow(ni_z,2.0) + Foam::pow(ni_y,2.0) + Foam::pow(ni_x,2.0))
            * (Foam::pow(ni_z,3.0)*ni_zz + Foam::pow(ni_y,3.0)*ni_yz + Foam::pow(ni_x,3.0)*ni_xz))
            - 16*ni_z*(Foam::pow(ni_z,2.0) + Foam::pow(ni_y,2.0) + Foam::pow(ni_x,2.0))
            * (-Foam::pow(ni_y,4.0) - Foam::pow(ni_x,4.0) + Foam::pow(ni_z,2.0)*(Foam::pow(ni_y,2.0) + Foam::pow(ni_x,2.0)))
            * (B*nj*ni_z + B*ni*nj_z + A*(nj_z*ni_zz + ni_z*nj_zz + nj_y*ni_yz + ni_y*nj_yz + nj_x*ni_xz + ni_x*nj_xz))
            + 4*A*nj_y*(Foam::pow(ni_z,2.0) + Foam::pow(ni_y,2.0) + Foam::pow(ni_x,2.0))
            * (4*(Foam::pow(ni_z,4.0) + Foam::pow(ni_y,4.0) + Foam::pow(ni_x,4.0))*(ni_z*ni_yz + ni_y*ni_yy + ni_x*ni_xy)
            - 4*(Foam::pow(ni_z,2.0) + Foam::pow(ni_y,2.0) + Foam::pow(ni_x,2.0))
            * (Foam::pow(ni_z,3.0)*ni_yz + Foam::pow(ni_y,3.0)*ni_yy + Foam::pow(ni_x,3.0)*ni_xy))
            - 16*(B*ni*nj + A*(ni_z*nj_z + ni_y*nj_y + ni_x*nj_x))
            * (3*Foam::pow(ni_y,2.0)*ni_yy*Foam::pow((Foam::pow(ni_z,2.0) + Foam::pow(ni_y,2.0) + Foam::pow(ni_x,2.0)),2.0)
            - ni_yy*(Foam::pow(ni_z,2.0) + Foam::pow(ni_y,2.0) + Foam::pow(ni_x,2.0))
            * (Foam::pow(ni_z,4.0) + Foam::pow(ni_y,4.0) + Foam::pow(ni_x,4.0))
            - 4*Foam::pow(ni_y,3.0)*(Foam::pow(ni_z,2.0) + Foam::pow(ni_y,2.0) + Foam::pow(ni_x,2.0))
            * (ni_z*ni_yz + ni_y*ni_yy + ni_x*ni_xy)
            + 6*ni_y*(Foam::pow(ni_z,4.0) + Foam::pow(ni_y,4.0) + Foam::pow(ni_x,4.0))
            * (ni_z*ni_yz + ni_y*ni_yy + ni_x*ni_xy)
            - 4*ni_y*(Foam::pow(ni_z,2.0) + Foam::pow(ni_y,2.0) + Foam::pow(ni_x,2.0))
            * (Foam::pow(ni_z,3.0)*ni_yz + Foam::pow(ni_y,3.0)*ni_yy + Foam::pow(ni_x,3.0)*ni_xy))
            - 16*ni_y*(Foam::pow(ni_z,2.0) + Foam::pow(ni_y,2.0) + Foam::pow(ni_x,2.0))
            * (-Foam::pow(ni_z,4.0) + Foam::pow(ni_z,2.0)*Foam::pow(ni_y,2.0) + Foam::pow(ni_y,2.0)*Foam::pow(ni_x,2.0) - Foam::pow(ni_x,4.0))
            * (B*nj*ni_y + B*ni*nj_y + A*(nj_z*ni_yz + ni_z*nj_yz + nj_y*ni_yy + ni_y*nj_yy + nj_x*ni_xy + ni_x*nj_xy))
            + 4*A*(Foam::pow(ni_z,2.0) + Foam::pow(ni_y,2.0) + Foam::pow(ni_x,2.0))*nj_x
            * (4*(Foam::pow(ni_z,4.0) + Foam::pow(ni_y,4.0) + Foam::pow(ni_x,4.0))*(ni_z*ni_xz + ni_y*ni_xy + ni_x*ni_xx)
            - 4*(Foam::pow(ni_z,2.0) + Foam::pow(ni_y,2.0) + Foam::pow(ni_x,2.0))
            * (Foam::pow(ni_z,3.0)*ni_xz + Foam::pow(ni_y,3.0)*ni_xy + Foam::pow(ni_x,3.0)*ni_xx))
            - 16*(B*ni*nj + A*(ni_z*nj_z + ni_y*nj_y + ni_x*nj_x))
            * (3*Foam::pow(ni_x,2.0)*Foam::pow((Foam::pow(ni_z,2.0) + Foam::pow(ni_y,2.0) + Foam::pow(ni_x,2.0)),2.0)*ni_xx
            - (Foam::pow(ni_z,2.0) + Foam::pow(ni_y,2.0) + Foam::pow(ni_x,2.0))
            * (Foam::pow(ni_z,4.0) + Foam::pow(ni_y,4.0) + Foam::pow(ni_x,4.0))*ni_xx
            - 4*Foam::pow(ni_x,3.0)*(Foam::pow(ni_z,2.0) + Foam::pow(ni_y,2.0) + Foam::pow(ni_x,2.0))
            * (ni_z*ni_xz + ni_y*ni_xy + ni_x*ni_xx)
            + 6*ni_x*(Foam::pow(ni_z,4.0) + Foam::pow(ni_y,4.0) + Foam::pow(ni_x,4.0))
            * (ni_z*ni_xz + ni_y*ni_xy + ni_x*ni_xx)
            - 4*ni_x*(Foam::pow(ni_z,2.0) + Foam::pow(ni_y,2.0) + Foam::pow(ni_x,2.0))
            * (Foam::pow(ni_z,3.0)*ni_xz + Foam::pow(ni_y,3.0)*ni_xy + Foam::pow(ni_x,3.0)*ni_xx))
            + A*Foam::pow((Foam::pow(ni_z,2.0) + Foam::pow(ni_y,2.0) + Foam::pow(ni_x,2.0)),2.0)
            * (3*Foam::pow((Foam::pow(ni_z,2.0) + Foam::pow(ni_y,2.0) + Foam::pow(ni_x,2.0)),2.0)
            - 4*(Foam::pow(ni_z,4.0) + Foam::pow(ni_y,4.0) + Foam::pow(ni_x,4.0)))*nj_xx
            - 16*ni_x*(Foam::pow(ni_z,2.0) + Foam::pow(ni_y,2.0) + Foam::pow(ni_x,2.0))
            * (-Foam::pow(ni_z,4.0) - Foam::pow(ni_y,4.0) + Foam::pow(ni_z,2.0)*Foam::pow(ni_x,2.0) + Foam::pow(ni_y,2.0)*Foam::pow(ni_x,2.0))
            * (B*nj*ni_x + B*ni*nj_x + A*(nj_z*ni_xz + ni_z*nj_xz + nj_y*ni_xy + ni_y*nj_xy + nj_x*ni_xx + ni_x*nj_xx)))
            / (Foam::pow((Foam::pow(ni_z,2.0) + Foam::pow(ni_y,2.0) + Foam::pow(ni_x,2.0)),4.0) + deltaN)
    );

    return Aniso_term;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //