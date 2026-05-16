/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019-2021 OpenCFD Ltd.
    Copyright (C) YEAR AUTHOR, AFFILIATION
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

\*---------------------------------------------------------------------------*/

#include "fixedValueFvPatchFieldTemplate.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "unitConversion.H"
#include "PatchFunction1.H"

//{{{ begin codeInclude

//}}} end codeInclude


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

//{{{ begin localCode

//}}} end localCode


// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

// dynamicCode:
// SHA1 = fb2363259c38fa83a1026fb454f75a72fd2635a4
//
// unique function name that can be checked if the correct library version
// has been loaded
extern "C" void planeStrainFront_fb2363259c38fa83a1026fb454f75a72fd2635a4(bool load)
{
    if (load)
    {
        // Code that can be explicitly executed after loading
    }
    else
    {
        // Code that can be explicitly executed before unloading
    }
}

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

makeRemovablePatchTypeField
(
    fvPatchVectorField,
    planeStrainFrontFixedValueFvPatchVectorField
);

} // End namespace Foam


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::
planeStrainFrontFixedValueFvPatchVectorField::
planeStrainFrontFixedValueFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    parent_bctype(p, iF)
{
    if (false)
    {
        printMessage("Construct planeStrainFront : patch/DimensionedField");
    }
}


Foam::
planeStrainFrontFixedValueFvPatchVectorField::
planeStrainFrontFixedValueFvPatchVectorField
(
    const planeStrainFrontFixedValueFvPatchVectorField& rhs,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    parent_bctype(rhs, p, iF, mapper)
{
    if (false)
    {
        printMessage("Construct planeStrainFront : patch/DimensionedField/mapper");
    }
}


Foam::
planeStrainFrontFixedValueFvPatchVectorField::
planeStrainFrontFixedValueFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    parent_bctype(p, iF, dict)
{
    if (false)
    {
        printMessage("Construct planeStrainFront : patch/dictionary");
    }
}


Foam::
planeStrainFrontFixedValueFvPatchVectorField::
planeStrainFrontFixedValueFvPatchVectorField
(
    const planeStrainFrontFixedValueFvPatchVectorField& rhs
)
:
    parent_bctype(rhs),
    dictionaryContent(rhs)
{
    if (false)
    {
        printMessage("Copy construct planeStrainFront");
    }
}


Foam::
planeStrainFrontFixedValueFvPatchVectorField::
planeStrainFrontFixedValueFvPatchVectorField
(
    const planeStrainFrontFixedValueFvPatchVectorField& rhs,
    const DimensionedField<vector, volMesh>& iF
)
:
    parent_bctype(rhs, iF)
{
    if (false)
    {
        printMessage("Construct planeStrainFront : copy/DimensionedField");
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::
planeStrainFrontFixedValueFvPatchVectorField::
~planeStrainFrontFixedValueFvPatchVectorField()
{
    if (false)
    {
        printMessage("Destroy planeStrainFront");
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void
Foam::
planeStrainFrontFixedValueFvPatchVectorField::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    if (false)
    {
        printMessage("updateCoeffs planeStrainFront");
    }

//{{{ begin code
    #line 43 "/home/tomflint/OpenFOAM/custom/Euler_Solid_Mech/tutorials/Validation/Raabe_Validation/0/U/boundaryField/Front"
const scalar edot  = 250000.0;
            const scalar t_ramp = 5e-9;
            const scalar t = this->db().time().value();
            const scalar ramp = Foam::min(t / t_ramp, 1.0);

            vectorField veloc(patch().size(), vector(0, 0, 0));
            const vectorField& cf = patch().Cf();
            forAll(cf, i)
            {
                veloc[i].x() = 0.0;
                veloc[i].y() = ramp * edot * cf[i].y();
                veloc[i].z() = 0.0;  // z = 0 on this face
            }
            operator==(veloc);
//}}} end code

    this->parent_bctype::updateCoeffs();
}


// ************************************************************************* //

