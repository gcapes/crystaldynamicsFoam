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
// SHA1 = 39da1269a5aae6c2138cc97c04608cf394e29b45
//
// unique function name that can be checked if the correct library version
// has been loaded
extern "C" void velocityProfileback_39da1269a5aae6c2138cc97c04608cf394e29b45(bool load)
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
    velocityProfilebackFixedValueFvPatchVectorField
);

} // End namespace Foam


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::
velocityProfilebackFixedValueFvPatchVectorField::
velocityProfilebackFixedValueFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    parent_bctype(p, iF)
{
    if (false)
    {
        printMessage("Construct velocityProfileback : patch/DimensionedField");
    }
}


Foam::
velocityProfilebackFixedValueFvPatchVectorField::
velocityProfilebackFixedValueFvPatchVectorField
(
    const velocityProfilebackFixedValueFvPatchVectorField& rhs,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    parent_bctype(rhs, p, iF, mapper)
{
    if (false)
    {
        printMessage("Construct velocityProfileback : patch/DimensionedField/mapper");
    }
}


Foam::
velocityProfilebackFixedValueFvPatchVectorField::
velocityProfilebackFixedValueFvPatchVectorField
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
        printMessage("Construct velocityProfileback : patch/dictionary");
    }
}


Foam::
velocityProfilebackFixedValueFvPatchVectorField::
velocityProfilebackFixedValueFvPatchVectorField
(
    const velocityProfilebackFixedValueFvPatchVectorField& rhs
)
:
    parent_bctype(rhs),
    dictionaryContent(rhs)
{
    if (false)
    {
        printMessage("Copy construct velocityProfileback");
    }
}


Foam::
velocityProfilebackFixedValueFvPatchVectorField::
velocityProfilebackFixedValueFvPatchVectorField
(
    const velocityProfilebackFixedValueFvPatchVectorField& rhs,
    const DimensionedField<vector, volMesh>& iF
)
:
    parent_bctype(rhs, iF)
{
    if (false)
    {
        printMessage("Construct velocityProfileback : copy/DimensionedField");
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::
velocityProfilebackFixedValueFvPatchVectorField::
~velocityProfilebackFixedValueFvPatchVectorField()
{
    if (false)
    {
        printMessage("Destroy velocityProfileback");
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void
Foam::
velocityProfilebackFixedValueFvPatchVectorField::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    if (false)
    {
        printMessage("updateCoeffs velocityProfileback");
    }

//{{{ begin code
    #line 59 "/home/tomflint/OpenFOAM/custom/Euler_Solid_Mech/tutorials/Fracture_Test_Cases/0/U/boundaryField/Back"
vectorField veloc(patch().size(), vector(0, 0, 0));
scalar radius=0.05e-3;
scalar Peak=10.0;
scalar P =4.0;

            const vectorField& centre(patch().Cf());
            forAll(centre, i)
            {
                veloc[i].z() =((Peak))*Foam::exp(-Foam::pow(1.0*Foam::pow((centre[i].y()-0.1e-3),2.0)/(2.0*radius*radius),P));
                if(centre[i].y()<0.005e-3||centre[i].y()>0.195e-3)
            {
            veloc[i].z() =0.0;
            }
            }
            

operator==(veloc);
//}}} end code

    this->parent_bctype::updateCoeffs();
}


// ************************************************************************* //

