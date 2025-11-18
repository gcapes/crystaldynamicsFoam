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
// SHA1 = 5d6d7c8ef316988a815e7b7e20bc21ec8713e372
//
// unique function name that can be checked if the correct library version
// has been loaded
extern "C" void velocityProfilefront_5d6d7c8ef316988a815e7b7e20bc21ec8713e372(bool load)
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
    velocityProfilefrontFixedValueFvPatchVectorField
);

} // End namespace Foam


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::
velocityProfilefrontFixedValueFvPatchVectorField::
velocityProfilefrontFixedValueFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    parent_bctype(p, iF)
{
    if (false)
    {
        printMessage("Construct velocityProfilefront : patch/DimensionedField");
    }
}


Foam::
velocityProfilefrontFixedValueFvPatchVectorField::
velocityProfilefrontFixedValueFvPatchVectorField
(
    const velocityProfilefrontFixedValueFvPatchVectorField& rhs,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    parent_bctype(rhs, p, iF, mapper)
{
    if (false)
    {
        printMessage("Construct velocityProfilefront : patch/DimensionedField/mapper");
    }
}


Foam::
velocityProfilefrontFixedValueFvPatchVectorField::
velocityProfilefrontFixedValueFvPatchVectorField
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
        printMessage("Construct velocityProfilefront : patch/dictionary");
    }
}


Foam::
velocityProfilefrontFixedValueFvPatchVectorField::
velocityProfilefrontFixedValueFvPatchVectorField
(
    const velocityProfilefrontFixedValueFvPatchVectorField& rhs
)
:
    parent_bctype(rhs),
    dictionaryContent(rhs)
{
    if (false)
    {
        printMessage("Copy construct velocityProfilefront");
    }
}


Foam::
velocityProfilefrontFixedValueFvPatchVectorField::
velocityProfilefrontFixedValueFvPatchVectorField
(
    const velocityProfilefrontFixedValueFvPatchVectorField& rhs,
    const DimensionedField<vector, volMesh>& iF
)
:
    parent_bctype(rhs, iF)
{
    if (false)
    {
        printMessage("Construct velocityProfilefront : copy/DimensionedField");
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::
velocityProfilefrontFixedValueFvPatchVectorField::
~velocityProfilefrontFixedValueFvPatchVectorField()
{
    if (false)
    {
        printMessage("Destroy velocityProfilefront");
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void
Foam::
velocityProfilefrontFixedValueFvPatchVectorField::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    if (false)
    {
        printMessage("updateCoeffs velocityProfilefront");
    }

//{{{ begin code
    #line 32 "/home/tomflint/OpenFOAM/custom/Euler_Solid_Mech/tutorials/MPMFoam/Case4/0/gridVelocity_/boundaryField/Front"
vectorField veloc(patch().size(), vector(0, 0, 0));
scalar radius=0.05e-3;
scalar Peak=-10.0;
scalar P =5.0;

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

