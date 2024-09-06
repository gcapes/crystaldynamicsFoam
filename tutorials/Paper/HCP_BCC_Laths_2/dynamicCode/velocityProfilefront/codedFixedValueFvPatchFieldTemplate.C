/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) YEAR OpenFOAM Foundation
     \\/     M anipulation  |
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

#include "codedFixedValueFvPatchFieldTemplate.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "unitConversion.H"
//{{{ begin codeInclude

//}}} end codeInclude


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

//{{{ begin localCode

//}}} end localCode


// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

extern "C"
{
    // dynamicCode:
    // SHA1 = 9b422f8cadbc73c636689cd219d93c7810e6102b
    //
    // unique function name that can be checked if the correct library version
    // has been loaded
    void velocityProfilefront_9b422f8cadbc73c636689cd219d93c7810e6102b(bool load)
    {
        if (load)
        {
            // code that can be explicitly executed after loading
        }
        else
        {
            // code that can be explicitly executed before unloading
        }
    }
}

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

makeRemovablePatchTypeField
(
    fvPatchVectorField,
    velocityProfilefrontFixedValueFvPatchVectorField
);


const char* const velocityProfilefrontFixedValueFvPatchVectorField::SHA1sum =
    "9b422f8cadbc73c636689cd219d93c7810e6102b";


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

velocityProfilefrontFixedValueFvPatchVectorField::
velocityProfilefrontFixedValueFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchField<vector>(p, iF)
{
    if (false)
    {
        Info<<"construct velocityProfilefront sha1: 9b422f8cadbc73c636689cd219d93c7810e6102b"
            " from patch/DimensionedField\n";
    }
}


velocityProfilefrontFixedValueFvPatchVectorField::
velocityProfilefrontFixedValueFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchField<vector>(p, iF, dict)
{
    if (false)
    {
        Info<<"construct velocityProfilefront sha1: 9b422f8cadbc73c636689cd219d93c7810e6102b"
            " from patch/dictionary\n";
    }
}


velocityProfilefrontFixedValueFvPatchVectorField::
velocityProfilefrontFixedValueFvPatchVectorField
(
    const velocityProfilefrontFixedValueFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchField<vector>(ptf, p, iF, mapper)
{
    if (false)
    {
        Info<<"construct velocityProfilefront sha1: 9b422f8cadbc73c636689cd219d93c7810e6102b"
            " from patch/DimensionedField/mapper\n";
    }
}


velocityProfilefrontFixedValueFvPatchVectorField::
velocityProfilefrontFixedValueFvPatchVectorField
(
    const velocityProfilefrontFixedValueFvPatchVectorField& ptf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchField<vector>(ptf, iF)
{
    if (false)
    {
        Info<<"construct velocityProfilefront sha1: 9b422f8cadbc73c636689cd219d93c7810e6102b "
            "as copy/DimensionedField\n";
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

velocityProfilefrontFixedValueFvPatchVectorField::
~velocityProfilefrontFixedValueFvPatchVectorField()
{
    if (false)
    {
        Info<<"destroy velocityProfilefront sha1: 9b422f8cadbc73c636689cd219d93c7810e6102b\n";
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void velocityProfilefrontFixedValueFvPatchVectorField::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    if (false)
    {
        Info<<"updateCoeffs velocityProfilefront sha1: 9b422f8cadbc73c636689cd219d93c7810e6102b\n";
    }

//{{{ begin code
    #line 33 "/home/tomflint/OpenFOAM/custom/Euler_Solid_Mech/tutorials/Paper/HCP_BCC_Laths_2/0/U/boundaryField/Front"
vectorField veloc(patch().size(), vector(0, 0, 0));
scalar radius=0.025e-3;
scalar Peak=-10.0;

            const vectorField& centre(patch().Cf());
            forAll(centre, i)
            {
                veloc[i].z() =((Peak))*Foam::exp(-1.0*Foam::pow((centre[i].y()-0.1e-3),2.0)/(2.0*radius*radius));
            }

operator==(veloc);
//}}} end code

    this->fixedValueFvPatchField<vector>::updateCoeffs();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //

