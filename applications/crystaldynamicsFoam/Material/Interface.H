/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | CIPHERFoam - Interface Class
     \\/     M anipulation  |
-------------------------------------------------------------------------------
Description
    Interface class for phase field interactions.
    Handles interface energies, mobilities, and anisotropy calculations.

\*---------------------------------------------------------------------------*/

#ifndef Interface_H
#define Interface_H

// OpenFOAM includes
#include "fvCFD.H"
#include "IOdictionary.H"
#include "volFields.H"
#include "fvMesh.H"
#include "Time.H"
#include "Pair.H"
#include "Hash.H"
#include "HashTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class Interface Declaration
\*---------------------------------------------------------------------------*/

class Interface
:
    public IOdictionary
{
public:

    // =========================================================================
    // Nested class interfacePair
    // =========================================================================
    class interfacePair
    :
        public Pair<word>
    {
    public:

        // Hash functor for interfacePair
        class hash
{
public:

    label operator()(const interfacePair& key) const
    {
        return word::hash()(key.first())
             ^ word::hash()(key.second());
    }
};

        // Constructors
        interfacePair()
        {}

        interfacePair(const word& alpha1Name, const word& alpha2Name)
        :
            Pair<word>(alpha1Name, alpha2Name)
        {}

        // Friend Operators
        friend bool operator==
        (
            const interfacePair& a,
            const interfacePair& b
        )
        {
            return
            (
                ((a.first() == b.first()) && (a.second() == b.second()))
            );
        }

        friend bool operator!=
        (
            const interfacePair& a,
            const interfacePair& b
        )
        {
            return (!(a == b));
        }
    };


private:

    // Private data members
    const fvMesh& mesh_;
    const volScalarField& T_;
    word Interface_type_;


public:

    // Constructors
    Interface(word interfaceType, const volScalarField& T);

    // Destructor
    virtual ~Interface()
    {}

    // Member Functions
    
    //- Assign interface energies from dictionary
    void Assign_Energies(scalar** NRG_arr);
    
    //- Assign interface mobilities from dictionary
    void Assign_Mobilities(scalar** Mob_arr);
    
    //- Assign anisotropy strengths from dictionary
    void Assign_Anisotropy_Strength(scalar** Aniso_strength_arr);
    
    //- Calculate anisotropic energy contribution
    scalar aniso_energy(
        scalar* tensor_i,
        scalar* tensor_j,
        scalar* vec_i,
        scalar* vec_j,
        const scalar& phi_i,
        const scalar& phi_j,
        const scalar& int_width,
        const scalar& deltaN
    );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //