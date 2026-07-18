/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2026 Aalto University
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is originating from OpenFOAM but modified by the authors
    described below.

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

Class
    Foam::Function1s::coolprop

Authors
    Stanislau Stasheuski, Aalto University, 2026.
    stanislau.stasheuski[at]aalto.fi

Description
    Function1 evaluating a CoolProp saturation-curve property, in SI units
    (active unit conversions are refused).

    Outside the saturation curve the value is clamped to the nearest
    endpoint (e.g. Tsat returns Tt below the triple pressure and Tc above
    the critical pressure), and the numerical derivative loses accuracy
    within ~1e-4 (relative) of those endpoints.

    Provides the CoolProp analogue of the saturation Function1s (Antoine,
    ArdenBuck, ...) and of the removed liquidProperties surface tension:
    with \c property \c sigma it slots into the temperatureDependent
    surface tension model, with \c Tsat or \c pSat into phase-change model
    dictionaries (e.g. the cavitation models' \c pSat and \c Tsat entries).

Usage
    \table
        Property | Description                                  | Required
        fluid    | Fluid specification with an optional backend
                   prefix, e.g. "H2O" or "REFPROP::H2O"         | yes
        property | sigma (of T), pSat (of T), Tsat (of p)
                   or hl (of T)                                 | yes
    \endtable

    Example, CoolProp surface tension of water in phaseProperties:
    \verbatim
    sigma
    {
        type            temperatureDependent;
        sigma
        {
            type            coolprop;
            fluid           H2O;
            property        sigma;
        }
    }
    \endverbatim

SourceFiles
    coolprop.C

See also
    Foam::coolpropProperties

\*---------------------------------------------------------------------------*/

#ifndef coolprop_H
#define coolprop_H

#include "coolpropProperties.H"
#include "Function1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace Function1s
{

/*---------------------------------------------------------------------------*\
                          Class coolprop Declaration
\*---------------------------------------------------------------------------*/

class coolprop
:
    public FieldFunction1<scalar, coolprop>
{
public:

    //- The evaluated property
    enum class propertyType
    {
        sigma,
        pSat,
        Tsat,
        hl
    };

    //- Property names, in propertyType enumeration order
    static const NamedEnum<propertyType, 4> propertyNames;


private:

    // Private Data

        //- The evaluated property
        const propertyType property_;

        //- The CoolProp-backed properties evaluator
        autoPtr<coolpropProperties> properties_;


public:

    //- Runtime type information
    TypeName("coolprop");


    // Constructors

        //- Construct from name, units and dictionary
        coolprop
        (
            const word& name,
            const unitSets& units,
            const dictionary& dict
        );

        //- Copy constructor, re-creating the CoolProp states
        coolprop(const coolprop& cp);


    //- Destructor
    virtual ~coolprop();


    // Member Functions

        //- Return the value at a scalar coordinate
        virtual scalar value(const scalar x) const;

        //- Return the derivative at a scalar coordinate
        virtual scalar derivative(const scalar x) const;

        //- Return the integral between two scalar coordinates
        virtual scalar integral(const scalar x1, const scalar x2) const;

        //- Write in dictionary format
        virtual void write(Ostream& os, const unitSets& units) const;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const coolprop&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Function1s
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
