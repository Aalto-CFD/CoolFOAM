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
    Foam::CoolpropTransport

Authors
    Stanislau Stasheuski, Aalto University, 2026.
    stanislau.stasheuski[at]aalto.fi

Description
    CoolProp-backed transport properties package, evaluating the dynamic
    viscosity and thermal conductivity of the state at (p, T) with the
    imposed phase through the CoolProp low-level interface.

Usage
    Both the \c mu and \c kappa sub-dictionaries take the same entries:

    \table
        Property | Description                                  | Required
        type     | Model type, coolprop                         | no
        fluid    | Fluid specification with an optional backend
                   prefix, e.g. "H2O" or "REFPROP::H2O"         | yes
        phase    | Phase of the evaluated state, liquid or gas  | no (liquid)
    \endtable

    Example specification of CoolProp transport for water:
    \verbatim
    transport
    {
        mu
        {
            type            coolprop;
            fluid           BICUBIC&HEOS::H2O;
        }
        kappa
        {
            type            coolprop;
            fluid           BICUBIC&HEOS::H2O;
        }
    }
    \endverbatim

SourceFiles
    CoolpropTransportI.H
    CoolpropTransport.C

See also
    Foam::coolpropProperties

\*---------------------------------------------------------------------------*/

#ifndef CoolpropTransport_H
#define CoolpropTransport_H

#include "coolpropProperties.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

template<class Thermo>
class CoolpropTransport;

template<class Thermo>
Ostream& operator<<(Ostream&, const CoolpropTransport<Thermo>&);


/*---------------------------------------------------------------------------*\
                      Class CoolpropTransport Declaration
\*---------------------------------------------------------------------------*/

template<class Thermo>
class CoolpropTransport
:
    public Thermo
{
    // Private Data

        //- The CoolProp-backed dynamic viscosity evaluator
        autoPtr<coolpropProperties> muProperties_;

        //- The CoolProp-backed thermal conductivity evaluator
        autoPtr<coolpropProperties> kappaProperties_;


    // Private Constructors

        //- Construct from name and parent dictionary and dictionaries
        CoolpropTransport
        (
            const word& name,
            const dictionary& dict,
            const dictionary& muDict,
            const dictionary& kappaDict
        );


public:

    // Constructors

        //- Construct from name and dictionary
        CoolpropTransport(const word& name, const dictionary& dict);

        //- Copy constructor, re-creating the CoolProp states
        inline CoolpropTransport(const CoolpropTransport&);

        //- Construct as named copy
        inline CoolpropTransport(const word&, const CoolpropTransport&);

        //- Construct and return a clone
        inline autoPtr<CoolpropTransport> clone() const;


    // Member Functions

        //- Return the instantiated type name
        static word typeName()
        {
            return "coolprop<" + Thermo::typeName() + '>';
        }

        //- Dynamic viscosity [kg/m/s]
        inline scalar mu(const scalar p, const scalar T) const;

        //- Thermal conductivity [W/m/K]
        inline scalar kappa(const scalar p, const scalar T) const;

        //- Write to Ostream
        void write(Ostream& os) const;


    // Member Operators

        //- Copy assignment, re-creating the CoolProp states
        inline void operator=(const CoolpropTransport&);


    // Ostream Operator

        friend Ostream& operator<< <Thermo>
        (
            Ostream&,
            const CoolpropTransport&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "CoolpropTransportI.H"

#ifdef NoRepository
    #include "CoolpropTransport.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
