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
    Foam::Coolprop

Authors
    Stanislau Stasheuski, Aalto University, 2026.
    stanislau.stasheuski[at]aalto.fi

Description
    CoolProp-backed equation of state, evaluating the state at (p, T) with
    the imposed phase through the CoolProp low-level interface.

    The specie sub-dictionary is honoured when present (as in the CoolFOAM
    reference implementation, Fadiga et al. 2020) but is not required: by
    default the molecular weight is obtained from CoolProp.

    The h, e, Cp, Cv, sp and sv "contributions" of the layered
    equation-of-state interface are zero: the CoolPropThermo layer this
    equation of state is intended to be stacked with evaluates the complete
    thermodynamic properties through CoolProp rather than assembling them
    from equation-of-state contributions.

Usage
    \table
        Property | Description                                  | Required
        type     | Model type, coolprop                         | no
        fluid    | Fluid specification with an optional backend
                   prefix, e.g. "H2O" or "REFPROP::H2O"         | yes
        phase    | Phase of the evaluated state, liquid or gas  | no (liquid)
    \endtable

    Example specification of the CoolProp equation of state for water:
    \verbatim
    equationOfState
    {
        type            coolprop;
        fluid           BICUBIC&HEOS::H2O;
    }
    \endverbatim

SourceFiles
    CoolpropI.H
    Coolprop.C

See also
    Foam::coolpropProperties

\*---------------------------------------------------------------------------*/

#ifndef Coolprop_H
#define Coolprop_H

#include "coolpropProperties.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

template<class Specie>
class Coolprop;

template<class Specie>
Ostream& operator<<(Ostream&, const Coolprop<Specie>&);


/*---------------------------------------------------------------------------*\
                          Class Coolprop Declaration
\*---------------------------------------------------------------------------*/

template<class Specie>
class Coolprop
:
    public Specie
{
    // Private Data

        //- The CoolProp-backed properties evaluator
        autoPtr<coolpropProperties> properties_;


    // Private Constructors

        //- Construct from a name, a parent dictionary and constructed
        //  properties
        inline Coolprop
        (
            const word& name,
            const dictionary& dict,
            autoPtr<coolpropProperties> properties
        );


public:

    // Constructors

        //- Construct from name and dictionary
        Coolprop(const word& name, const dictionary& dict);

        //- Copy constructor, re-creating the CoolProp states
        inline Coolprop(const Coolprop&);

        //- Construct as named copy
        inline Coolprop(const word& name, const Coolprop&);

        //- Construct and return a clone
        inline autoPtr<Coolprop> clone() const;


    // Member Functions

        //- Return the instantiated type name
        static word typeName()
        {
            return "coolprop<" + word(Specie::typeName_()) + '>';
        }

        //- Return the CoolProp-backed properties evaluator
        inline const coolpropProperties& properties() const;


        // Fundamental properties

            //- Is the equation of state is incompressible i.e. rho != f(p)
            static const bool incompressible = false;

            //- Is the equation of state is isochoric i.e. rho = const
            static const bool isochoric = false;

            //- Return density [kg/m^3]
            inline scalar rho(const scalar p, const scalar T) const;

            //- Return enthalpy contribution [J/kg]
            inline scalar h(const scalar p, const scalar T) const;

            //- Return Cp contribution [J/(kg K]
            inline scalar Cp(const scalar p, const scalar T) const;

            //- Return internal energy contribution [J/kg]
            inline scalar e(const scalar p, const scalar T) const;

            //- Return Cv contribution [J/(kg K]
            inline scalar Cv(const scalar p, const scalar T) const;

            //- Return entropy contribution to the integral of Cp/T [J/kg/K]
            inline scalar sp(const scalar p, const scalar T) const;

            //- Return entropy contribution to the integral of Cv/T [J/kg/K]
            inline scalar sv(const scalar p, const scalar T) const;

            //- Return compressibility [s^2/m^2]
            inline scalar psi(const scalar p, const scalar T) const;

            //- Return compression factor []
            inline scalar Z(const scalar p, const scalar T) const;

            //- Return (Cp - Cv) [J/(kg K]
            inline scalar CpMCv(const scalar p, const scalar T) const;

            //- Return volumetric coefficient of thermal expansion [1/T]
            inline scalar alphav(const scalar p, const scalar T) const;


        // I-O

            //- Write to Ostream
            void write(Ostream& os) const;


    // Member Operators

        //- Copy assignment, re-creating the CoolProp states
        inline void operator=(const Coolprop&);


    // Ostream Operator

        friend Ostream& operator<< <Specie>
        (
            Ostream&,
            const Coolprop&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "CoolpropI.H"

#ifdef NoRepository
    #include "Coolprop.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
