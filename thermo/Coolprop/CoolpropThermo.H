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
    Foam::CoolpropThermo

Authors
    Stanislau Stasheuski, Aalto University, 2026.
    stanislau.stasheuski[at]aalto.fi

Description
    CoolProp-backed thermodynamics package, evaluating the complete
    enthalpy, entropy and heat capacity of the state at (p, T) with the
    imposed phase through the CoolProp low-level interface (rather than
    assembling them from equation-of-state contributions, which are
    therefore zero in the accompanying CoolProp equation of state).

    The absolute enthalpy and entropy are relative to the CoolProp
    reference state of the fluid, optionally shifted by the hf and sf
    offsets (as in the CoolFOAM reference implementation, Fadiga et al.
    2020, and the built-in NSRDS thermodynamics); the "enthalpy of
    formation" is the absolute enthalpy at standard conditions plus hf, so
    the sensible energies are zero-referenced at standard conditions
    independently of the offsets.

Usage
    \table
        Property | Description                                  | Required
        type     | Model type, coolprop                         | no
        fluid    | Fluid specification with an optional backend
                   prefix, e.g. "H2O" or "REFPROP::H2O"         | yes
        phase    | Phase of the evaluated state, liquid or gas  | no (liquid)
        hf       | Offset added to the CoolProp reference
                   absolute enthalpy [J/kg]                     | no (0)
        sf       | Offset added to the CoolProp reference
                   entropy [J/kg/K]                             | no (0)
    \endtable

    Example specification of CoolProp thermodynamics for water:
    \verbatim
    thermodynamics
    {
        type            coolprop;
        fluid           BICUBIC&HEOS::H2O;
        hf              0;
        sf              0;
    }
    \endverbatim

SourceFiles
    CoolpropThermoI.H
    CoolpropThermo.C

See also
    Foam::coolpropProperties

\*---------------------------------------------------------------------------*/

#ifndef CoolpropThermo_H
#define CoolpropThermo_H

#include "coolpropProperties.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

template<class EquationOfState>
class CoolpropThermo;

template<class EquationOfState>
Ostream& operator<<(Ostream&, const CoolpropThermo<EquationOfState>&);


/*---------------------------------------------------------------------------*\
                       Class CoolpropThermo Declaration
\*---------------------------------------------------------------------------*/

template<class EquationOfState>
class CoolpropThermo
:
    public EquationOfState
{
    // Private Data

        //- The CoolProp-backed properties evaluator
        autoPtr<coolpropProperties> properties_;

        //- Offset added to the CoolProp reference absolute enthalpy [J/kg]
        scalar hf_;

        //- Offset added to the CoolProp reference entropy [J/kg/K]
        scalar sf_;


    // Private Constructors

        //- Construct from name and parent dictionary and dictionary
        CoolpropThermo
        (
            const word& name,
            const dictionary& dict,
            const dictionary& subDict
        );


public:

    // Constructors

        //- Construct from name and dictionary
        CoolpropThermo(const word& name, const dictionary& dict);

        //- Copy constructor, re-creating the CoolProp states
        inline CoolpropThermo(const CoolpropThermo&);

        //- Construct as named copy
        inline CoolpropThermo(const word&, const CoolpropThermo&);

        //- Construct and return a clone
        inline autoPtr<CoolpropThermo> clone() const;


    // Member Functions

        //- Return the instantiated type name
        static word typeName()
        {
            return "coolprop<" + EquationOfState::typeName() + '>';
        }

        //- Limit the temperature to the range where the backend evaluates
        //  a genuine state, so the energy -> temperature inversions cannot
        //  wander onto the clamped-saturation plateau and stall
        inline scalar limit(const scalar T) const;


        // Fundamental properties

            //- Heat capacity at constant pressure [J/kg/K]
            inline scalar Cp(const scalar p, const scalar T) const;

            //- Absolute enthalpy [J/kg]
            inline scalar ha(const scalar p, const scalar T) const;

            //- Sensible enthalpy [J/kg]
            inline scalar hs(const scalar p, const scalar T) const;

            //- Enthalpy of formation [J/kg]
            inline scalar hf() const;

            //- Entropy [J/kg/K]
            inline scalar s(const scalar p, const scalar T) const;

            //- Gibbs free energy of the mixture in the standard state [J/kg]
            inline scalar gStd(const scalar T) const;

            //- Temperature derivative of heat capacity at constant pressure
            inline scalar dCpdT(const scalar p, const scalar T) const;

            //- Heat capacity at constant volume [J/kg/K]
            inline scalar Cv(const scalar p, const scalar T) const;

            //- Sensible internal energy [J/kg]
            inline scalar es(const scalar p, const scalar T) const;

            //- Absolute internal energy [J/kg]
            inline scalar ea(const scalar p, const scalar T) const;


        // I-O

            //- Write to Ostream
            void write(Ostream& os) const;


    // Member Operators

        //- Copy assignment, re-creating the CoolProp states
        inline void operator=(const CoolpropThermo&);


    // Ostream Operator

        friend Ostream& operator<< <EquationOfState>
        (
            Ostream&,
            const CoolpropThermo&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "CoolpropThermoI.H"

#ifdef NoRepository
    #include "CoolpropThermo.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
