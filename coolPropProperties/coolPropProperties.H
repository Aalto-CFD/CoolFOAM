/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2025-2026 Aalto University
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
    Foam::coolPropProperties

Authors
    Stanislau Stasheuski, Aalto University, 2025-2026.
    stanislau.stasheuski[at]aalto.fi

Description
    The thermophysical properties of a liquid evaluated by the CoolProp
    library (https://coolprop.org) through the low-level interface of its
    shared library, for use with the standard liquid thermodynamics packages:

    \verbatim
    thermoType
    {
        type            heRhoThermo;
        properties      coolprop;
        mixture         pureMixture;
        energy          sensibleInternalEnergy;
    }

    mixture
    {
        // Any CoolProp fluid name, see
        // https://coolprop.org/fluid_properties/PurePseudoPure.html
        H2O;

        // REFPROP::H2O;  // Optional explicit backend prefix
    }
    \endverbatim

    Liquid-phase properties are evaluated from a pressure-temperature flash
    with an imposed liquid phase, which extends smoothly into the metastable
    region beyond saturation in the same spirit as the extrapolated built-in
    liquid correlations; if the flash fails the saturated-liquid state at the
    given temperature is used instead. Vapour properties (Cpg, mug, kappag)
    use the corresponding gas-phase flash, and saturation properties (pv, hl,
    sigma) a quality-temperature flash.

    The absolute enthalpy ha is relative to the CoolProp reference state of
    the fluid and hf is defined as ha at standard conditions, so the sensible
    energies are reference-independent; the absolute energy types are
    therefore not meaningful with this class. The vapour diffusivity D is not
    provided by CoolProp and falls back to the generic API vapour mass
    diffusivity function with the coefficients used by the built-in liquids.

    Limitations, shared with the built-in liquidProperties contract:
      - The equation of state is reported as incompressible (psi = 0,
        CpMCv = 0) although the density returned by the pressure-temperature
        flash does depend on pressure; the pressure dependence is therefore
        not seen by the compressibility terms, which is accurate away from
        the critical region and degrades as the critical point is approached.
      - The saturated-liquid fallback evaluates derivative outputs (alphav,
        Cp) at the Q = 0 endpoint; a backend that refuses first derivatives
        exactly on the saturation line would raise a fatal error there.
      - In the region where the single-phase flash fails, every evaluation
        pays a failed pressure-temperature flash before the saturated
        fallback; this is the cost of the metastable extension and matters
        only for cases sitting deep beyond saturation.

    Reference:
    \verbatim
        Fadiga, E., Casari, N., Suman, A., & Pinelli, M. (2020).
        CoolFOAM: The CoolProp wrapper for OpenFOAM.
        Computer Physics Communications, 250, 107047.
    \endverbatim

SourceFiles
    coolPropProperties.C

\*---------------------------------------------------------------------------*/

#ifndef coolPropProperties_H
#define coolPropProperties_H

#include "liquidProperties.H"
#include "APIdiffusionCoefficient.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class coolPropProperties Declaration
\*---------------------------------------------------------------------------*/

class coolPropProperties
:
    public liquidProperties
{
    // Private Data

        //- Fluid specification, e.g. "H2O" or "REFPROP::H2O"
        const word fluid_;

        //- Handle of the state evaluated at (p, T) with an imposed
        //  liquid phase
        const long liquidState_;

        //- Handle of the state evaluated at (p, T) with an imposed
        //  gas phase
        const long vapourState_;

        //- Handle of the state evaluated on the saturation curve
        const long saturationState_;

        // Cached conditions of the last state updates

            mutable scalar liquidP_, liquidT_;
            mutable scalar vapourP_, vapourT_;
            mutable scalar pvT_, pv_;
            mutable scalar hlT_, hl_;
            mutable scalar sigmaT_, sigma_;

        //- Vapour diffusivity function (not provided by CoolProp)
        Function2s::APIdiffusionCoefficient D_;

        //- Enthalpy at standard conditions [J/kg],
        //  the offset between the absolute and sensible enthalpies
        scalar hf_;


    // Private Member Functions

        //- Update the given memoised state to (p, T) with the imposed phase,
        //  falling back to its saturated state at quality Q and temperature T
        //  limited to [Tt, Tc] if the single-phase flash fails
        long flash
        (
            const long state,
            const char* phase,
            const scalar Q,
            scalar& pMemo,
            scalar& TMemo,
            const scalar p,
            const scalar T
        ) const;

        //- Update and return the liquid state at (p, T),
        //  falling back to the saturated liquid at T
        long liquid(scalar p, scalar T) const;

        //- Update and return the vapour state at (p, T),
        //  falling back to the saturated vapour at T
        long vapour(scalar p, scalar T) const;

        //- Update and return the saturation state at quality Q and
        //  temperature T limited to the [Tt, Tc] interval
        long saturation(scalar Q, scalar T) const;

        //- Construct from the fluid specification and a constructed state
        coolPropProperties(const word& fluid, const long state);


public:

    //- Runtime type information
    TypeName("coolprop");


    // Constructors

        //- Construct from the fluid specification
        coolPropProperties(const word& fluid);

        //- Construct from dictionary
        coolPropProperties(const dictionary& dict);

        //- Copy constructor
        coolPropProperties(const coolPropProperties& cpp);

        //- Construct and return clone
        virtual autoPtr<liquidProperties> clone() const
        {
            return autoPtr<liquidProperties>(new coolPropProperties(*this));
        }


    // Selectors

        //- Return a pointer to a new coolPropProperties created from
        //  the fluid specification
        static autoPtr<coolPropProperties> New(const word& fluid);

        //- Return a pointer to a new coolPropProperties created from
        //  dictionary
        static autoPtr<coolPropProperties> New(const dictionary& dict);


    //- Destructor
    virtual ~coolPropProperties();


    // Member Functions

        //- Return the fluid specification
        const word& fluid() const
        {
            return fluid_;
        }


        // Fundamental equation of state properties

            //- Liquid density [kg/m^3]
            virtual scalar rho(scalar p, scalar T) const;

            //- Volumetric coefficient of thermal expansion [1/T]
            virtual scalar alphav(scalar p, scalar T) const;


        // Fundamental thermodynamic properties

            //- Liquid heat capacity [J/kg/K]
            virtual scalar Cp(scalar p, scalar T) const;

            //- Liquid sensible enthalpy [J/kg]
            virtual scalar hs(scalar p, scalar T) const;

            //- Liquid enthalpy at standard conditions [J/kg]
            //  Note: this is the offset of the CoolProp reference enthalpy
            //  rather than the heat of formation
            virtual scalar hf() const;

            //- Liquid absolute enthalpy [J/kg]
            //  relative to the CoolProp reference state of the fluid
            virtual scalar ha(scalar p, scalar T) const;

            //- Liquid entropy [J/kg/K]
            //  relative to the CoolProp reference state of the fluid
            virtual scalar s(scalar p, scalar T) const;


        // Physical properties

            //- Vapour pressure [Pa]
            virtual scalar pv(scalar p, scalar T) const;

            //- Heat of vapourisation [J/kg]
            virtual scalar hl(scalar p, scalar T) const;

            //- Ideal gas heat capacity [J/kg/K]
            virtual scalar Cpg(scalar p, scalar T) const;

            //- Liquid viscosity [Pa s]
            virtual scalar mu(scalar p, scalar T) const;

            //- Vapour viscosity [Pa s]
            virtual scalar mug(scalar p, scalar T) const;

            //- Liquid thermal conductivity [W/m/K]
            virtual scalar kappa(scalar p, scalar T) const;

            //- Vapour thermal conductivity [W/m/K]
            virtual scalar kappag(scalar p, scalar T) const;

            //- Surface tension [N/m]
            virtual scalar sigma(scalar p, scalar T) const;

            //- Vapour diffusivity [m^2/s]
            virtual scalar D(scalar p, scalar T) const;

            //- Vapour diffusivity [m^2/s] with specified binary pair
            virtual scalar D(scalar p, scalar T, scalar Wb) const;

            //- Invert the vapour pressure relationship to retrieve the
            //  boiling temperature as a function of pressure
            virtual scalar pvInvert(scalar p) const;


    // I-O

        //- Write the fluid specification and the fluid constants
        virtual void write(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
