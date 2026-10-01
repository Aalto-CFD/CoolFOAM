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
    Foam::coolpropProperties

Authors
    Stanislau Stasheuski, Aalto University, 2025-2026.
    stanislau.stasheuski[at]aalto.fi

Description
    Thermophysical properties evaluator backing the CoolProp thermodynamics
    stack (Foam::Coolprop, Foam::CoolpropThermo and Foam::CoolpropTransport)
    through the CoolProp low-level interface. Selected by

    \verbatim
    thermoType
    {
        type            heRhoThermo;
        mixture         pureMixture;
        transport       coolprop;
        thermo          coolprop;
        equationOfState coolprop;
        specie          specie;
        energy          sensibleInternalEnergy;
    }

    mixture
    {
        // Any CoolProp fluid name with optional explicit backend prefix
        // https://coolprop.org/fluid_properties/PurePseudoPure.html
        equationOfState
        {
            type            coolprop;
            fluid           BICUBIC&HEOS::H2O;

            // Optional; "liquid" (default) or "gas": the imposed phase of
            // the flash. Needed for fluids at conditions where the liquid
            // role is meaningless, e.g. a permanent gas such as Air at
            // typical ambient temperature (T >> Tc)
            phase           gas;
        }
        thermodynamics  { $equationOfState; }
        transport
        {
            mu              { $equationOfState; }
            kappa           { $equationOfState; }
        }
    }
    \endverbatim

    It also remains a Foam::liquidProperties, usable wherever those are
    (e.g. Lagrangian liquid thermodynamics).

    Reference:
    \verbatim
        Fadiga, E., Casari, N., Suman, A., & Pinelli, M. (2020).
        CoolFOAM: The CoolProp wrapper for OpenFOAM.
        Computer Physics Communications, 250, 107047.
    \endverbatim

SourceFiles
    coolpropProperties.C

\*---------------------------------------------------------------------------*/

#ifndef coolpropProperties_H
#define coolpropProperties_H

#include "liquidProperties.H"
#include "APIdiffusionCoefficient.H"
#include "NamedEnum.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class coolpropProperties Declaration
\*---------------------------------------------------------------------------*/

class coolpropProperties
:
    public liquidProperties
{
public:

    //- Which phase this instance evaluates the primary properties as;
    //  both roles report the real CoolProp compressibility/energy
    enum class phaseType
    {
        liquid,
        gas
    };

    //- Phase type names, in phaseType enumeration order
    static const NamedEnum<phaseType, 2> phaseTypeNames;


private:

    // Private Data

        //- Fluid specification, e.g. "H2O" or "REFPROP::H2O"
        const word fluid_;

        //- Is the backend an SVDSBTL table, e.g. "SVDSBTL&REFPROP::H2O",
        //  which neither imposes the phase nor provides the heat capacities,
        //  derivatives, surface tension, acentric factor or saturated
        //  endpoints of a flash
        const bool tabulated_;

        //- Handle of the state evaluated at (p, T) with an imposed
        //  liquid phase
        const long liquidState_;

        //- Handle of the state evaluated at (p, T) with an imposed
        //  gas phase
        const long vapourState_;

        //- Handle of the state evaluated on the saturation curve
        const long saturationState_;

        //- Handles of the states of the source backend of a tabulated
        //  backend, e.g. "REFPROP::H2O", for the properties the table does
        //  not provide (-1 otherwise): the derivatives of the primary
        //  phase, and the other (surface tension, ideal-gas heat capacity)
        const long derivState_;
        const long sourceState_;

        // Cached conditions of the last state updates

            mutable scalar liquidP_, liquidT_;
            mutable scalar vapourP_, vapourT_;
            mutable scalar pvT_, pv_;
            mutable scalar hlT_, hl_;
            mutable scalar sigmaT_, sigma_;
            mutable scalar derivP_, derivT_;

        // Cached derivatives of a tabulated backend

            mutable scalar Cp_, alphav_, psi_, CpMCv_;

        //- Vapour diffusivity function (not provided by CoolProp)
        Function2s::APIdiffusionCoefficient D_;

        //- Which phase this instance evaluates the primary
        //  (rho, alphav, Cp, ha, s, mu, kappa, psi, CpMCv, ea) properties as
        const phaseType phase_;

        //- Enthalpy at standard conditions [J/kg],
        //  the offset between the absolute and sensible enthalpies
        scalar hf_;

        //- Internal energy at standard conditions [J/kg],
        //  the offset between the absolute and sensible internal energies
        scalar ef_;


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

        //- Update and return the given state on the saturation curve at
        //  quality Q and temperature T limited to the [Tt, Tc] interval
        long saturation(const long state, scalar Q, scalar T) const;

        //- Update and return the saturation state at quality Q and
        //  temperature T limited to the [Tt, Tc] interval
        long saturation(scalar Q, scalar T) const;

        //- Is the updated state of a tabulated backend, which does not
        //  impose the phase, within the table and on the side of the
        //  saturation curve of the phase with the saturated quality Q at
        //  temperature T
        bool imposedSide(const long state, const scalar Q, const scalar T)
            const;

        //- Evaluate and cache the Cp, alphav, psi and CpMCv of a tabulated
        //  backend at (p, T), which the table does not provide, from its
        //  source backend at the density and temperature of the resolved
        //  table state with the phase imposed
        void derivatives(const scalar p, const scalar T) const;

        //- Update and return the state for the imposed primary phase_
        //  (liquid() or vapour()) at (p, T)
        long primaryState(scalar p, scalar T) const;

        //- Return the pressure difference between p and the given primary
        //  state at temperature T of a tabulated backend, non-zero for its
        //  fallback next to the saturation curve, from which rho, ha, ea
        //  and s are extrapolated linearly in p (a first-order metastable
        //  extension consistent with psi, which the table does not
        //  provide), limited to a small relative density change, except
        //  for the density of the gas (not energy); zero otherwise
        scalar metastableDp
        (
            const long state,
            const scalar p,
            const scalar T,
            const bool energy
        ) const;

        //- Construct from the fluid specification, imposed phase and a
        //  constructed state
        coolpropProperties
        (
            const word& fluid,
            const long state,
            const phaseType phase
        );


public:

    //- Runtime type information
    TypeName("coolprop");


    // Constructors

        //- Construct from the fluid specification, optionally imposing
        //  the gas phase (defaults to the liquid phase)
        coolpropProperties
        (
            const word& fluid,
            const phaseType phase = phaseType::liquid
        );

        //- Construct from dictionary
        coolpropProperties(const dictionary& dict);

        //- Copy constructor
        coolpropProperties(const coolpropProperties& cpp);

        //- Construct and return clone
        virtual autoPtr<liquidProperties> clone() const
        {
            return autoPtr<liquidProperties>(new coolpropProperties(*this));
        }


    // Selectors

        //- Return a pointer to a new coolpropProperties created from
        //  the fluid specification, optionally imposing the gas phase
        static autoPtr<coolpropProperties> New
        (
            const word& fluid,
            const phaseType phase = phaseType::liquid
        );

        //- Return a pointer to a new coolpropProperties created from
        //  dictionary
        static autoPtr<coolpropProperties> New(const dictionary& dict);


    //- Destructor
    virtual ~coolpropProperties();


    // Member Functions

        //- Return the fluid specification
        const word& fluid() const
        {
            return fluid_;
        }

        //- Return which phase the primary properties evaluate as
        phaseType phase() const
        {
            return phase_;
        }


        // Fundamental equation of state properties

            //- Density of the primary phase [kg/m^3]
            virtual scalar rho(scalar p, scalar T) const;

            //- Volumetric coefficient of thermal expansion [1/T]
            virtual scalar alphav(scalar p, scalar T) const;

            //- Compressibility, drho/dp at constant T [s^2/m^2]
            //  from the real CoolProp isothermal compressibility.
            //  Hides (does not override) liquidProperties::psi, which is
            //  non-virtual upstream; see the ea() note on why hiding is
            //  safe here
            scalar psi(scalar p, scalar T) const;

            //- Cp - Cv [J/kg/K]
            //  Cp minus the real CoolProp Cv.
            //  Hides (does not override) liquidProperties::CpMCv; see the
            //  ea() note
            scalar CpMCv(scalar p, scalar T) const;


        // Fundamental thermodynamic properties

            //- Heat capacity of the primary phase [J/kg/K]
            virtual scalar Cp(scalar p, scalar T) const;

            //- Sensible enthalpy of the primary phase [J/kg]
            virtual scalar hs(scalar p, scalar T) const;

            //- Enthalpy at standard conditions [J/kg]
            //  Note: this is the offset of the CoolProp reference enthalpy
            //  rather than the heat of formation
            virtual scalar hf() const;

            //- Absolute enthalpy of the primary phase [J/kg]
            //  relative to the CoolProp reference state of the fluid
            virtual scalar ha(scalar p, scalar T) const;

            //- Absolute internal energy of the primary phase [J/kg]
            //  relative to the CoolProp reference state of the fluid.
            //  Unlike the liquidProperties contract (ea = ha), the real
            //  CoolProp internal energy, as are psi and CpMCv above.
            //  Hides (does not override) liquidProperties::ea, which is
            //  fine here: this class is only ever accessed concretely,
            //  through the layers of the Coolprop thermodynamics stack,
            //  never a liquidProperties base pointer, so the static (not
            //  dynamic) type always resolves to this definition
            scalar ea(scalar p, scalar T) const;

            //- Sensible internal energy of the primary phase [J/kg]
            //  See ea() for the non-virtual hiding note
            scalar es(scalar p, scalar T) const;

            //- Entropy of the primary phase [J/kg/K]
            //  relative to the CoolProp reference state of the fluid
            virtual scalar s(scalar p, scalar T) const;


        // Physical properties

            //- Vapour pressure [Pa]
            virtual scalar pv(scalar p, scalar T) const;

            //- Heat of vapourisation [J/kg]
            virtual scalar hl(scalar p, scalar T) const;

            //- Ideal gas heat capacity [J/kg/K]
            virtual scalar Cpg(scalar p, scalar T) const;

            //- Viscosity of the primary phase [Pa s]
            virtual scalar mu(scalar p, scalar T) const;

            //- Vapour viscosity [Pa s]
            virtual scalar mug(scalar p, scalar T) const;

            //- Thermal conductivity of the primary phase [W/m/K]
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
            //  boiling temperature as a function of pressure;
            //  -1 below the triple pressure, as the base implementation
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
