/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2025-2026 Aalto University
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is originating from OpenFOAM but modified by authors described
    in the according header file.

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

#include "coolpropProperties.H"
#include "addToRunTimeSelectionTable.H"

#include "thermodynamicConstants.H"
using namespace Foam::constant::thermodynamic;

#include "CoolPropLib.h"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(coolpropProperties, 0);
    addToRunTimeSelectionTable(liquidProperties, coolpropProperties, dictionary);

    const NamedEnum<coolpropProperties::phaseType, 2>
    coolpropProperties::phaseTypeNames
    {
        "liquid",
        "gas"
    };
}

// * * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * //

namespace Foam
{

//- Length of the message buffers of the CoolProp C interface
static const long coolpropMsgLen = 1000;


//- Input-pair and parameter indices of the CoolProp C interface,
//  resolved once on first use
struct coolpropIndices
{
    const long PT, QT, PQ;
    const long T, p, rho, Cp, Cp0, Cv, u, h, s, mu, kappa, sigma, alphav,
        kappaT;
    const long Tc, pc, rhoc, Tt, pt, W, omega;

    coolpropIndices()
    :
        PT(get_input_pair_index("PT_INPUTS")),
        QT(get_input_pair_index("QT_INPUTS")),
        PQ(get_input_pair_index("PQ_INPUTS")),
        T(get_param_index("T")),
        p(get_param_index("P")),
        rho(get_param_index("Dmass")),
        Cp(get_param_index("Cpmass")),
        Cp0(get_param_index("Cp0mass")),
        Cv(get_param_index("Cvmass")),
        u(get_param_index("Umass")),
        h(get_param_index("Hmass")),
        s(get_param_index("Smass")),
        mu(get_param_index("viscosity")),
        kappa(get_param_index("conductivity")),
        sigma(get_param_index("surface_tension")),
        alphav(get_param_index("isobaric_expansion_coefficient")),
        kappaT(get_param_index("isothermal_compressibility")),
        Tc(get_param_index("T_critical")),
        pc(get_param_index("p_critical")),
        rhoc(get_param_index("rhomolar_critical")),
        Tt(get_param_index("T_triple")),
        pt(get_param_index("p_triple")),
        W(get_param_index("molar_mass")),
        omega(get_param_index("acentric"))
    {}
};


static const coolpropIndices& coolprop()
{
    static const coolpropIndices indices;
    return indices;
}


//- Construct a new state for the given fluid specification, splitting an
//  optional "Backend::" prefix, e.g. "REFPROP::H2O"; defaults to "HEOS"
static long newCoolpropState(const word& fluid)
{
    const std::string::size_type sep = fluid.find("::");

    const std::string backend =
        sep == std::string::npos ? std::string("HEOS") : fluid.substr(0, sep);

    const std::string name =
        sep == std::string::npos
      ? static_cast<const std::string&>(fluid)
      : fluid.substr(sep + 2);

    long errcode = 0;
    char msg[coolpropMsgLen];

    const long state = AbstractState_factory
    (
        backend.c_str(),
        name.c_str(),
        &errcode,
        msg,
        coolpropMsgLen
    );

    if (errcode)
    {
        static const int fluidsLen = 16384;
        char fluids[fluidsLen] = "<unavailable>";
        get_global_param_string("FluidsList", fluids, fluidsLen);

        FatalErrorInFunction
            << "Cannot construct CoolProp fluid " << fluid
            << ": " << msg << nl << nl
            << "Valid CoolProp fluid names are:" << nl << fluids
            << exit(FatalError);
    }

    return state;
}


//- Free the state with the given handle
static void freeCoolpropState(const long state)
{
    long errcode = 0;
    char msg[coolpropMsgLen];

    AbstractState_free(state, &errcode, msg, coolpropMsgLen);
}


//- Update the given state, returning success and the error message
static bool coolpropUpdate
(
    const long state,
    const long inputPair,
    const scalar value1,
    const scalar value2,
    char* msg
)
{
    long errcode = 0;

    AbstractState_update
    (
        state,
        inputPair,
        value1,
        value2,
        &errcode,
        msg,
        coolpropMsgLen
    );

    return errcode == 0;
}


//- Impose the given phase, or lift the imposition for nullptr;
//  ignored for backends without phase imposition support
static void coolpropPhase(const long state, const char* phase)
{
    long errcode = 0;
    char msg[coolpropMsgLen];

    if (phase)
    {
        AbstractState_specify_phase
        (
            state,
            phase,
            &errcode,
            msg,
            coolpropMsgLen
        );
    }
    else
    {
        AbstractState_unspecify_phase(state, &errcode, msg, coolpropMsgLen);
    }
}


//- Return the keyed output of the given state of the given fluid
static scalar coolpropOutput
(
    const long state,
    const long key,
    const word& fluid
)
{
    long errcode = 0;
    char msg[coolpropMsgLen];

    const scalar value =
        AbstractState_keyed_output(state, key, &errcode, msg, coolpropMsgLen);

    if (errcode)
    {
        FatalErrorInFunction
            << "CoolProp error for fluid " << fluid << ": " << msg
            << exit(FatalError);
    }

    return value;
}


//- Return the keyed output of the given state,
//  or the given default if it is not available for the fluid
static scalar coolpropOutput
(
    const long state,
    const long key,
    const scalar deflt
)
{
    long errcode = 0;
    char msg[coolpropMsgLen];

    const scalar value =
        AbstractState_keyed_output(state, key, &errcode, msg, coolpropMsgLen);

    return errcode == 0 ? value : deflt;
}


//- Return the keyed output of the saturated phase ("liquid" or "gas") of the
//  given two-phase state of the given fluid, avoiding a second flash
static scalar coolpropOutputSat
(
    const long state,
    const char* saturatedPhase,
    const long key,
    const word& fluid
)
{
    long errcode = 0;
    char msg[coolpropMsgLen];

    const scalar value = AbstractState_keyed_output_satState
    (
        state,
        saturatedPhase,
        key,
        &errcode,
        msg,
        coolpropMsgLen
    );

    if (errcode)
    {
        FatalErrorInFunction
            << "CoolProp error for fluid " << fluid << ": " << msg
            << exit(FatalError);
    }

    return value;
}


//- Normal boiling temperature, or the triple-point temperature for fluids
//  which do not boil at atmospheric pressure (e.g. CO2)
static scalar coolpropTb(const long state, const word& fluid)
{
    char msg[coolpropMsgLen];

    // Normal boiling point is defined at one standard atmosphere, matching
    // the built-in liquids (e.g. Tb = 373.15 K for water), rather than at the
    // OpenFOAM standard pressure pStd used for the enthalpy reference
    const scalar pAtm = 101325;

    if (coolpropUpdate(state, coolprop().PQ, pAtm, 0, msg))
    {
        return coolpropOutput(state, coolprop().T, fluid);
    }
    else
    {
        return coolpropOutput(state, coolprop().Tt, scalar(0));
    }
}


} // End namespace Foam


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

long Foam::coolpropProperties::flash
(
    const long state,
    const char* phase,
    const scalar Q,
    scalar& pMemo,
    scalar& TMemo,
    const scalar p,
    const scalar T
) const
{
    if (p == pMemo && T == TMemo)
    {
        return state;
    }

    char msg[coolpropMsgLen];

    coolpropPhase(state, phase);

    if (!coolpropUpdate(state, coolprop().PT, p, T, msg))
    {
        // Beyond the (metastable) single-phase region: fall back to the
        // saturated state of the imposed phase at the given temperature
        coolpropPhase(state, nullptr);

        if
        (
           !coolpropUpdate
            (
                state,
                coolprop().QT,
                Q,
                min(max(T, Tt()), 0.9999*Tc()),
                msg
            )
        )
        {
            FatalErrorInFunction
                << "CoolProp error for fluid " << fluid_
                << " at p = " << p << " Pa, T = " << T << " K: " << msg
                << exit(FatalError);
        }
    }

    pMemo = p;
    TMemo = T;

    return state;
}


long Foam::coolpropProperties::liquid(scalar p, scalar T) const
{
    return flash(liquidState_, "phase_liquid", 0, liquidP_, liquidT_, p, T);
}


long Foam::coolpropProperties::vapour(scalar p, scalar T) const
{
    return flash(vapourState_, "phase_gas", 1, vapourP_, vapourT_, p, T);
}


long Foam::coolpropProperties::saturation(scalar Q, scalar T) const
{
    char msg[coolpropMsgLen];

    const scalar Tsat = min(max(T, Tt()), 0.9999*Tc());

    if (!coolpropUpdate(saturationState_, coolprop().QT, Q, Tsat, msg))
    {
        FatalErrorInFunction
            << "CoolProp error for fluid " << fluid_
            << " on the saturation curve at T = " << Tsat << " K: " << msg
            << exit(FatalError);
    }

    return saturationState_;
}


long Foam::coolpropProperties::primaryState(scalar p, scalar T) const
{
    return phase_ == phaseType::gas ? vapour(p, T) : liquid(p, T);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::coolpropProperties::coolpropProperties
(
    const word& fluid,
    const long state,
    const phaseType phase
)
:
    liquidProperties
    (
        fluid,
        1000*coolpropOutput(state, coolprop().W, fluid),
        coolpropOutput(state, coolprop().Tc, fluid),
        coolpropOutput(state, coolprop().pc, fluid),
        1000/coolpropOutput(state, coolprop().rhoc, fluid),
        1000*coolpropOutput(state, coolprop().pc, fluid)
       /(
            coolpropOutput(state, coolprop().rhoc, fluid)
           *RR
           *coolpropOutput(state, coolprop().Tc, fluid)
        ),
        coolpropOutput(state, coolprop().Tt, scalar(0)),
        coolpropOutput(state, coolprop().pt, scalar(0)),
        coolpropTb(state, fluid),
        0,
        coolpropOutput(state, coolprop().omega, scalar(0)),
        0
    ),
    fluid_(fluid),
    liquidState_(state),
    vapourState_(newCoolpropState(fluid)),
    saturationState_(newCoolpropState(fluid)),
    liquidP_(-vGreat),
    liquidT_(-vGreat),
    vapourP_(-vGreat),
    vapourT_(-vGreat),
    pvT_(-vGreat),
    pv_(0),
    hlT_(-vGreat),
    hl_(0),
    sigmaT_(-vGreat),
    sigma_(0),
    D_("D", 15.0, 15.0, W(), 28),
    phase_(phase),
    hf_(0),
    ef_(0)
{
    // The offset between the absolute enthalpy/internal energy, which are
    // relative to the CoolProp reference state of the fluid, and the
    // sensible enthalpy/internal energy
    hf_ = ha(pStd, Tstd);
    ef_ = ea(pStd, Tstd);
}


Foam::coolpropProperties::coolpropProperties
(
    const word& fluid,
    const phaseType phase
)
:
    coolpropProperties(fluid, newCoolpropState(fluid), phase)
{}


Foam::coolpropProperties::coolpropProperties(const dictionary& dict)
:
    coolpropProperties
    (
        dict.lookupOrDefault<word>("fluid", dict.dictName()),
        phaseTypeNames.lookupOrDefault("phase", dict, phaseType::liquid)
    )
{}


Foam::coolpropProperties::coolpropProperties(const coolpropProperties& cpp)
:
    liquidProperties(cpp),
    fluid_(cpp.fluid_),
    liquidState_(newCoolpropState(fluid_)),
    vapourState_(newCoolpropState(fluid_)),
    saturationState_(newCoolpropState(fluid_)),
    liquidP_(-vGreat),
    liquidT_(-vGreat),
    vapourP_(-vGreat),
    vapourT_(-vGreat),
    pvT_(-vGreat),
    pv_(0),
    hlT_(-vGreat),
    hl_(0),
    sigmaT_(-vGreat),
    sigma_(0),
    D_(cpp.D_),
    phase_(cpp.phase_),
    hf_(cpp.hf_),
    ef_(cpp.ef_)
{}


// * * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::coolpropProperties> Foam::coolpropProperties::New
(
    const word& fluid,
    const phaseType phase
)
{
    if (debug)
    {
        InfoInFunction << "Constructing coolpropProperties" << endl;
    }

    return autoPtr<coolpropProperties>(new coolpropProperties(fluid, phase));
}


Foam::autoPtr<Foam::coolpropProperties> Foam::coolpropProperties::New
(
    const dictionary& dict
)
{
    if (debug)
    {
        InfoInFunction << "Constructing coolpropProperties" << endl;
    }

    return autoPtr<coolpropProperties>(new coolpropProperties(dict));
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::coolpropProperties::~coolpropProperties()
{
    freeCoolpropState(liquidState_);
    freeCoolpropState(vapourState_);
    freeCoolpropState(saturationState_);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::coolpropProperties::rho(scalar p, scalar T) const
{
    return coolpropOutput(primaryState(p, T), coolprop().rho, fluid_);
}


Foam::scalar Foam::coolpropProperties::alphav(scalar p, scalar T) const
{
    return coolpropOutput(primaryState(p, T), coolprop().alphav, fluid_);
}


Foam::scalar Foam::coolpropProperties::psi(scalar p, scalar T) const
{
    if (phase_ == phaseType::liquid)
    {
        return 0;
    }

    const long state = primaryState(p, T);

    return
        coolpropOutput(state, coolprop().rho, fluid_)
       *coolpropOutput(state, coolprop().kappaT, fluid_);
}


Foam::scalar Foam::coolpropProperties::CpMCv(scalar p, scalar T) const
{
    if (phase_ == phaseType::liquid)
    {
        return 0;
    }

    const long state = primaryState(p, T);

    return
        coolpropOutput(state, coolprop().Cp, fluid_)
      - coolpropOutput(state, coolprop().Cv, fluid_);
}


Foam::scalar Foam::coolpropProperties::Cp(scalar p, scalar T) const
{
    return coolpropOutput(primaryState(p, T), coolprop().Cp, fluid_);
}


Foam::scalar Foam::coolpropProperties::hs(scalar p, scalar T) const
{
    return ha(p, T) - hf();
}


Foam::scalar Foam::coolpropProperties::hf() const
{
    return hf_;
}


Foam::scalar Foam::coolpropProperties::ha(scalar p, scalar T) const
{
    return coolpropOutput(primaryState(p, T), coolprop().h, fluid_);
}


Foam::scalar Foam::coolpropProperties::ea(scalar p, scalar T) const
{
    if (phase_ == phaseType::liquid)
    {
        return ha(p, T);
    }

    return coolpropOutput(primaryState(p, T), coolprop().u, fluid_);
}


Foam::scalar Foam::coolpropProperties::es(scalar p, scalar T) const
{
    return ea(p, T) - ef_;
}


Foam::scalar Foam::coolpropProperties::s(scalar p, scalar T) const
{
    return coolpropOutput(primaryState(p, T), coolprop().s, fluid_);
}


Foam::scalar Foam::coolpropProperties::pv(scalar p, scalar T) const
{
    if (T != pvT_)
    {
        pvT_ = T;
        pv_ =
            T < Tc()
          ? coolpropOutput(saturation(0, T), coolprop().p, fluid_)
          : Pc();
    }

    return pv_;
}


Foam::scalar Foam::coolpropProperties::hl(scalar p, scalar T) const
{
    if (T != hlT_)
    {
        hlT_ = T;

        if (T < Tc())
        {
            // Both saturation endpoints come from the single Q = 0 flash
            const long state = saturation(0, T);

            hl_ =
                coolpropOutputSat(state, "gas", coolprop().h, fluid_)
              - coolpropOutputSat(state, "liquid", coolprop().h, fluid_);
        }
        else
        {
            hl_ = 0;
        }
    }

    return hl_;
}


Foam::scalar Foam::coolpropProperties::Cpg(scalar p, scalar T) const
{
    return coolpropOutput(vapour(p, T), coolprop().Cp0, fluid_);
}


Foam::scalar Foam::coolpropProperties::mu(scalar p, scalar T) const
{
    return coolpropOutput(primaryState(p, T), coolprop().mu, fluid_);
}


Foam::scalar Foam::coolpropProperties::mug(scalar p, scalar T) const
{
    return coolpropOutput(vapour(p, T), coolprop().mu, fluid_);
}


Foam::scalar Foam::coolpropProperties::kappa(scalar p, scalar T) const
{
    return coolpropOutput(primaryState(p, T), coolprop().kappa, fluid_);
}


Foam::scalar Foam::coolpropProperties::kappag(scalar p, scalar T) const
{
    return coolpropOutput(vapour(p, T), coolprop().kappa, fluid_);
}


Foam::scalar Foam::coolpropProperties::sigma(scalar p, scalar T) const
{
    if (T != sigmaT_)
    {
        sigmaT_ = T;
        sigma_ =
            T < Tc()
          ? coolpropOutput(saturation(0, T), coolprop().sigma, fluid_)
          : 0;
    }

    return sigma_;
}


Foam::scalar Foam::coolpropProperties::D(scalar p, scalar T) const
{
    return D_.value(p, T);
}


Foam::scalar Foam::coolpropProperties::D(scalar p, scalar T, scalar Wb) const
{
    return D_.value(p, T, Wb);
}


Foam::scalar Foam::coolpropProperties::pvInvert(scalar p) const
{
    if (p >= Pc())
    {
        return Tc();
    }

    char msg[coolpropMsgLen];

    if (!coolpropUpdate(saturationState_, coolprop().PQ, p, 0, msg))
    {
        if (debug)
        {
            WarningInFunction
                << "CoolProp saturation temperature error for fluid "
                << fluid_ << " at p = " << p << " Pa: " << msg
                << nl << endl;
        }

        // The -1 sentinel below the triple pressure mirrors the upstream
        // liquidProperties::pvInvert contract, which the Lagrangian
        // phase-change models are written against; consumers that must
        // not receive it clamp at their boundary (see the Tsat case of
        // Function1s::coolprop::value)
        return -1;
    }

    return coolpropOutput(saturationState_, coolprop().T, fluid_);
}


// * * * * * * * * * * * * * * * * * * I-O  * * * * * * * * * * * * * * * * //

void Foam::coolpropProperties::write(Ostream& os) const
{
    writeEntry(os, "fluid", fluid_);

    if (phase_ == phaseType::gas)
    {
        writeEntry(os, "phase", phaseTypeNames[phase_]);
    }

    liquidProperties::write(os);
}


// ************************************************************************* //
