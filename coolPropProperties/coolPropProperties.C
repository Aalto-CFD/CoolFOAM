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

#include "coolPropProperties.H"

#include "thermodynamicConstants.H"
using namespace Foam::constant::thermodynamic;

#include "CoolPropLib.h"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(coolPropProperties, 0);

    const NamedEnum<coolPropProperties::phaseType, 2>
    coolPropProperties::phaseTypeNames
    {
        "liquid",
        "gas"
    };
}

// * * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * //

namespace Foam
{

//- Length of the message buffers of the CoolProp C interface
static const long coolPropMsgLen = 1000;


//- Input-pair and parameter indices of the CoolProp C interface,
//  resolved once on first use
struct coolPropIndices
{
    const long PT, QT, PQ;
    const long T, p, rho, Cp, Cp0, Cv, u, h, s, mu, kappa, sigma, alphav,
        kappaT;
    const long Tc, pc, rhoc, Tt, pt, W, omega;

    coolPropIndices()
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


static const coolPropIndices& coolProp()
{
    static const coolPropIndices indices;
    return indices;
}


//- Construct a new state for the given fluid specification, splitting an
//  optional "Backend::" prefix, e.g. "REFPROP::H2O"; defaults to "HEOS"
static long newCoolPropState(const word& fluid)
{
    const std::string::size_type sep = fluid.find("::");

    const std::string backend =
        sep == std::string::npos ? std::string("HEOS") : fluid.substr(0, sep);

    const std::string name =
        sep == std::string::npos
      ? static_cast<const std::string&>(fluid)
      : fluid.substr(sep + 2);

    long errcode = 0;
    char msg[coolPropMsgLen];

    const long state = AbstractState_factory
    (
        backend.c_str(),
        name.c_str(),
        &errcode,
        msg,
        coolPropMsgLen
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
static void freeCoolPropState(const long state)
{
    long errcode = 0;
    char msg[coolPropMsgLen];

    AbstractState_free(state, &errcode, msg, coolPropMsgLen);
}


//- Update the given state, returning success and the error message
static bool coolPropUpdate
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
        coolPropMsgLen
    );

    return errcode == 0;
}


//- Impose the given phase, or lift the imposition for nullptr;
//  ignored for backends without phase imposition support
static void coolPropPhase(const long state, const char* phase)
{
    long errcode = 0;
    char msg[coolPropMsgLen];

    if (phase)
    {
        AbstractState_specify_phase
        (
            state,
            phase,
            &errcode,
            msg,
            coolPropMsgLen
        );
    }
    else
    {
        AbstractState_unspecify_phase(state, &errcode, msg, coolPropMsgLen);
    }
}


//- Return the keyed output of the given state of the given fluid
static scalar coolPropOutput
(
    const long state,
    const long key,
    const word& fluid
)
{
    long errcode = 0;
    char msg[coolPropMsgLen];

    const scalar value =
        AbstractState_keyed_output(state, key, &errcode, msg, coolPropMsgLen);

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
static scalar coolPropOutput
(
    const long state,
    const long key,
    const scalar deflt
)
{
    long errcode = 0;
    char msg[coolPropMsgLen];

    const scalar value =
        AbstractState_keyed_output(state, key, &errcode, msg, coolPropMsgLen);

    return errcode == 0 ? value : deflt;
}


//- Return the keyed output of the saturated phase ("liquid" or "gas") of the
//  given two-phase state of the given fluid, avoiding a second flash
static scalar coolPropOutputSat
(
    const long state,
    const char* saturatedPhase,
    const long key,
    const word& fluid
)
{
    long errcode = 0;
    char msg[coolPropMsgLen];

    const scalar value = AbstractState_keyed_output_satState
    (
        state,
        saturatedPhase,
        key,
        &errcode,
        msg,
        coolPropMsgLen
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
static scalar coolPropTb(const long state, const word& fluid)
{
    char msg[coolPropMsgLen];

    // Normal boiling point is defined at one standard atmosphere, matching
    // the built-in liquids (e.g. Tb = 373.15 K for water), rather than at the
    // OpenFOAM standard pressure pStd used for the enthalpy reference
    const scalar pAtm = 101325;

    if (coolPropUpdate(state, coolProp().PQ, pAtm, 0, msg))
    {
        return coolPropOutput(state, coolProp().T, fluid);
    }
    else
    {
        return coolPropOutput(state, coolProp().Tt, scalar(0));
    }
}


} // End namespace Foam


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

long Foam::coolPropProperties::flash
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

    char msg[coolPropMsgLen];

    coolPropPhase(state, phase);

    if (!coolPropUpdate(state, coolProp().PT, p, T, msg))
    {
        // Beyond the (metastable) single-phase region: fall back to the
        // saturated state of the imposed phase at the given temperature
        coolPropPhase(state, nullptr);

        if
        (
           !coolPropUpdate
            (
                state,
                coolProp().QT,
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


long Foam::coolPropProperties::liquid(scalar p, scalar T) const
{
    return flash(liquidState_, "phase_liquid", 0, liquidP_, liquidT_, p, T);
}


long Foam::coolPropProperties::vapour(scalar p, scalar T) const
{
    return flash(vapourState_, "phase_gas", 1, vapourP_, vapourT_, p, T);
}


long Foam::coolPropProperties::saturation(scalar Q, scalar T) const
{
    char msg[coolPropMsgLen];

    const scalar Tsat = min(max(T, Tt()), 0.9999*Tc());

    if (!coolPropUpdate(saturationState_, coolProp().QT, Q, Tsat, msg))
    {
        FatalErrorInFunction
            << "CoolProp error for fluid " << fluid_
            << " on the saturation curve at T = " << Tsat << " K: " << msg
            << exit(FatalError);
    }

    return saturationState_;
}


long Foam::coolPropProperties::primaryState(scalar p, scalar T) const
{
    return phase_ == phaseType::gas ? vapour(p, T) : liquid(p, T);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::coolPropProperties::coolPropProperties
(
    const word& fluid,
    const long state,
    const phaseType phase
)
:
    liquidProperties
    (
        fluid,
        1000*coolPropOutput(state, coolProp().W, fluid),
        coolPropOutput(state, coolProp().Tc, fluid),
        coolPropOutput(state, coolProp().pc, fluid),
        1000/coolPropOutput(state, coolProp().rhoc, fluid),
        1000*coolPropOutput(state, coolProp().pc, fluid)
       /(
            coolPropOutput(state, coolProp().rhoc, fluid)
           *RR
           *coolPropOutput(state, coolProp().Tc, fluid)
        ),
        coolPropOutput(state, coolProp().Tt, scalar(0)),
        coolPropOutput(state, coolProp().pt, scalar(0)),
        coolPropTb(state, fluid),
        0,
        coolPropOutput(state, coolProp().omega, scalar(0)),
        0
    ),
    fluid_(fluid),
    liquidState_(state),
    vapourState_(newCoolPropState(fluid)),
    saturationState_(newCoolPropState(fluid)),
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


Foam::coolPropProperties::coolPropProperties
(
    const word& fluid,
    const phaseType phase
)
:
    coolPropProperties(fluid, newCoolPropState(fluid), phase)
{}


Foam::coolPropProperties::coolPropProperties(const dictionary& dict)
:
    coolPropProperties
    (
        dict.lookupOrDefault<word>("fluid", dict.dictName()),
        phaseTypeNames.lookupOrDefault("phase", dict, phaseType::liquid)
    )
{}


Foam::coolPropProperties::coolPropProperties(const coolPropProperties& cpp)
:
    liquidProperties(cpp),
    fluid_(cpp.fluid_),
    liquidState_(newCoolPropState(fluid_)),
    vapourState_(newCoolPropState(fluid_)),
    saturationState_(newCoolPropState(fluid_)),
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

Foam::autoPtr<Foam::coolPropProperties> Foam::coolPropProperties::New
(
    const word& fluid,
    const phaseType phase
)
{
    if (debug)
    {
        InfoInFunction << "Constructing coolPropProperties" << endl;
    }

    return autoPtr<coolPropProperties>(new coolPropProperties(fluid, phase));
}


Foam::autoPtr<Foam::coolPropProperties> Foam::coolPropProperties::New
(
    const dictionary& dict
)
{
    if (debug)
    {
        InfoInFunction << "Constructing coolPropProperties" << endl;
    }

    return autoPtr<coolPropProperties>(new coolPropProperties(dict));
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::coolPropProperties::~coolPropProperties()
{
    freeCoolPropState(liquidState_);
    freeCoolPropState(vapourState_);
    freeCoolPropState(saturationState_);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::coolPropProperties::rho(scalar p, scalar T) const
{
    return coolPropOutput(primaryState(p, T), coolProp().rho, fluid_);
}


Foam::scalar Foam::coolPropProperties::alphav(scalar p, scalar T) const
{
    return coolPropOutput(primaryState(p, T), coolProp().alphav, fluid_);
}


Foam::scalar Foam::coolPropProperties::psi(scalar p, scalar T) const
{
    if (phase_ == phaseType::liquid)
    {
        return 0;
    }

    const long state = primaryState(p, T);

    return
        coolPropOutput(state, coolProp().rho, fluid_)
       *coolPropOutput(state, coolProp().kappaT, fluid_);
}


Foam::scalar Foam::coolPropProperties::CpMCv(scalar p, scalar T) const
{
    if (phase_ == phaseType::liquid)
    {
        return 0;
    }

    const long state = primaryState(p, T);

    return
        coolPropOutput(state, coolProp().Cp, fluid_)
      - coolPropOutput(state, coolProp().Cv, fluid_);
}


Foam::scalar Foam::coolPropProperties::Cp(scalar p, scalar T) const
{
    return coolPropOutput(primaryState(p, T), coolProp().Cp, fluid_);
}


Foam::scalar Foam::coolPropProperties::hs(scalar p, scalar T) const
{
    return ha(p, T) - hf();
}


Foam::scalar Foam::coolPropProperties::hf() const
{
    return hf_;
}


Foam::scalar Foam::coolPropProperties::ha(scalar p, scalar T) const
{
    return coolPropOutput(primaryState(p, T), coolProp().h, fluid_);
}


Foam::scalar Foam::coolPropProperties::ea(scalar p, scalar T) const
{
    if (phase_ == phaseType::liquid)
    {
        return ha(p, T);
    }

    return coolPropOutput(primaryState(p, T), coolProp().u, fluid_);
}


Foam::scalar Foam::coolPropProperties::es(scalar p, scalar T) const
{
    return ea(p, T) - ef_;
}


Foam::scalar Foam::coolPropProperties::s(scalar p, scalar T) const
{
    return coolPropOutput(primaryState(p, T), coolProp().s, fluid_);
}


Foam::scalar Foam::coolPropProperties::pv(scalar p, scalar T) const
{
    if (T != pvT_)
    {
        pvT_ = T;
        pv_ =
            T < Tc()
          ? coolPropOutput(saturation(0, T), coolProp().p, fluid_)
          : Pc();
    }

    return pv_;
}


Foam::scalar Foam::coolPropProperties::hl(scalar p, scalar T) const
{
    if (T != hlT_)
    {
        hlT_ = T;

        if (T < Tc())
        {
            // Both saturation endpoints come from the single Q = 0 flash
            const long state = saturation(0, T);

            hl_ =
                coolPropOutputSat(state, "gas", coolProp().h, fluid_)
              - coolPropOutputSat(state, "liquid", coolProp().h, fluid_);
        }
        else
        {
            hl_ = 0;
        }
    }

    return hl_;
}


Foam::scalar Foam::coolPropProperties::Cpg(scalar p, scalar T) const
{
    return coolPropOutput(vapour(p, T), coolProp().Cp0, fluid_);
}


Foam::scalar Foam::coolPropProperties::mu(scalar p, scalar T) const
{
    return coolPropOutput(primaryState(p, T), coolProp().mu, fluid_);
}


Foam::scalar Foam::coolPropProperties::mug(scalar p, scalar T) const
{
    return coolPropOutput(vapour(p, T), coolProp().mu, fluid_);
}


Foam::scalar Foam::coolPropProperties::kappa(scalar p, scalar T) const
{
    return coolPropOutput(primaryState(p, T), coolProp().kappa, fluid_);
}


Foam::scalar Foam::coolPropProperties::kappag(scalar p, scalar T) const
{
    return coolPropOutput(vapour(p, T), coolProp().kappa, fluid_);
}


Foam::scalar Foam::coolPropProperties::sigma(scalar p, scalar T) const
{
    if (T != sigmaT_)
    {
        sigmaT_ = T;
        sigma_ =
            T < Tc()
          ? coolPropOutput(saturation(0, T), coolProp().sigma, fluid_)
          : 0;
    }

    return sigma_;
}


Foam::scalar Foam::coolPropProperties::D(scalar p, scalar T) const
{
    return D_.value(p, T);
}


Foam::scalar Foam::coolPropProperties::D(scalar p, scalar T, scalar Wb) const
{
    return D_.value(p, T, Wb);
}


Foam::scalar Foam::coolPropProperties::pvInvert(scalar p) const
{
    if (p >= Pc())
    {
        return Tc();
    }

    char msg[coolPropMsgLen];

    if (!coolPropUpdate(saturationState_, coolProp().PQ, p, 0, msg))
    {
        if (debug)
        {
            WarningInFunction
                << "CoolProp saturation temperature error for fluid "
                << fluid_ << " at p = " << p << " Pa: " << msg
                << nl << endl;
        }

        return -1;
    }

    return coolPropOutput(saturationState_, coolProp().T, fluid_);
}


// * * * * * * * * * * * * * * * * * * I-O  * * * * * * * * * * * * * * * * //

void Foam::coolPropProperties::write(Ostream& os) const
{
    writeEntry(os, "fluid", fluid_);

    if (phase_ == phaseType::gas)
    {
        writeEntry(os, "phase", phaseTypeNames[phase_]);
    }

    liquidProperties::write(os);
}


// ************************************************************************* //
