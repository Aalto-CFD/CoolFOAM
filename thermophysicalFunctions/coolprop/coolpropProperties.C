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


//- Minimum pressure above saturation of the liquid states of a tabulated
//  backend [Pa]: closer to the saturation curve the table can be wrong
//  (e.g. by up to 100% in density for SVDSBTL&REFPROP::H2O within 10% of
//  psat below 275 K), so there the liquid is extrapolated from this far
//  above it, which is accurate to psi*dp ~ 1e-5 relative for a liquid
static const scalar tabulatedLiquidDp = 1e4;

//- Relative pressure below saturation of the vapour fallback state of a
//  tabulated backend
static const scalar tabulatedVapourDp = 1e-4;

//- Maximum relative density change psi*dp/rho of the linear extrapolation
//  of a tabulated backend beyond saturation, beyond which it is held.
//  Without it, the energy of the vapour of the liquid cells of a VoF case
//  runs away (e.g. by -1.3e7 J/kg for water at 300 K compressed to 30 MPa,
//  extrapolated with the derivatives of the dilute saturated vapour), and
//  so does the liquid next to the critical point (e.g. to 2.9e7 J/kg from
//  psat = 22 MPa down to 1 kPa). The vapour density is not held though, to
//  stay consistent with psi, which the pressure equation relies on.
static const scalar tabulatedExtrapolation = 0.1;


//- Input-pair and parameter indices of the CoolProp C interface,
//  resolved once on first use
struct coolpropIndices
{
    const long PT, QT, PQ, DT;
    const long T, p, rho, Cp, Cp0, Cv, u, h, s, mu, kappa, sigma, alphav,
        kappaT;
    const long Tc, pc, rhoc, Tt, pt, W, omega;

    coolpropIndices()
    :
        PT(get_input_pair_index("PT_INPUTS")),
        QT(get_input_pair_index("QT_INPUTS")),
        PQ(get_input_pair_index("PQ_INPUTS")),
        DT(get_input_pair_index("DmassT_INPUTS")),
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


//- Return the backend of the given fluid specification, i.e. its optional
//  "Backend::" prefix, e.g. "REFPROP" for "REFPROP::H2O"; defaults to "HEOS"
static std::string coolpropBackend(const word& fluid)
{
    const std::string::size_type sep = fluid.find("::");

    return
        sep == std::string::npos ? std::string("HEOS") : fluid.substr(0, sep);
}


//- Return the fluid name of the given fluid specification, i.e. without its
//  optional "Backend::" prefix
static std::string coolpropName(const word& fluid)
{
    const std::string::size_type sep = fluid.find("::");

    return
        sep == std::string::npos
      ? static_cast<const std::string&>(fluid)
      : fluid.substr(sep + 2);
}


//- Is the backend of the given fluid specification an SVDSBTL table,
//  e.g. "SVDSBTL&REFPROP::H2O"
static bool coolpropTabulated(const word& fluid)
{
    return coolpropBackend(fluid).compare(0, 7, "SVDSBTL") == 0;
}


//- Return the fluid specification of the source backend the table of the
//  given tabulated fluid specification is built from, e.g. "REFPROP::H2O"
//  for "SVDSBTL&REFPROP::H2O"; defaults to "HEOS"
static word coolpropSource(const word& fluid)
{
    const std::string backend = coolpropBackend(fluid);
    const std::string::size_type sep = backend.find('&');

    const std::string source =
        sep == std::string::npos
      ? std::string("HEOS")
      : backend.substr(sep + 1);

    return source + "::" + coolpropName(fluid);
}


//- Construct a new state for the given fluid specification, splitting an
//  optional "Backend::" prefix, e.g. "REFPROP::H2O"; defaults to "HEOS"
static long newCoolpropState(const word& fluid)
{
    const std::string backend = coolpropBackend(fluid);
    const std::string name = coolpropName(fluid);

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
        // The conditions of the state, where available, to locate the error
        char buf[coolpropMsgLen];
        const scalar p =
            AbstractState_keyed_output
            (
                state, coolprop().p, &errcode, buf, coolpropMsgLen
            );
        const scalar T =
            AbstractState_keyed_output
            (
                state, coolprop().T, &errcode, buf, coolpropMsgLen
            );

        FatalErrorInFunction
            << "CoolProp error for fluid " << fluid
            << " at p = " << p << " Pa, T = " << T << " K: " << msg
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


//- Acentric factor, from the source backend for a tabulated fluid (whose
//  table does not provide it), or zero if it is not available
static scalar coolpropAcentric(const long state, const word& fluid)
{
    if (coolpropTabulated(fluid))
    {
        const long source = newCoolpropState(coolpropSource(fluid));
        const scalar omega =
            coolpropOutput(source, coolprop().omega, scalar(0));
        freeCoolpropState(source);

        return omega;
    }

    return coolpropOutput(state, coolprop().omega, scalar(0));
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

    bool single =
        (
            !tabulated_
         || Q != 0
         || T >= Tc()
         || p >= pv(p, T) + tabulatedLiquidDp
        )
     && coolpropUpdate(state, coolprop().PT, p, T, msg)
     && (!tabulated_ || imposedSide(state, Q, T));

    if (!single && tabulated_ && T >= Tc())
    {
        // Beyond the table at a supercritical temperature, i.e. below its
        // lowest pressure, the triple pressure (e.g. air below 5.3 kPa):
        // there is no saturation curve to fall back to, so fall back to the
        // table state at its lowest pressure, from which the extrapolation
        // continues
        single =
            coolpropUpdate(state, coolprop().PT, (1 + 1e-3)*Pt(), T, msg)
         && imposedSide(state, Q, T);
    }

    if (!single && tabulated_)
    {
        // A table has neither a metastable extension nor the single-phase
        // properties (e.g. viscosity) of the saturation curve itself: fall
        // back to its state off the saturation curve on the side of the
        // imposed phase, at the temperature limited to [Tt, Tc], and just
        // above Tt so that the vapour state stays above the triple
        // pressure, the lower limit of the table
        //
        // If the table misclassifies that state, step further off the
        // saturation curve (the extrapolation from it accounts for the
        // distance); a two-phase state would lack the transport properties
        const scalar Tsat = min(max(T, (1 + 1e-4)*Tt()), 0.9999*Tc());
        const scalar psat = pv(p, Tsat);

        scalar pa = psat;

        for (label i = 0; !single && i < 3; i++)
        {
            const scalar scale = pow(scalar(10), i);

            pa =
                Q == 0
              ? psat + scale*tabulatedLiquidDp
              : psat*(1 - scale*tabulatedVapourDp);

            single =
                coolpropUpdate(state, coolprop().PT, pa, Tsat, msg)
             && imposedSide(state, Q, Tsat);
        }

        if (!single)
        {
            FatalErrorInFunction
                << "No " << (Q == 0 ? "liquid" : "vapour")
                << " state of the tabulated fluid " << fluid_
                << " at p = " << p << " Pa, T = " << T
                << " K, nor next to the saturation curve at T = " << Tsat
                << " K, psat = " << psat << " Pa (last tried p = " << pa
                << " Pa: " << msg << ")"
                << exit(FatalError);
        }
    }

    if (!single)
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


long Foam::coolpropProperties::saturation
(
    const long state,
    scalar Q,
    scalar T
) const
{
    char msg[coolpropMsgLen];

    const scalar Tsat = min(max(T, Tt()), 0.9999*Tc());

    if (!coolpropUpdate(state, coolprop().QT, Q, Tsat, msg))
    {
        FatalErrorInFunction
            << "CoolProp error for fluid " << fluid_
            << " on the saturation curve at T = " << Tsat << " K: " << msg
            << exit(FatalError);
    }

    return state;
}


long Foam::coolpropProperties::saturation(scalar Q, scalar T) const
{
    return saturation(saturationState_, Q, T);
}


bool Foam::coolpropProperties::imposedSide
(
    const long state,
    const scalar Q,
    const scalar T
) const
{
    // Beyond its domain (e.g. below the triple point) the table reports a
    // successful update of a NaN state
    const scalar rho = coolpropOutput(state, coolprop().rho, fluid_);

    if (!std::isfinite(rho))
    {
        return false;
    }

    if (T >= Tc())
    {
        return true;
    }

    // Stable (not metastable) states are on the liquid side of the
    // saturation curve if and only if denser than the critical density
    const scalar rhoc = W()/Vc();

    return Q == 0 ? rho > rhoc : rho < rhoc;
}


void Foam::coolpropProperties::derivatives
(
    const scalar p,
    const scalar T
) const
{
    if (p == derivP_ && T == derivT_)
    {
        return;
    }

    // The resolved table state: (p, T) itself, or the fallback state next
    // to the saturation curve at the limited temperature
    const long state = primaryState(p, T);
    const scalar rho = coolpropOutput(state, coolprop().rho, fluid_);
    const scalar Ts = coolpropOutput(state, coolprop().T, fluid_);

    // Evaluate the source backend at the table density and temperature
    // with the imposed phase, which is explicit (no flash iteration) for a
    // Helmholtz-energy backend; a source without density-temperature
    // inputs (e.g. IF97) is evaluated at the table pressure instead
    char msg[coolpropMsgLen];

    coolpropPhase
    (
        derivState_,
        phase_ == phaseType::gas ? "phase_gas" : "phase_liquid"
    );

    if
    (
       !coolpropUpdate(derivState_, coolprop().DT, rho, Ts, msg)
    && !coolpropUpdate
        (
            derivState_,
            coolprop().PT,
            coolpropOutput(state, coolprop().p, fluid_),
            Ts,
            msg
        )
    )
    {
        FatalErrorInFunction
            << "CoolProp error for fluid " << coolpropSource(fluid_)
            << " at rho = " << rho << " kg/m^3, T = " << Ts << " K: " << msg
            << exit(FatalError);
    }

    const word source(coolpropSource(fluid_));

    Cp_ = coolpropOutput(derivState_, coolprop().Cp, source);
    alphav_ = coolpropOutput(derivState_, coolprop().alphav, source);
    psi_ =
        coolpropOutput(derivState_, coolprop().rho, source)
       *coolpropOutput(derivState_, coolprop().kappaT, source);
    CpMCv_ = Cp_ - coolpropOutput(derivState_, coolprop().Cv, source);

    derivP_ = p;
    derivT_ = T;
}


Foam::scalar Foam::coolpropProperties::metastableDp
(
    const long state,
    const scalar p,
    const scalar T,
    const bool energy
) const
{
    if (!tabulated_)
    {
        return 0;
    }

    const scalar dp = p - coolpropOutput(state, coolprop().p, fluid_);

    if (dp == 0 || (!energy && phase_ == phaseType::gas))
    {
        return dp;
    }

    derivatives(p, T);

    const scalar dpMax =
        tabulatedExtrapolation
       *coolpropOutput(state, coolprop().rho, fluid_)/max(psi_, small);

    return min(max(dp, -dpMax), dpMax);
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
        coolpropAcentric(state, fluid),
        0
    ),
    fluid_(fluid),
    tabulated_(coolpropTabulated(fluid)),
    liquidState_(state),
    vapourState_(newCoolpropState(fluid)),
    saturationState_(newCoolpropState(fluid)),
    derivState_
    (
        tabulated_ ? newCoolpropState(coolpropSource(fluid)) : -1
    ),
    sourceState_(tabulated_ ? newCoolpropState(coolpropSource(fluid)) : -1),
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
    derivP_(-vGreat),
    derivT_(-vGreat),
    Cp_(0),
    alphav_(0),
    psi_(0),
    CpMCv_(0),
    D_("D", 15.0, 15.0, W(), 28),
    phase_(phase),
    hf_(0),
    ef_(0)
{
    // The offset between the absolute enthalpy/internal energy, which are
    // relative to the CoolProp reference state of the fluid, and the
    // sensible enthalpy/internal energy. Taken from the resolved state
    // without the metastable extrapolation of a tabulated backend, so that
    // for a role beyond saturation at standard conditions (e.g. the gas) the
    // reference is the saturated state at Tstd, as for the other backends.
    const long stdState = primaryState(pStd, Tstd);
    hf_ = coolpropOutput(stdState, coolprop().h, fluid_);
    ef_ = coolpropOutput(stdState, coolprop().u, fluid_);
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
    tabulated_(cpp.tabulated_),
    liquidState_(newCoolpropState(fluid_)),
    vapourState_(newCoolpropState(fluid_)),
    saturationState_(newCoolpropState(fluid_)),
    derivState_
    (
        tabulated_ ? newCoolpropState(coolpropSource(fluid_)) : -1
    ),
    sourceState_(tabulated_ ? newCoolpropState(coolpropSource(fluid_)) : -1),
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
    derivP_(-vGreat),
    derivT_(-vGreat),
    Cp_(0),
    alphav_(0),
    psi_(0),
    CpMCv_(0),
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

    if (tabulated_)
    {
        freeCoolpropState(derivState_);
        freeCoolpropState(sourceState_);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::coolpropProperties::rho(scalar p, scalar T) const
{
    const long state = primaryState(p, T);
    const scalar rho = coolpropOutput(state, coolprop().rho, fluid_);
    const scalar dp = metastableDp(state, p, T, false);

    if (dp == 0)
    {
        return rho;
    }

    // Floored against the extrapolation of the gas to (unphysically) low
    // or negative pressures
    derivatives(p, T);
    return max(rho + psi_*dp, 1e-3*rho);
}


Foam::scalar Foam::coolpropProperties::alphav(scalar p, scalar T) const
{
    if (tabulated_)
    {
        derivatives(p, T);
        return alphav_;
    }

    return coolpropOutput(primaryState(p, T), coolprop().alphav, fluid_);
}


Foam::scalar Foam::coolpropProperties::psi(scalar p, scalar T) const
{
    if (tabulated_)
    {
        derivatives(p, T);
        return psi_;
    }

    const long state = primaryState(p, T);

    return
        coolpropOutput(state, coolprop().rho, fluid_)
       *coolpropOutput(state, coolprop().kappaT, fluid_);
}


Foam::scalar Foam::coolpropProperties::CpMCv(scalar p, scalar T) const
{
    if (tabulated_)
    {
        derivatives(p, T);
        return CpMCv_;
    }

    const long state = primaryState(p, T);

    return
        coolpropOutput(state, coolprop().Cp, fluid_)
      - coolpropOutput(state, coolprop().Cv, fluid_);
}


Foam::scalar Foam::coolpropProperties::Cp(scalar p, scalar T) const
{
    if (tabulated_)
    {
        derivatives(p, T);
        return Cp_;
    }

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
    const long state = primaryState(p, T);
    const scalar h = coolpropOutput(state, coolprop().h, fluid_);
    const scalar dp = metastableDp(state, p, T, true);

    if (dp == 0)
    {
        return h;
    }

    // (dh/dp)_T = (1 - T*alphav)/rho
    derivatives(p, T);
    return
        h
      + (1 - coolpropOutput(state, coolprop().T, fluid_)*alphav_)*dp
       /coolpropOutput(state, coolprop().rho, fluid_);
}


Foam::scalar Foam::coolpropProperties::ea(scalar p, scalar T) const
{
    const long state = primaryState(p, T);
    const scalar e = coolpropOutput(state, coolprop().u, fluid_);
    const scalar dp = metastableDp(state, p, T, true);

    if (dp == 0)
    {
        return e;
    }

    // (de/dp)_T = (p*kappaT - T*alphav)/rho, with kappaT = psi/rho
    derivatives(p, T);
    const scalar rhoa = coolpropOutput(state, coolprop().rho, fluid_);

    return
        e
      + (
            coolpropOutput(state, coolprop().p, fluid_)*psi_/rhoa
          - coolpropOutput(state, coolprop().T, fluid_)*alphav_
        )*dp/rhoa;
}


Foam::scalar Foam::coolpropProperties::es(scalar p, scalar T) const
{
    return ea(p, T) - ef_;
}


Foam::scalar Foam::coolpropProperties::s(scalar p, scalar T) const
{
    const long state = primaryState(p, T);
    const scalar s = coolpropOutput(state, coolprop().s, fluid_);
    const scalar dp = metastableDp(state, p, T, true);

    if (dp == 0)
    {
        return s;
    }

    // (ds/dp)_T = -alphav/rho
    derivatives(p, T);
    return s - alphav_*dp/coolpropOutput(state, coolprop().rho, fluid_);
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

        if (T < Tc() && tabulated_)
        {
            // A tabulated backend does not provide the saturation endpoints
            // of a flash, so evaluate each from its own
            const scalar hv =
                coolpropOutput(saturation(1, T), coolprop().h, fluid_);

            hl_ = hv - coolpropOutput(saturation(0, T), coolprop().h, fluid_);
        }
        else if (T < Tc())
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
    if (tabulated_)
    {
        // The ideal-gas heat capacity depends on T only: evaluate it from
        // the source backend in the dilute-gas limit, which exists at any T
        char msg[coolpropMsgLen];

        if (!coolpropUpdate(sourceState_, coolprop().DT, 1e-6, T, msg))
        {
            FatalErrorInFunction
                << "CoolProp error for fluid " << coolpropSource(fluid_)
                << " at T = " << T << " K: " << msg
                << exit(FatalError);
        }

        return coolpropOutput(sourceState_, coolprop().Cp0, fluid_);
    }

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
          ? coolpropOutput
            (
                saturation(tabulated_ ? sourceState_ : saturationState_, 0, T),
                coolprop().sigma,
                fluid_
            )
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
