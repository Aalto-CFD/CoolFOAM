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

Application
    Test-coolpropThermo

Authors
    Stanislau Stasheuski, Aalto University, 2026.
    stanislau.stasheuski[at]aalto.fi

Description
    Tests the Coolprop specie-layer thermodynamics stack: the composed
    type name, dictionary construction (including the $equationOfState
    macro), agreement with a directly-constructed coolpropProperties,
    the write/re-read round trip, and the Coolprop Function1.

\*---------------------------------------------------------------------------*/

#include "coolpropProperties.H"
#include "Coolprop.H"
#include "CoolpropThermo.H"
#include "CoolpropTransport.H"
#include "coolprop.H"

#include "specie.H"
#include "sensibleInternalEnergy.H"
#include "thermo.H"
#include "units.H"

#include "IStringStream.H"
#include "OStringStream.H"

using namespace Foam;

typedef
    CoolpropTransport
    <
        species::thermo
        <
            CoolpropThermo<Coolprop<specie>>,
            sensibleInternalEnergy
        >
    >
    thermoPhysics;


static int nFail = 0;

void check(const word& what, const scalar actual, const scalar expected)
{
    const scalar tol = 1e-10*max(mag(expected), small);

    if (mag(actual - expected) > tol)
    {
        Info<< "FAIL " << what << ": " << actual
            << " != " << expected << endl;
        nFail++;
    }
    else
    {
        Info<< "pass " << what << " = " << actual << endl;
    }
}


scalar function1Value(const dictionary& dict, const word& name, const scalar x)
{
    return
        Function1<scalar>::New
        (
            name,
            {units::any, units::any},
            dict
        )().value(x);
}


int main(int argc, char *argv[])
{
    const string mixtureString
    (
        "mixture"
        "{"
        "    equationOfState"
        "    {"
        "        type            coolprop;"
        "        fluid           H2O;"
        "    }"
        "    thermodynamics  { $equationOfState; }"
        "    transport"
        "    {"
        "        mu              { $equationOfState; }"
        "        kappa           { $equationOfState; }"
        "    }"
        "}"
    );

    IStringStream mixtureStream(mixtureString);
    const dictionary dict(mixtureStream);

    // Composed type name must match what basicThermo::New assembles from
    // the thermoType entries
    check
    (
        "typeName",
        thermoPhysics::typeName()
     == "coolprop<coolprop<coolprop<specie>>,sensibleInternalEnergy>",
        1
    );

    const thermoPhysics mix("mixture", dict.subDict("mixture"));
    const coolpropProperties props("H2O");

    const scalar p = 1e5, T = 300;

    // Molecular weight auto-derived from Coolprop, no specie sub-dict
    check("W", mix.W(), props.W());

    // Every layer must agree with the direct evaluator
    check("rho", mix.rho(p, T), props.rho(p, T));
    check("psi", mix.psi(p, T), props.psi(p, T));
    check("CpMCv", mix.CpMCv(p, T), props.CpMCv(p, T));
    check("alphav", mix.alphav(p, T), props.alphav(p, T));
    check("Cp", mix.Cp(p, T), props.Cp(p, T));
    check("Cv", mix.Cv(p, T), props.Cp(p, T) - props.CpMCv(p, T));
    check("ha", mix.ha(p, T), props.ha(p, T));
    check("hs", mix.hs(p, T), props.hs(p, T));
    check("hf", mix.hf(), props.hf());
    check("es", mix.es(p, T), props.es(p, T));
    check("ea", mix.ea(p, T), props.ea(p, T));
    check("s", mix.s(p, T), props.s(p, T));
    check("mu", mix.mu(p, T), props.mu(p, T));
    check("kappa", mix.kappa(p, T), props.kappa(p, T));

    // limit() clamps to the backend's genuine-state range (liquid role)
    check("limit(min)", mix.limit(1), props.Tt());
    check("limit(mid)", mix.limit(T), T);
    check("limit(max)", mix.limit(1e5), 0.9999*props.Tc());

    // Write/re-read round trip reconstructs the same physics
    {
        OStringStream os;
        mix.write(os);

        IStringStream restream(os.str());
        const dictionary redict(restream);
        const thermoPhysics remix("mixture", redict);

        check("round-trip rho", remix.rho(p, T), mix.rho(p, T));
        check("round-trip Cp", remix.Cp(p, T), mix.Cp(p, T));
        check("round-trip mu", remix.mu(p, T), mix.mu(p, T));
    }

    // Reference-style options: an honoured specie sub-dictionary and the
    // hf/sf reference offsets of the thermodynamics model
    {
        IStringStream offsetStream
        (
            "mixture"
            "{"
            "    specie          { molWeight 42; }"
            "    equationOfState"
            "    {"
            "        type            coolprop;"
            "        fluid           H2O;"
            "    }"
            "    thermodynamics  { $equationOfState; hf 1000; sf 10; }"
            "    transport"
            "    {"
            "        mu              { $equationOfState; }"
            "        kappa           { $equationOfState; }"
            "    }"
            "}"
        );
        const dictionary offsetDict(offsetStream);
        const thermoPhysics offmix("mixture", offsetDict.subDict("mixture"));

        check("specie W", offmix.W(), 42);
        check("hf offset ha", offmix.ha(p, T), props.ha(p, T) + 1000);
        check("hf offset hf", offmix.hf(), props.hf() + 1000);
        check("hf offset hs", offmix.hs(p, T), props.hs(p, T));
        check("hf offset ea", offmix.ea(p, T), props.ea(p, T) + 1000);
        check("hf offset es", offmix.es(p, T), props.es(p, T));
        check("sf offset s", offmix.s(p, T), props.s(p, T) + 10);
    }

    // The Coolprop Function1 tracks the evaluator's saturation properties
    {
        IStringStream f1Stream
        (
            "sigma { type coolprop; fluid H2O; property sigma; }"
            "pSat  { type coolprop; fluid H2O; property pSat; }"
            "Tsat  { type coolprop; fluid H2O; property Tsat; }"
            "hl    { type coolprop; fluid H2O; property hl; }"
        );
        const dictionary f1Dict(f1Stream);

        check("Function1 sigma", function1Value(f1Dict, "sigma", T), props.sigma(p, T));
        check("Function1 pSat", function1Value(f1Dict, "pSat", T), props.pv(p, T));
        check("Function1 Tsat", function1Value(f1Dict, "Tsat", p), props.pvInvert(p));
        check("Function1 hl", function1Value(f1Dict, "hl", T), props.hl(p, T));

        // Tsat clamps to the ends of the saturation curve instead of
        // returning pvInvert's -1 failure sentinel below the triple
        // pressure or an unphysical value above the critical pressure
        check("Function1 Tsat(p < pt)", function1Value(f1Dict, "Tsat", 100), props.Tt());
        check("Function1 Tsat(p > Pc)", function1Value(f1Dict, "Tsat", 3e7), props.Tc());
    }

    // Runtime-selectable through the liquidProperties dictionary table
    // (e.g. from a Lagrangian liquid specification)
    {
        IStringStream liquidStream("H2O { type coolprop; fluid H2O; }");
        const dictionary liquidDict(liquidStream);

        autoPtr<liquidProperties> liquidPtr
        (
            liquidProperties::New(liquidDict.subDict("H2O"))
        );

        check("liquidProperties::New W", liquidPtr->W(), props.W());
        check("liquidProperties::New rho", liquidPtr->rho(p, T), props.rho(p, T));
    }

    if (nFail)
    {
        Info<< nl << nFail << " checks FAILED" << endl;
        return 1;
    }

    Info<< nl << "All checks passed" << nl << "End" << endl;
    return 0;
}


// ************************************************************************* //
