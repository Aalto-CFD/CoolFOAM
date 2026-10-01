/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2026 Aalto University
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is originating from OpenFOAM but modified by authors described
    below.

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
    Test-coolpropProperties

Description
    Prints the coolpropProperties of the fluids given as arguments (H2O and
    Air if none are given) over a few liquid states and, when a built-in
    OpenFOAM liquid of the same name exists, the built-in properties for
    comparison.

    Also constructs each fluid a second time in the gas role ("phase gas;")
    and prints its properties over a p/T grid that spans both the
    metastable-fallback region (same points as the liquid role, deliberately
    subcooled from the gas role's point of view) and a genuinely superheated
    point, cross-checking against physical identities and the liquid role's
    vapour-companion properties (mug, kappag, Cpg). The liquid role's
    compressibility (psi, CpMCv, ea) is cross-checked the same way.

    Finally compares the tabulated SVDSBTL&HEOS::H2O with HEOS::H2O in both
    roles, returning non-zero if any property differs beyond its tolerance.

Authors
    Stanislau Stasheuski, Aalto University, 2026.
    stanislau.stasheuski[at]aalto.fi

\*---------------------------------------------------------------------------*/

#include "coolpropProperties.H"
#include "OSspecific.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

using namespace Foam;

void printConstants(const liquidProperties& l)
{
    Info<< "    W = " << l.W() << " kg/kmol"
        << ", Tc = " << l.Tc() << " K"
        << ", Pc = " << l.Pc() << " Pa"
        << ", Zc = " << l.Zc()
        << ", Tt = " << l.Tt() << " K"
        << ", Pt = " << l.Pt() << " Pa"
        << ", Tb = " << l.Tb() << " K"
        << ", omega = " << l.omega()
        << endl;
}


void printProperties(const liquidProperties& l, scalar p, scalar T)
{
    Info<< "    (p = " << p << " Pa, T = " << T << " K)"
        << " rho = " << l.rho(p, T)
        << " Cp = " << l.Cp(p, T)
        << " mu = " << l.mu(p, T)
        << " kappa = " << l.kappa(p, T)
        << " hs = " << l.hs(p, T)
        << " alphav = " << l.alphav(p, T)
        << " pv = " << l.pv(p, T)
        << " hl = " << l.hl(p, T)
        << " sigma = " << l.sigma(p, T)
        << " Cpg = " << l.Cpg(p, T)
        << " mug = " << l.mug(p, T)
        << " kappag = " << l.kappag(p, T)
        << endl;
}


// Prints the gas-role properties at (p, T), alongside the physical/cross
// checks against the liquid-role instance of the same fluid: mu/kappa
// should match the liquid role's mug/kappag exactly (both ultimately query
// the same vapour state), and ha - ea should match p/rho (h = u + p/rho,
// referenced to the same CoolProp state so offsets cancel)
void printGasProperties
(
    const coolpropProperties& g,
    const coolpropProperties& l,
    scalar p,
    scalar T
)
{
    Info<< "    (p = " << p << " Pa, T = " << T << " K)"
        << " rho = " << g.rho(p, T)
        << " Cp = " << g.Cp(p, T)
        << " Cpg(ideal) = " << g.Cpg(p, T)
        << " mu = " << g.mu(p, T) << " (mug = " << l.mug(p, T) << ")"
        << " kappa = " << g.kappa(p, T) << " (kappag = " << l.kappag(p, T)
        << ")"
        << " ha = " << g.ha(p, T)
        << " ea = " << g.ea(p, T)
        << " ha-ea = " << g.ha(p, T) - g.ea(p, T)
        << " (p/rho = " << p/g.rho(p, T) << ")"
        << " psi = " << g.psi(p, T)
        << " CpMCv = " << g.CpMCv(p, T)
        << endl;
}


// Prints psi, CpMCv and ha - ea at (p, T) alongside their physical
// cross-checks: psi against a central difference of rho in p,
// CpMCv against T*alphav^2/psi (Cp - Cv = T*alphav^2/(rho*kappaT)) and
// ha - ea against p/rho. The identities hold where the direct flash
// succeeds, not in the saturated fallback (e.g. liquid Air at 300 K)
void printCompressibility
(
    const coolpropProperties& l,
    scalar p,
    scalar T
)
{
    const scalar dp = 1e-4*p;
    const scalar psi = l.psi(p, T);

    Info<< "    (p = " << p << " Pa, T = " << T << " K)"
        << " psi = " << psi
        << " (drho/dp = "
        << (l.rho(p + dp, T) - l.rho(p - dp, T))/(2*dp) << ")"
        << " CpMCv = " << l.CpMCv(p, T)
        << " (T*alphav^2/psi = " << T*sqr(l.alphav(p, T))/psi << ")"
        << " ha-ea = " << l.ha(p, T) - l.ea(p, T)
        << " (p/rho = " << p/l.rho(p, T) << ")"
        << endl;
}


// Compares the given property of a tabulated (SVDSBTL) fluid with its
// reference backend, counting the differences relative to the reference,
// or to the given scale where the reference is smaller, beyond tol
label compare
(
    const word& what,
    const scalar tabulated,
    const scalar reference,
    const scalar tol,
    const scalar scale = rootVSmall
)
{
    const scalar err =
        mag(tabulated - reference)/max(mag(reference), scale);

    Info<< "    " << (err > tol ? "FAIL " : "pass ") << what
        << ": " << tabulated << " (reference " << reference
        << ", relative difference " << err << ")" << endl;

    return err > tol;
}


// Compares a tabulated (SVDSBTL) fluid in the given role with its
// reference backend: the tabulated backend does not impose the phase, so
// this also checks that each role stays on its side of the saturation
// curve, and it provides no derivatives, so this checks their
// finite-difference approximations
label compareTabulated
(
    const word& tabulatedFluid,
    const word& referenceFluid,
    const coolpropProperties::phaseType phase,
    const scalarList& ps,
    const scalarList& Ts
)
{
    autoPtr<coolpropProperties> tabPtr
    (
        coolpropProperties::New(tabulatedFluid, phase)
    );
    autoPtr<coolpropProperties> refPtr
    (
        coolpropProperties::New(referenceFluid, phase)
    );
    const coolpropProperties& tab = tabPtr();
    const coolpropProperties& ref = refPtr();

    Info<< nl << "CoolProp " << tabulatedFluid << " vs " << referenceFluid
        << " (" << coolpropProperties::phaseTypeNames[phase] << " role):"
        << endl;

    label nFail = 0;

    // The source backend correlations (e.g. REFPROP's) of the properties
    // the table does not provide may differ slightly from the reference's
    nFail += compare("omega", tab.omega(), ref.omega(), 1e-4);

    forAll(ps, i)
    {
        const scalar p = ps[i], T = Ts[i];

        Info<< "  (p = " << p << " Pa, T = " << T << " K)" << endl;

        nFail += compare("rho", tab.rho(p, T), ref.rho(p, T), 1e-2);
        nFail += compare("hs", tab.hs(p, T), ref.hs(p, T), 1e-2);
        nFail += compare("es", tab.es(p, T), ref.es(p, T), 1e-2);
        nFail += compare("Cp", tab.Cp(p, T), ref.Cp(p, T), 1e-2);
        nFail += compare
        (
            "alphav", tab.alphav(p, T), ref.alphav(p, T), 2e-2, 1e-4
        );
        nFail += compare("psi", tab.psi(p, T), ref.psi(p, T), 2e-2);
        nFail += compare("CpMCv", tab.CpMCv(p, T), ref.CpMCv(p, T), 2e-2, 10);
        nFail += compare("mu", tab.mu(p, T), ref.mu(p, T), 1e-2);
        nFail += compare("kappa", tab.kappa(p, T), ref.kappa(p, T), 1e-2);
        nFail += compare("pv", tab.pv(p, T), ref.pv(p, T), 1e-3);
        nFail += compare("sigma", tab.sigma(p, T), ref.sigma(p, T), 1e-2);
        nFail += compare("Cpg", tab.Cpg(p, T), ref.Cpg(p, T), 1e-6);
        nFail += compare("hl", tab.hl(p, T), ref.hl(p, T), 1e-2);

        // The density must respond to pressure as psi says, also where
        // extrapolated beyond saturation
        const scalar dp = 1e-3*p;
        nFail += compare
        (
            "drho/dp",
            (tab.rho(p + dp, T) - tab.rho(p - dp, T))/(2*dp),
            tab.psi(p, T),
            1e-2
        );
    }

    // Beyond the table (below the triple point pressure or temperature, as
    // reached transiently by a solver), where the table returns NaN states:
    // the properties must be those of the saturated fallback, finite and on
    // the side of the role
    const scalarList psOut({100, 1e5, -1e5});
    const scalarList TsOut({300, 260, 300});
    const scalar rhoc = tab.W()/tab.Vc();

    forAll(psOut, i)
    {
        const scalar p = psOut[i], T = TsOut[i];
        const scalar rho = tab.rho(p, T);

        const bool pass =
            std::isfinite(rho)
         && std::isfinite(tab.Cp(p, T))
         && std::isfinite(tab.psi(p, T))
         && std::isfinite(tab.mu(p, T))
         && (phase == coolpropProperties::phaseType::liquid)
         == (rho > rhoc);

        Info<< "  " << (pass ? "pass" : "FAIL") << " beyond the table (p = "
            << p << " Pa, T = " << T << " K): rho = " << rho
            << " Cp = " << tab.Cp(p, T) << " psi = " << tab.psi(p, T)
            << " mu = " << tab.mu(p, T) << endl;

        nFail += !pass;
    }

    // Liquid role next to the critical point at a low pressure, as seen by
    // the liquid thermo of an overheated vapour cell: the extrapolation
    // from psat ~ Pc must be held rather than let the energy run away
    // (the internal energy of the saturated liquid at Tc is ~ 2e6 J/kg)
    if (phase == coolpropProperties::phaseType::liquid)
    {
        const scalar p = 1e3, T = 0.9999*tab.Tc();
        const scalar ea = tab.ea(p, T), rho = tab.rho(p, T);

        const bool pass =
            std::isfinite(ea) && ea > 0 && ea < 2.5e6
         && std::isfinite(rho) && rho > rhoc;

        Info<< "  " << (pass ? "pass" : "FAIL")
            << " near-critical liquid held (p = " << p << " Pa, T = " << T
            << " K): ea = " << ea << " rho = " << rho << endl;

        nFail += !pass;
    }

    // Gas role in the liquid cells of a VoF case, compressed far beyond
    // saturation: the energy extrapolation must be held rather than run
    // away with the pressure, while the density keeps responding to it
    if (phase == coolpropProperties::phaseType::gas && 300 < tab.Tc())
    {
        const scalar T = 300;
        const scalar es1 = tab.es(1e6, T), es2 = tab.es(3e7, T);
        const scalar rho1 = tab.rho(1e6, T), rho2 = tab.rho(3e7, T);

        const bool pass =
            es1 == es2 && std::isfinite(rho2) && rho2 > rho1;

        Info<< "  " << (pass ? "pass" : "FAIL")
            << " compressed vapour energy held (T = " << T
            << " K): es(1e6) = " << es1 << " es(3e7) = " << es2
            << " rho(1e6) = " << rho1 << " rho(3e7) = " << rho2 << endl;

        nFail += !pass;
    }

    return nFail;
}


int main(int argc, char *argv[])
{
    wordList fluids;

    if (argc > 1)
    {
        for (int argi = 1; argi < argc; argi++)
        {
            fluids.append(argv[argi]);
        }
    }
    else
    {
        fluids.append("H2O");
        fluids.append("Air");
    }

    const scalarList ps({1e5, 1e6});
    const scalarList Ts({300, 400});

    // Gas-role grid: the two points above (deliberately subcooled from the
    // gas role's point of view, exercising the metastable Q = 1 fallback)
    // plus a genuinely superheated point (direct PT flash, no fallback)
    const scalarList gasPs({1e5, 1e6, 1e5});
    const scalarList gasTs({300, 400, 500});

    forAll(fluids, fluidi)
    {
        const word& fluid = fluids[fluidi];

        autoPtr<coolpropProperties> coolPtr(coolpropProperties::New(fluid));

        Info<< nl << "CoolProp " << fluid << ":" << endl;
        printConstants(coolPtr());

        forAll(ps, i)
        {
            printProperties(coolPtr(), ps[i], Ts[i]);
        }

        // Compare with the built-in liquid of the same name, if available
        if (liquidProperties::ConstructorTablePtr_->found(fluid))
        {
            autoPtr<liquidProperties> builtInPtr(liquidProperties::New(fluid));

            Info<< nl << "Built-in " << fluid << ":" << endl;
            printConstants(builtInPtr());

            forAll(ps, i)
            {
                printProperties(builtInPtr(), ps[i], Ts[i]);
            }
        }

        // Gas role
        autoPtr<coolpropProperties> gasPtr
        (
            coolpropProperties::New(fluid, coolpropProperties::phaseType::gas)
        );

        Info<< nl << "CoolProp " << fluid << " (gas role):" << endl;
        printConstants(gasPtr());

        forAll(gasPs, i)
        {
            printGasProperties(gasPtr(), coolPtr(), gasPs[i], gasTs[i]);
        }

        // The liquid role is compressible too
        Info<< nl << "CoolProp " << fluid << " (liquid role compressibility):"
            << endl;

        forAll(ps, i)
        {
            printCompressibility(coolPtr(), ps[i], Ts[i]);
        }
    }

    // Tabulated backend: the first run builds the tables into
    // ~/.CoolProp/SVDTables (or ALTERNATIVE_SVDTABLES_DIRECTORY), which
    // takes about a minute per input pair
    label nFail = 0;

    // The REFPROP-sourced table only if REFPROP is available
    wordList tabulatedFluids({"SVDSBTL&HEOS::H2O"});

    if (!getEnv("COOLPROP_REFPROP_ROOT").empty())
    {
        tabulatedFluids.append("SVDSBTL&REFPROP::H2O");
    }

    forAll(tabulatedFluids, i)
    {
        // Liquid role: subcooled, compressed, flashing (p < psat(T), where
        // the table extrapolates from off the saturation curve and HEOS
        // evaluates the metastable liquid), and next to the triple point
        // (where the REFPROP-sourced table is wrong close to saturation)
        nFail += compareTabulated
        (
            tabulatedFluids[i],
            "HEOS::H2O",
            coolpropProperties::phaseType::liquid,
            scalarList({1e5, 1e6, 2e6, 1e4, 1e3, 650, 1e3, 3e7}),
            scalarList({300, 400, 450, 350, 300, 273.3, 274, 273.2})
        );

        // Gas role: superheated, and slightly subcooled (p > psat(T), where
        // the table extrapolates and HEOS evaluates the metastable vapour)
        nFail += compareTabulated
        (
            tabulatedFluids[i],
            "HEOS::H2O",
            coolpropProperties::phaseType::gas,
            scalarList({1e5, 1e4, 1e6, 3.6e3, 640}),
            scalarList({400, 350, 500, 300, 274})
        );
    }

    // Supercritical gas (air): no saturation curve, but a table lower
    // pressure limit, the triple pressure (5.3 kPa), below which the table
    // is extrapolated
    nFail += compareTabulated
    (
        "SVDSBTL&HEOS::Air",
        "HEOS::Air",
        coolpropProperties::phaseType::gas,
        scalarList({1e5, 1e6, 1e4, 5e6, 2e3}),
        scalarList({300, 578, 250, 300, 300})
    );

    Info<< nl << (nFail ? "Some" : "All") << " tabulated checks "
        << (nFail ? "failed" : "passed") << nl;

    Info<< nl << "End" << nl << endl;

    return nFail ? 1 : 0;
}


// ************************************************************************* //
