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
    vapour-companion properties (mug, kappag, Cpg).

Authors
    Stanislau Stasheuski, Aalto University, 2026.
    stanislau.stasheuski[at]aalto.fi

\*---------------------------------------------------------------------------*/

#include "coolpropProperties.H"

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

        // The liquid role must stay exactly incompressible, unaffected by
        // the gas role's psi/CpMCv overrides
        Info<< "    Liquid role: psi = " << coolPtr().psi(ps[0], Ts[0])
            << " CpMCv = " << coolPtr().CpMCv(ps[0], Ts[0])
            << " (both should be 0)" << endl;
    }

    Info<< nl << "End" << nl << endl;

    return 0;
}


// ************************************************************************* //
