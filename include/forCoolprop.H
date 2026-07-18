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

Authors
    Stanislau Stasheuski, Aalto University, 2026.
    stanislau.stasheuski[at]aalto.fi

Description
    Instantiation macros for the CoolProp thermodynamics stack, mirroring
    the forNSRDSLiquids macros of $FOAM_SRC/thermophysicalModels/specie/
    include/forLiquids.H.

\*---------------------------------------------------------------------------*/

#ifndef forCoolprop_H
#define forCoolprop_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "specie.H"

#include "Coolprop.H"

#include "CoolpropThermo.H"

#include "sensibleEnthalpy.H"
#include "sensibleInternalEnergy.H"

#include "CoolpropTransport.H"

#include "thermo.H"

#include "forThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define forCoolpropEqns(Mu, He, Cp, Macro, Args...)                            \
    forThermo(Mu, He, Cp, Coolprop, specie, Macro, Args)

#define forCoolpropThermos(Mu, He, Macro, Args...)                             \
    forCoolpropEqns(Mu, He, CoolpropThermo, Macro, Args)

#define forCoolpropEnergies(Mu, Macro, Args...)                                \
    forCoolpropThermos(Mu, sensibleEnthalpy, Macro, Args);                     \
    forCoolpropThermos(Mu, sensibleInternalEnergy, Macro, Args)

#define forCoolprop(Macro, Args...)                                            \
    forCoolpropEnergies(CoolpropTransport, Macro, Args)

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
