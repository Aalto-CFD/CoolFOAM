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

Authors
    Henry Weller, CFD-Direct, 2017-2023.

    Stanislau Stasheuski, Aalto University, 2026.
    stanislau.stasheuski[at]aalto.fi

Description
    Instantiation of the fluid thermodynamics packages for the Coolprop
    properties, selected by

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
    \endverbatim

    Mirrors the instantiations of the built-in NSRDS liquid properties in
    $FOAM_SRC/thermophysicalModels/basic/rhoFluidThermo/rhoFluidThermos.C.

\*---------------------------------------------------------------------------*/

#include "rhoFluidThermo.H"

#include "pureMixture.H"

#include "forCoolprop.H"

#include "makeFluidThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    forCoolprop(makeFluidThermo, rhoFluidThermo, pureMixture);
}

// ************************************************************************* //
