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
    Helpers shared by the layers of the coolprop thermodynamics stack
    (Foam::Coolprop, Foam::CoolpropThermo and Foam::CoolpropTransport) for
    constructing their coolpropProperties evaluators from, and writing them
    back to, the model sub-dictionaries of the mixture specification.

\*---------------------------------------------------------------------------*/

#ifndef coolpropPropertiesModel_H
#define coolpropPropertiesModel_H

#include "coolpropProperties.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

//- Construct a coolpropProperties evaluator from a model sub-dictionary
//  (e.g. "equationOfState", "thermodynamics" or "transport/mu"),
//  validating the optional "type" entry
inline autoPtr<coolpropProperties> newCoolpropProperties
(
    const dictionary& dict
)
{
    if (dict.found("type"))
    {
        const word modelType = dict.lookup<word>("type");

        if (modelType != "coolprop")
        {
            FatalIOErrorInFunction(dict)
                << "Model for the CoolProp thermodynamics stack "
                << "should be of type coolprop. "
                << "Model type given is " << modelType << "."
                << exit(FatalIOError);
        }
    }

    return coolpropProperties::New
    (
        dict.lookup<word>("fluid"),
        coolpropProperties::phaseTypeNames.lookupOrDefault
        (
            "phase",
            dict,
            coolpropProperties::phaseType::liquid
        )
    );
}


//- Warn if a model's fluid/phase specification differs from the equation
//  of state's — a desynchronised sub-dictionary (e.g. the
//  $equationOfState macro not used after an edit) silently produces
//  hybrid physics otherwise
inline void warnCoolpropPropertiesMismatch
(
    const word& modelName,
    const coolpropProperties& properties,
    const coolpropProperties& eosProperties
)
{
    if
    (
        properties.fluid() != eosProperties.fluid()
     || properties.phase() != eosProperties.phase()
    )
    {
        WarningInFunction
            << "The " << modelName << " CoolProp specification (fluid "
            << properties.fluid() << ", phase "
            << coolpropProperties::phaseTypeNames[properties.phase()]
            << ") differs from the equationOfState specification (fluid "
            << eosProperties.fluid() << ", phase "
            << coolpropProperties::phaseTypeNames[eosProperties.phase()]
            << ")" << endl;
    }
}


//- Write the model entries that reconstruct the given coolpropProperties
//  evaluator, for use within a delimitDictionary scope
inline void writeCoolpropProperties
(
    Ostream& os,
    const coolpropProperties& properties
)
{
    writeEntry(os, "type", word("coolprop"));
    writeEntry(os, "fluid", properties.fluid());

    if (properties.phase() != coolpropProperties::phaseType::liquid)
    {
        writeEntry
        (
            os,
            "phase",
            coolpropProperties::phaseTypeNames[properties.phase()]
        );
    }
}

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
