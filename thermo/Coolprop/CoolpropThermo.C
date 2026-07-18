/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2026 Aalto University
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is originating from OpenFOAM but modified by the authors
    described in the according header file.

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

#include "CoolpropThermo.H"
#include "coolpropPropertiesModel.H"
#include "dictionary.H"
#include "delimitDictionary.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class EquationOfState>
Foam::CoolpropThermo<EquationOfState>::CoolpropThermo
(
    const word& name,
    const dictionary& dict,
    const dictionary& subDict
)
:
    EquationOfState(name, dict),
    properties_(newCoolpropProperties(subDict).ptr()),
    hf_
    (
        subDict.lookupOrDefaultBackwardsCompatible<scalar>
        (
            {"hf", "Hf"},
            dimensions::specificEnergy,
            0
        )
    ),
    sf_
    (
        subDict.lookupOrDefaultBackwardsCompatible<scalar>
        (
            {"sf", "Sf"},
            dimensions::specificEntropy,
            0
        )
    )
{
    warnCoolpropPropertiesMismatch
    (
        "thermodynamics",
        properties_(),
        EquationOfState::properties()
    );
}


template<class EquationOfState>
Foam::CoolpropThermo<EquationOfState>::CoolpropThermo
(
    const word& name,
    const dictionary& dict
)
:
    CoolpropThermo(name, dict, dict.subDict("thermodynamics"))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class EquationOfState>
void Foam::CoolpropThermo<EquationOfState>::write(Ostream& os) const
{
    EquationOfState::write(os);

    const delimitDictionary delimit(os, "thermodynamics");
    writeCoolpropProperties(os, properties_());
    writeEntry(os, "hf", hf_);
    writeEntry(os, "sf", sf_);
}


// * * * * * * * * * * * * * * * Ostream Operator  * * * * * * * * * * * * * //

template<class EquationOfState>
Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const CoolpropThermo<EquationOfState>& cpt
)
{
    cpt.write(os);
    return os;
}


// ************************************************************************* //
