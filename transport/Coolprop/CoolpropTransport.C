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

#include "CoolpropTransport.H"
#include "coolpropPropertiesModel.H"
#include "dictionary.H"
#include "delimitDictionary.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Thermo>
Foam::CoolpropTransport<Thermo>::CoolpropTransport
(
    const word& name,
    const dictionary& dict,
    const dictionary& muDict,
    const dictionary& kappaDict
)
:
    Thermo(name, dict),
    muProperties_(newCoolpropProperties(muDict).ptr()),
    kappaProperties_(newCoolpropProperties(kappaDict).ptr())
{
    warnCoolpropPropertiesMismatch
    (
        "transport/mu",
        muProperties_(),
        this->properties()
    );
    warnCoolpropPropertiesMismatch
    (
        "transport/kappa",
        kappaProperties_(),
        this->properties()
    );
}


template<class Thermo>
Foam::CoolpropTransport<Thermo>::CoolpropTransport
(
    const word& name,
    const dictionary& dict
)
:
    CoolpropTransport
    (
        name,
        dict,
        dict.subDict("transport").subDict("mu"),
        dict.subDict("transport").subDict("kappa")
    )
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Thermo>
void Foam::CoolpropTransport<Thermo>::write(Ostream& os) const
{
    Thermo::write(os);

    const delimitDictionary delimit(os, "transport");
    {
        const delimitDictionary delimitMu(os, "mu");
        writeCoolpropProperties(os, muProperties_());
    }
    {
        const delimitDictionary delimitKappa(os, "kappa");
        writeCoolpropProperties(os, kappaProperties_());
    }
}


// * * * * * * * * * * * * * * * Ostream Operator  * * * * * * * * * * * * * //

template<class Thermo>
Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const CoolpropTransport<Thermo>& cpt
)
{
    cpt.write(os);
    return os;
}


// ************************************************************************* //
