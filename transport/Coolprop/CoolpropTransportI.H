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

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Thermo>
inline Foam::CoolpropTransport<Thermo>::CoolpropTransport
(
    const CoolpropTransport& cpt
)
:
    Thermo(cpt),
    muProperties_(new coolpropProperties(cpt.muProperties_())),
    kappaProperties_(new coolpropProperties(cpt.kappaProperties_()))
{}


template<class Thermo>
inline Foam::CoolpropTransport<Thermo>::CoolpropTransport
(
    const word& name,
    const CoolpropTransport& cpt
)
:
    Thermo(name, cpt),
    muProperties_(new coolpropProperties(cpt.muProperties_())),
    kappaProperties_(new coolpropProperties(cpt.kappaProperties_()))
{}


template<class Thermo>
inline Foam::autoPtr<Foam::CoolpropTransport<Thermo>>
Foam::CoolpropTransport<Thermo>::clone() const
{
    return autoPtr<CoolpropTransport<Thermo>>
    (
        new CoolpropTransport<Thermo>(*this)
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Thermo>
inline Foam::scalar Foam::CoolpropTransport<Thermo>::mu
(
    const scalar p,
    const scalar T
) const
{
    return muProperties_().mu(p, T);
}


template<class Thermo>
inline Foam::scalar Foam::CoolpropTransport<Thermo>::kappa
(
    const scalar p,
    const scalar T
) const
{
    return kappaProperties_().kappa(p, T);
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Thermo>
inline void Foam::CoolpropTransport<Thermo>::operator=
(
    const CoolpropTransport& cpt
)
{
    Thermo::operator=(cpt);
    muProperties_.reset(new coolpropProperties(cpt.muProperties_()));
    kappaProperties_.reset(new coolpropProperties(cpt.kappaProperties_()));
}


// ************************************************************************* //
