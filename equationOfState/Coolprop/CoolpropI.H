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

// * * * * * * * * * * * * * Private Constructors  * * * * * * * * * * * * * //

template<class Specie>
inline Foam::Coolprop<Specie>::Coolprop
(
    const word& name,
    const dictionary& dict,
    autoPtr<coolpropProperties> properties
)
:
    Specie
    (
        dict.isDict("specie")
      ? Specie(name, dict)
      : Specie(name, 1, properties->W())
    ),
    properties_(properties.ptr())
{}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Specie>
inline Foam::Coolprop<Specie>::Coolprop(const Coolprop& cp)
:
    Specie(cp),
    properties_(new coolpropProperties(cp.properties_()))
{}


template<class Specie>
inline Foam::Coolprop<Specie>::Coolprop
(
    const word& name,
    const Coolprop& cp
)
:
    Specie(name, cp),
    properties_(new coolpropProperties(cp.properties_()))
{}


template<class Specie>
inline Foam::autoPtr<Foam::Coolprop<Specie>>
Foam::Coolprop<Specie>::clone() const
{
    return autoPtr<Coolprop<Specie>>(new Coolprop<Specie>(*this));
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Specie>
inline const Foam::coolpropProperties&
Foam::Coolprop<Specie>::properties() const
{
    return properties_();
}


template<class Specie>
inline Foam::scalar Foam::Coolprop<Specie>::rho
(
    const scalar p,
    const scalar T
) const
{
    return properties_().rho(p, T);
}


template<class Specie>
inline Foam::scalar Foam::Coolprop<Specie>::h
(
    const scalar p,
    const scalar T
) const
{
    return 0;
}


template<class Specie>
inline Foam::scalar Foam::Coolprop<Specie>::Cp
(
    const scalar p,
    const scalar T
) const
{
    return 0;
}


template<class Specie>
inline Foam::scalar Foam::Coolprop<Specie>::e
(
    const scalar p,
    const scalar T
) const
{
    return 0;
}


template<class Specie>
inline Foam::scalar Foam::Coolprop<Specie>::Cv
(
    const scalar p,
    const scalar T
) const
{
    return 0;
}


template<class Specie>
inline Foam::scalar Foam::Coolprop<Specie>::sp
(
    const scalar p,
    const scalar T
) const
{
    return 0;
}


template<class Specie>
inline Foam::scalar Foam::Coolprop<Specie>::sv
(
    const scalar p,
    const scalar T
) const
{
    return 0;
}


template<class Specie>
inline Foam::scalar Foam::Coolprop<Specie>::psi
(
    const scalar p,
    const scalar T
) const
{
    return properties_().psi(p, T);
}


template<class Specie>
inline Foam::scalar Foam::Coolprop<Specie>::Z
(
    const scalar p,
    const scalar T
) const
{
    return p/(properties_().rho(p, T)*this->R()*T);
}


template<class Specie>
inline Foam::scalar Foam::Coolprop<Specie>::CpMCv
(
    const scalar p,
    const scalar T
) const
{
    return properties_().CpMCv(p, T);
}


template<class Specie>
inline Foam::scalar Foam::Coolprop<Specie>::alphav
(
    const scalar p,
    const scalar T
) const
{
    return properties_().alphav(p, T);
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Specie>
inline void Foam::Coolprop<Specie>::operator=(const Coolprop& cp)
{
    Specie::operator=(cp);
    properties_.reset(new coolpropProperties(cp.properties_()));
}


// ************************************************************************* //
