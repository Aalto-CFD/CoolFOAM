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

#include "thermodynamicConstants.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class EquationOfState>
inline Foam::CoolpropThermo<EquationOfState>::CoolpropThermo
(
    const CoolpropThermo& cpt
)
:
    EquationOfState(cpt),
    properties_(new coolpropProperties(cpt.properties_())),
    hf_(cpt.hf_),
    sf_(cpt.sf_)
{}


template<class EquationOfState>
inline Foam::CoolpropThermo<EquationOfState>::CoolpropThermo
(
    const word& name,
    const CoolpropThermo& cpt
)
:
    EquationOfState(name, cpt),
    properties_(new coolpropProperties(cpt.properties_())),
    hf_(cpt.hf_),
    sf_(cpt.sf_)
{}


template<class EquationOfState>
inline Foam::autoPtr<Foam::CoolpropThermo<EquationOfState>>
Foam::CoolpropThermo<EquationOfState>::clone() const
{
    return autoPtr<CoolpropThermo<EquationOfState>>
    (
        new CoolpropThermo<EquationOfState>(*this)
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class EquationOfState>
inline Foam::scalar Foam::CoolpropThermo<EquationOfState>::limit
(
    const scalar T
) const
{
    // Beyond the triple point (and, for the liquid role, the critical
    // point) the backend falls back to the same clamped saturated state,
    // so the property curves plateau and a Newton inversion iterating out
    // there cannot converge. Match the backend's own clamp instead.
    const coolpropProperties& properties = properties_();

    // Tt is zero if CoolProp does not provide a triple point for the
    // fluid; keep the floor positive so iterates cannot reach exactly 0 K
    const scalar Tmin = max(properties.Tt(), small);

    if (properties.phase() == coolpropProperties::phaseType::liquid)
    {
        return min(max(T, Tmin), 0.9999*properties.Tc());
    }

    return max(T, Tmin);
}


template<class EquationOfState>
inline Foam::scalar Foam::CoolpropThermo<EquationOfState>::Cp
(
    const scalar p,
    const scalar T
) const
{
    return properties_().Cp(p, T);
}


template<class EquationOfState>
inline Foam::scalar Foam::CoolpropThermo<EquationOfState>::ha
(
    const scalar p,
    const scalar T
) const
{
    return properties_().ha(p, T) + hf_;
}


template<class EquationOfState>
inline Foam::scalar Foam::CoolpropThermo<EquationOfState>::hs
(
    const scalar p,
    const scalar T
) const
{
    return properties_().hs(p, T);
}


template<class EquationOfState>
inline Foam::scalar Foam::CoolpropThermo<EquationOfState>::hf() const
{
    return properties_().hf() + hf_;
}


template<class EquationOfState>
inline Foam::scalar Foam::CoolpropThermo<EquationOfState>::s
(
    const scalar p,
    const scalar T
) const
{
    return properties_().s(p, T) + sf_;
}


template<class EquationOfState>
inline Foam::scalar Foam::CoolpropThermo<EquationOfState>::gStd
(
    const scalar T
) const
{
    using namespace constant::thermodynamic;

    // For the liquid role above the boiling point at pStd the backend
    // falls back to the saturated state at T, so this is then the Gibbs
    // free energy on the saturation curve rather than at pStd — the
    // natural liquid analogue of the pressure-independent standard-state
    // correlations of the built-in thermodynamics packages
    return
        properties_().ha(pStd, T) + hf_
      - T*(properties_().s(pStd, T) + sf_);
}


template<class EquationOfState>
inline Foam::scalar Foam::CoolpropThermo<EquationOfState>::dCpdT
(
    const scalar p,
    const scalar T
) const
{
    // The relative step balances truncation against the CoolProp
    // evaluation noise (flash tolerances, table interpolation), which is
    // far larger than machine epsilon — a machine-epsilon-based step
    // would return mostly cancellation noise. Floored to protect the
    // division at (unphysical) zero temperature.
    const scalar dT = max(1e-4*T, rootSmall);

    return
        (properties_().Cp(p, T + dT) - properties_().Cp(p, T - dT))/(2*dT);
}


template<class EquationOfState>
inline Foam::scalar Foam::CoolpropThermo<EquationOfState>::Cv
(
    const scalar p,
    const scalar T
) const
{
    return properties_().Cp(p, T) - properties_().CpMCv(p, T);
}


template<class EquationOfState>
inline Foam::scalar Foam::CoolpropThermo<EquationOfState>::es
(
    const scalar p,
    const scalar T
) const
{
    return properties_().es(p, T);
}


template<class EquationOfState>
inline Foam::scalar Foam::CoolpropThermo<EquationOfState>::ea
(
    const scalar p,
    const scalar T
) const
{
    return properties_().ea(p, T) + hf_;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class EquationOfState>
inline void Foam::CoolpropThermo<EquationOfState>::operator=
(
    const CoolpropThermo& cpt
)
{
    EquationOfState::operator=(cpt);
    properties_.reset(new coolpropProperties(cpt.properties_()));
    hf_ = cpt.hf_;
    sf_ = cpt.sf_;
}


// ************************************************************************* //
