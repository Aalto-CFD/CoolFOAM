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

#include "coolprop.H"
#include "addToRunTimeSelectionTable.H"

#include "thermodynamicConstants.H"
using namespace Foam::constant::thermodynamic;

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace Function1s
{
    addScalarFunction1(coolprop);
}
}

const Foam::NamedEnum<Foam::Function1s::coolprop::propertyType, 4>
Foam::Function1s::coolprop::propertyNames
{
    "sigma",
    "pSat",
    "Tsat",
    "hl"
};


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::Function1s::coolprop::coolprop
(
    const word& name,
    const unitSets& units,
    const dictionary& dict
)
:
    FieldFunction1<scalar, coolprop>(name),
    property_(propertyNames[dict.lookup<word>("property")]),
    properties_(coolpropProperties::New(dict.lookup<word>("fluid")))
{
    // The values are evaluated in SI; refuse active unit conversions
    // rather than return silently mis-scaled results
    assertNoConvertUnits(typeName, units, dict);
}


Foam::Function1s::coolprop::coolprop(const coolprop& cp)
:
    FieldFunction1<scalar, coolprop>(cp),
    property_(cp.property_),
    properties_(new coolpropProperties(cp.properties_()))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::Function1s::coolprop::~coolprop()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::Function1s::coolprop::value(const scalar x) const
{
    // The saturation-curve properties are independent of the pressure
    // argument of the liquidProperties interface; pStd is passed for form
    switch (property_)
    {
        case propertyType::sigma:
            return properties_().sigma(pStd, x);
        case propertyType::pSat:
            return properties_().pv(pStd, x);
        case propertyType::Tsat:
        {
            // pvInvert returns -1 when the pressure is below the triple
            // pressure; clamp to the end of the saturation curve instead
            // of returning the sentinel as a temperature
            const scalar Tsat = properties_().pvInvert(x);
            const scalar Tt = properties_().Tt();

            if (Tsat < 0 && Tt < small)
            {
                FatalErrorInFunction
                    << "Saturation temperature of fluid "
                    << properties_().fluid() << " is not available at p = "
                    << x << " Pa, and the fluid provides no triple point "
                    << "to clamp to"
                    << exit(FatalError);
            }

            return max(Tsat, Tt);
        }
        case propertyType::hl:
            return properties_().hl(pStd, x);
    }

    return NaN;
}


Foam::scalar Foam::Function1s::coolprop::derivative(const scalar x) const
{
    // The relative step balances truncation against the CoolProp
    // evaluation noise, which is far larger than machine epsilon.
    // Accuracy degrades within ~1e-4*x of the ends of the saturation
    // curve (Tc, Pc, the triple point), where the stencil straddles the
    // clamped-constant branch of value().
    const scalar dx = 1e-4*max(mag(x), small);

    return (value(x + dx) - value(x - dx))/(2*dx);
}


Foam::scalar Foam::Function1s::coolprop::integral
(
    const scalar x1,
    const scalar x2
) const
{
    NotImplemented;
    return NaN;
}


void Foam::Function1s::coolprop::write
(
    Ostream& os,
    const unitSets& units
) const
{
    writeEntry(os, "fluid", properties_().fluid());
    writeEntry(os, "property", propertyNames[property_]);
}


// ************************************************************************* //
