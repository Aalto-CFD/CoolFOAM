/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2026 Aalto University
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is originating from OpenFOAM but modified by authors described
    in the according header file.

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

#include "coolPropPropertiesSelector.H"

// * * * * * * * * * * * * * * * Local Functions  * * * * * * * * * * * * * //

namespace Foam
{

//- Normalise a flat mixture dict ("fluid; key val; ...") into the nested
//  sub-dict form ("fluid { key val; ... }") expected by
//  thermophysicalPropertiesSelector, which only inspects the first entry.
//  Left unchanged if already nested, or if there is nothing to nest
//  (the plain "fluid;" form), so the base selector takes the same word-New
//  path it always has in those cases.
static dictionary normaliseMixtureDict(const dictionary& dict)
{
    const word fluidName(dict.first()->keyword());

    if (dict.isDict(fluidName) || dict.size() == 1)
    {
        return dict;
    }

    // "phase" (the only option key so far) placed before the fluid entry
    // would otherwise be mistaken for the fluid name below
    if (fluidName == "phase")
    {
        FatalIOErrorInFunction(dict)
            << "Expected the fluid specification as the first entry of "
            << "the mixture dict, found the option entry " << fluidName
            << " instead" << nl
            << "Move the fluid specification (e.g. \"H2O;\") to the start "
            << "of the mixture dict"
            << exit(FatalIOError);
    }

    dictionary fluidDict(dict);
    fluidDict.remove(fluidName);

    // Keep the original dict's name so a FatalIOError while parsing the
    // (e.g. "phase") entries below still points back at the case file
    dictionary wrapperDict(dict.name());
    wrapperDict.add(fluidName, fluidDict);

    return wrapperDict;
}

} // End namespace Foam


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::coolPropPropertiesSelector::coolPropPropertiesSelector
(
    const word& name,
    const dictionary& dict
)
:
    thermophysicalPropertiesSelector<coolPropProperties>
    (
        name,
        normaliseMixtureDict(dict)
    )
{}


// ************************************************************************* //
