/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2026 Aalto University
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is originating from OpenFOAM and modified by the authors
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
    Henry Weller, CFD-Direct, 2017-2026.

    Stanislau Stasheuski, Aalto University, 2026.
    stanislau.stasheuski[at]aalto.fi

Description
    Wrapper class providing run-time selection of coolPropProperties
    for the templated thermodynamics packages.

    Mirrors Foam::liquidPropertiesSelector.

SourceFiles
    coolPropPropertiesSelector.C
    coolPropPropertiesSelectorI.H

\*---------------------------------------------------------------------------*/

#ifndef coolPropPropertiesSelector_H
#define coolPropPropertiesSelector_H

#include "thermophysicalPropertiesSelector.H"
#include "coolPropProperties.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                 Class coolPropPropertiesSelector Declaration
\*---------------------------------------------------------------------------*/

class coolPropPropertiesSelector
:
    public thermophysicalPropertiesSelector<coolPropProperties>
{
public:

    // Constructors

        //- Inherit constructors
        using
            thermophysicalPropertiesSelector<coolPropProperties>::
            thermophysicalPropertiesSelector;

        //- Construct from name and dictionary, normalising the flat
        //  mixture form (e.g. "H2O; phase gas;") into the nested sub-dict
        //  form ("H2O { phase gas; }") expected by the base selector.
        //  Shadows (rather than being brought in by) the inherited
        //  constructor of the same signature.
        coolPropPropertiesSelector(const word& name, const dictionary& dict);


    //- Destructor
    virtual ~coolPropPropertiesSelector()
    {}


    // Member Functions

        // Physical properties

            //- Surface tension [N/m]
            virtual scalar sigma(const scalar p, const scalar T) const;

            //- Viscosity [Pa s]
            inline scalar mu(const scalar p, const scalar T) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "coolPropPropertiesSelectorI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
