/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2015 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

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

#include "UTwente.H"
#include "addToRunTimeSelectionTable.H"
#include "surfaceFields.H"
#include "volFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace viscosityModels
{
    defineTypeNameAndDebug(UTwente, 0);

    addToRunTimeSelectionTable
    (
        viscosityModel,
        UTwente,
        dictionary
    );
}
}


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::viscosityModels::UTwente::calcNu() const
{
    return max
    (
        nuMin_,
        min
        (
            nuMax_,
            k_*pow
            (
                max
                (
                    dimensionedScalar("one", dimTime, 1.0)*strainRate(),
                    dimensionedScalar("VSMALL", dimless, VSMALL)
                ),
                mu0_.value() - scalar(1.0)
            )
        )
    );
}
/*
Foam::tmp<Foam::volScalarField>
Foam::viscosityModels::UTwente::pressStatic() const
{
  
}
*/


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::viscosityModels::UTwente::UTwente
(
    const word& name,
    const dictionary& viscosityProperties,
    const volVectorField& U,
    const surfaceScalarField& phi
)
:
    viscosityModel(name, viscosityProperties, U, phi),
    UTwenteCoeffs_(viscosityProperties.subDict(typeName + "Coeffs")),
    k_("k", dimViscosity, UTwenteCoeffs_),
    mu0_("mu0", dimless, UTwenteCoeffs_),
    nuMin_("nuMin", dimViscosity, UTwenteCoeffs_),
    nuMax_("nuMax", dimViscosity, UTwenteCoeffs_),
    nu_
    (
        IOobject
        (
            name,
            U_.time().timeName(),
            U_.db(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        calcNu()
    ),
    strainRate_
    (
        IOobject
        (
            "strainRate",
            U_.time().timeName(),
            U_.db(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        strainRate()
    ),
    pressStatic_
    (
        IOobject
        (
            "pressStatic",
            U_.time().timeName(),
            U_.db(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        strainRate()
    )
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

bool Foam::viscosityModels::UTwente::read
(
    const dictionary& viscosityProperties
)
{
    viscosityModel::read(viscosityProperties);

    UTwenteCoeffs_ = viscosityProperties.subDict(typeName + "Coeffs");

    UTwenteCoeffs_.lookup("k") >> k_;
    UTwenteCoeffs_.lookup("mu0") >> mu0_;
    UTwenteCoeffs_.lookup("nuMin") >> nuMin_;
    UTwenteCoeffs_.lookup("nuMax") >> nuMax_;

    return true;
}


// ************************************************************************* //
