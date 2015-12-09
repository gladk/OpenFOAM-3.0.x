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
      max
      (
        pressStatic()*mu0_.value()/
        max
         (
             strainRate(),
             (dimensionedScalar("VSMALL", dimless/dimTime, 0.00001))
         ),
        dimensionedScalar("VSMALL", dimMass/dimLength/dimTime, 0.00001)
      )/dimensionedScalar("VSMALL", dimDensity, 1.)
    )
  );
}


Foam::tmp<Foam::volScalarField>
Foam::viscosityModels::UTwente::pressStatic() const {
  return (h0_ - U_.mesh().C().component(vector::Z))*rho0_*(9.80665*dimensionedScalar("VSMALL", dimLength/dimTime/dimTime, 1.0));
}

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
    mu0_("mu0", dimless, UTwenteCoeffs_),
    nuMin_("nuMin", dimViscosity, UTwenteCoeffs_),
    nuMax_("nuMax", dimViscosity, UTwenteCoeffs_),
    h0_("h0", dimLength, UTwenteCoeffs_),
    rho0_("rho0", dimMass/dimLength/dimLength/dimLength, UTwenteCoeffs_),
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
        pressStatic()
    )
{};


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

bool Foam::viscosityModels::UTwente::read
(
    const dictionary& viscosityProperties
)
{
    viscosityModel::read(viscosityProperties);

    UTwenteCoeffs_ = viscosityProperties.subDict(typeName + "Coeffs");

    UTwenteCoeffs_.lookup("mu0") >> mu0_;
    UTwenteCoeffs_.lookup("nuMin") >> nuMin_;
    UTwenteCoeffs_.lookup("nuMax") >> nuMax_;
    UTwenteCoeffs_.lookup("h0") >> h0_;
    UTwenteCoeffs_.lookup("rho0") >> rho0_;

    return true;
}


// ************************************************************************* //
