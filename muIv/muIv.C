/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
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

Author
    Alexander Jarosch research@alexj.at

\*---------------------------------------------------------------------------*/

#include "muIv.H"
#include "addToRunTimeSelectionTable.H"
#include "surfaceFields.H"
#include "fvcGrad.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace viscosityModels
{
    defineTypeNameAndDebug(muIv, 0);

    addToRunTimeSelectionTable
    (
        viscosityModel,
        muIv,
        dictionary
    );
}
}


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::viscosityModels::muIv::calcNu() const
{
    const objectRegistry& db = U_.db();
    const volScalarField& alphag = U_.mesh().lookupObject<volScalarField>("alpha.granul");

    if (db.foundObject<volScalarField>("p")) {
        Info<< "Calculate mu(Iv) based on pressure" << endl;
        return
        (
            max(
                min(
                    (mu_*peff_)/(2.0*rhog_*normD_*max(alphag, alphaSmall_)), nuMax_
                ), nuMin_
            )
        );
    } else{
        Info<< "Pressure not found for mu(Iv), return zero" << endl;
        return  tmp<volScalarField>
        (
            new volScalarField
            (
                IOobject
                (
                    "nuis0",
                    U_.time().timeName(),
                    U_.db(),
                    IOobject::NO_READ,
                    IOobject::NO_WRITE,
                    false
                ),
               U_.mesh(),
               dimensionedScalar("nuis0", dimViscosity, 0.0)
           )
        );
    }

}

Foam::tmp<Foam::volScalarField>
Foam::viscosityModels::muIv::calcMu() const
{
    const objectRegistry& db = U_.db();
    if (db.foundObject<volScalarField>("p")) {
        return
        (
           mus_ +  (mud_ - mus_)/
            (1 + I0_/Iv_)
        );
    } else {
        Info<< "Pressure not found for mu(Iv), return mus" << endl;
        return mus_*calcIv();
    }
}

Foam::tmp<Foam::volScalarField>
Foam::viscosityModels::muIv::calcIv() const
{
    const objectRegistry& db = U_.db();
    if (db.foundObject<volScalarField>("p")) {
        // Info<< "Calculate I based on pressure" << endl;
        return
        (
            normD_*nuWater_*rhoWater_/peff_
        );
    } else {
        Info<< "Pressure not found for I, return zero" << endl;
        return  tmp<volScalarField>
        (
            new volScalarField
            (
                IOobject
                (
                    "II0",
                    U_.time().timeName(),
                    U_.db(),
                    IOobject::NO_READ,
                    IOobject::NO_WRITE,
                    false
                ),
               U_.mesh(),
               dimensionedScalar("II0", dimless, 0.0)
           )
        );
    }
}

Foam::tmp<Foam::volScalarField>
Foam::viscosityModels::muIv::calcPeff() const
{
    const objectRegistry& db = U_.db();
    const Time& runTime= db.time();
    if (db.foundObject<volScalarField>("p") && runTime.timeIndex() > 1) {
        // Info<< "Calculate I based on pressure" << endl;
        const volScalarField& ptot = U_.mesh().lookupObject<volScalarField>("p");
        const volScalarField& gh = U_.mesh().lookupObject<volScalarField>("gh");
        if (rmHydWaterP_) {
            Info<< "Hydrostatic pressure of air phase removed" << endl;
            return max(ptot - rhoWater_*gh, pMin_);
        } else {
            return max(ptot, pMin_);
        }

    } else {
        Info<< "Effective pressure not calculated, return zero" << endl;
        return  tmp<volScalarField>
        (
            new volScalarField
            (
                IOobject
                (
                    "peff0",
                    U_.time().timeName(),
                    U_.db(),
                    IOobject::NO_READ,
                    IOobject::NO_WRITE,
                    false
                ),
               U_.mesh(),
               dimensionedScalar("peff0", dimPressure, SMALL)
           )
        );
    }
}

Foam::tmp<Foam::volScalarField>
Foam::viscosityModels::muIv::calcNormD() const
{
    // note this is different than the classical OpenFOAM strainRate
    return max(mag(symm(fvc::grad(U_)))/sqrt(2.0), dimensionedScalar ("vSmall", dimless/dimTime, VSMALL));
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::viscosityModels::muIv::muIv
(
    const word& name,
    const dictionary& viscosityProperties,
    const volVectorField& U,
    const surfaceScalarField& phi
)
:
    viscosityModel(name, viscosityProperties, U, phi),
    muIvCoeffs_(viscosityProperties.optionalSubDict(typeName + "Coeffs")),
    mus_("mus", dimless, muIvCoeffs_),
    mud_("mud", dimless, muIvCoeffs_),
    I0_("I0", dimless, muIvCoeffs_),
    dg_("dg", dimLength, muIvCoeffs_),
    rhog_("rhog", dimDensity, muIvCoeffs_),
    nuMax_("nuMax", dimViscosity, muIvCoeffs_),
    nuMin_("nuMin", dimViscosity, muIvCoeffs_),
    pMin_("pMin", dimPressure, muIvCoeffs_),
    nuWater_("nuWater", dimViscosity, muIvCoeffs_),
    rhoWater_("rhoWater", dimDensity, muIvCoeffs_),
    alphaSmall_("alphaSmall", dimless, muIvCoeffs_),
    rmHydWaterP_(muIvCoeffs_.lookupOrDefault<Switch>("rmHydWaterP", false)),
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
    mu_
    (
        IOobject
        (
            "mu",
            U_.time().timeName(),
            U_.db(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        calcMu()
    ),
    Iv_
    (
        IOobject
        (
            "Iv",
            U_.time().timeName(),
            U_.db(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        calcIv()
    ),
    peff_
    (
        IOobject
        (
            "peff",
            U_.time().timeName(),
            U_.db(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        calcPeff()
    ),
    normD_
    (
        IOobject
        (
            "normD",
            U_.time().timeName(),
            U_.db(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        calcNormD()
    )
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

bool Foam::viscosityModels::muIv::read
(
    const dictionary& viscosityProperties
)
{
    viscosityModel::read(viscosityProperties);

    muIvCoeffs_ = viscosityProperties.optionalSubDict(typeName + "Coeffs");

    muIvCoeffs_.lookup("mus") >> mus_;
    muIvCoeffs_.lookup("mud") >> mud_;
    muIvCoeffs_.lookup("I0") >> I0_;
    muIvCoeffs_.lookup("dg") >> dg_;
    muIvCoeffs_.lookup("rhog") >> rhog_;
    muIvCoeffs_.lookup("nuMax") >> nuMax_;
    muIvCoeffs_.lookup("nuMin") >> nuMin_;
    muIvCoeffs_.lookup("pMin") >> pMin_;
    muIvCoeffs_.lookup("rmHydWaterP_") >> rmHydWaterP_;
    muIvCoeffs_.lookup("rhoWater") >> rhoWater_;
    muIvCoeffs_.lookup("nuWater") >> nuWater_;  
  muIvCoeffs_.lookup("alphaSmall") >> alphaSmall_;
    return true;
}


// ************************************************************************* //
