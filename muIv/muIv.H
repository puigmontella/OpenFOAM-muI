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

Class
    Foam::viscosityModels::muIv

Description
    Standard mu(I) granular rheology. Following
    Barker, T., & Gray, J. (2017).
    Partial regularisation of the incompressible 𝜇(I)-rheology for granular
    flow. Journal of Fluid Mechanics, 828, 5-32. doi: 10.1017/jfm.2017.428

Author
    Alexander Jarosch research@alexj.at

SourceFiles
    muIv.C

\*---------------------------------------------------------------------------*/

#ifndef muIv_H
#define muIv_H

#include "viscosityModel.H"
#include "dimensionedScalar.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace viscosityModels
{

/*---------------------------------------------------------------------------*\
                           Class muI Declaration
\*---------------------------------------------------------------------------*/

class muIv
:
    public viscosityModel
{
    // Private Data

        dictionary muIvCoeffs_;

        dimensionedScalar mus_;
        dimensionedScalar mud_;
        dimensionedScalar I0_;
        dimensionedScalar dg_;
        dimensionedScalar rhog_;
        dimensionedScalar nuMax_;
        dimensionedScalar nuMin_;
        dimensionedScalar pMin_;
        dimensionedScalar rhoWater_;
        dimensionedScalar nuWater_;
        dimensionedScalar alphaSmall_;

        //- Switch to remove hydrostatic air pressure
        Switch rmHydWaterP_;

        volScalarField nu_;
        mutable volScalarField mu_;
        mutable volScalarField Iv_;
        mutable volScalarField peff_;
        mutable volScalarField normD_;


    // Private Member Functions

        //- Calculate and return the laminar viscosity
        tmp<volScalarField> calcNu() const;
        //- Calculate and return the mu(I) rheology number
        tmp<volScalarField> calcMu() const;
        //- Calculate and return the inertial I number
        tmp<volScalarField> calcIv() const;
        //- Calculate and return the effective pressure
        tmp<volScalarField> calcPeff() const;
        //- Calculate and return the effective pressure
        tmp<volScalarField> calcNormD() const;


public:

    //- Runtime type information
    TypeName("muIv");


    // Constructors

        //- Construct from components
        muIv
        (
            const word& name,
            const dictionary& viscosityProperties,
            const volVectorField& U,
            const surfaceScalarField& phi
        );


    //- Destructor
    virtual ~muIv()
    {}


    // Member Functions

        //- Return the laminar viscosity
        virtual tmp<volScalarField> nu() const
        {
            return nu_;
        }

        //- Return the mu(J) friction
        virtual tmp<volScalarField> mu() const
        {
            return mu_;
        }

        //- Return the viscous number Iv
        virtual tmp<volScalarField> Iv() const
        {
            return Iv_;
        }

        //- Return the effective pressure
        virtual tmp<volScalarField> peff() const
        {
            return peff_;
        }

        //- Return the norm of the deviatoric strain rate
        virtual tmp<volScalarField> normD() const
        {
            return normD_;
        }

        //- Return the laminar viscosity for patch
        virtual tmp<scalarField> nu(const label patchi) const
        {
            return nu_.boundaryField()[patchi];
        }

        //- Correct the laminar viscosity
        virtual void correct()
        {
            normD_ = calcNormD();
            peff_ = calcPeff();
            Iv_ = calcIv();
            mu_ = calcMu();
            nu_ = calcNu();
        }

        //- Read transportProperties dictionary
        virtual bool read(const dictionary& viscosityProperties);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace viscosityModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
