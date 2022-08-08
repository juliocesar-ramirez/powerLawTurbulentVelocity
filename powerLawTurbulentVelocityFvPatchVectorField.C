/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2022 OpenFOAM Foundation
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

#include "DimensionedScalarField.H"
#include "GeometricScalarField.H"
#include "UList.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "longDoubleScalar.H"
#include "mathematicalConstants.H"
#include "powerLawTurbulentVelocityFvPatchVectorField.H"
#include "primitiveFieldsFwd.H"
#include "scalar.H"
#include "scalarField.H"
#include "surfaceFields.H"
#include "vector.H"
#include "vectorField.H"
#include "volFields.H"
#include <cmath>
#include <math.h>

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::scalar Foam::powerLawTurbulentVelocityFvPatchVectorField::t() const
{
    return db().time().timeOutputValue();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::powerLawTurbulentVelocityFvPatchVectorField::
    powerLawTurbulentVelocityFvPatchVectorField(
        const fvPatch &p, const DimensionedField<vector, volMesh> &iF)
    : fixedValueFvPatchVectorField(p, iF), a_(0.0), alpha_(0.0), zref_(0.0),
      N_(0.0), n_(Zero), y_(Zero), wordData_("wordDefault"), labelData_(-1),
      boolData_(false) {}

Foam::powerLawTurbulentVelocityFvPatchVectorField::
    powerLawTurbulentVelocityFvPatchVectorField(
        const fvPatch &p, const DimensionedField<vector, volMesh> &iF,
        const dictionary &dict)
    : fixedValueFvPatchVectorField(p, iF), a_(dict.lookup<scalar>("a")),
      alpha_(dict.lookup<scalar>("alpha")), zref_(dict.lookup<scalar>("zref")),
      N_(dict.lookup<scalar>("N")),
      n_(dict.lookup<vector>("n")), y_(dict.lookup<vector>("y")),
      wordData_(dict.lookupOrDefault<word>("wordName", "wordDefault")),
      labelData_(-1), boolData_(false) {

  fixedValueFvPatchVectorField::evaluate();

  /*
  // Initialise with the value entry if evaluation is not possible
  fvPatchVectorField::operator=
  (
      vectorField("value", dict, p.size())
  );
  */
}

Foam::powerLawTurbulentVelocityFvPatchVectorField::
    powerLawTurbulentVelocityFvPatchVectorField(
        const powerLawTurbulentVelocityFvPatchVectorField &ptf,
        const fvPatch &p, const DimensionedField<vector, volMesh> &iF,
        const fvPatchFieldMapper &mapper)
    : fixedValueFvPatchVectorField(ptf, p, iF, mapper), a_(ptf.a_),
      alpha_(ptf.alpha_), zref_(ptf.zref_),N_(ptf.N_), n_(ptf.n_), y_(ptf.y_) ,
      wordData_(ptf.wordData_), labelData_(-1), boolData_(ptf.boolData_) {}

Foam::powerLawTurbulentVelocityFvPatchVectorField::
    powerLawTurbulentVelocityFvPatchVectorField(
        const powerLawTurbulentVelocityFvPatchVectorField &ptf,
        const DimensionedField<vector, volMesh> &iF)
    : fixedValueFvPatchVectorField(ptf, iF), a_(ptf.a_), alpha_(ptf.alpha_),
      zref_(ptf.zref_),N_(ptf.N_), n_(ptf.n_), y_(ptf.y_), wordData_(ptf.wordData_),
      labelData_(-1), boolData_(ptf.boolData_) {}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::powerLawTurbulentVelocityFvPatchVectorField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    fixedValueFvPatchVectorField::autoMap(m);
}


void Foam::powerLawTurbulentVelocityFvPatchVectorField::rmap
(
    const fvPatchVectorField& ptf,
    const labelList& addr
)
{
    fixedValueFvPatchVectorField::rmap(ptf, addr);


}


void Foam::powerLawTurbulentVelocityFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    vectorField cellxyz = patch().Cf();
    scalarField celly = cellxyz & y_;
    // power law
    vectorField Upl = n_ * a_ * pow(celly, alpha_);
    vectorField Ut = Upl * 0;
    // velocidad turbulenta
    scalarField Af = patch().magSf();
    scalar A = gSum(Af);

    boundBox bb(patch().patch().localPoints(), true);
    vector bmax = bb.max();
    scalar ly = bmax[1];
    scalar lz = bmax[2];
    // turbulent lenght
    scalar L = (2 * ly * lz) / (ly + ly);
    scalar sigma = 0.1 * L / 2;
    scalar sigma2 = 2 * sqr(sigma);

    scalar z0 = zref_ / exp(1 / alpha_);
    scalarField Iu = 1 / (log(celly / z0));
    scalarField k = 1.5 * pow(Upl.component(0), 2) * pow(Iu, 2);
    scalar a = 0.0;
    vector z1 = vector(1, 0, 0);

    forAll(patch(), facei) {
      vector cell = cellxyz[facei];
      scalar celly = cell[1];
      scalar cellz = cell[2];

      srand(time(NULL));
      vector Utx = vector(0, 0, 0);
      for (int i = 0; i < N_; i++) {
          // valores indices del patch random
        int random = rand() % cellxyz.size();

        vector xi = cellxyz[random];
        if (xi == cell) {
          scalar lxi = sqrt(patch().magSf()[i]);
          scalar add = 0.1 * lxi;
          xi[1] = xi[1] + add;
          xi[2] = xi[2] + add;
        }
        scalar xiy = xi[1];
        scalar xiz = xi[2];
        scalar dist = sqr(celly - xiy) + sqr(cellz - xiz);
        scalar ks = k[random];
        scalar circulacion =
            4 * sqrt(Foam::constant::mathematical::pi * A * ks) /
            (3 * N_ * (2 * log(3.0) - 3 * log(2.0)));

	scalar delta = sqrt(patch().magSf()[random]); //Delta calculated from face area
	if (sigma <= delta)
	{
	sigma = delta; 
	sigma2 = 2*sqr(sigma);
	}

        Utx += (1 / (Foam::constant::mathematical::twoPi)) * circulacion *
            (((xi - cell) ^ z1)/dist) * (1 - exp(-dist / sigma2)) *
            exp(-dist / sigma2);
      }
      Ut[facei] =  Utx;
    }

    fixedValueFvPatchVectorField::operator==(Upl + Ut);

    fixedValueFvPatchVectorField::updateCoeffs();
    }

void Foam::powerLawTurbulentVelocityFvPatchVectorField::write
(
    Ostream& os
) const
{
    fvPatchVectorField::write(os);
    writeEntry(os, "a", a_);
    writeEntry(os, "alpha", alpha_);
    writeEntry(os, "zref", zref_);
    writeEntry(os, "N", N_);
    writeEntry(os, "n", n_);
    writeEntry(os, "y", y_);
    writeEntry(os, "wordData", wordData_);
    writeEntry(os, "value", *this);
}


// * * * * * * * * * * * * * * Build Macro Function  * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchVectorField,
        powerLawTurbulentVelocityFvPatchVectorField
    );
}

// ************************************************************************* //
