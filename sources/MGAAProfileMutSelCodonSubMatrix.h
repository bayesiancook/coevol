
#ifndef MGAAPROFILEMUTSELCODONSUBMATRIX_H
#define MGAAPROFILEMUTSELCODONSUBMATRIX_H

#include "CodonSubMatrix.h"
// #include "CovMatrix.h"
#include "Var.h"

static const double TOOSMALL = 1e-20;
static const double TOOLARGE = 50;
static const double TOOLARGENEGATIVE = -50;

class MGAAProfileMutSelCodonSubMatrix : public NucCodonSubMatrix   {

        public:

        MGAAProfileMutSelCodonSubMatrix(CodonStateSpace* instatespace, SubMatrix* inNucMatrix, Var<Profile>* inAAProfile, Var<PosReal>* inNeff, bool innormalise = false) :
                        SubMatrix(instatespace->GetNstate(), innormalise),
                        CodonSubMatrix(instatespace, innormalise),
                        NucCodonSubMatrix(instatespace, inNucMatrix, innormalise),
                        AAProfile(inAAProfile),
			Neff(inNeff) {}

        Var<Profile>*   GetAAProfile() {return AAProfile;}

                protected:

        void                    ComputeArray(int state);
        void                    ComputeStationary();
        void                    SetAAProfile(Var<Profile>* inAAProfile) {AAProfile = inAAProfile;}
	void                    SetNeff(Var<PosReal>* inNeff)   {Neff = inNeff;}

        // data members
       
        Var<Profile>* AAProfile;
	Var<PosReal>* Neff;
};


class RandomMGAAProfileMutSelCodonSubMatrix : public MGAAProfileMutSelCodonSubMatrix, public RandomCodonSubMatrix   {

        public:

        RandomMGAAProfileMutSelCodonSubMatrix(CodonStateSpace* instatespace, RandomSubMatrix* inmatrix, Var<Profile>* inAAProfile, Var<PosReal>* inNeff, bool innormalise = false) :
                        SubMatrix(instatespace->GetNstate(), innormalise),
                        CodonSubMatrix(instatespace, innormalise),
			RandomSubMatrix(instatespace->GetNstate(),innormalise),
                        MGAAProfileMutSelCodonSubMatrix(instatespace,inmatrix,inAAProfile, inNeff, innormalise) ,
                        RandomCodonSubMatrix(instatespace, innormalise),
                        matrix(inmatrix) ,
                        RandomAAProfile(inAAProfile),
			RandomNeff(inNeff) {

                Register(matrix);
                Register(RandomAAProfile);
		Register(RandomNeff);
		specialUpdate();

        }

        void SetParameters()    {
                SetNucMatrix(matrix);
                SetAAProfile(RandomAAProfile);
		SetNeff(RandomNeff);
        }

	void ToStream(ostream& os)	{
		MGAAProfileMutSelCodonSubMatrix::ToStream(os);
		os << *RandomNeff << '\n';
	}

        protected:

        RandomSubMatrix* matrix;
        Var<Profile>* RandomAAProfile;
	Var<PosReal>* RandomNeff;
};

#endif

