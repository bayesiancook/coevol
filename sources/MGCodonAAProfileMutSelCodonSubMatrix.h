
#ifndef MGCODONAAPROFILEMUTSELCODONSUBMATRIX_H
#define MGCODONAAPROFILEMUTSELCODONSUBMATRIX_H

#include "CodonSubMatrix.h"
// #include "CovMatrix.h"
#include "Var.h"

class MGCodonAAProfileMutSelCodonSubMatrix : public virtual NucCodonSubMatrix   {

        public:

        MGCodonAAProfileMutSelCodonSubMatrix(CodonStateSpace* instatespace, SubMatrix* inNucMatrix, Var<Profile>* inAAProfile, Var<Profile>* inCodonProfile, Var<PosReal>* inNeff, bool innormalise = false) :
                        SubMatrix(instatespace->GetNstate(), innormalise),
                        CodonSubMatrix(instatespace,innormalise),
                        NucCodonSubMatrix(instatespace, inNucMatrix, innormalise),
                        CodonProfile(inCodonProfile),
			AAProfile(inAAProfile),
			Neff(inNeff) {}

        Var<Profile>*   GetCodonProfile() {return CodonProfile;}
        Var<Profile>*   GetAAProfile() {return AAProfile;}

                protected:

        void                    ComputeArray(int state);
        void                    ComputeStationary();
        void                    SetCodonProfile(Var<Profile>* inCodonProfile) {CodonProfile = inCodonProfile;}
        void                    SetAAProfile(Var<Profile>* inAAProfile) {AAProfile = inAAProfile;}
	void                    SetNeff(Var<PosReal>* inNeff)   {Neff = inNeff;}

        // data members
       
        Var<Profile>* CodonProfile;
        Var<Profile>* AAProfile;
	Var<PosReal>* Neff;
};


class RandomMGCodonAAProfileMutSelCodonSubMatrix : public virtual MGCodonAAProfileMutSelCodonSubMatrix, public virtual RandomCodonSubMatrix   {

        public:

        RandomMGCodonAAProfileMutSelCodonSubMatrix(CodonStateSpace* instatespace, RandomSubMatrix* inmatrix, Var<Profile>* inAAProfile, Var<Profile>* inCodonProfile, Var<PosReal>* inNeff, bool innormalise = false) :
                        SubMatrix(instatespace->GetNstate(), innormalise),
                        CodonSubMatrix(instatespace, innormalise),
                        NucCodonSubMatrix(instatespace, inmatrix, innormalise),
                        MGCodonAAProfileMutSelCodonSubMatrix(instatespace,inmatrix,inAAProfile,inCodonProfile,inNeff,innormalise),
			RandomSubMatrix(instatespace->GetNstate(),innormalise),
                        RandomCodonSubMatrix(instatespace, innormalise),
                        matrix(inmatrix) ,
                        RandomCodonProfile(inCodonProfile),
                        RandomAAProfile(inAAProfile),
			RandomNeff(inNeff) {

                Register(matrix);
                Register(RandomCodonProfile);
                Register(RandomAAProfile);
		Register(RandomNeff);
		specialUpdate();

        }

        void SetParameters()    {
                SetNucMatrix(matrix);
                SetCodonProfile(RandomCodonProfile);
                SetAAProfile(RandomAAProfile);
		SetNeff(RandomNeff);
        }

	void ToStream(ostream& os)	{
		MGCodonAAProfileMutSelCodonSubMatrix::ToStream(os);
		os << *RandomNeff << '\n';
	}

        protected:

        RandomSubMatrix* matrix;
        Var<Profile>* RandomCodonProfile;
        Var<Profile>* RandomAAProfile;
	Var<PosReal>* RandomNeff;
};

#endif

