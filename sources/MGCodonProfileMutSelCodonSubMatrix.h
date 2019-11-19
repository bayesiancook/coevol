
#ifndef MGCODONPROFILEMUTSELCODONSUBMATRIX_H
#define MGCODONPROFILEMUTSELCODONSUBMATRIX_H

#include "CodonSubMatrix.h"
#include "Var.h"

class MGCodonProfileMutSelCodonSubMatrix : public NucCodonSubMatrix   {

        public:

        MGCodonProfileMutSelCodonSubMatrix(CodonStateSpace* instatespace, SubMatrix* inNucMatrix, Var<Profile>* inCodonProfile, Var<PosReal>* inNeff, bool innormalise = false) :
                        SubMatrix(instatespace->GetNstate(), innormalise),
                        CodonSubMatrix(instatespace, innormalise),
                        NucCodonSubMatrix(instatespace, inNucMatrix, innormalise),
                        CodonProfile(inCodonProfile),
			Neff(inNeff) {}

        Var<Profile>*   GetCodonProfile() {return CodonProfile;}

                protected:

        void                    ComputeArray(int state);
        void                    ComputeStationary();
        void                    SetCodonProfile(Var<Profile>* inCodonProfile) {CodonProfile = inCodonProfile;}
	void                    SetNeff(Var<PosReal>* inNeff)   {Neff = inNeff;}

        // data members
       
        Var<Profile>* CodonProfile;
	Var<PosReal>* Neff;
};


class RandomMGCodonProfileMutSelCodonSubMatrix : public MGCodonProfileMutSelCodonSubMatrix, public RandomCodonSubMatrix   {

        public:

        RandomMGCodonProfileMutSelCodonSubMatrix(CodonStateSpace* instatespace, RandomSubMatrix* inmatrix, Var<Profile>* inCodonProfile, Var<PosReal>* inNeff, bool innormalise = false) :
                        SubMatrix(instatespace->GetNstate(), innormalise),
                        CodonSubMatrix(instatespace,innormalise) ,
			RandomSubMatrix(instatespace->GetNstate(),innormalise),
                        MGCodonProfileMutSelCodonSubMatrix(instatespace,inmatrix,inCodonProfile, inNeff, innormalise) ,
                        RandomCodonSubMatrix(instatespace, innormalise),
                        matrix(inmatrix) ,
                        RandomCodonProfile(inCodonProfile),
			RandomNeff(inNeff) {

                Register(matrix);
                Register(RandomCodonProfile);
		Register(RandomNeff);
		specialUpdate();

        }

        void SetParameters()    {
                SetNucMatrix(matrix);
                SetCodonProfile(RandomCodonProfile);
		SetNeff(RandomNeff);
        }

	void ToStream(ostream& os)	{
		MGCodonProfileMutSelCodonSubMatrix::ToStream(os);
		os << *RandomNeff << '\n';
	}

        protected:

        RandomSubMatrix* matrix;
        Var<Profile>* RandomCodonProfile;
	Var<PosReal>* RandomNeff;
};

#endif

