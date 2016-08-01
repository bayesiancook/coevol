#ifndef MVRELRATE_H
#define MVRELRATE_H

#include "RandomTypes.h"
#include "SimilarityMatrix.h"
#include "MultiVariateTreeProcess.h"
#include "LinRegContSub.h"

class MultiVariateRelRateCompensatoryMove : public MCUpdate, public Mnode {

	MultiVariateTreeProcess* tree;
	Rvar<PosRealVector>* relrate;
	SimilarityMatrix* simmat;
	double tuning;
	int index;

	public:

	MultiVariateRelRateCompensatoryMove(MultiVariateTreeProcess* intree, Rvar<PosRealVector>* inrelrate, SimilarityMatrix* insimmat, double intuning, int inindex){
		tree = intree;
		relrate = inrelrate;
		int Nrr = insimmat->GetDim() * (insimmat->GetDim() - 1) / 2;
		if (relrate->GetDim() != Nrr)	{
			cerr << "error in multivar relrate move\n";
			cerr << Nrr << '\t' << relrate->GetDim() << '\n';
			exit(1);
		}
		simmat = insimmat;
		tuning = intuning;
		index = inindex;
		tree->RecursiveRegister(this,tree->GetRoot());
		relrate->Register(this);
	}

	double Move(double tuning_modulator){
		Corrupt(true);
		double u = tuning * tuning_modulator * (Random::Uniform() - 0.5);
		double m = exp(-u);
		tree->PiecewiseTranslation(u,index,1);
		int k = 0;
		double radk = 0;
		for (int i=0; i<simmat->GetDim(); i++)	{
			for (int j=i+1; j<simmat->GetDim(); j++)	{
				if (simmat->isRadical(i,j))	{
					(*relrate)[k] *= m;
					radk++;
				}
				k++;
			}
		}
		double logratio = -radk*u + Update();
		// cerr << u << '\t' << logratio << '\n';
		bool accepted = (log(Random::Uniform()) < logratio);
		if (! accepted)	{
			Corrupt(false);
			Restore();
		}
		return (double) accepted;
	}
};

class MultiVariateRelRateMutCompensatoryMove : public MCUpdate, public Mnode {

	MultiVariateTreeProcess* tree;
	Rvar<PosRealVector>* relrate;
	double tuning;
	int tstvindex;
	int tvgcindex;

	public:

	MultiVariateRelRateMutCompensatoryMove(MultiVariateTreeProcess* intree, Rvar<PosRealVector>* inrelrate, double intuning, int intstvindex, int intvgcindex){
		tree = intree;
		relrate = inrelrate;
		int Nrr = Nnuc * (Nnuc - 1) / 2;
		if (relrate->GetDim() != Nrr)	{
			cerr << "error in multivar relrate move\n";
			cerr << Nrr << '\t' << relrate->GetDim() << '\n';
			exit(1);
		}
		tuning = intuning;
		tstvindex = intstvindex;
		tvgcindex = intvgcindex;
		tree->RecursiveRegister(this,tree->GetRoot());
		relrate->Register(this);
	}

	double Move(double tuning_modulator){
		Corrupt(true);
		double utstv = tuning * tuning_modulator * (Random::Uniform() - 0.5);
		double mtstv = exp(-utstv);
		double utvgc = tuning * tuning_modulator * (Random::Uniform() - 0.5);
		double mtvgc = exp(-utvgc);
		tree->PiecewiseTranslation(utstv,tstvindex,1);
		if (tvgcindex != -1)	{
			tree->PiecewiseTranslation(utvgc,tvgcindex,1);
		}

		/*
		0 ac tvgc
		1 ag ts
		2 at tv0
		3 cg tv0
		4 ct ts
		5 gt tvgc
		*/

		int ktstv = 0;
		int ktvgc = 0;
		if (tvgcindex != -1)	{
			(*relrate)[0] *= mtvgc;
			(*relrate)[1] *= mtstv;
			(*relrate)[4] *= mtstv;
			(*relrate)[5] *= mtvgc;
			ktstv = 2;
			ktvgc = 2;
		}
		else	{
			(*relrate)[1] *= mtstv;
			(*relrate)[4] *= mtstv;
			ktstv = 2;
			ktvgc = 0;
		}
		double logratio = -ktstv*utstv - ktvgc*utvgc + Update();

		bool accepted = (log(Random::Uniform()) < logratio);
		if (! accepted)	{
			Corrupt(false);
			Restore();
		}
		return (double) accepted;
	}
};

class SplitMultiVariateRelRateCompensatoryMove : public MCUpdate, public Mnode {

	MultiVariateTreeProcess* tree;
	Rvar<PosRealVector>* relrate;
	SimilarityMatrix* simmat;
	SplitAAMatrix* splitmat;
	int splittype;
	double tuning;
	int tstvindex;
	int omegatsindex;
	int omegatvindex;

	public:

	SplitMultiVariateRelRateCompensatoryMove(MultiVariateTreeProcess* intree, Rvar<PosRealVector>* inrelrate, SplitAAMatrix* insplitmat, int insplittype, SimilarityMatrix* insimmat, double intuning, int intstvindex, int inomegatsindex, int inomegatvindex){
		tree = intree;
		relrate = inrelrate;
		int Nrr = insimmat->GetDim() * (insimmat->GetDim() - 1) / 2;
		if (relrate->GetDim() != Nrr)	{
			cerr << "error in multivar relrate move\n";
			cerr << Nrr << '\t' << relrate->GetDim() << '\n';
			exit(1);
		}
		simmat = insimmat;
		splitmat = insplitmat;
		splittype = insplittype;
		tuning = intuning;
		tstvindex = intstvindex;
		omegatsindex = inomegatsindex;
		omegatvindex = inomegatvindex;
		tree->RecursiveRegister(this,tree->GetRoot());
		relrate->Register(this);
	}

	double Move(double tuning_modulator){
		Corrupt(true);
		double utstv = tuning * tuning_modulator * (Random::Uniform() - 0.5);
		double mtstv = exp(-utstv);
		double uomts = tuning * tuning_modulator * (Random::Uniform() - 0.5);
		double momts = exp(-uomts);
		double uomtv = tuning * tuning_modulator * (Random::Uniform() - 0.5);
		double momtv = exp(-uomtv);
		tree->PiecewiseTranslation(utstv,tstvindex,1);
		if (omegatsindex != -1)	{
			tree->PiecewiseTranslation(uomts,omegatsindex,1);
		}
		tree->PiecewiseTranslation(uomtv,omegatvindex,1);
		int k = 0;
		double ktstv = 0;
		double komts = 0;
		double komtv = 0;
		for (int i=0; i<simmat->GetDim(); i++)	{
			for (int j=i+1; j<simmat->GetDim(); j++)	{
				if (splitmat->IsNonCTNearest(i,j) == 1)	{
					if (simmat->isRadical(i,j))	{
						(*relrate)[k] *= momtv;
						komtv++;
					}
				}
				else	{
					if ((!splittype) || (splitmat->IsNonCTNearest(i,j) != -1))	{
						(*relrate)[k] *= mtstv;
						ktstv++;
						if (omegatsindex != -1)	{
							if (simmat->isRadical(i,j))	{
								(*relrate)[k] *= momts;
								komts++;
							}
						}
					}
				}
				k++;
			}
		}
		double logratio = -ktstv*utstv - komts*uomts - komtv*uomtv + Update();
		// cerr << u << '\t' << logratio << '\n';
		bool accepted = (log(Random::Uniform()) < logratio);
		if (! accepted)	{
			Corrupt(false);
			Restore();
		}
		return (double) accepted;
	}
};

class LinRegRelRateCompensatoryMove : public MCUpdate, public Mnode {

	LinRegNormalProcess* tree;
	Rvar<PosRealVector>* relrate;
	SimilarityMatrix* simmat;
	double tuning;
	int index;

	public:

	LinRegRelRateCompensatoryMove(LinRegNormalProcess* intree, Rvar<PosRealVector>* inrelrate, SimilarityMatrix* insimmat, double intuning, int inindex){
		tree = intree;
		relrate = inrelrate;
		int Nrr = insimmat->GetDim() * (insimmat->GetDim() - 1) / 2;
		if (relrate->GetDim() != Nrr)	{
			cerr << "error in multivar relrate move\n";
			cerr << Nrr << '\t' << relrate->GetDim() << '\n';
			exit(1);
		}
		simmat = insimmat;
		tuning = intuning;
		index = inindex;
		tree->RecursiveRegister(this,tree->GetRoot());
		relrate->Register(this);
	}

	double Move(double tuning_modulator){
		Corrupt(true);
		double u = tuning * tuning_modulator * (Random::Uniform() - 0.5);
		double m = exp(-u);
		tree->PiecewiseTranslation(u,index,1);
		int k = 0;
		double radk = 0;
		for (int i=0; i<simmat->GetDim(); i++)	{
			for (int j=i+1; j<simmat->GetDim(); j++)	{
				if (simmat->isRadical(i,j))	{
					(*relrate)[k] *= m;
					radk++;
				}
				k++;
			}
		}
		double logratio = -radk*u + Update();
		// cerr << u << '\t' << logratio << '\n';
		bool accepted = (log(Random::Uniform()) < logratio);
		if (! accepted)	{
			Corrupt(false);
			Restore();
		}
		return (double) accepted;
	}
};

#endif

