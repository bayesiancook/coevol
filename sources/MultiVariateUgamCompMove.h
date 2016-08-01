#ifndef MVUGAM_H
#define MVUGAM_H

#include "RandomTypes.h"
#include "MultiVariateTreeProcess.h"
#include "BranchProcess.h"

class MultiVariateUgamCompMove : public MCUpdate, public Mnode {

	MultiVariateTreeProcess* tree;
	BranchProcess<PosReal>* ugamtree;
	double tuning;
	int index;

	public:

	MultiVariateUgamCompMove(BranchProcess<PosReal>* inugamtree, MultiVariateTreeProcess* intree, int inindex, double intuning){
		tree = intree;
		ugamtree = inugamtree;
		tuning = intuning;
		index = inindex;
		tree->RecursiveRegister(this,tree->GetRoot());
		ugamtree->RecursiveRegister(this,tree->GetRoot());
	}

	double Move(double tuning_modulator){
		Corrupt(true);
		double u = tuning * tuning_modulator * (Random::Uniform() - 0.5);
		double m = exp(-u);
		tree->PiecewiseTranslation(u,index,1);
		int n = ugamtree->ScalarMultiplication(m);
		double loghastings = -n * u;
		double logratio = loghastings + Update();
		// cerr << u << '\t' << logratio << '\n';
		bool accepted = (log(Random::Uniform()) < logratio);
		if (! accepted)	{
			Corrupt(false);
			Restore();
		}
		return (double) accepted;
	}
};

class UniVariateUgamCompMove : public MCUpdate, public Mnode {

	LogNormalTreeProcess* tree;
	BranchProcess<PosReal>* ugamtree;
	double tuning;

	public:

	UniVariateUgamCompMove(BranchProcess<PosReal>* inugamtree, LogNormalTreeProcess* intree, double intuning){
		tree = intree;
		ugamtree = inugamtree;
		tuning = intuning;
		tree->RecursiveRegister(this,tree->GetRoot());
		ugamtree->RecursiveRegister(this,tree->GetRoot());
	}

	double Move(double tuning_modulator){
		Corrupt(true);
		double u = tuning * tuning_modulator * (Random::Uniform() - 0.5);
		double m = exp(-u);
		tree->Translation(u);
		int n = ugamtree->ScalarMultiplication(m);
		double loghastings = -n * u;
		double logratio = loghastings + Update();
		bool accepted = (log(Random::Uniform()) < logratio);
		if (! accepted)	{
			Corrupt(false);
			Restore();
		}
		return (double) accepted;
	}
};

/*
class MultiVariateUgamLocalCompMove : public MCUpdate {

	MultiVariateTreeProcess* tree;
	BranchProcess<PosReal>* ugamtree;
	double tuning;
	int index;

	public:

	MultiVariateUgamCompMove(BranchProcess<PosReal>* inugamtree, MultiVariateTreeProcess* intree, int inindex, double intuning){
		tree = intree;
		ugamtree = inugamtree;
		tuning = intuning;
		index = inindex;
	}

	double Move(double tuning_modulator){
		Corrupt(true);
		double u = tuning * tuning_modulator * (Random::Uniform() - 0.5);
		double m = exp(-u);
		tree->PiecewiseTranslation(u,index,1);
		int n = ugamtree->ScalarMultiplication(m);
		double loghastings = -n * u;
		double logratio = loghastings + Update();
		// cerr << u << '\t' << logratio << '\n';
		bool accepted = (log(Random::Uniform()) < logratio);
		if (! accepted)	{
			Corrupt(false);
			Restore();
		}
		return (double) accepted;
	}
};
*/

#endif

