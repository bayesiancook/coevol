
#ifndef CHRONOCOMP_H
#define CHRONOCOMP_H

#include "PrecisionNormalTreeProcess.h"
#include "PartitionMultiVariateTreeProcess.h"
#include "Chronogram.h"
#include "Move.h"

class MultiVariateChronoCompMove : public MCUpdate {

	Chronogram* chrono;
	MultiVariateTreeProcess* process;
	int index;
	double tuning;

	public:

	MultiVariateChronoCompMove(Chronogram* inchrono, MultiVariateTreeProcess* inprocess, int inindex, double intuning)	{
		chrono = inchrono;
		process = inprocess;
		index = inindex;
		if (index < 0)	{
			cerr << "error : index out of bound\n";
			exit(1);
		}
		tuning = intuning;
	}

	Tree* GetTree() {return chrono->GetTree();}

	const Link* GetRoot() {return GetTree()->GetRoot();}

	double Move(double tuning_modulator){

		bool accepted = false;
		Mnode* mnode = new Mnode(false);

		const Link* from = 0;
		do	{
			from = GetTree()->ChooseInternalNode();
		} while ((! from) || (from == GetRoot()));

		double t = tuning * tuning_modulator;
		double m = t * (Random::Uniform() - 0.5);
		double e = exp(m);
		chrono->SetMinMax(from);
		double min = chrono->GetNodeDate(from->GetNode())->GetMin();
		double max = chrono->GetNodeDate(from->GetNode())->GetMax();
		double time = chrono->GetAbsoluteTime(from);
		if ((e*time > min) && (e*time < max))	{
			chrono->RecursiveRegister(mnode,from);
			process->RecursiveRegister(mnode,from);
			mnode->Corrupt(true);
			chrono->MultiplyTimes(from,e);
			process->RecursivePiecewiseTranslation(from,-m,index,1);
			int nnode = chrono->RecursiveGetNinternalNode(from);
			double logratio = mnode->Update();
			logratio += nnode * m;
			bool accepted = (log(Random::Uniform()) < logratio);
			if (! accepted)	{
				mnode->Corrupt(false);
				mnode->Restore();
			}
			delete mnode;
		}

		return (double) accepted;
	}

};


class MultiVariateScaleCompMove : public MCUpdate, public Mnode {

	Rvar<PosReal>* scale;
	MultiVariateTreeProcess* process;
	int index;
	double tuning;

	public:

	MultiVariateScaleCompMove(Rvar<PosReal>* inscale, MultiVariateTreeProcess* inprocess, int inindex, double intuning)	{
		scale = inscale;
		process = inprocess;
		scale->Register(this);
		process->RecursiveRegister(this,GetRoot());
		index = inindex;
		tuning = intuning;
	}

	const Tree* GetTree() {return process->GetTree();}

	const Link* GetRoot() {return GetTree()->GetRoot();}

	double Move(double tuning_modulator){

		Corrupt(true);
		double t = tuning * tuning_modulator;
		double m = t * (Random::Uniform() - 0.5);
		double e = exp(m);
		scale->ScalarMultiplication(e);
		process->RecursivePiecewiseTranslation(GetRoot(),-m,index,1);
		process->RootNeighborPiecewiseTranslation(2*m,index,1);
		double logratio = Update();
		logratio += m;
		bool accepted = (log(Random::Uniform()) < logratio);
		if (! accepted)	{
			Corrupt(false);
			Restore();
		}
		return (double) accepted;
	}

};


class RootSigmaMove : public MCUpdate, public Mnode	{

	double tuning;
	CalibratedChronogram* chrono;
	// ConjugatePartitionMultiVariateTreeProcess* process;
	Var<CovMatrix>* mat;

	public:

	RootSigmaMove(CalibratedChronogram* inchrono, Var<CovMatrix>* inmat, double intuning)	{
	// RootSigmaMove(CalibratedChronogram* inchrono, ConjugatePartitionMultiVariateTreeProcess* inprocess, double intuning)	{
		tuning = intuning;
		chrono = inchrono;
		mat = inmat;
		// process = inprocess;
		chrono->RecursiveRegister(this,chrono->GetRoot());
		chrono->GetScale()->Register(this);
		mat->Register(this);
	}

	double Move(double tuning_modulator){
		Corrupt(true);
		double m = tuning * tuning_modulator * (Random::Uniform() - 0.5);
		double e = exp(m);
		double f = exp(-m);
		chrono->GetScale()->ScalarMultiplication(e);
		chrono->MultiplyTimes(chrono->GetRoot(),f);
		chrono->GetNodeDate(chrono->GetRoot()->GetNode())->setval(1.0);
		mat->ScalarMultiplication(e);
		double logratio = Update();
		logratio -= (chrono->GetNinternalNode() - 2) * m;
		// logratio += mat->GetDim() * (mat->GetDim() + 1) / 2 * m;
		logratio += mat->GetDim() * m;
		bool accepted = (log(Random::Uniform()) < logratio);
		if (! accepted)	{
			Corrupt(false);
			Restore();
		}
		return ((double) accepted);
	}
};

class RootSigmaProcessMove : public MCUpdate, public Mnode	{

	double tuning;
	CalibratedChronogram* chrono;
	MultiVariateTreeProcess* process;
	int index;
	Var<CovMatrix>* mat;

	public:

	RootSigmaProcessMove(CalibratedChronogram* inchrono, Var<CovMatrix>* inmat, MultiVariateTreeProcess* inprocess, int inindex, double intuning)	{
		tuning = intuning;
		chrono = inchrono;
		mat = inmat;
		process = inprocess;
		index = inindex;
		chrono->RecursiveRegister(this,chrono->GetRoot());
		chrono->GetScale()->Register(this);
		process->RecursiveRegister(this,GetRoot());
		mat->Register(this);
	}

	const Tree* GetTree() {return process->GetTree();}

	const Link* GetRoot() {return GetTree()->GetRoot();}

	double Move(double tuning_modulator){
		Corrupt(true);
		double m = tuning * tuning_modulator * (Random::Uniform() - 0.5);
		double e = exp(m);
		double f = exp(-m);
		chrono->GetScale()->ScalarMultiplication(e);
		chrono->MultiplyTimes(chrono->GetRoot(),f);
		chrono->GetNodeDate(chrono->GetRoot()->GetNode())->setval(1.0);
		process->RecursivePiecewiseTranslation(GetRoot(),-m,index,1);
		mat->ScalarMultiplication(e);
		double logratio = Update();
		logratio -= (chrono->GetNinternalNode() - 2) * m;
		logratio += mat->GetDim() * m;
		bool accepted = (log(Random::Uniform()) < logratio);
		if (! accepted)	{
			Corrupt(false);
			Restore();
		}
		return ((double) accepted);
	}
};


class LogNormalChronoCompMove : public MCUpdate {

	Chronogram* chrono;
	LogNormalTreeProcess* process;
	double tuning;

	public:

	LogNormalChronoCompMove(Chronogram* inchrono, LogNormalTreeProcess* inprocess, double intuning)	{
		chrono = inchrono;
		process = inprocess;
		tuning = intuning;
	}

	Tree* GetTree() {return chrono->GetTree();}

	const Link* GetRoot() {return GetTree()->GetRoot();}

	double Move(double tuning_modulator){

		Mnode* mnode = new Mnode(false);

		const Link* from = 0;
		do	{
			from = GetTree()->ChooseInternalNode();
		} while ((! from) || (from == GetRoot()));

		chrono->RecursiveRegister(mnode,from);
		process->RecursiveRegister(mnode,from);

		mnode->Corrupt(true);

		double t = tuning * tuning_modulator;
		double m = t * (Random::Uniform() - 0.5);
		double e = exp(m);
		chrono->MultiplyTimes(from,e);
		process->RecursiveTranslation(from,-m);
		int nnode = chrono->RecursiveGetNinternalNode(from);
		double logratio = mnode->Update();
		logratio += nnode * m;
		bool accepted = (log(Random::Uniform()) < logratio);
		if (! accepted)	{
			mnode->Corrupt(false);
			mnode->Restore();
		}
		delete mnode;

		return (double) accepted;
	}

};


class LogNormalScaleCompMove : public MCUpdate, public Mnode {

	Rvar<PosReal>* scale;
	LogNormalTreeProcess* process;
	double tuning;

	public:

	LogNormalScaleCompMove(Rvar<PosReal>* inscale, LogNormalTreeProcess* inprocess, double intuning)	{
		scale = inscale;
		process = inprocess;
		scale->Register(this);
		process->RecursiveRegister(this,GetRoot());
		tuning = intuning;
	}

	const Tree* GetTree() {return process->GetTree();}

	const Link* GetRoot() {return GetTree()->GetRoot();}

	double Move(double tuning_modulator){

		Corrupt(true);
		double t = tuning * tuning_modulator;
		double m = t * (Random::Uniform() - 0.5);
		double e = exp(m);
		scale->ScalarMultiplication(e);
		process->RecursiveTranslation(GetRoot(),-m);
		double logratio = Update();
		logratio += m;
		bool accepted = (log(Random::Uniform()) < logratio);
		if (! accepted)	{
			Corrupt(false);
			Restore();
		}
		return (double) accepted;
	}

};

#endif
