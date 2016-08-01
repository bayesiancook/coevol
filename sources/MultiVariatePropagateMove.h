
#ifndef MVMULTIPROP_H
#define MVMULTIPROP_H

#include "MultiVariateTreeProcess.h"


class MultiVariatePropagateMove : public MCUpdate {

	public:

	MultiVariatePropagateMove(MultiVariateTreeProcess* inprocess, double intuning, double ins01, double ins10){
		process= inprocess;
		tuning = intuning;
		s01 = ins01;
		s10 = ins10;
	}

	int GetDim()	{
		return process->GetDim();
	}

	const Link* GetRoot() {return process->GetRoot();}

	double Move(double tuning_modulator){

		double t = tuning * tuning_modulator;
		double* delta = new double[GetDim()];
		for (int i=0; i<GetDim(); i++)	{
			delta[i] = t * (Random::Uniform() - 0.5);
		}

		Mnode* mnode = new Mnode(false);
		ChooseNodeSetAndMove(GetRoot(),mnode,delta,false);

		double logratio = mnode->Update();
		bool accepted = (log(Random::Uniform()) < logratio);
		if (! accepted)	{
			mnode->Corrupt(false);
			mnode->Restore();
		}
		delete mnode;
		delete[] delta;
		return (double) accepted;
	}

	private:

	bool ChooseNodeSetAndMove(const Link* from, Mnode* mnode, double* delta, bool sw)	{
		bool ret = false;
		if (sw)	{
			if (Random::Uniform() < s10)	{
				sw = false;
			}
		}
		else	{
			if (Random::Uniform() < s01)	{
				sw = true;
				ret = true;
			}
		}
		if (sw)	{
			process->GetNodeVal(from->GetNode())->Register(mnode);
			process->GetNodeVal(from->GetNode())->Corrupt(true);
			process->GetMultiNormal(from->GetNode())->Shift(delta);
		}
		const Link* link = from->Next();
		while ((!ret) && (link!=from))	{
			ret |= ChooseNodeSetAndMove(link->Out(),mnode,delta,sw);
			link = link->Next();
		}
		return ret;
	}

	MultiVariateTreeProcess* process;
	double tuning;
	double s01;
	double s10;
};


#endif

