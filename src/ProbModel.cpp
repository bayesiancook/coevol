
#include <cstdlib>

#include "ProbModel.h"

ProbModel::ProbModel() : scheduler(this) {}

ProbModel::~ProbModel() {}

bool ProbModel::CheckUpdateFlags()	{
	bool ret = true;
	for(clit i=root.begin(); i!=root.end(); ++i)	{
		ret &= (*i)->CheckUpdateFlags();
	}
	return ret;
}

void ProbModel::Register(DAGnode* var)	{
	state.insert(var);
}

void ProbModel::RootRegister(DAGnode* var)	{
	root.insert(var);
}

void ProbModel::Corrupt()	{
	map<DAGnode*,int> m;
	for(clit i=root.begin(); i!=root.end(); ++i)	{
		(*i)->FullCorrupt(m);
	}
}

double ProbModel::Update(bool check)	{
	Corrupt();
	double total = 0;
	for(clit i=root.begin(); i!=root.end(); ++i)	{
		total += (*i)->FullUpdate(check);
	}
	DAGnode::initmode = false;
	return total;
}

void ProbModel::Register()	{

	/*
	if (! state.empty())	{
		cerr << "error : state is not empty\n";
		exit(1);
	}
	Corrupt();
	for(clit i=root.begin(); i!=root.end(); ++i)	{
		(*i)->RecursiveRegister(this);
	}
	cerr << "model size : " << state.size() << '\n';
	*/
}


double ProbModel::Move(double tuning_modulator, int ncycle, bool verbose, bool check)	{
	scheduler.Cycle(tuning_modulator,ncycle,verbose,check);
	return 1;
}

void ProbModel::Monitor(ostream& os, ostream& osdetail)	{
	scheduler.ToStream(os, osdetail);
		Details(osdetail);
}

/*
void ProbModel::Monitor(ostream& os)	{
	scheduler.ToStream(os);
}
*/

