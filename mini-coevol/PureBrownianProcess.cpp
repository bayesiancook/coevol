#include "PureBrownianProcess.h"

PureBrownianProcess::PureBrownianProcess(Chronogram *intree, Var<CovMatrix> *insigma, Var<PosReal>* inagescale, GlobalScalingFunction* inscalefunction)
{
	SetWithRoot(false);
	tree = intree;
	sigma = insigma;
	agescale = inagescale;
	scalefunction = inscalefunction;


	RecursiveCreate(GetRoot());

}

PureBrownianProcess::~PureBrownianProcess()
{
	RecursiveDelete(GetRoot());
}



Rvar<BrownianBridge> * PureBrownianProcess::CreateBranchVal(const Link* link)	{

	if (link->isRoot())	{
		cerr << "error : create branch val called on root\n";
		exit(1);
	}

	// grep the time associated to the branch
	//Var<PosReal>* time = GetLengthTree()->GetBranchLength(link->GetBranch());

	// grep the age associated to the nodes
	Var<PosReal>* up = GetLengthTree()->GetNodeDate(link->GetNode());
	Var<PosReal>* down = GetLengthTree()->GetNodeDate(link->Out()->GetNode());


	// make the new brownian path, and return the pointer
	return new RandomBrownianPath(up, down, sigma);

}


double PureBrownianProcess::Move(double tuning)	{
	int n = 0;
	double tot = Move(this->GetRoot(),tuning,n);
	return tot / n;
}

double PureBrownianProcess::Move(Link* from, double tuning, int& count)	{
	double total = 0;
	if(!from->isRoot()) {
		total = this->GetBranchVal(from->GetBranch())->Move(tuning);
		count++;
	}
	for(Link* link=from->Next(); link!=from; link=link->Next())	{
		total += Move(link->Out(),tuning,count);
	}
	return total;
}

double PureBrownianProcess::GetLogProb()	{
	return GetLogProb(this->GetRoot());
}

double PureBrownianProcess::GetLogProb(Link* from)	{
	double total = 0;
	if(!from->isRoot())
		total = this->GetBranchVal(from->GetBranch())->GetLogProb();

	for(Link* link=from->Next(); link!=from; link=link->Next())	{
		total += GetLogProb(link->Out());
	}
	return total;
}

void PureBrownianProcess::Clamp() {
	Clamp(GetRoot());
}

void PureBrownianProcess::Clamp(Link* from)	{
	if(!from->isRoot())
		this->GetBranchVal(from->GetBranch())->Clamp();
	for(Link* link=from->Next(); link!=from; link=link->Next())	{
		Clamp(link->Out());
	}
}

void PureBrownianProcess::drawSample()	{
	SampleBranch(this->GetRoot());
}

void PureBrownianProcess::SampleBranch(Link* from)	{

	if(!from->isRoot())
		this->GetBranchVal(from->GetBranch())->Sample();
	for(Link* link=from->Next(); link!=from; link=link->Next())	{
		SampleBranch(link->Out());
	}
}


ConjugatePureBrownianProcess::ConjugatePureBrownianProcess(Chronogram *intree, ConjugateInverseWishart *insigma, Var<PosReal>* inagescale, GlobalScalingFunction* inscalefunction)	{

	SetWithRoot(false);
	tree = intree;
	sigma = insigma;
	agescale = inagescale;
	scalefunction = inscalefunction;
	PureBrownianProcess::sigma = insigma;

	RecursiveCreate(GetRoot());

}


Rvar<BrownianBridge> * ConjugatePureBrownianProcess::CreateBranchVal(const Link* link)	{

	if (link->isRoot())	{
		cerr << "error : create branch val called on root\n";
		exit(1);
	}

	// grep the time associated to the branch
	//Var<PosReal>* time = GetLengthTree()->GetBranchLength(link->GetBranch());

	// grep the age associated to the nodes
	Var<PosReal>* up = GetLengthTree()->GetNodeDate(link->GetNode());
	Var<PosReal>* down = GetLengthTree()->GetNodeDate(link->Out()->GetNode());


	// make the new brownian path, and return the pointer
	return new ConjugateRandomBrownianPath(up, down, sigma);
}

