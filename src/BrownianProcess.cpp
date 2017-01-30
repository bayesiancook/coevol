#include "BrownianProcess.h"

BrownianProcess::BrownianProcess(Chronogram *intree, Var<CovMatrix>* insigma, Var<PosReal>* inagescale, GlobalScalingFunction* inscalefunction) {

	tree = intree;
	sigma = insigma;
	agescale = inagescale;
	scalefunction = inscalefunction;

	pureBrownianProcess = new PureBrownianProcess(tree, sigma,agescale,scalefunction);
	instantProcess = new MultiVariateTreeProcess(sigma,tree,agescale,scalefunction);
	instantProcess->Reset();

}


BrownianProcess::~BrownianProcess() {
	delete pureBrownianProcess;
	delete instantProcess;
}

void BrownianProcess::drawSample() {
	pureBrownianProcess->Sample();
	instantProcess->Sample();
}

double BrownianProcess::GetLogProb() {
	return pureBrownianProcess->GetLogProb() + instantProcess->GetLogProb();
}

double BrownianProcess::Move(double tuning) {
	double t1 = pureBrownianProcess->Move(tuning);
	double t2 = instantProcess->Move(tuning);
	return t1+t2/2;
}

double BrownianProcess::getTotalLength() {
	return getTotalLength(GetTree()->GetRoot());
}

double BrownianProcess::getTotalLength(const Link* from) {
	double r = 0;
	if(!from->isRoot()) {
		r = pureBrownianProcess->GetRandomBrownianPath(from)->getLength();
	}

	for(Link* link=from->Next(); link!=from; link=link->Next())
		r+=getTotalLength(link->Out());

	return r;
}

void BrownianProcess::PrintLengths(ostream& os)	{
	PrintLengths(os,GetTree()->GetRoot());
}

void BrownianProcess::PrintLengths(ostream& os, const Link* from)	{

	os << GetTree()->GetLeftMost(from) << '\t' << GetTree()->GetRightMost(from) << '\t';
	os << instantProcess->GetNodeVal(from->GetNode())->val() << '\t';
	if (! from->isRoot())	{
		os << pureBrownianProcess->GetRandomBrownianPath(from)->getLength() << '\n';
	}
	else	{
		os << 0 << '\n';
	}
	for(Link* link=from->Next(); link!=from; link=link->Next())	{
		PrintLengths(os,link->Out());
	}
}

void BrownianProcess::PrintNodeVals(ostream& os)	{
	PrintNodeVals(os,GetTree()->GetRoot());
}

void BrownianProcess::PrintNodeVals(ostream& os, const Link* from)	{

	os << GetTree()->GetLeftMost(from) << '\t' << GetTree()->GetRightMost(from);
	Rvar<RealVector>* tmp = instantProcess->GetNodeVal(from->GetNode());
	for (int i=0; i<tmp->GetDim(); i++)	{
		os << '\t' << (*tmp)[i];
	}
	os << '\n';
	for(Link* link=from->Next(); link!=from; link=link->Next())	{
		PrintNodeVals(os,link->Out());
	}
}

double** BrownianProcess::GetIncrementSum() {
	int dim = sigma->GetDim();

	double** A = new double*[dim];
	for(int i=0; i<dim;i++) {
		A[i] = new double[dim];
		for(int j=0; j<sigma->GetDim(); j++)
			A[i][j] = 0;
	}

	GetIncrementSum(A, GetTree()->GetRoot());

	return A;

}

void BrownianProcess::GetIncrementSum(double** A, const Link* from) {
	int dim = sigma->GetDim();

	for(const Link* link = from->Next(); link!=from; link = link->Next()) {
		GetIncrementSum(A, link->Out());

		const RealVector initValue = instantProcess->GetNodeVal(link->GetNode())->val();
		const RealVector finalValue = instantProcess->GetNodeVal(link->Out()->GetNode())->val();
		BrownianBridge* bb = pureBrownianProcess->GetRandomBrownianPath(link);
			  
		int N = bb->getNSegments();
		double* dValue = new double[dim];
		for(int i=0; i<dim; i++)
			// scalefunction 
			dValue[i] = (finalValue[i]-initValue[i])/bb->getLength();		//Part of the transition associated to the ramp

		for(int k=0; k<N; k++) {
			double* y = new double[dim];
			double segmentLength = bb->getSegmentLength(k);
			for(int i=0; i<dim; i++)
				y[i] = dValue[i]*segmentLength + bb->getValue(k+1, i) - bb->getValue(k, i);
			for(int i=0; i<dim; i++)
				for(int j=0; j<dim; j++)
					A[i][j] += y[i]*y[j]/segmentLength;
			delete[] y;
		}
		delete[] dValue;
	}
		  
}

double BrownianProcess::GetMeanNSubBranch() {
	return ((double) GetNSubBranch(GetTree()->GetRoot())) / GetNBranch(GetTree()->GetRoot()) ;
}

int BrownianProcess::GetNBranch()	{
	return GetNBranch(GetTree()->GetRoot());
}

int BrownianProcess::GetNBranch(const Link* from) {
	int i=0;
	for(const Link* link = from->Next(); link!=from; link = link->Next()) {
		i++;
		i+=GetNBranch(link->Out());
	}
	return i;
}

int BrownianProcess::GetMinNSubBranch() {
	return GetMinNSubBranch(GetTree()->GetRoot());
}

int BrownianProcess::GetMinNSubBranch(const Link* from) {
	int min = -1;
	for(const Link* link = from->Next(); link!=from; link = link->Next()) {
		int tmp = pureBrownianProcess->GetRandomBrownianPath(link)->getNSegments();
		if ((min == -1) || (min > tmp))	{
			min = tmp;
		}
		tmp = GetMinNSubBranch(link->Out());
		if ((min == -1) || (min > tmp))	{
			if (tmp != -1)	{
				min = tmp;
			}
		}
	}
	return min;
}

int BrownianProcess::GetNOneSubBranch() {
	return GetNOneSubBranch(GetTree()->GetRoot());
}

int BrownianProcess::GetNOneSubBranch(const Link* from) {
	int tot = 0;
	for(const Link* link = from->Next(); link!=from; link = link->Next()) {
		int tmp = pureBrownianProcess->GetRandomBrownianPath(link)->getNSegments();
		if (tmp == 1)	{
			tot++;
		}
		tot += GetMinNSubBranch(link->Out());
	}
	return tot;
}

int BrownianProcess::GetMaxNSubBranch() {
	return GetMaxNSubBranch(GetTree()->GetRoot());
}

int BrownianProcess::GetMaxNSubBranch(const Link* from) {
	int max = 0;
	for(const Link* link = from->Next(); link!=from; link = link->Next()) {
		int tmp = pureBrownianProcess->GetRandomBrownianPath(link)->getNSegments();
		if (max < tmp)	{
			max = tmp;
		}
		tmp = GetMaxNSubBranch(link->Out());
		if (max < tmp)	{
			max = tmp;
		}
	}
	return max;
}

int BrownianProcess::GetNSubBranch() {

	return GetNSubBranch(GetTree()->GetRoot());
}

int BrownianProcess::GetNSubBranch(const Link* from) {
	int i=0;
	for(const Link* link = from->Next(); link!=from; link = link->Next()) {
		i+= pureBrownianProcess->GetRandomBrownianPath(link)->getNSegments();
		i+=GetNSubBranch(link->Out());
	}
	return i;
}

 double BrownianProcess::GetMeanRate(int v){
	int n = 0;
	double total = GetTotalRate(v, GetTree()->GetRoot(),n);
	return total / n;
}

  double BrownianProcess::GetIntegralRate(int v){

	return GetIntegralRate(v, GetTree()->GetRoot());
}

 double BrownianProcess::GetVarRate(int v){
	int n = 0;
	double total1 = GetTotalRate(v, GetTree()->GetRoot(),n);
	n = 0;
	double total2 = GetTotalSquareRate(v, GetTree()->GetRoot(),n);
	total1 /= n;
	total2 /= n;
	total2 -= total1 * total1;
	return total2;
}

 double BrownianProcess::GetTotalRate(int v, const Link* from, int& n)	{
	double total = exp(instantProcess->GetNodeVal(from->GetNode())->val().GetArray()[v]);
	n++;
	for (const Link* link=from->Next(); link!=from; link=link->Next())	{
			total += GetTotalRate(v, link->Out(),n);
	}
	return total;
}

double BrownianProcess::GetTotalSquareRate(int v, const Link* from, int& n)	{
	double temp =  exp(instantProcess->GetNodeVal(from->GetNode())->val().GetArray()[v]);
	double total = temp * temp;
	n++;
	for (const Link* link=from->Next(); link!=from; link=link->Next())	{
			total += GetTotalSquareRate(v, link->Out(),n);
	}
	return total;
}

 double BrownianProcess::GetIntegralRate(int v, const Link* from)	{
	 double total = 0;
	for (const Link* link=from->Next(); link!=from; link=link->Next())	{
			total += GetIntegralRate(v, link->Out());
			double init = instantProcess->GetNodeVal(link->GetNode())->val().GetArray()[v];
			double final = instantProcess->GetNodeVal(link->Out()->GetNode())->val().GetArray()[v];
			total += pureBrownianProcess->GetBranchVal(link->GetBranch())->getIntegral(v, init, final);
	}
	return total;
}



 double BrownianProcess::GetMeanGC(int v){
	int n = 0;
	double total = GetTotalGC(v, GetTree()->GetRoot(),n);
	return total / n;
}

 double BrownianProcess::GetVarGC(int v){
	int n = 0;
	double total1 = GetTotalGC(v, GetTree()->GetRoot(),n);
	n = 0;
	double total2 = GetTotalSquareGC(v, GetTree()->GetRoot(),n);
	total1 /= n;
	total2 /= n;
	total2 -= total1 * total1;
	return total2;
}

 double BrownianProcess::GetTotalGC(int v, const Link* from, int& n)	{
	double total = 1.0/(1.0+exp(-instantProcess->GetNodeVal(from->GetNode())->val().GetArray()[v]));
	n++;
	for (const Link* link=from->Next(); link!=from; link=link->Next())	{
			total += GetTotalGC(v, link->Out(),n);
	}
	return total;
}

double BrownianProcess::GetTotalSquareGC(int v, const Link* from, int& n)	{
	double temp =  1.0/(1.0+exp(-instantProcess->GetNodeVal(from->GetNode())->val().GetArray()[v]));
	double total = temp * temp;
	n++;
	for (const Link* link=from->Next(); link!=from; link=link->Next())	{
			total += GetTotalSquareGC(v, link->Out(),n);
	}
	return total;
}


ConjugateBrownianProcess::ConjugateBrownianProcess(Chronogram *intree, ConjugateInverseWishart* insigma) {

	tree = intree;
	sigma = insigma;
	conjugatesigma = insigma;

	conjugatePureBrownianProcess = new ConjugatePureBrownianProcess(tree, conjugatesigma);
	pureBrownianProcess = conjugatePureBrownianProcess;

	conjugateInstantProcess = new ConjugateMultiVariateTreeProcess(conjugatesigma,tree);
	instantProcess = conjugateInstantProcess;
	instantProcess->Reset();

}
