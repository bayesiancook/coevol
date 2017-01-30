#ifndef BROWNIANPROCESS_H
#define	BROWNIANPROCESS_H

#include "MultiVariateTreeProcess.h"
#include "PureBrownianProcess.h"
#include "ConjugateInverseWishart.h"
#include "ConjugateMultiVariateTreeProcess.h"
#include "Move.h"

class BrownianProcess : public MCMC {

	public :
		BrownianProcess() {}
		BrownianProcess(Chronogram *intree, Var<CovMatrix>* insigma, Var<PosReal>* inagescale = 0, GlobalScalingFunction* inscalefunction = 0);
		virtual ~BrownianProcess();

		virtual void drawSample();
		virtual double GetLogProb();

		virtual double Move(double tuning);

		//Accessors
		virtual Tree* GetTree() { return GetLengthTree()->GetTree();}
		virtual Chronogram* GetLengthTree() { return tree; }

		PureBrownianProcess* GetPureBrownianProcess() {return pureBrownianProcess;}
		MultiVariateTreeProcess* GetInstantProcess() {return instantProcess; }

		//Statistics methods

		void PrintLengths(ostream& os);
		void PrintLengths(ostream& os, const Link* from);

		void PrintNodeVals(ostream& os);
		void PrintNodeVals(ostream& os, const Link* from);

		double getTotalLength();
		double getTotalLength(const Link* from);
		double GetMeanRate(int v);
		double GetIntegralRate(int v);
		double GetVarRate(int v);
		double GetTotalRate(int v, const Link* from, int& n);
		double GetTotalSquareRate(int v, const Link* from, int& n);
		double GetIntegralRate(int v, const Link* from);
		double GetMeanGC(int v);
		double GetVarGC(int v);
		double GetTotalGC(int v, const Link* from, int& n);
		double GetTotalSquareGC(int v, const Link* from, int& n);
	  
		double** GetIncrementSum();
		void GetIncrementSum(double** A, const Link* from);
		int GetNBranch();
		int GetNSubBranch();
		int GetNOneSubBranch();
		int GetMinNSubBranch();
		int GetMaxNSubBranch();
		double GetMeanNSubBranch();
		int GetNBranch(const Link* from);
		int GetNSubBranch(const Link* from);
		int GetNOneSubBranch(const Link* from);
		int GetMinNSubBranch(const Link* from);
		int GetMaxNSubBranch(const Link* from);

		protected:

		Chronogram *tree;		//The chronogram of the tree
		Var<CovMatrix> *sigma;		//The covariance matrix
		Var<PosReal>* agescale;
		GlobalScalingFunction* scalefunction;
		
		PureBrownianProcess *pureBrownianProcess;
		MultiVariateTreeProcess *instantProcess;
};

class BrownianSigmaMove : public MCUpdate, public Mnode	{

	public :

	BrownianSigmaMove(InverseWishartMatrix* insigma, BrownianProcess* inprocess, double intuning)	{
		sigma = insigma;
		process = inprocess;
		brownianprocess = inprocess->GetPureBrownianProcess();
		instantprocess = inprocess->GetInstantProcess();
		tuning = intuning;
		RecursiveRegister(GetRoot());
		sigma->Register(this);
	}

	const Link* GetRoot() {
		return instantprocess->GetRoot();
	}

	int GetDim()	{
		return sigma->GetDim();
	}

	void RecursiveRegister(const Link* from)	{

		// instantprocess->GetMultiNormal(from)->Register(this);
		if (! from->isRoot())	{
			brownianprocess->GetRandomBrownianPath(from)->Register(this);
		}

		for(const Link* link=from->Next(); link!=from; link=link->Next())	{
			RecursiveRegister(link->Out());
		}
	}

	void RecursiveMultiply(const Link* from, double** Q)	{

		// instantprocess->GetMultiNormal(from)->LeftMultiply(Q);
		if (! from->isRoot())	{
			brownianprocess->GetRandomBrownianPath(from)->LeftMultiply(Q);
		}

		for(const Link* link=from->Next(); link!=from; link=link->Next())	{
			RecursiveMultiply(link->Out(),Q);
		}
	}

	void ComputeExponential(double** S)	{

		int nstep = 5;

		double** Q = new double*[GetDim()];
		for (int i=0; i<GetDim(); i++)	{
			Q[i] = new double[GetDim()];
		}

		double z = exp(log(2.0) * 2 * nstep);
		for (int i=0; i<GetDim(); i++)	{
			for (int j=0; j<GetDim(); j++)	{
				S[i][j] /= z;
			}
		}
		for (int i=0; i<GetDim(); i++)	{
			S[i][i] += 1.0;
		}

		for (int n=0; n<nstep; n++)	{

			for (int i=0; i<GetDim(); i++)	{
				for (int j=0; j<GetDim(); j++)	{
					double tmp = 0;
					for (int k=0; k<GetDim(); k++)	{
						tmp += S[i][k] * S[k][j];
					}
					Q[i][j] = tmp;
				}
			}

			for (int i=0; i<GetDim(); i++)	{
				for (int j=0; j<GetDim(); j++)	{
					double tmp = 0;
					for (int k=0; k<GetDim(); k++)	{
						tmp += Q[i][k] * Q[k][j];
					}
					S[i][j] = tmp;
				}
			}
		}
		
		for (int i=0; i<GetDim(); i++)	{
			delete[] Q[i];
		}
		delete[] Q;
	}

	double Move(double tuning_mod) {

		double logd1 = sigma->GetLogDeterminant();

		// double logs1 = sigma->GetLogProb();
		// double logq1 = brownianprocess->GetLogProb() + instantprocess->GetLogProb();

		Corrupt(true);

		// draw symmetric matrix
		double** S = new double*[GetDim()];
		for (int i=0; i<GetDim(); i++)	{
			S[i] = new double[GetDim()];
		}

		for (int i=0; i<GetDim(); i++)	{
			for (int j=0; j<GetDim(); j++)	{
				S[i][j] = tuning * tuning_mod * Random::sNormal();
			}
		}
		
		// make exponential
		ComputeExponential(S);

		// apply matrix recursively too all instant processes and all brownian paths
		RecursiveMultiply(GetRoot(),S);

		// apply matrix to sigma
		sigma->LeftRightMultiply(S);

		double deltalogp = Update();

		double logd2 = sigma->GetLogDeterminant();
 
		// double logs2 = sigma->GetLogProb();
		// double logq2 = brownianprocess->GetLogProb() + instantprocess->GetLogProb();

		double logdetp = 0.5 * (logd2 - logd1);
		// + 2 : for sigma
		double logh = (process->GetNSubBranch() - process->GetNBranch() + 2) * logdetp;

		// cerr << tuning << '\t' << logh << '\t' << deltalogp << '\t' << deltalogp + logh << '\n';
		// cerr << tuning << '\t' << (process->GetNSubBranch()) * logdetp + (logq2 - logq1) << '\n';

		deltalogp += logh;

		int accept = (log(Random::Uniform()) < deltalogp);

		if (! accept)	{
			Corrupt(false);
			Restore();
		}

		for (int i=0; i<GetDim(); i++)	{
			delete[] S[i];
		}
		delete[] S;

		return accept;	
		
	}

	private:

	
	PureBrownianProcess* brownianprocess;
	MultiVariateTreeProcess* instantprocess;
	InverseWishartMatrix *sigma;
	BrownianProcess *process;
	double tuning;

};

/*
class BrownianTrinodeMove : public MCUpdate, public Mnode	{

	public :

	BrownianTrinodeMove(LengthTree* intimetree, BrownianProcess* inprocess, double intuning)	{
		process = inprocess;
		brownianprocess = inprocess->GetPureBrownianProcess();
		instantprocess = inprocess->GetInstantProcess();
		tuning = intuning;
	}

	const Link* GetRoot() {
		return instantprocess->GetRoot();
	}

	double Move(double tuning_mod)	{

		double ntot = 0;
		double nacc = 0;

		RecursiveMove(GetRoot(), nacc, ntot, tuning*tuning_mod);

		return nacc / ntot;
	}

	double RecursiveMove(const Link* from, double& naccept, double& ntot, double tuning)	{

		if (! from->isRoot())	{
			naccept += LocalMove(from,tuning);
			ntot++;
		}
		for(const Link* link=from->Next(); link!=from; link=link->Next())	{
			RecursiveMove(link->Out(),naccept,ntot,tuning);
		}
	}

	double LocalMove(const Link* from, double tuning)	{

		Mnode* mnode = new Mnode;

		instantprocess->GetMultiNormal(from)->Register(mnode);
		brownianprocess->GetRandomBrownianPath(from)->Register(mnode);
		double mean = 1.0 / brownianprocess->GetRandomBrownianPath(from)->getLength();
		for(const Link* link=from->Next(); link!=from; link=link->Next())	{
			brownianprocess->GetRandomBrownianPath(link->Out())->Register(mnode);
			mean += 1.0 / brownianprocess->GetRandomBrownianPath(link->Out())->getLength();
		}
		mean = 1.0 / mean;
		
		mnode->Corrupt(true);
		
		instantprocess->GetMultiNormal(from)->MatrixProposeMove();
		// not finished;

		for(const Link* link=from->Next(); link!=from; link=link->Next())	{
		double logratio = mnode->Update();

		int accepted = (log(Random::Uniform()) < logratio);
		if (! accepted)	{
			mnode->Corrupt(false);
			mnode->Restore();
		}
		delete mnode;
		return accepted;
	}

	private:
	
	PureBrownianProcess* brownianprocess;
	MultiVariateTreeProcess* instantprocess;
	BrownianProcess *process;
	double tuning;

};
*/

class ConjugateBrownianSigmaMove : public MCUpdate {

public:
	ConjugateBrownianSigmaMove(ConjugateInverseWishart* insigma, BrownianProcess *inprocess) {
		sigma = insigma;
		process = inprocess;
	}

	double Move(double tuning) {
		double **A = process->GetIncrementSum();
		int N = process->GetNSubBranch();
		sigma->Integrate();
		sigma->AddToScale(A);
		sigma->AddToShape(N);
		sigma->Resample();
		for(int i=0; i<sigma->GetDim(); i++)
			delete[] A[i];
		delete[] A;
		return 1;
	}

private:
	ConjugateInverseWishart *sigma;
	BrownianProcess *process;
};


class ConjugateBrownianProcess : public BrownianProcess	{

	public:

	ConjugateBrownianProcess(Chronogram *intree, ConjugateInverseWishart* insigma);

	ConjugatePureBrownianProcess* GetConjugatePureBrownianProcess() {return conjugatePureBrownianProcess;}
	ConjugateMultiVariateTreeProcess* GetConjugateInstantProcess() {return conjugateInstantProcess; }
	ConjugateInverseWishart* GetConjugateInverseWishart()	{return conjugatesigma;}

	protected:

	ConjugateInverseWishart* conjugatesigma;
	ConjugatePureBrownianProcess *conjugatePureBrownianProcess;
	ConjugateMultiVariateTreeProcess *conjugateInstantProcess;
};


class FullConjugateBrownianSigmaMove : public MCUpdate {

	public:

	FullConjugateBrownianSigmaMove(ConjugateInverseWishart* insigma, ConjugateBrownianProcess *inprocess, double intuning, int innrep) {
		sigma = insigma;
		process = inprocess;
		tuning = intuning;
		nrep = innrep;
	}

	double Move(double tuning_modulator) {

		sigma->Integrate();

		double r = 0;
		for(int i=0; i<nrep; i++){
			r += process->Move(tuning * tuning_modulator);
		}

		sigma->Resample();
		return nrep ? r / nrep : 1.0;
	}

	private:

	ConjugateInverseWishart *sigma;
	ConjugateBrownianProcess *process;
	double tuning;
	int nrep;
};

/*
class ConjugateBrownianExternalMove : public MCUpdate {

	ConjugateBrownianProcess *process;
	ConjugateInverseWishart* truemat;
	InverseWishartMatrix* mat;
	ConjugatePureBrownianProcess* brownianprocess;
	ConjugateMultiVariateTreeProcess* instantprocess;

	double tuning;
	int nrep;
	int dim;

	int shapestat;
	CovMatrix scalestat;
	int bkshapestat;
	CovMatrix bkscalestat;

	public:

	ConjugateBrownianExternalMove(ConjugateInverseWishart* inmat, ConjugateBrownianProcess* inprocess, double intuning, int innrep) : scalestat(inmat->GetDim()), bkscalestat(inmat->GetDim())	{
		process = inprocess;
		brownianprocess = inprocess->GetConjugatePureBrownianProcess();
		instantprocess = inprocess->GetConjugateInstantProcess();
		truemat = inmat;
		tuning = intuning;
		nrep = innrep;
		dim = truemat->GetDim();
	}

	const Link* GetRoot() {
		return instantprocess->GetRoot();
	}

	double Move(double tuning_modulator){

		double naccept = 0;
		double ntot = 0;

		double naccept2 = 0;
		double ntot2 = 0;

		mat = new InverseWishartMatrix(truemat->GetDiagonalMatrix(),truemat->GetP());
		mat->SetAtSigmaZero();
		mat->Corrupt(false);
		mat->Update();

		PrepareSufficientStatistic();
		for (int rep=0; rep<nrep; rep++)	{
			RecursiveNodeConjugateMove(GetRoot(),naccept,ntot,tuning*tuning_modulator);
			RecursiveBranchConjugateMove(GetRoot(),naccept2,ntot2,tuning*tuning_modulator);
		}

		truemat->Integrate();
		truemat->Resample();

		delete mat;

		return naccept / ntot;
	}

	void ResetSufficientStatistic()	{
		shapestat = 0;
		for( int i=0; i<dim; i++){
			for( int j=0; j<dim; j++){
				scalestat[i][j] = 0;
			}
		}
	}

	void SaveSufficientStatistic()	{
		bkshapestat = shapestat;
		for( int i=0; i<dim; i++){
			bkscalestat[i][i] = scalestat[i][i];
			for( int j=0; j<i; j++){
				bkscalestat[i][j] = scalestat[i][j];
			}
		}
	}

	void RestoreSufficientStatistic()	{
		bkshapestat = shapestat;
		for( int i=0; i<dim; i++){
			scalestat[i][i] = bkscalestat[i][i];
			for( int j=0; j<i; j++){
				scalestat[i][j] = bkscalestat[i][j];
				scalestat[j][i] = bkscalestat[i][j];
			}
		}
	}

	void AddToScale(const double* in, bool pos = true)	{
		double epsilon = pos ? 1.0 : -1.0;
		for( int i=0; i<dim; i++){
			scalestat[i][i] += epsilon * in[i] * in[i];
			for( int j=0; j<i; j++){
				scalestat[i][j] += epsilon * in[i] * in[j];
				scalestat[j][i] = scalestat[i][j];
			}
		}
	}

	void AddToScale(double** in, bool pos = true)	{
		double epsilon = pos ? 1.0 : -1.0;
		for( int i=0; i<dim; i++){
			for (int j=0; j<dim; j++)	{
				scalestat[i][j] += epsilon * in[i][j];
			}
		}
	}

	void PrepareSufficientStatistic()	{
		ResetSufficientStatistic();
		RecursiveAddSufficientStatistic(GetRoot());
	}

	void RecursiveAddSufficientStatistic(const Link* from)	{

		if (! from->isRoot())	{

			ConjugateMultiNormal* target = instantprocess->GetConjugateMultiNormal(from);
			target->ComputeContrast();
			AddToScale(target->GetContrast(),true);

			ConjugateRandomBrownianPath* path = brownianprocess->GetConjugateRandomBrownianPath(from->Out());
			path->ComputeContrast();
			AddToScale(path->GetContrast(),true);

			shapestat += path->getNSegments();
		}
		for(const Link* link=from->Next(); link!=from; link=link->Next())	{
			RecursiveAddSufficientStatistic(link->Out());
		}
	}

	double RecursiveNodeConjugateMove(const Link* from, double& naccept, double& ntot, double tuning)	{

		naccept += LocalNodeConjugateMove(from,tuning);
		ntot++;
		for(const Link* link=from->Next(); link!=from; link=link->Next())	{
			RecursiveNodeConjugateMove(link->Out(),naccept,ntot,tuning);
		}
	}

	double RecursiveBranchConjugateMove(const Link* from, double& naccept, double& ntot, double tuning)	{

		if (! from->isRoot())	{
			naccept += LocalBranchConjugateMove(from,tuning);
			ntot++;
		}
		for(const Link* link=from->Next(); link!=from; link=link->Next())	{
			RecursiveBranchConjugateMove(link->Out(),naccept,ntot,tuning);
		}
	}

	double LocalNodeConjugateMove(const Link* from, double tuning)	{

		ConjugateMultiNormal* target = instantprocess->GetConjugateMultiNormal(from);

		// integrated log p(X)
		double logp1 = truemat->MarginalLogProb(scalestat,shapestat);

		// non integrated log p(X | sigma) + log p(sigma)
		double logq1 = 0;
		// double logq1 = mat->GetLogProb();
		logq1 += instantprocess->GetConjugateMultiNormal(from)->GetLogProb();
		for(const Link* link=from->Next(); link!=from; link=link->Next())	{
			logq1 += instantprocess->GetConjugateMultiNormal(link->Out())->GetLogProb();
		}

		SaveSufficientStatistic();

		if (! from->isRoot())	{
			instantprocess->GetConjugateMultiNormal(from)->ComputeContrast();
			AddToScale(instantprocess->GetConjugateMultiNormal(from)->GetContrast(),false);
		}
		for(const Link* link=from->Next(); link!=from; link=link->Next())	{
			instantprocess->GetConjugateMultiNormal(link->Out())->ComputeContrast();
			AddToScale(instantprocess->GetConjugateMultiNormal(link->Out())->GetContrast(),false);
		}


		target->Corrupt(true);

		double* deltax = new double[dim];
		mat->drawVal(deltax);

		// add delta to target
		instantprocess->GetConjugateMultiNormal(from)->Shift(deltax,tuning);

		if (! from->isRoot())	{
			instantprocess->GetConjugateMultiNormal(from)->ComputeContrast();
			AddToScale(instantprocess->GetConjugateMultiNormal(from)->GetContrast(),true);
		}
		for(const Link* link=from->Next(); link!=from; link=link->Next())	{
			instantprocess->GetConjugateMultiNormal(link->Out())->ComputeContrast();
			AddToScale(instantprocess->GetConjugateMultiNormal(link->Out())->GetContrast(),true);
		}


		// integrated log p(X')
		double logp2 = truemat->MarginalLogProb(scalestat,shapestat);

		double deltalogp = target->Update();

		// non integrated log p(X' | sigma') + log p(sigma')
		double logq2 = 0;
		// double logq2 = mat->GetLogProb();
		logq2 += instantprocess->GetConjugateMultiNormal(from)->GetLogProb();
		for(const Link* link=from->Next(); link!=from; link=link->Next())	{
			logq2 += instantprocess->GetConjugateMultiNormal(link->Out())->GetLogProb();
		}

		deltalogp -= logq2 - logq1;
		deltalogp += logp2 - logp1;

		int accept = (log(Random::Uniform()) < deltalogp);
		if (! accept)	{
			RestoreSufficientStatistic();
			target->Corrupt(false);
			target->Restore();
		}

		delete[] deltax;

		return accept;	
	}

	double LocalBranchConjugateMove(const Link* from, double tuning)	{

		ConjugateRandomBrownianPath* target = brownianprocess->GetConjugateRandomBrownianPath(from);

		// integrated log p(X)
		double logp1 = truemat->MarginalLogProb(scalestat,shapestat);

		// non integrated log p(X | sigma) + log p(sigma)
		double logq1 = 0;
		// double logq1 = mat->GetLogProb();
		logq1 += brownianprocess->GetConjugateRandomBrownianPath(from)->GetLogProb();
		for(const Link* link=from->Next(); link!=from; link=link->Next())	{
			logq1 += brownianprocess->GetConjugateRandomBrownianPath(link->Out())->GetLogProb();
		}

		SaveSufficientStatistic();

		if (! from->isRoot())	{
			brownianprocess->GetConjugateRandomBrownianPath(from)->ComputeContrast();
			AddToScale(brownianprocess->GetConjugateRandomBrownianPath(from)->GetContrast(),false);
		}
		for(const Link* link=from->Next(); link!=from; link=link->Next())	{
			brownianprocess->GetConjugateRandomBrownianPath(link->Out())->ComputeContrast();
			AddToScale(brownianprocess->GetConjugateRandomBrownianPath(link->Out())->GetContrast(),false);
		}

		target->Corrupt(true);

		BrownianBridge *other = new BrownianBridge(target->getInitAge(), target->getFinalAge(), *mat , target->getNSegments(), target->getSegmentLength());
		other->generateBridge();
		target->Shift(other,tuning);

		if (! from->isRoot())	{
			brownianprocess->GetConjugateRandomBrownianPath(from)->ComputeContrast();
			AddToScale(brownianprocess->GetConjugateRandomBrownianPath(from)->GetContrast(),true);
		}
		for(const Link* link=from->Next(); link!=from; link=link->Next())	{
			brownianprocess->GetConjugateRandomBrownianPath(link->Out())->ComputeContrast();
			AddToScale(brownianprocess->GetConjugateRandomBrownianPath(link->Out())->GetContrast(),true);
		}

		// integrated log p(X')
		double logp2 = truemat->MarginalLogProb(scalestat,shapestat);

		double deltalogp = target->Update();

		// non integrated log p(X' | sigma') + log p(sigma')
		double logq2 = 0;
		// double logq2 = mat->GetLogProb();
		logq2 += brownianprocess->GetConjugateRandomBrownianPath(from)->GetLogProb();
		for(const Link* link=from->Next(); link!=from; link=link->Next())	{
			logq2 += brownianprocess->GetConjugateRandomBrownianPath(link->Out())->GetLogProb();
		}

		deltalogp -= logq2 - logq1;
		deltalogp += logp2 - logp1;

		int accept = (log(Random::Uniform()) < deltalogp);
		if (! accept)	{
			RestoreSufficientStatistic();
			target->Corrupt(false);
			target->Restore();
		}

		delete other;

		return accept;	
	}
};
*/

class ConjugateBrownianExternalMove : public MCUpdate {

	ConjugateBrownianProcess *process;
	ConjugateInverseWishart* truemat;
	InverseWishartMatrix* mat;
	ConjugatePureBrownianProcess* brownianprocess;
	ConjugateMultiVariateTreeProcess* instantprocess;

	double tuning;
	int nrep;
	int dim;

	int shapestat;
	CovMatrix scalestat;
	int bkshapestat;
	CovMatrix bkscalestat;

	public:

	ConjugateBrownianExternalMove(ConjugateInverseWishart* inmat, ConjugateBrownianProcess* inprocess, double intuning, int innrep) : scalestat(inmat->GetDim()), bkscalestat(inmat->GetDim())	{
		process = inprocess;
		brownianprocess = inprocess->GetConjugatePureBrownianProcess();
		instantprocess = inprocess->GetConjugateInstantProcess();
		truemat = inmat;
		tuning = intuning;
		nrep = innrep;
		dim = truemat->GetDim();
	}

	const Link* GetRoot() {
		return instantprocess->GetRoot();
	}

	double Move(double tuning_modulator){

		double naccept = 0;
		double ntot = 0;

		mat = new InverseWishartMatrix(truemat->GetDiagonalMatrix(),truemat->GetP());
		mat->setval(truemat->val());
		mat->Corrupt(false);
		mat->Update();

		PrepareSufficientStatistic();
		for (int rep=0; rep<nrep; rep++)	{
			RecursiveNodeConjugateMove(GetRoot(),naccept,ntot,tuning*tuning_modulator);
			RecursiveBranchConjugateMove(GetRoot(),naccept,ntot,tuning*tuning_modulator);
		}

		truemat->setval(mat->val());
		truemat->Corrupt(false);
		truemat->Update();
		delete mat;

		return naccept / ntot;
	}

	void ResetSufficientStatistic()	{
		shapestat = 0;
		for( int i=0; i<dim; i++){
			for( int j=0; j<dim; j++){
				scalestat[i][j] = 0;
			}
		}
	}

	void SaveSufficientStatistic()	{
		bkshapestat = shapestat;
		for( int i=0; i<dim; i++){
			bkscalestat[i][i] = scalestat[i][i];
			for( int j=0; j<i; j++){
				bkscalestat[i][j] = scalestat[i][j];
			}
		}
	}

	void RestoreSufficientStatistic()	{
		bkshapestat = shapestat;
		for( int i=0; i<dim; i++){
			scalestat[i][i] = bkscalestat[i][i];
			for( int j=0; j<i; j++){
				scalestat[i][j] = bkscalestat[i][j];
				scalestat[j][i] = bkscalestat[i][j];
			}
		}
	}

	void AddToScale(const double* in, bool pos = true)	{
		double epsilon = pos ? 1.0 : -1.0;
		for( int i=0; i<dim; i++){
			scalestat[i][i] += epsilon * in[i] * in[i];
			for( int j=0; j<i; j++){
				scalestat[i][j] += epsilon * in[i] * in[j];
				scalestat[j][i] = scalestat[i][j];
			}
		}
	}

	void AddToScale(double** in, bool pos = true)	{
		double epsilon = pos ? 1.0 : -1.0;
		for( int i=0; i<dim; i++){
			for (int j=0; j<dim; j++)	{
				scalestat[i][j] += epsilon * in[i][j];
			}
		}
	}

	void PrepareSufficientStatistic()	{
		ResetSufficientStatistic();
		RecursiveAddSufficientStatistic(GetRoot());
	}

	void RecursiveAddSufficientStatistic(const Link* from)	{

		if (! from->isRoot())	{

			ConjugateMultiNormal* target = instantprocess->GetConjugateMultiNormal(from);
			target->ComputeContrast();
			AddToScale(target->GetContrast(),true);

			ConjugateRandomBrownianPath* path = brownianprocess->GetConjugateRandomBrownianPath(from->Out());
			path->ComputeContrast();
			AddToScale(path->GetContrast(),true);

			shapestat += path->getNSegments();
		}
		for(const Link* link=from->Next(); link!=from; link=link->Next())	{
			RecursiveAddSufficientStatistic(link->Out());
		}
	}

	double RecursiveNodeConjugateMove(const Link* from, double& naccept, double& ntot, double tuning)	{

		naccept += LocalNodeConjugateMove(from,tuning);
		ntot++;
		for(const Link* link=from->Next(); link!=from; link=link->Next())	{
			RecursiveNodeConjugateMove(link->Out(),naccept,ntot,tuning);
		}
	}

	double RecursiveBranchConjugateMove(const Link* from, double& naccept, double& ntot, double tuning)	{

		if (! from->isRoot())	{
			naccept += LocalBranchConjugateMove(from,tuning);
			ntot++;
		}
		for(const Link* link=from->Next(); link!=from; link=link->Next())	{
			RecursiveBranchConjugateMove(link->Out(),naccept,ntot,tuning);
		}
	}

	double LocalNodeConjugateMove(const Link* from, double tuning)	{

		ConjugateMultiNormal* target = instantprocess->GetConjugateMultiNormal(from);

		// integrated log p(X)
		double logp1 = truemat->MarginalLogProb(scalestat,shapestat);

		// non integrated log p(X | sigma) + log p(sigma)
		double logq1 = 0;
		// double logq1 = mat->GetLogProb();
		logq1 += instantprocess->GetConjugateMultiNormal(from)->GetLogProb();
		for(const Link* link=from->Next(); link!=from; link=link->Next())	{
			logq1 += instantprocess->GetConjugateMultiNormal(link->Out())->GetLogProb();
		}

		SaveSufficientStatistic();

		if (! from->isRoot())	{
			instantprocess->GetConjugateMultiNormal(from)->ComputeContrast();
			AddToScale(instantprocess->GetConjugateMultiNormal(from)->GetContrast(),false);
		}
		for(const Link* link=from->Next(); link!=from; link=link->Next())	{
			instantprocess->GetConjugateMultiNormal(link->Out())->ComputeContrast();
			AddToScale(instantprocess->GetConjugateMultiNormal(link->Out())->GetContrast(),false);
		}


		target->Corrupt(true);
		mat->Corrupt(true);

		double* deltax = new double[dim];
		mat->drawVal(deltax);
		double logh1 = mat->logValProb(deltax);

		// add delta to target
		instantprocess->GetConjugateMultiNormal(from)->Shift(deltax,tuning);

		if (! from->isRoot())	{
			instantprocess->GetConjugateMultiNormal(from)->ComputeContrast();
			AddToScale(instantprocess->GetConjugateMultiNormal(from)->GetContrast(),true);
		}
		for(const Link* link=from->Next(); link!=from; link=link->Next())	{
			instantprocess->GetConjugateMultiNormal(link->Out())->ComputeContrast();
			AddToScale(instantprocess->GetConjugateMultiNormal(link->Out())->GetContrast(),true);
		}


		// integrated log p(X')
		double logp2 = truemat->MarginalLogProb(scalestat,shapestat);

		mat->GibbsResample(scalestat,shapestat);
		mat->Update();

		double deltalogp = target->Update();

		// logh2 = log p(delta |  sigma')
		double logh2 = mat->logValProb(deltax);
		
		// non integrated log p(X' | sigma') + log p(sigma')
		double logq2 = 0;
		// double logq2 = mat->GetLogProb();
		logq2 += instantprocess->GetConjugateMultiNormal(from)->GetLogProb();
		for(const Link* link=from->Next(); link!=from; link=link->Next())	{
			logq2 += instantprocess->GetConjugateMultiNormal(link->Out())->GetLogProb();
		}

		deltalogp -= logq2 - logq1;
		deltalogp += logp2 - logp1;

		int accept = (log(Random::Uniform()) < deltalogp);
		if (! accept)	{
			RestoreSufficientStatistic();
			target->Corrupt(false);
			target->Restore();
			mat->Corrupt(false);
			mat->Restore();
		}

		delete[] deltax;

		return accept;	
	}

	double LocalBranchConjugateMove(const Link* from, double tuning)	{

		ConjugateRandomBrownianPath* target = brownianprocess->GetConjugateRandomBrownianPath(from);

		// integrated log p(X)
		double logp1 = truemat->MarginalLogProb(scalestat,shapestat);

		// non integrated log p(X | sigma) + log p(sigma)
		double logq1 = 0;
		// double logq1 = mat->GetLogProb();
		logq1 += brownianprocess->GetConjugateRandomBrownianPath(from)->GetLogProb();
		for(const Link* link=from->Next(); link!=from; link=link->Next())	{
			logq1 += brownianprocess->GetConjugateRandomBrownianPath(link->Out())->GetLogProb();
		}

		SaveSufficientStatistic();

		if (! from->isRoot())	{
			brownianprocess->GetConjugateRandomBrownianPath(from)->ComputeContrast();
			AddToScale(brownianprocess->GetConjugateRandomBrownianPath(from)->GetContrast(),false);
		}
		for(const Link* link=from->Next(); link!=from; link=link->Next())	{
			brownianprocess->GetConjugateRandomBrownianPath(link->Out())->ComputeContrast();
			AddToScale(brownianprocess->GetConjugateRandomBrownianPath(link->Out())->GetContrast(),false);
		}

		target->Corrupt(true);
		mat->Corrupt(true);

		BrownianBridge *other = new BrownianBridge(target->getInitAge(), target->getFinalAge(), *mat , target->getNSegments(), target->getSegmentLength());
		other->generateBridge();
		target->Shift(other,tuning);

		double logh1 = other->getLogProb();

		if (! from->isRoot())	{
			brownianprocess->GetConjugateRandomBrownianPath(from)->ComputeContrast();
			AddToScale(brownianprocess->GetConjugateRandomBrownianPath(from)->GetContrast(),true);
		}
		for(const Link* link=from->Next(); link!=from; link=link->Next())	{
			brownianprocess->GetConjugateRandomBrownianPath(link->Out())->ComputeContrast();
			AddToScale(brownianprocess->GetConjugateRandomBrownianPath(link->Out())->GetContrast(),true);
		}

		// integrated log p(X')
		double logp2 = truemat->MarginalLogProb(scalestat,shapestat);

		mat->GibbsResample(scalestat,shapestat);
		mat->Update();

		double deltalogp = target->Update();

		// logh2 = log p(delta |  sigma')
		double logh2 = other->getLogProb(*mat);
		
		// non integrated log p(X' | sigma') + log p(sigma')
		double logq2 = 0;
		// double logq2 = mat->GetLogProb();
		logq2 += brownianprocess->GetConjugateRandomBrownianPath(from)->GetLogProb();
		for(const Link* link=from->Next(); link!=from; link=link->Next())	{
			logq2 += brownianprocess->GetConjugateRandomBrownianPath(link->Out())->GetLogProb();
		}

		deltalogp -= logq2 - logq1;
		deltalogp += logp2 - logp1;

		int accept = (log(Random::Uniform()) < deltalogp);
		if (! accept)	{
			RestoreSufficientStatistic();
			target->Corrupt(false);
			target->Restore();
			mat->Corrupt(false);
			mat->Restore();
		}

		delete other;

		return accept;	
	}
};

class ConjugateBrownianExternalAllBranchMove : public MCUpdate {

	ConjugateBrownianProcess *process;
	ConjugateInverseWishart* truemat;
	InverseWishartMatrix* mat;
	ConjugatePureBrownianProcess* brownianprocess;
	ConjugateMultiVariateTreeProcess* instantprocess;

	double tuning;
	int nrep;
	int dim;

	int shapestat;
	CovMatrix scalestat;
	int bkshapestat;
	CovMatrix bkscalestat;

	Mnode* mnode;

	public:

	ConjugateBrownianExternalAllBranchMove(ConjugateInverseWishart* inmat, ConjugateBrownianProcess* inprocess, double intuning, int innrep) : scalestat(inmat->GetDim()), bkscalestat(inmat->GetDim())	{
		process = inprocess;
		brownianprocess = inprocess->GetConjugatePureBrownianProcess();
		instantprocess = inprocess->GetConjugateInstantProcess();
		truemat = inmat;
		tuning = intuning;
		nrep = innrep;
		dim = truemat->GetDim();
		mnode = new Mnode;
		RecursiveRegister(GetRoot());
	}

	const Link* GetRoot() {
		return instantprocess->GetRoot();
	}

	void RecursiveRegister(const Link* from)	{

		if (! from->isRoot())	{
			brownianprocess->GetConjugateRandomBrownianPath(from)->Register(mnode);
		}

		for(const Link* link=from->Next(); link!=from; link=link->Next())	{
			RecursiveRegister(link->Out());
		}
	}

	double Move(double tuning_modulator){

		double naccept = 0;
		double ntot = 0;


		PrepareBranchSufficientStatistic();
		mat = new InverseWishartMatrix(truemat->GetDiagonalMatrix(),truemat->GetP());
		mat->Corrupt(false);
		mat->Update();
		mat->GibbsResample(scalestat,shapestat);

		PrepareSufficientStatistic();
		for (int rep=0; rep<nrep; rep++)	{
			naccept += AllBranchConjugateMove(tuning*tuning_modulator);
			ntot++;
		}

		truemat->Integrate();
		truemat->Resample();

		delete mat;

		return naccept / ntot;
	}

	void ResetSufficientStatistic()	{
		shapestat = 0;
		for( int i=0; i<dim; i++){
			for( int j=0; j<dim; j++){
				scalestat[i][j] = 0;
			}
		}
	}

	void SaveSufficientStatistic()	{
		bkshapestat = shapestat;
		for( int i=0; i<dim; i++){
			bkscalestat[i][i] = scalestat[i][i];
			for( int j=0; j<i; j++){
				bkscalestat[i][j] = scalestat[i][j];
			}
		}
	}

	void RestoreSufficientStatistic()	{
		bkshapestat = shapestat;
		for( int i=0; i<dim; i++){
			scalestat[i][i] = bkscalestat[i][i];
			for( int j=0; j<i; j++){
				scalestat[i][j] = bkscalestat[i][j];
				scalestat[j][i] = bkscalestat[i][j];
			}
		}
	}

	void AddToScale(const double* in, bool pos = true)	{
		double epsilon = pos ? 1.0 : -1.0;
		for( int i=0; i<dim; i++){
			scalestat[i][i] += epsilon * in[i] * in[i];
			for( int j=0; j<i; j++){
				scalestat[i][j] += epsilon * in[i] * in[j];
				scalestat[j][i] = scalestat[i][j];
			}
		}
	}

	void AddToScale(double** in, bool pos = true)	{
		double epsilon = pos ? 1.0 : -1.0;
		for( int i=0; i<dim; i++){
			for (int j=0; j<dim; j++)	{
				scalestat[i][j] += epsilon * in[i][j];
			}
		}
	}

	void PrepareSufficientStatistic()	{
		ResetSufficientStatistic();
		RecursiveAddSufficientStatistic(GetRoot());
	}

	void RecursiveAddSufficientStatistic(const Link* from)	{

		if (! from->isRoot())	{

			ConjugateMultiNormal* target = instantprocess->GetConjugateMultiNormal(from);
			target->ComputeContrast();
			AddToScale(target->GetContrast(),true);

			ConjugateRandomBrownianPath* path = brownianprocess->GetConjugateRandomBrownianPath(from->Out());
			path->ComputeContrast();
			AddToScale(path->GetContrast(),true);

			shapestat += path->getNSegments();
		}
		for(const Link* link=from->Next(); link!=from; link=link->Next())	{
			RecursiveAddSufficientStatistic(link->Out());
		}
	}

	void PrepareBranchSufficientStatistic()	{
		ResetSufficientStatistic();
		RecursiveAddBranchSufficientStatistic(GetRoot());
	}

	void RecursiveAddBranchSufficientStatistic(const Link* from)	{

		if (! from->isRoot())	{

			ConjugateMultiNormal* target = instantprocess->GetConjugateMultiNormal(from);
			target->ComputeContrast();
			AddToScale(target->GetContrast(),true);

			shapestat ++;
		}
		for(const Link* link=from->Next(); link!=from; link=link->Next())	{
			RecursiveAddBranchSufficientStatistic(link->Out());
		}
	}

	double AllBranchConjugateMove(double tuning)	{

		double logp1 = truemat->MarginalLogProb(scalestat,shapestat);
		double logq1 = brownianprocess->GetLogProb();

		SaveSufficientStatistic();

		mnode->Corrupt(true);

		brownianprocess->RecursiveProposeMove(GetRoot(),tuning,*mat);

		PrepareSufficientStatistic();
		
		// integrated log p(X')
		double logp2 = truemat->MarginalLogProb(scalestat,shapestat);
		double logq2 = brownianprocess->GetLogProb();

		double deltalogp = mnode->Update();
		deltalogp -= logq2 - logq1;
		deltalogp += logp2 - logp1;

		int accept = (log(Random::Uniform()) < deltalogp);
		if (! accept)	{
			RestoreSufficientStatistic();
			mnode->Corrupt(false);
			mnode->Restore();
		}

		return accept;	
	}

};


#endif	/* BROWNIANPROCESS_H */

