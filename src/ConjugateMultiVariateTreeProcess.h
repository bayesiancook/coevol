#ifndef CONJUGATEMULTIVARIATETREEPROCESS
#define CONJUGATEMULTIVARIATETREEPROCESS

#include "MultiVariateTreeProcess.h"
#include "ConjugateInverseWishart.h"
#include "ContinuousData.h"
#include "AutoRegressiveMultiVariateTreeProcess.h"


class ConjugateMultiVariateTreeProcess : public virtual MultiVariateTreeProcess	{

	protected:

	ConjugateInverseWishart* sigma;

	public:

	ConjugateMultiVariateTreeProcess() {}

	ConjugateMultiVariateTreeProcess(ConjugateInverseWishart* insigma, NodeBranchVarTree<PosReal,PosReal>* inchrono, Var<PosReal>* inagescale, GlobalScalingFunction* inscalefunction, Var<RealVector>* inrootmean = 0, Var<PosRealVector>* inrootvar = 0){

		sigma = insigma;
		MultiVariateTreeProcess::sigma = insigma;
		tree = inchrono;
		chrono = inchrono;
		agescale = inagescale;
		scalefunction = inscalefunction;
		scaletree = 0;
		drift = 0;
		rootmean = inrootmean;
		rootvar = inrootvar;
		RecursiveCreate(GetRoot());
	}

	ConjugateMultiVariateTreeProcess(ConjugateInverseWishart* insigma, LengthTree* intree, LengthTree* inscaletree = 0, Var<RealVector>* indrift = 0, Var<RealVector>* inrootmean = 0, Var<PosRealVector>* inrootvar = 0){

		sigma = insigma;
		MultiVariateTreeProcess::sigma = insigma;
		tree = intree;
		scaletree = inscaletree;
		drift = indrift;
		rootmean = inrootmean;
		rootvar = inrootvar;
		scalefunction = 0;
		RecursiveCreate(GetRoot());
	}

	~ConjugateMultiVariateTreeProcess(){
		RecursiveDelete(GetRoot());
	}


	ConjugateMultiNormal* GetConjugateMultiNormal(const Link* link)	{
		ConjugateMultiNormal* m = dynamic_cast<ConjugateMultiNormal*> (GetNodeVal(link->GetNode()));
		return m;
	}

	void DAgostinosKandZ(double* KK, double* ZZ1, double* ZZ2)	{

		double* m1 = new double[GetDim()];
		double* m2 = new double[GetDim()];
		double* m3 = new double[GetDim()];
		double* m4 = new double[GetDim()];

		double n = GetContrastMean(m1);
		GetContrastCentralMoment(m2,m1,2);
		GetContrastCentralMoment(m3,m1,3);
		GetContrastCentralMoment(m4,m1,4);

		for (int i=0; i<GetDim(); i++)	{
			double g1 = m3[i] / sqrt(m2[i] * m2[i] * m2[i]);
			double g2 = m4[i] / m2[i] / m2[i] - 3;

			double mu2 = 6.0 * (n-2) / (n+1) / (n+3);
			double gamma2 = 36.0 * (n-7) * (n*n + 2*n -5) / (n-2) / (n+5) / (n+7) / (n+9);
			double W2 = sqrt(2*gamma2 +4)- 1;
			double W = sqrt(W2);
			double delta = 1.0 / sqrt(log(W));
			double alpha = sqrt(2 / (W2 - 1));
			double t = g1 / alpha / sqrt(mu2);
			double Z1 = delta * log(t + sqrt(t*t+1));

			double mmu1 = -6.0 / (n+1);
			double mmu2 = 24.0 * n * (n-2) * (n-3) / (n+1) / (n+1) / (n+3) / (n+5);
			double ggamma1 = 6.0 * (n*n -5*n + 2) / (n+7) / (n+9) * sqrt(6.0 * (n+3) * (n+5) / n / (n-2) / (n-3));
			double A = 6.0 + 8.0 / ggamma1* (2.0 / ggamma1 + sqrt(1 + 4.0 / ggamma1 / ggamma1));
			double Z2 = sqrt(9.0 * A / 2) * (1 - 2.0 / 9 / A - exp( 1.0 / 3 * log( (1 - 2.0 / A) / (1.0 + (g2 - mmu1) / sqrt(mmu2) * sqrt(2.0 / (A - 4))))));

			ZZ1[i] = Z1 * Z1;
			ZZ2[i] = Z2 * Z2;
			KK[i] = ZZ1[i] + ZZ2[i];
		}

		delete[] m1;
		delete[] m2;
		delete[] m3;
		delete[] m4;
	}

	int GetContrastMean(double* m)	{
		for (int i=0; i<GetDim(); i++)	{
			m[i] = 0;
		}
		int n = 0;
		RecursiveGetContrastMean(GetRoot(),m,n);
		for (int i=0; i<GetDim(); i++)	{
			m[i] /= n;
		}
		return n;
	}

	int GetContrastCentralMoment(double* mk, double* m1, int k)	{
		for (int i=0; i<GetDim(); i++)	{
			mk[i] = 0;
		}
		int n = 0;
		RecursiveGetContrastCentralMoment(GetRoot(),mk,m1,k,n);
		for (int i=0; i<GetDim(); i++)	{
			mk[i] /= n;
		}
		return n;
	}

	void RecursiveGetContrastMean(const Link* from, double* m, int& n)	{
		if (! from->isRoot())	{
			double* c = GetConjugateMultiNormal(from)->GetContrast();
			for (int i=0; i<GetDim(); i++)	{
				m[i] += c[i];
			}
			n++;
		}
		for (const Link* link=from->Next(); link!=from; link=link->Next())	{
			RecursiveGetContrastMean(link->Out(),m,n);
		}
	}

	void RecursiveGetContrastCentralMoment(const Link* from, double* mk, const double* m1, int k, int& n)	{
		if (! from->isRoot())	{
			double* c = GetConjugateMultiNormal(from)->GetContrast();
			for (int i=0; i<GetDim(); i++)	{
				mk[i] += pow((c[i] - m1[i]),k);
			}
			n++;
		}
		for (const Link* link=from->Next(); link!=from; link=link->Next())	{
			RecursiveGetContrastCentralMoment(link->Out(),mk,m1,k,n);
		}
	}


	protected :

	Rvar<RealVector>* CreateNodeVal(const Link* link){
		if(!link->isRoot()){
			if (scalefunction)	{
				return new ConjugateMultiNormal(sigma, GetMultiNormal(link->Out()), tree->GetBranchVal(link->GetBranch()), chrono->GetNodeVal(link->GetNode()), agescale, scalefunction);
			}
			return new ConjugateMultiNormal(sigma, GetMultiNormal(link->Out()), tree->GetBranchVal(link->GetBranch()), (scaletree ? scaletree->GetBranchVal(link->GetBranch()) : 0), drift);
		}
		else{
 			return new ConjugateMultiNormal(sigma);
		}
	}

};

class ConjugateMultiVariateMove : public MCUpdate {

	ConjugateInverseWishart* sigma;
	ConjugateMultiVariateTreeProcess* tree;
	VarArray<PosReal>* priorOnSigmaZero;
	int nmove;
	double tuning;

	public:

	ConjugateMultiVariateMove(ConjugateInverseWishart* insigma, ConjugateMultiVariateTreeProcess* intree, double intuning, int innmove){
	// ConjugateMultiVariateMove(ConjugateInverseWishart* insigma, ConjugateMultiVariateTreeProcess* intree, RvarVec* inpriorOnSigmaZero, int innmove){
		sigma = insigma;
		tree = intree;
		tuning = intuning;
		nmove =innmove;
		// priorOnSigmaZero = inpriorOnSigmaZero;
	}

	double Move(double tuning_modulator){
		sigma->Integrate();
		// sigma->ActivateSufficientStatistic();
		double r = 0;
		// double p = 0;
		for(int i=0; i<nmove; i++){
			r += tree->Move(tuning * tuning_modulator);
			/*
			for(int j=0; j<10; j++){
				p += priorOnSigmaZero->Move(tuning);
			}
			*/
		}
		sigma->Resample();
		// sigma->InactivateSufficientStatistic();
		return nmove ? r/nmove : 1.0;
	}
};

class ConjugateInverseWishartPriorMove : public MCUpdate {

	ConjugateInverseWishart* sigma;
	RandomVarArray<PosReal>* priorOnSigmaZero;
	int nmove;
	double tuning;

	public:

	ConjugateInverseWishartPriorMove(ConjugateInverseWishart* insigma, RandomVarArray<PosReal>* inpriorOnSigmaZero, double intuning, int innmove){
		sigma = insigma;
		priorOnSigmaZero = inpriorOnSigmaZero;
		tuning = intuning;
		nmove =innmove;
	}

	double Move(double tuning_modulator){
		sigma->Integrate();
		// sigma->ActivateSufficientStatistic();
		double r = 0;
		// double p = 0;
		for(int i=0; i<nmove; i++){
			r += priorOnSigmaZero->Move(tuning * tuning_modulator);
		}
		sigma->Resample();
		// sigma->InactivateSufficientStatistic();
		return r/nmove;
	}
};

class ConjugateMultiVariatePiecewiseTranslationMove : public MCUpdate {

	ConjugateInverseWishart* sigma;
	ConjugateMultiVariateTreeProcess* tree;
	RandomVarArray<PosReal>* priorOnSigmaZero;
	int nmove;
	double tuning;
	int index;
	int k;

	public:

	ConjugateMultiVariatePiecewiseTranslationMove(ConjugateInverseWishart* insigma, ConjugateMultiVariateTreeProcess* intree, double intuning, int inindex, int ink,int innmove){
	// ConjugateMultiVariateMove(ConjugateInverseWishart* insigma, ConjugateMultiVariateTreeProcess* intree, RvarVec* inpriorOnSigmaZero, int innmove){
		sigma = insigma;
		tree = intree;
		tuning = intuning;
		nmove =innmove;
		index = inindex;
		k = ink;
		// priorOnSigmaZero = inpriorOnSigmaZero;
	}

	double Move(double tuning_modulator){
		sigma->Integrate();
		// sigma->ActivateSufficientStatistic();
		double r = 0;
		// double p = 0;
		for(int i=0; i<nmove; i++){
			r += tree->PiecewiseTranslationMove(tuning * tuning_modulator,index,k);
			/*
			for(int j=0; j<10; j++){
				p += priorOnSigmaZero->Move(tuning);
			}
			*/
		}
		sigma->Resample();
		// sigma->InactivateSufficientStatistic();
		return r/nmove;
	}
};

class ConjugateMultiVariateWholeTreePiecewiseTranslationMove : public MCUpdate, public Mnode {

	ConjugateInverseWishart* sigma;
	ConjugateMultiVariateTreeProcess* tree;
	RandomVarArray<PosReal>* priorOnSigmaZero;
	int nmove;
	double tuning;
	int index;
	int k;

	public:

	ConjugateMultiVariateWholeTreePiecewiseTranslationMove(ConjugateInverseWishart* insigma, ConjugateMultiVariateTreeProcess* intree, double intuning, int inindex, int ink,int innmove){
		sigma = insigma;
		tree = intree;
		tuning = intuning;
		nmove =innmove;
		index = inindex;
		k = ink;
		tree->RecursiveRegister(this,tree->GetRoot());
	}

	double Move(double tuning_modulator){
		double tot = 0;
		sigma->Integrate();
		for(int i=0; i<nmove; i++){
			Corrupt(true);
			double u = tuning * tuning_modulator * (Random::Uniform() - 0.5);
			tree->PiecewiseTranslation(u,index,k);
			double logratio = Update();
			bool accepted = (log(Random::Uniform()) < logratio);
			if (! accepted)	{
				Corrupt(false);
				Restore();
			}
			else 	{
				tot++;
			}
		}
		sigma->Resample();
		// sigma->FullCorrupt();
		// sigma->FullUpdate();
		return tot/nmove;
	}
};



class ConjugateAutoRegressiveMultiNormal : public ConjugateSampling<RealVector>, public AutoRegressiveMultiNormal	{

	public:

	ConjugateAutoRegressiveMultiNormal(MultiNormalSemiConjugate* insigma, Var<RealVector>* inmean, Var<PosReal>* inphi, Var<RealVector>* inup = 0, Var<PosReal>* intime = 0) :
		Rvar<RealVector>(),
		ConjugateSampling<RealVector>(),
		AutoRegressiveMultiNormal(insigma, inmean, inphi, inup, intime)
	{
		conjugate_up.insert(insigma);
		contrast = 0;
	}

	~ConjugateAutoRegressiveMultiNormal() {}

	void AddSufficientStatistic(SemiConjPrior* parent)	{
		MultiNormalSemiConjugate* prior = dynamic_cast<MultiNormalSemiConjugate*>(parent);
		if (! prior)	{
			cout << "cast error in ConjugateMultiNormal::AddSuffStat\n";
			exit(1);
		}
		prior->AddToShape();
		ComputeContrast();
		prior->AddToScale(contrast);
	}

	void ComputeContrast()	{
		if (! contrast)	{
			contrast = new double[GetDim()];
		}
		if(isRoot()){
			double mu = 1.0 / sqrt(2 * phi->val());
			for ( int i=0; i <GetDim(); i++){
				contrast[i] = ((*this)[i] -(*mean)[i]) / mu;
			}
		}
		else{
			double expo = exp(-phi->val() * time->val());
			double mu = sqrt ((1 - exp(-2 * phi->val() * time->val())) / 2 / phi->val());
			for ( int i=0; i <GetDim(); i++){
				contrast[i] = ((*this)[i] - (*up)[i] * expo - (*mean)[i] * (1 - expo)) / mu;
			}
		}
	}

	private:
	double* contrast;

};

class ConjugateAutoRegressiveMultiVariateTreeProcess : public AutoRegressiveMultiVariateTreeProcess, ConjugateMultiVariateTreeProcess	{

	public:

	ConjugateAutoRegressiveMultiVariateTreeProcess(ConjugateInverseWishart* insigma, Var<RealVector>* inmean, Var<PosReal>* inphi, LengthTree* intree){
		sigma = insigma;
		phi = inphi;
		mean = inmean;
		AutoRegressiveMultiVariateTreeProcess::sigma = insigma;
		tree = intree;
		RecursiveCreate(GetRoot());
	}

	~ConjugateAutoRegressiveMultiVariateTreeProcess(){
		RecursiveDelete(GetRoot());
	}


	protected :

	Rvar<RealVector>* CreateNodeVal(const Link* link){
		if(!link->isRoot()){
			return new ConjugateAutoRegressiveMultiNormal(sigma, mean, phi, GetMultiNormal(link->Out()), tree->GetBranchVal(link->GetBranch()));
		}
		else{
 			return new ConjugateAutoRegressiveMultiNormal(sigma, mean, phi);
		}
	}

};


class ConjugateMultiVariateExternalMove : public MCUpdate {

	ConjugateMultiVariateTreeProcess* tree;
	ConjugateInverseWishart* truemat;
	InverseWishartMatrix* mat;
	double tuning;
	int nrep;
	int dim;

	int shapestat;
	CovMatrix scalestat;
	int bkshapestat;
	CovMatrix bkscalestat;

	public:

	ConjugateMultiVariateExternalMove(ConjugateInverseWishart* inmat, ConjugateMultiVariateTreeProcess* intree, double intuning, int innrep) : scalestat(inmat->GetDim()), bkscalestat(inmat->GetDim())	{
		tree = intree;
		truemat = inmat;
		tuning = intuning;
		nrep = innrep;
		dim = truemat->GetDim();
	}

	const Link* GetRoot() {
		return tree->GetRoot();
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
			RecursiveConjugateMove(GetRoot(),naccept,ntot,tuning*tuning_modulator);
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

	void PrepareSufficientStatistic()	{
		ResetSufficientStatistic();
		RecursiveAddSufficientStatistic(GetRoot());
	}

	void RecursiveAddSufficientStatistic(const Link* from)	{

		if (! from->isRoot())	{
			ConjugateMultiNormal* target = tree->GetConjugateMultiNormal(from);
			target->ComputeContrast();
			AddToScale(target->GetContrast());
			shapestat++;
		}
		for(const Link* link=from->Next(); link!=from; link=link->Next())	{
			RecursiveAddSufficientStatistic(link->Out());
		}
	}

	double RecursiveConjugateMove(const Link* from, double& naccept, double& ntot, double tuning)	{

		naccept += LocalConjugateMove(from,tuning);
		ntot++;
		for(const Link* link=from->Next(); link!=from; link=link->Next())	{
			RecursiveConjugateMove(link->Out(),naccept,ntot,tuning);
		}
	}

	double LocalConjugateMove(const Link* from, double tuning)	{

		ConjugateMultiNormal* target = tree->GetConjugateMultiNormal(from);

		// integrated log p(X)
		double logp1 = truemat->MarginalLogProb(scalestat,shapestat);

		// non integrated log p(X | sigma) + log p(sigma)
		double logq1 = 0;
		// double logq1 = mat->GetLogProb();
		logq1 += tree->GetConjugateMultiNormal(from)->GetLogProb();
		for(const Link* link=from->Next(); link!=from; link=link->Next())	{
			logq1 += tree->GetConjugateMultiNormal(link->Out())->GetLogProb();
		}

		SaveSufficientStatistic();
		if (! from->isRoot())	{
			tree->GetConjugateMultiNormal(from)->ComputeContrast();
			AddToScale(tree->GetConjugateMultiNormal(from)->GetContrast(),false);
		}
		for(const Link* link=from->Next(); link!=from; link=link->Next())	{
			tree->GetConjugateMultiNormal(link->Out())->ComputeContrast();
			AddToScale(tree->GetConjugateMultiNormal(link->Out())->GetContrast(),false);
		}

		target->Corrupt(true);
		mat->Corrupt(true);

		double* deltax = new double[dim];
		mat->drawVal(deltax);
		double logh1 = mat->logValProb(deltax);

		// add delta to target
		tree->GetConjugateMultiNormal(from)->Shift(deltax,tuning);

		if (! from->isRoot())	{
			tree->GetConjugateMultiNormal(from)->ComputeContrast();
			AddToScale(tree->GetConjugateMultiNormal(from)->GetContrast(),true);
		}
		for(const Link* link=from->Next(); link!=from; link=link->Next())	{
			tree->GetConjugateMultiNormal(link->Out())->ComputeContrast();
			AddToScale(tree->GetConjugateMultiNormal(link->Out())->GetContrast(),true);
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
		logq2 += tree->GetConjugateMultiNormal(from)->GetLogProb();
		for(const Link* link=from->Next(); link!=from; link=link->Next())	{
			logq2 += tree->GetConjugateMultiNormal(link->Out())->GetLogProb();
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
};

#endif

