#ifndef CONJUGATEGTRMODEL_H
#define CONJUGATEGTRMODEL_H

#include "ProbModel.h"
#include "ConjugatePath.h"
#include "ConjugateGammaTree.h"
#include "ConjugateOneMatrixPhyloProcess.h"

#include "IID.h"


//
// A Conjugate Model uses Conjugate Gibbs Sampling methods such as described in Lartillot, 2006, J Comput Biol, 13:1701
//
// it can be conjugate with respect to branch lengths, site specific rates and relative exchange rates of the matrix
// it is also "semi-conjugate" w.r.t. the stationary probabilities of the susbtitution matrix
// (in the loose sense that it implements Metropolis based on data augmentation and sufficient statistic, but not Gibbs sampling)
//


class ConjugateGTRModel : public ProbModel	{

	public:

	ConjugateGTRModel(string datafile, string treefile);

	~ConjugateGTRModel();

	double GetLogPrior();
	double GetLogLikelihood();
	double GetLogProb();

	void MakeScheduler();

	double Move(double tuning = 1);

	void drawSample();

	double GetMeanRate();
	double GetVarRate();
	double GetLength();

	void TraceHeader(ostream& os);
	void Trace(ostream& os);

	void PrintTree(ostream& os);
	void TreeToPS(string filename);

	void ToStream(ostream& os);
	void FromStream(istream& is);

	void PostMapping(ostream& os, bool red);
	void PostMapping(int site, ostream& os, bool red);
	void PostPredMapping(ostream& os, bool red);
	void PostPredMapping(int site, ostream& os, bool red);

	ConjugateGammaTree* gamtree;

	private:
	Tree* tree;
	SequenceAlignment* data;
	TaxonSet* taxonset;

	int Nsite;
	int Nstate;

	// rates across sites
	Dvar<PosReal>* PriorVarRate;
	Exponential* alpha;
	ConjugateGammaIIDArray* rate;

	// branch lengths
	Dvar<PosReal>* PriorLambda;
	Dvar<PosReal>* mu;
	Exponential* lambda;

	// relative exchange rates of the matrix
	IIDConjugateGamma* relrate;

	// equilibrium frequencies of the matrix
	SemiConjugateDirichlet* stationary;

	GTRRandomSubMatrix* matrix;

	ConjugateOneMatrixRASPhyloProcess* phyloprocess;

};

#endif // CONJUGATEGTRMODEL_H
