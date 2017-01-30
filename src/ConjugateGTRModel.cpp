
#include "DrawTree.h"
#include "ConjugateGTRModel.h"
#include "IID.h"

ConjugateGTRModel::ConjugateGTRModel(string datafile, string treefile)	{

	// fetch data from file
	data = new FileSequenceAlignment(datafile);
	Nsite = data->GetNsite();	// # columns
	Nstate = data->GetNstate();	// # states (20 for amino acids)

	taxonset = data->GetTaxonSet();

	// get tree from file (newick format)
	tree = new Tree(treefile);

	// check whether tree and data fits together
	tree->RegisterWith(taxonset);

	cerr << "tree and data ok\n";

	PriorLambda = new Const<PosReal>(10);
	mu = new Const<PosReal>(1);
	lambda = new Exponential(PriorLambda,Exponential::MEAN);
	gamtree = new ConjugateGammaTree(tree,mu,lambda);

	PriorVarRate = new Const<PosReal>(1);
	alpha = new Exponential(PriorVarRate,Exponential::MEAN);
	rate = new ConjugateGammaIIDArray(Nsite,alpha,alpha);

	relrate = new IIDConjugateGamma(Nstate*(Nstate-1)/2);
	stationary = new SemiConjugateDirichlet(Nstate);
	matrix = new GTRRandomSubMatrix(relrate, stationary);

	/*
	alpha->ClampAt(1);
	for (int i=0; i<Nsite; i++)	{
		(*rate)[i]->ClampAt(1);
	}
	relrate->setall(1);
	relrate->Clamp();
	stationary->setuniform();
	stationary->Clamp();
	*/

	cerr << "new phyloprocess\n";
	phyloprocess = new ConjugateOneMatrixRASPhyloProcess(gamtree, rate, matrix, data);

	cerr << "unfold\n";
	phyloprocess->Unfold();

	cerr << "register\n";
	RootRegister(PriorVarRate);
	RootRegister(PriorLambda);
	RootRegister(mu);
	RootRegister(relrate);
	RootRegister(stationary);
	Register();

		cerr << "initialise\n";
		Sample();
		Trace(cerr);

	cerr << "model created\n";
}

ConjugateGTRModel::~ConjugateGTRModel()	{

}

double ConjugateGTRModel::GetLogPrior()	{
	double total = 0;
	total += alpha->GetLogProb();
	total += rate->GetLogProb();
	total += lambda->GetLogProb();
	total += gamtree->GetLogProb();
	total += stationary->GetLogProb();
	total += relrate->GetLogProb();
	return total;
}

double ConjugateGTRModel::GetLogLikelihood()	{
	return phyloprocess->GetLogProb();
}

double ConjugateGTRModel::GetLogProb()	{
	return GetLogPrior() + GetLogLikelihood();
}

double ConjugateGTRModel::GetMeanRate()	{
	double total = 0;
	for (int i=0; i<Nsite; i++)	{
		total += (*rate)[i]->val();
	}
	total /= Nsite;
	return total;
}

double ConjugateGTRModel::GetVarRate()	{
	double mean = 0;
	double var = 0;
	for (int i=0; i<Nsite; i++)	{
		double tmp = (*rate)[i]->val();
		mean += tmp;
		var += tmp * tmp;
	}
	mean /= Nsite;
	var /= Nsite;
	var -= mean * mean;
	return var;
}

double ConjugateGTRModel::GetLength()	{
	return gamtree->GetTotalLength() * phyloprocess->GetMatrix()->GetRate() * GetMeanRate();
}

void ConjugateGTRModel::TraceHeader(ostream& os)	{
	os << "#logprior\tlnL\tlength\tlambda\talpha\tmeanrate\tvarrate\tstatent\tmeanrr\trrent\n";
}

void ConjugateGTRModel::Trace(ostream& os)	{
	os << GetLogPrior() << '\t' << GetLogLikelihood();
	os << '\t' << lambda->val();
	os << '\t' << alpha->val();
	os << '\t' << GetMeanRate();
	os << '\t' << GetVarRate();
	os << '\t' << stationary->val().GetEntropy();
	os << '\t' << relrate->val().GetMean() << '\t' << relrate->val().GetEntropy();
	os << '\n';
}

void ConjugateGTRModel::MakeScheduler()	{
	cerr << "scheduler is inactivated\n";
}

double ConjugateGTRModel::Move(double tuning)	{

	for (int i=0; i<Nsite; i++)	{
		(*rate)[i]->Integrate();
	}
	for (int j=0; j<10; j++)	{
		alpha->Move(tuning);
	}
	for (int i=0; i<Nsite; i++)	{
		(*rate)[i]->Resample();
	}
	// phyloprocess->Move(tuning);

	gamtree->Integrate();
	for (int j=0; j<10; j++)	{
		lambda->Move(tuning);
	}
	gamtree->Resample();
	// phyloprocess->Move(tuning);

	relrate->Integrate();
	relrate->Resample();
	// phyloprocess->Move(tuning);

	stationary->ActivateSufficientStatistic();

	for (int j=0; j<100; j++)	{
		stationary->Move(tuning,1);
		stationary->Move(tuning/10,1);
		stationary->Move(tuning/50,2);
		stationary->Move(0.1*tuning);
		stationary->Move(0.01*tuning);
		stationary->Move(0.001*tuning);

	}
	stationary->InactivateSufficientStatistic();

	phyloprocess->Move(tuning);

	return 1;
}

void ConjugateGTRModel::drawSample()	{
	// draw rates
	alpha->Sample();
	rate->Sample();

	// draw gamma tree
	lambda->Sample();
	gamtree->Sample();

	// draw matrix
	stationary->Sample();
	relrate->Sample();

	// make simulation along tree
	phyloprocess->Sample();
}

void ConjugateGTRModel::ToStream(ostream& os)	{

	os << *alpha << '\n';
	os << *rate << '\n';
	os << *lambda << '\n';
	os << *gamtree << '\n';
	os << *stationary << '\n';
	os << *relrate << '\n';
}

void ConjugateGTRModel::FromStream(istream& is)	{

	is >> *alpha;
	is >> *rate;
	is >> *lambda;
	is >> *gamtree;
	is >> *stationary;
	is >> *relrate;
}



