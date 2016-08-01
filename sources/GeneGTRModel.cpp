
#include "BaseType.h"
#include "RandomTypes.h"
#include "ProbModel.h"
#include "IID.h"
#include "OneMatrixPhyloProcess.h"
#include "Move.h"

#include "Normal.h"

// this is a model that runs on several genes in parallel
// for each gene, one has a tree and a multiple alignment
// the evolution of each gene is described by a GTR matrix:
//
// for gene g, and states i and j, the rate of substitution is:
//
// Q_{gij} = \rho_{ij} \pi_{gj}
//
// for g = 1.. Ngene
// and i,j=1..Nstate (Nstate = 4 or 20)
//
// thus, this rate is the combination of
// - a set of global relative rates
// - a gene-specific profile of stationary probabilities (\pi_{gi})_{i=1..Nstate}
//
// the gene-specific profiles are themselves a combination of a global "selection profile" and a gene-specific "selection profile"
// I call them selection profiles, although in reality, they might include mutational aspects as well
// but this is because the ultimate objective is to reformulate this model in a codon framework
//
// the "selection profiles" are iid normal random variables (thus, RealVector variables, of dimension Nstate)
// of mean 0, and variance V0, which is itself endowed with an exponential prior of mean 1.
// positive values mean that the amino-acid or nucleotide is favored, negative values that it is not favored
//
// the global selection profile is:
// (h_{0i})_{i=1..Nstate}
//
// the gene-specific profiles, for g=1..Ngene:
// (f_{gi})_{i=1..Nstate}
//
// and the stationary probability profile for gene g is obtained by
// - taking the sum of the global and gene-specific profiles
// - taking the exponential of that sum
// - renormalizing
//
// \pi_{gi} = 1/Z_g * exp ( h_{0i} + f_{gi} )
//
// Z_g = \sum_{j=1..Nstate} exp( h_{0j} + f_{gj} )
//
// otherwise, each gene tree is endowed with iid gamma (actually, exponentially) distributed branch lengths
//


//
// first, we need a class that takes two Var<RealVector> variables a and b of dimension N (these will be selection profiles)
// and computes a Profile (which will be the stationary probability vector)
//
// pi[i] = 1/Z * exp(a[i] + b[i]), for all i=1..N,
// where Z is such that the c[i] sum to 1 : Z = \sum_i c[i]
//

class SimpleRenormalizedIIDExp : public Dvar<Profile> {

	public:

	SimpleRenormalizedIIDExp(Var<RealVector>* ina)	{
		setval(Profile(ina->GetDim()));
		bkvalue = Profile(ina->GetDim());
		a = ina;
		Register(a);
	}

	// the only function that we need to define is the one that makes the sum,exponential and renormalization
	void specialUpdate()	{

		// compute the exponential of sums
		double total = 0;
		for (int i=0; i<GetDim(); i++)	{
			(*this)[i] = exp( (*a)[i] );
			total += (*this)[i];
		}

		// renormalize
		for (int i=0; i<GetDim(); i++)	{
			(*this)[i] /= total;
		}
	}

	private:
	Var<RealVector>* a;

};
// this class is therefore a Dvar<Profile>
// Dvar because it is a Deterministic function of its parents
// and Profile because it is a vector of positive real numbers that sum to 1
//
// a Dvar should essentially provide
// - a constructor (to register the parents)
// - a specialUpdate function, which specifies how to compute the value of the Dvar given the parents.
// this specialUpdate method is called each time the parents have been corrupted.


class RenormalizedIIDExp : public Dvar<Profile> {

	public:

	// the constructor takes 2 arguments: pointers to tha 2 variables a and b
	// which we want to sum, exponentiate, and renormalized, all this in a componentwise fashion.
	//
	RenormalizedIIDExp(Var<RealVector>* ina, Var<RealVector>* inb)	{
		if (ina->GetDim() != inb->GetDim())	{
			cerr << "error in RenormalizedIIDExp : non matching dimension\n";
			throw;
		}
		setval(Profile(ina->GetDim()));
		bkvalue = Profile(ina->GetDim());
		a = ina;
		b = inb;
		Register(a);
		Register(b);
	}

	// the only function that we need to define is the one that makes the sum,exponential and renormalization
	void specialUpdate()	{

		// compute the exponential of sums
		double total = 0;
		for (int i=0; i<GetDim(); i++)	{
			(*this)[i] = exp( (*a)[i] + (*b)[i] );
			total += (*this)[i];
		}

		// renormalize
		for (int i=0; i<GetDim(); i++)	{
			(*this)[i] /= total;
		}
	}

	private:
	Var<RealVector>* a;
	Var<RealVector>* b;

};


// the model
// all gene-specific variables begin with "Gene"
// all other variables are global variables (same for all genes)
//
class GeneGTRModel : public ProbModel	{

	// number of genes
	int Ngene;

	Tree** GeneTree;
	FileSequenceAlignment** GeneData;
	TaxonSet** GeneTaxonSet;

	int* Nsite;

	// number of states (4 for nucleic acids, 20 for amino-acids)
	// same for all genes
	int Nstate;

	// ---------
	// the random variables of the model
	// ---------

	Const<PosReal>* PriorLambda;
	Exponential* lambda;
	Const<PosReal>* mu;

	// one tree for each gene
	GammaTree** GeneGamTree;

	// relative exchange rates of the matrix
	// (same for all genes)
	Dvar<PosReal>* PriorMeanRelRate;
	IIDExp* relrate;

	// the equilibrium frequencies of the matrix are modelled as
	// proportional to exp(f_g + f_i)
	//
	Dvar<Real>* MeanPropensity;
	Dvar<PosReal>* VarPropensityPrior;
	Exponential* VarPropensity;

	// IIDNormal* BasalSelectionProfile;

	IIDNormal** GeneSelectionProfile;

	SimpleRenormalizedIIDExp** GeneStationary;

	GTRRandomSubMatrix** GeneMatrix;

	Dvar<PosReal>* PriorVarRate;
	Exponential* alpha;
	GammaIIDArray** GeneRate;

	OneMatrixRASPhyloProcess** GenePhyloProcess;

	public:

	// constructor
	// this is where the entire graph structure of the model is created

	GeneGTRModel(int inNgene, string* datafile, string* treefile) {

		// number of genes
		Ngene = inNgene;

		// for each gene,
		// one tree and one dataset (alignment)
		GeneData = new FileSequenceAlignment*[Ngene];
		GeneTree = new Tree*[Ngene];
		GeneTaxonSet = new TaxonSet*[Ngene];
		Nsite = new int[Ngene];

		for (int gene=0; gene<Ngene; gene++)	{
			GeneData[gene] = new FileSequenceAlignment(datafile[gene]);
			GeneTree[gene] = new Tree(treefile[gene]);
			GeneTaxonSet[gene] = GeneData[gene]->GetTaxonSet();
			GeneTree[gene]->RegisterWith(GeneTaxonSet[gene]);
			Nsite[gene] = GeneData[gene]->GetNsite();
		}

		cerr << "tree and data ok\n";

		// fix the number of states
		// should perhaps check that all genes indeed have the same number of states...
		Nstate = GeneData[0]->GetNstate();	// # states (4 for nucleotides, 20 for amino acids)

		// ----------
		// construction of the model
		// ----------

		// branch lengths (for all genes)
		// are iid gamma of shape parameter Mu=1, and scale parameter Lamda ~ Exp(PriorLambda=1)
		PriorLambda = new Const<PosReal>(10);
		mu = new Const<PosReal>(1);
		lambda = new Exponential(PriorLambda,Exponential::MEAN);

		GeneGamTree = new GammaTree*[Ngene];
		for (int gene=0; gene<Ngene; gene++)	{
			GeneGamTree[gene] = new GammaTree(GeneTree[gene],mu,lambda);
		}

		// constructing the general time reversible (GTR) substitution matrices of each gene

		// the global relative rates
		relrate = new IIDExp(Nstate*(Nstate-1)/2);

		// all selection profiles are IIDGamma of mean MeanPropensity=0,
		// and variance VarPropensity ~ Exponential of mean VarPropensityPrior = 1;
		MeanPropensity = new Const<Real>(0);
		VarPropensityPrior = new Const<PosReal>(1);
		VarPropensity = new Exponential(VarPropensityPrior,Exponential::MEAN);

		// the global selection profile for all genes (h_0 above)
		// BasalSelectionProfile = new IIDNormal(Nstate,MeanPropensity,VarPropensity);

		// the gene selection profiles (f_g above)
		GeneSelectionProfile = new IIDNormal*[Ngene];
		for (int gene=0; gene<Ngene; gene++)	{
			GeneSelectionProfile[gene] = new IIDNormal(Nstate,MeanPropensity,VarPropensity);
		}

		// the stationary probabilities of each gene's substitution process
		// here we use the class that we have created above
		GeneStationary = new SimpleRenormalizedIIDExp*[Ngene];
		for (int gene=0; gene<Ngene; gene++)	{
			GeneStationary[gene] = new SimpleRenormalizedIIDExp(GeneSelectionProfile[gene]);
		}

		// gene specific gtr matrices, combining the global relative rates with the gene-specific stationary probability profiles
		GeneMatrix = new GTRRandomSubMatrix*[Ngene];
		for (int gene=0; gene<Ngene; gene++)	{
			GeneMatrix[gene] = new GTRRandomSubMatrix(relrate,GeneStationary[gene]);
		}

		PriorVarRate = new Const<PosReal>(PosReal(1));
		alpha = new Exponential(PriorVarRate,Exponential::MEAN);
		GeneRate = new GammaIIDArray*[Ngene];

		// for each gene, a simple OneMatrixPhyloProcess, with the relevant tree, alignment, and matrix
		GenePhyloProcess = new OneMatrixRASPhyloProcess*[Ngene];
		for (int gene=0; gene<Ngene; gene++)	{
			GeneRate[gene] = new GammaIIDArray(Nsite[gene],alpha,alpha);
			cerr << "create phylo process\n";
			GenePhyloProcess[gene] = new OneMatrixRASPhyloProcess(GeneGamTree[gene],GeneRate[gene],GeneMatrix[gene],GeneData[gene]);
			cerr << "unfold\n";
			GenePhyloProcess[gene]->Unfold();

		}

		// to REGISTER a model: gives all ROOTS of the graph
		// (all nodes that do not have parents)
		// and call Register()
		// Register() will make a recursive enumeration of all nodes of the graph
		// starting from the roots
		cerr << "register\n";
		RootRegister(PriorLambda);
		RootRegister(mu);
		RootRegister(relrate);
		RootRegister(PriorVarRate);
		RootRegister(MeanPropensity);
		RootRegister(VarPropensityPrior);
		Register();

		cerr << "initialise\n";
		Sample();
		Update();
		cerr << "ok\n";
		Trace(cerr);

		cerr << "scheduler\n";
		MakeScheduler();

		cerr << "model created\n";
	}

	// destructor
	// deallocations should normally be done here
	// but in general, the model is deleted just before the program exits, so we can dispense with it for the moment

	~GeneGTRModel() {}

	// log probability of the model is the sum of the log prior and the log likelihood

	double GetLogProb()	{
		return GetLogPrior() + GetLogLikelihood();
	}
	double GetLogPrior()	{
		double total = 0;
		total += lambda->GetLogProb();
		total += alpha->GetLogProb();
		total += relrate->GetLogProb();
		total += VarPropensity->GetLogProb();
		// total += BasalSelectionProfile->GetLogProb();
		for (int gene=0; gene<Ngene; gene++)	{
			total += GeneRate[gene]->GetLogProb();
			total += GeneGamTree[gene]->GetLogProb();
			total += GeneSelectionProfile[gene]->GetLogProb();
		}
		return total;
	}

	double GetLogLikelihood()	{
		double total = 0;
		for (int gene=0; gene<Ngene; gene++)	{
			total += GenePhyloProcess[gene]->GetLogProb();
		}
		return total;
	}

	void MakeScheduler()	{

		// first register the moves on the global variables
		scheduler.Register(new SimpleMove(lambda,1),10,"lambda");
		scheduler.Register(new SimpleMove(lambda,0.1),10,"lambda");

		scheduler.Register(new SimpleMove(relrate,0.1),10,"relrates");
		scheduler.Register(new SimpleMove(relrate,0.01),10,"relrates");

		scheduler.Register(new SimpleMove(alpha,1),10,"alpha");
		scheduler.Register(new SimpleMove(alpha,0.1),10,"alpha");
		scheduler.Register(new SimpleMove(alpha,0.01),10,"alpha");

	/*
		scheduler.Register(new SimpleMove(BasalSelectionProfile,0.1),10,"basal aa");
		scheduler.Register(new SimpleMove(BasalSelectionProfile,0.01),10,"basal aa");
		scheduler.Register(new SimpleMove(BasalSelectionProfile,0.001),10,"basal aa");

		scheduler.Register(new RealVectorMove(BasalSelectionProfile,1,1),10,"basal aa2");
		scheduler.Register(new RealVectorMove(BasalSelectionProfile,0.3,2),10,"basal aa2");
		scheduler.Register(new RealVectorMove(BasalSelectionProfile,0.1,4),10,"basal aa2");
	*/

		scheduler.Register(new SimpleMove(VarPropensity,1),10,"var");
		scheduler.Register(new SimpleMove(VarPropensity,0.1),10,"var");
		scheduler.Register(new SimpleMove(VarPropensity,0.01),10,"var");

		// second register the moves on all gene specific variables
		for (int gene=0; gene<Ngene; gene++)	{

			scheduler.Register(new SimpleMove(GeneGamTree[gene],1),10,"branch lengths");
			scheduler.Register(new SimpleMove(GeneGamTree[gene],0.1),10,"branch lengths");
			scheduler.Register(new SimpleMove(GeneGamTree[gene],0.01),10,"branch lengths");

			scheduler.Register(new SimpleMove(GeneRate[gene],1),10,"rates");
			scheduler.Register(new SimpleMove(GeneRate[gene],0.1),10,"rates");
			scheduler.Register(new SimpleMove(GeneRate[gene],0.01),10,"rates");

			scheduler.Register(new SimpleMove(GeneSelectionProfile[gene],0.1),10,"gene aa");
			scheduler.Register(new SimpleMove(GeneSelectionProfile[gene],0.01),10,"gene aa");
			scheduler.Register(new SimpleMove(GeneSelectionProfile[gene],0.001),10,"gene aa");

			scheduler.Register(new RealVectorMove(GeneSelectionProfile[gene],1,1),10,"gene aa2");
			scheduler.Register(new RealVectorMove(GeneSelectionProfile[gene],0.3,2),10,"gene aa2");
			scheduler.Register(new RealVectorMove(GeneSelectionProfile[gene],0.1,4),10,"gene aa2");


			scheduler.Register(new SimpleMove(GenePhyloProcess[gene],1),1,"mapping");
		}
	/*
		scheduler.Register(new OneToManyRealVectorComponentwiseCompensatoryMove(BasalSelectionProfile,GeneSelectionProfile,Ngene,1),10,"compensate");
		scheduler.Register(new OneToManyRealVectorComponentwiseCompensatoryMove(BasalSelectionProfile,GeneSelectionProfile,Ngene,0.1),10,"compensate");
		scheduler.Register(new OneToManyRealVectorComponentwiseCompensatoryMove(BasalSelectionProfile,GeneSelectionProfile,Ngene,0.01),10,"compensate");
	*/
	}

	// Draw a sample from the prior

	/*
	double Move(double tuning)	{
		scheduler.Cycle(1,1,true,true);
		return 1;
	}
	*/

	void drawSample()	{
		lambda->Sample();
		alpha->Sample();
		relrate->Sample();
		VarPropensity->Sample();
		// BasalSelectionProfile->Sample();
		for (int gene=0; gene<Ngene; gene++)	{
			GeneRate[gene]->Sample();
			GeneGamTree[gene]->Sample();
			GeneSelectionProfile[gene]->Sample();
			GenePhyloProcess[gene]->Sample();
		}
	}

	double GetMeanLength()	{
		double total = 0;
		for (int gene=0; gene<Ngene; gene++)	{
			total += GeneGamTree[gene]->GetTotalLength() * GeneMatrix[gene]->GetRate();
		}
		total /= Ngene;
		return total;
	}

	double GetMeanStationaryEntropy()	{
		double total = 0;
		for (int gene=0; gene<Ngene; gene++)	{
			total += GeneStationary[gene]->GetEntropy();
		}
		total /= Ngene;
		return total;
	}

	// creates the header of the <model_name>.trace file
	void TraceHeader(ostream& os)	{
		os << "#logprior\tlnL\tlength\tlalpha\ttstatent\tmeanrr\trrent\n";
	}

	// writes all summary statistics on one single line
	// in the same order as that provided by the header
	void Trace(ostream& os)	{
	os << GetLogPrior() << '\t' << GetLogLikelihood();
		os << '\t' << GetMeanLength();
		os << '\t' << alpha->val();
		os << '\t' << GetMeanStationaryEntropy();
		os << '\t' << relrate->val().GetMean() << '\t' << relrate->val().GetEntropy();
		os << '\n';
		os.flush();
	}

	void ToStream(ostream& os)	{
		// to be implemented
	}

	void FromStream(istream& is)	{
		// to be implemented
	}

};


int main(int argc, char* argv[])	{

	string datalist = argv[1];
	string treelist = argv[2];
	string name = argv[3];

	ifstream dis(datalist.c_str());
	ifstream tis(treelist.c_str());
	int Ngene;
	dis >> Ngene;
	int Ntree;
	tis >> Ntree;
	if (Ngene != Ntree)	{
		cerr << "error : non matching number of gene alignments and trees\n";
		exit(1);
	}

	string* datafile = new string[Ngene];
	string* treefile = new string[Ngene];
	for (int gene=0; gene<Ngene; gene++)	{
		dis >> datafile[gene];
		tis >> treefile[gene];
	}

	cerr << "create model\n";
	GeneGTRModel* model = new GeneGTRModel(Ngene, datafile, treefile);

	cerr << "start\n";

	ofstream os((name + ".trace").c_str());
	model->TraceHeader(os);

	while (1)	{
		model->Move(1);
		model->Move(0.1);
		model->Trace(os);
		ofstream mon_os((name + ".monitor").c_str());
		model->Monitor(mon_os);
	}
}

