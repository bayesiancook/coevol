
#include "BaseType.h"
#include "RandomTypes.h"
#include "ProbModel.h"
#include "IID.h"
#include "OneMatrixPhyloProcess.h"
#include "Move.h"

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
	Dirichlet** GeneStationary;

	GTRRandomSubMatrix** GeneMatrix;

	Dvar<PosReal>* PriorVarRate;
	Exponential* alpha;
	GammaIIDArray** GeneRate;

	OneMatrixPhyloProcess** GenePhyloProcess;

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

		GeneStationary = new Dirichlet*[Ngene];
		for (int gene=0; gene<Ngene; gene++)	{
			GeneStationary[gene] = new Dirichlet(Nstate);
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
		GenePhyloProcess = new OneMatrixPhyloProcess*[Ngene];
		for (int gene=0; gene<Ngene; gene++)	{
			GeneRate[gene] = new GammaIIDArray(Nsite[gene],alpha,alpha);
			cerr << "create phylo process\n";
			// GenePhyloProcess[gene] = new OneMatrixPhyloProcess(GeneGamTree[gene],GeneMatrix[gene],GeneData[gene]);
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
		RootRegister(PriorVarRate);
		RootRegister(relrate);
		for (int gene=0; gene<Ngene; gene++)	{
			RootRegister(GeneStationary[gene]);
		}
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
		// total += BasalSelectionProfile->GetLogProb();
		for (int gene=0; gene<Ngene; gene++)	{
			total += GeneRate[gene]->GetLogProb();
			total += GeneStationary[gene]->GetLogProb();
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

		// second register the moves on all gene specific variables
		for (int gene=0; gene<Ngene; gene++)	{

			scheduler.Register(new SimpleMove(GeneGamTree[gene],1),10,"branch lengths");
			scheduler.Register(new SimpleMove(GeneGamTree[gene],0.1),10,"branch lengths");
			scheduler.Register(new SimpleMove(GeneGamTree[gene],0.01),10,"branch lengths");

			scheduler.Register(new SimpleMove(GeneRate[gene],1),10,"rates");
			scheduler.Register(new SimpleMove(GeneRate[gene],0.1),10,"rates");
			scheduler.Register(new SimpleMove(GeneRate[gene],0.01),10,"rates");

			scheduler.Register(new ProfileMove(GeneStationary[gene],0.3,1),10,"gene stat2");
			scheduler.Register(new ProfileMove(GeneStationary[gene],0.1,2),10,"gene stat4");
			scheduler.Register(new ProfileMove(GeneStationary[gene],0.1,4),10,"gene stat8");

			scheduler.Register(new SimpleMove(GenePhyloProcess[gene],1),1,"mapping");
		}
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
		for (int gene=0; gene<Ngene; gene++)	{
			GeneRate[gene]->Sample();
			GeneGamTree[gene]->Sample();
			GeneStationary[gene]->Sample();
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
		os << "#logprior\tlnL\tlength\talpha\ttstatent\tmeanrr\trrent\n";
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

