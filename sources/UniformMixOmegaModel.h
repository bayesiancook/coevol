
#ifndef MIXOMEGA
#define MIXOMEGA

#include "MeanValTree.h"

#include "BaseType.h"
#include "RandomTypes.h"
#include "ProbModel.h"
#include "IID.h"
#include "PrecisionNormalTreeProcess.h"
#include "AutoRegressiveNormalTreeProcess.h"
#include "CodonSequenceAlignment.h"
#include "BDCalibratedChronogram.h"
// #include "CoalCalibratedChronogram.h"

#include "BranchProcess.h"
#include "NodeProcess.h"
#include "MatrixTree.h"
#include "OneMatrixPhyloProcess.h"
#include "BranchMatrixPhyloProcess.h"
#include "ContinuousData.h"
#include "MeanExpTree.h"
#include "Normal.h"
#include "MultiVarNormal.h"

#include "AutoRegressiveMultiVariateTreeProcess.h"

#include "GCProcess.h"

#include "WhiteNoise.h"

#include "GeneralConjugatePath.h"
#include "CodonConjugatePath.h"

#include "Jeffreys.h"

#include "Partition.h"
#include "SplitPartition.h"

#include "SplitLengthTree.h"
#include "SplitChronogram.h"
#include "SplitMultiVariateMove.h"
#include "MultiVariatePropagateMove.h"

#include "PartitionMultiVariateTreeProcess.h"
#include "KalmanMultiVariateTreeProcess.h"

#include "LinRegOmega.h"

#include "Parallel.h"

const double minjeff = 0.001;
const double maxjeff = 1000;
const double minuni = -1000;
const double maxuni = 1000;

class MixOmegaModel : public ProbModel {

	public:

	// data fields

	// ---------
	// the fixed parameters of the model
	// ---------

	// a fixed tree (read from file)
	Tree* tree;
	FileSequenceAlignment** nucdata;
	CodonSequenceAlignment** codondata;
	string* genename;
	int Ngene;
	int MaxNgene;
	int* procngene;
	int* genesize;
	int* genealloc;
	double* genelnl;
	double* tmpgenelnl;
	int* procnsite;
	int totnsite;

	map<const Link*,int>* missingmap;

	ContinuousData* contdata;
	TaxonSet* taxonset;
	int Ncont;
	CodonStateSpace* codonstatespace;

	// ---------
	// the random variables of the model
	// ---------

	Const<Real>* Zero;
	Const<PosReal>* One;

	Const<PosReal>* MeanChi;
	Const<PosReal>* MeanChi2;
	Exponential* Chi;
	Exponential* Chi2;

	string calibfile;
	string omegafile;

	int chronoprior;
	double meanchi;
	double meanchi2;
	// 0 : uniform;
	// 1 : bd;
	// 2 : bd with cauchy proper lower bounds
	// 3 : rbd with cauchy proper lower bounds

	bool iscalib;

	// chronogram
	Chronogram* chronogram;

	JeffreysIIDArray* DiagArray;
	SigmaZero* sigmaZero;
	Rvar<CovMatrix>* sigma;

	Const<RealVector>* rootmean;
	Const<PosRealVector>* rootvar;

	MultiVariateTreeProcess* process;

	Jeffreys* syntaushape;
	Jeffreys* syntauscale;
	GammaIIDArray* genesyntau;

	Uniform* rootsynmean;
	Jeffreys* rootsynvar;
	NormalIIDArray* generootsyn;

	BidimIIDUniform* synregarray;

	LinRegNormalProcess** genesynprocess;
	MeanExpTree** genesyntree;

	Jeffreys* taushape;
	Jeffreys* tauscale;
	GammaIIDArray* genetau;

	Uniform* rootomegamean;
	Jeffreys* rootomegavar;
	NormalIIDArray* generootomega;

	BidimIIDUniform* regarray;

	LinRegNormalProcess** geneprocess;
	MeanExpTree** geneomegatree;

	Dirichlet* relratecenter;
	Dirichlet* stationarycenter;
	Jeffreys* relrateconcentration;
	Jeffreys* stationaryconcentration;
	DirichletIIDArray* generelrate;
	DirichletIIDArray* genestationary;

	GTRRandomSubMatrixWithNormRates** nucmatrix;

	MatrixTree** genecodonmatrixtree;

	// phylo process
	PathConjugateTree** pathconjtree;
	PhyloProcess** phyloprocess;

	bool meanexp;

	bool conjpath;
	bool normalise;

	GeneticCodeType type;

	int nrep;

	int df;

	int clampsuffstat;
	string suffstatfile;

	bool uniform;
	bool clamptree;
	bool clampdiag;
	bool clampreg;
	bool clamptoggle;
	bool contml;

	bool priorsampling;

	double mappingfreq;
	int mappingevery;

	bool fixomega;
	int myid;
	int nprocs;
	Mnode* mpinode;
	int mpiparamdim;
	double* mpiparamarray;
	int mpiprocessdim;
	double* mpiprocessarray;

	public:

	CodonSequenceAlignment* GetData(int gene)	{
		return codondata[gene];
	}

	ConjugateInverseWishart* GetConjugateInverseWishart() {
		ConjugateInverseWishart* tmp = dynamic_cast<ConjugateInverseWishart*>(sigma);
		if (! tmp)	{
			cerr << "error : dynamic cast of conjugate inverse wishart : " << sigma << '\t' << tmp << '\n';
			exit(1);
		}
		return tmp;
	}

	ConjugateMultiVariateTreeProcess* GetConjugateMultiVariateTreeProcess() {
		ConjugateMultiVariateTreeProcess* tmp = dynamic_cast<ConjugateMultiVariateTreeProcess*>(process);
		if (! tmp)	{
			cerr << "error : dynamic cast of multivariate tree process : " << process << '\t' << tmp << '\n';
			exit(1);
		}
		return tmp;
	}

	ConditionalExternalKalmanMultiVariateTreeProcess* GetConditionalExternalKalmanMultiVariateTreeProcess() {
		ConditionalExternalKalmanMultiVariateTreeProcess* tmp = dynamic_cast<ConditionalExternalKalmanMultiVariateTreeProcess*>(process);
		if (! tmp)	{
			cerr << "error : dynamic cast of multivariate tree process into kalman : " << process << '\t' << tmp << '\n';
			exit(1);
		}
		return tmp;
	}

	Link* GetRoot() {return process->GetRoot();}

	bool isCalibrated()	{
		return iscalib;
	}

	bool isConstrained()	{
		return clamptree;
	}

	MixOmegaModel(string datafile, string treefile, string contdatafile, string incalibfile, double rootage, double rootstdev, int inchronoprior, double inmeanchi, double inmeanchi2, int indf, int inclampdiag, int inconjpath, int inclamptree, int inmeanexp, int innormalise, int innrep, string insuffstatfile, string rootfile, int inclampreg, double inmappingfreq, int inmappingevery, string omegafile, GeneticCodeType intype, int inmyid, int innprocs, bool sample=true)	{

		myid = inmyid;
		nprocs = innprocs;

		if (omegafile != "None")	{
			fixomega = true;
		}
		else	{
			fixomega = false;
		}

		calibfile = incalibfile;
		if (inclampreg == 2)	{
			clampreg = true;
		}
		else	{
			clampreg = false;
		}
		mappingfreq = inmappingfreq;
		mappingevery = inmappingevery;
		type = intype;

		clamptree = inclamptree;

		suffstatfile = insuffstatfile;
		clampsuffstat = (suffstatfile != "None");

		df = indf;

		chronoprior = inchronoprior;
		iscalib = false;
		meanchi = inmeanchi;
		meanchi2 = inmeanchi2;

		if (inclampdiag == 2)	{
			clampdiag = false;
			contml = true;
		}
		else if (inclampdiag == 1)	{
			clampdiag = true;
			contml = false;
		}
		else	{
			clampdiag = false;
			contml = false;
		}

		meanexp = inmeanexp;

		priorsampling = false;

		if (inconjpath == -1)	{
			conjpath = true;
		}
		else if (inconjpath == 2)	{
			conjpath = false;
			priorsampling = true;
		}
		else	{
			conjpath = inconjpath;
		}
		nrep = innrep;
		/*
		if (nrep == 0)	{
			nrep = conjpath ? 10 : 1;
		}
		*/
		normalise = innormalise;
		if (clampsuffstat && ! conjpath)	{
			cerr << "error : suffstat requires conjugate path sampling\n";
			exit(1);
		}

		AllocateAlignments(datafile);

		taxonset = nucdata[0]->GetTaxonSet();
		codonstatespace = codondata[0]->GetCodonStateSpace();

		// get tree from file (newick format)
		tree = new Tree(treefile);
		// check whether tree and data fit together
		tree->RegisterWith(taxonset);

		CreateMissingMap(tree);

		// get continuous data from file
		if (contdatafile != "None")	{
			contdata = new FileContinuousData(contdatafile);
			Ncont = contdata->GetNsite();
		}
		else	{
			cerr << "error : life history model requires continuous data\n";
			exit(1);
			contdata = 0;
			Ncont = 0;
		}

		if (!myid)	{
			MasterCreate();
			CreateMPIParamArrays();
			CreateMPIGeneProcessesArrays();
			MasterRegisterMPINode();
			MasterMakeScheduler();
			if (sample)	{
				Update();
				TraceHeader(cerr);
				Trace(cerr);
			}
		}
		else	{
			SlaveCreate(sample);
			CreateMPIParamArrays();
			CreateMPIGeneProcessesArrays();
			SlaveRegisterMPINode();
			SlaveMakeScheduler();
			if (sample)	{
				Update();
				SlaveTrace();
			}
		}

	}

	// destructor
	// deallocations should normally be done here
	// but in general, the model is deleted just before the program exits, so we can dispense with it for the moment
	~MixOmegaModel() {}

	Tree* GetTree() {return tree;}

	// BranchVarTree<UnitReal>* GetGCTree() {return gcprocess;}

	MultiVariateTreeProcess* GetProcess() {return process;}

	int GetNcont()	{
		return Ncont;
	}

	int GetNgene()	{
		return Ngene;
	}

	string GetGeneName(int gene)	{
		return genename[gene];
	}

	Chronogram* GetChronogram() {
		return chronogram;
	}

	LengthTree* GetLengthTree()	{
		return chronogram;
	}


	ContinuousData* GetContinuousData() {return contdata;}

	CovMatrix* GetSigma() {return sigma;}

	CalibratedChronogram* GetCalibratedChronogram()	{
		return dynamic_cast<CalibratedChronogram*>(chronogram);
	}

	double Update(bool check = false)	{
		double ret = ProbModel::Update();
		/*
		if (phyloprocess)	{
			phyloprocess->Sample();
			ret = ProbModel::Update();
		}
		*/
		return ret;
	}

	double GetLogProb()	{
		if (myid)	{
			cerr << "error: slave in get log prob\n";
			exit(1);
		}
		return GetLogPrior() + GetLogLikelihood();
	}

	double GetLogPrior()	{
		if (myid)	{
			cerr << "error: slave in get log prior\n";
			exit(1);
		}
		double total = 0;

		if (! clamptree)	{
			if ((chronoprior >= 1) && (chronoprior <= 3))	{
				total += Chi->GetLogProb();
				total += Chi2->GetLogProb();
			}
			total += chronogram->GetLogProb();
		}

		total += DiagArray->GetLogProb();
		total += sigma->GetLogProb();
		total += process->GetLogProb();

		if (! clampreg)	{
			total += synregarray->GetLogProb();
			total += regarray->GetLogProb();
		}

		total += syntaushape->GetLogProb();
		total += syntauscale->GetLogProb();
		total += genesyntau->GetLogProb();

		if (! fixomega)	{
			total += rootsynmean->GetLogProb();
			total += rootsynvar->GetLogProb();
			total += generootsyn->GetLogProb();
		}

		for (int gene=0; gene<Ngene; gene++)	{
			total += genesynprocess[gene]->GetLogProb();
		}

		total += taushape->GetLogProb();
		total += tauscale->GetLogProb();
		total += genetau->GetLogProb();

		if (! fixomega)	{
			total += rootomegamean->GetLogProb();
			total += rootomegavar->GetLogProb();
			total += generootomega->GetLogProb();
		}

		for (int gene=0; gene<Ngene; gene++)	{
			total += geneprocess[gene]->GetLogProb();
		}

		if (! fixomega)	{
			total += relratecenter->GetLogProb();
			total += stationarycenter->GetLogProb();
			total += relrateconcentration->GetLogProb();
			total += stationaryconcentration->GetLogProb();
			total += generelrate->GetLogProb();
			total += genestationary->GetLogProb();
		}

		return total;
	}

	double GetLogLikelihood()	{
		if (myid)	{
			cerr << "error: get likelihood called by slave\n";
			exit(1);
		}
		double ret = 0;
		if (priorsampling || fixomega)	{
			return 0;
		}
		MasterCollectLogLikelihood();
		double total = 0;
		for (int gene=0; gene<Ngene; gene++)	{
			total += genelnl[gene];
		}
		return total;
	}

	void MasterCollectLogLikelihood()	{
		MPI_Status stat;
		for (int j=1; j<nprocs; j++)	{
			MPI_Recv(tmpgenelnl,Ngene,MPI_DOUBLE,j,TAG1,MPI_COMM_WORLD,&stat);
			for (int gene=0; gene<Ngene; gene++)	{
				if (genealloc[gene] == j)	{
					genelnl[gene] = tmpgenelnl[gene];
				}
			}
		}
	}

	void SlaveSendLogLikelihood()	{
		for (int gene=0; gene<Ngene; gene++)	{
			if (genealloc[gene] == myid)	{
				if (clampsuffstat)	{
					genelnl[gene] = pathconjtree[gene]->GetLogProb();
				}
				else	{
					genelnl[gene] = phyloprocess[gene]->GetLogProb();
				}
			}
		}
		MPI_Send(genelnl,Ngene,MPI_DOUBLE,0,TAG1,MPI_COMM_WORLD);
	}

	void MakeScheduler()	{
		cerr << "error: in make scheduler\n";
		exit(1);
	}

	void MasterMakeScheduler();
	void SlaveMakeScheduler();

	void AllocateAlignments(string datafile);
	void CreateMissingMap(Tree* intree);
	bool RecursiveFillMissingMap(const Link* from, int gene);

	void MasterCreate();
	void SlaveCreate(bool sample);

	void MasterSendParameters();
	void SlaveReceiveParameters();

	void SlaveSendGeneProcesses();
	void MasterReceiveGeneProcesses();

	void MasterSendGeneProcesses();
	void SlaveReceiveGeneProcesses();

	void MasterRegisterMPINode();
	void SlaveRegisterMPINode();

	void CreateMPIParamArrays();
	void CreateMPIGeneProcessesArrays();

	/*
	double Move(double tuning = 1)	{
		if (! myid)	{
			scheduler.Cycle(1,1,true,true);
		}
		else	{
			scheduler.Cycle(1,1,false,false);
		}
		// Cycle(1,1,verbose,check)
		// scheduler.Cycle(1,1,true,true);
		return 1;
	}
	*/

	void drawSample()	{
		cerr << "in sample\n";
		exit(1);
	}

	Var<PosReal>* GetScale()	{
		if (isCalibrated())	{
			return GetCalibratedChronogram()->GetScale();
		}
		return 0;
	}

	double GetRootAge()	{
		if (isCalibrated())	{
			return GetCalibratedChronogram()->GetRootAge();
			// return GetCalibratedChronogram()->GetScale()->val();
		}
		return 1;
	}

	LengthTree* GetGeneOmegaTree(int gene)	{
		return geneomegatree[gene];
	}

	void UpdateGeneOmegaTree(int gene)	{
		geneomegatree[gene]->specialUpdate();
	}

	void UpdateGeneOmegaTrees()	{
		for (int gene=0; gene<Ngene; gene++)	{
			UpdateGeneOmegaTree(gene);
		}
	}

	double GetTotalTime()	{
		return chronogram->GetTotalTime();
	}

	double GetRegCoeff(int i, int j)	{
		return regarray->GetCell(i,j)->val();
	}

	double GetSynRegCoeff(int i, int j)	{
		return synregarray->GetCell(i,j)->val();
	}

	void SlaveTrace()	{
		SlaveSendLogLikelihood();
	}

	void TraceHeader(ostream& os)	{
		os << "#logprior\tlnL";
		int n = 8;
		if (Ngene < n)	{
			n = Ngene;
		}
		for (int gene=0; gene<n; gene++)	{
			os << '\t' << "syn" << gene;
			os << '\t' << "om" << gene;
		}
		os << "\tsyntaumean\tvar";
		os << "\ttaumean\tvar";
		for (int i=0; i<Ncont; i++)	{
			os << '\t' << "synregvar" << i;
		}
		if (! fixomega)	{
			os << "\tsynrootmean";
			os << "\tsynrootvar";
		}
		for (int i=0; i<Ncont; i++)	{
			os << '\t' << "regvar" << i;
		}
		if (! fixomega)	{
			os << "\trootmean";
			os << "\trootvar";
		}

		for (int k=0; k<Ncont+2; k++)	{
			for (int l=k+1; l<Ncont+2; l++)	{
				os << '\t' << "sigma" << k << '_' << l;
			}
		}
		for (int k=0; k<Ncont+2; k++)	{
			os << '\t' << "sigma" << k << '_' << k;
		}
		if (! fixomega)	{
			os << "\trelrate";
			os << "\trelrateconc";
			os << "\tstatent";
			os << "\tstatconc";
		}
		if (isCalibrated())	{
			os << "\trootage";
		}
		if ((chronoprior >= 1) && (chronoprior <= 3))	{
			os << "\tp1\tp2";
		}
		for (int k=0; k<Ncont; k++)	{
			os << '\t' << "root_" << k;
		}
		os << '\n';
	}

	// writes all summary statistics on one single line
	// in the same order as that provided by the header
	void Trace(ostream& os)	{

		os << GetLogPrior() << '\t' << GetLogLikelihood();
		int n = 8;
		if (Ngene < n)	{
			n = Ngene;
		}
		for (int gene=0; gene<n; gene++)	{
			os << '\t' << genesynprocess[gene]->GetMean();
			os << '\t' << geneprocess[gene]->GetMean();
		}

		os << '\t' << genesyntau->GetMean();
		os << '\t' << sqrt(genesyntau->GetVar()) / genesyntau->GetMean();
		os << '\t' << genetau->GetMean();
		os << '\t' << sqrt(genetau->GetVar()) / genetau->GetMean();

		for (int i=0; i<Ncont; i++)	{
			os << '\t' << synregarray->GetVar(i);
		}
		if (! fixomega)	{
			os << '\t' << *rootsynmean;
			os << '\t' << *rootsynvar;
		}

		for (int i=0; i<Ncont; i++)	{
			os << '\t' << regarray->GetVar(i);
		}
		if (! fixomega)	{
			os << '\t' << *rootomegamean;
			os << '\t' << *rootomegavar;
		}

		for (int k=0; k<Ncont+2; k++)	{
			for (int l=k+1; l<Ncont+2; l++)	{
				os << '\t' << (*sigma)[k][l];
			}
		}
		for (int k=0; k<Ncont+2; k++)	{
			os << '\t' << (*sigma)[k][k];
		}

		if (! fixomega)	{
			os << '\t' << relratecenter->GetEntropy();
			os << '\t' << *relrateconcentration;
			os << '\t' << stationarycenter->GetEntropy();
			os << '\t' << *stationaryconcentration;
		}

		if (isCalibrated())	{
			os << '\t' << GetRootAge();
		}
		if ((chronoprior >= 1) && (chronoprior <= 3))	{
			os << '\t' << *Chi << '\t' << *Chi2;
		}

		for (int k=0; k<Ncont; k++)	{
			os << '\t' << GetProcess()->GetMultiNormal(GetTree()->GetRoot())->val()[k];
		}

		os << '\n';
		os.flush();
	}

	void ToStream(ostream& os)	{
		os << *chronogram << '\n';
		if (isCalibrated())	{
			os << *GetCalibratedChronogram()->GetScale() << '\n';
		}
		if ((chronoprior >= 1) && (chronoprior <= 3))	{
			os << *Chi << '\t' << *Chi2 << '\n';
		}
		os << *DiagArray << '\n';
		os << *sigma << '\n';
		os << *process << '\n';

		os << *syntaushape << '\n';
		os << *syntauscale << '\n';
		os << *genesyntau << '\n';

		synregarray->ToStream(os);
		os << '\n';

		os << *rootsynmean << '\n';
		os << *rootsynvar << '\n';
		os << *generootsyn << '\n';
		for (int gene = 0; gene<Ngene; gene++)	{
			os << *genesynprocess[gene] << '\n';
		}

		os << *taushape << '\n';
		os << *tauscale << '\n';
		os << *genetau << '\n';

		regarray->ToStream(os);
		os << '\n';

		os << *rootomegamean << '\n';
		os << *rootomegavar << '\n';
		os << *generootomega << '\n';
		for (int gene = 0; gene<Ngene; gene++)	{
			os << *geneprocess[gene] << '\n';
		}

		if (! fixomega)	{
			os << *relratecenter << '\n';
			os << *stationarycenter << '\n';
			os << *relrateconcentration << '\n';
			os << *stationaryconcentration << '\n';
			os << *generelrate << '\n';
			os << *genestationary << '\n';
		}
	}

	void FromStream(istream& is)	{
		is >> *chronogram;
		if (isCalibrated())	{
			is >> *GetCalibratedChronogram()->GetScale();
		}
		if ((chronoprior >= 1) && (chronoprior <= 3))	{
			is >> *Chi >> *Chi2;
		}
		is >> *DiagArray;
		is >> *sigma;
		is >> *process;

		is >> *syntaushape;
		is >> *syntauscale;
		is >> *genesyntau;

		synregarray->FromStream(is);
		is >> *rootsynmean;
		is >> *rootsynvar;
		is >> *generootsyn;
		for (int gene = 0; gene<Ngene; gene++)	{
			is >> *genesynprocess[gene];
		}

		is >> *taushape;
		is >> *tauscale;
		is >> *genetau;

		regarray->FromStream(is);
		is >> *rootomegamean;
		is >> *rootomegavar;
		is >> *generootomega;
		for (int gene = 0; gene<Ngene; gene++)	{
			is >> *geneprocess[gene];
		}

		if (! fixomega)	{
			is >> *relratecenter;
			is >> *stationarycenter;
			is >> *relrateconcentration;
			is >> *stationaryconcentration;
			is >> *generelrate;
			is >> *genestationary;
		}
	}

	void WriteOmega(ostream& os)	{
		for (int gene = 0; gene<Ngene; gene++)	{
			os << *genesynprocess[gene] << '\n';
		}
		for (int gene = 0; gene<Ngene; gene++)	{
			os << *geneprocess[gene] << '\n';
		}
	}

	void ReadOmega(istream& is)	{
		for (int gene = 0; gene<Ngene; gene++)	{
			is >> *genesynprocess[gene];
		}
		for (int gene = 0; gene<Ngene; gene++)	{
			is >> *geneprocess[gene];
		}

		RecursiveSetToMean(GetTree()->GetRoot());
	}

	void RecursiveSetToMean(const Link* from)	{
		double meansyn = 0;
		for (int gene=0; gene<Ngene; gene++)	{
			meansyn += genesynprocess[gene]->GetNodeVal(from->GetNode())->val();
		}
		meansyn /= Ngene;
		(*process->GetNodeVal(from->GetNode()))[Ncont] = meansyn;

		double meanom = 0;
		for (int gene=0; gene<Ngene; gene++)	{
			meanom += genesynprocess[gene]->GetNodeVal(from->GetNode())->val();
		}
		meanom /= Ngene;
		(*process->GetNodeVal(from->GetNode()))[Ncont+1] = meanom;

		for (const Link* link=from->Next(); link!=from; link=link->Next())	{
			RecursiveSetToMean(link->Out());
		}
	}

	void ComputeKalmanSuffStat(ConditionalExternalKalmanMultiVariateTreeProcess* kalman);

	map<const Link*, double**>& GetKy() {return Ky;}
	double** GetW() {return W;}

	void RecursiveComputeKalmanSuffStat(const Link* from, ConditionalExternalKalmanMultiVariateTreeProcess* kalman);

	map<const Link*, double**> Ky;
	double** W;
};

void MixOmegaModel::AllocateAlignments(string datafile)	{

	ifstream is(datafile.c_str());
	is >> Ngene;
	genename = new string[Ngene];
	genesize = new int[Ngene];
	genealloc = new int[Ngene];
	genelnl = new double[Ngene];
	tmpgenelnl = new double[Ngene];
	nucdata = new FileSequenceAlignment*[Ngene];
	codondata = new CodonSequenceAlignment*[Ngene];
	int* geneweight = genesize;

	totnsite = 0;
	for (int gene=0; gene<Ngene; gene++)	{
		if (! myid)	{
			cerr << gene << '\t' << Ngene << '\t' << genename[gene] << '\n';
		}
		is >> genename[gene];
		nucdata[gene] = new FileSequenceAlignment(genename[gene]);
		codondata[gene] = new CodonSequenceAlignment(nucdata[gene],type);

		if (! gene)	{
			codonstatespace = codondata[gene]->GetCodonStateSpace();
		}
		genesize[gene] = codondata[gene]->GetNsite();
		totnsite += genesize[gene];
	}

	if ((! fixomega) && (nprocs > 1))	{
		// sort alignments by decreasing size
		int permut[Ngene];
		for (int gene=0; gene<Ngene; gene++)	{
			permut[gene] = gene;
		}
		for (int i=0; i<Ngene; i++)	{
			for (int j=Ngene-1; j>i; j--)	{
				if (geneweight[permut[i]] < geneweight[permut[j]])	{
				// if (genesize[permut[i]] < genesize[permut[j]])	{
					int tmp = permut[i];
					permut[i] = permut[j];
					permut[j] = tmp;
				}
			}
		}

		int totsize[nprocs];
		for (int i=0; i<nprocs; i++)	{
			totsize[i] = 0;
		}

		for (int i=0; i<Ngene; i++)	{
			int gene = permut[i];
			int size = geneweight[gene];
			// int size = genesize[gene];

			int min = 0;
			int jmin = 0;
			for (int j=1; j<nprocs; j++)	{
				if ((j==1) || (min > totsize[j]))	{
					min = totsize[j];
					jmin = j;
				}
			}
			genealloc[gene] = jmin;
			totsize[jmin] += size;
		}

		if (totsize[0])	{
			cerr << "error in alloc\n";
			exit(1);
		}

		procngene = new int[nprocs];
		for (int i=0; i<nprocs; i++)	{
			procngene[i] = 0;
		}

		int total = 0;
		for (int i=1; i<nprocs; i++)	{
			int tot = 0;
			for (int gene=0; gene<Ngene; gene++)	{
				if (genealloc[gene] == i)	{
					tot += geneweight[gene];
					// tot += genesize[gene];
					total++;
					procngene[genealloc[gene]]++;
				}
			}
			if (tot != totsize[i])	{
				cerr << "error in allocation\n";
				cerr << tot << '\t' << totsize[i] << '\n';
				exit(1);
			}
		}
		if (total != Ngene)	{
			cerr << "error in total allocation\n";
			exit(1);
		}
		MaxNgene = 0;
		for (int i=1; i<nprocs; i++)	{
			if (MaxNgene < procngene[i])	{
				MaxNgene = procngene[i];
			}
		}

		procnsite = new int[nprocs];
		for (int i=0; i<nprocs; i++)	{
			procnsite[i] = 0;
		}
		for (int gene=0; gene<Ngene; gene++)	{
			if ((genealloc[gene] < 0) || (genealloc[gene] >= nprocs))	{
				cerr << "alloc : " << genealloc[gene] << '\t' << gene << '\n';
				exit(1);
			}
			procnsite[0] += genesize[gene];
			procnsite[genealloc[gene]] += genesize[gene];
		}
	}
}

void MixOmegaModel::CreateMissingMap(Tree* intree)	{

	missingmap = new map<const Link*, int>[Ngene];
	for (int gene=0; gene<Ngene; gene++)	{
		RecursiveFillMissingMap(intree->GetRoot(),gene);
	}
}

bool MixOmegaModel::RecursiveFillMissingMap(const Link* from, int gene)	{

	bool some = false;
	if (from->isLeaf())	{
		bool miss = codondata[gene]->AllMissingTaxon(from->GetNode()->GetName());
		missingmap[gene][from] = miss;
		some = ! miss;
	}
	else	{
		for (const Link* link=from->Next(); link!=from; link=link->Next())	{
			some |= RecursiveFillMissingMap(link->Out(),gene);
		}
		missingmap[gene][from] = !some;
	}
	return some;
}

void MixOmegaModel::MasterCreate()	{

	// ----------
	// construction of the graph
	// ----------

	Zero = new Const<Real>(0);
	One = new Const<PosReal>(1);

	MeanChi = 0;
	MeanChi2 = 0;
	Chi = 0;
	Chi2 = 0;

	if (calibfile != "None")	{
		cerr << "calib\n";
		exit(1);
		/*
		iscalib = true;

		double a = rootage * rootage / rootstdev / rootstdev;
		double b = rootage / rootstdev / rootstdev;
		if (rootage == -1)	{
			a = b = -1;
		}
		CalibrationSet* calibset = new FileCalibrationSet(calibfile, tree);

		if (chronoprior == 0)	{
			chronogram = new CalibratedChronogram(tree,One,a,b,calibset);
		}
		else {
			cerr << "BD\n";
			MeanChi = new Const<PosReal>(meanchi);
			MeanChi2 = new Const<PosReal>(meanchi2);
			Chi = new Exponential(MeanChi,Exponential::MEAN);
			Chi2 = new Exponential(MeanChi2,Exponential::MEAN);
			chronogram = new BDCalibratedChronogram(tree,One,Chi,Chi2,a,b,calibset,chronoprior);
		}
		*/
	}
	else	{
		chronogram = new Chronogram(tree,One);
	}

	if (clamptree)	{
		chronogram->Clamp();
	}

	DiagArray = new JeffreysIIDArray(Ncont+2,minjeff,maxjeff,Zero);
	DiagArray->setval(0.1);

	sigmaZero = new SigmaZero(DiagArray);

	if (clampdiag)	{
		sigma = new DiagonalCovMatrix(sigmaZero,Ncont+2+df);
	}
	else	{
		// sigma = new InverseWishartMatrix(sigmaZero,Ncont+2+df);
		sigma = new ConjugateInverseWishart(sigmaZero,Ncont+2+df);
	}

	rootmean = 0;
	rootvar = 0;
	/*
	if (rootfile != "None")	{
		cerr << "check if root mean and var compatible with conjugate multivariate process\n";
		exit(1);
		rootmean = new Const<RealVector>(RealVector(Ncont+2));
		rootvar = new Const<PosRealVector>(PosRealVector(Ncont+2));
		(*rootmean)[Ncont] = 0;
		(*rootvar)[Ncont] = 0;
		(*rootmean)[Ncont+1] = 0;
		(*rootvar)[Ncont+1] = 0;
		ifstream is(rootfile.c_str());
		for (int i=0; i<Ncont; i++)	{
			double mean, var;
			is >> mean >> var;
			(*rootmean)[i] = mean;
			(*rootvar)[i] = var;
			cerr << mean << '\t' << var << '\n';
		}
	}
	*/

	if (clampdiag)	{
		process = new MultiVariateTreeProcess(sigma,chronogram,0,0,rootmean,rootvar);
	}
	else	{
		// process = new MultiVariateTreeProcess(sigma,chronogram,0,0,rootmean, rootvar);
		// process = new ConjugateExternalKalmanMultiVariateTreeProcess(GetConjugateInverseWishart(),chronogram,0,0,rootmean, rootvar);
		process = new ConjugateConditionalExternalKalmanMultiVariateTreeProcess(GetConjugateInverseWishart(),chronogram,0,0,rootmean, rootvar);
	}

	process->Reset();
	process->GetMultiNormal(process->GetRoot())->ClampAt(0,Ncont);
	process->GetMultiNormal(process->GetRoot())->ClampAt(0,Ncont+1);

	// here make the regcoef array

	synregarray = new BidimIIDUniform(Ngene,Ncont+2,minuni,maxuni,Zero);
	synregarray->SetAt(1.0,Ncont);
	synregarray->Clamp(Ncont);
	synregarray->SetAt(0,Ncont+1);
	synregarray->Clamp(Ncont+1);

	if (clampreg)	{
		cerr << "clamp reg\n";
		for (int i=0; i<Ncont; i++)	{
			synregarray->SetAt(0,i);
			synregarray->Clamp(i);
		}
	}
	else	{
		for (int i=0; i<Ncont; i++)	{
			synregarray->SetAtRandom(0,0.1,i);
		}
	}

	cerr << "tau\n";
	syntaushape = new Jeffreys(minjeff,maxjeff,Zero);
	syntauscale = new Jeffreys(minjeff,maxjeff,Zero);
	syntaushape->setval(10.0);
	syntauscale->setval(100.0);
	genesyntau = new GammaIIDArray(Ngene,syntaushape,syntauscale);
	for (int gene=0; gene<Ngene; gene++)	{
		genesyntau->GetVal(gene)->setval(0.1 + 0.01 * Random::Uniform());
		// genesyntau->GetVal(gene)->ClampAt(1.0);
	}

	cerr << "gene process\n";

	genesynprocess = new LinRegNormalProcess*[Ngene];
	for (int gene=0; gene<Ngene; gene++)	{
		genesynprocess[gene] = new LinRegNormalProcess(chronogram,process,synregarray,0,gene,genesyntau->GetVal(gene));
		genesynprocess[gene]->Reset();
		genesynprocess[gene]->ClampRoot(0);
	}

	cerr << "root vals\n";
	rootsynmean = new Uniform(minuni,maxuni,Zero);
	rootsynmean->setval(-1.0);
	rootsynvar = new Jeffreys(minjeff,maxjeff,Zero);
	rootsynvar->setval(0.1);
	generootsyn = new NormalIIDArray(Ngene,rootsynmean,rootsynvar);
	// generootsyn->Reset();
	// does not work anyway (clamps the Array as a MCMC, not its entries)
	// generootsyn->Clamp();

	genesyntree = 0;
	/*
	if (! fixomega)	{
		cerr << "syn trees\n";
		genesyntree = new MeanExpTree*[Ngene];
		for (int gene=0; gene<Ngene; gene++)	{
			genesyntree[gene] = new MeanExpTree(genesynprocess[gene],chronogram,INTEGRAL,generootsyn->GetVal(gene));
		}
	}
	*/

	regarray = new BidimIIDUniform(Ngene,Ncont+2,minuni,maxuni,Zero);
	regarray->SetAt(0,Ncont);
	regarray->Clamp(Ncont);
	regarray->SetAt(1.0,Ncont+1);
	regarray->Clamp(Ncont+1);

	if (clampreg)	{
		cerr << "clamp reg\n";
		for (int i=0; i<Ncont; i++)	{
			regarray->SetAt(0,i);
			regarray->Clamp(i);
		}
	}
	else	{
		for (int i=0; i<Ncont; i++)	{
			regarray->SetAtRandom(0,0.1,i);
		}
	}


	cerr << "tau\n";
	taushape = new Jeffreys(minjeff,maxjeff,Zero);
	tauscale = new Jeffreys(minjeff,maxjeff,Zero);
	taushape->setval(10.0);
	tauscale->setval(100.0);
	genetau = new GammaIIDArray(Ngene,taushape,tauscale);
	for (int gene=0; gene<Ngene; gene++)	{
		// genetau->GetVal(gene)->setval(0.1);
		genetau->GetVal(gene)->setval(0.1 + 0.01 * Random::Uniform());
		// genetau->GetVal(gene)->ClampAt(1.0);
	}

	cerr << "gene process\n";

	geneprocess = new LinRegNormalProcess*[Ngene];
	for (int gene=0; gene<Ngene; gene++)	{
		geneprocess[gene] = new LinRegNormalProcess(chronogram,process,regarray,0,gene,genetau->GetVal(gene));
		geneprocess[gene]->Reset();
		geneprocess[gene]->ClampRoot(0);
	}

	cerr << "root vals\n";
	rootomegamean = new Uniform(minuni,maxuni,Zero);
	rootomegamean->setval(-1.0);
	rootomegavar = new Jeffreys(minjeff,maxjeff,Zero);
	rootomegavar->setval(0.1);
	generootomega = new NormalIIDArray(Ngene,rootomegamean,rootomegavar);
	// generootomega->Reset();
	// does not work anyway (clamps the Array as a MCMC, not its entries)
	// generootomega->Clamp();

	geneomegatree = 0;
	/*
	if (! fixomega)	{
		cerr << "omega trees\n";
		geneomegatree = new MeanExpTree*[Ngene];
		for (int gene=0; gene<Ngene; gene++)	{
			geneomegatree[gene] = new MeanExpTree(geneprocess[gene],chronogram,MEAN,generootomega->GetVal(gene));
		}
	}
	*/

	if (fixomega)	{
		ifstream is(omegafile.c_str());
		ReadOmega(is);
	}

	if (! fixomega)	{
		cerr << "matrix\n";
		// create a GTR nucleotide matrix
		relratecenter = new Dirichlet(Nnuc*(Nnuc-1)/2);
		stationarycenter = new Dirichlet(Nnuc);
		relratecenter->setuniform();
		stationarycenter->setuniform();
		relrateconcentration = new Jeffreys(minjeff,maxjeff,Zero);
		relrateconcentration->setval(100);
		stationaryconcentration = new Jeffreys(minjeff,maxjeff,Zero);
		stationaryconcentration->setval(100);

		generelrate = new DirichletIIDArray(Ngene,relratecenter,relrateconcentration);
		genestationary = new DirichletIIDArray(Ngene,stationarycenter,stationaryconcentration);

		/*
		nucmatrix = new GTRRandomSubMatrixWithNormRates*[Ngene];
		for (int gene=0; gene < Ngene; gene++)	{
			nucmatrix[gene] = new GTRRandomSubMatrixWithNormRates(generelrate->GetVal(gene),genestationary->GetVal(gene),true);
		}

		cerr << "codon matrix trees\n";
		genecodonmatrixtree = new MatrixTree*[Ngene];
		for (int gene=0; gene<Ngene; gene++)	{
			genecodonmatrixtree[gene] = new MatrixTree(codondata[gene]->GetCodonStateSpace(), nucmatrix[gene], geneomegatree[gene], One);
		}
		*/
	}

	cerr << "set and clamp\n";
	if (contdata)	{
		for (int i=0; i<Ncont; i++)	{
			process->SetAndClamp(contdata,i,i);
		}

		if (contml)	{
			for (int k=0; k<Ncont; k++)	{
				process->ML(k);
			}
		}

	}
	/*
	for (int i=0; i<Ncont+2; i++)	{
		process->ClampLeaves(i);
	}
	*/

	/*
	if (! fixomega)	{
		cerr << "phyloprocess\n";
		// make substitution mappings
		if (conjpath)	{
			pathconjtree = new PathConjugateTree*[Ngene];
			phyloprocess = new PhyloProcess*[Ngene];
			for (int gene=0; gene<Ngene; gene++)	{
				pathconjtree[gene] = new MGCodonBranchMatrixPathConjugateTree(genesyntree[gene], genecodonmatrixtree[gene], GetData(gene));
				phyloprocess[gene] = new PathConjugatePhyloProcess(pathconjtree[gene]);
			}
		}
		else	{
			pathconjtree = 0;
			if (priorsampling)	{
				phyloprocess = 0;
			}
			phyloprocess = new PhyloProcess*[Ngene];
			for (int gene=0; gene<Ngene; gene++)	{
				phyloprocess[gene] = new BranchMatrixPhyloProcess(genesyntree[gene], genecodonmatrixtree[gene], GetData(gene));
			}
		}
	}
	else	{
		phyloprocess = 0;
	}

	if (phyloprocess)	{
		cerr << "unfold\n";
		for (int gene=0; gene<Ngene; gene++)	{
			phyloprocess[gene]->Unfold();
		}
	}
	if (sample)	{
		cerr << "sample\n";
		if (phyloprocess)	{
			for (int gene=0; gene<Ngene; gene++)	{
				phyloprocess[gene]->Sample();
			}
		}
	}
	*/

	cerr << "register\n";
	// register model
	RootRegister(Zero);
	RootRegister(One);
	if (chronoprior)	{
		RootRegister(MeanChi);
		RootRegister(MeanChi2);
	}
	if (rootmean)	{
		RootRegister(rootmean);
		RootRegister(rootvar);
	}
	if (! fixomega)	{
		RootRegister(relratecenter);
		RootRegister(stationarycenter);
	}
	Register();

	cerr << "model created\n";
}

void MixOmegaModel::SlaveCreate(bool sample)	{

	// ----------
	// construction of the graph
	// ----------

	Zero = new Const<Real>(0);
	One = new Const<PosReal>(1);

	MeanChi = 0;
	MeanChi2 = 0;
	Chi = 0;
	Chi2 = 0;

	if (calibfile != "None")	{
		iscalib = true;
		/*
		double a = rootage * rootage / rootstdev / rootstdev;
		double b = rootage / rootstdev / rootstdev;
		if (rootage == -1)	{
			a = b = -1;
		}
		CalibrationSet* calibset = new FileCalibrationSet(calibfile, tree);

		if (chronoprior == 0)	{
			chronogram = new CalibratedChronogram(tree,One,a,b,calibset);
		}
		else {
			cerr << "BD\n";
			MeanChi = new Const<PosReal>(meanchi);
			MeanChi2 = new Const<PosReal>(meanchi2);
			Chi = new Exponential(MeanChi,Exponential::MEAN);
			Chi2 = new Exponential(MeanChi2,Exponential::MEAN);
			chronogram = new BDCalibratedChronogram(tree,One,Chi,Chi2,a,b,calibset,chronoprior);
		}
		*/
	}
	else	{
		chronogram = new Chronogram(tree,One);
	}

	if (clamptree)	{
		chronogram->Clamp();
	}

	DiagArray = new JeffreysIIDArray(Ncont+2,minjeff,maxjeff,Zero);
	DiagArray->setval(0.1);

	sigmaZero = new SigmaZero(DiagArray);

	if (clampdiag)	{
		sigma = new DiagonalCovMatrix(sigmaZero,Ncont+2+df);
	}
	else	{
		// sigma = new InverseWishartMatrix(sigmaZero,Ncont+2+df);
		sigma = new ConjugateInverseWishart(sigmaZero,Ncont+2+df);
	}

	rootmean = 0;
	rootvar = 0;
	/*
	if (rootfile != "None")	{
		cerr << "check if root mean and var compatible with conjugate multivariate process\n";
		exit(1);
		rootmean = new Const<RealVector>(RealVector(Ncont+2));
		rootvar = new Const<PosRealVector>(PosRealVector(Ncont+2));
		(*rootmean)[Ncont] = 0;
		(*rootvar)[Ncont] = 0;
		(*rootmean)[Ncont+1] = 0;
		(*rootvar)[Ncont+1] = 0;
		ifstream is(rootfile.c_str());
		for (int i=0; i<Ncont; i++)	{
			double mean, var;
			is >> mean >> var;
			(*rootmean)[i] = mean;
			(*rootvar)[i] = var;
			cerr << mean << '\t' << var << '\n';
		}
	}
	*/

	if (clampdiag)	{
		process = new MultiVariateTreeProcess(sigma,chronogram,0,0,rootmean,rootvar);
	}
	else	{
		// process = new MultiVariateTreeProcess(sigma,chronogram,0,0,rootmean, rootvar);
		process = new ConjugateMultiVariateTreeProcess(GetConjugateInverseWishart(),chronogram,0,0,rootmean, rootvar);
	}

	process->Reset();
	process->GetMultiNormal(process->GetRoot())->ClampAt(0,Ncont);
	process->GetMultiNormal(process->GetRoot())->ClampAt(0,Ncont+1);

	// here make the regcoef array

	synregarray = new BidimIIDUniform(Ngene,Ncont+2,minuni,maxuni,Zero);
	synregarray->SetAt(1.0,Ncont);
	synregarray->Clamp(Ncont);
	synregarray->SetAt(0,Ncont+1);
	synregarray->Clamp(Ncont+1);

	if (clampreg)	{
		cerr << "clamp reg\n";
		for (int i=0; i<Ncont; i++)	{
			synregarray->SetAt(0,i);
			synregarray->Clamp(i);
		}
	}
	else	{
		for (int i=0; i<Ncont; i++)	{
			synregarray->SetAtRandom(0,0.1,i);
		}
	}


	cerr << "tau\n";
	syntaushape = new Jeffreys(minjeff,maxjeff,Zero);
	syntauscale = new Jeffreys(minjeff,maxjeff,Zero);
	syntaushape->setval(10.0);
	syntauscale->setval(100.0);
	genesyntau = new GammaIIDArray(Ngene,syntaushape,syntauscale);
	for (int gene=0; gene<Ngene; gene++)	{
		// genesyntau->GetVal(gene)->setval(0.1);
		genesyntau->GetVal(gene)->setval(0.1 + 0.01 * Random::Uniform());
		// genesyntau->GetVal(gene)->ClampAt(1.0);
	}

	cerr << "gene process\n";

	genesynprocess = new LinRegNormalProcess*[Ngene];
	for (int gene=0; gene<Ngene; gene++)	{
		if (genealloc[gene] == myid)	{
			genesynprocess[gene] = new LinRegNormalProcess(chronogram,process,synregarray,0,gene,genesyntau->GetVal(gene));
			genesynprocess[gene]->Reset();
			genesynprocess[gene]->ClampRoot(0);
		}
		else	{
			genesynprocess[gene] = 0;
		}
	}

	cerr << "root vals\n";
	rootsynmean = new Uniform(minuni,maxuni,Zero);
	rootsynmean->setval(-1.0);
	rootsynvar = new Jeffreys(minjeff,maxjeff,Zero);
	rootsynvar->setval(0.1);
	generootsyn = new NormalIIDArray(Ngene,rootsynmean,rootsynvar);
	// generootsyn->Reset();
	// does not work anyway (clamps the Array as a MCMC, not its entries)
	// generootsyn->Clamp();

	if (! fixomega)	{
		cerr << "syn trees\n";
		genesyntree = new MeanExpTree*[Ngene];
		for (int gene=0; gene<Ngene; gene++)	{
			if (genealloc[gene] == myid)	{
				genesyntree[gene] = new MeanExpTree(genesynprocess[gene],chronogram,INTEGRAL,generootsyn->GetVal(gene));
			}
		}
	}

	regarray = new BidimIIDUniform(Ngene,Ncont+2,minuni,maxuni,Zero);
	regarray->SetAt(0,Ncont);
	regarray->Clamp(Ncont);
	regarray->SetAt(1.0,Ncont+1);
	regarray->Clamp(Ncont+1);

	if (clampreg)	{
		cerr << "clamp reg\n";
		for (int i=0; i<Ncont; i++)	{
			regarray->SetAt(0,i);
			regarray->Clamp(i);
		}
	}
	else	{
		for (int i=0; i<Ncont; i++)	{
			regarray->SetAtRandom(0,0.1,i);
		}
	}


	cerr << "tau\n";
	taushape = new Jeffreys(minjeff,maxjeff,Zero);
	tauscale = new Jeffreys(minjeff,maxjeff,Zero);
	taushape->setval(10.0);
	tauscale->setval(100.0);
	genetau = new GammaIIDArray(Ngene,taushape,tauscale);
	for (int gene=0; gene<Ngene; gene++)	{
		// genetau->GetVal(gene)->setval(0.1);
		genetau->GetVal(gene)->setval(0.1 + 0.01 * Random::Uniform());
		// genetau->GetVal(gene)->ClampAt(1.0);
	}

	cerr << "gene process\n";

	geneprocess = new LinRegNormalProcess*[Ngene];
	for (int gene=0; gene<Ngene; gene++)	{
		if (genealloc[gene] == myid)	{
			geneprocess[gene] = new LinRegNormalProcess(chronogram,process,regarray,0,gene,genetau->GetVal(gene));
			geneprocess[gene]->Reset();
			geneprocess[gene]->ClampRoot(0);
		}
	}

	cerr << "root vals\n";
	rootomegamean = new Uniform(minuni,maxuni,Zero);
	rootomegamean->setval(-1.0);
	rootomegavar = new Jeffreys(minjeff,maxjeff,Zero);
	rootomegavar->setval(0.1);
	generootomega = new NormalIIDArray(Ngene,rootomegamean,rootomegavar);
	// generootomega->Reset();
	// does not work anyway (clamps the Array as a MCMC, not its entries)
	// generootomega->Clamp();

	if (! fixomega)	{
		cerr << "omega trees\n";
		geneomegatree = new MeanExpTree*[Ngene];
		for (int gene=0; gene<Ngene; gene++)	{
			if (genealloc[gene] == myid)	{
				geneomegatree[gene] = new MeanExpTree(geneprocess[gene],chronogram,MEAN,generootomega->GetVal(gene));
			}
		}
	}

	if (fixomega)	{
		cerr << "in slave: fix omega\n";
		exit(1);
		/*
		ifstream is(omegafile.c_str());
		ReadOmega(is);
		*/
	}

	if (! fixomega)	{
		cerr << "matrix\n";
		// create a GTR nucleotide matrix
		relratecenter = new Dirichlet(Nnuc*(Nnuc-1)/2);
		stationarycenter = new Dirichlet(Nnuc);
		relratecenter->setuniform();
		stationarycenter->setuniform();

		relrateconcentration = new Jeffreys(minjeff,maxjeff,Zero);
		relrateconcentration->setval(100);
		stationaryconcentration = new Jeffreys(minjeff,maxjeff,Zero);
		stationaryconcentration->setval(100);

		generelrate = new DirichletIIDArray(Ngene,relratecenter,relrateconcentration);
		genestationary = new DirichletIIDArray(Ngene,stationarycenter,stationaryconcentration);

		nucmatrix = new GTRRandomSubMatrixWithNormRates*[Ngene];

		for (int gene=0; gene < Ngene; gene++)	{
			if (genealloc[gene] == myid)	{
				nucmatrix[gene] = new GTRRandomSubMatrixWithNormRates(generelrate->GetVal(gene),genestationary->GetVal(gene),true);
			}
		}

		cerr << "codon matrix trees\n";
		genecodonmatrixtree = new MatrixTree*[Ngene];
		for (int gene=0; gene<Ngene; gene++)	{
			if (genealloc[gene] == myid)	{
				genecodonmatrixtree[gene] = new MatrixTree(codondata[gene]->GetCodonStateSpace(), nucmatrix[gene], geneomegatree[gene], One, &missingmap[gene]);
			}
		}
	}

	cerr << "set and clamp\n";
	if (contdata)	{
		for (int i=0; i<Ncont; i++)	{
			process->SetAndClamp(contdata,i,i);
		}

		if (contml)	{
			for (int k=0; k<Ncont; k++)	{
				process->ML(k);
			}
		}

	}
	/*
	for (int i=0; i<Ncont+2; i++)	{
		process->ClampLeaves(i);
	}
	*/

	if (! fixomega)	{
		cerr << "phyloprocess\n";
		// make substitution mappings
		if (conjpath)	{
			pathconjtree = new PathConjugateTree*[Ngene];
			phyloprocess = new PhyloProcess*[Ngene];
			for (int gene=0; gene<Ngene; gene++)	{
				if (genealloc[gene] == myid)	{
					pathconjtree[gene] = new MGCodonBranchMatrixPathConjugateTree(genesyntree[gene], genecodonmatrixtree[gene], GetData(gene));
					phyloprocess[gene] = new PathConjugatePhyloProcess(pathconjtree[gene]);
				}
			}
		}
		else	{
			pathconjtree = 0;
			if (priorsampling)	{
				phyloprocess = 0;
			}
			phyloprocess = new PhyloProcess*[Ngene];
			for (int gene=0; gene<Ngene; gene++)	{
				if (genealloc[gene] == myid)	{
					phyloprocess[gene] = new BranchMatrixPhyloProcess(genesyntree[gene], genecodonmatrixtree[gene], GetData(gene));
				}
			}
		}
	}
	else	{
		phyloprocess = 0;
	}

	if (phyloprocess)	{
		cerr << "unfold\n";
		for (int gene=0; gene<Ngene; gene++)	{
			if (genealloc[gene] == myid)	{
				phyloprocess[gene]->Unfold();
			}
		}
	}
	if (sample)	{
		cerr << "sample\n";
		if (phyloprocess)	{
			for (int gene=0; gene<Ngene; gene++)	{
				if (genealloc[gene] == myid)	{
					phyloprocess[gene]->Sample();
				}
			}
		}
	}

	cerr << "register\n";
	// register model
	RootRegister(Zero);
	RootRegister(One);
	if (chronoprior)	{
		RootRegister(MeanChi);
		RootRegister(MeanChi2);
	}
	if (rootmean)	{
		RootRegister(rootmean);
		RootRegister(rootvar);
	}
	if (! fixomega)	{
		RootRegister(relratecenter);
		RootRegister(stationarycenter);
	}
	Register();

	cerr << "model created\n";
}

void MixOmegaModel::CreateMPIParamArrays()	{

	mpiparamdim = chronogram->GetNinternalNode() + (Ncont + 2) * process->GetNnode() + 4 + Nnuc*(Nnuc-1)/2 + 1 + Nnuc + 1 + 4;
	mpiparamarray = new double[mpiparamdim];
}

void MixOmegaModel::CreateMPIGeneProcessesArrays()	{

	mpiprocessdim = 2 * (process->GetNnode() + 1) + Nnuc*(Nnuc-1)/2 + Nnuc + Ncont * 2 + 2;
	mpiprocessarray = new double[MaxNgene * mpiprocessdim];
}

void MixOmegaModel::MasterSendGeneProcesses()	{

	for (int j=1; j<nprocs; j++)	{

		double* ptr = mpiprocessarray;

		for (int gene=0; gene<Ngene; gene++)	{
			if (genealloc[gene] == j)	{
				ptr = genesynprocess[gene]->GetNodeVals(ptr);
				(*ptr++) = generootsyn->GetVal(gene)->val();
				ptr = geneprocess[gene]->GetNodeVals(ptr);
				(*ptr++) = generootomega->GetVal(gene)->val();
				for (int i=0; i<Nnuc*(Nnuc-1)/2; i++)	{
					(*ptr++) = (*generelrate->GetVal(gene))[i];
				}
				for (int i=0; i<Nnuc; i++)	{
					(*ptr++) = (*genestationary->GetVal(gene))[i];
				}
				for (int i=0; i<Ncont; i++)	{
					(*ptr++) = synregarray->GetCell(gene,i)->val();
					(*ptr++) = regarray->GetCell(gene,i)->val();
				}
				(*ptr++) = genesyntau->GetVal(gene)->val();
				(*ptr++) = genetau->GetVal(gene)->val();
			}
		}

		if (ptr - mpiprocessarray != procngene[j] * mpiprocessdim)	{
			cerr << "error in send processes: non matching dim\n";
			exit(1);
		}

		MPI_Send(mpiprocessarray,procngene[j]*mpiprocessdim,MPI_DOUBLE,j,TAG1,MPI_COMM_WORLD);
	}
}

void MixOmegaModel::SlaveReceiveGeneProcesses()	{

	// return;
	mpinode->Corrupt(false);

	MPI_Status stat;
	MPI_Recv(mpiprocessarray,procngene[myid]*mpiprocessdim,MPI_DOUBLE,0,TAG1,MPI_COMM_WORLD,&stat);

	double* ptr = mpiprocessarray;

	for (int gene=0; gene<Ngene; gene++)	{
		if (genealloc[gene] == myid)	{
			ptr = genesynprocess[gene]->SetNodeVals(ptr);
			generootsyn->GetVal(gene)->setval(*ptr++);
			ptr = geneprocess[gene]->SetNodeVals(ptr);
			generootomega->GetVal(gene)->setval(*ptr++);
			for (int i=0; i<Nnuc*(Nnuc-1)/2; i++)	{
				(*generelrate->GetVal(gene))[i] = (*ptr++);
			}
			for (int i=0; i<Nnuc; i++)	{
				(*genestationary->GetVal(gene))[i] = (*ptr++);
			}
			for (int i=0; i<Ncont; i++)	{
				synregarray->GetCell(gene,i)->setval(*ptr++);
				regarray->GetCell(gene,i)->setval(*ptr++);
			}
			genesyntau->GetVal(gene)->setval(*ptr++);
			genetau->GetVal(gene)->setval(*ptr++);
		}
	}

	if (ptr - mpiprocessarray != procngene[myid] * mpiprocessdim)	{
		cerr << "error in send processes: non matching dim\n";
		exit(1);
	}

	mpinode->Update();
}

void MixOmegaModel::SlaveSendGeneProcesses()	{

	// return;
	double* ptr = mpiprocessarray;
	for (int gene=0; gene<Ngene; gene++)	{
		if (genealloc[gene] == myid)	{
			ptr = genesynprocess[gene]->GetNodeVals(ptr);
			(*ptr++) = generootsyn->GetVal(gene)->val();
			ptr = geneprocess[gene]->GetNodeVals(ptr);
			(*ptr++) = generootomega->GetVal(gene)->val();
			for (int i=0; i<Nnuc*(Nnuc-1)/2; i++)	{
				(*ptr++) = (*generelrate->GetVal(gene))[i];
			}
			for (int i=0; i<Nnuc; i++)	{
				(*ptr++) = (*genestationary->GetVal(gene))[i];
			}
			for (int i=0; i<Ncont; i++)	{
				(*ptr++) = synregarray->GetCell(gene,i)->val();
				(*ptr++) = regarray->GetCell(gene,i)->val();
			}
			(*ptr++) = genesyntau->GetVal(gene)->val();
			(*ptr++) = genetau->GetVal(gene)->val();
		}
	}

	if (ptr - mpiprocessarray != procngene[myid] * mpiprocessdim)	{
		cerr << "error in send processes: non matching dim\n";
		exit(1);
	}

	MPI_Send(mpiprocessarray,procngene[myid]*mpiprocessdim,MPI_DOUBLE,0,TAG1,MPI_COMM_WORLD);
}

void MixOmegaModel::MasterReceiveGeneProcesses()	{

	// return;
	mpinode->Corrupt(false);

	MPI_Status stat;
	for (int j=1; j<nprocs; j++)	{

		MPI_Recv(mpiprocessarray,procngene[j]*mpiprocessdim,MPI_DOUBLE,j,TAG1,MPI_COMM_WORLD,&stat);

		double* ptr = mpiprocessarray;

		for (int gene=0; gene<Ngene; gene++)	{
			if (genealloc[gene] == j)	{
				ptr = genesynprocess[gene]->SetNodeVals(ptr);
				generootsyn->GetVal(gene)->setval(*ptr++);
				ptr = geneprocess[gene]->SetNodeVals(ptr);
				generootomega->GetVal(gene)->setval(*ptr++);
				for (int i=0; i<Nnuc*(Nnuc-1)/2; i++)	{
					(*generelrate->GetVal(gene))[i] = (*ptr++);
				}
				for (int i=0; i<Nnuc; i++)	{
					(*genestationary->GetVal(gene))[i] = (*ptr++);
				}
				for (int i=0; i<Ncont; i++)	{
					synregarray->GetCell(gene,i)->setval(*ptr++);
					regarray->GetCell(gene,i)->setval(*ptr++);
				}
				genesyntau->GetVal(gene)->setval(*ptr++);
				genetau->GetVal(gene)->setval(*ptr++);
			}
		}

		if (ptr - mpiprocessarray != procngene[j] * mpiprocessdim)	{
			cerr << "error in send processes: non matching dim\n";
			exit(1);
		}
	}

	mpinode->Update();
}

void MixOmegaModel::MasterRegisterMPINode()	{

	mpinode = new Mnode;
	for (int gene=0; gene<Ngene; gene++)	{
		genesynprocess[gene]->RegisterNodeTree(mpinode);
		generootsyn->GetVal(gene)->Register(mpinode);
		geneprocess[gene]->RegisterNodeTree(mpinode);
		generootomega->GetVal(gene)->Register(mpinode);
		generelrate->GetVal(gene)->Register(mpinode);
		genestationary->GetVal(gene)->Register(mpinode);
	}
	synregarray->RegisterArray(mpinode);
	genesyntau->RegisterArray(mpinode);
	regarray->RegisterArray(mpinode);
	genetau->RegisterArray(mpinode);
}

void MixOmegaModel::SlaveRegisterMPINode()	{

	mpinode = new Mnode;
	// chronogram??
	process->RegisterNodeTree(mpinode);
	rootsynmean->Register(mpinode);
	rootsynvar->Register(mpinode);
	rootomegamean->Register(mpinode);
	rootomegavar->Register(mpinode);
	relratecenter->Register(mpinode);
	relrateconcentration->Register(mpinode);
	stationarycenter->Register(mpinode);
	stationaryconcentration->Register(mpinode);
	syntaushape->Register(mpinode);
	syntauscale->Register(mpinode);
	taushape->Register(mpinode);
	tauscale->Register(mpinode);
}

void MixOmegaModel::MasterSendParameters()	{

	// return;
	double* ptr = mpiparamarray;
	// send
	// chronogram ?
	// process
	// synregarray
	// genesyntau
	// regarray
	// genetau
	ptr = chronogram->GetDates(ptr);
	ptr = process->GetNodeVals(ptr);

	// rootsynmean and var
	// root mean and var
	(*ptr++) = rootsynmean->val();
	(*ptr++) = rootsynvar->val();
	(*ptr++) = rootomegamean->val();
	(*ptr++) = rootomegavar->val();

	(*ptr++) = syntaushape->val();
	(*ptr++) = syntauscale->val();
	(*ptr++) = taushape->val();
	(*ptr++) = tauscale->val();


	// relrate center and concentration
	for (int i=0; i<Nnuc*(Nnuc-1)/2; i++)	{
		(*ptr++) = (*relratecenter)[i];
	}
	(*ptr++) = relrateconcentration->val();

	// stat center and concentration
	for (int i=0; i<Nnuc; i++)	{
		(*ptr++) = (*stationarycenter)[i];
	}
	(*ptr++) = stationaryconcentration->val();

	if (ptr - mpiparamarray != mpiparamdim)	{
		cerr << "error in send param: non matching dimension\n";
		exit(1);
	}

	MPI_Bcast(mpiparamarray,mpiparamdim,MPI_DOUBLE,0,MPI_COMM_WORLD);
}

void MixOmegaModel::SlaveReceiveParameters()	{

	// return;
	MPI_Bcast(mpiparamarray,mpiparamdim,MPI_DOUBLE,0,MPI_COMM_WORLD);

	mpinode->Corrupt(false);

	double* ptr = mpiparamarray;
	// send
	// chronogram ?
	// process
	// synregarray
	// genesyntau
	// regarray
	// genetau
	ptr = chronogram->SetDates(ptr);
	ptr = process->SetNodeVals(ptr);

	// rootsynmean and var
	// root mean and var
	rootsynmean->setval(*ptr++);
	rootsynvar->setval(*ptr++);
	rootomegamean->setval(*ptr++);
	rootomegavar->setval(*ptr++);

	syntaushape->setval(*ptr++);
	syntauscale->setval(*ptr++);
	taushape->setval(*ptr++);
	tauscale->setval(*ptr++);

	// relrate center and concentration
	for (int i=0; i<Nnuc*(Nnuc-1)/2; i++)	{
		(*relratecenter)[i] = (*ptr++);
	}
	relrateconcentration->setval(*ptr++);

	// stat center and concentration
	for (int i=0; i<Nnuc; i++)	{
		(*stationarycenter)[i] = (*ptr++);
	}
	stationaryconcentration->setval(*ptr++);

	if (ptr - mpiparamarray != mpiparamdim)	{
		cerr << "error in receive param: non matching dimension\n";
		exit(1);
	}
	mpinode->Update();
}


class MasterSendParametersMove : public MCUpdate	{

	public:

	MasterSendParametersMove(MixOmegaModel* inmodel)	{
		model = inmodel;
	}

	double Move(double tuning = 1)	{
		model->MasterSendParameters();
		return 1;
	}

	private:

	MixOmegaModel* model;
};

class SlaveReceiveParametersMove : public MCUpdate	{

	public:

	SlaveReceiveParametersMove(MixOmegaModel* inmodel)	{
		model = inmodel;
	}

	double Move(double tuning = 1)	{
		model->SlaveReceiveParameters();
		return 1;
	}

	private:

	MixOmegaModel* model;
};

class SlaveSendGeneProcessesMove : public MCUpdate	{

	public:

	SlaveSendGeneProcessesMove(MixOmegaModel* inmodel)	{
		model = inmodel;
	}

	double Move(double tuning = 1)	{
		model->SlaveSendGeneProcesses();
		return 1;
	}

	private:

	MixOmegaModel* model;
};

class MasterReceiveGeneProcessesMove : public MCUpdate	{

	public:

	MasterReceiveGeneProcessesMove(MixOmegaModel* inmodel)	{
		model = inmodel;
	}

	double Move(double tuning = 1)	{
		model->MasterReceiveGeneProcesses();
		return 1;
	}

	private:

	MixOmegaModel* model;
};

void MixOmegaModel::ComputeKalmanSuffStat(ConditionalExternalKalmanMultiVariateTreeProcess* kalman)	{

	double** W = kalman->GetW();
	for (int i=0; i<Ncont+2; i++)	{
		for (int j=0; j<Ncont+2; j++)	{
			W[i][j] = 0;
			for (int gene=0; gene<Ngene; gene++)	{
				W[i][j] += genesyntau->GetVal(gene)->val() * GetSynRegCoeff(gene,i) * GetSynRegCoeff(gene,j);
				W[i][j] += genetau->GetVal(gene)->val() * GetRegCoeff(gene,i) * GetRegCoeff(gene,j);
			}
		}
	}

	RecursiveComputeKalmanSuffStat(GetRoot(), kalman);

}

void MixOmegaModel::RecursiveComputeKalmanSuffStat(const Link* from, ConditionalExternalKalmanMultiVariateTreeProcess* kalman)	{

	double** tmp = kalman->GetKa(from);
	for (int i=0; i<Ncont+2; i++)	{
		tmp[i][0] = 0;
	}
	if (! from->isRoot())	{
		for (int gene=0; gene<Ngene; gene++)	{
			double syntemp = genesyntau->GetVal(gene)->val() * (genesynprocess[gene]->GetNodeVal(from->GetNode())->val() - genesynprocess[gene]->GetNodeVal(from->Out()->GetNode())->val());
			double temp = genetau->GetVal(gene)->val() * (geneprocess[gene]->GetNodeVal(from->GetNode())->val() - geneprocess[gene]->GetNodeVal(from->Out()->GetNode())->val());
			for (int i=0; i<Ncont+2; i++)	{
				tmp[i][0] += GetSynRegCoeff(gene,i) * syntemp;
				tmp[i][0] += GetRegCoeff(gene,i) * temp;
			}
		}
	}
	for (const Link* link=from->Next(); link!=from; link=link->Next())	{
		RecursiveComputeKalmanSuffStat(link->Out(),kalman);
	}

}

class ConditionalExternalMultiVariateKalmanMove : public MCUpdate, public Mnode {

	ConditionalExternalKalmanMultiVariateTreeProcess* tree;
	MixOmegaModel* model;
	int* index;

	public:

	ConditionalExternalMultiVariateKalmanMove(ConditionalExternalKalmanMultiVariateTreeProcess* intree, MixOmegaModel* inmodel, int* inindex)	{
		tree = intree;
		model = inmodel;
		index = new int[tree->GetDim()];
		for (int i=0; i<tree->GetDim(); i++)	{
			index[i] = inindex[i];
		}
		tree->RecursiveRegister(this,tree->GetRoot());
	}

	double Move(double tuning_modulator){
		Corrupt(false);
		tree->KalmanCreate(index);
		model->ComputeKalmanSuffStat(tree);
		tree->KalmanInit();
		tree->KalmanMove();
		tree->KalmanDelete();
		Update();
		return 1;
	}
};

class GammaIIDArrayHyperMove : public MCUpdate, public Mnode	{

	GammaIIDArray* array;
	Jeffreys* shape;
	Jeffreys* scale;
	double tuning;
	int n;

	public:

	GammaIIDArrayHyperMove(GammaIIDArray* inarray, Jeffreys* inshape, Jeffreys* inscale, double intuning, int inn)	{
		array = inarray;
		shape = inshape;
		scale = inscale;
		tuning = intuning;
		n = inn;
		shape->Register(this);
		scale->Register(this);
		// array->RegisterArray(this);
	}

	double Move(double tuning_mod = 1)	{

		Corrupt(false);

		ComputeSuffStat();

		// make series of MH moves
		double acc = 0;
		double tot = 0;
		for (int rep=0; rep<n; rep++)	{
			acc += moveshape(tuning_mod * tuning);
			tot++;
			acc += moveshape(tuning_mod * tuning * 0.3);
			tot++;
			acc += moveshape(tuning_mod * tuning * 0.1);
			tot++;
			acc += movescale(tuning_mod * tuning);
			tot++;
			acc += movescale(tuning_mod * tuning * 0.3);
			tot++;
			acc += movescale(tuning_mod * tuning * 0.1);
			tot++;
		}

		// clean up
		Update();
		return acc / tot;
	}

	void ComputeSuffStat()	{
		sum = 0;
		sumlog = 0;
		for (int i=0; i<array->GetSize(); i++)	{
			double tmp = array->GetVal(i)->val();
			sum += tmp;
			sumlog += log(tmp);
		}
	}

	double logprob()	{
		int N = array->GetSize();
		double alpha = shape->val();
		double beta = scale->val();
		// jeffreys prior for alpha and beta included here
		return alpha * N * log(beta) - N * Random::logGamma(alpha) + (alpha - 1) * sumlog - beta * sum - log(alpha) - log(beta);
	}

	int moveshape(double tuning)	{

		double bk = shape->val();
		double m = tuning * (Random::Uniform() - 0.5);
		double e = exp(m);
		double logratio = -logprob();
		double newval = bk*e;
		shape->setval(newval);
		logratio += logprob();
		logratio += m;
		int acc = ((newval > minjeff) && (newval < maxjeff) && (log(Random::Uniform()) < logratio));
		if (! acc)	{
			shape->setval(bk);
		}
		return acc;
	}

	int movescale(double tuning)	{

		double bk = scale->val();
		double m = tuning * (Random::Uniform() - 0.5);
		double e = exp(m);
		double logratio = -logprob();
		double newval = bk * e;
		scale->setval(newval);
		logratio += logprob();
		logratio += m;
		int acc = ((newval > minjeff) && (newval < maxjeff) && (log(Random::Uniform()) < logratio));
		if (! acc)	{
			scale->setval(bk);
		}
		return acc;
	}

	double sum;
	double sumlog;
};

class DirichletIIDArrayHyperMove : public MCUpdate, public Mnode	{

	DirichletIIDArray* array;
	Jeffreys* concentration;
	Dirichlet* center;
	int dim;
	double* sumlog;
	double* bkcenter;
	int* indices;
	int nrep;

	public:

	DirichletIIDArrayHyperMove(DirichletIIDArray* inarray, Dirichlet* incenter, Jeffreys* inconcentration, int innrep)	{
		array = inarray;
		center = incenter;
		concentration = inconcentration;
		nrep = innrep;
		concentration->Register(this);
		center->Register(this);
		// array->RegisterArray(this);
		dim = center->GetDim();
		sumlog = new double[dim];
		indices = new int[dim];
		bkcenter = new double[dim];
	}

	double Move(double tuning_mod = 1)	{

		Corrupt(false);

		ComputeSuffStat();

		// make series of MH moves
		double acc = 0;
		double tot = 0;
		for (int rep=0; rep<nrep; rep++)	{
			acc += moveconcentration(3);
			tot++;
			acc += moveconcentration(1);
			tot++;
			acc += moveconcentration(0.3);
			tot++;
			acc += moveconcentration(0.1);
			tot++;
			acc += movecenter(0.3,1);
			tot++;
			acc += movecenter(0.1,1);
			tot++;
			acc += movecenter(0.1,2);
			tot++;
			acc += movecenter(0.03,2);
			tot++;
			acc += movecenter(0.01,5);
			tot++;
		}

		// clean up
		double total = 0;
		for (int k=0; k<dim; k++)	{
			total += (*center)[k];
		}
		for (int k=0; k<dim; k++)	{
			(*center)[k] /= total;
		}
		Update();
		return acc / tot;
	}

	void ComputeSuffStat()	{
		for (int k=0; k<dim; k++)	{
			sumlog[k] = 0;
		}
		for (int i=0; i<array->GetSize(); i++)	{
			for (int k=0; k<dim; k++)	{
				sumlog[k] += log((*array->GetVal(i))[k]);
			}
		}
	}

	double logprob()	{

		int N = array->GetSize();

		double ret = N * Random::logGamma(concentration->val());
		for (int k=0; k<dim; k++)	{
			ret -= N * Random::logGamma(concentration->val() * (*center)[k]);
			ret += (concentration->val() * (*center)[k] - 1) * sumlog[k];
		}
		// jeffreys prior on concentration
		ret -= log(concentration->val());
		return ret;
	}

	int moveconcentration(double tuning)	{

		double bk = concentration->val();
		double m = tuning * (Random::Uniform() - 0.5);
		double e = exp(m);
		double logratio = -logprob();
		double newval = bk*e;
		concentration->setval(newval);
		logratio += logprob();
		logratio += m;
		int acc = ((newval > minjeff) && (newval < maxjeff) && (log(Random::Uniform()) < logratio));
		if (! acc)	{
			concentration->setval(bk);
		}
		return acc;
	}

	int movecenter(double tuning, int n)	{

		double logratio = -logprob();

		for (int k=0; k<dim; k++)	{
			bkcenter[k] = (*center)[k];
		}
		if (2*n > dim)	{
			n = dim / 2;
		}
		Random::DrawFromUrn(indices,2*n,dim);
		for (int i=0; i<n; i++)	{
			int i1 = indices[2*i];
			int i2 = indices[2*i+1];
			double tot = (*center)[i1] + (*center)[i2];
			double x = (*center)[i1];

			// double h = tuning * (Random::Uniform() - 0.5);
			double h = tot * tuning * (Random::Uniform() - 0.5);
			x += h;
			while ((x<0) || (x>tot))	{
				if (x<0)	{
					x = -x;
				}
				if (x>tot)	{
					x = 2*tot - x;
				}
			}
			(*center)[i1] = x;
			(*center)[i2] = tot - x;
		}

		logratio += logprob();
		int acc = (log(Random::Uniform()) < logratio);
		if (! acc)	{
			for (int k=0; k<dim; k++)	{
				(*center)[k] = bkcenter[k];
			}
		}
		return acc;
	}
};

class NormalIIDArrayHyperMove : public MCUpdate, public Mnode	{

	NormalIIDArray* array;
	Uniform* mean;
	Jeffreys* var;

	public:

	NormalIIDArrayHyperMove(NormalIIDArray* inarray, Uniform* inmean, Jeffreys* invar)	{
		array = inarray;
		mean = inmean;
		var = invar;
		mean->Register(this);
		var->Register(this);
		// array->RegisterArray(this);
	}

	double Move(double tuning = 1)	{

		Corrupt(false);
		MoveMean();
		MoveVar();
		Update();
		return 1;
	}

	void MoveMean()	{

		double m1 = 0;
		for (int i=0; i<array->GetSize(); i++)	{
			m1 += array->GetVal(i)->val();
		}
		m1 /= array->GetSize();
		double x = Random::sNormal() * sqrt(var->val() / array->GetSize()) + m1;
		if ((x > minuni) && (x < maxuni))	{
			mean->setval(x);
		}
	}

	void MoveVar()	{

		double m2 = 0;
		double m = mean->val();
		for (int i=0; i<array->GetSize(); i++)	{
			double tmp = array->GetVal(i)->val() - m;
			m2 += tmp * tmp;
		}
		double prec = Random::Gamma(0.5 * array->GetSize(), 0.5 * m2);
		if ((prec > minjeff) && (prec < maxjeff))	{
			var->setval(1.0 / prec);
		}
	}
};

template<class T> class ArrayMove : public MCUpdate	{

	T* array;
	int N;
	int* alloc;
	int myid;
	double tuning;

	public:

	ArrayMove(T* inarray, int inN, int* inalloc, int inmyid, double intuning)	{
		array = inarray;
		N = inN;
		alloc = inalloc;
		myid = inmyid;
		tuning = intuning;
	}

	double Move(double tuning_mod = 1)	{
		double tot = 0;
		double acc = 0;
		for (int i=0; i<N; i++)	{
			if (alloc[i] == myid)	{
				acc += array->GetVal(i)->Move(tuning_mod * tuning);
				tot ++;
			}
		}
		return acc / tot;
	}
};

template<class T> class BidimArrayMove : public MCUpdate	{

	BidimArray<T>* array;
	int N;
	int* alloc;
	int myid;
	double tuning;

	public:

	BidimArrayMove(BidimArray<T>* inarray, int inN, int* inalloc, int inmyid, double intuning)	{
		array = inarray;
		N = inN;
		alloc = inalloc;
		myid = inmyid;
		tuning = intuning;
	}

	double Move(double tuning_mod = 1)	{
		return array->ArrayMove(alloc,myid,tuning_mod * tuning);
	}
};

void MixOmegaModel::MasterMakeScheduler()	{

	if (! fixomega)	{
		scheduler.Register(new MasterSendParametersMove(this),1,"send param");
		scheduler.Register(new MasterReceiveGeneProcessesMove(this),1,"receive gene processes");
	}

	for (int i=0; i<nrep; i++)	{
		if (! clamptree)	{
			if ((chronoprior >= 1) && (chronoprior <= 3))	{
				scheduler.Register(new SimpleMove(Chi,1),10,"bd hyper");
				scheduler.Register(new SimpleMove(Chi,0.1),10,"bd hyper");
				scheduler.Register(new SimpleMove(Chi2,1),10,"bd hyper");
				scheduler.Register(new SimpleMove(Chi2,0.1),10,"bd hyper");
			}
			scheduler.Register(new SimpleMove(chronogram,1),10,"chrono");
			scheduler.Register(new SimpleMove(chronogram,0.1),10,"chrono");
			scheduler.Register(new SimpleMove(chronogram,0.01),10,"chrono");
		}

		if (isCalibrated() && (! clamptree))	{
			scheduler.Register(new SimpleMove(GetCalibratedChronogram()->GetScale(),1),10,"root age");
			scheduler.Register(new SimpleMove(GetCalibratedChronogram()->GetScale(),0.1),10,"root age");
			scheduler.Register(new SimpleMove(GetCalibratedChronogram()->GetScale(),0.01),10,"root age");
			scheduler.Register(new SimpleMove(GetCalibratedChronogram()->GetScale(),0.001),10,"root age");
		}

		int* index = new int[Ncont+2];
		for (int j=0; j<Ncont; j++)	{
			index[j] = 0;
		}

		index[Ncont] = 1;
		index[Ncont+1] = 1;
		scheduler.Register(new ConditionalExternalMultiVariateKalmanMove(GetConditionalExternalKalmanMultiVariateTreeProcess(),this,index),1,"kalman");

		if (! contml)	{
			index[Ncont] = 0;
			index[Ncont+1] = 0;
			for (int j=0; j<Ncont; j++)	{
				index[j] = 1;
				scheduler.Register(new ConditionalExternalMultiVariateKalmanMove(GetConditionalExternalKalmanMultiVariateTreeProcess(),this,index),1,"kalman");
				index[j] = 0;
			}
		}

		delete[] index;

		scheduler.Register(new ConjugateMultiVariateMove(GetConjugateInverseWishart(),GetConjugateMultiVariateTreeProcess(),1,0),1,"conjugate sigma - process");

		scheduler.Register(new SimpleMove(DiagArray,3),100,"diag array");
		scheduler.Register(new SimpleMove(DiagArray,1),100,"diag array");
		scheduler.Register(new SimpleMove(DiagArray,0.1),100,"diag array");

		// CONJUGATE
		scheduler.Register(new GammaIIDArrayHyperMove(genesyntau,syntaushape,syntauscale,1,100),1,"suffstat syn tau shape and scale");
		scheduler.Register(new GammaIIDArrayHyperMove(genetau,taushape,tauscale,1,100),1,"suffstat tau shape and scale");

		/*
		scheduler.Register(new SimpleMove(syntaushape,1),10,"syn tau shape");
		scheduler.Register(new SimpleMove(syntaushape,0.1),10,"syn tau shape");
		scheduler.Register(new SimpleMove(syntaushape,0.03),10,"syn tau shape");

		scheduler.Register(new SimpleMove(syntauscale,1),10,"syn tau scale");
		scheduler.Register(new SimpleMove(syntauscale,0.1),10,"syn tau scale");
		scheduler.Register(new SimpleMove(syntauscale,0.03),10,"syn tau scale");

		scheduler.Register(new SimpleMove(taushape,1),10,"tau shape");
		scheduler.Register(new SimpleMove(taushape,0.1),10,"tau shape");
		scheduler.Register(new SimpleMove(taushape,0.03),10,"tau shape");

		scheduler.Register(new SimpleMove(tauscale,1),10,"tau scale");
		scheduler.Register(new SimpleMove(tauscale,0.1),10,"tau scale");
		scheduler.Register(new SimpleMove(tauscale,0.03),10,"tau scale");
		*/

		if (! fixomega)	{
			scheduler.Register(new NormalIIDArrayHyperMove(generootsyn,rootsynmean,rootsynvar),1,"conjugate root syn mean and var");
			scheduler.Register(new NormalIIDArrayHyperMove(generootomega,rootomegamean,rootomegavar),1,"conjugate root omega mean and var");
		}

		/*
		if (! fixomega)	{
			scheduler.Register(new SimpleMove(rootsynmean,3),10,"root syn mean");
			scheduler.Register(new SimpleMove(rootsynmean,1),10,"root syn mean");
			scheduler.Register(new SimpleMove(rootsynmean,0.1),10,"root syn mean");

			scheduler.Register(new SimpleMove(rootsynvar,3),10,"root syn var");
			scheduler.Register(new SimpleMove(rootsynvar,1),10,"root syn var");
			scheduler.Register(new SimpleMove(rootsynvar,0.1),10,"root syn var");
		}


		if (! fixomega)	{
			scheduler.Register(new SimpleMove(rootomegamean,3),10,"root omega mean");
			scheduler.Register(new SimpleMove(rootomegamean,1),10,"root omega mean");
			scheduler.Register(new SimpleMove(rootomegamean,0.1),10,"root omega mean");

			scheduler.Register(new SimpleMove(rootomegavar,3),10,"root omega var");
			scheduler.Register(new SimpleMove(rootomegavar,1),10,"root omega var");
			scheduler.Register(new SimpleMove(rootomegavar,0.1),10,"root omega var");
		}
		*/

		if (! fixomega)	{
			scheduler.Register(new DirichletIIDArrayHyperMove(generelrate,relratecenter,relrateconcentration,100),1,"suffstat hyper relrate");
			scheduler.Register(new DirichletIIDArrayHyperMove(genestationary,stationarycenter,stationaryconcentration,100),1,"suffstat hyper stationary");
		}

		/*
		if (! fixomega)	{
			scheduler.Register(new ProfileMove(relratecenter,0.1,1),10,"relrates center");
			scheduler.Register(new ProfileMove(relratecenter,0.03,2),10,"relrates center");
			scheduler.Register(new SimpleMove(relratecenter,0.01),10,"relrates center");

			scheduler.Register(new ProfileMove(stationarycenter,0.01,2),10,"stat4 center");
			scheduler.Register(new ProfileMove(stationarycenter,0.03,2),10,"stat4 center");
			scheduler.Register(new ProfileMove(stationarycenter,0.01,5),10,"stat10 center");
			scheduler.Register(new SimpleMove(stationarycenter,0.001),10,"stat center");

			scheduler.Register(new SimpleMove(relrateconcentration,10),10,"relrates conc");
			scheduler.Register(new SimpleMove(relrateconcentration,1),10,"relrates conc");
			scheduler.Register(new SimpleMove(relrateconcentration,0.1),10,"relrates conc");
			scheduler.Register(new SimpleMove(relrateconcentration,0.01),10,"relrates conc");

			scheduler.Register(new SimpleMove(stationaryconcentration,10),10,"stat conc");
			scheduler.Register(new SimpleMove(stationaryconcentration,1),10,"stat conc");
			scheduler.Register(new SimpleMove(stationaryconcentration,0.1),10,"stat conc");
			scheduler.Register(new SimpleMove(stationaryconcentration,0.01),10,"stat conc");
		}
		*/
	}

}

class RootProcessCompensatoryMove : public MCUpdate, public Mnode	{

	LinRegNormalProcess* process;
	Normal* rootval;
	double tuning;

	public:

	RootProcessCompensatoryMove(LinRegNormalProcess* inprocess, Normal* inrootval, double intuning)	{
		process = inprocess;
		rootval = inrootval;
		tuning = intuning;
		rootval->Register(this);
		process->RecursiveRegister(this,process->GetRoot());
	}

	double Move(double tuning_mod = 1)	{

		Corrupt(true);
		double m = tuning_mod * tuning * (Random::Uniform() - 0.5);
		rootval->setval(rootval->val() + m);
		process->Translation(-m);
		double logprob = Update();
		int acc = (log(Random::Uniform()) < logprob);
		if (! acc)	{
			Corrupt(false);
			Restore();
		}
		return (double) acc;
	}
};

void MixOmegaModel::SlaveMakeScheduler()	{


	scheduler.Register(new SlaveReceiveParametersMove(this),1,"receive params");

	if (! fixomega)	{
		for (int gene=0; gene<Ngene; gene++)	{
			if (genealloc[gene] == myid)	{
				if (conjpath)	{
					if (! clampsuffstat)	{
						scheduler.Register(new DSemiConjugateMappingMove(phyloprocess[gene],pathconjtree[gene],mappingfreq),1,"mapping + sufficient stat",mappingevery);
					}
				}
				else	{
					if (phyloprocess)	{
						scheduler.Register(new SimpleMove(phyloprocess[gene],mappingfreq),1,"mapping",mappingevery);
					}
				}
				scheduler.Register(new SimpleMove(genesynprocess[gene],3),10,"gene syn process");
				scheduler.Register(new SimpleMove(genesynprocess[gene],1),10,"gene syn process");
				scheduler.Register(new SimpleMove(genesynprocess[gene],0.1),10,"gene syn process");

				scheduler.Register(new RootProcessCompensatoryMove(genesynprocess[gene],generootsyn->GetNormal(gene),1),10,"gene syn process root comp move");
				scheduler.Register(new RootProcessCompensatoryMove(genesynprocess[gene],generootsyn->GetNormal(gene),0.1),10,"gene syn process root comp move");

				scheduler.Register(new SimpleMove(geneprocess[gene],3),10,"gene process");
				scheduler.Register(new SimpleMove(geneprocess[gene],1),10,"gene process");
				scheduler.Register(new SimpleMove(geneprocess[gene],0.1),10,"gene process");

				scheduler.Register(new RootProcessCompensatoryMove(geneprocess[gene],generootomega->GetNormal(gene),1),10,"gene process root comp move");
				scheduler.Register(new RootProcessCompensatoryMove(geneprocess[gene],generootomega->GetNormal(gene),0.1),10,"gene process root comp move");

				scheduler.Register(new ProfileMove(generelrate->GetVal(gene),0.1,1),10,"relrates");
				scheduler.Register(new ProfileMove(generelrate->GetVal(gene),0.03,2),10,"relrates");
				scheduler.Register(new SimpleMove(generelrate->GetVal(gene),0.01),10,"relrates");

				scheduler.Register(new ProfileMove(genestationary->GetVal(gene),0.01,2),10,"stat4");
				scheduler.Register(new ProfileMove(genestationary->GetVal(gene),0.03,2),10,"stat4");
				scheduler.Register(new ProfileMove(genestationary->GetVal(gene),0.01,5),10,"stat10");
				scheduler.Register(new SimpleMove(genestationary->GetVal(gene),0.001),10,"stat");
			}
		}

		scheduler.Register(new ArrayMove<NormalIIDArray>(generootsyn,Ngene,genealloc,myid,1),10,"gene root syn");
		scheduler.Register(new ArrayMove<NormalIIDArray>(generootsyn,Ngene,genealloc,myid,0.1),10,"gene root syn");

		scheduler.Register(new ArrayMove<NormalIIDArray>(generootomega,Ngene,genealloc,myid,1),10,"gene root omega");
		scheduler.Register(new ArrayMove<NormalIIDArray>(generootomega,Ngene,genealloc,myid,0.1),10,"gene root omega");

	}

	// for (int rep=0; rep<nrep; rep++)	{

	scheduler.Register(new BidimArrayMove<Real>(synregarray,Ngene,genealloc,myid,3),10,"syn reg");
	scheduler.Register(new BidimArrayMove<Real>(synregarray,Ngene,genealloc,myid,1),10,"syn reg");
	scheduler.Register(new BidimArrayMove<Real>(synregarray,Ngene,genealloc,myid,0.1),10,"syn reg");

	scheduler.Register(new ArrayMove<GammaIIDArray>(genesyntau,Ngene,genealloc,myid,3),10,"genesyntau");
	scheduler.Register(new ArrayMove<GammaIIDArray>(genesyntau,Ngene,genealloc,myid,1),10,"genesyntau");
	scheduler.Register(new ArrayMove<GammaIIDArray>(genesyntau,Ngene,genealloc,myid,0.1),10,"genesyntau");

	scheduler.Register(new BidimArrayMove<Real>(regarray,Ngene,genealloc,myid,3),10,"reg");
	scheduler.Register(new BidimArrayMove<Real>(regarray,Ngene,genealloc,myid,1),10,"reg");
	scheduler.Register(new BidimArrayMove<Real>(regarray,Ngene,genealloc,myid,0.1),10,"reg");

	scheduler.Register(new ArrayMove<GammaIIDArray>(genetau,Ngene,genealloc,myid,3),10,"genetau");
	scheduler.Register(new ArrayMove<GammaIIDArray>(genetau,Ngene,genealloc,myid,1),10,"genetau");
	scheduler.Register(new ArrayMove<GammaIIDArray>(genetau,Ngene,genealloc,myid,0.1),10,"genetau");
	// }

	scheduler.Register(new SlaveSendGeneProcessesMove(this),1,"send gene processes");

}

#endif
