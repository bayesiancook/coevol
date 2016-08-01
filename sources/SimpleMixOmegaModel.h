
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

#include "LinRegOmega.h"

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
	ContinuousData* contdata;
	TaxonSet* taxonset;

	// number of columns
	int Ngene;
	int* Nsite;
	// number of states (4 for nucleic acids, 20 for amino-acids. 61 for codons)
	// int Nstate;

	int Ncont;

	// ---------
	// the random variables of the model
	// ---------

	Const<Real>* Zero;
	Const<PosReal>* One;

	Const<PosReal>* MeanChi;
	Const<PosReal>* MeanChi2;
	Exponential* Chi;
	Exponential* Chi2;

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
	LengthTree* lengthtree;

	Jeffreys* synsigma;
	LogNormalTreeProcess* synratetree;

	JeffreysIIDArray* DiagArray;
	SigmaZero* sigmaZero;
	Rvar<CovMatrix>* sigma;

	Const<RealVector>* rootmean;
	Const<PosRealVector>* rootvar;

	MultiVariateTreeProcess* process;

	Jeffreys* taushape;
	Jeffreys* tauscale;
	GammaIIDArray* genetau;

	Uniform* rootomegamean;
	Jeffreys* rootomegavar;
	NormalIIDArray* generootomega;

	Jeffreys** jeffregvar;
	Var<Real>** regmean;
	Var<PosReal>** regvar;

	BidimIIDNormal* regarray;

	LinRegNormalProcess** geneprocess;
	MeanExpTree** geneomegatree;

	Dirichlet* relrate;
	Dirichlet* stationary;
	GTRRandomSubMatrixWithNormRates* nucmatrix;

	MatrixTree** genecodonmatrixtree;

	// phylo process
	PathConjugateTree** pathconjtree;
	PhyloProcess** phyloprocess;

	bool meanexp;

	bool conjpath;
	bool normalise;

	int nrep;

	int df;

	int clampsuffstat;
	string suffstatfile;

	bool clamptree;
	bool clampdiag;
	bool clampreg;

	bool priorsampling;

	double mappingfreq;

	public:

	CodonSequenceAlignment* GetData(int gene)	{
		return codondata[gene];
	}

	bool isCalibrated()	{
		return iscalib;
	}

	bool isConstrained()	{
		return clamptree;
	}

	MixOmegaModel(string datafile, string treefile, string contdatafile, string calibfile, double rootage, double rootstdev, int inchronoprior, double inmeanchi, double inmeanchi2, int indf, bool inclampdiag, int inconjpath, bool inclamptree, bool inmeanexp, bool innormalise, int innrep, string insuffstatfile, string rootfile, bool inclampreg, double inmappingfreq, GeneticCodeType type = Universal, bool sample=true)	{

		clampreg = inclampreg;
		mappingfreq = inmappingfreq;

		clamptree = inclamptree;

		suffstatfile = insuffstatfile;
		clampsuffstat = (suffstatfile != "None");

		df = indf;

		chronoprior = inchronoprior;
		iscalib = false;
		meanchi = inmeanchi;
		meanchi2 = inmeanchi2;

		clampdiag = inclampdiag;
		meanexp = inmeanexp;

		// get data from file

		ifstream is(datafile.c_str());
		is >> Ngene;
		Nsite = new int[Ngene];
		nucdata = new FileSequenceAlignment*[Ngene];
		codondata = new CodonSequenceAlignment*[Ngene];
		for (int gene=0; gene<Ngene; gene++)	{
			string filename;
			is >> filename;
			nucdata[gene] = new FileSequenceAlignment(filename);
			codondata[gene] = new CodonSequenceAlignment(nucdata[gene],type);
			Nsite[gene] = GetData(gene)->GetNsite();
			// Nstate = GetData(gene)->GetNstate();
			cerr << filename << '\t' << Nsite[gene] << '\n';
		}

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
		if (nrep == 0)	{
			nrep = conjpath ? 10 : 1;
		}
		cerr << "nrep : " << nrep << '\n';
		cerr << conjpath << '\n';
		normalise = innormalise;
		if (clampsuffstat && ! conjpath)	{
			cerr << "error : suffstat requires conjugate path sampling\n";
			exit(1);
		}

		// asssume all datasets have same taxonset!
		// or at least first dataset
		taxonset = nucdata[0]->GetTaxonSet();

		// get tree from file (newick format)
		tree = new Tree(treefile);
		// check whether tree and data fit together
		tree->RegisterWith(taxonset);

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

		cerr << "tree and data ok\n";
		cerr << '\n';

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
		}
		else	{
			chronogram = new Chronogram(tree,One);
		}

		if (clamptree)	{
			chronogram->Clamp();
		}
		lengthtree = chronogram;

		double minjeff = 0.001;
		double maxjeff = 1000;

		synsigma = new Jeffreys(minjeff,maxjeff,Zero);
		synsigma->setval(1);
		synratetree = new LogNormalTreeProcess(lengthtree,synsigma,INTEGRAL);
		synratetree->Reset();

		DiagArray = new JeffreysIIDArray(Ncont+1,minjeff,maxjeff,Zero);
		DiagArray->setval(1.0);

		sigmaZero = new SigmaZero(DiagArray);

		if (clampdiag)	{
			sigma = new DiagonalCovMatrix(sigmaZero,Ncont+1+df);
		}
		else	{
			sigma = new InverseWishartMatrix(sigmaZero,Ncont+1+df);
		}

		rootmean = 0;
		rootvar = 0;
		if (rootfile != "None")	{
			rootmean = new Const<RealVector>(RealVector(Ncont+1));
			rootvar = new Const<PosRealVector>(PosRealVector(Ncont+1));
			(*rootmean)[0] = 0;
			(*rootvar)[0] = 0;
			ifstream is(rootfile.c_str());
			for (int i=0; i<Ncont; i++)	{
				double mean, var;
				is >> mean >> var;
				(*rootmean)[i+1] = mean;
				(*rootvar)[i+1] = var;
			}
		}

		if (clampdiag)	{
			process = new MultiVariateTreeProcess(sigma,lengthtree,0,0,rootmean,rootvar);
		}
		else	{
			process = new MultiVariateTreeProcess(sigma,lengthtree,0,0,rootmean, rootvar);
		}

		process->Reset();

		// here make the regcoef array

		cerr << "reg mean and var\n";
		jeffregvar = new Jeffreys*[Ncont];
		regmean = new Var<Real>*[Ncont+1];
		regvar = new Var<PosReal>*[Ncont+1];
		regmean[0] = Zero;
		regvar[0] = One;
		for (int i=0; i<Ncont; i++)	{
			jeffregvar[i] = new Jeffreys(minjeff,maxjeff,Zero);
			regvar[i+1] = jeffregvar[i];
			regmean[i+1] = Zero;
		}

		cerr << "regarray\n";
		cerr << Ngene << '\t' << Ncont + 1 << '\n';
		regarray = new BidimIIDNormal(Ngene,Ncont+1,regmean,regvar);
		regarray->SetAt(1.0,0);
		for (int i=0; i<Ncont; i++)	{
			cerr << "set \n";
			regarray->SetAt(0,i+1);
		}
		if (clampreg)	{
			cerr << "clamp\n";
			for (int i=0; i<Ncont; i++)	{
				regarray->Clamp(i+1);
			}
		}
		// clamp first colum at one

		cerr << "tau\n";
		taushape = new Jeffreys(minjeff,maxjeff,Zero);
		tauscale = new Jeffreys(minjeff,maxjeff,Zero);
		taushape->setval(10.0);
		tauscale->setval(100.0);
		genetau = new GammaIIDArray(Ngene,taushape,tauscale);

		cerr << "gene process\n";

		geneprocess = new LinRegNormalProcess*[Ngene];
		for (int gene=0; gene<Ngene; gene++)	{
			geneprocess[gene] = new LinRegNormalProcess(lengthtree,process,regarray,gene,genetau->GetVal(gene));
			geneprocess[gene]->Reset();
		}

		double minuni = -1000;
		double maxuni = 1000;

		cerr << "root vals\n";
		rootomegamean = new Uniform(minuni,maxuni,Zero);
		rootomegamean->setval(-1.0);
		rootomegavar = new Jeffreys(minjeff,maxjeff,Zero);
		rootomegavar->setval(0.1);
		generootomega = new NormalIIDArray(Ngene,rootomegamean,rootomegavar);

		cerr << "omega trees\n";
		geneomegatree = new MeanExpTree*[Ngene];
		for (int gene=0; gene<Ngene; gene++)	{
			geneomegatree[gene] = new MeanExpTree(geneprocess[gene],lengthtree,MEAN,generootomega->GetVal(gene));
		}

		cerr << "matrix\n";
		// create a GTR nucleotide matrix
		relrate = new Dirichlet(Nnuc*(Nnuc-1)/2);
		stationary = new Dirichlet(Nnuc);
		nucmatrix = new GTRRandomSubMatrixWithNormRates(relrate,stationary,true);

		cerr << "codon matrix trees\n";
		genecodonmatrixtree = new MatrixTree*[Ngene];
		for (int gene=0; gene<Ngene; gene++)	{
			genecodonmatrixtree[gene] = new MatrixTree(codondata[gene]->GetCodonStateSpace(), nucmatrix, geneomegatree[gene], One);
		}

		cerr << "set and clamp\n";
		if (contdata)	{
			for (int i=0; i<Ncont; i++)	{
				process->SetAndClamp(contdata,i+1,i);
			}
		}

		cerr << "phyloprocess\n";
		// make substitution mappings
		if (conjpath)	{
			pathconjtree = new PathConjugateTree*[Ngene];
			phyloprocess = new PhyloProcess*[Ngene];
			for (int gene=0; gene<Ngene; gene++)	{
				pathconjtree[gene] = new MGCodonBranchMatrixPathConjugateTree(synratetree, genecodonmatrixtree[gene], GetData(gene));
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
				phyloprocess[gene] = new BranchMatrixPhyloProcess(synratetree, genecodonmatrixtree[gene], GetData(gene));
			}
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
		RootRegister(relrate);
		RootRegister(stationary);
		Register();

		cerr << "make scheduler\n";
		MakeScheduler();
		if (sample)	{
			cerr << "update\n";
			Update();
			TraceHeader(cerr);
			Trace(cerr);
		}
		cerr << "model created\n";
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

	Chronogram* GetChronogram() {
		return chronogram;
	}
	LengthTree* GetLengthTree() {return lengthtree;}

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
		return GetLogPrior() + GetLogLikelihood();
	}

	double GetLogPrior()	{
		double total = 0;

		if ((chronoprior >= 1) && (chronoprior <= 3))	{
			total += Chi->GetLogProb();
			total += Chi2->GetLogProb();
		}
		total += chronogram->GetLogProb();

		total += synsigma->GetLogProb();
		total += synratetree->GetLogProb();

		total += DiagArray->GetLogProb();
		total += sigma->GetLogProb();
		total += process->GetLogProb();

		for (int i=0; i<Ncont; i++)	{
			total += jeffregvar[i]->GetLogProb();
		}

		total += regarray->GetLogProb();

		total += taushape->GetLogProb();
		total += tauscale->GetLogProb();
		total += genetau->GetLogProb();

		total += rootomegamean->GetLogProb();
		total += rootomegavar->GetLogProb();
		total += generootomega->GetLogProb();

		for (int gene=0; gene<Ngene; gene++)	{
			total += geneprocess[gene]->GetLogProb();
		}

		total += relrate->GetLogProb();
		total += stationary->GetLogProb();

		return total;
	}

	double GetLogLikelihood()	{
		double ret = 0;
		if (priorsampling)	{
			return 0;
		}
		else if (clampsuffstat)	{
			for (int gene=0; gene<Ngene; gene++)	{
				ret += pathconjtree[gene]->GetLogProb();
			}
		}
		else	{
			for (int gene=0; gene<Ngene; gene++)	{
				ret += phyloprocess[gene]->GetLogProb();
			}
		}
		return ret;
	}

	virtual void MakeScheduler()	{

		if (conjpath)	{
			if (! clampsuffstat)	{
				for (int gene=0; gene<Ngene; gene++)	{
					scheduler.Register(new DSemiConjugateMappingMove(phyloprocess[gene],pathconjtree[gene],mappingfreq),1,"mapping + sufficient stat");
				}
			}
		}
		else	{
			if (phyloprocess)	{
				for (int gene=0; gene<Ngene; gene++)	{
					scheduler.Register(new SimpleMove(phyloprocess[gene],1),1,"mapping");
				}
			}
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

			scheduler.Register(new SimpleMove(synratetree,10),10,"synratetree");
			scheduler.Register(new SimpleMove(synratetree,1),10,"synratetree");
			scheduler.Register(new SimpleMove(synratetree,0.1),10,"synratetree");
			scheduler.Register(new SimpleMove(synratetree,0.01),10,"synratetree");

			scheduler.Register(new SimpleMove(synsigma,10),100,"syn sigma");
			scheduler.Register(new SimpleMove(synsigma,1),100,"syn sigma");
			scheduler.Register(new SimpleMove(synsigma,0.1),100,"syn sigma");
			scheduler.Register(new SimpleMove(synsigma,0.01),100,"syn sigma");

			scheduler.Register(new SimpleMove(process,10),10,"multinormal");
			scheduler.Register(new SimpleMove(process,1),10,"multinormal");
			scheduler.Register(new SimpleMove(process,0.1),10,"multinormal");
			scheduler.Register(new SimpleMove(process,0.01),10,"multinormal");

			scheduler.Register(new SimpleMove(sigma,10),100,"sigma");
			scheduler.Register(new SimpleMove(sigma,1),100,"sigma");
			scheduler.Register(new SimpleMove(sigma,0.1),100,"sigma");
			scheduler.Register(new SimpleMove(sigma,0.01),100,"sigma");

			scheduler.Register(new SimpleMove(DiagArray,10),10,"theta");
			scheduler.Register(new SimpleMove(DiagArray,1),10,"theta");
			scheduler.Register(new SimpleMove(DiagArray,0.1),10,"theta");

			scheduler.Register(new SimpleMove(regarray,10),100,"reg");
			scheduler.Register(new SimpleMove(regarray,1),100,"reg");
			scheduler.Register(new SimpleMove(regarray,0.1),100,"reg");
			scheduler.Register(new SimpleMove(regarray,0.01),100,"reg");

			for (int i=0; i<Ncont; i++)	{
				scheduler.Register(new SimpleMove(jeffregvar[i],10),100,"reg var");
				scheduler.Register(new SimpleMove(jeffregvar[i],1),100,"reg var");
				scheduler.Register(new SimpleMove(jeffregvar[i],0.1),100,"reg var");
				scheduler.Register(new SimpleMove(jeffregvar[i],0.01),100,"reg var");
			}

			scheduler.Register(new SimpleMove(genetau,10),100,"genetau");
			scheduler.Register(new SimpleMove(genetau,1),100,"genetau");
			scheduler.Register(new SimpleMove(genetau,0.1),100,"genetau");
			scheduler.Register(new SimpleMove(genetau,0.01),100,"genetau");

			scheduler.Register(new SimpleMove(taushape,10),100,"tau shape");
			scheduler.Register(new SimpleMove(taushape,1),100,"tau shape");
			scheduler.Register(new SimpleMove(taushape,0.1),100,"tau shape");
			scheduler.Register(new SimpleMove(taushape,0.01),100,"tau shape");

			scheduler.Register(new SimpleMove(tauscale,10),100,"tau scale");
			scheduler.Register(new SimpleMove(tauscale,1),100,"tau scale");
			scheduler.Register(new SimpleMove(tauscale,0.1),100,"tau scale");
			scheduler.Register(new SimpleMove(tauscale,0.01),100,"tau scale");

			for (int gene=0; gene<Ngene; gene++)	{
				scheduler.Register(new SimpleMove(geneprocess[gene],1),10,"gene process");
				scheduler.Register(new SimpleMove(geneprocess[gene],0.1),10,"gene process");
				scheduler.Register(new SimpleMove(geneprocess[gene],0.01),10,"gene process");
			}

			scheduler.Register(new SimpleMove(generootomega,1),10,"gene root omega");
			scheduler.Register(new SimpleMove(generootomega,0.1),10,"gene root omega");

			scheduler.Register(new SimpleMove(rootomegamean,10),100,"root omega mean");
			scheduler.Register(new SimpleMove(rootomegamean,1),100,"root omega mean");
			scheduler.Register(new SimpleMove(rootomegamean,0.1),100,"root omega mean");

			scheduler.Register(new SimpleMove(rootomegavar,10),100,"root omega var");
			scheduler.Register(new SimpleMove(rootomegavar,1),100,"root omega var");
			scheduler.Register(new SimpleMove(rootomegavar,0.1),100,"root omega var");

			scheduler.Register(new ProfileMove(relrate,0.1,1),10,"relrates");
			scheduler.Register(new ProfileMove(relrate,0.03,2),10,"relrates");
			scheduler.Register(new SimpleMove(relrate,0.01),10,"relrates");

			scheduler.Register(new ProfileMove(stationary,0.01,2),10,"stat4");
			scheduler.Register(new ProfileMove(stationary,0.03,2),10,"stat4");
			scheduler.Register(new ProfileMove(stationary,0.01,5),10,"stat10");
			scheduler.Register(new SimpleMove(stationary,0.001),10,"stat");
		}
	}

	/*
	double Move(double tuning = 1)	{
		// Cycle(1,1,verbose,check)
		scheduler.Cycle(1,1,true,true);
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

	int GetNgene()	{
		return Ngene;
	}

	double GetRootAge()	{
		if (isCalibrated())	{
			return GetCalibratedChronogram()->GetRootAge();
			// return GetCalibratedChronogram()->GetScale()->val();
		}
		return 1;
	}

	LogNormalTreeProcess* GetSynRateTree()	{
		return synratetree;
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

	double GetMeanSynRate()	{
		return GetSynRateTree()->GetMeanRate();
	}

	double GetTotalTime()	{
		return chronogram->GetTotalTime();
	}

	double GetRegCoef(int i, int j)	{
		return regarray->GetCell(i,j)->val();
	}

	void TraceHeader(ostream& os)	{
		os << "#logprior\tlnL";
		os << "\tsynrate";
		for (int i=0; i<Ncont; i++)	{
			os << '\t' << "regvar" << i;
		}
		os << "\trootmean";
		os << "\trootvar";
		for (int k=0; k<Ncont+1; k++)	{
			for (int l=k+1; l<Ncont+1; l++)	{
				os << '\t' << "cont_" << k << '_' << l;
			}
		}
		for (int k=0; k<Ncont+1; k++)	{
			os << '\t' << "cont_" << k << '_' << k;
		}
		os << "\tsynsigma";
		os << "\trelrate";
		os << "\tstatent";
		if (isCalibrated())	{
			os << "\trootage";
		}
		if ((chronoprior >= 1) && (chronoprior <= 3))	{
			os << "\tp1\tp2";
		}
		os << "\tdim";
		for (int k=0; k<Ncont; k++)	{
			os << '\t' << "root_" << k;
		}
		os << '\n';
	}

	// writes all summary statistics on one single line
	// in the same order as that provided by the header
	void Trace(ostream& os)	{

		os << GetLogPrior() << '\t' << GetLogLikelihood();
		os << '\t' << GetMeanSynRate();
		for (int i=0; i<Ncont; i++)	{
			os << '\t' << regarray->GetVar(i+1);
		}
		os << '\t' << *rootomegamean;
		os << '\t' << *rootomegavar;

		for (int k=0; k<Ncont+1; k++)	{
			for (int l=k+1; l<Ncont+1; l++)	{
				os << '\t' << (*sigma)[k][l];
			}
		}
		for (int k=0; k<Ncont+1; k++)	{
			os << '\t' << (*sigma)[k][k];
		}

		os << '\t' << synsigma->val();
		os << '\t' << relrate->GetEntropy();
		os << '\t' << stationary->GetEntropy();

		if (isCalibrated())	{
			os << '\t' << GetRootAge();
		}
		if ((chronoprior >= 1) && (chronoprior <= 3))	{
			os << '\t' << *Chi << '\t' << *Chi2;
		}

		os << '\t' << GetProcess()->GetMultiNormal(GetTree()->GetRoot())->val();
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
		os << *synsigma << '\n';
		os << *synratetree << '\n';
		os << *DiagArray << '\n';
		os << *sigma << '\n';
		os << *process << '\n';
		for (int i=0; i<Ncont; i++)	{
			os << *jeffregvar[i] << '\n';
		}
		regarray->ToStream(os);
		os << '\n';
		os << *rootomegamean << '\n';
		os << *rootomegavar << '\n';
		os << *generootomega << '\n';
		os << *relrate << '\n';
		os << *stationary << '\n';

		for (int gene = 0; gene<Ngene; gene++)	{
			os << *geneomegatree[gene] << '\n';
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
		is >> *synsigma;
		is >> *synratetree;
		is >> *DiagArray;
		is >> *sigma;
		is >> *process;
		for (int i=0; i<Ncont; i++)	{
			is >> *jeffregvar[i];
		}
		regarray->FromStream(is);
		is >> *rootomegamean;
		is >> *rootomegavar;
		is >> *generootomega;
		is >> *relrate;
		is >> *stationary;

		for (int gene = 0; gene<Ngene; gene++)	{
			is >> *geneomegatree[gene];
		}
	}
};

#endif
