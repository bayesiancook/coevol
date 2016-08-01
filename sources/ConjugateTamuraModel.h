
#include "ConjugateMultiVariateTreeProcess.h"

#include "BaseType.h"
#include "RandomTypes.h"
#include "ProbModel.h"
#include "IID.h"
#include "PrecisionNormalTreeProcess.h"
#include "CodonSubMatrix.h"
#include "CodonSequenceAlignment.h"
#include "BranchProcess.h"
#include "OneMatrixPhyloProcess.h"
#include "ContinuousData.h"
#include "MeanExpTree.h"

#include "TamuraModel.h"

class ConjugateTamuraModel : public TamuraModel {

	public:

	ConjugateTamuraModel(string datafile, string treefile, string contdatafile, string calibfile, double rootage, double rootstdev, int inchronoprior, double inmeanchi, double inmeanchi2, double priorsigma, bool inautoregressive, int inconjpath, int contdatatype, bool inomegaratiotree, bool inclamproot, bool inmeanexp, bool innormalise, int innrep, bool sample=true, GeneticCodeType type=Universal)	{

		normalise = innormalise;

		chronoprior = inchronoprior;
		meanchi = inmeanchi;
		meanchi2 = inmeanchi2;

		autoregressive = inautoregressive;
		clampdiag = false;
		clamproot = inclamproot;
		meanexp = inmeanexp;
		omegaratiotree = inomegaratiotree;

		if (omegaratiotree)	{
			L = 4;
		}
		else	{
			L = 3;
		}

		if (inconjpath == -1)	{
			conjpath = true;
		}
		else	{
			conjpath = inconjpath;
		}

		nrep = innrep;
		if (nrep == 0)	{
			nrep = conjpath ? 30 : 1;
		}

		// get data from file
		nucdata = new FileSequenceAlignment(datafile);
		if (omegaratiotree)	{
			codondata = new CodonSequenceAlignment(nucdata, true, type);
		}
		else	{
			codondata = 0;
		}

		Nsite = GetData()->GetNsite();
		Nstate = GetData()->GetNstate();

		taxonset = nucdata->GetTaxonSet();

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

		RootAlpha = 0;
		RootBeta = 0;

		MeanChi = 0;
		MeanChi2 = 0;
		Chi = 0;
		Chi2 = 0;

		cerr << "chronoprior : " << chronoprior << '\n';
		if (calibfile != "None")	{
			double a = rootage * rootage / rootstdev / rootstdev;
			double b = rootage / rootstdev / rootstdev;
			RootAlpha = new Const<PosReal>(a);
			RootBeta = new Const<PosReal>(b);
			CalibrationSet* calibset = new FileCalibrationSet(calibfile, tree);

			if (chronoprior == 0)	{
				chronogram = new CalibratedChronogram(tree,One,RootAlpha,RootBeta,calibset);
			}
			else if (chronoprior == 1)	{
				cerr << "BD\n";
				MeanChi = new Const<PosReal>(meanchi);
				MeanChi2 = new Const<PosReal>(meanchi2);
				Chi = new Exponential(MeanChi,Exponential::MEAN);
				Chi2 = new Exponential(MeanChi2,Exponential::MEAN);
				chronogram = new BDCalibratedChronogram(tree,One,Chi,Chi2,RootAlpha,RootBeta,calibset);
			}
			else	{
				cerr << "error : chronoprior\n";
				exit(1);
			}
		}
		else	{
			chronogram = new Chronogram(tree,One);
		}

		GammaDiag = new Gamma*[Ncont+L];
		diag = new Rvar<PosReal>*[Ncont+L];
		for (int k=0; k<Ncont+L; k++)	{
			GammaDiag[k] = new Gamma(One,One);
			diag[k] = GammaDiag[k];
			GammaDiag[k]->ClampAt(priorsigma);
		}
		priorOnSigmaZero = new RvarVec(diag,Ncont+L);
		sigmaZero = new SigmaZero(priorOnSigmaZero);
		sigma = new ConjugateInverseWishart(sigmaZero,Ncont+L+1);
		sigma->SetIdentity();

		if (autoregressive)	{
			phi = new Gamma(One,One);
			mean = new IIDUniform(Ncont+L,100);
			process = new ConjugateAutoRegressiveMultiVariateTreeProcess(GetConjugateInverseWishart(),mean,phi,chronogram);
		}
		else	{
			phi = 0;
			mean = 0;
			process = new ConjugateMultiVariateTreeProcess(GetConjugateInverseWishart(),chronogram);
		}

		// process->GetMultiNormal(tree->GetRoot())->ClampAt(log(100.0),3);

		if (clamproot)	{
			process->ClampRoot();
		}

		if (contdata)	{
			for (int i=0; i<Ncont; i++)	{
				process->SetAndClamp(contdata,L+i,i,contdatatype);
			}
		}

		// process->Reset();

		// cut off to avoid numerical errors
		/*
		for (int l=0; l<L; l++)	{
			process->CutOff(1,l);
		}
		*/

		cerr << "subs process\n";

		rootgc = new Beta(One,One);
		gcsynratetree = new MeanExpTreeFromMultiVariate(process,0,INTEGRAL,false,meanexp);
		atsynratetree = new MeanExpTreeFromMultiVariate(process,1,INTEGRAL,false,meanexp);
		tstvtree = new MeanExpTreeFromMultiVariate(process,2,MEAN,false,meanexp);
		nonsynratetree = 0;
		if (omegaratiotree == 0)	{
			omegatree = 0;
		}
		else if (omegaratiotree == 1)	{
			omegatree = new MeanExpTreeFromMultiVariate(process,3,MEAN,false,meanexp);
		}
		else	{
			nonsynratetree = new MeanExpTreeFromMultiVariate(process,3,INTEGRAL,false,meanexp);
			omegatree = new RatioTree(nonsynratetree,gcsynratetree);
		}

		cerr << "nucmatrixtree\n";
		nucmatrixtree = new TamuraMatrixTree(tstvtree,gcsynratetree,atsynratetree,rootgc,normalise);

		cerr << "phyloprocess\n";
		if (omegaratiotree)	{
			matrixtree = new GCMatrixTree((CodonStateSpace*) codondata->GetStateSpace(), nucmatrixtree, omegatree, One);
			if (conjpath)	{
				pathconjtree = new BranchMatrixPathConjugateTree(gcsynratetree, matrixtree, codondata);
				phyloprocess = new PathConjugatePhyloProcess(pathconjtree);
			}
			else	{
				pathconjtree = 0;
				phyloprocess = new BranchMatrixPhyloProcess(gcsynratetree, matrixtree, codondata);
			}
		}
		else	{
			cerr << "here\n";
			matrixtree = 0;
			if (conjpath)	{
				pathconjtree = new BranchMatrixPathConjugateTree(gcsynratetree, nucmatrixtree, nucdata);
				phyloprocess = new PathConjugatePhyloProcess(pathconjtree);
			}
			else	{
				pathconjtree = 0;
				phyloprocess = new BranchMatrixPhyloProcess(gcsynratetree, nucmatrixtree, nucdata);
			}
		}

		phyloprocess->Unfold();
		if (sample)	{
			phyloprocess->Sample();
		}

		// register model
		RootRegister(One);
		if (RootAlpha)	{
			RootRegister(RootAlpha);
			RootRegister(RootBeta);
		}
		if (chronoprior == 1)	{
			RootRegister(MeanChi);
			RootRegister(MeanChi2);
		}
		if (autoregressive)	{
			RootRegister(Zero);
		}
		Register();

		MakeScheduler();
		if (sample)	{
			Update();
			if (RootAlpha)	{
				cerr << "starting chrono : " << GetCalibratedChronogram()->GetLogProb() << '\n';
				cerr << "scale progeny : " << GetCalibratedChronogram()->GetScale()->down.size() << '\n';
			}
		}
	}

	~ConjugateTamuraModel() {}

	ConjugateInverseWishart* GetConjugateInverseWishart() {
		ConjugateInverseWishart* tmp = dynamic_cast<ConjugateInverseWishart*>(sigma);
		if (! tmp)	{
			cerr << "error : dynamic castt of conjugate inverse wishart : " << sigma << '\t' << tmp << '\n';
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

	void MakeScheduler()	{

		if (conjpath)	{
			scheduler.Register(new DSemiConjugateMappingMove(phyloprocess,pathconjtree),1,"mapping + sufficient stat");
		}
		else	{
			scheduler.Register(new SimpleMove(phyloprocess,1),1,"mapping");
		}

		for (int i=0; i<nrep; i++)	{
			if (chronoprior == 1)	{
				scheduler.Register(new SimpleMove(Chi,1),10,"chi");
				scheduler.Register(new SimpleMove(Chi,0.1),10,"chi");
				scheduler.Register(new SimpleMove(Chi2,1),10,"chi2");
				scheduler.Register(new SimpleMove(Chi2,0.1),10,"chi2");
			}
			scheduler.Register(new SimpleMove(chronogram,1),10,"chrono");
			scheduler.Register(new SimpleMove(chronogram,0.1),10,"chrono");
			scheduler.Register(new SimpleMove(chronogram,0.01),10,"chrono");

			if (RootAlpha)	{
				scheduler.Register(new SimpleMove(GetCalibratedChronogram()->GetScale(),1),1,"root age");
				scheduler.Register(new SimpleMove(GetCalibratedChronogram()->GetScale(),0.1),1,"root age");
				scheduler.Register(new SimpleMove(GetCalibratedChronogram()->GetScale(),0.01),1,"root age");
			}

			scheduler.Register(new ConjugateMultiVariateMove(GetConjugateInverseWishart(),GetConjugateMultiVariateTreeProcess(),10,10),1,"conjugate sigma - process");
			scheduler.Register(new ConjugateMultiVariateMove(GetConjugateInverseWishart(),GetConjugateMultiVariateTreeProcess(),1,10),1,"conjugate sigma - process");
			scheduler.Register(new ConjugateMultiVariateMove(GetConjugateInverseWishart(),GetConjugateMultiVariateTreeProcess(),0.1,10),1,"conjugate sigma - process");
			scheduler.Register(new ConjugateMultiVariateMove(GetConjugateInverseWishart(),GetConjugateMultiVariateTreeProcess(),0.01,10),1,"conjugate sigma - process");

			if (autoregressive)	{
				scheduler.Register(new SimpleMove(phi,10),100,"phi");
				scheduler.Register(new SimpleMove(phi,1),100,"phi");
				scheduler.Register(new SimpleMove(phi,0.1),100,"phi");
				scheduler.Register(new SimpleMove(phi,0.01),100,"phi");

				scheduler.Register(new SimpleMove(mean,1),100,"mean");
				scheduler.Register(new SimpleMove(mean,0.1),100,"mean");
				scheduler.Register(new SimpleMove(mean,0.01),100,"mean");
			}

			scheduler.Register(new SimpleMove(rootgc,1),10,"root gc");
			scheduler.Register(new SimpleMove(rootgc,0.1),10,"root gc");
			scheduler.Register(new SimpleMove(rootgc,0.01),10,"root gc");
		}
	}

};

