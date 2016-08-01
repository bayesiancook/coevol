
#ifndef AMINOACIDMODEL_H
#define AMINOACIDMODEL_H

#include "MeanValTree.h"
#include "BaseType.h"
#include "RandomTypes.h"
#include "ProbModel.h"
#include "CalibratedChronogram.h"
#include "BranchProcess.h"
#include "OneMatrixPhyloProcess.h"
#include "ContinuousData.h"
#include "MeanExpTree.h"
#include "Normal.h"
#include "GeneralConjugatePath.h"
#include "Jeffreys.h"
#include "ConjugateMultiVariateTreeProcess.h"
#include "MultiVariatePropagateMove.h"
#include "MultiVariateRelRateCompensatoryMove.h"
#include "MeanChronogram.h"
#include "AuxCoevol.h"
#include "CodonSequenceAlignment.h"
#include "ProteinSequenceAlignment.h"
#include "AminoAcidOmegaSubMatrix.h"
#include "SimilarityMatrix.h"
#include "AminoAcidMatrixTree.h"
#include "BranchMatrixPhyloProcess.h"

class RenormalizedPosRealVector : public Dvar<Profile>	{

	public:

	RenormalizedPosRealVector(Var<PosRealVector>* inup)	{
		int dimension = inup->GetDim();
		setval(Profile(dimension));
		bkvalue = Profile(dimension);
		up = inup;
		Register(up);
		specialUpdate();
	}

	double GetEntropy()	{
		double tot = 0;
		for (int i=0; i<GetDim(); i++)	{
			tot -= (*this)[i] * log((*this)[i]);
		}
		return tot;
	}

	protected:

	void specialUpdate()	{
		double tot = 0;
		for (int i=0; i<GetDim(); i++)	{
			tot += (*up)[i];
		}
		for (int i=0; i<GetDim(); i++)	{
			(*this)[i] = (*up)[i] / tot;
		}
	}

	Var<PosRealVector>* up;

};


class AminoAcidOmegaModel : public ProbModel {

	public:

	// data fields

	// ---------
	// the fixed parameters of the model
	// ---------

	// a fixed tree (read from file)
	Tree* tree;
	SequenceAlignment* AAdata;
	CodonSequenceAlignment* codondata;
	ContinuousData* contdata;
	ContinuousData* nuccomp;
	TaxonSet* taxonset;

	// number of columns
	int Nsite;
	// number of states (4 for nucleic acids, 20 for amino-acids. 61 for codons)
	int Nstate;

	int Ncont;

	// ---------
	// the random variables of the model
	// ---------

	Const<Real>* Zero;
	Const<PosReal>* One;

	Const<PosReal>* RootAlpha;
	Const<PosReal>* RootBeta;
	Dvar<PosReal>* PriorMu;
	Gamma* mu;
	Chronogram* chronogram;

	LengthTree* lengthtree;
	LengthTree* synratetree;
	GammaTree* syngammatree;
	MeanExpTreeFromMultiVariate* meanexptree;

	MeanExpTreeFromMultiVariate* omegatree;

	JeffreysIIDArray* DiagArray;
	SigmaZero* sigmaZero;
	ConjugateInverseWishart* sigma;


	ConjugateMultiVariateTreeProcess* process;


	// nucleotide mutation matrix is relrate * stationary
	IIDExp* exprelrate;
	RenormalizedPosRealVector* relrate;
	// Dirichlet* relrate;
	Dirichlet* stationary;
	AminoAcidModelType type;

	//Similarity Matrix
	SimilarityMatrix * aasimilarityMatrix ;
	RandomSubMatrix* AAmatrix;
	SubMatrix* aasubmatrix;

	BranchValPtrTree<RandomSubMatrix>* matrixtree;

	// phylo process
	PathConjugateTree* pathconjtree;
	PhyloProcess* phyloprocess;

	// if true: covariances are all set equal to 0
	bool clampdiag;

	bool clamptree;
	bool meanexp;

	int omegaindex;

	// total number of substitution parameters modelled as non homogeneous
	int L;

	int nrep;

	int df;

	Const<RealVector>* rootmean;
	Const<PosRealVector>* rootvar;

	bool Unconstrained()	{
		return (syngammatree != 0);
	}


	AminoAcidOmegaModel(string datafile, string treefile, string contdatafile, string calibfile, double rootage, double rootstdev, double priorsigma, int indf, int contdatatype, bool inclamptree, bool inmeanexp, int innrep, bool sample, AminoAcidModelType inmodeltype, GeneticCodeType geneticcodetype, string rootfile= "NONE")	{

		clamptree = inclamptree;
		meanexp = inmeanexp;
		L = 2;

		cerr << "aadata\n";
		AAdata = new FileSequenceAlignment(datafile);

		Nsite = GetData()->GetNsite();
		Nstate = GetData()->GetNstate();
		cerr << "data read : " << Nsite << '\t' << Nstate << '\n';

		if(Nstate == Nnuc){
			cerr << "Converting the data into amino acid alignement "<< endl;
			codondata = new CodonSequenceAlignment(AAdata,true, geneticcodetype);
			AAdata = new ProteinSequenceAlignment(codondata);
		}

		Nsite = GetData()->GetNsite();
		Nstate = GetData()->GetNstate();

		cerr << "data read : " << Nsite << '\t' << Nstate << '\n';

		// verifier Nsate
		if(Nstate != Naa){
			cerr << "Data File Problems: Not a Protein Alignement"<< endl;
			exit(1);
		}

		StateSpace * dataStateSpace = GetData()->GetStateSpace();
		if (dynamic_cast<ProteinStateSpace*> (dataStateSpace) == 0 ){
			cerr <<"Data File Problem: Not a Protein Alignement "<< endl;
			exit(1);
		}

		nrep = innrep;
		  if (nrep == 0){
			nrep = 30;
		}

		taxonset = AAdata->GetTaxonSet();

		cerr << "tree\n";
		cerr << treefile << '\n';
		// get tree from file (newick format)
		tree = new Tree(treefile);
		cerr << "tree ok\n";
		// check whether tree and data fit together
		tree->RegisterWith(taxonset);
		cerr << "register ok\n";

		cerr << "cont : " << contdatafile << '\n';
		// get continuous data from file
		if (contdatafile != "None")	{
			contdata = new FileContinuousData(contdatafile);
			Ncont = contdata->GetNsite();
			cerr << "Continious caracters number: "<<Ncont << endl;
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

		if(calibfile == "Unconstrained"){
			cerr << "unconstrained"<<endl;
			L--;
			PriorMu = new Const<PosReal>(0.1);
			mu = new Gamma(One,PriorMu);
			chronogram = 0;
			meanexptree = 0;
			syngammatree = new GammaTree(tree,One,mu);
		}
		else {
			PriorMu = new Const<PosReal>(1);
			mu = new Gamma(One,PriorMu);
			mu->ClampAt(1);
			syngammatree = 0;

			if (calibfile != "None"){
				double a = rootage * rootage / rootstdev / rootstdev;
				double b = rootage / rootstdev / rootstdev;
				RootAlpha = new Const<PosReal>(a);
				RootBeta = new Const<PosReal>(b);
				CalibrationSet* calibset = new FileCalibrationSet(calibfile, tree);

				chronogram = new CalibratedChronogram(tree,mu,RootAlpha,RootBeta,calibset);
			} 
			else	{
				cerr << "\nnew chrono\n";
				RootAlpha = 0;
				RootBeta = 0;
				chronogram = new Chronogram(tree,One);
				cerr << "ok"<<endl;
			}

			if (clamptree)	{
				chronogram->Clamp();
			}
		}

		if(Unconstrained())	{
			lengthtree = syngammatree;
		}
		else	{
			lengthtree = chronogram;
		}

		rootmean = 0;
		rootvar = 0;
		if (rootfile != "None") {
			cerr << "Prior on root:" << endl;
			rootmean = new Const<RealVector>(RealVector(Ncont+L));
			rootvar = new Const<PosRealVector>(PosRealVector(Ncont+L));
			ifstream is(rootfile.c_str());
			for (int i=0; i<Ncont+L; i++)   {
				double mean, var;
				is >> mean >> var;
				cerr << mean << "\t"<< var << endl;
				(*rootmean)[i] = mean;
				(*rootvar)[i] = var;
			}
		}


		double mindiag = 0.001;
		double maxdiag = 1000;
		DiagArray = new JeffreysIIDArray(Ncont+L,mindiag,maxdiag,Zero);
		if (priorsigma == -1)	{
			DiagArray->setval(1.0);
		}
		else	{
			DiagArray->ClampAt(priorsigma);
		}
		sigmaZero = new SigmaZero(DiagArray);

		df = Ncont + L + indf;
		cerr << "sigma\n";
		sigma = new ConjugateInverseWishart(sigmaZero, df);

		cerr << "process\n";

		if(rootmean)	{
			process = new ConjugateMultiVariateTreeProcess(sigma,lengthtree, 0 , 0 , rootmean, rootvar);
		}
		else	{
			process = new ConjugateMultiVariateTreeProcess(sigma,lengthtree);
		}

		for (int i=0; i<Ncont; i++)	{
			process->SetAndClamp(contdata,L+i,i,contdatatype);
		}

		for (int l=0; l<L; l++)	{
			process->CutOff(1,l);
		}

		cerr << "syn and omega\n";

		if(Unconstrained()){
			process->GetMultiNormal(process->GetRoot())->ClampAt(0 , 0);
			synratetree = syngammatree;
			omegatree = new MeanExpTreeFromMultiVariate(process,0,MEAN,false,meanexp);
			omegaindex = 0;
		}
		else{
			process->GetMultiNormal(process->GetRoot())->ClampAt(0 , 1);
			synratetree = new MeanExpTreeFromMultiVariate(process,0,INTEGRAL,false,meanexp);
			omegatree = new MeanExpTreeFromMultiVariate(process,1,MEAN,false,meanexp);
			omegaindex = 1;
		}

		type = inmodeltype;
		cerr << "Using the model based on " << type << endl;
		switch(type){
			case 0 : aasimilarityMatrix = new PolarityBasedMatrix();break;
			case 1 : aasimilarityMatrix = new VolumeBasedMatrix(); break;
			case 2 : aasimilarityMatrix = new ChargeBasedMatrix(); break;
			case 3 : aasimilarityMatrix = new PolarityAndVolumeBasedMatrix(); break;
			case 4 : cerr<< " NO Model Specfified" << endl; exit(0);
		}

		aasimilarityMatrix->Affiche();
		cerr << "matrix\n";
		exprelrate = new IIDExp(Naa*(Naa-1)/2);
		relrate = new RenormalizedPosRealVector(exprelrate);
		// relrate = new Dirichlet(Naa*(Naa-1)/2);
		stationary = new Dirichlet(Naa);
		matrixtree = new AminoAcidOmegaMatrixTree(aasimilarityMatrix,relrate, stationary, omegatree, One);

		// substitution process
		pathconjtree = new BranchMatrixPathConjugateTree (synratetree, matrixtree, AAdata);
		phyloprocess = new PathConjugatePhyloProcess(pathconjtree);
		//phyloprocess = new BranchMatrixPhyloProcess(synratetree,matrixtree,AAdata);

		cerr << "unfold\n";
		phyloprocess->Unfold();

		cerr << "sample\n";
		if (sample)	{
			if (phyloprocess)	{
				phyloprocess->Sample();
			}
		}

		cerr << "register\n";
		// register model
		RootRegister(PriorMu);
		RootRegister(Zero);
		RootRegister(One);
		if (RootAlpha)	{
			RootRegister(RootAlpha);
			RootRegister(RootBeta);
		}
		RootRegister(exprelrate);
		RootRegister(stationary);
		if(rootmean){
			RootRegister(rootmean);
			RootRegister(rootvar);
		}
		Register();

		MakeScheduler();
		if (sample)	{
			Update();
		}
		if (RootAlpha)	{
			cerr << "starting chrono : " << GetCalibratedChronogram()->GetLogProb() << '\n';
			cerr << "scale progeny : " << GetCalibratedChronogram()->GetScale()->down.size() << '\n';
			cerr << "starting chrono : " << GetCalibratedChronogram()->GetScale()->val() << '\n';
		}
	}

	// destructor
	// deallocations should normally be done here
	// but in general, the model is deleted just before the program exits, so we can dispense with it for the moment
	~AminoAcidOmegaModel() {}


	// accessors
	Tree* GetTree() {return tree;}

	SequenceAlignment* GetData()	{
		return AAdata;
	}

	MeanExpTreeFromMultiVariate* GetSynRateTree() {
		if (Unconstrained())	{
			cerr << "error : unconstrained model\n";
			exit(1);
		}
		MeanExpTreeFromMultiVariate* tmp = dynamic_cast<MeanExpTreeFromMultiVariate*>(synratetree);
		if (! tmp)	{
			cerr << "error in get synratetree: null pointer\n";
			exit(1);
		}
		return tmp;
	}

	MeanExpTreeFromMultiVariate* GetOmegaTree() {return omegatree;}

	MultiVariateTreeProcess* GetMultiVariateProcess() {return process;}
	Chronogram* GetChronogram() {return chronogram;}
	LengthTree * GetLengthTree () {return lengthtree;}


	ContinuousData* GetContinuousData() {return contdata;}

	int GetL() {return L;}

	CovMatrix* GetCovMatrix() {return sigma;}

	// probability computation

	double GetLogProb()	{
		return GetLogPrior() + GetLogLikelihood();
	}

	double GetLogPrior()	{
		double total = 0;

		if (Unconstrained())	{
			total += mu->GetLogProb();
			total += syngammatree->GetLogProb();
		}
		else 	{
			total += chronogram->GetLogProb();
		}

		total += DiagArray->GetLogProb();
		total += sigma->GetLogProb();
		total += process->GetLogProb();

		// total += exprelrate->GetLogProb();
		total += stationary->GetLogProb();
		return total;
	}

	double GetLogLikelihood()	{
		double ret = pathconjtree->GetLogProb();
		// double ret = phyloprocess->GetLogProb();
		return ret;
	}

	// summary statistics
	double GetTotalLength()	{
		return GetSynRateTree()->GetTotal();
	}

	double GetMeanOmega()	{
		return omegatree->GetMean();
	}

	double GetMeanSynRate()	{
		if (Unconstrained())	{
			return syngammatree->GetMean();
		}
		return GetSynRateTree()->GetTotal();
	}
	RandomSubMatrix* GetSubMatrix(){
		return matrixtree->GetBranchVal(matrixtree->GetRoot()->Next()->GetBranch());
	}

	/*
	double Move(double tuning = 1)	{
		// Cycle(1,1,verbose,check)
		scheduler.Cycle(1,1,true,true);
		return 1;
	}
	*/

	// MCMC schedule
	virtual void MakeScheduler()	{

		// scheduler.Register(new SimpleMove(phyloprocess,1),1,"mapping");
		scheduler.Register(new DSemiConjugateMappingMove(phyloprocess,pathconjtree),1,"mapping + sufficient stat");

		for (int i=0; i<nrep; i++)	{
			if (! clamptree)	{
				if (Unconstrained())	{
					scheduler.Register(new SimpleMove(mu,1),10,"mu");
					scheduler.Register(new SimpleMove(mu,0.1),10,"mu");
					scheduler.Register(new SimpleMove(syngammatree,1),10,"syngamtree");
					scheduler.Register(new SimpleMove(syngammatree,0.1),10,"syngamtree");
					scheduler.Register(new SimpleMove(syngammatree,0.01),10,"syngamtree");
				}
				else {
					scheduler.Register(new SimpleMove(chronogram,1),10,"chrono");
					scheduler.Register(new SimpleMove(chronogram,0.1),10,"chrono");
					scheduler.Register(new SimpleMove(chronogram,0.01),10,"chrono");
					if (RootAlpha)	{
						scheduler.Register(new SimpleMove(GetCalibratedChronogram()->GetScale(),1),10,"root age");
						scheduler.Register(new SimpleMove(GetCalibratedChronogram()->GetScale(),0.1),10,"root age");
						scheduler.Register(new SimpleMove(GetCalibratedChronogram()->GetScale(),0.01),10,"root age");
					}
				}
			}
	 
			scheduler.Register(new SimpleMove(process,10),30,"multinormal");
			scheduler.Register(new SimpleMove(process,1),30,"multinormal");
			scheduler.Register(new SimpleMove(process,0.1),30,"multinormal");
			scheduler.Register(new SimpleMove(process,0.01),30,"multinormal");

			int n = taxonset->GetNtaxa() * 10;
			scheduler.Register(new MultiVariatePropagateMove(process,1,0.1,0.1),n,"propmove");
			scheduler.Register(new MultiVariatePropagateMove(process,0.1,0.5,0.5),n,"propmove");
			scheduler.Register(new MultiVariatePropagateMove(process,0.1,0.9,0.9),n,"propmove");
			scheduler.Register(new MultiVariatePropagateMove(process,0.01,0.99,0.99),n,"propmove");

			scheduler.Register(new MultiVariateRelRateCompensatoryMove(process,exprelrate,aasimilarityMatrix,1,omegaindex),10,"process relrate comp move");
			scheduler.Register(new MultiVariateRelRateCompensatoryMove(process,exprelrate,aasimilarityMatrix,0.1,omegaindex),10,"process relrate comp move");
			scheduler.Register(new MultiVariateRelRateCompensatoryMove(process,exprelrate,aasimilarityMatrix,0.01,omegaindex),10,"process relrate comp move");

			scheduler.Register(new ConjugateMultiVariateMove(sigma,process,10,10),1,"conjugate sigma - process");
			scheduler.Register(new ConjugateMultiVariateMove(sigma,process,1,10),1,"conjugate sigma - process");
			scheduler.Register(new ConjugateMultiVariateMove(sigma,process,0.1,10),1,"conjugate sigma - process");
			scheduler.Register(new ConjugateMultiVariateMove(sigma,process,0.01,10),1,"conjugate sigma - process");
			scheduler.Register(new ConjugateMultiVariateMove(sigma,process,0.001,10),1,"conjugate sigma - process");

			scheduler.Register(new SimpleMove(DiagArray,10),10,"theta");
			scheduler.Register(new SimpleMove(DiagArray,1),10,"theta");
			scheduler.Register(new SimpleMove(DiagArray,0.1),10,"theta");

			scheduler.Register(new PosRealVectorMove(exprelrate,1,1),10,"relrates");
			scheduler.Register(new PosRealVectorMove(exprelrate,0.1,2),10,"relrates");
			scheduler.Register(new PosRealVectorMove(exprelrate,0.03,4),10,"relrates");
			scheduler.Register(new SimpleMove(exprelrate,0.01),10,"relrates");
			scheduler.Register(new PosRealVectorTranslationMove(exprelrate,1),3,"relrates");
			scheduler.Register(new PosRealVectorTranslationMove(exprelrate,0.1),3,"relrates");
			scheduler.Register(new PosRealVectorTranslationMove(exprelrate,0.01),3,"relrates");

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

		if(Unconstrained()){
			mu->setval(10);
			syngammatree->Sample();
		}
		else	{
			chronogram->Sample();
		}

		DiagArray->Sample();
		sigma->Sample();
		// sigma->SetIdentity();
		process->Sample();
		exprelrate->Sample();
		stationary->Sample();
	}


	void PrintEntries(ostream& os, int* array = 0)	{

		int cumul = 0;
		if(!Unconstrained()){
			os << "Kc\n";
			cumul++;
		}

		os << "omega (Kr/Kc)\n";
		cumul++;

		for (int k=0; k<Ncont; k++)	{
			if ((! array) || (array[cumul] == 1))	{
				os << "character " << k + 1 << '\n';
			}
			cumul++;
		}
	}

	CalibratedChronogram* GetCalibratedChronogram()	{
		if (Unconstrained())	{
			cerr << "error : calibrated chronogram does not exist under unconstrained model\n";
			exit(1);
		}
		return dynamic_cast<CalibratedChronogram*>(chronogram);
	}

	Var<PosReal>* GetScale()	{
		if (RootAlpha)	{
			return GetCalibratedChronogram()->GetScale();
		}
		return 0;
	}

	double GetRootAge()	{
		if (RootAlpha)	{
			return GetCalibratedChronogram()->GetScale()->val();
		}
		return 1;
	}


	// trace
	void TraceHeader(ostream& os)	{
		os << "#logprior\tlnL";
		if(!Unconstrained())	{
			os<<"\tKc";
		}
		os << "\tomega";

		for (int k=0; k<Ncont+L; k++)	{
			for (int l=k+1; l<Ncont+L; l++)	{
				os << '\t' << "sigma_" << k << '_' << l;
			}
		}
		for (int k=0; k<Ncont+L; k++)	{
			os << '\t' << "sigma_" << k << '_' << k;
		}

		if (RootAlpha)	{
			os << '\t' << "rootage";
		}
		if(rootmean)	{
			for(int i=0; i<Ncont+L; i++){
				os << '\t' << "root_" << i;
			}
		}
		os << "\tstatent";
		os << "\trrent";
		if(Unconstrained()){
			os << "\tsyngamma";
			os << "\tsyngammaVar";
		}
		os << '\n';
	}

	// writes all summary statistics on one single line
	// in the same order as that provided by the header
	void Trace(ostream& os)	{
		os << fixed << setprecision(6) << GetLogPrior() << '\t' << GetLogLikelihood();
		if(! Unconstrained())	{
			os << '\t' << GetMeanSynRate();
		}
		os << '\t' << GetMeanOmega();

		for (int k=0; k<Ncont+L; k++)	{
			for (int l=k+1; l<Ncont+L; l++)	{
				os << '\t' << (*sigma)[k][l];
			}
		}
		for (int k=0; k<Ncont+L; k++)	{
			os << '\t' << (*sigma)[k][k];
		}

		if (RootAlpha)	{
			os << '\t' << GetRootAge();
		}
		if(rootmean){
			MultiNormal* m = process->GetMultiNormal(process->GetRoot());  
			for(int i=0; i<Ncont+L; i++)	{
				os << '\t' << exp((*m)[i]);
			}
		}

		os << '\t' << stationary->val().GetEntropy();
		os << '\t' << relrate->GetEntropy();
		os << '\t' << exprelrate->val().GetMean();
		os << '\t' << exprelrate->val().GetVar();
		if(Unconstrained()){
			os << '\t' << syngammatree->GetMean();
			os << '\t' << syngammatree->GetVar();
		}

		os << '\n';
		os.flush();
	}

	// save current state
	void ToStream(ostream& os){
		os << *mu << '\n';
		if (Unconstrained())	{
			os << *syngammatree << '\n';
		}
		else	{
			os << *chronogram << '\n';
			if (RootAlpha)	{
				os << *GetCalibratedChronogram()->GetScale() << '\n';
			}
		}

		os << *DiagArray << '\n';
		os << *sigma << '\n';
		os << *process << '\n';
		os << *exprelrate << '\n';
		os << *stationary << '\n';
	}

	void FromStream(istream& is){
		is >> *mu;
		if (Unconstrained())	{
			is >> *syngammatree;
		}
		else	{
			is >> *chronogram;
			if (RootAlpha)	{
				is >> *GetCalibratedChronogram()->GetScale();
			}
		}

		is >> *DiagArray;
		is >> *sigma;
		is >> *process;
		is >> *exprelrate;
		is >> *stationary;
	}
};

#endif
