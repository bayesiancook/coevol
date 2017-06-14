
#include "Sample.h"
#include "NearlyNeutralModel1.h"
#include "MeanValTree.h"
#include "MeanICTree.h"
#include "MeanCovMatrix.h"
#include "StringStreamUtils.h"
#include "MeanChronogram.h"
#include "MeanChronoBubbleTree.h"

class BranchOmegaMultivariateSample : public Sample	{

	private:
	string modeltype;
	string datafile;
	string treefile;
	string contdatafile;
	string calibfile;

	double rootage;
	double rootstdev;

	bool withNe;
	bool clamptree;
	bool meanexp;
	GeneticCodeType type;

	int contdatatype;

	int nrep;

	int df;

	public:

	string GetModelType() {return modeltype;}

	BranchOmegaMultivariateModel* GetModel() {return (BranchOmegaMultivariateModel*) model;}

	BranchOmegaMultivariateSample(string filename, int inburnin, int inevery, int inuntil) : Sample(filename,inburnin,inevery,inuntil) {	
		Open();
	}


	void Open()	{

		ifstream is((name + ".param").c_str());
		// check that file exists
		if (!is)	{
			cerr << "error : cannot find file : " << name << ".param\n";
			exit(1);
		}

		double priorsigma = 1;

		is >> modeltype;
		is >> type;
		is >> datafile >> treefile >> contdatafile;
		is >> calibfile >> rootage >> rootstdev;
		is >> contdatatype;
		is >> withNe;
		is >> meanexp;
		is >> nrep;
		is >> df;
		is >> priorsigma;
		is >> clamptree;

		int check;
		is >> check;
		if (check)	{
			cerr << "error when reading model\n";
			exit(1);
		}

		is >> chainevery >> chainuntil >> chainsize;

		// the chain's saving frequency, upper limit and current size
		// not to be confused with the sample's subsampling frequency, upper limit and size

		// make a new model depending on the type obtained from the file
		if (modeltype == "CONJUGATEBRANCHOMEGAMULTIVARIATE")	{
			model = new BranchOmegaMultivariateModel(datafile,treefile,contdatafile,calibfile,rootage,rootstdev,priorsigma,df,contdatatype,withNe,clamptree,meanexp,nrep,false,type);
		}
		else	{
			cerr << "error when opening file "  << name << " : does not recognise model type : " << modeltype << '\n';
			exit(1);
		}

		// read model (i.e. chain's last point) from <name>.param
		model->FromStream(is);
		// model->Update();

		// open <name>.chain, and prepare stream and stream iterator
		OpenChainFile();
		// now, size is defined (it is the total number of points with which this Sample object will make all its various posterior averages)
		// all these points can be accessed to (only once) by repeated calls to GetNextPoint()

		cerr << "number of points to read : " << size << '\n';
		cerr << '\n';
	}

	void ReadIC()	{

		MeanICTree* ictree = new MeanICTree(GetModel()->GetMultiVariateProcess());
		ictree->Reset();
		for (int i=0; i<size; i++)	{
			cerr << '.';
			GetNextPoint();
			ictree->Add();
		}
		cerr << '\n';
		ictree->Normalise();
		ofstream os((name + ".postmeanic").c_str());
		ictree->Tabulate(os,false);
		ofstream los((name + ".postmeanleafic").c_str());
		ictree->Tabulate(los,true);
		cerr << "mean independent contrasts in " << name  << ".postmeanic\n";
		cerr << "for terminal branches only in " << name  << ".postmeanleafic\n";
	}

	void Read(bool printlog, bool printmean, bool printci, bool printstdev, bool withleaf, bool withinternal, string mulreg, bool tex, double xscale, double yscale, double nodescale, double nodepower, double barwidth, int fontsize, bool bubbletext, double meanreg, double stdevreg)	{

		int Ncont = GetModel()->Ncont;
		int dim = GetModel()->GetCovMatrix()->GetDim();

		MeanChronogram* meanchrono = new MeanChronogram(GetModel()->GetTree());
		MeanExpNormTree* meansynrate = new MeanExpNormTree(GetModel()->GetTree(),false,printlog,printmean,printci,printstdev,withleaf,withinternal,meanreg,stdevreg);
		MeanExpNormTree* meanomega = new MeanExpNormTree(GetModel()->GetTree(),false,printlog,printmean,printci,printstdev,withleaf,withinternal);
		MeanExpNormTree* meanmutrate = new MeanExpNormTree(GetModel()->GetTree(),false,printlog,printmean,printci,printstdev,withleaf,withinternal);
		MeanExpNormTree* meanu_Ne = new MeanExpNormTree(GetModel()->GetTree(),false,printlog,printmean,printci,printstdev,withleaf,withinternal);
		MeanExpNormTree* meanNe = new MeanExpNormTree(GetModel()->GetTree(),false,printlog,printmean,printci,printstdev,withleaf,withinternal);	
	
		meanu_Ne->SetLogScale(10.0);
		meanNe->SetLogScale(10.0);

		double alpha[dim];
		MeanExpNormTree** tree = new MeanExpNormTree*[Ncont];
		for (int k=0; k<Ncont; k++)	{
			tree[k] = new MeanExpNormTree(GetModel()->GetTree(),false,printlog,printmean,printci,printstdev,withleaf,withinternal);
		}
		
		double* synrateslope(GetModel()->GetSynrateSlope());
		double* omegaslope(GetModel()->GetOmegaSlope());
		double* u_Neslope(GetModel()->GetU_NeSlope());
		
		
		int indice1(0);
		int indice2(0);
		int indice3(0);
		int indice4(0);
		
		for (int k=0; k<Ncont; k++)	{
			if (GetModel()->GetContinuousData()->GetCharacterName(k) == "maturity") {
				indice1 = k+dim-Ncont;
			}	
			else if (GetModel()->GetContinuousData()->GetCharacterName(k) == "mass") {
				indice2 = k+dim-Ncont;
			}
			else if (GetModel()->GetContinuousData()->GetCharacterName(k) == "longevity") {
				indice3 = k+dim-Ncont;
			}
			else if (GetModel()->GetContinuousData()->GetCharacterName(k) == "piNpiS") {
				indice4 = k+dim-Ncont;
			}	
		}
		
		/*
		MeanCovMatrix*  maty1 = new MeanCovMatrix(dim);
		MeanCovMatrix*  maty2 = new MeanCovMatrix(dim);
		MeanCovMatrix*  maty3 = new MeanCovMatrix(dim);
		MeanCovMatrix*  maty4 = new MeanCovMatrix(dim);
		*/
		MeanCovMatrix*  mat = new MeanCovMatrix(dim);
		

		for (int i=0; i<size; i++)	{
			cerr << '.';
	
			GetNextPoint();

			GetModel()->GetSynRateTree()->specialUpdate();
			GetModel()->GetChronogram()->specialUpdate();
			GetModel()->GetSynrateNodeTree()->specialUpdate();
			GetModel()->GetOmegaNodeTree()->specialUpdate();
			GetModel()->GetNeNodeTree()->specialUpdate();
			
			Var<Real>* gamma(GetModel()->GetGamma());

			meanchrono->Add(GetModel()->GetChronogram());

			if (!withNe) {meanmutrate->Add(GetModel()->GetMultiVariateProcess(), GetModel()->GetChronogram(), 0);}
			if (withNe) {meanNe->Add(GetModel()->GetMultiVariateProcess(), GetModel()->GetChronogram(), 0);}
			meansynrate->Add(GetModel()->GetSynrateNodeTree(), GetModel()->GetChronogram());
			meanomega->Add(GetModel()->GetOmegaNodeTree(), GetModel()->GetChronogram());
			meanu_Ne->Add(GetModel()->GetNeNodeTree(), GetModel()->GetChronogram());
			
				
			for (int k=0; k<Ncont; k++)	{
				tree[k]->Add(GetModel()->GetMultiVariateProcess(), GetModel()->GetChronogram(), GetModel()->GetL()+k);
			}

			CovMatrix& m = *(GetModel()->GetCovMatrix());
			mat->Add(&m);
			
			/*
			double mas1_1[dim][dim];
			double mas1_2[dim][dim];
			double mas1_3[dim][dim];
			double mas1_4[dim][dim];
			double mas1_5[dim][dim];
			CovMatrix my1(dim);

			for (int i = 0; i<dim; i++) {
				for (int j = 0; j<dim; j++) {
					if (i == indice1) {
						mas1_1[i][j]=0;
						for (int k=0; k<dim; k++) {
							mas1_1[i][j] += synrateslope[k] * m[k][j];
						}
					}		
					else {	
					mas1_1[i][j] = m[i][j];
					}
				}
			}
			
			for (int i = 0; i<dim; i++) {
				for (int j = 0; j<dim; j++) {
					if (j == indice1) {
						mas1_2[i][j]=0;
						for (int k=0; k<dim; k++) {
							mas1_2[i][j] += synrateslope[k] * mas1_1[i][k];
						}
					}	
					else {	
					mas1_2[i][j] = mas1_1[i][j];
					}
				}
			}
			
			for (int i = 0; i<dim; i++) {
				for (int j = 0; j<dim; j++) {
					if (i == indice2) {
						mas1_3[i][j]=0;
						for (int k=0; k<dim; k++) {
							mas1_3[i][j] += omegaslope[k] * mas1_2[k][j] * *gamma;
						}
					}		
					else {	
					mas1_3[i][j] = mas1_2[i][j];
					}
				}
			}
			
			for (int i = 0; i<dim; i++) {
				for (int j = 0; j<dim; j++) {
					if (j == indice2) {
						mas1_4[i][j]=0;
						for (int k=0; k<dim; k++) {
							mas1_4[i][j] += omegaslope[k] * mas1_3[i][k] * *gamma;
						}
					}	
					else {	
					mas1_4[i][j] = mas1_3[i][j];
					}
				}
			}
			
			for (int i = 0; i<dim; i++) {
				for (int j = 0; j<dim; j++) {
					if (i == indice3) {
						mas1_5[i][j]=0;
						for (int k=0; k<dim; k++) {
							mas1_5[i][j] += u_Neslope[k] * mas1_4[k][j];
						}
					}		
					else {	
					mas1_5[i][j] = mas1_4[i][j];
					}
				}
			}
			
			for (int i = 0; i<dim; i++) {
				for (int j = 0; j<dim; j++) {
					my1[i][j]=0;
					if (j == indice3) {
						for (int k=0; k<dim; k++) {
							my1[i][j] += u_Neslope[k] * mas1_5[i][k];
						}
					}	
					else {	
					my1[i][j] = mas1_5[i][j];
					}
				}
			}
			
			
			
			
			double mas2[dim][dim];
			CovMatrix my2(dim);
			
			
			for (int i = 0; i<dim; i++) {
				for (int j = 0; j<dim; j++) {
					if (i == indice4) {
						mas2[i][j] = 0;
						for (int k=0; k<dim; k++) {
							mas2[i][j] += synrateslope[k] * m[k][j];
						}
					}		
					else {	
					mas2[i][j] = m[i][j];
					}
				}
			}
			
			for (int i = 0; i<dim; i++) {
				for (int j = 0; j<dim; j++) {
					if (j == indice4) {
						my2[i][j] = 0;
						for (int k=0; k<dim; k++) {
							my2[i][j] += synrateslope[k] * mas2[i][k];
						}
					}	
					else {	
					my2[i][j] = mas2[i][j];
					}
				}
			}
			
			
			
			double mas3[dim][dim];
			CovMatrix my3(dim);
			
			
			for (int i = 0; i<dim; i++) {
				for (int j = 0; j<dim; j++) {
					if (i == indice4) {
						mas3[i][j] = 0;
						for (int k=0; k<dim; k++) {
							mas3[i][j] += omegaslope[k] * m[k][j] * *gamma;
						}
					}		
					else {	
					mas3[i][j] = m[i][j];
					}
				}
			}
			
			for (int i = 0; i<dim; i++) {
				for (int j = 0; j<dim; j++) {
					if (j == indice4) {
						my3[i][j] = 0;
						for (int k=0; k<dim; k++) {
							my3[i][j] += omegaslope[k] * mas3[i][k] * *gamma;
						}
					}	
					else {	
					my3[i][j] = mas3[i][j];
					}
				}
			}
			
			
			
			double mas4[dim][dim];
			CovMatrix my4(dim);
			
			
			for (int i = 0; i<dim; i++) {
				for (int j = 0; j<dim; j++) {
					if (i == indice4) {
						mas4[i][j] = 0;
						for (int k=0; k<dim; k++) {
							mas4[i][j] += u_Neslope[k] * m[k][j];
						}
					}		
					else {	
					mas4[i][j] = m[i][j];
					}
				}
			}
			
			for (int i = 0; i<dim; i++) {
				for (int j = 0; j<dim; j++) {
					if (j == indice4) {
						my4[i][j] = 0;
						for (int k=0; k<dim; k++) {
							my4[i][j] += u_Neslope[k] * mas4[i][k];
						}
					}	
					else {	
					my4[i][j] = mas4[i][j];
					}
				}
			}
			
			maty1->Add(&my1);
			maty2->Add(&my2);
			maty3->Add(&my3);
			maty4->Add(&my4);
			*/
			
		}
		cerr << '\n';
		cerr << "normalise\n";

		mat->Normalize();
		ofstream cout((GetName() + ".cov").c_str());

		cout << *mat;

		cerr << "covariance matrix in " << name << ".cov\n";
		cerr << '\n';
		
		mat->PrintSlopes(cout);
		mat->PrintSlopes2(cout);
		
		/*
		maty1->Normalize();
		ofstream cout1((GetName() + ".cov_all").c_str());

		cout1 << *maty1;

		cerr << "covariance matrix in " << name << ".cov_all\n";
		cerr << '\n';
		
		maty1->PrintSlopes(cout1);
		maty1->PrintSlopes2(cout1);

		
		maty2->Normalize();
		ofstream cout2((GetName() + ".cov_dS->piNpiS").c_str());

		cout2 << *maty2;

		cerr << "covariance matrix in " << name << ".cov_dS->piNpiS\n";
		cerr << '\n';
		
		maty2->PrintSlopes(cout2);
		maty2->PrintSlopes2(cout2);
		
		
		maty3->Normalize();
		ofstream cout3((GetName() + ".cov_omega->piNpiS").c_str());

		cout3 << *maty3;

		cerr << "covariance matrix in " << name << ".cov_omega->piNpiS\n";
		cerr << '\n';
		
		maty3->PrintSlopes(cout3);
		maty3->PrintSlopes2(cout3);
		
		
		if (!withNe) {
			maty4->Normalize();
			ofstream cout4((GetName() + ".cov_Ne->piNpiS").c_str());
	
			cout4 << *maty4;

			cerr << "covariance matrix in " << name << ".cov_Ne->piNpiS\n";
			cerr << '\n';
			
			maty4->PrintSlopes(cout4);
			maty4->PrintSlopes2(cout4);
		}

		if (withNe) {
			maty4->Normalize();
			ofstream cout4((GetName() + ".cov_u->piNpiS").c_str());
	
			cout4 << *maty4;

			cerr << "covariance matrix in " << name << ".cov_u->piNpiS\n";
			cerr << '\n';
			
			maty4->PrintSlopes(cout4);
			maty4->PrintSlopes2(cout4);
		}

		*/

		meanchrono->Normalise();
		ofstream chos((GetName() + ".postmeandates.tre").c_str());
		meanchrono->ToStream(chos);

		ofstream cchos((GetName() + ".postmeandates.tab").c_str());
		meanchrono->Tabulate(cchos);


		if (!withNe) {
			meanmutrate->Normalise();
			ofstream mos((GetName() + ".postmeanmutrate.tre").c_str());
			meanmutrate->ToStream(mos);
			cerr << "reconstructed variations of mutrate in " << name << ".postmeanmutrate.tre\n";
			cerr << "pp of mean leaf values > root value : " << meanmutrate->GetPPLeafRoot() << '\n';
		}	


		meanomega->Normalise();
		ofstream oos((GetName() + ".postmeanomega.tre").c_str());
		meanomega->ToStream(oos);
		cerr << "reconstructed variations of omega in " << name << ".postmeanomega.tre\n";
		cerr << "pp of mean leaf values > root value : " << meanomega->GetPPLeafRoot() << '\n';


		if (withNe) {
			meanNe->Normalise();
			ofstream Neos((GetName() + ".postmeanNe.tre").c_str());
			meanNe->ToStream(Neos);
			cerr << "reconstructed variations of Ne in " << name << ".postmeanNe.tre\n";
			cerr << "pp of mean leaf values > root value : " << meanNe->GetPPLeafRoot() << '\n';
		}
		
		meansynrate->Normalise();
		ofstream sos((GetName() + ".postmeansynrate.tre").c_str());
		meansynrate->ToStream(sos);
		cerr << "reconstructed variations of Ks in " << name << ".postmeansynrate.tre\n";
		cerr << "pp of mean leaf values > root value : " << meansynrate->GetPPLeafRoot() << '\n';



		meanu_Ne->Normalise();
		if (!withNe) {
			ofstream mnos((GetName() + ".postmeanNe.tre").c_str());
			meanu_Ne->ToStream(mnos);
			cerr << "reconstructed variations of Ne in " << name << ".postmeanNe.tre\n";
			cerr << "pp of mean leaf values > root value : " << meanu_Ne->GetPPLeafRoot() << '\n';
		}
		if (withNe) {
			ofstream mnos((GetName() + ".postmeanmutrate.tre").c_str());
			meanu_Ne->ToStream(mnos);
			cerr << "reconstructed variations of mutrate in " << name << ".postmeanmutrate.tre\n";
			cerr << "pp of mean leaf values > root value : " << meanu_Ne->GetPPLeafRoot() << '\n';
		}
		

		for (int k=0; k<Ncont; k++)	{
			tree[k]->Normalise();
			ostringstream s;
			s << GetName() << ".postmean" << k+1 << ".tre";
			ofstream os(s.str().c_str());
			tree[k]->ToStream(os);
			cerr << "reconstructed variations of continuous character # " << k+1 << " in "  << name << ".postmean" << k+1 << ".tre\n";
			cerr << "pp of mean leaf values > root value : " << tree[k]->GetPPLeafRoot() << '\n';
		}

		if (!withNe) {
		ofstream mmos((GetName() + ".postmeanmutrate.tab").c_str());
		meanmutrate->Tabulate(mmos);
		mmos.close();
		}
		
		ofstream ssos((GetName() + ".postmeansynrate.tab").c_str());
		meansynrate->Tabulate(ssos);
		ssos.close();

		ofstream ooos((GetName() + ".postmeanomega.tab").c_str());
		meanomega->Tabulate(ooos);
		ooos.close();

		if (withNe) {
			ofstream Neoos((GetName() + ".postmeanNe.tab").c_str());
			meanNe->Tabulate(Neoos);
			Neoos.close();
		}	
		
		if (!withNe) {
			ofstream mnoos((GetName() + ".postmeanNe.tab").c_str());
			meanu_Ne->Tabulate(mnoos);
			mnoos.close();
		}
		if (withNe) {
			ofstream mnoos((GetName() + ".postmeanmutrate.tab").c_str());
			meanu_Ne->Tabulate(mnoos);
			mnoos.close();
		}	
		
		

		for (int k=0; k<Ncont; k++)	{
			ostringstream s;
			s << GetName() << ".postmean" << k+1 << ".tab";
			ofstream os(s.str().c_str());
			tree[k]->Tabulate(os);
		}

		cerr << '\n';
	}
};


int main(int argc, char* argv[])	{

	int burnin = 0;
	int every = 1;
	int until = -1;
	string name; 

	bool ic = false;

	bool printlog = false;
	bool printmean = false;
	bool printci = true;
	bool printstdev = false;
	bool withleaf = true;
	bool withinternal = true;

	string mulreg = "";

	bool tex = false;
	double nodescale = 5;
	double nodepower = 1;
	double barwidth = 0.04;
	double x = 6;
	double y = 8;
	// bool withnodebubbles = true;
	// bool withtimeintervals = true;
	int fontsize = 4;

	bool bubbletext = false;

	double meanreg = 0;
	double stdevreg = 0;

	try	{

		if (argc == 1)	{
			throw(0);
		}

		int i = 1;
		while (i < argc)	{
			string s = argv[i];

			if (s == "-tex")	{
				tex = true;
			}
			else if (s == "-ns")	{
				i++;
				nodescale = atof(argv[i]);
			}
			else if (s == "-np")	{
				i++;
				nodepower= atof(argv[i]);
			}
			else if (s == "-bw")	{
				i++;
				barwidth = atof(argv[i]);
			}
			else if (s == "-fs")	{
				i++;
				fontsize = atoi(argv[i]);
			}
			else if (s == "-xs")	{
				i++;
				x = atof(argv[i]);
			}
			else if (s == "-ys")	{
				i++;
				y = atof(argv[i]);
			}
			else if (s == "+bt")	{
				bubbletext = true;
			}
			else if (s == "-bt")	{
				bubbletext = false;
			}
			else if (s == "-pointwise")	{
				i++;
				meanreg = atof(argv[i]);
				i++;
				stdevreg = atof(argv[i]);
			}
			/*
			else if (s == "+nb")	{
				withnodebubbles = true;
			}
			else if (s == "-nb")	{
				withnodebubbles = false;
			}
			*/
			else if (s == "-ic")	{
				ic = true;
			}
			else if (s == "+log")	{
				printlog = true;
			}
			else if (s == "-log")	{
				printlog = false;
			}
			else if (s == "+mean")	{
				printmean = true;
			}
			else if (s == "-mean")	{
				printmean = false;
			}
			else if (s == "+stdev")	{
				printstdev = true;
			}
			else if (s == "-stdev")	{
				printstdev = false;
			}
			else if (s == "+ci")	{
				printci = true;
			}
			else if (s == "-ci")	{
				printci = false;
			}
			else if (s == "+leaf")	{
				withleaf = true;
			}
			else if (s == "-leaf")	{
				withleaf = false;
			}
			else if (s == "+internal")	{
				withinternal = true;
			}
			else if (s == "-internal")	{
				withinternal = false;
			}
			else if (s == "-mulreg")	{
				i++;
				mulreg = argv[i];
			}
			else if (s == "-partial")	{
				i++;
				mulreg = argv[i];
			}
			else if ( (s == "-x") || (s == "-extract") )	{
				i++;
				if (i == argc) throw(0);
				s = argv[i];
				if (! IsInt(s))	{
					throw(0);
				}
				burnin = atoi(argv[i]);
				i++;
				if (i == argc) throw(0);
				s = argv[i];
				if (IsInt(s))	{
					every = atoi(argv[i]);
					i++;
					if (i == argc) throw(0);
					s = argv[i];
					if (IsInt(s))	{
						until = atoi(argv[i]);
					}
					else	{
						i--;
					}
				}
				else	{
					i--;
				}
			}
			else	{
					if (i != (argc -1))	{
					throw(0);
				}
				name = argv[i];
			}
			i++;
		}
		if (name == "")	{
			throw(0);
		}
	}
	catch(...)	{
		cerr << "readcoevol [-x <burnin> <every> <until>] <chainname> \n";
		cerr << '\n';
		exit(1);
	}

	BranchOmegaMultivariateSample sample(name,burnin,every,until);

	if (ic)	{
		sample.ReadIC();
		exit(1);
	}

	sample.Read(printlog,printmean,printci,printstdev,withleaf,withinternal,mulreg,tex,x,y,nodescale,nodepower,barwidth,fontsize,bubbletext,meanreg,stdevreg);

}



