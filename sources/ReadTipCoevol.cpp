
#include "Sample.h"
#include "TipCoevolModel.h"
#include "MeanValTree.h"
#include "StringStreamUtils.h"
#include "MeanChronogram.h"
#include "MeanCovMatrix.h"
#include "MeanICTree.h"

class LogNSample : public Sample	{

	private:
	string modeltype;
	string treefile;
	string nucdatafile;
	string morpho2datafile;
	string morpho3datafile;
	string morpho4datafile;
	string contdatafile;
	int contdatatype;
	string calibfile;
	double rootage;
	double rootstdev;

	double divrate;
	double extrate;
	double massext;
	double K0;
	double K1;
	double divratestdev;
	double extratestdev;
	double massextstdev;
	double K0stdev;
	double K1stdev;
	double T1;

	double bddivratemean;
	double mumean;
	double psimean;
	double rhomean;
	double bddivratestdev;
	double mustdev;
	double psistdev;
	double rhostdev;

	double divcutoff;
	int Nextant;
	int chronoprior;
	int conjpath;
	int clampdiag;
	int clamptree;
	int meanexp;
	int df;

	int brownian;
	int nSegments;
	Segmentation segm;

	int clockprior;
	double ratemean;
	double ratestdev;
	double morphoratemean;
	double morphoratestdev;

	int scaling;
	double scalet0val;
	double scalefactorval;
	double scalefactorstdev;
	double scalerateval;
	double scaleratestdev;

	int prior;

	public:

	string GetModelType() {return modeltype;}

	GTRLogNormalModel* GetModel() {return (GTRLogNormalModel*) model;}

	LogNSample(string filename, int inburnin, int inevery, int inuntil) : Sample(filename,inburnin,inevery,inuntil)	{
		Open();
	}

	void Open()	{

		ifstream is((name + ".param").c_str());

		// check that file exists
		if (!is)	{
			cerr << "error : cannot find file : " << name << ".param\n";
			exit(1);
		}

		brownian = 0;
		clockprior = 0;

		is >> modeltype;
		is >> treefile >> nucdatafile;
		is >> morpho2datafile;
		is >> morpho3datafile;
		is >> morpho4datafile;
		is >> contdatafile;
		is >> contdatatype;
		is >> calibfile;

		is >> rootage >> rootstdev;
		is >> divrate >> extrate >> massext >> K0 >> K1;
		is >> divratestdev >> extratestdev >> massextstdev >> K0stdev >> K1stdev;
		is >> T1;
		is >> bddivratemean >> mumean >> psimean >> rhomean >> bddivratestdev >> mustdev >> psistdev >> rhostdev;

		is >> divcutoff;
		is >> Nextant;

		is >> chronoprior;
		is >> conjpath;
		is >> clampdiag >> clamptree >> meanexp >> df;

		brownian = 0;
		clockprior = 0;
		scaling = 0;
		prior = 0;

		int check;
		is >> check;
		if (check)	{
			is >> brownian >> nSegments;
			int tmp;
			is >> tmp;
			if (tmp == 0)	{
				segm = SEGM_REGULAR;
			}
			else	{
				segm = SEGM_ABSOLUTE;
			}
			is >> check;
			if (check)	{
				is >> clockprior;
				is >> ratemean >> ratestdev;
				is >> morphoratemean >> morphoratestdev;
				is >> check;
				if (check)	{
					is >> scaling;
					is >> scalet0val;
					is >> scalefactorval >> scalefactorstdev;
					is >> scalerateval >> scaleratestdev;
					is >> check;
					if (check)	{
						is >> prior;
						is >> check;
						if (check)	{
							cerr << "error when reading model\n";
							exit(1);
						}
					}
				}
			}
		}

		is >> chainevery >> chainuntil >> chainsize;

		if (modeltype == "TIPCOEVOL")	{
			model = new GTRLogNormalModel(nucdatafile,morpho2datafile,morpho3datafile,morpho4datafile,contdatafile,contdatatype,treefile,calibfile,rootage,rootstdev,divrate,extrate,massext,K0,K1,divratestdev,extratestdev,massextstdev,K0stdev,K1stdev,T1,bddivratemean,mumean,psimean,rhomean,bddivratestdev,mustdev,psistdev,rhostdev,divcutoff,Nextant,chronoprior,clampdiag,clamptree,meanexp,df,brownian,nSegments,segm,clockprior,ratemean,ratestdev,morphoratemean,morphoratestdev,scaling,scalet0val,scalefactorval,scalefactorstdev,scalerateval,scaleratestdev,prior,false);
		}
		else	{
			cerr << "error when opening file "  << name << " : does not recognise model type : " << modeltype << '\n';
			exit(1);
		}

		model->FromStream(is);
		// model->Update();

		OpenChainFile();

		cerr << "number of points to read : " << size << '\n';
		cerr << '\n';
	}

	void GetTrees(string name)	{

		for (int i=0; i<size; i++)	{
			cerr << '.';
			GetNextPoint();

			ostringstream c;
			c << name << i << ".treeparam";
			ofstream cos(c.str().c_str());

			cos << *(GetModel()->GetChronogram());
			cos << *(GetModel()->GetChronogram()->GetScale());

			ostringstream s;
			s << name << i << ".tree";
			ofstream os(s.str().c_str());

			MeanChronogram* meanchrono = new MeanChronogram(GetModel()->GetTree(),false,false,false,true);
			meanchrono->Add(GetModel()->GetChronogram());
			meanchrono->Normalise();
			meanchrono->ToStream(os);

			ostringstream t;
			t << name << i << ".dates.tab";
			ofstream tos(t.str().c_str());
			meanchrono->Tabulate(tos);

			delete meanchrono;
		}
		cerr << '\n';
		cerr << "samples in " << name << "X.tree\n";
		cerr << '\n';
	}

	void Simulate(string name, double minlength, double maxlength, double mintotallength, double maxtotallength, string treefile, string paramfile)	{

		for (int i=0; i<size; i++)	{
			cerr << '.';
			GetNextPoint();
			ostringstream s;
			s << name << i;
			if ((treefile != "null") && (paramfile != "null"))	{
				GetModel()->SimulateFromFixedParam(s.str(),minlength,maxlength,mintotallength,maxtotallength,treefile,paramfile,i);
			}
			else	{
				GetModel()->Simulate(s.str(),minlength,maxlength,mintotallength,maxtotallength);
			}
			/*
			for (int j=0; j<nrep; j++)	{
				ostringstream s;
				s << name << j;
				s << name << i << j;
				GetModel()->Simulate(s.str());
			}
			*/
		}
		cerr << '\n';
	}

	void ReadIC()	{

		int dim = GetModel()->GetMultiVariateProcess()->GetDim();

		MeanICTree* ictree = new MeanICTree(GetModel()->GetMultiVariateProcess());
		ictree->Reset();
		
		/*
		MeanICTree* ppictree = new MeanICTree(GetModel()->GetMultiVariateProcess());
		ppictree->Reset();
		double* pp = new double[dim];
		for (int k=0; k<dim; k++)	{
			pp[k] = 0;
		}
		*/

		for (int i=0; i<size; i++)	{
			cerr << '.';
			GetNextPoint();
			ictree->Add();
			// GetModel()->Update();
			/*
			GetModel()->SampleProcess();
			ppictree->Add();
			for (int k=0; k<dim; k++)	{
				if (ppictree->GetGmean(k) > ictree->GetGmean(k))	{
					pp[k]++;
				}
			}
			*/
		}
		cerr << '\n';
		ictree->Normalise();
		ofstream os((name + ".postmeanic").c_str());
		ictree->Tabulate(os,false);
		ofstream los((name + ".postmeanleafic").c_str());
		ictree->Tabulate(los,true);
		cerr << "mean independent contrasts in " << name  << ".postmeanic\n";
		cerr << "for terminal branches only in " << name  << ".postmeanleafic\n";
		/*
		ppictree->Normalise();
		cerr << "pp positive trend : \n";
		for (int k=0; k<dim; k++)	{
			cerr << k << '\t' << pp[k] / size << '\n';
		}
		cerr << '\n';
		*/
		ictree->DAgostinosKandZ();
	}

	void ReadAcrossKT()	{

		double ncross = 0;

		for (int i=0; i<size; i++)	{
			cerr << '.';
			// get next point -> will be stored into "model", and thus, will be accessible through GetModel()
			GetNextPoint();
			GetModel()->FastUpdate();
			ncross += GetModel()->GetChronogram()->GetHowManyCross(65.0);
		}
			
		ncross /= size;
		cerr << '\n';
		cerr << "number of lineages crossing KT: " << ncross << '\n';
		cerr << '\n';
	}

	void ReadRate()	{

		double totvar = 0;
		double lnvar = 0;
		double wnvar = 0;

		double totmvar = 0;
		double lnmvar = 0;
		double wnmvar = 0;

		for (int i=0; i<size; i++)	{
			cerr << '.';

			GetNextPoint();
			GetModel()->FastUpdate();

			totvar += GetModel()->GetTotalRateVariance();
			lnvar += GetModel()->GetLnRateVariance();
			wnvar += GetModel()->GetWnRateVariance();
			/*
			double tmp1 = GetModel()->GetLnRateVariance();
			double tmp2 = GetModel()->GetWnRateVariance();
			cout << tmp1 / (tmp1 + tmp2) << '\t' << tmp2 / (tmp1 + tmp2) << '\n';
			*/

			totmvar += GetModel()->GetTotalMorphoRateVariance();
			lnmvar += GetModel()->GetLnMorphoRateVariance();
			wnmvar += GetModel()->GetWnMorphoRateVariance();

		}
		cerr << '\n';

		totvar /= size;
		lnvar /= size;
		wnvar /= size;
		cerr << "nuc rate variance\n";
		cerr << "ln    : " << lnvar << '\t' << lnvar / (lnvar + wnvar) << '\n';
		cerr << "wn    : " << wnvar << '\t' << wnvar / (lnvar + wnvar) << '\n';
		cerr << "total : " << totvar << '\t' << lnvar + wnvar << '\n';
		cerr << '\n';

		ofstream os((name + ".ratevar").c_str());
		os << lnvar / (lnvar + wnvar) << '\t' << wnvar / (lnvar + wnvar) << '\n';

		totmvar /= size;
		lnmvar /= size;
		wnmvar /= size;
		cerr << "morpho rate variance\n";
		cerr << "ln    : " << lnmvar << '\t' << lnmvar / (lnmvar + wnmvar) << '\n';
		cerr << "wn    : " << wnmvar << '\t' << wnmvar / (lnmvar + wnmvar) << '\n';
		cerr << "total : " << totmvar << '\t' << lnmvar + wnmvar << '\n';
		cerr << '\n';

	}

	void ReadRateTrees(bool printlog, bool printmean, bool printmed, bool printci, bool printstdev, bool withleaf, bool withinternal)	{

		if (prior == 1)	{
			cerr << "error in read rate trees: chain under the prior\n";
			exit(1);
		}

		MeanExpNormTree* meannucratetree =  new MeanExpNormTree(GetModel()->GetTree(),false,printlog,printmean,printmed,printci,printstdev,withleaf,withinternal);
		MeanExpNormTree* meanwnnucratetree =  new MeanExpNormTree(GetModel()->GetTree(),false,printlog,printmean,printmed,printci,printstdev,withleaf,withinternal);
		MeanExpNormTree* meanmorphoratetree = 0;
		if (GetModel()->Morpho())	{
			meanmorphoratetree =  new MeanExpNormTree(GetModel()->GetTree(),false,printlog,printmean,printmed,printci,printstdev,withleaf,withinternal);
		}

		int L = GetModel()->GetL();
		int Ncont = GetModel()->GetNcont();
		int dim = 0;
		if (! prior)	{
			dim = GetModel()->GetDim();
		}
		MeanExpNormTree** tree = new MeanExpNormTree*[dim];
		for (int k=0; k<dim; k++)	{
			tree[k] = new MeanExpNormTree(GetModel()->GetTree(),false,printlog,printmean,printmed,printci,printstdev,withleaf,withinternal);
		}


		for (int i=0; i<size; i++)	{

			cerr << '.';

			// get next point -> will be stored into "model", and thus, will be accessible through GetModel()

			GetNextPoint();
			GetModel()->FastUpdate();

			meannucratetree->Add((BranchVarTree<PosReal>*) GetModel()->GetNucRateTree(), GetModel()->GetCalibratedChronogram(),true);
			if (clockprior)	{
				if (clockprior == 1)	{
					meanwnnucratetree->Add((BranchVarTree<PosReal>*) GetModel()->GetWnNucRateTree(), GetModel()->GetCalibratedChronogram(),true);
				}
				else	{
					meanwnnucratetree->Add((BranchVarTree<PosReal>*) GetModel()->GetWnNucRateTree(), GetModel()->GetCalibratedChronogram(),false);
				}
			}
				

			if (GetModel()->Morpho())	{
				meanmorphoratetree->Add((BranchVarTree<PosReal>*) GetModel()->GetMorphoRateTree(), GetModel()->GetCalibratedChronogram(),true);
			}

			for (int k=0; k<dim; k++)	{
				tree[k]->Add(GetModel()->GetMultiVariateProcess(),GetModel()->GetChronogram(),k);
			}

		}

		cerr << '\n';

		meannucratetree->Normalise();
		ofstream ros((name + ".nucrates").c_str());
		meannucratetree->ToStream(ros);
		cerr << "nuc rates in " << name << ".nucrates\n";

		if (GetModel()->GetWnNucRateTree())	{
			meanwnnucratetree->Normalise();
			meanwnnucratetree->CutoffFromBelow(1e-4);
			ofstream ros((name + ".wnnucrates").c_str());
			meanwnnucratetree->ToStream(ros);
			cerr << "wn nuc rates in " << name << ".wnnucrates\n";
		}

		if (GetModel()->Morpho())	{
			meanmorphoratetree->Normalise();
			ofstream ros((name + ".morphorates").c_str());
			meanmorphoratetree->ToStream(ros);
			cerr << "morpho rates in " << name << ".morphorates\n";
		}

		for (int k=0; k<dim; k++)	{
			tree[k]->Normalise();
			ostringstream s;
			s << GetName() << ".postmean" << k+1 << ".tre";
			ofstream os(s.str().c_str());
			tree[k]->ToStream(os);
		}

	}

	void Read(bool printlog, bool printmean, bool printmed, bool printci, bool printstdev, bool withleaf, bool withinternal, string mulreg)	{

		cerr << "read\n";

		MeanChronogram* meanchrono = new MeanChronogram(GetModel()->GetTree(),true,false,false,true);
		// meanchrono->SetWithLeafDates(true);

		int L = GetModel()->GetL();
		int Ncont = GetModel()->GetNcont();
		int dim = 0;
		if (! prior)	{
			dim = GetModel()->GetDim();
		}
		MeanExpNormTree** tree = new MeanExpNormTree*[dim];
		for (int k=0; k<dim; k++)	{
			tree[k] = new MeanExpNormTree(GetModel()->GetTree(),false,printlog,printmean,printmed,printci,printstdev,withleaf,withinternal);
		}

		MeanCovMatrix*  mat = 0;
		if (dim)	{
			mat = new MeanCovMatrix(dim);
		}

		double meanwnsigma = 0;
		int nstate = GetModel()->GetNstate();
		int nrr = nstate * (nstate -1) / 2;
		double* meanrelrate = new double[nrr];
		for (int i=0; i<nrr; i++)	{
			meanrelrate[i] = 0;
		}
		double* meanstat = new double[nstate];
		for (int i=0; i<nstate; i++)	{
			meanstat[i] = 0;
		}

		// cycle over the sample
		for (int i=0; i<size; i++)	{
			cerr << '.';

			// get next point -> will be stored into "model", and thus, will be accessible through GetModel()

			GetNextPoint();
			GetModel()->FastUpdate();

			if (clockprior)	{
				meanwnsigma += GetModel()->GetWnNucSigma();
			}

			if (! prior)	{
				double* relrate = GetModel()->GetRelRate();
				for (int j=0; j<nrr; j++)	{
					meanrelrate[j] += relrate[j];
				}
				double* stat = GetModel()->GetStationary();
				for (int j=0; j<nstate; j++)	{
					meanstat[j] += stat[j];
				}
			}

			// GetModel()->GetChronogram()->specialUpdate();
			meanchrono->Add(GetModel()->GetCalibratedChronogram());

			if (! prior)	{

				for (int k=0; k<dim; k++)	{
					tree[k]->Add(GetModel()->GetMultiVariateProcess(),GetModel()->GetChronogram(),k);
				}

				if (dim)	{
					CovMatrix& m = *(GetModel()->GetCovMatrix());
					mat->Add(&m);
				}
			}

		}
		cerr << '\n';
		cerr << "normalise\n";

		meanwnsigma /= size;
		for (int j=0; j<nrr; j++)	{
			meanrelrate[j] /= size;
		}
		for (int j=0; j<nstate; j++)	{
			meanstat[j] /= size;
		}

		meanchrono->Normalise();
		ofstream dos((name + ".dates.tre").c_str());
		meanchrono->ToStream(dos);
		cerr << "dates in " << name << ".dates.tre\n";
		ofstream ddos((GetName() + ".dates.tab").c_str());
		meanchrono->Tabulate(ddos);

		if (! prior)	{

		for (int k=0; k<dim; k++)	{
			tree[k]->Normalise();
			ostringstream s;
			s << GetName() << ".postmean" << k+1 << ".tre";
			ofstream os(s.str().c_str());
			tree[k]->ToStream(os);
			// cerr << "reconstructed variation in " << GetModel()->GetContinuousData()->GetCharacterName(k) << " in "  << name << ".postmean" << k+1 << ".tre\n";
		}

		ofstream pos((name + ".meanparam").c_str());

		if (dim)	{
		mat->Normalize();
		if (mulreg != "")	{
			ofstream cout((GetName() + ".marginalcov").c_str());
			if (mulreg.size() != ((unsigned int) mat->GetDim()))	{
				cerr << "error when specifying multiple regression : " << mulreg << '\n';
				exit(1);
			}

			int* indexarray = new int[mat->GetDim()];
			int dim = 0;
			for (int l=0; l<mat->GetDim(); l++)	{
				if (mulreg[l] == '1')	{
					indexarray[l] = 1;
					dim++;
				}
				else if (mulreg[l] == '0')	{
					indexarray[l] = 0;
				}
				else	{
					cerr << "error when specifying multiple regression : " << mulreg << '\n';
					exit(1);
				}
			}

			cout << "entries are in the following order:\n";
			cout << "nucrate\n";
			if (GetModel()->Morpho())	{
				cout << "morphorate\n";
			}
			for (int k=0; k<Ncont; k++)	{
				cout << "trait " << k << '\n';
			}

			cout << '\n';
			PartialMeanCovMatrix* partmat = new PartialMeanCovMatrix(mat,indexarray,dim);

			ReducedMeanCovMatrix* redmat = new ReducedMeanCovMatrix(mat,indexarray,dim);
			cout << '\n';
			cout << "partial correl\n";
			cout << '\n';
			cout << *partmat;

			ofstream cout2((GetName() + ".controlcov").c_str());
			cout2 << "entries are in the following order:\n";
			cout2 << "nucrate\n";
			if (GetModel()->Morpho())	{
				cout2 << "morphorate\n";
			}
			for (int k=0; k<Ncont; k++)	{
				cout2 << "trait " << k << '\n';
			}

			cout2 << '\n';
			cout2 << "reduced correl\n";
			cout2 << '\n';
			redmat->PrintCovariances(cout2);
			redmat->PrintCorrel(cout2);
			redmat->PrintPosteriorProbs(cout2);
			// cout2 << *redmat;
			cout2 << '\n';
			delete partmat;
			delete redmat;
		}
		ofstream cout((GetName() + ".cov").c_str());
		cout << "entries are in the following order:\n";
		cout << "nucrate\n";
		if (GetModel()->Morpho())	{
			cout << "morphorate\n";
		}
		for (int k=0; k<Ncont; k++)	{
			cout << "trait " << k << '\n';
		}
		cout << '\n';
		// mat->SetLatex(tex);
		cout << *mat;

		cerr << "covariance matrix in " << name << ".cov\n";
		cerr << '\n';

		mat->PrintCovariances(pos);
		}

		if (clockprior)	{
			pos << meanwnsigma << '\n';
		}
		pos << nrr << '\t';
		for (int j=0; j<nrr; j++)	{
			pos << meanrelrate[j] << '\t';
		}
		pos << '\n';
		pos << nstate << '\t';
		for (int j=0; j<nstate; j++)	{
			pos << meanstat[j] << '\t';
		}
		pos << '\n';

		}

		cerr << '\n';
	}
};


int main(int argc, char* argv[])	{

	int burnin = 0;
	int every = 1;
	int until = -1;
	string name;

	string mulreg = "";

	bool printlog = false;
	bool printmean = false;
    bool printmed = false;
	bool printci = true;
	bool printstdev = false;
	bool withleaf = true;
	bool withinternal = true;

	int ic = 0;
	int simu = 0;
	string simuname = "";
	double mintotallength = 0;
	double maxtotallength = 50;
	double minlength = 0;
	double maxlength = 2;
	string treefile = "null";
	string paramfile = "null";

	int trees = 0;
	string treename = "";

	int nrep = 1;

	int rate = 0;
	int ratetree = 0;

	int kt = 0;

	try	{

		if (argc == 1)	{
			throw(0);
		}

		int i = 1;
		while (i < argc)	{
			string s = argv[i];

			if (s == "-mulreg")	{
				i++;
				mulreg = argv[i];
			}
			else if (s == "-partial")	{
				i++;
				mulreg = argv[i];
			}
			else if (s == "+log")	{
				printlog = true;
			}
			else if (s == "-log")	{
				printlog = false;
			}
			else if (s == "+mean")	{
				printmean = true;
                printmed = false;
			}
			else if (s == "-mean")	{
				printmean = false;
			}
			else if ((s == "+med") || (s == "+median"))	{
				printmed = true;
                printmean = false;
			}
			else if ((s == "-med") || (s == "-median"))	{
				printmed = false;
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
			else if (s == "-ic")	{
				ic = 1;
			}
			else if (s == "-simu")	{
				simu = 1;
				i++;
				simuname = argv[i];
			}
			else if (s == "-nrep")	{
				i++;
				nrep = atoi(argv[i]);
			}
			else if (s == "-minmax")	{
				i++;
				minlength = atof(argv[i]);
				i++;
				maxlength = atof(argv[i]);
				i++;
				mintotallength = atof(argv[i]);
				i++;
				maxtotallength = atof(argv[i]);
			}
			else if (s == "-t")	{
				i++;
				treefile = argv[i];
			}
			else if (s == "-p")	{
				i++;
				paramfile = argv[i];
			}
			else if (s == "-trees")	{
				trees = 1;
				i++;
				treename = argv[i];
			}
			else if (s == "-rate")	{
				rate = 1;
			}
			else if (s == "-ratetree")	{
				ratetree = 1;
			}
			else if (s == "-kt")	{
				kt = 1;
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
		cerr << "readlogn [-x <burnin> <every> <until>] <chainname> \n";
		cerr << '\n';
		exit(1);
	}

	LogNSample sample(name,burnin,every,until);

	if (ic)	{
		sample.ReadIC();
	}
	else if (trees)	{
		sample.GetTrees(treename);
	}
	else if (simu)	{
		sample.Simulate(simuname,minlength,maxlength,mintotallength,maxtotallength,treefile,paramfile);
	}
	else if (rate)	{
		sample.ReadRate();
	}
	else if (kt)	{
		sample.ReadAcrossKT();
	}
	else if (ratetree)	{
		sample.ReadRateTrees(printlog,printmean,printmed,printci,printstdev,withleaf,withinternal);
	}
	else	{
		sample.Read(printlog,printmean,printmed,printci,printstdev,withleaf,withinternal,mulreg);
	}

}



