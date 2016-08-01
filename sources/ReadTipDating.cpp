
#include "Sample.h"
#include "TipDatingModel.h"
#include "MeanValTree.h"
#include "StringStreamUtils.h"
#include "MeanChronogram.h"
#include "MeanICTree.h"

class LogNSample : public Sample	{

	private:
	string modeltype;
	string treefile;
	string nucdatafile;
	string morpho2datafile;
	string morpho3datafile;
	string morpho4datafile;
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
	int clockprior;
	double ratemean;
	double ratestdev;
	double morphoratemean;
	double morphoratestdev;
	int conjpath;

	int clamptree;
	int prior;
	int nSegments;
	double fixalpha;

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

		is >> modeltype;
		is >> treefile >> nucdatafile;
		is >> morpho2datafile;
		is >> morpho3datafile;
		is >> morpho4datafile;
		is >> calibfile;
		is >> rootage >> rootstdev;
		is >> divrate >> extrate >> massext >> K0 >> K1;
		is >> divratestdev >> extratestdev >> massextstdev >> K0stdev >> K1stdev;
		is >> T1;
		is >> bddivratemean >> mumean >> psimean >> rhomean >> bddivratestdev >> mustdev >> psistdev >> rhostdev;
		is >> divcutoff;
		is >> Nextant;
		is >> chronoprior;
		is >> clockprior;
		is >> ratemean >> ratestdev;
		is >> morphoratemean >> morphoratestdev;
		is >> conjpath;

		clamptree = 0;
		prior = 0;
		nSegments = 100;
		fixalpha = 0;

		int check;
		is >> check;
		if (check)	{
			is >> clamptree;
			is >> prior;
			is >> nSegments;
			is >> fixalpha;
			is >> check;
			if (check)	{
				cerr << "error when reading model\n";
				cerr << divcutoff << '\t' << chronoprior << '\t' << conjpath << '\n';
				exit(1);
			}
		}

		is >> chainevery >> chainuntil >> chainsize;

		if (modeltype == "TIP")	{
			model = new GTRLogNormalModel(nucdatafile,morpho2datafile,morpho3datafile,morpho4datafile,treefile,calibfile,rootage,rootstdev,divrate,extrate,massext,K0,K1,divratestdev,extratestdev,massextstdev,K0stdev,K1stdev,T1,bddivratemean,mumean,psimean,rhomean,bddivratestdev,mustdev,psistdev,rhostdev,divcutoff,Nextant,chronoprior,clockprior,ratemean,ratestdev,morphoratemean,morphoratestdev,clamptree,prior,nSegments,fixalpha,false);
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

	void Simulate(string name, double minlength, double maxlength)	{

		for (int i=0; i<size; i++)	{
			cerr << '.';
			GetNextPoint();
			ostringstream s;
			s << name << i;
			GetModel()->Simulate(s.str(),minlength,maxlength);
		}
		cerr << '\n';
	}

	void ReadIC()	{

		MeanICTree* ictree = 0;
		// MeanICTree* ictree = new MeanICTree(GetModel()->GetAlphaStableProcess());
		// MeanICTree* ppictree = new MeanICTree(GetModel()->GetMultiVariateProcess());
		ictree->Reset();
		// ppictree->Reset();
		int dim = 1;
		/*
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

	void ReadRateTrees(bool printlog, bool printmean, bool printci, bool printstdev, bool withleaf, bool withinternal)	{

		MeanExpNormTree* meannucratetree =  new MeanExpNormTree(GetModel()->GetTree(),false,printlog,printmean,printci,printstdev,withleaf,withinternal);
		MeanExpNormTree* meanwnnucratetree =  new MeanExpNormTree(GetModel()->GetTree(),false,printlog,printmean,printci,printstdev,withleaf,withinternal);
		for (int i=0; i<size; i++)	{

			cerr << '.';

			GetNextPoint();
			GetModel()->FastUpdate();

			meannucratetree->Add((BranchVarTree<PosReal>*) GetModel()->GetNucRateTree(), GetModel()->GetCalibratedChronogram(),true);
			meanwnnucratetree->Add((BranchVarTree<PosReal>*) GetModel()->GetWnNucRateTree(), GetModel()->GetCalibratedChronogram(),false);
		}

		cerr << '\n';

		meannucratetree->Normalise();
		ofstream ros((name + ".nucrates").c_str());
		meannucratetree->ToStream(ros);
		cerr << "nuc rates in " << name << ".nucrates\n";

		meanwnnucratetree->Normalise();
		meanwnnucratetree->CutoffFromBelow(1e-4);
		ofstream wos((name + ".wnnucrates").c_str());
		meanwnnucratetree->ToStream(wos);
		cerr << "wn nuc rates in " << name << ".wnnucrates\n";

	}

	void ReadRate()	{

		double totvar = 0;
		double lnvar = 0;
		double wnvar = 0;

		/*
		double totmvar = 0;
		double lnmvar = 0;
		double wnmvar = 0;
		*/

		for (int i=0; i<size; i++)	{
			cerr << '.';

			GetNextPoint();
			GetModel()->FastUpdate();

			totvar += GetModel()->GetTotalRateVariance();
			lnvar += GetModel()->GetLnRateVariance();
			wnvar += GetModel()->GetWnRateVariance();

			/*
			totmvar += GetModel()->GetTotalMorphoRateVariance();
			lnmvar += GetModel()->GetLnMorphoRateVariance();
			wnmvar += GetModel()->GetWnMorphoRateVariance();
			*/

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

		/*
		totmvar /= size;
		lnmvar /= size;
		wnmvar /= size;
		cerr << "morpho rate variance\n";
		cerr << "ln    : " << lnmvar << '\t' << lnmvar / (lnmvar + wnmvar) << '\n';
		cerr << "wn    : " << wnmvar << '\t' << wnmvar / (lnmvar + wnmvar) << '\n';
		cerr << "total : " << totmvar << '\t' << lnmvar + wnmvar << '\n';
		cerr << '\n';
		*/

	}

	void Read(bool printlog, bool printmean, bool printci, bool printstdev, bool withleaf, bool withinternal)	{

		cerr << "read\n";

		// MeanExpNormTree* meanratetree =  new MeanExpNormTree(GetModel()->GetTree(),false,printlog,printmean,printci,printstdev,withleaf,withinternal);
		MeanChronogram* meanchrono = new MeanChronogram(GetModel()->GetTree(),true,false,false,true);
		// meanchrono->SetWithLeafDates(true);

		double ncross = 0;

		// cycle over the sample
		for (int i=0; i<size; i++)	{
			cerr << '.';

			// get next point -> will be stored into "model", and thus, will be accessible through GetModel()

			GetNextPoint();

			// meanratetree->Add(GetModel()->GetLengthTree(), GetModel()->GetCalibratedChronogram());
			// meanratetree->Add((BranchVarTree<PosReal>*) GetModel()->GetLogNormalProcess(), GetModel()->GetCalibratedChronogram());
			meanchrono->Add(GetModel()->GetCalibratedChronogram());

			ncross += GetModel()->GetCalibratedChronogram()->GetHowManyCross(65.0);
			
		}
		cerr << '\n';
		cerr << "normalise\n";

		/*
		meanratetree->Normalise();
		ofstream ros((name + ".rates").c_str());
		meanratetree->ToStream(ros);
		cerr << "rates in " << name << ".rates\n";
		*/

		meanchrono->Normalise();
		ofstream dos((name + ".dates.tre").c_str());
		meanchrono->ToStream(dos);
		cerr << "dates in " << name << ".dates.tre\n";
		ofstream ddos((GetName() + ".dates.tab").c_str());
		meanchrono->Tabulate(ddos);

		ncross /= size;
		cerr << '\n';
		cerr << "number of lineages crossing KT: " << ncross << '\n';
		cerr << '\n';

		cerr << '\n';
	}
};


int main(int argc, char* argv[])	{

	int burnin = 0;
	int every = 1;
	int until = -1;
	string name;

	bool printlog = false;
	bool printmean = false;
	bool printci = true;
	bool printstdev = false;
	bool withleaf = true;
	bool withinternal = true;

	int simu = 0;
	string simuname = "";
	double minlength = 0;
	double maxlength = 50;

	int ratetree = 0;
	int rate = 0;

	try	{

		if (argc == 1)	{
			throw(0);
		}

		int i = 1;
		while (i < argc)	{
			string s = argv[i];

			if (s == "+log")	{
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
			else if (s == "-simu")	{
				simu = 1;
				i++;
				simuname = argv[i];
			}
			else if (s == "-minmax")	{
				i++;
				minlength = atof(argv[i]);
				i++;
				maxlength = atof(argv[i]);
			}
			else if (s == "-ratetree")	{
				ratetree = 1;
			}
			else if (s == "-rate")	{
				rate = 1;
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

	if (simu)	{
		sample.Simulate(simuname,minlength,maxlength);
	}
	else if (ratetree)	{
		sample.ReadRateTrees(printlog,printmean,printci,printstdev,withleaf,withinternal);
	}
	else if (rate)	{
		sample.ReadRate();
	}
	else	{
		sample.Read(printlog,printmean,printci,printstdev,withleaf,withinternal);
	}
}



