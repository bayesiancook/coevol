

#include "Chain.h"
#include "TipDatingModel.h"
#include "StringStreamUtils.h"

class LogNChain : public Chain	{

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

	LogNChain(string intree, string innucdata, string inmorpho2data, string inmorpho3data, string inmorpho4data, string incalibfile, double inrootage, double inrootstdev, double indivrate, double inextrate, double inmassext, double inK0, double inK1, double indivratestdev, double inextratestdev, double inmassextstdev, double inK0stdev, double inK1stdev, double inT1, double inbddivratemean, double inmumean, double inpsimean, double inrhomean, double inbddivratestdev, double inmustdev, double inpsistdev, double inrhostdev, double indivcutoff, int inNextant, int inchronoprior, int inclockprior, double inratemean, double inratestdev, double inmorphoratemean, double inmorphoratestdev, bool inconjpath, int inclamptree, int inprior, int innSegments, double infixalpha, string filename, int force = 1)	{
		modeltype = "TIP";
		treefile = intree;
		nucdatafile = innucdata;
		morpho2datafile = inmorpho2data;
		morpho3datafile = inmorpho3data;
		morpho4datafile = inmorpho4data;
		rootage = inrootage;
		rootstdev = inrootstdev;

		divrate = indivrate;
		extrate = inextrate;
		massext = inmassext;
		K0 = inK0;
		K1 = inK1;
		divratestdev = indivratestdev;
		extratestdev = inextratestdev;
		massextstdev = inmassextstdev;
		K0stdev = inK0stdev;
		K1stdev = inK1stdev;
		T1 = inT1;

		bddivratemean = inbddivratemean;
		mumean = inmumean;
		psimean = inpsimean;
		rhomean = inrhomean;

		bddivratestdev = inbddivratestdev;
		mustdev = inmustdev;
		psistdev = inpsistdev;
		rhostdev = inrhostdev;

		divcutoff = indivcutoff;
		Nextant = inNextant;
		chronoprior = inchronoprior;
		clockprior = inclockprior;
		ratemean = inratemean;
		ratestdev = inratestdev;
		morphoratemean = inmorphoratemean;
		morphoratestdev = inmorphoratestdev;

		calibfile = incalibfile;

		conjpath = inconjpath;

		clamptree = inclamptree;
		prior = inprior;

		nSegments = innSegments;
		fixalpha = infixalpha;

		name = filename;
		New(force);
	}

	LogNChain(string filename)	{
		name = filename;
		Open();
	}

	void New(int force)	{
		if (modeltype == "TIP")	{
			model = new GTRLogNormalModel(nucdatafile,morpho2datafile,morpho3datafile,morpho4datafile,treefile,calibfile,rootage,rootstdev,divrate,extrate,massext,K0,K1,divratestdev,extratestdev,massextstdev,K0stdev,K1stdev,T1,bddivratemean,mumean,psimean,rhomean,bddivratestdev,mustdev,psistdev,rhostdev,divcutoff,Nextant,chronoprior,clockprior,ratemean,ratestdev,morphoratemean,morphoratestdev,clamptree,prior,nSegments,fixalpha);
		}
		else	{
			cerr << "error, does not recognise model type : " << modeltype << '\n';
			exit(1);
		}
		Reset(force);
		// model->Update();
		cerr << "ln prob = " << model->GetLogProb() << "\n";
	}

	void Open()	{
		ifstream is((name + ".param").c_str());
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
				exit(1);
			}
		}

		is >> every >> until >> size;

		if (modeltype == "TIP")	{
			model = new GTRLogNormalModel(nucdatafile,morpho2datafile,morpho3datafile,morpho4datafile,treefile,calibfile,rootage,rootstdev,divrate,extrate,massext,K0,K1,divratestdev,extratestdev,massextstdev,K0stdev,K1stdev,T1,bddivratemean,mumean,psimean,rhomean,bddivratestdev,mustdev,psistdev,rhostdev,divcutoff,Nextant,chronoprior,clockprior,ratemean,ratestdev,morphoratemean,morphoratestdev,clamptree,prior,nSegments,fixalpha);
		}
		else	{
			cerr << "error when opening file "  << name << " : does not recognise model type : " << modeltype << '\n';
			exit(1);
		}
		model->FromStream(is);
		cerr << "update\n";
		model->Update();
		cerr << size << " points saved, current ln prob = " << GetModel()->GetLogProb() << "\n";
	}

	void Save()	{
		ofstream param_os((name + ".param").c_str());
		param_os << GetModelType() << '\n';
		param_os << treefile << '\t' << nucdatafile << '\n';
		param_os << morpho2datafile << '\n';
		param_os << morpho3datafile << '\n';
		param_os << morpho4datafile << '\n';
		param_os << calibfile << '\n';

		param_os << rootage << '\t' << rootstdev << '\n';
		param_os << divrate << '\t' << extrate << '\t' << massext << '\t' << K0 << '\t' << K1 << '\n';
		param_os << divratestdev << '\t' << extratestdev << '\t' << massextstdev << '\t' << K0stdev << '\t' << K1stdev << '\n';
		param_os << T1 << '\n';
		param_os << bddivratemean << '\t' << mumean << '\t' << psimean << '\t' << rhomean << '\t';
		param_os << bddivratestdev << '\t' << mustdev << '\t' << psistdev << '\t' << rhostdev << '\n';

		param_os << divcutoff << '\n';
		param_os << Nextant << '\n';

		param_os << chronoprior << '\n';
		param_os << clockprior << '\n';
		param_os << ratemean << '\t' << ratestdev << '\t' << morphoratemean << '\t' << morphoratestdev << '\n';
		param_os << conjpath << '\n';
		param_os << 1 << '\n';
		param_os << clamptree << '\n';
		param_os << prior << '\n';
		param_os << nSegments << '\n';
		param_os << fixalpha << '\n';
		param_os << 0 << '\n';
		param_os << every << '\t' << until << '\t' << size << '\n';
		model->ToStream(param_os);
	}

	void Move()	{
		for (int i=0; i<every; i++)	{
			model->Move(1);
		}
		SavePoint();
		Save();
		Monitor();
	}
};

int main(int argc, char* argv[])	{

	if (argc == 2)	{
		string name = argv[1];
		LogNChain* chain = new LogNChain(name);
		cerr << "start\n";
		cerr << '\n';
		chain->Start();
		cerr << '\n';
		cerr << "exit\n";
		cerr << chain->GetSize() << " points saved, current ln prob = " << chain->GetModel()->GetLogProb() << "\n";
		cerr << '\n';
	}
	else	{

		string treefile = "";
		string nucdatafile = "";
		string morpho2datafile = "null";
		string morpho3datafile = "null";
		string morpho4datafile = "null";
		string calibfile = "null";

		double rootmean = 1;
		double rootstdev = 1;;

		int chronoprior = 0;
		double divcutoff = 0;
		int Nextant = 0;

		double divrate = 0;
		double extrate = 0;
		double massext = 0;
		double K0 = 0;
		double K1 = 0;
		double divratestdev = 0;
		double extratestdev = 0;
		double massextstdev = 0;
		double K0stdev = 0;
		double K1stdev = 0;
		double T1 = 0;

		double bddivratemean = 0;
		double mumean = 0;
		double psimean = 0;
		double rhomean = 0;
		double bddivratestdev = 0;
		double mustdev = 0;
		double psistdev = 0;
		double rhostdev = 0;

		int clockprior = 0;
		double ratemean = 1.0;
		double ratestdev = 1.0;
		double morphoratemean = 1.0;
		double morphoratestdev = 1.0;

		string name = "";
		bool conjpath = true;

		int clamptree = 0;
		int prior = 0;

		int nSegments = 100;
		double fixalpha = 0;

		int every = 1;
		int until = -1;

		int force = 0;

		try	{

			if (argc == 1)	{
				throw(0);
			}

			int i = 1;
			while (i < argc)	{
				string s = argv[i];

				if (s == "-nuc")	{
					i++;
					nucdatafile = argv[i];
				}
				else if ((s == "-t") || (s == "-T"))	{
					i++;
					treefile = argv[i];
				}
				else if (s == "-morpho")	{
					i++;
					morpho2datafile = argv[i];
					i++;
					morpho3datafile = argv[i];
					i++;
					morpho4datafile = argv[i];
				}
				else if (s == "-fixbl")	{
					clamptree = 1;
				}
				else if (s == "-prior")	{
					prior = 1;
				}
				else if (s == "-cal")	{
					i++;
					calibfile = argv[i];
				}
				else if (s == "-rp")	{
					i++;
					rootmean = atof(argv[i]);
					i++;
					rootstdev = atof(argv[i]);
				}
				else if (s == "-serialcoal")	{
					chronoprior = 1;
					i++;
					rootmean = atof(argv[i]);
					i++;
					rootstdev = atof(argv[i]);
				}
				else if (s == "-logisticcoal")	{
					chronoprior = 2;
					i++;
					rootmean = atof(argv[i]);
					i++;
					rootstdev = atof(argv[i]);
					i++;
					T1 = atof(argv[i]);
					i++;
					divrate = atof(argv[i]);
					i++;
					divratestdev = atof(argv[i]);
					i++;
					extrate = atof(argv[i]);
					i++;
					extratestdev = atof(argv[i]);
					i++;
					K0 = atof(argv[i]);
					i++;
					K0stdev = atof(argv[i]);
					i++;
					massext = atof(argv[i]);
					i++;
					massextstdev = atof(argv[i]);
					i++;
					K1 = atof(argv[i]);
					i++;
					K1stdev = atof(argv[i]);
				}
				else if (s == "-serialBD")	{
					chronoprior = 3;
					i++;
					rootmean = atof(argv[i]);
					i++;
					rootstdev = atof(argv[i]);
					i++;
					bddivratemean = atof(argv[i]);
					i++;
					bddivratestdev = atof(argv[i]);
					i++;
					mumean = atof(argv[i]);
					i++;
					mustdev = atof(argv[i]);
					i++;
					psimean = atof(argv[i]);
					i++;
					psistdev = atof(argv[i]);
					i++;
					rhomean = atof(argv[i]);
					i++;
					rhostdev = atof(argv[i]);
				}
				else if (s == "-logisticBD")	{
					chronoprior = 4;
					i++;
					rootmean = atof(argv[i]);
					i++;
					rootstdev = atof(argv[i]);
					i++;
					T1 = atof(argv[i]);
					i++;
					divrate = atof(argv[i]);
					i++;
					divratestdev = atof(argv[i]);
					i++;
					extrate = atof(argv[i]);
					i++;
					extratestdev = atof(argv[i]);
					i++;
					K0 = atof(argv[i]);
					i++;
					K0stdev = atof(argv[i]);
					i++;
					massext = atof(argv[i]);
					i++;
					massextstdev = atof(argv[i]);
					i++;
					K1 = atof(argv[i]);
					i++;
					K1stdev = atof(argv[i]);
					i++;
					psimean = atof(argv[i]);
					i++;
					psistdev = atof(argv[i]);
				}
				else if (s == "-divsampling")	{
					i++;
					divcutoff = atof(argv[i]);
					i++;
					Nextant = atof(argv[i]);
				}
				else if (s == "-nodebd")	{
					chronoprior = 4;
					i++;
					rootmean = atof(argv[i]);
					i++;
					rootstdev = atof(argv[i]);
					i++;
					bddivratemean = atof(argv[i]);
					i++;
					bddivratestdev = atof(argv[i]);
					i++;
					mumean = atof(argv[i]);
					i++;
					mustdev = atof(argv[i]);
				}
				else if (s == "-cauchy")	{
					if (chronoprior != 4)	{
						throw;
					}
					chronoprior = 5;
				}
				else if (s == "-rootbounded")	{
					if (chronoprior != 4)	{
						throw;
					}
					chronoprior = 6;
				}
				else if (s == "-ln")	{
					clockprior = 0;
				}
				else if (s == "-wn")	{
					clockprior = 1;
					i++;
					ratemean = atof(argv[i]);
					i++;
					ratestdev = atof(argv[i]);
					i++;
					morphoratemean = atof(argv[i]);
					i++;
					morphoratestdev = atof(argv[i]);
				}
				else if (s == "-lnwn")	{
					clockprior = 2;
				}
				else if (s == "-asjump")	{
					clockprior = 3;
				}
				else if (s == "-as")	{
					clockprior = 4;
				}
				else if (s == "-aswn")	{
					clockprior = 5;
				}
				else if (s == "-brjump")	{
					clockprior = 6;
				}
				else if (s == "-alpha")	{
					i++;
					fixalpha = atof(argv[i]);
				}
				else if (s == "-nseg")	{
					i++;
					nSegments = atoi(argv[i]);
				}
				else if (s == "-conjpath")	{
					conjpath = 1;
				}
				else if (s == "-nonconjpath")	{
					conjpath = 0;
				}
				else if (s == "-f")	{
					force = 1;
				}
				else if (s == "-x")	{
					i++;
					if (i == argc) throw(0);
					every = atoi(argv[i]);
					i++;
					if (i == argc) throw(0);
					until = atoi(argv[i]);
				}
				else	{
					if (i != (argc -1))	{
						throw(0);
					}
					name = argv[i];
				}
				i++;
			}
			if ((nucdatafile == "") || (treefile == "") || (name == ""))	{
				throw(0);
			}
		}
		catch(...)	{
			cerr << "logn -d <datafile> -morpho <m2 m3 m4> -t <tree> [-x <every> <until>] <chainname> \n";
			cerr << '\n';
			exit(1);
		}

		LogNChain* chain = new LogNChain(treefile,nucdatafile,morpho2datafile,morpho3datafile,morpho4datafile,calibfile,rootmean,rootstdev,divrate,extrate,massext,K0,K1,divratestdev,extratestdev,massextstdev,K0stdev,K1stdev,T1,bddivratemean,mumean,psimean,rhomean,bddivratestdev,mustdev,psistdev,rhostdev,divcutoff,Nextant,chronoprior,clockprior,ratemean,ratestdev,morphoratemean,morphoratestdev,conjpath,clamptree,prior,nSegments,fixalpha,name,force);
		chain->SetEvery(every);
		chain->SetUntil(until);
		cerr << "start\n";
		cerr << '\n';
		chain->Start();
		cerr << '\n';
		cerr << "exit\n";
		cerr << chain->GetSize() << " points saved, current ln prob = " << chain->GetModel()->GetLogProb() << "\n";
		cerr << '\n';
	}
}

