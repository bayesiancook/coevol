
#include "Chain.h"
#include "BranchOmegaMultivariateModel.h"
#include "ConjugateBranchOmegaMultivariateModel.h"


int main(int argc, char* argv[])	{

	string datafile = "";
	string treefile = "";
	string contdatafile = "None";

	string calibfile = "None";
	double rootage = 0;
	double rootstdev = 0;
	int chronoprior = 0;
	double meanchi = 1;
	double meanchi2 = 1;

	double priorsigma = 0.1;
	int df = 2;
	int mutmodel = 0;
	int gc = 0;

	string name = "";

	bool clampdiag = false;
	bool autoregressive = false;

	GeneticCodeType type = Universal;

	int conjpath = 0;
	int contdatatype = 0;
	int omegaratiotree = 1;

	bool clamproot = true;
	bool meanexp = false;
	bool normalise = true;

	string bounds = "None";

	int from = 0;
	bool postpred = false;

	int nrep;

	string nucmatrixfile = "None";
	try	{

		if (argc == 1)	{
			throw(0);
		}

		int i = 1;
		while (i < argc)	{
			string s = argv[i];

			if (s == "-d")	{
				i++;
				datafile = argv[i];
			}
			else if (s == "-ppred")	{
				postpred = true;
			}
			else if ((s == "-t") || (s == "-T"))	{
				i++;
				treefile = argv[i];
			}
			else if (s == "-from")	{
				i++;
				from = atoi(argv[i]);
			}
			else if (s == "-nuc")	{
				i++;
				nucmatrixfile = argv[i];
			}
			else if (s == "-nrep")	{
				i++;
				nrep = atoi(argv[i]);
			}
			else if (s == "-c")	{
				i++;
				contdatafile = argv[i];
			}
			else if (s == "-gc")	{
				gc = true;
			}
			else if (s == "-oup")	{
				autoregressive = true;
			}
			else if (s == "-bp")	{
				autoregressive = false;
			}
			else if (s == "-omegaratiotree")	{
				omegaratiotree = true;
			}
			else if (s == "-meanexp")	{
				meanexp = true;
			}
			else if (s == "-priornu")	{
				i++;
				priorsigma = atof(argv[i]);
			}
			else if (s == "-diag")	{
				clampdiag = true;
			}
			else if (s == "-mtmam")	{
				type = MtMam;
			}
			else if (s == "-mtinv")	{
				type = MtInv;
			}
			else if (s == "-uni")	{
				type = Universal;
			}
			else if (s == "-cal")	{
				i++;
				calibfile = argv[i];
				i++;
				rootage = atof(argv[i]);
				i++;
				rootstdev = atof(argv[i]);
			}
			else	{
				if (i != (argc -1))	{
					throw(0);
				}
				name = argv[i];
			}
			i++;
		}
		if ((datafile == "") || (treefile == "") || (name == ""))	{
			throw(0);
		}
	}
	catch(...)	{
		cerr << "mulomega -d <alignment> -t <tree> [-c <phenodata>] [-conj -nonconj -bp -oup -diag] [-uni -mtmam -mtinv] [-x <every> <until>] <chainname> \n";
		cerr << '\n';
		exit(1);
	}

	BranchOmegaMultivariateModel* model = new BranchOmegaMultivariateModel(datafile,treefile,contdatafile,calibfile,rootage,rootstdev,chronoprior,meanchi,meanchi2,priorsigma,df,mutmodel,gc,clampdiag,autoregressive,conjpath,contdatatype,omegaratiotree,clamproot,meanexp,normalise,0 /* nrep */ ,bounds,false,type);

	/*
	if (postpred)	{
		model->PostPredSimu(nrep,name);
		exit(1);
	}
	*/

	cerr << "read nucmatrix\n";
	ifstream nis(nucmatrixfile.c_str());
	model->GetNucMatrix(nis);

	int rep = from;
	while (rep < nrep)	{
	// for (int rep=from; rep<nrep; rep++)	{
		ostringstream s;
		s << name << rep;
		cerr << s.str() << '\n';
		model->Simulate(s.str());
		ofstream os((s.str() + ".").c_str());
		cerr << "COVARIANCE : \n";
		cerr << *(model->GetCovMatrix()) << '\n';

		ofstream cos((s.str() + ".cov").c_str());
		cos << *(model->GetCovMatrix()) << '\n';
		double c1 = model->GetMaxdS();
		double c2 = 0;
		cerr << "MAX dS : " << c1 << '\n';
		if (omegaratiotree == 2)	{
			c2 = model->GetMaxdN();
			cerr << "MAX dN : " << c2 << '\n';
		}
		else	{
			c2 = model->GetMaxOmega();
			cerr << "Max omega : " << c2 << '\n';
		}
		if (autoregressive)	{
			cerr << "MEAN : \n";
			cerr << *(model->GetMeanVector()) << '\n';
		}
		if ((c1 < 5) && (c2 < 5))	{
			rep++;
		}
	}

	/*
	ifstream is(treefile.c_str());
	for (int rep=0; rep<from; rep++)	{
		string tree;
		is >> tree;
	}

	for (int rep=from; rep<nrep; rep++)	{
		ostringstream s;
		s << name << rep;
		cerr << s.str() << '\n';

		string tree;
		is >> tree;
		ofstream tos((s.str() + ".tree").c_str());
		tos << tree << '\n';
		tos.close();

		double priorsigma = 1;

		BranchOmegaMultivariateModel* model = new BranchOmegaMultivariateModel(datafile,s.str() + ".tree",contdatafile,calibfile,rootage,rootstdev,priorsigma,gc,clampdiag,autoregressive,false,true,type);
		model->Simulate(s.str());
		ofstream os((s.str() + ".").c_str());
		cerr << "COVARIANCE : \n";
		cerr << *(model->GetCovMatrix()) << '\n';
		if (autoregressive)	{
			cerr << "MEAN : \n";
			cerr << *(model->GetMeanVector()) << '\n';
		}
	}
	*/
}

