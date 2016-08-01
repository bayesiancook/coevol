
#include "MeanValTree.h"

#include "BaseType.h"
#include "RandomTypes.h"
#include "ProbModel.h"
#include "IID.h"
#include "PrecisionNormalTreeProcess.h"
#include "AutoRegressiveMultiVariateTreeProcess.h"
#include "ConjugateMultiVariateTreeProcess.h"
#include "BranchProcess.h"
#include "ContinuousData.h"
#include "AncestralData.h"
#include "MeanExpTree.h"
#include "Normal.h"

#include "Jeffreys.h"
#include "Partition.h"
#include "MultiVarNormal.h"

#include "FixedLengthTree.h"


class Simu	{

	Tree* tree;
	int Ncont;
	int Nanc;
	ContinuousData* contdata;
	AncestralData* ancdata;
	AncestralData* contancdata;
	const TaxonSet* taxset;

	Const<Real>* Zero;
	Const<PosReal>* One;

	FixedLengthTree* fixtree;

	JeffreysIIDArray* DiagArray;
	SigmaZero* sigmaZero;
	Rvar<CovMatrix>* sigma;

	MultiVariateTreeProcess* process;

	public:

	Simu(string treefile, int inNcont, int inNanc, string ancdatafile, double incorrel, double incov1, double incov2, string name, int df, int nrep)	{

		cerr << "tree\n";
		tree = new Tree(treefile);
		Ncont = inNcont;
		Nanc = inNanc;

		cerr << "anc\n";
		ancdata = new FileAncestralData(tree,ancdatafile,Nanc);
		contancdata = new FileAncestralData(tree,ancdatafile,Ncont);
		taxset = ancdata->GetTaxonSet();

		cerr << "cont\n";
		contdata = new ContinuousData(taxset,Ncont);

		cerr << "model\n";
		Zero = new Const<Real>(0);
		One = new Const<PosReal>(1);

		fixtree = new FixedLengthTree(tree,One);

		double mindiag = 0.001;
		double maxdiag = 1000;
		DiagArray = new JeffreysIIDArray(Ncont+Nanc,mindiag,maxdiag,Zero);
		DiagArray->ClampAt(1.0);
		sigmaZero = new SigmaZero(DiagArray);
		sigma = new InverseWishartMatrix(sigmaZero, Ncont+Nanc+df);
		if (incorrel != -1)	{
			if (sigma->GetDim() != 2)	{
				cerr << "error when setting correlation coef\n";
				exit(1);
			}
			double cov = incorrel * sqrt(incov1 * incov2);
			(*sigma)[0][0] = incov1;
			(*sigma)[1][1] = incov2;
			(*sigma)[0][1] = cov;
			(*sigma)[1][0] = cov;
			sigma->Corrupt(true);
		}

		// sigma->SetIdentity();

		ofstream sos((name + ".sigma").c_str());
		sos << *sigma << '\n';

		process = new MultiVariateTreeProcess(sigma,fixtree,0,0);

		if (nrep == 1)	{
			process->GetLeafData(contdata,0);
			process->GetAncData(ancdata,Ncont);
			process->GetAncData(contancdata,0);

			ofstream cos((name + ".cont").c_str());
			contdata->ToStream(cos);

			ofstream aos((name + ".anc").c_str());
			ancdata->ToStream(aos);

			ofstream caos((name + ".fullcont").c_str());
			contancdata->ToStream(caos);

			double anc1 = (*process->GetNodeVal(process->GetRoot()->GetNode()))[0];
			double anc2 = (*process->GetNodeVal(process->GetRoot()->Next()->Out()->GetNode()))[0];
			double anc3 = (*process->GetNodeVal(process->GetRoot()->Next()->Next()->Out()->GetNode()))[0];
			ofstream sos((name + ".deltaroot").c_str());
			sos << anc1 << '\t' << anc2 << '\t' << anc3 << '\n';
		}
		else	{
			for (int rep=0; rep<nrep; rep++)	{

				cerr << '.';
				process->Sample();

				ostringstream s;
				s << name << rep;
				process->GetLeafData(contdata,0);
				process->GetAncData(ancdata,Ncont);
				process->GetAncData(contancdata,0);

				ofstream cos((s.str() + ".cont").c_str());
				contdata->ToStream(cos);

				ofstream aos((s.str() + ".anc").c_str());
				ancdata->ToStream(aos);

				ofstream caos((s.str() + ".fullcont").c_str());
				contancdata->ToStream(caos);

				double anc1 = (*process->GetNodeVal(process->GetRoot()->GetNode()))[0];
				double anc2 = (*process->GetNodeVal(process->GetRoot()->Next()->Out()->GetNode()))[0];
				double anc3 = (*process->GetNodeVal(process->GetRoot()->Next()->Next()->Out()->GetNode()))[0];
				ofstream sos((s.str() + ".deltaroot").c_str());
				sos << anc1 << '\t' << anc2 << '\t' << anc3 << '\n';

			}
			cerr << '\n';
		}
	}

	~Simu()	{
	}
};

int main(int argc, char* argv[])	{

	string treefile = "";
	string ancdatafile = "None";
	string name = "none";
	int Ncont = 1;
	int Nanc = 1;
	int df = 0;
	int nrep = 1;
	double correl = -1;
	double cov1 = -1;
	double cov2 = -1;

	try	{

		if (argc == 1)	{
			throw(0);
		}

		int i = 1;
		while (i < argc)	{
			string s = argv[i];

			if (s == "-anc")	{
				i++;
				ancdatafile = argv[i];
			}
			else if ((s == "-t") || (s == "-T"))	{
				i++;
				treefile = argv[i];
			}
			else if (s == "-ncont")	{
				i++;
				Ncont = atoi(argv[i]);
			}
			else if (s == "-nanc")	{
				i++;
				Nanc = atoi(argv[i]);
			}
			else if (s == "-correl")	{
				i++;
				cov1 = atof(argv[i]);
				i++;
				cov2 = atof(argv[i]);
				i++;
				correl = atof(argv[i]);
			}
			else if (s == "-df")	{
				i++;
				df = atoi(argv[i]);
			}
			else if (s == "-nrep")	{
				i++;
				nrep = atoi(argv[i]);
			}
			else	{
				if (i != (argc -1))	{
					throw(0);
				}
				name = argv[i];
			}
			i++;
		}
		if ((ancdatafile == "") || (treefile == "") || (name == ""))	{
			throw(0);
		}
	}
	catch(...)	{
		cerr << "error in command\n";
		cerr << '\n';
		exit(1);
	}

	Simu* simu = new Simu(treefile,Ncont,Nanc,ancdatafile,correl,cov1,cov2,name,df,nrep);
	cerr << simu << '\n';
}


