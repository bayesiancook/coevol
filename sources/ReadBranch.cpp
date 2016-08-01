
#include "Sample.h"
#include "ConjugateBranchModel.h"
// #include "BranchModel.h"
#include "MeanValTree.h"
#include "MeanICTree.h"
#include "MeanCovMatrix.h"
#include "StringStreamUtils.h"
#include "MeanChronogram.h"

class BranchSample : public Sample	{

	private:
	string modeltype;
	string datafile;
	string treefile;
	string calibfile;

	string suffstatfile;

	double rootage;
	double rootstdev;

	int chronoprior;
	double meanchi;
	double meanchi2;

	int gc;
	int codonmodel;
	bool clamptree;
	GeneticCodeType type;

	int conjpath;
	int fullconj;
	double mappingfreq;

	bool normalise;
	int nrep;


	public:

	string GetModelType() {return modeltype;}

	BranchModel* GetModel() {return (BranchModel*) model;}

	BranchSample(string filename, int inburnin, int inevery, int inuntil) : Sample(filename,inburnin,inevery,inuntil)	{
		Open();
	}

	void Open()	{

		codonmodel = 0;
		gc = 0;
		clamptree = false;

		// open <name>.param
		ifstream is((name + ".param").c_str());

		// check that file exists
		if (!is)	{
			cerr << "error : cannot find file : " << name << ".param\n";
			exit(1);
		}

		// read model type, and other standard fields
		is >> modeltype;
		is >> type;
		is >> datafile >> treefile;
		is >> calibfile >> rootage >> rootstdev;
		is >> chronoprior >> meanchi >> meanchi2;
		is >> gc >> codonmodel;
		is >> conjpath;
		is >> normalise >> nrep;
		is >> clamptree;
		int check;
		is >> check;
		if (check)	{
			is >> suffstatfile;
			is >> check;
			if (check)	{
				is >> mappingfreq;
				is >> check;
				if (check)	{
					is >> fullconj;
					is >> check;
					if (check)	{
						cerr << "error when reading model\n";
						exit(1);
					}
				}
			}
		}

		is >> chainevery >> chainuntil >> chainsize;

		if (modeltype == "BRANCHWISE")	{
			model = new BranchModel(datafile,treefile,calibfile,rootage,rootstdev,chronoprior,meanchi,meanchi2,gc,codonmodel,conjpath,fullconj,mappingfreq,clamptree,normalise,nrep,suffstatfile,false,type);
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

	void ComputeMeanOMGC()	{

		double meanomts = 0;
		double meanomtv0 = 0;
		double meanomtvgc = 0;
		double meangc = 0;
		double varomts = 0;
		double varomtv0 = 0;
		double varomtvgc = 0;
		double vargc = 0;

		for (int i=0; i<size; i++)	{

			cerr << '.';
			GetNextPoint();

			double tmp = GetModel()->GetMeanOmegaTs();
			meanomts += tmp;
			varomts += tmp * tmp;
			tmp = GetModel()->GetMeanOmegaTv0();
			meanomtv0 += tmp;
			varomtv0 += tmp * tmp;
			tmp = GetModel()->GetMeanOmegaTvGC();
			meanomtvgc += tmp;
			varomtvgc += tmp * tmp;
			tmp = GetModel()->GetMeanGC();
			meangc += tmp;
			vargc += tmp * tmp;
		}

		cerr << '\n';

		meanomts /= size;
		varomts /= size;
		varomts -= meanomts * meanomts;
		meanomtv0 /= size;
		varomtv0 /= size;
		varomtv0 -= meanomtv0 * meanomtv0;
		meanomtvgc /= size;
		varomtvgc /= size;
		varomtvgc -= meanomtvgc * meanomtvgc;
		meangc /= size;
		vargc /= size;
		vargc -= meangc * meangc;

		ofstream os((GetName() + ".means").c_str());
		os << meanomts << '\t' << sqrt(varomts) << '\t';
		os << meanomtv0 << '\t' << sqrt(varomtv0) << '\t';
		os << meanomtvgc << '\t' << sqrt(varomtvgc) << '\t';
		os << meangc << '\t' << sqrt(vargc) << '\n';
	}

	void WriteSynNonSyn()	{

		GetNextPoint();
		cerr << "update\n";
		GetModel()->Update();
		cerr << "update ok\n";
		GetModel()->PrintSynNonSyn(name);
	}

	void Read(bool printlog, bool printmean, bool printci, bool printstdev, bool withleaf, bool withinternal)	{

		cerr << "read\n";

		MeanBranchTree* meantree = new MeanBranchTree(GetModel()->GetTree(),false);
		MeanChronogram* meanchrono = new MeanChronogram(GetModel()->GetTree());

		double meanom1 =0;
		double varom1 = 0;
		double meanom2 =0;
		double varom2 = 0;
		double meanom =0;
		double varom = 0;

		MeanExpNormTree* meansynrate = new MeanExpNormTree(GetModel()->GetTree(),false,printlog,printmean,printci,printstdev,withleaf,withinternal);

		MeanExpNormTree* meanomega = 0;
		MeanExpNormTree* meanomegats = 0;
		MeanExpNormTree* meanomegatv0 = 0;
		MeanExpNormTree* meanomegatvgc = 0;
		if (GetModel()->Codon3())	{
			meanomegats = new MeanExpNormTree(GetModel()->GetTree(),false,printlog,printmean,printci,printstdev,withleaf,withinternal);
			meanomegatv0 = new MeanExpNormTree(GetModel()->GetTree(),false,printlog,printmean,printci,printstdev,withleaf,withinternal);
			meanomegatvgc = new MeanExpNormTree(GetModel()->GetTree(),false,printlog,printmean,printci,printstdev,withleaf,withinternal);
		}
		else if (GetModel()->Codon())	{
			meanomega = new MeanExpNormTree(GetModel()->GetTree(),false,printlog,printmean,printci,printstdev,withleaf,withinternal);
			meanomega->ActivatePP(1.0);
		}

		MeanExpNormTree* meangc = 0;
		MeanExpNormTree* meangc1 = 0;
		MeanExpNormTree* meangc2 = 0;
		MeanExpNormTree* meangc3 = 0;
		if (GetModel()->isGCActivated())	{
			meangc = new MeanExpNormTree(GetModel()->GetTree(),true,printlog,printmean,printci,printstdev,withleaf,withinternal);
		}
		if (GetModel()->isGC3Activated())	{
			meangc1 = new MeanExpNormTree(GetModel()->GetTree(),true,printlog,printmean,printci,printstdev,withleaf,withinternal);
			meangc2 = new MeanExpNormTree(GetModel()->GetTree(),true,printlog,printmean,printci,printstdev,withleaf,withinternal);
			meangc3 = new MeanExpNormTree(GetModel()->GetTree(),true,printlog,printmean,printci,printstdev,withleaf,withinternal);
		}

		// cycle over the sample
		int halfsize = size / 2;
		for (int i=0; i<size; i++)	{
			cerr << '.';

			GetNextPoint();

			if (!GetModel()->Unconstrained())	{
				GetModel()->GetChronogram()->specialUpdate();
				meanchrono->Add(GetModel()->GetChronogram());
				// meantree->Add(GetModel()->GetSynRateTree());
				meansynrate->Add(GetModel()->GetSynRateTree(), GetModel()->GetLengthTree(),true);
				// meansynrate->Add((NodeVarTree<Real>* ) GetModel()->GetSynRateTree(), GetModel()->GetLengthTree());
			}
			else	{
				meantree->Add(GetModel()->GetSynGammaTree());
			}


			if (GetModel()->isGCActivated())	{
				meangc->Add(GetModel()->GetGCTree(), GetModel()->GetLengthTree());
			}
			if (GetModel()->isGC3Activated())	{
				meangc1->Add(GetModel()->GetGCTree1(), GetModel()->GetLengthTree());
				meangc2->Add(GetModel()->GetGCTree2(), GetModel()->GetLengthTree());
				meangc3->Add(GetModel()->GetGCTree3(), GetModel()->GetLengthTree());
			}
			if (GetModel()->Codon3())	{
				meanomegats->Add(GetModel()->GetOmegaTsTree(), GetModel()->GetLengthTree());
				meanomegatv0->Add(GetModel()->GetOmegaTv0Tree(), GetModel()->GetLengthTree());
				meanomegatvgc->Add(GetModel()->GetOmegaTvGCTree(), GetModel()->GetLengthTree());
			}
			else if (GetModel()->Codon())	{
				meanomega->Add(GetModel()->GetOmegaTree(), GetModel()->GetLengthTree());
				double tmp = GetModel()->GetMeanOmega();
				meanom += tmp;
				varom += tmp * tmp;
				if (i<halfsize)	{
					meanom1 += tmp;
					varom1 += tmp * tmp;
				}
				else	{
					meanom2 += tmp;
					varom2 += tmp * tmp;
				}
			}
		}

		cerr << '\n';
		cerr << "normalise\n";

		meantree->Normalise();
		ofstream tos((GetName() + ".postmean.tre").c_str());
		meantree->ToStream(tos);

		if (!GetModel()->Unconstrained())	{
			meanchrono->Normalise();
			ofstream chos((GetName() + ".postmeandates.tre").c_str());
			meanchrono->ToStream(chos);

			ofstream cchos((GetName() + ".postmeandates.tab").c_str());
			meanchrono->Tabulate(cchos);

			meansynrate->Normalise();
			ofstream sos((GetName() + ".postmeansynrate.tre").c_str());
			meansynrate->ToStream(sos);
			ofstream ssos((GetName() + ".postmeansynrate.tab").c_str());
			meansynrate->Tabulate(ssos);
			cerr << "reconstructed variation in Ks in " << name << ".postmeansynrate.tre\n";

		}

		if (GetModel()->Codon3())	{
			meanomegats->Normalise();
			ofstream tsoos((GetName() + ".postmeanomegats.tre").c_str());
			meanomegats->ToStream(tsoos);
			meanomegatv0->Normalise();
			ofstream tv0oos((GetName() + ".postmeanomegatv0.tre").c_str());
			meanomegatv0->ToStream(tv0oos);
			meanomegatvgc->Normalise();
			ofstream tvgcoos((GetName() + ".postmeanomegatvgc.tre").c_str());
			meanomegatvgc->ToStream(tvgcoos);
		}
		else if (GetModel()->Codon())	{
			meanomega->Normalise();
			ofstream oos((GetName() + ".postmeanomega.tre").c_str());
			meanomega->ToStream(oos);
			cerr << "reconstructed variation in omega in " << name << ".postmeanomega.tre\n";
			meanom /= size;
			varom /= size;
			varom -= meanom * meanom;
			meanom1 /= halfsize;
			varom1 /= halfsize;
			varom1 -= meanom1 * meanom1;
			meanom2 /= (size - halfsize);
			varom2 /= (size - halfsize);
			varom2 -= meanom2 * meanom2;
			double z = fabs(meanom2 - meanom1) / (sqrt(varom1) + sqrt(varom2));

			ofstream moos((GetName() + ".postmeanomega").c_str());
			moos << z << '\t' << meanom1 << '\t' << meanom2 << '\t' << sqrt(varom1) << '\t' << sqrt(varom2) << '\t' << meanom << '\t' << sqrt(varom) << '\n';
		}

		if (GetModel()->isGCActivated())	{
			meangc->Normalise();
			ofstream gcos((GetName() + ".postmeangc.tre").c_str());
			meangc->ToStream(gcos);
			cerr << "reconstructed variation in gc in " << name << ".postmeangc.tre\n";
			meangc->SetWithLeaf(true);
			meangc->SetWithInternal(true);
			ofstream ggcos((GetName() + ".postmeangc.tab").c_str());
			meangc->Tabulate(ggcos);
			ggcos.close();
		}
		if (GetModel()->isGC3Activated())	{
			meangc1->Normalise();
			meangc2->Normalise();
			meangc3->Normalise();
			ofstream gc1os((GetName() + ".postmeangc1.tre").c_str());
			ofstream gc2os((GetName() + ".postmeangc2.tre").c_str());
			ofstream gc3os((GetName() + ".postmeangc3.tre").c_str());
			meangc1->ToStream(gc1os);
			meangc2->ToStream(gc2os);
			meangc3->ToStream(gc3os);
			cerr << "reconstructed variation in gc1 in " << name << ".postmeangc1.tre\n";
			cerr << "reconstructed variation in gc2 in " << name << ".postmeangc2.tre\n";
			cerr << "reconstructed variation in gc3 in " << name << ".postmeangc3.tre\n";
		}

		if (GetModel()->Codon3())	{
			meanomegats->SetWithLeaf(true);
			meanomegats->SetWithInternal(true);
			ofstream tsos((GetName() + ".postmeanomegats.tab").c_str());
			meanomegats->Tabulate(tsos);
			tsos.close();
			meanomegatv0->SetWithLeaf(true);
			meanomegatv0->SetWithInternal(true);
			ofstream tv0os((GetName() + ".postmeanomegatv0.tab").c_str());
			meanomegatv0->Tabulate(tv0os);
			tv0os.close();
			meanomegatvgc->SetWithLeaf(true);
			meanomegatvgc->SetWithInternal(true);
			ofstream tvgcos((GetName() + ".postmeanomegatvgc.tab").c_str());
			meanomegatvgc->Tabulate(tvgcos);
			tvgcos.close();
		}
		else if (GetModel()->Codon())	{
			meanomega->SetWithLeaf(true);
			meanomega->SetWithInternal(true);
			ofstream ooos((GetName() + ".postmeanomega.tab").c_str());
			meanomega->Tabulate(ooos);
			ooos.close();
			ofstream cos((GetName() + ".postmeantime.tab").c_str());
			meanchrono->TabulateTimes(cos);
			cos.close();

			ofstream ppos((GetName() + ".ppomega.tab").c_str());
			meanomega->TabulatePP(ppos);
			ppos.close();
		}

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

	bool means = false;

	// bool nans = false;

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
			else if (s == "-m")	{
				means = true;
			}
			/*
			else if (s == "-nans")	{
				nans = true;
			}
			*/
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

	BranchSample sample(name,burnin,every,until);

	/*
	if (nans)	{
		sample.WriteSynNonSyn();
		exit(1);
	}
	*/

	if (means)	{
		sample.ComputeMeanOMGC();
		exit(1);
	}

	sample.Read(printlog,printmean,printci,printstdev,withleaf,withinternal);

}



