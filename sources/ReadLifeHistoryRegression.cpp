
#include "Sample.h"
#include "LifeHistoryRegressionModel.h"
#include "MeanValTree.h"
#include "MeanICTree.h"
#include "MeanCovMatrix.h"
#include "StringStreamUtils.h"
#include "MeanChronogram.h"
#include "MeanChronoBubbleTree.h"

#include "MeanTimeLine.h"


class LifeHistoryRegressionSample : public Sample	{

	private:
	string modeltype;
	string datafile;
	string treefile;
	string contdatafile;
	string calibfile;
	string suffstatfile;
	string rootfile;

	double rootage;
	double rootstdev;

	int chronoprior;
	double meanchi;
	double meanchi2;

	double softa;

	int mutmodel;
	int gc;
	int ry;
	bool clampdiag;
	int noreg;
	bool autoregressive;
	bool clamproot;
	bool clamptree;
	bool meanexp;
	GeneticCodeType type;

	int conjpath;
	int contdatatype;

	int omegaratiotree;

	bool normalise;
	int nrep;

	int df;

	int krkctype;

	string bounds;
	string mix;

	int nsplit;
	int withdrift;
	int uniformprior;
	bool withtimeline;

	public:

	string GetModelType() {return modeltype;}

	LifeHistoryRegressionModel* GetModel() {return (LifeHistoryRegressionModel*) model;}

	LifeHistoryRegressionSample(string filename, int inburnin, int inevery, int inuntil) : Sample(filename,inburnin,inevery,inuntil)	{
		bounds = "None";
		mix = "None";
		suffstatfile = "None";
		rootfile = "None";
		withtimeline = false;
		autoregressive = false;
		krkctype = 0;
		Open();
	}

	void Open()	{

		mutmodel = 0;

		// open <name>.param
		ifstream is((name + ".param").c_str());

		// check that file exists
		if (!is)	{
			cerr << "error : cannot find file : " << name << ".param\n";
			exit(1);
		}

		double priorsigma = 1;
		ry = 0;

		nsplit = 1;
		withdrift = 0;
		uniformprior = 0;

		// read model type, and other standard fields
		is >> modeltype;
		is >> type;
		is >> datafile >> treefile >> contdatafile;
		is >> calibfile >> rootage >> rootstdev;
		is >> chronoprior >> meanchi >> meanchi2;
		if (chronoprior == 4)	{
			is >> softa;
		}
		is >> clampdiag >> gc;
		is >> conjpath;
		is >> contdatatype;
		is >> omegaratiotree;
		is >> clamproot >> meanexp;
		is >> normalise >> nrep;
		int check;
		is >> check;
		if (check)	{
			is >> mutmodel;
			is >> check;
			if (check)	{
				is >> df;
				is >> check;
				if (check)	{
					is >> bounds;
					is >> check;
					if (check)	{
						is >> priorsigma;
						is >> check;
						if (check)	{
							is >> mix;
							is >> check;
							if (check)	{
								is >> nsplit;
								is >> check;
								if (check)	{
									is >> withdrift;
									is >> check;
									if (check)	{
										is >> clamptree;
										is >> check;
										if (check)	{
											is >> suffstatfile;
											is >> check;
											if (check)	{
												is >> autoregressive;
												is >> check;
												if (check)	{
													is >> rootfile;
													is >> check;
													if (check)	{
														is >> ry;
														is >> check;
														if (check)	{
															is >> noreg;
															is >> check;
															if (check)	{
																is >> uniformprior;
																is >> check;
																if (check)	{
																	is >> krkctype;
																	is >> check;
																	if (check)	{
																		cerr << "error when reading model\n";
																		exit(1);
																	}
																}
															}
														}
													}
												}
											}
										}
									}
								}
							}
						}
					}
				}
			}
		}
		is >> chainevery >> chainuntil >> chainsize;
		// the chain's saving frequency, upper limit and current size
		// not to be confused with the sample's subsampling frequency, upper limit and size

		// make a new model depending on the type obtained from the file
		if (modeltype == "LIFEHISTORYREGRESSION")	{
			model = new LifeHistoryRegressionModel(datafile,treefile,contdatafile,calibfile,rootage,rootstdev,chronoprior,meanchi,meanchi2,softa,priorsigma,df,mutmodel,gc,ry,clampdiag,noreg,autoregressive,conjpath,contdatatype,omegaratiotree,clamproot,clamptree,meanexp,normalise,nrep,bounds,mix,nsplit,withdrift,uniformprior,suffstatfile,rootfile,krkctype,false,type);
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

		/*
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
		*/
	}

	void MakeTimeLine()	{
		int N = 100;
		MeanTimeLine* synrate = new MeanTimeLine(N);
		MeanTimeLine* omega = new MeanTimeLine(N);

		int Ncont = GetModel()->Ncont;
		MeanTimeLine** tree = new MeanTimeLine*[Ncont];
		for (int k=0; k<Ncont; k++)	{
			tree[k] = new MeanTimeLine(N);
		}

		for (int j=0; j<size; j++)	{
			cerr << '.';
			GetNextPoint();
			GetModel()->GetSynRateTree()->specialUpdate();
			GetModel()->GetChronogram()->specialUpdate();
			synrate->Add(GetModel()->GetSynRateTree(),GetModel()->GetChronogram(),true);
			if (GetModel()->omegaratiotree)	{
				GetModel()->GetOmegaTree()->specialUpdate();
				omega->Add(GetModel()->GetOmegaTree(),GetModel()->GetChronogram(),false);
			}
			// GetModel()->GetContProcess()->specialUpdate();
			for (int k=0; k<Ncont; k++)	{
				tree[k]->Add(GetModel()->GetContProcess(),GetModel()->GetChronogram(),k);
			}
		}
		cerr << '\n';
		synrate->Normalize();
		ofstream sos((GetName() + ".synratetimeline").c_str());
		synrate->ToStream(sos);

		if (GetModel()->omegaratiotree)	{
			omega->Normalize();
			ofstream oos((GetName() + ".omegatimeline").c_str());
			omega->ToStream(oos);
		}

		for (int k=0; k<Ncont; k++)	{
			tree[k]->Normalize();
			ostringstream s;
			s << GetName() << "." << k << "timeline";
			ofstream os(s.str().c_str());
			tree[k]->ToStream(os);
		}
	}

	void ReadDrift()	{

		int Ncont = GetModel()->Ncont;
		int Nmat = GetModel()->Nmat;
		double mean[Nmat][Ncont];
		double var[Nmat][Ncont];
		double pp[Nmat][Ncont];
		for (int mat=0; mat<Nmat; mat++)	{
			for (int j=0; j<Ncont; j++)	{
				mean[mat][j] = 0;
				var[mat][j] = 0;
				pp[mat][j] = 0;
			}
		}

		// list<double>* sample = new list<double>[Ncont];

		for (int i=0; i<size; i++)	{
			cerr << '.';
			GetNextPoint();
			for (int mat=0; mat<Nmat; mat++)	{
				for (int j=0; j<Ncont; j++)	{
					double tmp = (*(GetModel()->GetDrift(mat)))[j];
					mean[mat][j] += tmp;
					var[mat][j] += tmp * tmp;
					if (tmp > 0)	{
						pp[j][mat]++;
					}
					// sample[j].push_back(tmp);
				}
			}
		}
		cerr << '\n';

		ofstream os((GetName() + ".postmeandrift").c_str());
		for (int mat=0; mat<Nmat; mat++)	{
			for (int j=0; j<Ncont; j++)	{
				mean[mat][j] /= size;
				var[mat][j] /= size;
				var[mat][j] -= mean[mat][j] * mean[mat][j];
				pp[mat][j] /= size;
				//sample[j].sort(sample[j].begin(), sample[j].end());
				os << mat << '\t' << j << '\t' << mean[mat][j] << '\t' << sqrt(var[mat][j]) << '\t' << pp[mat][j] << '\n';
			}
		}
	}

	void Read(bool printlog, bool printmean, bool printci, bool printstdev, bool withleaf, bool withinternal, string mulreg, bool tex, double xscale, double yscale, double nodescale, double nodepower, double barwidth, int fontsize, bool bubbletext)	{

		double meansyn = 0;
		double varsyn = 0;
		double meantime = 0;
		double vartime = 0;

		cerr << "read\n";
		cerr << "df : " << GetModel()->df << '\n';
		int Ncont = GetModel()->Ncont;
		if (GetModel()->Split())	{
			NewickTree::Simplify();
		}

		MeanBranchTree* meantree = new MeanBranchTree(GetModel()->GetFineGrainedTree(),false);

		MeanChronogram* meanchrono = new MeanChronogram(GetModel()->GetTree());

		MeanExpNormTree* meansynrate = new MeanExpNormTree(GetModel()->GetFineGrainedTree(),false,printlog,printmean,printci,printstdev,withleaf,withinternal);
		MeanExpNormTree* meanbranchsynrate = new MeanExpNormTree(GetModel()->GetFineGrainedTree(),false,printlog,printmean,printci,printstdev,withleaf,withinternal);

		MeanExpNormTree* meanomega = 0;
		if (omegaratiotree)	{
			meanomega = new MeanExpNormTree(GetModel()->GetFineGrainedTree(),false,printlog,printmean,printci,printstdev,withleaf,withinternal);
		}

		MeanExpNormTree* meangc = 0;
		MeanExpNormTree* meangc1 = 0;
		MeanExpNormTree* meangc2 = 0;
		MeanExpNormTree* meangc3 = 0;
		if (GetModel()->isGCActivated())	{
			meangc = new MeanExpNormTree(GetModel()->GetFineGrainedTree(),true,printlog,printmean,printci,printstdev,withleaf,withinternal);
		}
		if (GetModel()->isGC3Activated())	{
			meangc1 = new MeanExpNormTree(GetModel()->GetFineGrainedTree(),true,printlog,printmean,printci,printstdev,withleaf,withinternal);
			meangc2 = new MeanExpNormTree(GetModel()->GetFineGrainedTree(),true,printlog,printmean,printci,printstdev,withleaf,withinternal);
			meangc3 = new MeanExpNormTree(GetModel()->GetFineGrainedTree(),true,printlog,printmean,printci,printstdev,withleaf,withinternal);
		}

		MeanExpNormTree** tree = new MeanExpNormTree*[Ncont];
		for (int k=0; k<Ncont; k++)	{
			tree[k] = new MeanExpNormTree(GetModel()->GetFineGrainedTree(),false,printlog,printmean,printci,printstdev,withleaf,withinternal);
		}

		int dim = GetModel()->GetContMatrix()->GetDim();
		MeanCovMatrix*  mat = new MeanCovMatrix(dim);

		double logmixalpha = 0;
		ofstream los((GetName() + ".logmixalpha").c_str());

		// cycle over the sample
		for (int i=0; i<size; i++)	{
			cerr << '.';

			// get next point -> will be stored into "model", and thus, will be accessible through GetModel()

			GetNextPoint();

			GetModel()->GetSynRateTree()->specialUpdate();
			GetModel()->GetChronogram()->specialUpdate();
			double r = GetModel()->GetMeanSynRate();
			meansyn += r;
			varsyn += r*r;

			double t = GetModel()->GetTotalTime();
			meantime += t;
			vartime += t*t;

			double tmp = log(GetModel()->GetMixAlpha());
			logmixalpha += tmp;
			los << tmp << '\n';


			meanchrono->Add(GetModel()->GetChronogram());

			if (GetModel()->Split())	{
				GetModel()->GetSplitLengthTree()->specialUpdate();
			}

			// GetModel()->GetGCTree()->specialUpdate();
			// cout << GetModel()->GetVarGC() << '\n';
			meantree->Add(GetModel()->GetSynRateTree());

			meanbranchsynrate->Add(GetModel()->GetSynRateTree(), GetModel()->GetLengthTree(),true);
			meansynrate->Add(GetModel()->GetSubProcess(), GetModel()->GetLengthTree(), 0);
			if (GetModel()->isGCActivated())	{
				meangc->Add(GetModel()->GetSubProcess(), GetModel()->GetLengthTree(), GetModel()->GetGCIndex());
			}
			if (GetModel()->isGC3Activated())	{
				meangc1->Add(GetModel()->GetSubProcess(), GetModel()->GetLengthTree(), GetModel()->GetGCIndex());
				meangc2->Add(GetModel()->GetSubProcess(), GetModel()->GetLengthTree(), GetModel()->GetGCIndex()+1);
				meangc3->Add(GetModel()->GetSubProcess(), GetModel()->GetLengthTree(), GetModel()->GetGCIndex()+2);
			}
			if (omegaratiotree)	{
				meanomega->Add(GetModel()->GetSubProcess(), GetModel()->GetLengthTree(), 1);
			}
			for (int k=0; k<Ncont; k++)	{
				tree[k]->Add(GetModel()->GetContProcess(), GetModel()->GetLengthTree(), k);
			}

			CovMatrix& m = *(GetModel()->GetContMatrix());
			mat->Add(&m);
		}
		cerr << '\n';
		cerr << "normalise\n";

		meansyn /= size;
		varsyn /= size;
		varsyn -= meansyn * meansyn;
		cerr << '\n';
		cerr << "mean total syn : " << meansyn << " +/- " << sqrt(varsyn) << '\n';

		meantime /= size;
		vartime /= size;
		vartime -= meantime * meantime;
		cerr << "mean total time: " << meantime << " +/- " << sqrt(vartime) << '\n';
		cerr << "average syn rate per tree depth : " << meansyn / meantime << '\n';
		cerr << '\n';

		logmixalpha /= size;
		if (logmixalpha)	{
			cerr << "mean log mix alpha: " << logmixalpha << '\n';
		}

		mat->Normalize();
		ofstream cout((GetName() + ".cov").c_str());

		if (mulreg != "")	{
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
			GetModel()->PrintEntries(cout, indexarray);

			cout << '\n';
			cerr << dim << '\t' << mat->GetDim() << '\n';
			PartialMeanCovMatrix* partmat = new PartialMeanCovMatrix(mat,indexarray,dim);
			// ReducedMeanCovMatrix* redmat = new ReducedMeanCovMatrix(mat,indexarray,dim);
			cout << '\n';
			// cout << *redmat;
			cout << *partmat;
			delete partmat;
			// delete redmat;
		}
		else	{
			cout << "entries are in the following order:\n";
			GetModel()->PrintEntries(cout);

			cout << '\n';
			cout << *mat;
		}

		cerr << "covariance matrix in " << name << ".cov\n";
		cerr << '\n';

		meantree->Normalise();
		ofstream tos((GetName() + ".postmean.tre").c_str());
		meantree->ToStream(tos);

		meanchrono->Normalise();
		ofstream chos((GetName() + ".postmeandates.tre").c_str());
		meanchrono->ToStream(chos);

		ofstream cchos((GetName() + ".postmeandates.tab").c_str());
		meanchrono->Tabulate(cchos);

		if (omegaratiotree)	{
			meanomega->Normalise();
			ofstream oos((GetName() + ".postmeanomega.tre").c_str());
			meanomega->ToStream(oos);
			cerr << "reconstructed variation in omega in " << name << ".postmeanomega.tre\n";
			cerr << "pp of mean leaf values > root value : " << meanomega->GetPPLeafRoot() << '\n';
		}

		if (GetModel()->isGCActivated())	{
			meangc->Normalise();
			ofstream gcos((GetName() + ".postmeangc.tre").c_str());
			meangc->ToStream(gcos);
			cerr << "reconstructed variation in gc in " << name << ".postmeangc.tre\n";
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

		meansynrate->Normalise();
		ofstream sos((GetName() + ".postmeansynrate.tre").c_str());
		meansynrate->ToStream(sos);
		cerr << "reconstructed variation in Ks in " << name << ".postmeansynrate.tre\n";
		cerr << "pp of mean leaf values > root value : " << meansynrate->GetPPLeafRoot() << '\n';

		meanbranchsynrate->Normalise();
		ofstream bsos((GetName() + ".postmeanbranchsynrate.tre").c_str());
		meanbranchsynrate->ToStream(bsos);
		cerr << "reconstructed mean Ks per branch in " << name << ".postmeanbranchsynrate.tre\n";

		for (int k=0; k<Ncont; k++)	{
			tree[k]->Normalise();
			ostringstream s;
			s << GetName() << ".postmean" << k+1 << ".tre";
			ofstream os(s.str().c_str());
			tree[k]->ToStream(os);
			cerr << "reconstructed variation in continuous character # " << k+1 << " in "  << name << ".postmean" << k+1 << ".tre\n";
			cerr << "pp of mean leaf values > root value : " << tree[k]->GetPPLeafRoot() << '\n';
			if (tex)	{
				MeanChronoBubbleTree* textree = new MeanChronoBubbleTree(meanchrono,tree[k],xscale,yscale,nodescale,nodepower,barwidth,fontsize,bubbletext);
				textree->Draw((s.str() + ".tex").c_str());
			}
		}

		ofstream ssos((GetName() + ".postmeansynrate.tab").c_str());
		meansynrate->Tabulate(ssos);
		ssos.close();

		ofstream bssos((GetName() + ".postmeanbranchsynrate.tab").c_str());
		meanbranchsynrate->Tabulate(bssos);
		bssos.close();

		if (omegaratiotree)	{
			ofstream ooos((GetName() + ".postmeanomega.tab").c_str());
			meanomega->Tabulate(ooos);
			ooos.close();
		}

		for (int k=0; k<Ncont; k++)	{
			ostringstream s;
			s << GetName() << ".postmean" << k+1 << ".tab";
			ofstream os(s.str().c_str());
			tree[k]->Tabulate(os);
			// tree[k]->TabulateDistribution(os);
			// tree[k]->TabulateDistribution(os,false);
			ostringstream s2;
			s2 << GetName() << ".postmean" << k << ".loo";
			ofstream os2(s2.str().c_str());
			tree[k]->LooTabulate(os2, GetModel()->GetContinuousData());
		}

		/*
		ostringstream s;
		s << "paste " << GetName() << ".postmeansynrate.tab " << GetName() << ".postmeanomega.tab ";
		for (int k=0; k<Ncont; k++)	{
			s << GetName() << ".postmean" << k << ".tab ";
		}
		s << " > " << GetName() << ".tab";
		system(s.str().c_str());
		*/

		cerr << '\n';
	}
};


int main(int argc, char* argv[])	{

	int burnin = 0;
	int every = 1;
	int until = -1;
	string name;

	string truefile = "";
	bool ic = false;

	bool drift = false;

	bool printlog = false;
	bool printmean = false;
	bool printci = true;
	bool printstdev = false;
	bool withleaf = true;
	bool withinternal = true;

	string mulreg = "";

	bool timeline = false;

	bool tex = false;
	double nodescale = 5;
	double nodepower = 1;
	double barwidth = 0.04;
	double x = 6;
	double y = 8;
	bool withnodebubbles = true;
	bool withtimeintervals = true;
	int fontsize = 4;

	bool bubbletext = false;

	try	{

		if (argc == 1)	{
			throw(0);
		}

		int i = 1;
		while (i < argc)	{
			string s = argv[i];

			if (s == "-ic")	{
				ic = true;
			}
			else if (s == "-timeline")	{
				timeline = true;
			}
			else if (s == "-drift")	{
				drift = true;
			}
			else if (s == "-tex")	{
				tex = true;
			}
			else if (s == "-ns")	{
				i++;
				nodescale = atof(argv[i]);
			}
			else if (s == "-np")	{
				i++;
				nodepower = atof(argv[i]);
			}
			else if (s == "-bw")	{
				i++;
				barwidth = atof(argv[i]);
			}
			else if (s == "+bt")	{
				bubbletext = true;
			}
			else if (s == "-bt")	{
				bubbletext = false;
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
			else if (s == "+nb")	{
				withnodebubbles = true;
			}
			else if (s == "-nb")	{
				withnodebubbles = false;
			}
			else if (s == "+ti")	{
				withtimeintervals = true;
			}
			else if (s == "-ti")	{
				withtimeintervals = false;
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

	LifeHistoryRegressionSample sample(name,burnin,every,until);

	if (drift)	{
		sample.ReadDrift();
		exit(1);
	}

	if (timeline)	{
		sample.MakeTimeLine();
		exit(1);
	}
	if (ic)	{
		sample.ReadIC();
		exit(1);
	}

	sample.Read(printlog,printmean,printci,printstdev,withleaf,withinternal,mulreg,tex,x,y,nodescale,nodepower,barwidth,fontsize,bubbletext);

}



