


#include "Chain.h"
#include "ConjugateBranchOmegaMultivariateModel.h"
#include "StringStreamUtils.h"

class BranchOmegaMultivariateChain : public Chain	{

	private:
	string modeltype;
	string datafile;
	string treefile;
	string contdatafile;
	string calibfile;

	string suffstatfile;

	double rootage;
	double rootstdev;

	int chronoprior;
	double softa;

	double meanchi;
	double meanchi2;

	bool clampdiag;
	bool autoregressive;
	bool clamproot;
	bool clamptree;
	bool meanexp;
	int withdrift;
	int uniformprior;
	string rootfile;
	bool withtimeline;
	bool separatesyn;
	bool separateomega;

	GeneticCodeType type;

	double priorsigma;
	string priorsigmafile;
	int gc;

	int conjpath;
	double mappingfreq;
	int contdatatype;

	int omegaratiotree;

	bool normalise;
	int nrep;
	int ncycle;

	int nsplit;

	int mutmodel;

	int df;

	int krkctype;
	int jitter;

	string bounds;
	string mix;

	int nprocs;
	int myid;

	Chrono mappingchrono;
	Chrono paramchrono;
	Chrono chrono;

	public:

	string GetModelType() {return modeltype;}

	BranchOmegaMultivariateModel* GetModel() {return (BranchOmegaMultivariateModel*) model;}

	BranchOmegaMultivariateChain(string indata, string intree, string incontdata, string incalibfile, double inrootage, double inrootstdev, int inchronoprior, double insofta, double inmeanchi, double inmeanchi2, double inpriorsigma, string inpriorsigmafile, int indf, int inmutmodel, int ingc, string filename, bool inclampdiag, bool inautoregressive, GeneticCodeType intype, bool conjugate, int inconjpath, double inmappingfreq, int incontdatatype, int inomegaratiotree, bool inclamproot, bool inclamptree, bool inmeanexp, bool innormalise, int innrep, int inncycle, string inbounds, string inmix, int innsplit, int inwithdrift, int inuniformprior, string inrootfile, string insuffstatfile, bool inwithtimeline, bool inseparatesyn, bool inseparateomega, int inkrkctype, int injitter, int inmyid, int innprocs, int force)	{
		if (conjugate)	{
			modeltype = "CONJUGATEBRANCHOMEGAMULTIVARIATE";
		}
		else	{
			modeltype = "BRANCHOMEGAMULTIVARIATE";
		}

		myid = inmyid;
		nprocs = innprocs;

		type = intype;
		datafile = indata;
		treefile = intree;
		contdatafile = incontdata;
		calibfile = incalibfile;
		suffstatfile = insuffstatfile;
		rootage = inrootage;
		rootstdev = inrootstdev;
		chronoprior = inchronoprior;
		softa = insofta;
		meanchi = inmeanchi;
		meanchi2 = inmeanchi2;
		name = filename;
		clampdiag = inclampdiag;
		autoregressive = inautoregressive;
		priorsigma = inpriorsigma;
		priorsigmafile = inpriorsigmafile;
		df = indf;
		gc = ingc;
		conjpath = inconjpath;
		mappingfreq = inmappingfreq;
		contdatatype = incontdatatype;
		omegaratiotree = inomegaratiotree;
		clamproot = inclamproot;
		clamptree = inclamptree;
		meanexp = inmeanexp;
		normalise = innormalise;
		nrep = innrep;
		ncycle = inncycle;
		bounds = inbounds;
		mix = inmix;
		mutmodel = inmutmodel;
		nsplit = innsplit;
		withdrift = inwithdrift;
		uniformprior = inuniformprior;
		rootfile = inrootfile;
		withtimeline = inwithtimeline;
		separatesyn = inseparatesyn;
		separateomega = inseparateomega;
		krkctype = inkrkctype;
		jitter = injitter;
		New(force);
	}

	BranchOmegaMultivariateChain(string filename, int inmyid, int innprocs)	{

		myid = inmyid;
		nprocs = innprocs;

		conjpath = 1;
		mappingfreq = 1;
		contdatatype = 0;
		name = filename;
		priorsigma = 1;
		priorsigmafile = "None";
		df = 2;
		omegaratiotree = 0;
		clamproot = false;
		clamptree = false;
		meanexp = false;
		mutmodel = 0;
		nsplit = 1;
		bounds = "None";
		mix = "None";
		cerr << "open\n";
		withdrift = 0;
		uniformprior = 0;
		rootfile = "None";
        softa = 0;
		withtimeline = false;
		separatesyn = false;
		separateomega = false;
		krkctype = 0;
		jitter = 0;
		suffstatfile = "None";
		nrep = 0;
		ncycle = 1;
		Open();
	}

	void New(int force)	{

		if (modeltype == "BRANCHOMEGAMULTIVARIATE")	{
			model = new BranchOmegaMultivariateModel(datafile,treefile,contdatafile,calibfile,rootage,rootstdev,chronoprior,softa,meanchi,meanchi2,priorsigma,priorsigmafile,df,mutmodel,gc,clampdiag,autoregressive,conjpath,mappingfreq,contdatatype,omegaratiotree,clamproot,clamptree,meanexp,normalise,nrep,ncycle,bounds,mix,nsplit,withdrift,uniformprior,rootfile,suffstatfile,withtimeline,separatesyn,separateomega,krkctype,jitter,myid,nprocs,1,type);
		}
		else if (modeltype == "CONJUGATEBRANCHOMEGAMULTIVARIATE")	{
			model = new ConjugateBranchOmegaMultivariateModel(datafile,treefile,contdatafile,calibfile,rootage,rootstdev,chronoprior,softa,meanchi,meanchi2,priorsigma,priorsigmafile,df,mutmodel,gc,autoregressive,conjpath,mappingfreq,contdatatype,omegaratiotree,clamproot,clamptree,meanexp,normalise,nrep,ncycle,bounds,mix,nsplit,withdrift,uniformprior,rootfile,suffstatfile,withtimeline,separatesyn,separateomega,krkctype,jitter,myid,nprocs,1,type);
		}
		else	{
			cerr << "error, does not recognise model type : " << modeltype << '\n';
			exit(1);
		}
		if (! myid)	{
			Reset(force);
			// model->Update();
			cerr << "ln prob = " << model->GetLogProb() << "\n";
		}
	}

	void Open()	{

		ifstream is((name + ".param").c_str());
		if (!is)	{
			cerr << "error : cannot find file : " << name << ".param\n";
			exit(1);
		}
		is >> modeltype;
		is >> type;
		is >> datafile >> treefile >> contdatafile;
		is >> calibfile >> rootage >> rootstdev;
		is >> chronoprior >> meanchi >> meanchi2;
		if (chronoprior == 4)	{
			is >> softa;
		}
		is >> clampdiag >> autoregressive >> gc;
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
												is >> withtimeline;
												is >> check;
												if (check)	{
													is >> uniformprior;
													is >> rootfile;
													is >> check;
													if (check)	{
														is >> separatesyn;
														is >> separateomega;
														is >> check;
														if (check)	{
															is >> priorsigmafile;
															is >> check;
															if (check)	{
																is >> krkctype;
																is >> check;
																if (check)	{
																	is >> mappingfreq;
																	is >> check;
																	if (check)	{
																		is >> jitter;
																		is >> check;
																		if (check)	{
																			is >> ncycle;
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
			}
		}

		is >> every >> until >> size;

		if (modeltype == "BRANCHOMEGAMULTIVARIATE")	{
			model = new BranchOmegaMultivariateModel(datafile,treefile,contdatafile,calibfile,rootage,rootstdev,chronoprior,softa,meanchi,meanchi2,priorsigma,priorsigmafile,df,mutmodel,gc,clampdiag,autoregressive,conjpath,mappingfreq,contdatatype,omegaratiotree,clamproot,clamptree,meanexp,normalise,nrep,ncycle,bounds,mix,nsplit,withdrift,uniformprior,rootfile,suffstatfile,withtimeline,separatesyn,separateomega,krkctype,jitter,myid,nprocs,1,type);
		}
		else if (modeltype == "CONJUGATEBRANCHOMEGAMULTIVARIATE")	{
			model = new ConjugateBranchOmegaMultivariateModel(datafile,treefile,contdatafile,calibfile,rootage,rootstdev,chronoprior,softa,meanchi,meanchi2,priorsigma,priorsigmafile,df,mutmodel,gc,autoregressive,conjpath,mappingfreq,contdatatype,omegaratiotree,clamproot,clamptree,meanexp,normalise,nrep,ncycle,bounds,mix,nsplit,withdrift,uniformprior,rootfile,suffstatfile,withtimeline,separatesyn,separateomega,krkctype,jitter,myid,nprocs,1,type);
		}
		else	{
			cerr << "error when opening file "  << name << " : does not recognise model type : " << modeltype << '\n';
			exit(1);
		}
		if (! myid)	{
			model->FromStream(is);
			if (nprocs > 1)	{
				GetModel()->GlobalUpdateParameters();
			}
		}
		else	{
			GetModel()->SlaveUpdateParameters();
		}
		model->Update();
		if (! myid)	{
			cerr << size << " points saved, current ln prob = " << GetModel()->GetLogProb() << "\n";
		}
	}

	void Save()	{

		ofstream param_os((name + ".param").c_str());
		param_os << GetModelType() << '\n';
		param_os << type << '\n';
		param_os << datafile << '\t' << treefile << '\t' << contdatafile << '\n';
		param_os << calibfile << '\t' << rootage << '\t' << rootstdev << '\n';
		param_os << chronoprior << '\t' << meanchi << '\t' << meanchi2 << '\n';
		if (chronoprior == 4)	{
			param_os << softa << '\n';
		}
		param_os << clampdiag << '\t' << autoregressive << '\t' << gc << '\n';
		param_os << conjpath << '\n';
		param_os << contdatatype << '\n';
		param_os << omegaratiotree << '\n';
		param_os << clamproot << '\t' << meanexp << '\n';
		param_os << normalise << '\t' << nrep << '\n';
		param_os << 1 << '\n';
		param_os << mutmodel << '\n';
		param_os << 1 << '\n';
		param_os << df << '\n';
		param_os << 1 << '\n';
		param_os << bounds << '\n';
		param_os << 1 << '\n';
		param_os << priorsigma << '\n';
		param_os << 1 << '\n';
		param_os << mix << '\n';
		param_os << 1 << '\n';
		param_os << nsplit << '\n';
		param_os << 1 << '\n';
		param_os << withdrift << '\n';
		param_os << 1 << '\n';
		param_os << clamptree << '\n';
		param_os << 1 << '\n';
		param_os << suffstatfile << '\n';
		param_os << 1 << '\n';
		param_os << withtimeline << '\n';
		param_os << 1 << '\n';
		param_os << uniformprior << '\n';
		param_os << rootfile << '\n';
		param_os << 1 << '\n';
		param_os << separatesyn << '\n';
		param_os << separateomega << '\n';
		param_os << 1 << '\n';
		param_os << priorsigmafile << '\n';
		param_os << 1 << '\n';
		param_os << krkctype << '\n';
		param_os << 1 << '\n';
		param_os << mappingfreq << '\n';
		param_os << 1 << '\n';
		param_os << jitter << '\n';
		param_os << 1 << '\n';
		param_os << ncycle << '\n';
		param_os << 0 << '\n';
		param_os << every << '\t' << until << '\t' << size << '\n';

		model->ToStream(param_os);
	}

	void Monitor()	{
		Chain::Monitor();
		ofstream mon_os((name + ".chrono").c_str());
		mon_os << "total time per cycle : " << chrono.GetTime() / 1000 << '\n';
	}

	void MakeFiles(int force)	{
		if (! myid)	{
			Chain::MakeFiles(force);
		}
	}

	void Move()	{
		chrono.Reset();
		chrono.Start();

		for (int i=0; i<every; i++)	{
			model->Move(1);
		}

		chrono.Stop();
		if (! GetModel()->GetMyid())	{
			SavePoint();
			Save();
			Monitor();
		}
	}
};

int main(int argc, char* argv[])	{

	int myid = 0;
	int nprocs = 1;

#ifdef USE_MPI
	MPI_Init(&argc,&argv);
	MPI_Comm_rank(MPI_COMM_WORLD,&myid);
	MPI_Comm_size(MPI_COMM_WORLD,&nprocs);
#endif

	if (! myid)	{
		cerr << "Coevol version 1.6\n";
		cerr << '\n';
	}

	if (argc == 2)	{
		string name = argv[1];
		BranchOmegaMultivariateChain* chain = new BranchOmegaMultivariateChain(name,myid,nprocs);
		if (! myid)	{
			cerr << "start\n";
			cerr << '\n';
		}
		chain->Start();
		if (! myid)	{
			cerr << '\n';
			cerr << "exit\n";
			cerr << chain->GetSize() << " points saved, current ln prob = " << chain->GetModel()->GetLogProb() << "\n";
			cerr << '\n';
		}
	}
	else	{

		string datafile = "";
		string treefile = "";
		string contdatafile = "None";
		string calibfile = "None";
		string suffstatfile = "None";
		double rootage = 0;
		double rootstdev = 0;

		int chronoprior = 0;
		double softa = 0.025;
		double meanchi = -1;
		double meanchi2 = -1;

		string name = "";
		bool autoregressive = false;
		bool clampdiag = false;
		bool clamproot = false;
		bool clamptree = false;
		bool meanexp = false;
		bool conjugate = true;
		GeneticCodeType type = Universal;

		double priorsigma = -1;
		string priorsigmafile = "None";
		int df = 0;
		int gc = 0;
		int mutmodel = 0;

		int conjpath = -1;
		double mappingfreq = -1;
		int contdatatype = 0;

		int every = 1;
		int until = -1;

		int omegaratiotree = 0;

		int nrep = 0;
		int ncycle = 1;
		bool normalise = true;

		string bounds = "None";
		string mix = "None";

		int nsplit = 1;
		int withdrift = 0;
		int uniformprior = 0;
		string rootfile = "None";
		string sigmafile = "None";

		bool withtimeline = false;
		bool separatesyn = false;
		bool separateomega = false;

		int krkctype = -1;
		int splitkrkctype = 0;

		int jitter = 0;

		int force = 0;

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
				else if ((s == "-t") || (s == "-T"))	{
					i++;
					treefile = argv[i];
				}
				else if (s == "-c")	{
					i++;
					contdatafile = argv[i];
				}
                /*
				else if (s == "-nn")	{
					splitkrkctype = 1;
				}
				else if (s == "-krkc2")	{
					omegaratiotree = 5;
				}
				else if (s == "-krkc3")	{
					omegaratiotree = 4;
				}
                */
				else if (s == "-charge"){
					krkctype = CHARGE;
					// omegaratiotree = 1;
				}
				else if (s == "-pol"){
					krkctype = POLARITY;
					// omegaratiotree = 1;
				}
				else if (s == "-vol"){
					krkctype = VOLUME;
					// omegaratiotree = 1;
				}
				else if (s == "-polvol"){
					krkctype = POLARITYVOLUME;
					// omegaratiotree = 1;
				}
                /*
				else if (s == "-jitter")	{
					jitter = 1;
				}
				else if (s == "-wnjitter")	{
					jitter = 2;
				}
				else if (s == "-contjitter")	{
					jitter = 3;
				}
				else if (s == "-wncontjitter")	{
					jitter = 4;
				}
                */
				else if (s == "-suffstat")	{
					i++;
					suffstatfile = argv[i];
				}
				else if (s == "-f")	{
					force = 1;
				}
                /*
				else if (s == "-timeline")	{
					withtimeline = true;
				}
                */
				else if (s == "-sepsyn")	{
					separatesyn = true;
				}
				else if (s == "-sepomega")	{
					separateomega = true;
				}
				else if (s == "-linear")	{
					contdatatype = 2;
				}
				else if (s == "-logit")	{
					contdatatype = 1;
				}
                /*
				else if (s == "-drift")	{
					withdrift = 1;
				}
				else if (s == "-expdrift")	{
					withdrift = 2;
				}
				else if (s == "-doubleexpdrift")	{
					withdrift = 3;
				}
				else if (s == "-reldrift")	{
					withdrift = 4;
				}
				else if (s == "-unidrift")	{
					uniformprior = 1;
				}
                */
				else if (s == "-root")	{
					i++;
					rootfile = argv[i];
				}
                /*
				else if (s == "-split")	{
					i++;
					s = argv[i];
					i++;
					nsplit = atoi(argv[i]);
					if (s == "length")	{
						nsplit = -nsplit;
					}
				}
				else if (s == "-bounds")	{
					i++;
					bounds = argv[i];
				}
				else if (s == "-mix")	{
					i++;
					mix = argv[i];
				}
				else if (s == "-uncons")	{
					calibfile = "Unconstrained";
				}
                */
				else if (s == "-cal")	{
					i++;
					calibfile = argv[i];
					i++;
					if ((i == argc) || (! IsFloat(argv[i])))	{
						cerr << "error in command: -cal <mean> <stdev>\n";
						cerr << '\n';
						exit(1);
					}
					rootage = atof(argv[i]);
					i++;
					if ((i == argc) || (! IsFloat(argv[i])))	{
						cerr << "error in command: -cal <mean> <stdev>\n";
						cerr << '\n';
						exit(1);
					}
					rootstdev = atof(argv[i]);
				}
				else if (s == "-unif")	{
					chronoprior = 0;
				}
				else if (s == "-bd")	{
					chronoprior = 1;
				}
                /*
				else if (s == "-cbd")	{
					chronoprior = 2;
				}
				else if (s == "-rbd")	{
					chronoprior = 3;
				}
				else if (s == "-sbd")	{
					chronoprior = 4;
					i++;
					softa = atof(argv[i]);
				}
                */
				else if (s == "-bdhyperprior")	{
					i++;
					if ((i == argc) || (! IsFloat(argv[i])))	{
						cerr << "error in command: -bd <p1> <p2>\n";
						cerr << '\n';
						exit(1);
					}
					meanchi = atof(argv[i]);
					i++;
					if ((i == argc) || (! IsFloat(argv[i])))	{
						cerr << "error in command: -bd <p1> <p2>\n";
						cerr << '\n';
						exit(1);
					}
					meanchi2 = atof(argv[i]);
				}
				else if ((s == "-priorsigma") || (s == "-priornu"))	{
					i++;
					if ((i == argc) || (! IsFloat(argv[i])))	{
						cerr << "error in command: -priorsigma <priorsigma>\n";
						cerr << '\n';
						exit(1);
					}
					priorsigma = atof(argv[i]);
					if (priorsigma == -2)	{
						i++;
						priorsigmafile = argv[i];
					}
				}
				else if (s == "-df")	{
					i++;
					df = atoi(argv[i]);
				}
				else if (s == "-conj")	{
					conjugate = true;
				}
				else if (s == "-nonconj")	{
					conjugate = false;
				}
				else if (s == "-simpleconjpath")	{
					conjpath = 1;
				}
				else if (s == "-fastconjpath")	{
					conjpath = 2;
				}
				else if (s == "-nonconjpath")	{
					conjpath = 0;
				}
				else if (s == "-prior")	{
					conjpath = 3;
				}
				else if (s == "-mapfreq")	{
					i++;
					mappingfreq = atof(argv[i]);
				}
				else if (s == "-norm")	{
					normalise = true;
				}
				else if (s == "-nonnorm")	{
					normalise = false;
				}
				else if (s == "-nrep")	{
					i++;
					nrep = atoi(argv[i]);
				}
				else if (s == "-ncycle")	{
					i++;
					ncycle = atoi(argv[i]);
				}
				else if ((s == "-dSdN") || (s == "-dsdn"))	{
					omegaratiotree = 2;
				}
				else if ((s == "-dSomega") || (s == "-dsomega") || (s == "-dsom"))	{
					omegaratiotree = 1;
				}
				else if (s == "-dsom3")	{
					omegaratiotree = 3;
				}
				else if (s == "-dsom2")	{
					omegaratiotree = 4;
				}
				/*
				else if (s == "-hky")	{
					mutmodel = 1;
				}
				else if (s == "-hky2")	{
					mutmodel = 2;
				}
				else if (s == "-hky3")	{
					mutmodel = 3;
				}
				*/
				else if (s == "-tstv")	{
					mutmodel = 4;
				}
				else if (s == "-tstvgc")	{
					mutmodel = 5;
				}
                /*
				else if (s == "-gcstat")	{
					gc = -1;
				}
                */
				else if (s == "-gc")	{
					gc = 1;
				}
                /*
				else if (s == "-gc3")	{
					gc = 3;
				}
				else if (s == "-gc2")	{
					gc = 2;
				}
				else if (s == "-oup")	{
					autoregressive = true;
				}
				else if (s == "-bp")	{
					autoregressive = false;
				}
				else if (s == "-clamproot")	{
					clamproot = true;
				}
                */
				else if ((s == "-fixtimes") || (s == "-fixbl"))	{
					clamptree = true;
				}
				else if (s == "-meanexp")	{
					meanexp = true;
				}
				else if (s == "-arith")	{
					meanexp = false;
				}
				else if (s == "-geod")	{
					meanexp = true;
				}
				else if (s == "-diag")	{
					clampdiag = true;
					conjugate = false;
				}
				else if (s == "-mtvert")	{
					type = MtMam;
				}
				else if (s == "-mtmam")	{
					type = MtMam;
				}
				else if (s == "-mtinv")	{
					type = MtInv;
				}
				else if (s == "-univ")	{
					type = Universal;
				}
				else if ( (s == "-x") || (s == "-extract") )	{
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
			/*
			if ((datafile == "") || (treefile == "") || (name == ""))	{
				throw(0);
			}
			*/
			if (contdatafile == "")	{
				contdatafile = "None";
			}
			if (clampdiag && conjugate)	{
				cerr << "error : conjugate sampling and diagonal model are not compatible\n";
				exit(1);
			}
		}
		catch(...)	{
			cerr << "coevol -d <alignment> -t <tree> [-c <trait_data>]  <chainname> \n";
			cerr << "see manual for details\n";
			cerr << '\n';
			exit(1);
		}


		if ((krkctype != -1) && (! omegaratiotree))	{
			omegaratiotree = 1;
		}
		if (splitkrkctype)	{
			krkctype += 4;
		}
		if (chronoprior > 1)	{
			cerr << "only classical hard bounds\n";
			exit(1);
		}
		BranchOmegaMultivariateChain* chain = new BranchOmegaMultivariateChain(datafile,treefile,contdatafile,calibfile,rootage,rootstdev,chronoprior,softa,meanchi,meanchi2,priorsigma,priorsigmafile,df,mutmodel,gc,name,clampdiag,autoregressive,type,conjugate,conjpath,mappingfreq,contdatatype,omegaratiotree,clamproot,clamptree,meanexp,normalise,nrep,ncycle,bounds,mix,nsplit,withdrift,uniformprior,rootfile,suffstatfile,withtimeline,separatesyn,separateomega,krkctype,jitter,myid,nprocs,force);
		chain->SetEvery(every);
		chain->SetUntil(until);
		if (!myid)	{
			cerr << "start\n";
			cerr << '\n';
		}
		chain->Start();
		if (!myid)	{
			cerr << '\n';
			cerr << "exit\n";
			cerr << chain->GetSize() << " points saved, current ln prob = " << chain->GetModel()->GetLogProb() << "\n";
			cerr << '\n';
		}
	}
#ifdef USE_MPI
	MPI_Finalize();
#endif
	exit(0);
}

