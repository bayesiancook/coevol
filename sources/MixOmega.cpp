
#include "Chain.h"
#include "MixOmegaModel.h"
#include "StringStreamUtils.h"
#include "Chrono.h"

#include "Parallel.h"

// int BGCMutSelSubMatrix::bgccount = 0;

MPI_Datatype Propagate_arg;

class MixOmegaChain : public Chain	{

	private:

	int burnindone;
	int myid;
	int nprocs;

	string modeltype;
	string datafile;
	string treefile;
	string contdatafile;
	string calibfile;

	string ontologyfile;
	int synontology;
	int omegaontology;
	int clampbetatoggle;

	string suffstatfile;
	string rootfile;

	double rootage;
	double rootstdev;

	int chronoprior;
	double meanchi;
	double meanchi2;

	int clampdiag;
	int clamptree;
	int meanexp;

	GeneticCodeType type;

	int conjpath;

	int normalise;
	int nrep;

	int df;
	int clampreg;

	double mappingfreq;
	double mappingevery;
	string omegafile;

	public:

	MixOmegaModel* GetModel() {return (MixOmegaModel*) model;}

	string GetModelType() {return modeltype;}

	MixOmegaChain(string indata, string intree, string incontdata, string incalibfile, string inontologyfile, int insynontology, int inomegaontology, int inclampbetatoggle, double inrootage, double inrootstdev, int inchronoprior, double inmeanchi, double inmeanchi2, int indf, GeneticCodeType intype, string filename, int inclampdiag, int inconjpath, int inclamptree, int inmeanexp, int innormalise, int innrep, string insuffstatfile, string inrootfile, int inclampreg, double inmappingfreq, double inmappingevery, string inomegafile, int inmyid, int innprocs, int force = 1)	{
		burnindone = 0;
		myid = inmyid;
		nprocs = innprocs;
		modeltype = "MIXOMEGA";
		datafile = indata;
		treefile = intree;
		contdatafile = incontdata;
		calibfile = incalibfile;
		ontologyfile = inontologyfile;
		synontology = insynontology;
		omegaontology = inomegaontology;
		clampbetatoggle = inclampbetatoggle;
		suffstatfile = insuffstatfile;
		rootfile = inrootfile;
		rootage = inrootage;
		rootstdev = inrootstdev;
		chronoprior = inchronoprior;
		meanchi = inmeanchi;
		meanchi2 = inmeanchi2;
		name = filename;
		clampdiag = inclampdiag;
		df = indf;
		type = intype;
		conjpath = inconjpath;
		clamptree = inclamptree;
		meanexp = inmeanexp;
		normalise = innormalise;
		nrep = innrep;
		clampreg = inclampreg;
		mappingfreq = inmappingfreq;
		mappingevery = inmappingevery;
		omegafile = inomegafile;
		New(force);
	}

	MixOmegaChain(string filename, int inmyid, int innprocs)	{
		name = filename;
		myid = inmyid;
		nprocs = innprocs;
		Open();
	}

	void New(int force)	{
		if (modeltype == "MIXOMEGA")	{
			model = new MixOmegaModel(datafile,treefile,contdatafile,calibfile,ontologyfile,synontology,omegaontology,clampbetatoggle,rootage,rootstdev,chronoprior,meanchi,meanchi2,df,clampdiag,conjpath,clamptree,meanexp,normalise,nrep,suffstatfile,rootfile,clampreg,mappingfreq,mappingevery,omegafile,type,myid,nprocs,true);
		}
		else	{
			cerr << "error, does not recognise model type : " << modeltype << '\n';
			exit(1);
		}
		if (! myid)	{
			Reset(force);
		}
	}

	void Open()	{
		ifstream is((name + ".param").c_str());
		if (!is)	{
			cerr << "error : cannot find file : " << name << ".param\n";
			exit(1);
		}
		// burnindone = 1;

		ontologyfile = "None";
		synontology = 0;
		omegaontology = 0;
		clampbetatoggle = 1;

		is >> modeltype;
		is >> datafile >> treefile >> contdatafile;
		is >> calibfile >> rootage >> rootstdev;
		is >> chronoprior >> meanchi >> meanchi2;

		is >> clampdiag;
		is >> conjpath;
		is >> meanexp;
		is >> normalise >> nrep;
		is >> df;
		is >> type;
		is >> clamptree;
		is >> suffstatfile;
		is >> rootfile;
		is >> clampreg;
		is >> mappingfreq;
		is >> mappingevery;
		is >> omegafile;

		int check;
		is >> check;
		if (check)	{
			is >> ontologyfile;
			is >> synontology >> omegaontology >> clampbetatoggle;

			is >> check;
			if (check)	{
				cerr << "error when reading model\n";
				exit(1);
			}
		}

		is >> every >> until >> size;

		if (modeltype == "MIXOMEGA")	{
			model = new MixOmegaModel(datafile,treefile,contdatafile,calibfile,ontologyfile,synontology,omegaontology,clampbetatoggle,rootage,rootstdev,chronoprior,meanchi,meanchi2,df,clampdiag,conjpath,clamptree,meanexp,normalise,nrep,suffstatfile,rootfile,clampreg,mappingfreq,mappingevery,omegafile,type,myid,nprocs,false);
		}
		else	{
			cerr << "error when opening file "  << name << " : does not recognise model type : " << modeltype << '\n';
			exit(1);
		}
		// model->FromStream(is);
		if (! myid)	{
			GetModel()->FromStream(is);
			GetModel()->MasterSendParameters();
			GetModel()->MasterSendGeneProcesses();
			GetModel()->Update();
			GetModel()->TraceHeader(cerr);
			GetModel()->Trace(cerr);
		}
		else	{
			GetModel()->SlaveReceiveParameters();
			GetModel()->SlaveReceiveGeneProcesses();
			GetModel()->Update();
			GetModel()->SlaveTrace();
		}
		cerr << "model restarted\n";
		// cerr << size << " points saved, current ln prob = " << GetModel()->GetLogProb() << "\n";
	}

	void Save()	{
		ofstream param_os((name + ".param").c_str());
		param_os << GetModelType() << '\n';
		param_os << datafile << '\t' << treefile << '\t' << contdatafile << '\n';
		param_os << calibfile << '\t' << rootage << '\t' << rootstdev << '\n';
		param_os << chronoprior << '\t' << meanchi << '\t' << meanchi2 << '\n';
		param_os << clampdiag << '\n';
		param_os << conjpath << '\n';
		param_os << meanexp << '\n';
		param_os << normalise << '\t' << nrep << '\n';
		param_os << df << '\n';
		param_os << type << '\n';
		param_os << clamptree << '\n';
		param_os << suffstatfile << '\n';
		param_os << rootfile << '\n';
		param_os << clampreg << '\n';
		param_os << mappingfreq << '\n';
		param_os << mappingevery << '\n';
		param_os << omegafile << '\n';
		param_os << 1 << '\n';
		param_os << ontologyfile << '\n';
		param_os << synontology << '\t' << omegaontology << '\t' << clampbetatoggle << '\n';
		param_os << 0 << '\n';
		param_os << every << '\t' << until << '\t' << size << '\n';
		model->ToStream(param_os);
	}

	void Start()	{
		if (! myid)	{
			ofstream run_os((name + ".run").c_str());
			run_os << 1 << '\n';
			run_os.close();
		}
		Run();
	}

	void Run()	{
		while (GetRunningStatus() && ((until == -1) || (size <= until)))	{
			Chrono chrono;
			chrono.Reset();
			chrono.Start();
			Move();
			chrono.Stop();
			ofstream check_os((name + ".time").c_str());
			check_os << chrono.GetTime() / 1000 << '\n';
		}
		if (! myid)	{
			ofstream run_os((name + ".run").c_str());
			run_os << 0 << '\n';
		}
	}

	int GetRunningStatus()	{
		if (! myid)	{
			ifstream run_is((name + ".run").c_str());
			int run;
			run_is >> run;
			MPI_Bcast(&run,1,MPI_INT,0,MPI_COMM_WORLD);
			return run;
		}
		else	{
			int run;
			MPI_Bcast(&run,1,MPI_INT,0,MPI_COMM_WORLD);
			return run;
		}
	}

	void Move()	{

		/*
		if (! burnindone)	{
			GetModel()->PreMove(100);
			burnindone=1;
		}
		*/

		for (int i=0; i<every; i++)	{
			model->Move(1);
		}

		if (! myid)	{
			SavePoint();
			Save();
			Monitor();
		}
		else	{
			SlaveMonitor();
		}
	}

	void SlaveMonitor()	{
		GetModel()->SlaveTrace();
		ostringstream s;
		s << name << myid << ".monitor";
		ofstream mon_os(s.str().c_str());
		ostringstream s2;
		s2 << name << myid << ".detailedmonitor";
		ofstream mond_os(s2.str().c_str());
		model->Monitor(mon_os,mond_os);
	}
};

int main(int argc, char* argv[])	{

	int myid,nprocs;

	MPI_Init(&argc,&argv);
	MPI_Comm_rank(MPI_COMM_WORLD,&myid);
	MPI_Comm_size(MPI_COMM_WORLD,&nprocs);

	int blockcounts[2] = {1,3};
	MPI_Datatype types[2] = {MPI_DOUBLE,MPI_INT};
	MPI_Aint dtex,displacements[2];

	displacements[0] = (MPI_Aint) 0;
	MPI_Type_extent(MPI_DOUBLE,&dtex);
	displacements[1] = dtex;
	MPI_Type_struct(2,blockcounts,displacements,types,&Propagate_arg);
	MPI_Type_commit(&Propagate_arg);

	if (! myid)	{
		cerr << '\n';
	}

	if (argc == 2)	{

		string name = argv[1];
		MixOmegaChain* chain = new MixOmegaChain(name,myid,nprocs);
		cerr << "start\n";
		cerr << '\n';
		chain->Start();
		cerr << '\n';
		cerr << "exit\n";
		cerr << chain->GetSize() << " points saved, current ln prob = " << chain->GetModel()->GetLogProb() << "\n";
		cerr << '\n';
	}
	else	{

		GeneticCodeType type = Universal;

		string datafile = "";
		string treefile = "";
		string contdatafile = "None";
		string calibfile = "None";

		string ontologyfile = "None";
		int synontology = 0;
		int omegaontology = 0;
		int clampbetatoggle = 1;

		string suffstatfile = "None";
		string rootfile = "None";
		double rootage = 0;
		double rootstdev = 0;

		int chronoprior = 0;
		double meanchi = 1e-3;
		double meanchi2 = 1e-3;

		string name = "";
		int clampdiag = 2;
		int clamptree = 1;
		int meanexp = 0;

		int df = 0;

		int conjpath = -1;

		int every = 50;
		int until = -1;

		int nrep = 1;
		int normalise = 0;

		int force = 0;

		int clampreg = 1;

		double mappingfreq = 1;
		double mappingevery = 50;

		string omegafile = "None";

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
				else if (s == "-mapfreq")	{
					i++;
					mappingfreq = atof(argv[i]);
				}
				else if (s == "-mapevery")	{
					i++;
					mappingevery = atof(argv[i]);
				}
				else if (s == "-fixom")	{
					i++;
					omegafile = argv[i];
				}
				else if (s == "-suffstat")	{
					i++;
					suffstatfile = argv[i];
				}
				else if (s == "-clampreg")	{
					clampreg = 2;
				}
				else if (s == "-clamptoggle")	{
					clampreg = 1;
				}
				else if (s == "-toggle")	{
					clampreg = 0;
				}
				else if (s == "-o")	{
					i++;
					ontologyfile = argv[i];
					omegaontology = 1;
				}
				else if ((s == "-uni") || (s == "-uniform"))	{
					clampreg = 3;
				}
				else if (s == "-root")	{
					i++;
					rootfile = argv[i];
				}
				else if (s == "-f")	{
					force = 1;
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
				else if (s == "-bd")	{
					chronoprior = 1;
				}
				else if (s == "-cbd")	{
					chronoprior = 2;
				}
				else if (s == "-rbd")	{
					chronoprior = 3;
				}
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
				else if (s == "-df")	{
					i++;
					df = atoi(argv[i]);
				}
				else if (s == "-conjpath")	{
					conjpath = 1;
				}
				else if (s == "-nonconjpath")	{
					conjpath = 0;
				}
				else if (s == "-prior")	{
					conjpath = 2;
				}
				else if (s == "-norm")	{
					normalise = 1;
				}
				else if (s == "-nonnorm")	{
					normalise = 0;
				}
				else if (s == "-nrep")	{
					i++;
					nrep = atoi(argv[i]);
				}
				else if (s == "-fixbl")	{
					clamptree = 1;
				}
				else if (s == "-meanexp")	{
					meanexp = 1;
				}
				else if (s == "-arith")	{
					meanexp = 0;
				}
				else if (s == "-geod")	{
					meanexp = 1;
				}
				else if (s == "-diag")	{
					clampdiag = 1;
				}
				else if (s == "-contml")	{
					clampdiag = 2;
				}
				else if (s == "-contbayes")	{
					clampdiag = 0;
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
			if ((datafile == "") || (treefile == "") || (name == ""))	{
				throw(0);
			}
			if (contdatafile == "")	{
				contdatafile = "None";
			}
		}
		catch(...)	{
			cerr << "error in command\n";
			cerr << '\n';
			exit(1);
		}

		MixOmegaChain* chain = new MixOmegaChain(datafile,treefile,contdatafile,calibfile,ontologyfile,synontology,omegaontology,clampbetatoggle,rootage,rootstdev,chronoprior,meanchi,meanchi2,df,type,name,clampdiag,conjpath,clamptree,meanexp,normalise,nrep,suffstatfile,rootfile,clampreg,mappingfreq,mappingevery,omegafile,myid,nprocs,force);

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
			cerr << chain->GetSize() << " points saved\n";
			cerr << '\n';
		}
	}

	MPI_Finalize();
}

