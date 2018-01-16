
#include "Sample.h"
#include "DirichletCodonUsageSelectionModelMS.h"
#include "GTRSubMatrix.h"

class DirichletCodonUsageSelectionMSSample : public Sample	{

	private:
	string modeltype;
	string datafile;
    string contdatafile;
	string treefile;
    int ncond;
    int fixglob;
    int codonmodel;

    double alpha;

	public:

	string GetModelType() {return modeltype;}

	DirichletCodonUsageSelectionModelMS* GetModel() {return (DirichletCodonUsageSelectionModelMS*) model;}

	DirichletCodonUsageSelectionMSSample(string filename, int inburnin, int inevery, int inuntil, double inalpha = 0.1) : Sample(filename,inburnin,inevery,inuntil)	{
		alpha = inalpha;
		Open();
	}

	void Open()	{

		// open <name>.param
		ifstream is((name + ".param").c_str());

		// check that file exists
		if (!is)	{
			cerr << "error : cannot find file : " << name << ".param\n";
			exit(1);
		}

		// read model type, and other standard fields
		is >> modeltype;
		is >> datafile >> contdatafile >> treefile;
        is >> ncond;
        is >> fixglob >> codonmodel;
		is >> chainevery >> chainuntil >> chainsize;

		if (modeltype == "DIFFSELDIR")	{
			cerr << "CREATE\n";
			model = new DirichletCodonUsageSelectionModelMS(datafile,contdatafile,treefile,ncond,fixglob,codonmodel,false);
		}
		else	{
			cerr << "error when opening file "  << name << '\n';
			exit(1);
		}

		// read model (i.e. chain's last point) from <name>.param
		cerr << "from stream\n";
		model->FromStream(is);

		/*
		cerr << "UPDATE\n";
		model->Update();
		*/

		// open <name>.chain, and prepare stream and stream iterator
		OpenChainFile();
		// now, size is defined (it is the total number of points with which this Sample object will make all its various posterior averages)
		// all these points can be accessed to (only once) by repeated calls to GetNextPoint()
	}

	void ReadFalsePositives(int ncat)	{

		int Nsite = GetModel()->GetSite();
		int Naa = GetModel()->GetaaState();
		int K = GetModel()->GetCategory();

		double*** meansel = new double**[K];
		double*** ppsel = new double**[K];

		int kmin = 1;

		for (int k=kmin; k<K; k++)	{
			meansel[k] = new double*[Nsite];
			ppsel[k] = new double*[Nsite];
			for (int i=0; i<Nsite; i++)	{
				meansel[k][i] = new double[Naa];
				ppsel[k][i] = new double[Naa];
				for (int j=0; j<Naa; j++)	{
					meansel[k][i][j] = 0;
					ppsel[k][i][j] = 0;
				}
			}
		}

		// cycle over the sample
		for (int c=0; c<size; c++)	{
			GetNextPoint();
			cerr << '.';

			//get selection profile
			for(int k=kmin;k<K;k++)	{
				for(int i=0;i<Nsite;i++)	{
					double logsum = 0; 
					for(int j=0;j<Naa;j++)	{ 
						double tmp = log(GetModel()->GetSelectionProfile(k, i, j)) - log(GetModel()->GetSelectionProfile(k-1,i,j));
						logsum += tmp;
					}
					for(int j=0;j<Naa;j++)	{  
						double tmp = log(GetModel()->GetSelectionProfile(k, i, j)) - log(GetModel()->GetSelectionProfile(k-1,i,j)) - (1.0/20) * (logsum);
						meansel[k][i][j] += tmp;
						if(tmp > 0)	{
							ppsel[k][i][j]++;
						}
					}
				}
			}
		}
		cerr << '\n';

		double** falsepositives = new double*[K];
		for(int k=kmin;k<K;k++)	{
			falsepositives[k] = new double[ncat];
			for (int cat=0; cat<ncat; cat++)	{
				falsepositives[k][cat] = 0;
			}
		}
		for(int k=kmin;k<K;k++)	{
			for(int i=0;i<Nsite;i++)	{
				for(int j=0;j<Naa;j++)	{ 
					meansel[k][i][j] /= size;
					ppsel[k][i][j] /= size;
					int cat = (int) (10*(ppsel[k][i][j] - 0.5));
					if (cat < 0)	{
						cat = -cat;
					}
					if (cat == ncat)	{
						cat--;
					}
					falsepositives[k][cat]++;
				}
			}
		}

		for(int k=kmin;k<K;k++)	{
			for (int cat=0; cat<ncat; cat++)	{
				for (int cat2=cat+1; cat2<ncat; cat2++)	{
					falsepositives[k][cat] +=falsepositives[k][cat2];
				}
			}
		}

		for(int k=kmin;k<K;k++)	{
			for (int cat=0; cat<ncat; cat++)	{
				falsepositives[k][cat] /= Nsite*Naa;
			}
		}


		ofstream cout((GetName() + ".fp").c_str());
		cout << '\n';
		cout << "false positive rates:\n";
		cout << '\n';
		for(int k=kmin;k<K;k++)	{
			cout << "category : " << k << '\n';
			for (int cat=0; cat<ncat; cat++)	{
				cout << 0.5 + ((double) cat) / 2 / ncat << '\t' << 0.5 + ((double) (cat+1)) / 2 / ncat << '\t';
				cout << 100 * falsepositives[k][cat] << '\t' << Nsite*Naa * falsepositives[k][cat] << '\n';
			}
			cout << '\n';
		}
	}

	void NewRead(double cutoff)	{

		// prepare the mean and variance
		
		int Nsite = GetModel()->GetSite();
		int K = GetModel()->GetCategory();

		double*** meansel = new double**[K];
		for (int k=1; k<K; k++)	{
			meansel[k] = new double*[Nsite];
			for (int i=0; i<Nsite; i++)	{
				meansel[k][i] = new double[Naa];
				for (int a=0; a<Naa; a++)	{
					meansel[k][i][a] = 0;
				}
			}
		}

		double*** ppsel = new double**[K];
		for (int k=1; k<K; k++)	{
			ppsel[k] = new double*[Nsite];
			for (int i=0; i<Nsite; i++)	{
				ppsel[k][i] = new double[Naa];
				for (int a=0; a<Naa; a++)	{
					ppsel[k][i][a] = 0;
				}
			}
		}

		// cycle over the sample
		for (int c=0; c<size; c++)	{
			GetNextPoint();
			cerr << '.';

			for (int k=1; k<K; k++)	{
				for (int i=0; i<Nsite; i++)	{

					double delta[Naa];
					for(int j=0;j<Naa;j++)	{ 
						delta[j] = log(GetModel()->GetSelectionProfile(k, i, j)) - log(GetModel()->GetSelectionProfile(k-1,i,j));
					}

					double mean = 0;
					for (int a=0; a<Naa; a++)	{
						mean += delta[a];
					}
					mean /= Naa;

					for (int a=0; a<Naa; a++)	{
						double tmp = delta[a] - mean;
						meansel[k][i][a] += tmp;
						if (tmp > 0)	{
							ppsel[k][i][a] ++;
						}
					}
				}
			}
		}
		cerr << '\n';

		for (int k=1; k<K; k++)	{

			ostringstream s1,s2;
			s1 << GetName() << "_" << k << ".meandiffsel";
			s2 << GetName() << "_" << k << ".signdiffsel";
			ofstream os1(s1.str().c_str());
			ofstream os2(s2.str().c_str());

			for (int i=0; i<Nsite; i++)	{

				os1 << i;
				for (int a=0; a<Naa; a++)	{
					ppsel[k][i][a] /= size;
					os1 << '\t' << ppsel[k][i][a];
					
				}
				for (int a=0; a<Naa; a++)	{
					meansel[k][i][a] /= size;
					os1 << '\t' << meansel[k][i][a];
				}
				os1 << '\n';

				for (int a=0; a<Naa; a++)	{
                    if ((ppsel[k][i][a] > cutoff) || (ppsel[k][i][a] < 1-cutoff)) {
						os2 << i << '\t' << AminoAcids[a] << '\t' << ppsel[k][i][a] << '\t' << meansel[k][i][a] << '\n';
					}
				}
			}
		}

		ofstream os((GetName() + ".postmeansel").c_str());
		for (int i=0; i<Nsite; i++)	{
			os << i << '\t';
			int k = K-1;
			// for (int k=1; k<K; k++)	{
			for (int a=0; a<Naa; a++)	{
				os << meansel[k][i][a] << '\t' << ppsel[k][i][a] << '\t';
			}
			// }
			os << '\n';
		}
	}

	// a very simple (and quite uninteresting) method for obtaining
	// the posterior mean and variance of the total length of the tree
	void Read(double cutoff)	{

		// prepare the mean and variance
		
		int Nsite = GetModel()->GetSite();

		int K = GetModel()->GetCategory();

		double selectionarray[K][Nsite][Naa];
		double pp[K][Nsite][Naa];
		double absolutarray[K][Nsite][Naa];
		double differentialfitness[K-1];
		double c3c4array[Nsite][Naa];
		double pp1[Nsite][Naa];


		for(int k=0;k<K;k++)	{
			for(int i=0;i<Nsite;i++)	{
				for(int j=0;j<Naa;j++)	{
					selectionarray[k][i][j] = 0;
					absolutarray[k][i][j] = 0;
					pp[k][i][j] = 0;
					c3c4array[i][j] = 0;
					pp1[i][j] = 0;
				}
			}
		}

		for(int i=0;i<K-1;i++)	{
			differentialfitness[i] = 0;
		}
	
		ostringstream delta;
		delta << GetName() +  ".deltacondition";
		string deltaname = delta.str();
		ofstream osdelta(deltaname.c_str());
	

		// cycle over the sample
		for (int c=0; c<size; c++)	{
			GetNextPoint();
			cerr << '*';

			//for global aveselection
			for(int i=0;i<Nsite;i++)	{
				for(int j=0;j<Naa;j++)	{ 
					double tmp = GetModel()->GetSelectionProfile(0,i,j);
					if (tmp > 0.00001)	{
						selectionarray[0][i][j] += tmp;
					}
					pp[0][i][j]++;
				}
			}

			//for absolutarray for all conditions
			for(int k=0;k<K;k++)	{
				for(int i=0;i<Nsite;i++)	{
					for(int j=0;j<Naa;j++)	{ 
						double tmp = GetModel()->GetSelectionProfile(k,i,j);
						if (tmp > 0.00001)	{
							absolutarray[k][i][j] += tmp;
						}
					}
				}
			}
			

			for(int k=1;k<K;k++)	{
				for(int i=0;i<Nsite;i++)	{
					/*
					double logsum = 0; 
					for(int j=0;j<Naa;j++)	{ 
						logsum += log(GetModel()->GetSelectionProfile(k, i, j)) - log(GetModel()->GetSelectionProfile(k-1,i,j));
					}
					logsum /= Naa;
					*/
					double total = 0;
					for(int j=0;j<Naa;j++)	{  
						double tmp = log(GetModel()->GetSelectionProfile(k, i, j)) - log(GetModel()->GetSelectionProfile(k-1,i,j));
						// double tmp = GetModel()->GetSelectionProfile(k, i, j) - GetModel()->GetSelectionProfile(k-1,i,j);
						total += tmp;
					}
					total /= Naa;

					if(k>2)	{
						for(int j=0;j<Naa;j++)	{  
							double tmp = log(GetModel()->GetSelectionProfile(k, i, j)) - log(GetModel()->GetSelectionProfile(1,i,j));
							total += tmp;
						}
						total /= Naa;
					}

					for(int j=0;j<Naa;j++)	{  
						double tmp = log(GetModel()->GetSelectionProfile(k, i, j)) - log(GetModel()->GetSelectionProfile(k-1,i,j)) - total;
						if(k>2)	{
							tmp = log(GetModel()->GetSelectionProfile(k, i, j)) - log(GetModel()->GetSelectionProfile(1,i,j)) - total;
						}
						if (fabs(tmp) > 0.00001)	{
						// double tmp = log(GetModel()->GetSelectionProfile(k, i, j)) - log(GetModel()->GetSelectionProfile(k-1,i,j)) - logsum;
							selectionarray[k][i][j] += tmp;
						}
						if (tmp > 0)	{
							pp[k][i][j]++;
						}
					}
					
					if (k>1){
						for(int j=0;j<Naa;j++)	{  
							double tmp = log(GetModel()->GetSelectionProfile(k, i, j)) - log(GetModel()->GetSelectionProfile(k-1,i,j)) - total;
							if (fabs(tmp) > 0.00001)	{
								c3c4array[i][j] += tmp;
							}
							if (tmp > 0)	{
								pp1[i][j]++;
							}
						}
					  
					}  
				}
			}
			
			// get true fitness difference
			for(int k=0;k<K-1;k++)	{
				double total = 0;
				for(int i=0;i<Nsite;i++)	{
					GetModel()->Getsubmatrix(k,i)->UpdateMatrix();
					for(int j=0;j<Naa;j++)	{    
						total += ((GetModel()->Getsubmatrix(k+1,i)->Stationary(j)) - (GetModel()->Getsubmatrix(k,i)->Stationary(j))) * log(GetModel()->GetSelectionProfile(k+1, i, (GetModel()->GetCodonStateSpace()->Translation(j))));
					}
				}
				differentialfitness[k] += total;
				osdelta << total << '\t';
			}
			osdelta << '\n';
			
			
		}
		cerr << '\n';

		for(int k=0;k<K;k++)	{
			for(int i=0;i<Nsite;i++)	{
				for(int j=0;j<Naa;j++)	{ 
					selectionarray[k][i][j] /= size;
					absolutarray[k][i][j] /= size;
					pp[k][i][j] /= size;
				}
			}
		}
		
		for(int i=0;i<Nsite;i++)	{
			for(int j=0;j<Naa;j++)	{ 
				c3c4array[i][j] /= size;
				pp1[i][j] /= size;
			}
		}
		

		for(int k=0;k<K;k++)	{
			ostringstream s;
			s << GetName() + "_" << k << ".aveselection";
			string name = s.str();
			ofstream os(name.c_str());
					
			os << Nsite << '\n';
			for(int i=0;i<Nsite;i++)	{
				for(int j=0;j<Naa;j++)	{    
					os << selectionarray[k][i][j] << '\t';
				}
				os << '\n';
			}
		}
		
		for(int k=0;k<K;k++)	{
			ostringstream abs;
			abs << GetName() + "_" << k << ".abselection";
			string abname = abs.str();
			ofstream oabs(abname.c_str());
					
			oabs << Nsite << '\n';
			for(int i=0;i<Nsite;i++)	{
				for(int j=0;j<Naa;j++)	{    
					oabs << absolutarray[k][i][j] << '\t';
				}
				oabs << '\n';
			}
		}


		ostringstream cs;
		cs << GetName() << ".c3c4selection";
		string cname = cs.str();
		ofstream ocs(cname.c_str());

		ostringstream cps;
		cps << GetName() << ".c3c4pp";
		string cpname = cps.str();
		ofstream ocps(cpname.c_str());

		ocs << Nsite << '\n';

		for(int i=0;i<Nsite;i++)	{
			for(int j=0;j<Naa;j++)	{  
				ocs << c3c4array[i][j] << '\t';
 				if (pp1[i][j] > cutoff)	{
					ocps << i << '\t' << AminoAcids[j] << '\t' << pp1[i][j] << '\n';
				}
				if (pp1[i][j] < 1-cutoff)	{
					ocps << i << '\t' << AminoAcids[j] << '\t' << 1-pp1[i][j] << '\n';
				}
			}
			ocs << '\n';
			
		} 

		ostringstream outfile;
		outfile << GetName() << ".significance";
		string outfilename = outfile.str();
		ofstream osoutfile(outfilename.c_str());

		for(int k=1;k<K;k++)	{
			ostringstream s;
			s << GetName() + "_" << k << ".pos";
			string name = s.str();
			ofstream ppos(name.c_str());
					
			int n70 = 0;
			int n90 = 0;

			for(int i=0;i<Nsite;i++)	{
				for(int j=0;j<Naa;j++)	{    
					ppos << pp[k][i][j] << '\t';

					if (pp[k][i][j] < 1-cutoff)	{
						cerr << k << '\t' << i << '\t' << AminoAcids[j] << '\t' << 1-pp[k][i][j] << '\n';
						osoutfile << k << '\t' << i << '\t' << AminoAcids[j] << '\t' << 1-pp[k][i][j] << '\n';
						n70++;
					}
					if (pp[k][i][j] > cutoff)	{
						cerr << k << '\t' << i << '\t' << AminoAcids[j] << '\t' << pp[k][i][j] << '\n';
						osoutfile << k << '\t' << i << '\t' << AminoAcids[j] << '\t' << pp[k][i][j] << '\n';
						n70++;
					}
					if(pp[k][i][j] > 0.90 || pp[k][i][j] < 0.10)	{
						n90++;
					}

				}
				ppos << '\n';
			}
			osoutfile << k << '\t' << n70 << '\t' << n90 << '\t' << Nsite*Naa << '\n';
			osoutfile << '\n';

		}
	}

	void PostPred(string basename)	{

		for (int i=0; i<size; i++)	{
			cerr << '.';
			// get next point -> will be stored into "model", and thus, will be accessible through GetModel()
			GetNextPoint();

			GetModel()->Update();

			ostringstream s;
			s << basename << "_" << i;

			GetModel()->PostPredAli(s.str());
		}

		cerr << "simulated data in " << basename << "_[rep].ali\n";
	}
};


int main(int argc, char* argv[])	{

	int burnin = 0;
	int every = 1;
	int until = -1;
	string name;

	double alpha = 0.1;
	double cutoff = 0.9;

	int ppred = 0;
	string basename = "";

	int fp = 0;
	int ncat = 10;

	int comp = 1;

	try	{

		if (argc == 1)	{
			throw(0);
		}

		int i = 1;
		while (i < argc)	{
			string s = argv[i];
			if (s == "-a")	{
				i++;
				alpha = atof(argv[i]);
			}
			else if (s == "-c")	{
				i++;
				cutoff = atof(argv[i]);
			}
			else if (s == "-fp")	{
				fp = 1;
			}
			else if (s == "-ncat")	{
				i++;
				ncat = atoi(argv[i]);
			}
			else if (s == "-ppred")	{
				ppred = 1;
				i++;
				basename = argv[i];
			}
			else if ( (s == "-x") || (s == "-extract") )	{
				i++;
				if (i == argc) throw(0);
				s = argv[i];
				burnin = atoi(argv[i]);
				i++;
				if (i == argc) throw(0);
				s = argv[i];
				every = atoi(argv[i]);
				i++;
				if (i == argc) throw(0);
				s = argv[i];
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
		if (name == "")	{
			throw(0);
		}
	}
	catch(...)	{
		cerr << "readpb [-x <burnin> <every> <until>] <chainname> \n";
		cerr << '\n';
		exit(1);
	}

	DirichletCodonUsageSelectionMSSample sample(name,burnin,every,until,alpha);
	if (ppred)	{
		cerr << "posterior predictive simulation\n";
		sample.PostPred(basename);
	}
	else if (fp)	{
		sample.ReadFalsePositives(ncat);
	}
	else if (comp)	{
		cerr << "reading, new version\n";
		sample.NewRead(cutoff);
	}
	else	{
		cerr << "read\n";
		sample.Read(cutoff);
	}

}


