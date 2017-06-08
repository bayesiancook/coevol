
#include "Sample.h"
#include "DirichletNormalCodonUsageSelectionModelMS.h"
#include "GTRSubMatrix.h"
#include "CI.h"

class DirichletNormalCodonUsageSelectionMSSample : public Sample	{

	private:
	string modeltype;
	string datafile;
	string treefile;
	int category;
	int conjugate;
	string type;
	string mechanism;
	double alpha;
	double nucmutationrate[4][4];
	double codonusage[61];
	double cupparray[61][61];
	int pparray[4][20];

	public:

	string GetModelType() {return modeltype;}

	DirichletNormalCodonUsageSelectionModelMS* GetModel() {return (DirichletNormalCodonUsageSelectionModelMS*) model;}

	DirichletNormalCodonUsageSelectionMSSample(string filename, int inburnin, int inevery, int inuntil, double inalpha = 0.1) : Sample(filename,inburnin,inevery,inuntil)	{
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
		is >> datafile >> treefile >> category;
		is >> conjugate;
		is >> type >> mechanism;
		is >> chainevery >> chainuntil >> chainsize;
		// the chain's saving frequency, upper limit and current size
		// not to be confused with the sample's subsampling frequency, upper limit and size


		// make a new model depending on the type obtained from the file
		if (modeltype == "SELECTIONGTR")	{
			cerr << "CREATE\n";
			model = new DirichletNormalCodonUsageSelectionModelMS(datafile,treefile,category,type,conjugate,mechanism,false);
		}
		else	{
			cerr << "error when opening file "  << name << '\n';
			exit(1);
		}

		// read model (i.e. chain's last point) from <name>.param
		cerr << "from stream\n";
		model->FromStream(is);
		// cerr << "UPDATE\n";
		// model->Update();
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
						logsum += log(GetModel()->GetSelectionProfile(k, i, j));
					}
					for(int j=0;j<Naa;j++)	{  
						double tmp = (log(GetModel()->GetSelectionProfile(k, i, j))) - (1.0/20) * (logsum);
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
				cout << 100 * falsepositives[k][cat] << '\n';
			}
			cout << '\n';
		}
	}

	void ReadFDR(string truefile, int ncat)	{

		int nlevel = 4;
		double levelval[nlevel];
		levelval[0] = 1.0;
		levelval[1] = 2.0;
		levelval[2] = 3.0;
		levelval[3] = 4.0;

		int Nsite = GetModel()->GetSite();
		int Naa = GetModel()->GetaaState();
		int K = GetModel()->GetCategory();

		// the true selection profiles (that have been used to make the posterior predictive simulation)
		double*** truesel = new double**[K];

		// the mean selection profiles that will be estimated
		double*** meansel = new double**[K];
		list<double>*** listsel = new list<double>**[K];

		// the posterior probabilities associated to our estimation
		double*** ppsel = new double**[K];

		// make the array and get the true profiles from file
		ifstream is((truefile + ".trueprofiles").c_str());

		int kmin = 1;

		for (int k=kmin; k<K; k++)	{
			int tmpk;
			is >> tmpk;
			if (tmpk != k)	{
				cerr << "error when reading true profiles\n";
				exit(1);
			}
			truesel[k] = new double*[Nsite];
			meansel[k] = new double*[Nsite];
			listsel[k] = new list<double>*[Nsite];
			ppsel[k] = new double*[Nsite];
			for (int i=0; i<Nsite; i++)	{
				truesel[k][i] = new double[Naa];
				meansel[k][i] = new double[Naa];
				listsel[k][i] = new list<double>[Naa];
				ppsel[k][i] = new double[Naa];
				int tmpi;
				is >> tmpi;
				for (int j=0; j<Naa; j++)	{
					is >> truesel[k][i][j];
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
						logsum += log(GetModel()->GetSelectionProfile(k, i, j));
					}
					for(int j=0;j<Naa;j++)	{  
						double tmp = (log(GetModel()->GetSelectionProfile(k, i, j))) - (1.0/20) * (logsum);
						meansel[k][i][j] += tmp;
						listsel[k][i][j].push_back(tmp);
						if(tmp > 0)	{
							ppsel[k][i][j]++;
						}
					}
				}
			}
		}
		cerr << '\n';

		// for each differential selection category (k=1..K-1)
		// and for each possible cutoff c between 0.5 and 1
		// calculate the number of discoveries
		// (defined as those site/aa for which  pp < c or pp > 1 - c)
		// and the number of those discoveries that are false

		double* error = new double[K];
		double* coverage = new double[K];
		for(int k=kmin;k<K;k++)	{
			coverage[k] = 0;
			error[k] = 0;
		}

		double** discoveries = new double*[K];
		double** falsediscoveries = new double*[K];
		double*** ddiscoveries = new double**[K];
		double*** ffalsediscoveries = new double**[K];
		for(int k=kmin;k<K;k++)	{
			discoveries[k] = new double[ncat];
			falsediscoveries[k] = new double[ncat];
			ddiscoveries[k] = new double*[ncat];
			ffalsediscoveries[k] = new double*[ncat];
			for (int cat=0; cat<ncat; cat++)	{
				discoveries[k][cat] = 0;
				falsediscoveries[k][cat] = 0;
				ddiscoveries[k][cat] = new double[nlevel];
				ffalsediscoveries[k][cat] = new double[nlevel];
				for (int level=0; level<nlevel; level++)	{
					ddiscoveries[k][cat][level] = 0;
					ffalsediscoveries[k][cat][level] = 0;
				}
			}
		}

		double** levelcount = new double*[K];
		for(int k=kmin;k<K;k++)	{
			levelcount[k] = new double[nlevel];
			for (int level=0; level<nlevel; level++)	{
				levelcount[k][level] = 0;
			}
		}

		for(int k=kmin;k<K;k++)	{
			for(int i=0;i<Nsite;i++)	{
				for(int j=0;j<Naa;j++)	{ 
					meansel[k][i][j] /= size;
					ppsel[k][i][j] /= size;


					error[k] += (meansel[k][i][j] - truesel[k][i][j]) * (meansel[k][i][j] - truesel[k][i][j]);
					double median, min, max;
					GetCI(listsel[k][i][j],0.975,median,min,max);
					if ((truesel[k][i][j] > min) && (truesel[k][i][j] < max))	{
						coverage[k]++;
					}
					// get the smallest cutoff at which this pp will be considered as a discovery
					// thus, for instance, if ncat = 10
					// 0.45 < pp < 0.55: cat = 0
					// 0.40 < pp < 0.45 or 0.55 < pp < 0.60: cat = 1
					// ..
					// pp < 0.05 or pp > 0.95: cat = 9
					int cat = (int) (2*ncat*(ppsel[k][i][j] - 0.5));
					if (cat < 0)	{
						cat = -cat;
					}
					if (cat == ncat)	{
						cat--;
					}

					// increment the number of discoveries
					discoveries[k][cat]++;

					int level = nlevel-1;
					while (level && (fabs(truesel[k][i][j]) < levelval[level]))	{
						level--;
					}
					levelcount[k][level]++;
					ddiscoveries[k][cat][level]++;

					// this discovery will be a false discovery if	
					// either pp > 0.5 (we infer a positive effect) but in fact the true selection effect is negavite
					// or pp < 0.5 (we infer a negative effect) but in fact the true effect is positive
					if ((ppsel[k][i][j] - 0.5) *truesel[k][i][j] < 0)	{
						falsediscoveries[k][cat]++;
						ffalsediscoveries[k][cat][level]++;
					}
				}
			}

			// make an output file in which the mean selection effects are tabulated against the true selection effects
			// they can then be compared by plotting true versus estimated
			ostringstream s;
			s << GetName() << "_" << k << ".compsel";
			ofstream os(s.str().c_str());
			for(int i=0;i<Nsite;i++)	{
				for(int j=0;j<Naa;j++)	{ 
					os << meansel[k][i][j] << '\t' << truesel[k][i][j] << '\n';
				}
			}

			error[k] /= Nsite*Naa;
			coverage[k] /= Nsite*Naa;
		}

		// make the cumulated numbers of discoveries and false discoveries
		// (i.e.: if you decide for a threshold of 0.80, then, you should sum all discoveries for all categories corrsponding to pp > 0.8)
		for(int k=kmin;k<K;k++)	{
			for (int cat=0; cat<ncat; cat++)	{
				for (int cat2=cat+1; cat2<ncat; cat2++)	{
					discoveries[k][cat] += discoveries[k][cat2];
					falsediscoveries[k][cat] += falsediscoveries[k][cat2];
					for (int level=0; level<nlevel; level++)	{
						ddiscoveries[k][cat][level] += ddiscoveries[k][cat2][level];
						ffalsediscoveries[k][cat][level] += ffalsediscoveries[k][cat2][level];
					}
				}
			}
		}

		// error and coverage
		ofstream cos((GetName() + ".coverage").c_str());
		cos << "cat\terror\tcoverage\n";
		for(int k=kmin;k<K;k++)	{
			cos << k << '\t' << sqrt(error[k]) << '\t' << coverage[k] << '\n';
		}
		cos.close();

		// output estimated numbers of (false) disoveries and FDR
		ofstream cout((GetName() + ".fdr").c_str());
		cout.precision(4);
		cout << '\n';
		cout << "false discovery rates:\n";
		cout << '\n';
		for(int k=kmin;k<K;k++)	{
			cout << "category : " << k << '\n';
			cout << "c\t#disc\t#false\tFDR\t";
			for (int level=0; level<nlevel; level++)	{
				cout << ">" << levelval[level] << '\t' << '\t';
			}
			cout << '\n';
			cout << '\n';
			for (int cat=0; cat<ncat; cat++)	{
				cout << 0.5 + ((double) cat) / 2 / ncat << '\t';
				cout << discoveries[k][cat] << '\t' << falsediscoveries[k][cat] << '\t';
				if (discoveries[k][cat])	{
					cout << 100 * falsediscoveries[k][cat] / discoveries[k][cat] << '\t';
				}
				else	{
					cout << 0 << '\t';
				}

				for (int level=0; level<nlevel; level++)	{
					if (levelcount[k][level])	{
						cout << 100 * ddiscoveries[k][cat][level] / levelcount[k][level] << '\t';
					}
					else	{
						cout << 0 << '\t';
					}
					if (ddiscoveries[k][cat][level])	{
						cout << 100 * ffalsediscoveries[k][cat][level] / ddiscoveries[k][cat][level] << '\t';
					}
					else	{
						cout << 0 << '\t';
					}
				}
				cout << '\n';
			}
			cout << '\n';
		}

		/*
		for(int k=kmin;k<K;k++)	{
			cout << "category : " << k << '\n';
			cout << "c\t#disc\t#false\tFDR\n";
			for (int cat=0; cat<ncat; cat++)	{
				cout << 0.5 + ((double) cat) / 2 / ncat << '\t';
				cout << discoveries[k][cat] << '\t' << falsediscoveries[k][cat] << '\t';
				if (discoveries[k][cat])	{
					cout << 100 * falsediscoveries[k][cat] / discoveries[k][cat] << '\n';
				}
				else	{
					cout << 0 << '\n';
				}
			}
			cout << '\n';
		}
		*/
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
					const double* delta = GetModel()->GetDelta(k,i);

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
					if (ppsel[k][i][a] > cutoff)	{
						os2 << i << '\t' << a << '\t' << ppsel[k][i][a] << '\t' << meansel[k][i][a] << '\n';
					}
				}
			}
		}
	}

	void Read(double cutoff)	{

		// prepare the mean and variance
		
		int Nsite = GetModel()->GetSite();
		int Naastate = GetModel()->GetaaState();
		int Ncodonstate = GetModel()->GetCodonState();

		int K = GetModel()->GetCategory();

		list<double> listaveselectionarray[K][Nsite][Naastate];
		double aveselectionarray[K][Nsite][Naastate];
		double stationary[K][Nsite][Naastate];
		double possel[K][Nsite][Naastate];
		double multiplicativeselectionarray[K][Nsite][Naastate];
		double differentialfitness[K-1];


		for(int k=0;k<K;k++)	{
			for(int i=0;i<Nsite;i++)	{
				for(int j=0;j<Naastate;j++)	{
					aveselectionarray[k][i][j] = 0;
					multiplicativeselectionarray[k][i][j] = 0;
					possel[k][i][j] = 0;
					stationary[k][i][j] = 0;
				}
			}
		}


		for(int i=0;i<4;i++)	{
			for(int j=0;j<4;j++)	{
				nucmutationrate[i][j] = 0;
			}
		}


		for(int i=0;i<Ncodonstate;i++)	{
			codonusage[i] = 0;
		}

		for(int i=0;i<Ncodonstate;i++)	{
			for(int j=0;j<Ncodonstate;j++)	{
				cupparray[i][j] = 0;
			}
		}

		for(int i=0;i<K;i++)	{
			for(int j=0;j<20;j++)	{
				pparray[i][j] = 0;
			}
		}


		for(int i=0;i<K-1;i++)	{
			differentialfitness[i] = 0;
		}


		ostringstream s7;
		s7 << GetName() <<  ".B57";
		string name7 = s7.str();
		ofstream os7(name7.c_str());

		ostringstream s8;
		s8 << GetName() <<  ".B35";
		string name8 = s8.str();
		ofstream os8(name8.c_str());

		ostringstream s9;
		s9 << GetName() +  ".deltacondition";
		string name9 = s9.str();
		ofstream os9(name9.c_str());

		// cycle over the sample
		for (int c=0; c<size; c++)	{
			GetNextPoint();
			cerr << '.';

			//get selection profile
			for(int k=0;k<K;k++)	{
				ostringstream s10;
				s10 << GetName() + "_" << k << ".selectionfrq";
				string name10 = s10.str();
				ofstream os10(name10.c_str(), ios_base::app);

				for(int i=0;i<Nsite;i++)	{
					double logsum = 0; 
					for(int j=0;j<Naastate;j++)	{ 
						logsum += log(GetModel()->GetSelectionProfile(k, i, j));
					}
					for(int j=0;j<Naastate;j++)	{  
						double tmp = 0;
						if (k != 0)	{
							tmp = (log(GetModel()->GetSelectionProfile(k, i, j))) - (1.0/20) * (logsum);
						}

						else	{
							tmp = GetModel()->GetSelectionProfile(k,i,j);
						}
						aveselectionarray[k][i][j] += tmp;
						listaveselectionarray[k][i][j].push_back(tmp);

						if((((log(GetModel()->GetSelectionProfile(k, i, j))) - (1.0/20) * logsum)) > 0)	{
							possel[k][i][j]++;
						}

						os10 << log(GetModel()->GetSelectionProfile(k, i, j))- (1.0/20) * (logsum) << '\n';
					}
				}
			}
      

			//get stationary
			for(int k=0;k<K;k++)	{
				for(int i=0;i<Nsite;i++)	{
					for(int j=0;j<Naastate;j++)	{    
						stationary[k][i][j] += GetModel()->GetStationary(k, i, j);
					}
				}
			}


			//get nucmatrix which is the mutation rate 
			GetModel()->Getnucmatrix()->localUpdate();

			for(int i=0;i<4;i++)	{
				for(int j=0;j<4;j++)	{
					nucmutationrate[i][j] += (*GetModel()->Getnucmatrix())(i,j);
				}
			}


			//get codon usage posterior probability
			for(int i=0;i<Ncodonstate;i++)	{
				codonusage[i] += GetModel()->GetCodonUsage(i);
				for(int j=i;j<Ncodonstate;j++)	{
					if(GetModel()->GetCodonUsage(i) > GetModel()->GetCodonUsage(j))	{
						cupparray[i][j]++;
					}
				}
			}

			// get true fitness difference
			for(int k=0;k<K-1;k++)	{
				double total = 0;
				for(int i=0;i<Nsite;i++)	{
					GetModel()->Getsubmatrix(k,i)->UpdateMatrix();
					for(int j=0;j<Ncodonstate;j++)	{    
						total += ((GetModel()->Getsubmatrix(k+1,i)->Stationary(j)) - (GetModel()->Getsubmatrix(k,i)->Stationary(j))) * log(GetModel()->GetSelectionProfile(k+1, i, (GetModel()->GetCodonStateSpace()->Translation(j))));
					}
				}
				differentialfitness[k] += total;
				os9 << total << '\t';
			}
			os9 << '\n';
		}
		cerr << '\n';

		ostringstream outfile;
		outfile << GetName() << ".significance";
		string outfilename = outfile.str();
		ofstream osoutfile(outfilename.c_str());


		for(int k=0;k<K;k++)	{

			ostringstream s1;
			s1 << GetName() + "_" << k << ".aveselection";
			string name1 = s1.str();
			ofstream os1(name1.c_str());
					
			ostringstream s14;
			s14 << GetName() + "_" << k << ".multiselection";
			string name14 = s14.str();
			ofstream os14(name14.c_str());
					

			ostringstream s5;
			s5 << GetName() + "_" << k << ".pos";
			string name5 = s5.str();
			ofstream os5(name5.c_str());

			// os1 << Nsite << '\n';
			os14 << Nsite << '\n';

			int n70=0;
			int n90=0;

			for(int i=0;i<Nsite;i++)	{
				for(int j=0;j<Naastate;j++)	{    
					aveselectionarray[k][i][j] /= size;
					os1 << aveselectionarray[k][i][j] << '\t';

					if(k ==1)	{
						multiplicativeselectionarray[k][i][j] = aveselectionarray[k][i][j] * aveselectionarray[0][i][j];
					}
					if(k ==2)	{
						multiplicativeselectionarray[k][i][j] = aveselectionarray[k][i][j] * aveselectionarray[k-1][i][j] * aveselectionarray[0][i][j];
					}
					if(k ==3)	{
						multiplicativeselectionarray[k][i][j] = aveselectionarray[k][i][j] * aveselectionarray[k-2][i][j] * aveselectionarray[0][i][j];
					}

					os14 << multiplicativeselectionarray[k][i][j] << '\t';

					possel[k][i][j] /= size;
					os5 << possel[k][i][j] << '\t';
//cerr << "pos " << k << "," << i << "," << j << ": " << possel[k][i][j] << '\n';

					if((possel[k][i][j] > 0) &  (possel[k][i][j] <= 0.05))   {
						pparray[k][0]++;
					}
					if((possel[k][i][j] > 0.05) &  (possel[k][i][j] <= 0.1))   {
						pparray[k][1]++;
					}
					if((possel[k][i][j] > 0.1) &  (possel[k][i][j] <= 0.15))   {
						pparray[k][2]++;
					}
					if((possel[k][i][j] > 0.15) &  (possel[k][i][j] <= 0.2))   {
						pparray[k][3]++;
					}
					if((possel[k][i][j] > 0.2) &  (possel[k][i][j] <= 0.25))   {
						pparray[k][4]++;
					}
					if((possel[k][i][j] > 0.25) &  (possel[k][i][j] <= 0.3))   {
						pparray[k][5]++;
					}
					if((possel[k][i][j] > 0.3) &  (possel[k][i][j] <= 0.35))   {
						pparray[k][6]++;
					}
					if((possel[k][i][j] > 0.35) &  (possel[k][i][j] <= 0.4))   {
						pparray[k][7]++;
					}
					if((possel[k][i][j] > 0.4) &  (possel[k][i][j] <= 0.45))   {
						pparray[k][8]++;
					}
					if((possel[k][i][j] > 0.45) &  (possel[k][i][j] <= 0.5))   {
						pparray[k][9]++;
					}
					if((possel[k][i][j] > 0.5) &  (possel[k][i][j] <= 0.55))   {
						pparray[k][10]++;
					}
					if((possel[k][i][j] > 0.55) &  (possel[k][i][j] <= 0.6))   {
						pparray[k][11]++;
					}
					if((possel[k][i][j] > 0.6) &  (possel[k][i][j] <= 0.65))   {
						pparray[k][12]++;
					}
					if((possel[k][i][j] > 0.65) &  (possel[k][i][j] <= 0.7))   {
						pparray[k][13]++;
					}
					if((possel[k][i][j] > 0.7) &  (possel[k][i][j] <= 0.75))   {
						pparray[k][14]++;
					}
					if((possel[k][i][j] > 0.75) &  (possel[k][i][j] <= 0.8))   {
						pparray[k][15]++;
					}
					if((possel[k][i][j] > 0.8) &  (possel[k][i][j] <= 0.85))   {
						pparray[k][16]++;
					}
					if((possel[k][i][j] > 0.85) &  (possel[k][i][j] <= 0.9))   {
						pparray[k][17]++;
					}
					if((possel[k][i][j] > 0.9) &  (possel[k][i][j] <= 0.95))   {
						pparray[k][18]++;
					}
					if((possel[k][i][j] > 0.95) &  (possel[k][i][j] <= 1.0))   {
						pparray[k][19]++;
					}
  
					if((possel[k][i][j] > cutoff) || (possel[k][i][j] < 1 - cutoff))	{
						if (k > 1)	{
							double median, min, max;
							GetCI(listaveselectionarray[k][i][j],0.975,median,min,max);
							cerr << k << '\t' << i << '\t' << AminoAcids[j] << '\t' << possel[k][i][j] << '\t' << median << '\t' << min << '\t' << max << '\n';
							osoutfile << k << '\t' << i << '\t' << AminoAcids[j] << '\t' << possel[k][i][j] << '\t' << median << '\t' << min << '\t' << max << '\n';
							n70++;
						}
					}
					if(possel[k][i][j] > 0.90 || possel[k][i][j] < 0.10)	{
						if (k > 1)	{
							n90++;
						}
					}

				}
				os1 << '\n';
				os5 << '\n';
				os14 << '\n';
			}

			if (k)	{
				cerr << '\n';
				cerr << k << '\t' << n70 << '\t' << n90 << '\t' << Nsite*Naastate << '\n';

				osoutfile << '\n';
				osoutfile << k << '\t' << n70 << '\t' << n90 << '\t' << Nsite*Naastate << '\n';
			}

		}

		for(int k=0;k<K;k++)	{
			ostringstream t0;
			t0 << GetName() + "_" << k << ".pp";
			string namet0 = t0.str();
			ofstream out(namet0.c_str());

			for(int l=0;l<20;l++)	{
				out << pparray[k][l] << '\n';
			}
		}


		ostringstream s2;
		s2 << GetName() +  ".mutationrate";
		string name2 = s2.str();
		ofstream os2(name2.c_str());

		for(int i=0;i<4;i++)	{
			for(int j=0;j<4;j++)	{
				nucmutationrate[i][j] /= size;
				os2 << nucmutationrate[i][j] << '\t';
			}
			os2 << '\n';
		}


		ostringstream s3;
		s3 << GetName() +  ".codonusage";
		string name3 = s3.str();
		ofstream os3(name3.c_str());

		for(int i=0;i<Ncodonstate;i++)	{
			codonusage[i] /= size;
			os3 << codonusage[i] << '\t';
		}



		ostringstream s6;
		s6 << GetName() +  ".cupp";
		string name6 = s6.str();
		ofstream os6(name6.c_str());

		for(int i=0;i<Ncodonstate;i++)	{
			for(int j=0;j<Ncodonstate;j++)	{
				cupparray[i][j] /= size;
				os6 << cupparray[i][j] << '\t';
			}
			os6 << '\n';
		}

		ofstream os4((GetName() + ".fitness").c_str());
		os4 << Nsite << '\n';
		for(int k=0;k<K;k++)	{
			for(int i=0;i<Nsite;i++)	{
				for(int j=0;j<Naastate;j++)	{    
					stationary[k][i][j] /= size;
					os4 << stationary[k][i][j] << '\t';
				}
				os4 << '\n';
			}
		}

	}


	
	void PostPred(string basename, int sitemin, int sitemax)	{

		for (int i=0; i<size; i++)	{
			cerr << '.';
			// get next point -> will be stored into "model", and thus, will be accessible through GetModel()
			GetNextPoint();

			GetModel()->Update();

			ostringstream s;
			s << basename << "_" << i;

			GetModel()->OutputSelectionProfile(s.str(),sitemin,sitemax);
			GetModel()->PostPredAli(s.str(),sitemin,sitemax);
		}

		cerr << "simulated data in " << basename << "_[rep].ali\n";
	}
	
};


int main(int argc, char* argv[])	{

	int burnin = 0;
	int every = 1;
	int until = -1;
	string name;

	int ppred = 0;

	int sitemin = -1;
	int sitemax = -1;

	double alpha = 0.1;
	double cutoff = 0.7;
	string basename = "";

	int fdr = 0;
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
			else if (s == "-old")	{
				comp = 0;
			}
			else if (s == "-ppred")	{
				ppred = 1;
				i++;
				basename = argv[i];
			}
			else if (s == "-sub")	{
				i++;
				sitemin = atoi(argv[i]);
				i++;
				sitemax = atoi(argv[i]);
			}
			else if (s == "-fdr")	{
				fdr = 1;
				i++;
				basename = argv[i];
			}
			else if (s == "-fp")	{
				fp = 1;
			}
			else if (s == "-ncat")	{
				i++;
				ncat = atoi(argv[i]);
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


	DirichletNormalCodonUsageSelectionMSSample sample(name,burnin,every,until,alpha);
	
	if (ppred)	{
		cerr << "posterior predictive simulation\n";
		sample.PostPred(basename,sitemin,sitemax);
	}
	else if (fdr)	{
		sample.ReadFDR(basename,ncat);
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

