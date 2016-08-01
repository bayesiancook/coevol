
#include "Sample.h"
#include "NormalOntologyModel.h"
#include "StringStreamUtils.h"

class NormalOntologySample : public Sample	{

	private:
	string modeltype;
	string datafile;
	string ontologyfile;
	int mingene;
	int minconcept;
	int withtoggle;
	double logvar;

	public:

	string GetModelType() {return modeltype;}

	NormalOntologyModel* GetModel() {return (NormalOntologyModel*) model;}

	NormalOntologySample(string filename, int inburnin, int inevery, int inuntil) : Sample(filename,inburnin,inevery,inuntil)	{
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

		withtoggle = 0;
		logvar = -100;
		// read model type, and other standard fields
		is >> modeltype;
		is >> datafile >> ontologyfile;
		is >> mingene >> minconcept;
		int check;
		is >> check;
		if (check)	{
			is >> withtoggle;
			is >> check;
			if (check)	{
				is >> logvar;
				is >> check;
				if (check)	{
					cerr << "error when reading model\n";
					exit(1);
				}
			}
		}

		is >> chainevery >> chainuntil >> chainsize;

		if (modeltype == "NORMALONTOLOGY")	{
			model = new NormalOntologyModel(datafile,ontologyfile,mingene,minconcept,withtoggle,logvar);
		}
		else	{
			cerr << "error when opening file "  << name << " : does not recognise model type : " << modeltype << '\n';
			exit(1);
		}
		model->FromStream(is);

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

	void ReadHisto(int ncat)	{

		int Nconcept = GetModel()->GetNconcept();
		int Ncont = GetModel()->GetNcont();

		double kappahisto[Ncont][ncat];
		double meanprop[Ncont][ncat];
		double meanlogvaralpha[Ncont];
		double varlogvaralpha[Ncont];
		for (int cont=0; cont<Ncont; cont++)	{
			for (int cat=0; cat<ncat; cat++)	{
				kappahisto[cont][cat] = 0;
				meanprop[cont][cat] = 0;
			}
			meanlogvaralpha[cont] = 0;
			varlogvaralpha[cont] = 0;
		}

		double log10 = log(10.0);
		double min = log(minjeff) / log10;
		double max = log(maxjeff) / log10;

		for (int i=0; i<size; i++)	{
			cerr << '.';

			GetNextPoint();

			for (int cont=0; cont<Ncont; cont++)	{
				double prop = GetModel()->toggle->GetMean(cont);
				double tmp = log(GetModel()->jeffkappa->GetVal(cont)->val()) / log10;
				int bin = ncat * (tmp - min) / (max- min);
				if (bin == ncat) bin--;
				kappahisto[cont][bin]++;
				meanprop[cont][bin] += prop;

				double logvaralpha = log(GetModel()->jeffalphavar->GetVal(cont)->val()) / log10;
				meanlogvaralpha[cont] += logvaralpha;
				varlogvaralpha[cont] += logvaralpha * logvaralpha;
			}

		}
		cerr << '\n';

		ofstream os((GetName() + ".histo").c_str());
		for (int cat=0; cat<ncat; cat++)	{
			os << min + ((double) (cat + 0.5)) / ncat * (max-min);
			for (int cont=0; cont<Ncont; cont++)	{
				if (kappahisto[cont][cat])	{
					meanprop[cont][cat] /= kappahisto[cont][cat];
				}
				kappahisto[cont][cat] /= size;
				os << '\t' << kappahisto[cont][cat] << '\t' << meanprop[cont][cat];
			}
			os << '\n';
		}

		cerr << '\n';
		cerr << "log variance of alpha : \n";
		for (int cont=0; cont<Ncont; cont++)	{
			meanlogvaralpha[cont] /= size;
			varlogvaralpha[cont] /= size;
			varlogvaralpha[cont] -= meanlogvaralpha[cont] * meanlogvaralpha[cont];
			cerr << cont << '\t' << meanlogvaralpha[cont] << '\t' << sqrt(varlogvaralpha[cont]) << '\n';
		}
		cerr << '\n';
		cerr << "output in " << GetName() << ".histo\n";
		cerr << '\n';
	}

	void Read()	{

		int Nconcept = GetModel()->GetNconcept();
		int Ncont = GetModel()->GetNcont();

		double meanbeta[Nconcept][Ncont];
		double varbeta[Nconcept][Ncont];
		double ppbeta[Nconcept][Ncont];

		double meanprob[Ncont];
		double varprob[Ncont];
		for (int cont=0; cont<Ncont; cont++)	{
			meanprob[cont] = 0;
			varprob[cont] = 0;
		}

		for (int g=0; g<Nconcept; g++)	{
			for (int c=0; c<Ncont; c++)	{
				meanbeta[g][c] = 0;
				varbeta[g][c] = 0;
				ppbeta[g][c] = 0;
			}
		}

		ofstream tos((GetName() + ".beta").c_str());

		for (int i=0; i<size; i++)	{
			cerr << '.';

			GetNextPoint();

			int nactive[Ncont];
			for (int cont=0; cont<Ncont; cont++)	{
				nactive[cont] = 0;
			}
			for (int g=0; g<Nconcept; g++)	{
				for (int c=0; c<Ncont; c++)	{
					double tmp = GetModel()->GetBeta(g,c);
					tos << tmp << '\t';
					meanbeta[g][c] += tmp;
					varbeta[g][c] += tmp * tmp;
					if (GetModel()->withToggle())	{
						if (tmp)	{
							ppbeta[g][c] ++;
							nactive[c]++;
						}
					}
					else	{
						if (tmp > 0)	{
							ppbeta[g][c] ++;
						}
					}
				}
				tos << '\n';
			}
			for (int cont=0; cont<Ncont; cont++)	{
				double tmp = ((double) nactive[cont]) / Nconcept;
				meanprob[cont] += tmp;
				varprob[cont] += tmp * tmp;
			}
		}

		cerr << '\n';
		cerr << "normalise\n";

		for (int g=0; g<Nconcept; g++)	{
			for (int c=0; c<Ncont; c++)	{
				meanbeta[g][c] /= size;
				varbeta[g][c] /= size;
				varbeta[g][c] -= meanbeta[g][c] * meanbeta[g][c];
				ppbeta[g][c] /= size;
			}
		}

		ofstream os((GetName() + ".postmeanbeta").c_str());
		for (int g=0; g<Nconcept; g++)	{
			for (int c=0; c<Ncont; c++)	{
				os << GetModel()->GetOntology()->GetConceptName(g) << '\t' << meanbeta[g][c] << '\t' << varbeta[g][c] << '\t' << ppbeta[g][c] << '\n';
			}
		}

		if (GetModel()->withToggle())	{
			cerr << '\n';
			cerr << "proportion of active concepts\n";
			for (int cont=0; cont<Ncont; cont++)	{
				meanprob[cont] /= size;
				varprob[cont] /= size;
				varprob[cont] -= meanprob[cont] * meanprob[cont];
				cerr << cont << '\t' << meanprob[cont] << " +/- " << sqrt(varprob[cont]) << '\n';
			}
			cerr << '\n';
		}
		cerr << '\n';
		cerr << "output in " << GetName() << ".postmeanbeta\n";
		cerr << '\n';
	}
};


int main(int argc, char* argv[])	{

	int burnin = 0;
	int every = 1;
	int until = -1;

	int ncat = 20;
	int histo = 0;

	string name;

	try	{

		if (argc == 1)	{
			throw(0);
		}

		int i = 1;
		while (i < argc)	{
			string s = argv[i];

			if (s == "-histo")	{
				histo = 1;
			}
			else if (s == "-ncat")	{
				i++;
				ncat = atoi(argv[i]);
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
		cerr << "readnormalontology [-x <burnin> <every> <until>] <chainname> \n";
		cerr << '\n';
		exit(1);
	}

	NormalOntologySample sample(name,burnin,every,until);

	if (histo)	{
		sample.ReadHisto(ncat);
	}
	else	{
		sample.Read();
	}
}



