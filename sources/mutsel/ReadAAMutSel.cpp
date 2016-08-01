
#include "Sample.h"
#include "DirAAMutSelMatMixModel.h"
// #include "AAMutSelMatMixModel.h"
#include "MeanValTree.h"
#include "ProteinSequenceAlignment.h"

class AAMutSelSample : public Sample	{

	private:
	string modeltype;
	string datafile;
	string treefile;
	int P;
	GeneticCodeType type;

	public:

	string GetModelType() {return modeltype;}

	AAMutSelModel* GetModel() {return (AAMutSelModel*) model;}

	AAMutSelSample(string filename, int inburnin, int inevery, int inuntil) : Sample(filename,inburnin,inevery,inuntil)	{
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
		is >> type;
		is >> datafile >> treefile;
		is >> P;
		is >> chainevery >> chainuntil >> chainsize;
		// the chain's saving frequency, upper limit and current size
		// not to be confused with the sample's subsampling frequency, upper limit and size

		// make a new model depending on the type obtained from the file
		if (modeltype == "AAMUTSEL")	{
			model = new AAMutSelModel(datafile,treefile,P,true,type);
		}
		else	{
			cerr << "error when opening file "  << name << " : does not recognise model type : " << modeltype << '\n';
			exit(1);
		}

		// read model (i.e. chain's last point) from <name>.param
		model->FromStream(is);
		cerr << "UPDATE\n";
		// model->Update();

		// open <name>.chain, and prepare stream and stream iterator
		OpenChainFile();
		// now, size is defined (it is the total number of points with which this Sample object will make all its various posterior averages)
		// all these points can be accessed to (only once) by repeated calls to GetNextPoint()
	}

	// a very simple (and quite uninteresting) method for obtaining
	// the posterior mean and variance of the total length of the tree
	void Read()	{

		int Nsite = GetModel()->GetNsite();
		
		// compute mean fitness profiles
		double** fit = new double*[Nsite];
		for (int i=0; i<Nsite; i++)	{
			fit[i] = new double[Naa];
			for (int j=0; j<Naa; j++)	{
				fit[i][j] = 0;
			}
		}

		// and mean aa eq freq profiles
		/*
		double** stat = new double*[Nsite];
		for (int i=0; i<Nsite; i++)	{
			stat[i] = new double[Naa];
			for (int j=0; j<Naa; j++)	{
				stat[i][j] = 0;
			}
		}
		*/

		// to be compared with empirical profiles
		double** emp = new double*[Nsite];
		for (int i=0; i<Nsite; i++)	{
			emp[i] = new double[Naa];
		}
		ProteinSequenceAlignment* protali = new ProteinSequenceAlignment(GetModel()->GetCodonData());
		protali->GetSiteEmpiricalFreq(emp);

		ofstream empos((name + ".ssempfreq").c_str());
		empos << Nsite << '\n';
		for (int i=0; i<Nsite; i++)	{
			for (int j=0; j<Naa; j++)	{
				empos << emp[i][j] << '\t';
			}
			empos << '\n';
		}

		for (int i=0; i<size; i++)	{
			cerr << '.';
			GetNextPoint();

			for (int i=0; i<Nsite; i++)	{
				for (int j=0; j<Naa; j++)	{
					fit[i][j] += GetModel()->GetSiteStat(i)[j];
				}
			}
		}
		cerr << '\n';

		ofstream fitos((name + ".ssfitness").c_str());
		fitos << Nsite << '\n';
		for (int i=0; i<Nsite; i++)	{
			for (int j=0; j<Naa; j++)	{
				fit[i][j] /= size;
				fitos << fit[i][j] << '\t';
			}
			fitos << '\n';
		}
				
	}

	void PostPredCompositionalHeterogeneity()	{

		ofstream compos((name + ".comp").c_str());
		compos << size << '\t' << GetModel()->GetNtaxa() << '\t' << GetModel()->GetNstate() << '\n';
		compos << '\n';

		ofstream compaaos((name + ".compaa").c_str());
		compos << size << '\t' << GetModel()->GetNtaxa() << '\t' << Naa << '\n';
		compos << '\n';

		ofstream compnucos((name + ".compnuc").c_str());
		compos << size << '\t' << GetModel()->GetNtaxa() << '\t' << Nnuc << '\n';
		compos << '\n';

		double obs = GetModel()->GetCodonData()->CompositionalHeterogeneity(&compos);
		double obsaa = GetModel()->GetCodonData()->AminoAcidCompositionalHeterogeneity(&compaaos);
		double obsnuc = GetModel()->GetCodonData()->NucleotideCompositionalHeterogeneity(&compnucos);
		double obsnucdiff = GetModel()->GetCodonData()->Nucleotide123CompositionalHeterogeneity();
		double* obsnucpos = new double[3];
		for (int pos=0; pos<3; pos++) {
			obsnucpos[pos] = GetModel()->GetCodonData()->NucleotideCompositionalHeterogeneity(0,pos);
		}


		double mean = 0;
		double var = 0;
		double pp = 0;
		double meanaa = 0;
		double varaa = 0;
		double ppaa = 0;
		double meannuc = 0;
		double varnuc = 0;
		double ppnuc = 0;
		double meannucdiff = 0;
		double varnucdiff = 0;
		double ppnucdiff = 0;
		double* meannucpos = new double[3];
		double* varnucpos = new double[3];
		double* ppnucpos = new double[3];
		for (int pos=0; pos<3; pos++) {
			meannucpos[pos] = 0;
			varnucpos[pos] = 0;
			ppnucpos[pos] = 0;
		}
		
		// cycle over the sample
		cerr << "total number of points : " << size << '\n';
		for (int i=0; i<size; i++)	{
			cerr << '.';

			// get next point -> from now on, accessible through GetModel()
			GetNextPoint();
			
			// update point (to be sure that stochastic mapping will work based on correct model specs)
			GetModel()->Update();
			GetModel()->PostPredSample();

			double tmp = GetModel()->GetPostPredCodonData()->CompositionalHeterogeneity(&compos);
			mean += tmp;
			var += tmp * tmp;
			if (obs < tmp)	{
				pp++;
			}

			tmp = GetModel()->GetPostPredCodonData()->AminoAcidCompositionalHeterogeneity(&compaaos);
			meanaa += tmp;
			varaa += tmp * tmp;
			if (obsaa < tmp)	{
				ppaa++;

			}
			tmp = GetModel()->GetPostPredCodonData()->NucleotideCompositionalHeterogeneity(&compnucos);
			meannuc += tmp;
			varnuc += tmp * tmp;
			if (obsnuc < tmp)	{
				ppnuc++;
			}

			tmp = GetModel()->GetPostPredCodonData()->Nucleotide123CompositionalHeterogeneity();
			meannucdiff += tmp;
			varnucdiff += tmp * tmp;
			if (obsnucdiff < tmp)	{
				ppnucdiff++;
			}

			for (int pos=0; pos<3; pos++) {
				tmp = GetModel()->GetPostPredCodonData()->NucleotideCompositionalHeterogeneity(0,pos);
				meannucpos[pos] += tmp;
				varnucpos[pos] += tmp * tmp;
				if (obsnucpos[pos] < tmp)	{
					ppnucpos[pos] ++;
				}
			}
		}
		mean /= size;
		var /= size;
		var -= mean * mean;
		pp /= size;

		meanaa /= size;
		varaa /= size;
		varaa -= meanaa * meanaa;
		ppaa /= size;

		meannuc /= size;
		varnuc /= size;
		varnuc -= meannuc * meannuc;
		ppnuc /= size;

		meannucdiff /= size;
		varnucdiff /= size;
		varnucdiff -= meannucdiff * meannucdiff;
		ppnucdiff /= size;

		for (int pos=0; pos<3; pos++) {
			meannucpos[pos] /= size;
			varnucpos[pos] /= size;
			varnucpos[pos] -= meannucpos[pos] * meannucpos[pos];
			ppnucpos[pos] /= size;
		}

		ofstream os((name + ".ppsummary").c_str());
		os << '\n';
		os << "codon\n";
		os << "observed  : " << obs << '\n';
		os << "predictive: " << mean << '\t' << sqrt(var) << '\n';
		os << "pp        : " << pp << '\n';
		os << "z-score   : " << (obs - mean) / sqrt(var) << '\n';
		os << '\n';
		os << "aa\n";
		os << "observed  : " << obsaa << '\n';
		os << "predictive: " << meanaa << '\t' << sqrt(varaa) << '\n';
		os << "pp        : " << ppaa << '\n';
		os << "z-score   : " << (obsaa - meanaa) / sqrt(varaa) << '\n';
		os << '\n';
		os << "nuc\n";
		os << "observed  : " << obsnuc << '\n';
		os << "predictive: " << meannuc << '\t' << sqrt(varnuc) << '\n';
		os << "pp        : " << ppnuc << '\n';
		os << "z-score   : " << (obsnuc - meannuc) / sqrt(varnuc) << '\n';
		os << '\n';
		os << "nuc diff between (12) and (3)\n";
		os << "observed  : " << obsnucdiff << '\n';
		os << "predictive: " << meannucdiff << '\t' << sqrt(varnucdiff) << '\n';
		os << "pp        : " << ppnucdiff << '\n';
		os << "z-score   : " << (obsnucdiff - meannucdiff) / sqrt(varnucdiff) << '\n';
		os << '\n';
		
		for (int pos=0; pos<3; pos++) {
			os << "nuc pos " << pos << '\n';
			os << "observed  : " << obsnucpos[pos] << '\n';
			os << "predictive: " << meannucpos[pos] << '\t' << sqrt(varnucpos[pos]) << '\n';
			os << "pp        : " << ppnucpos[pos] << '\n';
			os << "z-score   : " << (obsnucpos[pos] - meannucpos[pos]) / sqrt(varnucpos[pos]) << '\n';
			os << '\n';
		}

		delete[] obsnucpos;
		delete[] meannucpos;
		delete[] varnucpos;
		delete[] ppnucpos;
	}
};


int main(int argc, char* argv[])	{

	int burnin = 0;
	int every = 1;
	int until = -1;
	string name;

	int postpred = 0;

	try	{

		if (argc == 1)	{
			throw(0);
		}

		int i = 1;
		while (i < argc)	{
			string s = argv[i];

			if ( (s == "-x") || (s == "-extract") )	{
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
			else if (s == "-pp")	{
				postpred = 1;
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

	AAMutSelSample sample(name,burnin,every,until);
	if (postpred)	{
		sample.PostPredCompositionalHeterogeneity();
	}
	else	{
		sample.Read();
	}

}



