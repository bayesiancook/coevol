
#include "Sample.h"
#include "AAProfileMutSelMatMixModel.h"
//#include "NeffAAProfileInfMixMutSelModel.h"
//#include "MeanValTree.h"

class AAProfileMutSelMixSample : public Sample	{

	private:
	string modeltype;
	string datafile;
	string treefile;
	string profilefile;
	int P;
	GeneticCodeType type;

	public:

	string GetModelType() {return modeltype;}

	AAProfileMutSelModel* GetModel() {return (AAProfileMutSelModel*) model;}
	//AAProfileInfMixMutSelModel* GetModel() {return (NeffAAProfileInfMixMutSelModel*) model;}

	AAProfileMutSelMixSample(string filename, int inburnin, int inevery, int inuntil) : Sample(filename,inburnin,inevery,inuntil)	{
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
		is >> profilefile;
		is >> chainevery >> chainuntil >> chainsize;
		// the chain's saving frequency, upper limit and current size
		// not to be confused with the sample's subsampling frequency, upper limit and size

		// make a new model depending on the type obtained from the file
		if (modeltype == "AAPROFILEMUTSEL")	{
			model = new AAProfileMutSelModel(datafile,treefile,P,true,type,profilefile);
		}
		else	{
			cerr << "error when opening file "  << name << " : does not recognise model type : " << modeltype << '\n';
			exit(1);
		}

		// read model (i.e. chain's last point) from <name>.param
		model->FromStream(is);
		cerr << "UPDATE\n";
		model->Update();

		// open <name>.chain, and prepare stream and stream iterator
		OpenChainFile();
		// now, size is defined (it is the total number of points with which this Sample object will make all its various posterior averages)
		// all these points can be accessed to (only once) by repeated calls to GetNextPoint()
	}

	void SiteSpecificSelectionCoeffs()	{

		string seqfile = ((name + ".testseq").c_str());
		SequenceAlignment* nucseqdata = new FileSequenceAlignment(seqfile); 
		GeneticCodeType type = Universal;
		CodonSequenceAlignment* codonseqdata = new CodonSequenceAlignment(nucseqdata, true, type);
		int Nstate = codonseqdata->GetNstate();
		int Nsite = codonseqdata->Nsite;
		double* tempAAProfile = new double[Naa];
		double* tempSelCoeffs = new double[Naa];
		for (int aa=0; aa<Naa; aa++)	{
			tempSelCoeffs[aa] = 0;
		}
		double** meanSelCoeffs = new double*[Nsite];
		for (int site=0; site<Nsite; site++)	{
			meanSelCoeffs[site] = new double[Naa];
			for (int aa=0; aa<Naa; aa++)	{
				meanSelCoeffs[site][aa] = 0;
			}
		}


		for (int i=0; i<size; i++)	{
			GetNextPoint();
			for (int site=0; site<Nsite; site++)	{
				for (int aa=0; aa<Naa; aa++)	{
					tempAAProfile[aa] = (*GetModel()->aaprofilemutselmix->GetRandomVariable(site))[aa];
					tempSelCoeffs[aa] = 0;
				}
				int codonFrom = codonseqdata->GetState(0,site);
				for (int codonTo=0; codonTo<Nstate; codonTo++)	{
					int posDiff = codonseqdata->GetCodonStateSpace()->GetDifferingPosition(codonFrom,codonTo);
					if (posDiff == 0 || posDiff==1 || posDiff==2)     {
						if (! codonseqdata->GetCodonStateSpace()->Synonymous(codonFrom,codonTo))	{
							int aaFrom = codonseqdata->GetCodonStateSpace()->Translation(codonFrom);
							int aaTo = codonseqdata->GetCodonStateSpace()->Translation(codonTo);
							tempSelCoeffs[aaTo] = log(tempAAProfile[aaTo]) - log(tempAAProfile[aaFrom]);
						}
					}
				}
				for (int aa=0; aa<Naa; aa++)	{
					meanSelCoeffs[site][aa] += tempSelCoeffs[aa];
				}
			}
		}

		ofstream os((name + ".siteselcoeffs").c_str());
		ofstream ostestseqaa((name + ".testseqaaprofile").c_str());
		os << Nsite << "\n";
		ostestseqaa << Nsite << "\n";
		for (int site=0; site<Nsite; site++)	{
			int state = codonseqdata->GetState(0,site);
			for (int aa=0; aa<Naa; aa++)	{
				meanSelCoeffs[site][aa] /= size;
				os << meanSelCoeffs[site][aa] << '\t';
				if (codonseqdata->GetCodonStateSpace()->Translation(state) == aa)	{
					ostestseqaa << 1.0 << '\t';
				}
				else	{
					ostestseqaa << 0.0 << '\t';
				}
			}
			os << '\n';
			ostestseqaa << '\n';
		}



		for (int site=0; site<Nsite; site++)	{
			delete[] meanSelCoeffs[site];
		}
		delete[] tempAAProfile;
		delete[] tempSelCoeffs;
		delete[] meanSelCoeffs;

	}

	void Read()	{


		int Nsite = GetModel()->aaprofilemutselmix->GetSize();
		ofstream os((name + ".weightedLogL").c_str());

		// cycle over the sample
		cout << "reading chain... ";
		cout.flush();
		for (int i=0; i<size; i++)	{
			cout << i+1 << " ";
			cout.flush();
			GetNextPoint();
			GetModel()->Update();
			double logL = 0;
			for (int site=0; site<Nsite; site++)	{
				double weightedSiteLikelihood = 0;
				for (int k = 0; k < P; k++)	{
					GetModel()->aaprofilemutselmix->SetAllocation(site, k);
					double w = GetModel()->aaprofilemutselmix->GetWeight(k);
					double l = GetModel()->GetSiteLogLikelihood(site);
					weightedSiteLikelihood += w * exp(l);
				}
				logL += log(weightedSiteLikelihood);
			}
			os << logL << '\n';
		}
		cout << " processed " << size << " sample points.\n";

		//ofstream osfreq((name + ".ssempfreqs").c_str());
		//GetModel()->EmpiricalSiteSpecificFrequencies(osfreq);

		/*

		cout << "IN READ\n";
		cout.flush();
		// create array in which to load profiles for every site
		int Nsite = GetModel()->aaprofilemutselmix->GetSize();
		double** meanAAProfile = new double*[Nsite];
		double* tempAAProfile = new double[Naa];
		for (int i=0; i<Nsite; i++)	{
			meanAAProfile[i] = new double[Naa];
			for (int j=0; j<Naa; j++)	{
				meanAAProfile[i][j] = 0;
			}
		}

		cout << "size is: " << size << "\n";


		cout << "reading chain...";
		// cycle over the sample
		for (int i=0; i<size; i++)	{
			//cout << "i: " << i << "\n";
			//cout.flush();
			GetNextPoint();

			//cout << GetModel()->aaprofilemutselmix->GetComponentNumber() << "\n";

			for (int site=0; site<Nsite; site++)	{
				for (int aa=0; aa<Naa; aa++)	{
					tempAAProfile[aa] = (*GetModel()->aaprofilemutselmix->GetRandomVariable(site))[aa];
					meanAAProfile[site][aa] += tempAAProfile[aa];
				}
			}
			//cout << (*GetModel()->aaprofilemutselmix->GetRandomVariable(0))[0] << "\n"; // site 0, just to check
			//tempAAProfile = GetModel()->aaprofilemutselmix->GetRandomVariable(0) << "\n"; // site 0, just to check

			//cout << *GetModel()->AAProfileConcentration << '\n';
			//cout << *GetModel()->
		}

		ofstream os((name + ".ssprofiles").c_str());
		os << Nsite << '\n';
		for (int site=0; site<Nsite; site++)	{
			for (int aa=0; aa<Naa; aa++)	{
				meanAAProfile[site][aa] /= size;
				os << meanAAProfile[site][aa] << '\t';
			}
			os << '\n';
		}
		cout << '\n';

		//ofstream osfreq((name + ".ssempfreqs").c_str());
		//GetModel()->EmpiricalSiteSpecificFrequencies(osfreq);
		*/
	}

	void PostPredNonsynSubs()	{
		cout << "Not yet availible for this model.\n";
		exit(1);
	}

};


int main(int argc, char* argv[])	{

	int burnin = 0;
	int every = 1;
	int until = -1;
	string name;

	int postpred = 0;
	int sssc = 0; // site-specific selection coefficients of a test sequence found in file <name>.testseq

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
			else if (s == "-sssc")	{
				sssc = 1;
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
		cerr << "readaaprofilemutsel [-x <burnin> <every> <until>] <chainname> \n";
		cerr << '\n';
		exit(1);
	}

	AAProfileMutSelMixSample sample(name,burnin,every,until);
	if (postpred)	{
		sample.PostPredNonsynSubs();
	}
	else	{
		if (sssc)	{
			sample.SiteSpecificSelectionCoeffs();
		}
		else {
			sample.Read();
		}
	}

}



