
#include "Sample.h"
#include "CodonProfileMutSelMatInfMixModel.h"
//#include "NeffAAProfileInfMixMutSelModel.h"
//#include "MeanValTree.h"

class CodonProfileMutSelInfMixSample : public Sample	{

	private:
	string modeltype;
	string datafile;
	string treefile;
	int P;
	GeneticCodeType type;

	public:

	string GetModelType() {return modeltype;}

	CodonProfileMutSelMatInfMixModel* GetModel() {return (CodonProfileMutSelMatInfMixModel*) model;}
	//AAProfileInfMixMutSelModel* GetModel() {return (NeffAAProfileInfMixMutSelModel*) model;}

	CodonProfileMutSelInfMixSample(string filename, int inburnin, int inevery, int inuntil) : Sample(filename,inburnin,inevery,inuntil)	{
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
		if (modeltype == "CODONPROFILEINFMIXMUTSEL")	{
			model = new CodonProfileMutSelMatInfMixModel(datafile,treefile,P,true,type);
			//model = new AAProfileMutSelMatInfMixModel(datafile,treefile,P,true,type); replace true with false to read chain without sampling mapping, and without update
		}
		else	{
			cerr << "error when opening file "  << name << " : does not recognise model type : " << modeltype << '\n';
			exit(1);
		}

		// read model (i.e. chain's last point) from <name>.param
		model->FromStream(is);
		cerr << "UPDATE\n";
		//model->Update();

		// open <name>.chain, and prepare stream and stream iterator
		OpenChainFile();
		// now, size is defined (it is the total number of points with which this Sample object will make all its various posterior averages)
		// all these points can be accessed to (only once) by repeated calls to GetNextPoint()
	}

	// a very simple (and quite uninteresting) method for obtaining
	// the posterior mean and variance of the total length of the tree

	void SiteSpecificSelectionCoeffs()	{
		cout << "Not yet implemented\n";
		cout.flush();

		/*

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

		*/

	}


	//-------------------
	// ReadClusterModes
	//-------------------

	void ReadClusterModes(int SizeThreshold = 1, double ClusteringRange = 0.1)	{

		cout << "not yet implemented\n";
		cout.flush();
		/*

		int NClusterMax = 1000;

		//string BaseName = SampleName + ".modeclusterlogo";
		string BaseName = name + ".modeclusterlogo";
		int pageNumber = 1;
		int lineNumber = 0;
		ostringstream appel;
		//appel <<  "cp " << auxpath + modeheader << " " << BaseName;
		appel <<  "cp modeheader " << BaseName;
		system(appel.str().c_str());

		//ofstream logo_os(BaseName.c_str(), IOS_APPEND);
		ofstream logo_os(BaseName.c_str(), std::ios_base::app);
		logo_os << "%%Page: " << 1 << ' ' << 1 << '\n';
		logo_os << "startpage\n";
		logo_os << "startline\n";

		int Nstate = Naa;
		//int Nstate = mParam->Nstate;
		//int size = GetSize();
		int Nsite = GetModel()->aaprofilemutselmix->GetSize();

		// ---------------------------------------
		// creating clustering related arrays :
		// ---------------------------------------

		// NCluster : number of identified clusters
		// clusterTable(i,j) : frequency of amino acid j in cluster i
		// ClusterWeight[i] : weight of cluster i (mean number of sites affiliated to modes aggregated to this cluster)
		// ClusterSize[i] : size of cluster i (number of modes aggregated to this cluster)

		double  clusterTable[NClusterMax][Nstate];
		for (int i=0; i<NClusterMax; i++)	{
			for (int j=0; j<Nstate; j++)	{
				clusterTable[i][j] = 0;
			}
		}
		int NCluster = 0;

		int* ClusterWeight = new int[NClusterMax];
		for (int i=0; i<NClusterMax; i++)	{
			ClusterWeight[i] = 0;
		}

		int* ClusterSize = new int[NClusterMax];
		for (int i=0; i<NClusterMax; i++)	{
			ClusterSize[i] = 0;
		}

		int** IsAffiliated = new int*[NClusterMax];
		for (int i=0; i<NClusterMax; i++)	{
			IsAffiliated[i] = new int[size];
			for (int j=0; j<size; j++)	{
				IsAffiliated[i][j] = 0;
			}
		}
		int ClusterStability[NClusterMax];

		double** SiteAff = new double*[Nsite];
		for (int i=0; i<Nsite; i++)	{
			SiteAff[i] = new double[NClusterMax];
			for (int j=0; j<NClusterMax; j++)	{
				SiteAff[i][j] = 0;
			}
		}


		double* stat = new double[Naa];

		// -------------------------
		// reading sample file
		// -------------------------

		for (int i=0; i<size; i++)	{
			cerr << '.';
			cerr.flush();
			//PhyloBayes* pb = GetNextPB();
			//PhyloBayes& PB = *pb;
			//PB.UpdateSiteNumber();

			GetNextPoint();
			int Ncomponent = GetModel()->aaprofilemutselmix->GetComponentNumber();
			//for (int j=0; j<PB.GetModeNumber(); j++)	{
			for (int j=0; j<Ncomponent; j++)	{

				//double* stat = PB.Stationary[j];
				//int sitenumber = PB.SiteNumber[j];
				for (int aa=0; aa<Naa; aa++)	{
					stat[aa] = (*GetModel()->aaprofilemutselmix->GetComponent(j))[aa];
					//cout << stat[aa] << "\n";
					//cout.flush();
				}
				int sitenumber = GetModel()->aaprofilemutselmix->GetAllocationStatistic(j);

				int affiliated = -1;

				for (int k=0; k<NCluster; k++)	{
					// compare current mode to kth cluster
					// if close enough :
					// 	if not yet affiliated, then affiliated to cluster k
					// 	else affiliation cluster and cluster k are merged
					// 		and cluster k is flagged for elimination

					double kl = 0;
					double quad = 0;
					if (ClusterSize[k] == 0)	{
						cerr << "error while aggregating clusters : found an empty cluster\n";
						exit(1);
					}

					for (int l=0; l<Nstate; l++)	{
						kl += clusterTable[k][l]/ClusterWeight[k] * (log(clusterTable[k][l])/ClusterWeight[k] - log(stat[l]));
						double temp = clusterTable[k][l]/ClusterWeight[k] - stat[l];
						quad += temp * temp;
					}

					if (kl < 0)	{
					//	cerr << "error : negative KL divergence\n";
					//	exit(1);
					}

					if (quad < 0)	{
						cerr << "error : negative quadratic distance\n";
						exit(1);
					}

					if (quad < ClusteringRange)	{

						if (affiliated == -1)	{
							for (int l=0; l<Nstate; l++)	{
								clusterTable[k][l] += sitenumber * stat[l];
							}
							ClusterSize[k] ++;
							ClusterWeight[k] += sitenumber;
							affiliated = k;
							IsAffiliated[k][i] = 1;
							for (int site=0; site<Nsite; site++)	{
								//if (PB.Mode[site] == j)	{
								if (GetModel()->aaprofilemutselmix->GetAllocation(site) == j)	{
									SiteAff[site][k]++;
								}
							}
						}
						else	{
							// merging clusters

							for (int l=0; l<Nstate; l++)	{
								clusterTable[affiliated][l] += clusterTable[k][l];
								clusterTable[k][l] = 0;
							}
							ClusterSize[affiliated] += ClusterSize[k];
							for (int ii=0; ii<i; ii++)	{
								IsAffiliated[affiliated][ii] |= IsAffiliated[k][ii];
								IsAffiliated[k][ii] = 0;
							}
							IsAffiliated[affiliated][i] = 1;
							IsAffiliated[k][i] = 0;
							ClusterSize[k] = 0;
							ClusterWeight[affiliated] += ClusterWeight[k];
							ClusterWeight[k] = 0;
							for (int site=0; site<Nsite; site++)	{
								SiteAff[site][affiliated] += SiteAff[site][k];
								SiteAff[site][k] = 0;
							}
						}
					}

				}

				// look whether should stack down clusters

				// scan for first empty cluster : k;
				// scan for next non empty cluster : l;
				// while l < NCluster
				// else
				//	tranfer l to k
				// and increment l until a new non empty cluster is found or l == NCluster

				int k=0;
				while ( (k<NCluster) && (ClusterSize[k]))	{
					k++;
				}

				int l = k+1;
				while (l < NCluster)	{

					while ( (l<NCluster) && (! ClusterSize[l]))	{
						l++;
					}
					if (l<NCluster)	{
						// tranfer l to k;

						ClusterSize[k] = ClusterSize[l];
						for (int m=0; m<Nstate; m++)	{
							clusterTable[k][m] = clusterTable[l][m];
							clusterTable[l][m] = 0;
						}
						for (int ii=0; ii<=i; ii++)	{
							IsAffiliated[k][ii] = IsAffiliated[l][ii];
							IsAffiliated[l][ii] = 0;
						}
						for (int site=0; site<Nsite; site++)	{
							SiteAff[site][k] = SiteAff[site][l];
							SiteAff[site][l] = 0;
						}
						ClusterSize[l] = 0;
						ClusterWeight[k] = ClusterWeight[l];
						ClusterWeight[l] = 0;
						k++;
						l++;
					}
				}

				NCluster = k;
				if (NCluster > NClusterMax)	{
					cerr << "error : NCluster overflow\n";
					exit(1);
				}

				// if not affiliated, create a new cluster

				if (affiliated == -1)	{
					ClusterSize[NCluster] = 1;
					for (int l=0; l<Nstate; l++)	{
						clusterTable[NCluster][l] = sitenumber * stat[l];
					}
					ClusterWeight[NCluster] = sitenumber;
					for (int site=0; site<Nsite; site++)	{
						//if (PB.Mode[site] == j)	{
						if (GetModel()->aaprofilemutselmix->GetAllocation(site) == j)	{
							SiteAff[site][NCluster]=1;
						}
					}
					NCluster++;
					if (NCluster >= NClusterMax)	{
						cerr << "error : NCluster overflow\n";
						exit(1);
					}

				}

			}
		}
		cerr << '\n';

		// ---------------------------------------
		// computing cluster related statistics
		// ---------------------------------------

		// compute cluster profiles

		for (int k=0; k<NCluster; k++)	{

			if (! ClusterSize[k])	{
				cerr << "error : void cluster\n";
				exit(1);
			}

			for (int m=0; m<Nstate; m++)	{
				clusterTable[k][m] /= ClusterWeight[k];
			}
		}

		int totalweight = 0;
		for (int l=0; l<NCluster; l++)	{
			totalweight += ClusterWeight[l];
		}

		// cluster stability
		for (int l=0; l<NCluster; l++)	{
			ClusterStability[l] = 0;
		}
		for (int k=0; k<NCluster; k++)	{
			for (int l=0; l<size; l++)	{
				ClusterStability[k] += IsAffiliated[k][l];
			}
		}
		// ---------
		// output
		// ---------

		// prepare file for output

		ofstream Cluster_os( (name + ".clustermodes").c_str() );
		ofstream ClusterStability_os( (name + ".clusterstability").c_str() );

		int MaxHeight = 3;
		double logoThreshold = 0.001;
		int LinePerPage = 8;
		//int SiteLinePerPage = 2;
		double SpaceBetweenLines = -50;
		//double SiteSpaceBetweenLines = -130;

		// write down cluster
		// and make logos

		double totalWeight = 0;
		double* array = new double[Nstate];

		Cluster_os << "profile\tweight";
		for (int k=0; k<Nstate; k++)	{
			//Cluster_os << '\t' << mParam->Alphabet[k];
			Cluster_os << '\t' << AminoAcids[k];
		}
		Cluster_os << '\n' << '\n';

		int* clusterlist = new int[NCluster];
		int FinalNCluster = 0;

		double garbageclusterstat[Nstate];
		double garbageclusterweight = 0;
		for (int k=0; k<Nstate; k++)	{
			garbageclusterstat[k] = 0;
		}

		for (int i=0; i<NCluster; i++)	{

			int max = 0;
			int cluster = 0;
			for (int j=0; j<NCluster; j++)	{
				if (max < ClusterWeight[j])	{
					max = ClusterWeight[j];
					cluster = j;
				}
			}


			clusterlist[i] = cluster;

			double weight = ((double) ClusterWeight[cluster]) / size;

			totalWeight += weight;
			ClusterWeight[cluster] =0;

			if (weight >= SizeThreshold)	{

				FinalNCluster++;

				double max = 0;
				double min = 1;

				for (int k=0; k<Nstate; k++)	{
					if (min > clusterTable[cluster][k])	{
						min = clusterTable[cluster][k];
					}
					if (max < clusterTable[cluster][k])	{
						max = clusterTable[cluster][k];
					}
				}

				for (int l=0; l<Nstate; l++)	{
					double importance = clusterTable[cluster][l]/max;
					if (importance > 0.6)	{
						Cluster_os << AminoAcids[l];
					}
					else if (importance > 0.30)	{
						Cluster_os << aminoacids[l];
					}
				}

				for (int l=0; l<Nstate; l++)	{
					double importance = clusterTable[cluster][l]/max;
					if (importance > 0.6)	{
						ClusterStability_os << AminoAcids[l];
					}
					else if (importance > 0.30)	{
						ClusterStability_os << aminoacids[l];
					}
				}

				ClusterStability_os << '\t' << ((int) ( ((double) (ClusterStability[cluster]) / size) * 100)) << '\n';

				Cluster_os << '\t' << weight ;
				for (int l=0; l<Nstate; l++)	{
					Cluster_os << '\t' << clusterTable[cluster][l];
				}
				Cluster_os << '\n';

				// logo
				for (int k=0; k<Nstate; k++)	{

					logo_os << "gsave\n";

					array[k] = clusterTable[cluster][k] * MaxHeight ;
					if (array[k] > logoThreshold)	{
						logo_os << array[k] << " (" << AminoAcids[k] << ") numchar\n";
					}

					logo_os << "grestore\n";
					logo_os << "shift\n";
				}

				logo_os << "endline\n";

				logo_os << "0 0 0 setrgbcolor\n";
				logo_os << "240 30 moveto\n";

				ostringstream s;
				s << '(' << ((int) weight) << ") show\n" ;
				logo_os << s.str();


				lineNumber ++;

				if (lineNumber == LinePerPage)	{
					lineNumber = 0;
					pageNumber++;
					logo_os << "endpage\n";

					logo_os << "%%Page: " << pageNumber << ' ' << pageNumber << '\n';
					logo_os << "startpage\n";
					logo_os << "startline\n";
				}
				else		{
					logo_os << " 0 " << SpaceBetweenLines << " translate\n";
					logo_os << "startline\n";
				}
			}

			else 	{

				for (int k=0; k<Nstate; k++)	{
					garbageclusterstat[k] += weight * clusterTable[cluster][k];
				}
				garbageclusterweight += weight;

			}
		}

		FinalNCluster++;
		for (int k=0; k<Nstate; k++)	{
			garbageclusterstat[k] /= garbageclusterweight;
		}

		Cluster_os << '\n';

		double max = 0;
		double min = 1;

		for (int k=0; k<Nstate; k++)	{
			if (min > garbageclusterstat[k])	{
				min = garbageclusterstat[k];
			}
			if (max < garbageclusterstat[k])	{
				max = garbageclusterstat[k];
			}
		}

		for (int l=0; l<Nstate; l++)	{
			double importance = garbageclusterstat[l]/max;
			if (importance > 0.6)	{
				Cluster_os << AminoAcids[l];
			}
			else if (importance > 0.30)	{
				Cluster_os << aminoacids[l];
			}
		}
		Cluster_os << '\t' << garbageclusterweight ;
		for (int l=0; l<Nstate; l++)	{
			Cluster_os << '\t' << garbageclusterstat[l];
		}
		Cluster_os << '\n';

			cout << "number of clusters : " << FinalNCluster << '\n';

		Cluster_os << '\n';
		Cluster_os << "total number of clusters (not including the garbage-cluster) : " << FinalNCluster - 1 << '\n';
		Cluster_os << "total weight of garbage cluster : " << ((int) (garbageclusterweight / totalWeight * 100)) << " \%" << '\n';
		// site affiliations
		ofstream SiteAff_os( (name + ".siteaff").c_str() );
		SiteAff_os << "site\taffiliation frequencies\n";

		double check[FinalNCluster];
		for (int i=0; i<FinalNCluster; i++)	{
			check[i] = 0;
		}
		for (int i=0; i<Nsite; i++)	{
			SiteAff_os << i;
			double total = 0;
			for (int j=0; j<NCluster; j++)	{
				total += SiteAff[i][clusterlist[j]];
			}
			double wgarb = 0;
			for (int j=0; j<FinalNCluster-1; j++)	{
				double tmp = SiteAff[i][clusterlist[j]] / total;
				SiteAff_os << '\t' << tmp;
				check[j] += tmp;
				wgarb += tmp;
			}
			check[FinalNCluster-1] += 1 - wgarb;
			SiteAff_os << '\t' << 1 - wgarb << '\n';
		}
		//
		//for (int i=0; i<FinalNCluster; i++)	{
		//	cerr << i << '\t' << check[i] << '\n';
		//}
		//

		// logo
		for (int k=0; k<Nstate; k++)	{

			logo_os << "gsave\n";

			array[k] = garbageclusterstat[k] * MaxHeight ;
			if (array[k] > logoThreshold)	{
				logo_os << array[k] << " (" << AminoAcids[k] << ") numchar\n";
			}

			logo_os << "grestore\n";
			logo_os << "shift\n";
		}

		logo_os << "endline\n";

		logo_os << "0 0 0 setrgbcolor\n";
		logo_os << "240 30 moveto\n";

		ostringstream s;
		s << '(' << ((int) garbageclusterweight) << ") show\n" ;
		logo_os << s.str();

		lineNumber ++;

		// finishing up logo work

		logo_os << "endline\n";
		logo_os << "endpage\n";

		logo_os << "%%Trailer\n";
		logo_os << "%%Pages: " << pageNumber << "\n";

		logo_os.close();

		cout << '\n';
		cout << "cluster logos in " << name << ".modeclusterlogo" << '\n';
		cout << '\n';

		*/
	}

	//-----------------
	// Read()
	//-----------------

	void Read()	{


		cout << "IN READ\n";
		cout.flush();
		// create array in which to load profiles for every site
		int Nsite = GetModel()->codonprofilemutselmix->GetSize();
		int Nstate = GetModel()->GetNstate();
		double** meanCodonProfile = new double*[Nsite];
		double* tempCodonProfile = new double[Nstate];
		for (int i=0; i<Nsite; i++)	{
			meanCodonProfile[i] = new double[Nstate];
			for (int j=0; j<Naa; j++)	{
				meanCodonProfile[i][j] = 0;
			}
		}

		cout << "size is: " << size << "\n";


		// cycle over the sample
		for (int i=0; i<size; i++)	{

			GetNextPoint();

			//cout << GetModel()->aaprofilemutselmix->GetComponentNumber() << "\n";

			for (int site=0; site<Nsite; site++)	{
				for (int codon=0; codon<Nstate; codon++)	{
					tempCodonProfile[codon] = (*GetModel()->codonprofilemutselmix->GetRandomVariable(site))[codon];
					meanCodonProfile[site][codon] += tempCodonProfile[codon];
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
			for (int codon=0; codon<Nstate; codon++)	{
				meanCodonProfile[site][codon] /= size;
				os << meanCodonProfile[site][codon] << '\t';
			}
			os << '\n';
		}
		cout << '\n';

		ofstream osfreq((name + ".ssempfreqs").c_str());
		GetModel()->EmpiricalSiteSpecificFrequencies(osfreq);
	}

	void PostPredNonsynSubs()	{

		double obsnsmean;
		double prednsmean;
		double obsnsvar;
		double prednsvar;
		double obsaadiv;
		double predaadiv;
		int countPredNSMeanGreaterThanObs = 0;
		int countPredNSVarGreaterThanObs = 0;
		int countPredAADivGreaterThanObs = 0;
		ofstream osnsmean((name + ".nsmean").c_str());
		ofstream osnsvar((name + ".nsvar").c_str());
		ofstream osaadiv((name + ".aadiversity").c_str());
		osnsmean << "#obs\tpred\n";
		osnsvar << "#obs\tpred\n";
		osaadiv << "#obs\tpred\n";

		// cycle over the sample
		for (int i=0; i<size; i++)	{
			GetNextPoint();
			GetModel()->Update();
			GetModel()->SamplePosteriorMapping();

			// posterior statistics

			obsnsmean = GetModel()->NonsynSubMean();
			osnsmean << obsnsmean << '\t';
			obsnsvar = GetModel()->NonsynSubVariance();
			osnsvar << obsnsvar << '\t';
			obsaadiv = GetModel()->MeanObservedAADiversity();
			osaadiv << obsaadiv << '\t';

			GetModel()->SamplePosteriorPredictiveMapping();
			GetModel()->PostPredSample();  // should eventually be changed...

			// posterior predictive statistics

			prednsmean = GetModel()->NonsynSubMean();
			osnsmean << prednsmean << '\n';
			prednsvar = GetModel()->NonsynSubVariance();
			osnsvar << prednsvar << '\n';
			predaadiv = GetModel()->MeanPredictiveAADiversity();
			osaadiv << predaadiv << '\n';

			// keep counts for p-values

			if (prednsmean >= obsnsmean)	{
				countPredNSMeanGreaterThanObs++;
			}
			if (prednsvar >= obsnsvar)	{
				countPredNSVarGreaterThanObs++;
			}
			if (predaadiv >= obsaadiv)	{
				countPredAADivGreaterThanObs++;
			}

		}

		// Output p-values at the tails of files.

		double pvalue;
		pvalue = (double)(countPredNSMeanGreaterThanObs)/size;
		osnsmean << "#pvalue: " << pvalue << "\n";
		cout << "nsmean p-value: " << pvalue << "\n";

		pvalue = (double)(countPredNSVarGreaterThanObs)/size;
		osnsvar << "#pvalue: " << pvalue << "\n";
		cout << "nsvar p-value: " << pvalue << "\n";

		pvalue = (double)(countPredAADivGreaterThanObs)/size;
		osaadiv << "#pvalue: " << pvalue << "\n";
		cout << "aadiv p-value: " << pvalue << "\n";

		cout << '\n';
	}

};


int main(int argc, char* argv[])	{

	int burnin = 0;
	int every = 1;
	int until = -1;
	string name;

	int postpred = 0;
	int sssc = 0; // site-specific selection coefficients of a test sequence found in file <name>.testseq
	int cm = 0;	// cluster modes
	int SizeThreshold = 1;
	double ClusteringRange = 0.05;

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
			else if (s == "-cm")	{
				cm = 1;
				i++;
				if (i == argc) throw(0);
				s = argv[i];
				SizeThreshold = atoi(argv[i]);
				i++;
				if (i == argc) throw(0);
				s = argv[i];
				ClusteringRange = atof(argv[i]);
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

	CodonProfileMutSelInfMixSample sample(name,burnin,every,until);
	if (postpred)	{
		sample.PostPredNonsynSubs();
	}
	else	{
		if (sssc)	{
			sample.SiteSpecificSelectionCoeffs();
		}
		else if (cm)	{
			sample.ReadClusterModes(SizeThreshold, ClusteringRange);
		}
		else {
			sample.Read();
		}
	}

}



