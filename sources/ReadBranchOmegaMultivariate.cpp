
#include "Sample.h"
#include "ConjugateBranchOmegaMultivariateModel.h"
#include "MeanValTree.h"
#include "MeanICTree.h"
#include "MeanCovMatrix.h"
#include "StringStreamUtils.h"
#include "MeanChronogram.h"
#include "MeanChronoBubbleTree.h"

#include "TexTab.h"
#include "MeanTimeLine.h"


class BranchOmegaMultivariateSample : public Sample	{

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

	int mutmodel;
	int gc;
	bool clampdiag;
	bool autoregressive;
	bool clamproot;
	bool clamptree;
	bool meanexp;
	GeneticCodeType type;

	int conjpath;
	double mappingfreq;
	int contdatatype;

	int omegaratiotree;

	bool normalise;
	int nrep;
	int ncycle;

	int df;

	string bounds;
	string mix;

	int nsplit;
	int withdrift;
	bool withtimeline;
	bool separatesyn;
	bool separateomega;

	int krkctype;
	int jitter;

	int uniformprior;
	string rootfile;

	int sample;

	public:

	string GetModelType() {return modeltype;}

	BranchOmegaMultivariateModel* GetModel() {return (BranchOmegaMultivariateModel*) model;}

	ConjugateBranchOmegaMultivariateModel* GetConjugateModel() {
		ConjugateBranchOmegaMultivariateModel* tmp = dynamic_cast<ConjugateBranchOmegaMultivariateModel*>(model);
		if (! tmp)	{
			cerr << "error in sample: dynamic cast to conjugate\n";
			exit(1);
		}
		return tmp;
	}

	BranchOmegaMultivariateSample(string filename, int inburnin, int inevery, int inuntil, int insample) : Sample(filename,inburnin,inevery,inuntil)	{
		sample = insample;

		bounds = "None";
		mix = "None";
		suffstatfile = "None";
		withtimeline = false;
		separatesyn = false;
		separateomega = false;
		krkctype = 0;
		jitter = 0;
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
		string priorsigmafile = "None";

		nsplit = 1;
		withdrift = 0;
		uniformprior = 0;
		rootfile = "None";
        softa = 0;

		// read model type, and other standard fields
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
		is >> chainevery >> chainuntil >> chainsize;
		// the chain's saving frequency, upper limit and current size
		// not to be confused with the sample's subsampling frequency, upper limit and size

		// make a new model depending on the type obtained from the file
		if (modeltype == "BRANCHOMEGAMULTIVARIATE")	{
			model = new BranchOmegaMultivariateModel(datafile,treefile,contdatafile,calibfile,rootage,rootstdev,chronoprior,softa,meanchi,meanchi2,priorsigma,priorsigmafile,df,mutmodel,gc,clampdiag,autoregressive,conjpath,mappingfreq,contdatatype,omegaratiotree,clamproot,clamptree,meanexp,normalise,nrep,ncycle,bounds,mix,nsplit,withdrift,uniformprior,rootfile,suffstatfile,withtimeline,separatesyn,separateomega,krkctype,jitter,0,1,sample,type);
		}
		else if (modeltype == "CONJUGATEBRANCHOMEGAMULTIVARIATE")	{
			model = new ConjugateBranchOmegaMultivariateModel(datafile,treefile,contdatafile,calibfile,rootage,rootstdev,chronoprior,softa,meanchi,meanchi2,priorsigma,priorsigmafile,df,mutmodel,gc,autoregressive,conjpath,mappingfreq,contdatatype,omegaratiotree,clamproot,clamptree,meanexp,normalise,nrep,ncycle,bounds,mix,nsplit,withdrift,uniformprior,rootfile,suffstatfile,withtimeline,separatesyn,separateomega,krkctype,jitter,0,1,sample,type);
		}
		else	{
			cerr << "error when opening file "  << name << " : does not recognise model type : " << modeltype << '\n';
			exit(1);
		}

		// read model (i.e. chain's last point) from <name>.param
		model->FromStream(is);
		// model->Update();
		// model->GetLogProb();

		// open <name>.chain, and prepare stream and stream iterator
		OpenChainFile();
		// now, size is defined (it is the total number of points with which this Sample object will make all its various posterior averages)
		// all these points can be accessed to (only once) by repeated calls to GetNextPoint()

		cerr << "number of points to read : " << size << '\n';
		cerr << '\n';
	}

	void PostPredNormalityTest()	{

		ConjugateMultiVariateTreeProcess* process = GetConjugateModel()->GetConjugateMultiVariateTreeProcess();
		int dim = process->GetDim();
		double* ppk = new double[dim];
		double* ppz1 = new double[dim];
		double* ppz2 = new double[dim];
		for (int k=0; k<dim; k++)	{
			ppk[k] = 0;
			ppz1[k] = 0;
			ppz2[k] = 0;
		}
		double* obsk = new double[dim];
		double* obsz1 = new double[dim];
		double* obsz2 = new double[dim];
		double* predk = new double[dim];
		double* predz1 = new double[dim];
		double* predz2 = new double[dim];

		for (int i=0; i<size; i++)	{
			cerr << '.';
			GetNextPoint();
			GetModel()->Update();

			GetConjugateModel()->GetConjugateInverseWishart(0)->ActivateSufficientStatistic();
			process->DAgostinosKandZ(obsk,obsz1,obsz2);
			GetConjugateModel()->GetConjugateInverseWishart(0)->InactivateSufficientStatistic();

			GetModel()->SampleProcess();

			GetConjugateModel()->GetConjugateInverseWishart(0)->ActivateSufficientStatistic();
			process->DAgostinosKandZ(predk,predz1,predz2);
			GetConjugateModel()->GetConjugateInverseWishart(0)->InactivateSufficientStatistic();

			for (int k=0; k<dim; k++)	{
				if (obsk[k] > predk[k])	{
					ppk[k] ++;
				}
				if (obsz1[k] > predz1[k])	{
					ppz1[k] ++;
				}
				if (obsz2[k] > predz2[k])	{
					ppz2[k] ++;
				}
			}
		}
		cerr << '\n';

		cout << "#\tK\tskewness\tkurtosis\n";
		cout << '\n';
		for (int k=0; k<dim; k++)	{
			ppk[k] /= size;
			ppz1[k] /= size;
			ppz2[k] /= size;
			cout << k << '\t' << ppk[k] << '\t' << ppz1[k] << '\t' << ppz2[k] << '\n';
		}


	}

	/*
	void PostStrandSymmetry()	{

		double ct2ga = 0;
		double cg2ta = 0;
		double gt2ca = 0;
		for (int i=0; i<size; i++)	{
			cerr << '.';
			GetNextPoint();
			const double* rr = GetModel()->GetRelativeRates();
			if (rr[0] == rr[5])	{
				gt2ca++;
			}
			if (rr[1] == rr[4])	{
				ct2ga++;
			}
			if (rr[2] == rr[3])	{
				cg2ta++;
			}
		}
		cerr << "gt2ca : " << gt2ca / size << '\n';
		cerr << "ct2ga : " << ct2ga / size << '\n';
		cerr << "cg2ta : " << cg2ta / size << '\n';

	}

	void ReadLogKappa()	{

		MultiVariateTreeProcess* process = GetModel()->GetMultiVariateProcess();
		int dim = process->GetDim();
		double* meanlogkappa = new double[dim];
		double* varlogkappa = new double[dim];
		for (int k=0; k<dim; k++)	{
			meanlogkappa[k] = 0;
			varlogkappa[k] = 0;
		}
		for (int i=0; i<size; i++)	{
			cerr << '.';
			GetNextPoint();
			for (int k=0; k<dim; k++)	{
				double tmp = log(GetModel()->DiagArray->GetVal(k)->val());
				meanlogkappa[k] += tmp;
				varlogkappa[k] += tmp * tmp;
			}
		}
		cerr << '\n';
		for (int k=0; k<dim; k++)	{
			meanlogkappa[k] /= size;
			varlogkappa[k] /= size;
			varlogkappa[k] -= meanlogkappa[k] * meanlogkappa[k];
			cout << meanlogkappa[k] << '\t' << sqrt(varlogkappa[k]) << '\t' << exp(meanlogkappa[k]) << '\t' << exp(sqrt(varlogkappa[k])) << '\n';
		}

	}


	void ReadBF()	{

		double* logpostarray = new double[size];
		double* logpriorarray = new double[size];
		double maxpost = 0;
		double maxprior = 0;
		for (int i=0; i<size; i++)	{
			cerr << '.';
			GetNextPoint();
			GetModel()->Update();
			GetConjugateModel()->GetConjugateInverseWishart(0)->ActivateSufficientStatistic();
			logpostarray[i] = GetConjugateModel()->GetConjugateInverseWishart(0)->DiagonalLogPosterior();
			logpriorarray[i] = GetConjugateModel()->GetConjugateInverseWishart(0)->DiagonalLogPrior();
			GetConjugateModel()->GetConjugateInverseWishart(0)->InactivateSufficientStatistic();
			if ((!i) || (maxpost < logpostarray[i]))	{
				maxpost = logpostarray[i];
			}
			if ((!i) || (maxprior < logpriorarray[i]))	{
				maxprior = logpriorarray[i];
			}
		}
		double totpost = 0;
		for (int i=0; i<size; i++)	{
			totpost += exp(logpostarray[i] - maxpost);
		}
		totpost /= size;
		double logpost = log(totpost) + maxpost;
		double totprior = 0;
		for (int i=0; i<size; i++)	{
			totprior += exp(logpriorarray[i] - maxprior);
		}
		totprior /= size;
		double logprior = log(totprior) + maxprior;
		// double logprior = GetConjugateModel()->GetConjugateInverseWishart(0)->DiagonalLogPrior(0.001,1000.0);

		cout << '\n';
		cout << "log prior : " << logprior << '\n';
		cout << "log post  : " << logpost << '\n';
		cout << "log BF    : " << logpost - logprior << '\n';
		cout << '\n';
	}
	*/

	void ReadIC()	{

		MeanICTree* ictree = new MeanICTree(GetModel()->GetMultiVariateProcess());
		MeanICTree* ppictree = new MeanICTree(GetModel()->GetMultiVariateProcess());
		ictree->Reset();
		ppictree->Reset();
		int dim = GetModel()->GetMultiVariateProcess()->GetDim();
		double* pp = new double[dim];
		for (int k=0; k<dim; k++)	{
			pp[k] = 0;
		}
		for (int i=0; i<size; i++)	{
			cerr << '.';
			GetNextPoint();
			ictree->Add();
			// GetModel()->Update();
			GetModel()->SampleProcess();
			ppictree->Add();
			for (int k=0; k<dim; k++)	{
				if (ppictree->GetGmean(k) > ictree->GetGmean(k))	{
					pp[k]++;
				}
			}
		}
		cerr << '\n';
		ictree->Normalise();
		ppictree->Normalise();
		ofstream os((name + ".postmeanic").c_str());
		ictree->Tabulate(os,false);
		ofstream los((name + ".postmeanleafic").c_str());
		ictree->Tabulate(los,true);
		cerr << "mean independent contrasts in " << name  << ".postmeanic\n";
		cerr << "for terminal branches only in " << name  << ".postmeanleafic\n";
		cerr << "pp positive trend : \n";
		for (int k=0; k<dim; k++)	{
			cerr << k << '\t' << pp[k] / size << '\n';
		}
		cerr << '\n';
		ictree->DAgostinosKandZ();
	}

	void PostPredTrend()	{

		MultiVariateTreeProcess* process = GetModel()->GetMultiVariateProcess();
		int dim = process->GetDim();
		double* obs = new double[dim];
		double* pred = new double[dim];
		double* mean = new double[dim];
		double* var = new double[dim];
		double* ppmean = new double[dim];
		double* ppvar = new double[dim];
		double* pp = new double[dim];
		for (int k=0; k<dim; k++)	{
			mean[k] = 0;
			var[k] = 0;
			ppmean[k] = 0;
			ppvar[k] = 0;
			pp[k] = 0;
		}

		for (int j=0; j<size; j++)	{
			cerr << '.';
			GetNextPoint();
			for (int k=0; k<dim; k++)	{
				obs[k] = process->GetLeafMean(k) - (*process->GetMultiNormal(GetModel()->GetTree()->GetRoot()))[k];
				mean[k] += obs[k];
				var[k] += obs[k] * obs[k];
			}
			GetModel()->SampleProcess();
			for (int k=0; k<dim; k++)	{
				pred[k] = process->GetLeafMean(k) - (*process->GetMultiNormal(GetModel()->GetTree()->GetRoot()))[k];
				ppmean[k] += pred[k];
				ppvar[k] += pred[k] * pred[k];
				if (pred[k] > obs[k])	{
					pp[k] ++;
				}
			}

		}
		cerr << '\n';
		ofstream os((GetName() + ".trends").c_str());
		for (int k=0; k<dim; k++)	{
			mean[k] /= size;
			var[k] /= size;
			var[k] -= mean[k] * mean[k];
			ppmean[k] /= size;
			ppvar[k] /= size;
			ppvar[k] -= ppmean[k] * ppmean[k];
			pp[k] /= size;
			os << pp[k] << '\t' << mean[k] << '\t' << sqrt(var[k]) << '\t' << ppmean[k] << '\t' << sqrt(ppvar[k]) << '\n';
		}
	}

	void MakeTimeLine()	{
		int N = 100;
		MeanTimeLine* synrate = new MeanTimeLine(N);
		MeanTimeLine* omega = 0;
		MeanTimeLine* omegats = 0;
		MeanTimeLine* omegatv0 = 0;
		MeanTimeLine* omegatvgc = 0;
		MeanTimeLine* gc = 0;
		if (GetModel()->Has3Omega())	{
			omegats = new MeanTimeLine(N);
			omegatv0 = new MeanTimeLine(N);
			omegatvgc = new MeanTimeLine(N);
		}
		else if (GetModel()->Has2Omega())	{
			omegats = new MeanTimeLine(N);
			omegatv0 = new MeanTimeLine(N);
		}
		else if (GetModel()->HasOmega())	{
			omega = new MeanTimeLine(N);
		}
		if (GetModel()->isGCActivated())	{
			gc = new MeanTimeLine(N);
		}

		int Ncont = GetModel()->Ncont;
		MeanTimeLine** tree = new MeanTimeLine*[Ncont];
		for (int k=0; k<Ncont; k++)	{
			tree[k] = new MeanTimeLine(N);
		}

		for (int j=0; j<size; j++)	{
			cerr << '.';
			GetNextPoint();
			if (!GetModel()->Unconstrained() && ! GetModel()->SeparateSyn())	{
				GetModel()->GetSynRateTree()->specialUpdate();
			}
			GetModel()->GetChronogram()->specialUpdate();
			if (!GetModel()->Unconstrained() && ! GetModel()->SeparateSyn())	{
				synrate->Add(GetModel()->GetSynRateTree(),GetModel()->GetChronogram(),true);
			}
			if (GetModel()->Has3Omega())	{
				GetModel()->GetOmegaTSTree()->specialUpdate();
				omegats->Add(GetModel()->GetOmegaTSTree(),GetModel()->GetChronogram(),false);
				GetModel()->GetOmegaTV0Tree()->specialUpdate();
				omegatv0->Add(GetModel()->GetOmegaTV0Tree(),GetModel()->GetChronogram(),false);
				GetModel()->GetOmegaTVGCTree()->specialUpdate();
				omegatvgc->Add(GetModel()->GetOmegaTVGCTree(),GetModel()->GetChronogram(),false);
			}
			else if (GetModel()->Has2Omega())	{
				GetModel()->GetOmegaTSTree()->specialUpdate();
				omegats->Add(GetModel()->GetOmegaTSTree(),GetModel()->GetChronogram(),false);
				GetModel()->GetOmegaTV0Tree()->specialUpdate();
				omegatv0->Add(GetModel()->GetOmegaTV0Tree(),GetModel()->GetChronogram(),false);
			}
			else if (GetModel()->HasOmega())	{
				GetModel()->UpdateOmegaTree();
				omega->Add(GetModel()->GetOmegaTree(),GetModel()->GetChronogram(),false);
			}
			if (GetModel()->isGCActivated())	{
				GetModel()->GetMeanLogitGCTree()->specialUpdate();
				gc->Add(GetModel()->GetGCTree(),GetModel()->GetChronogram(),false);
			}
			for (int k=0; k<Ncont; k++)	{
				tree[k]->Add(GetModel()->GetMultiVariateProcess(),GetModel()->GetChronogram(),GetModel()->GetL() + k);
			}
		}
		cerr << '\n';
		synrate->Normalize();
		ofstream sos((GetName() + ".synratetimeline").c_str());
		synrate->ToStream(sos);

		if (GetModel()->Has3Omega())	{
			omegats->Normalize();
			ofstream otsos((GetName() + ".omegatstimeline").c_str());
			omegats->ToStream(otsos);
			omegatv0->Normalize();
			ofstream otv0os((GetName() + ".omegatv0timeline").c_str());
			omegatv0->ToStream(otv0os);
			omegatvgc->Normalize();
			ofstream otvgcos((GetName() + ".omegatvgctimeline").c_str());
			omegatvgc->ToStream(otvgcos);
		}
		else if (GetModel()->Has2Omega())	{
			omegats->Normalize();
			ofstream otsos((GetName() + ".omegatstimeline").c_str());
			omegats->ToStream(otsos);
			omegatv0->Normalize();
			ofstream otv0os((GetName() + ".omegatv0timeline").c_str());
			omegatv0->ToStream(otv0os);
		}
		else if (GetModel()->HasOmega())	{
			omega->Normalize();
			ofstream oos((GetName() + ".omegatimeline").c_str());
			omega->ToStream(oos);
		}
		if (GetModel()->isGCActivated())	{
			gc->Normalize();
			ofstream gcos((GetName() + ".gctimeline").c_str());
			gc->ToStream(gcos);
		}

		for (int k=0; k<Ncont; k++)	{
			tree[k]->Normalize();
			ostringstream s;
			s << GetName() << "." << k << "timeline";
			ofstream os(s.str().c_str());
			tree[k]->ToStream(os);
		}
	}

	void ReadTimeLine()	{

		int timesize = GetModel()->GetTimeLine()->GetSize();
		int K = GetModel()->GetCovMatrix()->GetDim();
		double* date = new double[timesize];
		for (int i=0; i<timesize; i++)	{
			date[i] = 1 - GetModel()->GetTimeIntervals()->GetDate(i);
		}
		double** current = new double*[timesize];
		double** meandiff = new double*[timesize];
		double** vardiff = new double*[timesize];
		double** ppdiff = new double*[timesize];
		double** mean = new double*[timesize];
		double** var = new double*[timesize];
		for (int i=0; i<timesize; i++)	{
			mean[i] = new double[K];
			var[i] = new double[K];
			current[i] = new double[K];
			meandiff[i] = new double[K];
			vardiff[i] = new double[K];
			ppdiff[i] = new double[K];
			for (int k=0; k<K; k++)	{
				mean[i][k] = 0;
				var[i][k] = 0;
				current[i][k] = 0;
				meandiff[i][k] = 0;
				vardiff[i][k] = 0;
				ppdiff[i][k] = 0;
			}
		}

		for (int j=0; j<size; j++)	{
			cerr << '.';
			GetNextPoint();
			for (int i=0; i<timesize; i++)	{
				for (int k=0; k<K; k++)	{
					double tmp = (*GetModel()->GetTimeLine()->GetVal(i))[k];
					current[i][k] = tmp;
					double tmp2 = (*GetModel()->GetTimeLine()->GetVal(2))[k];
					tmp-=tmp2;
					mean[i][k] += tmp;
					var[i][k] += tmp * tmp;
				}
			}
			for (int i=1; i<timesize; i++)	{
				for (int k=0; k<K; k++)	{
					double tmp = (current[i][k] - current[i-1][k]) / sqrt(date[i] - date[i-1]);
					meandiff[i][k] += tmp;
					vardiff[i][k] += tmp * tmp;
					if (tmp > 0)	{
						ppdiff[i][k] ++;
					}
				}
			}
		}
		cerr << '\n';

		ofstream os((name + ".timeline").c_str());
		for (int i=0; i<timesize; i++)	{
			os << date[i];
			for (int k=0; k<K; k++)	{
				mean[i][k] /= size;
				var[i][k] /= size;
				var[i][k] -= mean[i][k] * mean[i][k];

				meandiff[i][k] /= size;
				vardiff[i][k] /= size;
				vardiff[i][k] -= meandiff[i][k] * meandiff[i][k];
				if (!i)	{
					vardiff[i][k] = 0;
				}
				ppdiff[i][k] /= size;
				os << '\t' << mean[i][k] << '\t' << sqrt(var[i][k]) << '\t';
				// os << '\t' << mean[i][k] << '\t' << sqrt(var[i][k]) << '\t' << meandiff[i][k] << '\t' << sqrt(vardiff[i][k]) << '\t' << ppdiff[i][k] << '\t';
			}
			os << '\n';
		}
	}

	void DrawDensities(string taxpairfile, int cont)	{

		MeanExpNormTree* tree = new MeanExpNormTree(GetModel()->GetFineGrainedTree(),false,true,true,false,false,true,true,true);

		// cycle over the sample
		for (int i=0; i<size; i++)	{
			cerr << '.';

			// get next point -> will be stored into "model", and thus, will be accessible through GetModel()

			GetNextPoint();
			tree->Add(GetModel()->GetMultiVariateProcess(), GetModel()->GetLengthTree(), GetModel()->GetL()+cont);
		}
		cerr << '\n';
		tree->Normalise();

		ifstream is(taxpairfile.c_str());
		int ntax;
		is >> ntax;
		map<pair<string,string>,int> m;
		for (int i=0; i<ntax; i++)	{
			string name1, name2, name3;
			is >> name1 >> name2 >> name3;
			m[pair<string,string>(name1,name2)] = 1;
		}
		ostringstream s;
		s << GetName() << ".postmean" << cont+1 << ".tab";
		ofstream os(s.str().c_str());
		// tree->TabulateDistribution(os,&m);

	}

	void Read(bool printlog, bool printmean, bool printmed, bool printci, bool printstdev, bool withleaf, bool withinternal, string mulreg, bool tex, double xscale, double yscale, double nodescale, double nodepower, double barwidth, int fontsize, bool bubbletext, bool withheader, double leafnameshift, double meanreg, double stdevreg, bool postdist)	{

		double meansyn = 0;
		double varsyn = 0;
		double meantime = 0;
		double vartime = 0;

		int Ncont = GetModel()->Ncont;
		if (GetModel()->Split())	{
			NewickTree::Simplify();
		}

		MeanBranchTree* meantree = new MeanBranchTree(GetModel()->GetFineGrainedTree(),false);

		MeanChronogram* meanchrono = 0;
		if (! GetModel()->Unconstrained())	{
			meanchrono = new MeanChronogram(GetModel()->GetTree());
		}

		MeanExpNormTree* meansynrate = new MeanExpNormTree(GetModel()->GetFineGrainedTree(),false,printlog,printmean,printmed,printci,printstdev,withleaf,withinternal,meanreg,stdevreg);
		MeanExpNormTree* meanbranchsynrate = new MeanExpNormTree(GetModel()->GetFineGrainedTree(),false,printlog,printmean,printmed,printci,printstdev,withleaf,withinternal);

		MeanExpNormTree* meanomega = 0;
		MeanExpNormTree* meanomegats = 0;
		MeanExpNormTree* meanomegatv0 = 0;
		MeanExpNormTree* meanomegatvgc = 0;
		MeanExpNormTree* meanbranchomega = 0;

		if (GetModel()->Has3Omega())	{
			// cerr << "3 omega : " << GetModel()->GetOmegaTsIndex() << '\t' << GetModel()->GetOmegaTv0Index() << '\t' << GetModel()->GetOmegaTvGCIndex() << '\n';
			meanomegats = new MeanExpNormTree(GetModel()->GetFineGrainedTree(),false,printlog,printmean,printmed,printci,printstdev,withleaf,withinternal);
			meanomegatv0 = new MeanExpNormTree(GetModel()->GetFineGrainedTree(),false,printlog,printmean,printmed,printci,printstdev,withleaf,withinternal);
			meanomegatvgc = new MeanExpNormTree(GetModel()->GetFineGrainedTree(),false,printlog,printmean,printmed,printci,printstdev,withleaf,withinternal);
		}
		else if (GetModel()->Has2Omega())	{
			// cerr << "2 omega : " << GetModel()->GetOmegaTsIndex() << '\t' << GetModel()->GetOmegaTv0Index() << '\n';
			meanomegats = new MeanExpNormTree(GetModel()->GetFineGrainedTree(),false,printlog,printmean,printmed,printci,printstdev,withleaf,withinternal);
			meanomegatv0 = new MeanExpNormTree(GetModel()->GetFineGrainedTree(),false,printlog,printmean,printmed,printci,printstdev,withleaf,withinternal);
		}
		else if (GetModel()->HasOmega())	{
			// cerr << "1 omega : " << GetModel()->GetOmegaIndex() << '\n';
			meanomega = new MeanExpNormTree(GetModel()->GetFineGrainedTree(),false,printlog,printmean,printmed,printci,printstdev,withleaf,withinternal);
			meanbranchomega = new MeanExpNormTree(GetModel()->GetFineGrainedTree(),false,printlog,printmean,printmed,printci,printstdev,withleaf,withinternal);
		}

		MeanExpNormTree* meantstv = 0;
		MeanExpNormTree* meantvgc = 0;

		if (GetModel()->HasTsTv())	{
			// cerr << "TsTv : " << GetModel()->GetTsTvIndex() << '\n';
			meantstv = new MeanExpNormTree(GetModel()->GetFineGrainedTree(),false,printlog,printmean,printmed,printci,printstdev,withleaf,withinternal);
		}
		if (GetModel()->HasTvGC())	{
			// cerr << "TvGC: " << GetModel()->GetTvGCIndex() << '\n';
			meantvgc = new MeanExpNormTree(GetModel()->GetFineGrainedTree(),false,printlog,printmean,printmed,printci,printstdev,withleaf,withinternal);
		}

		MeanExpNormTree* meangc = 0;
		MeanExpNormTree* meangc1 = 0;
		MeanExpNormTree* meangc2 = 0;
		MeanExpNormTree* meangc3 = 0;
		if (GetModel()->isGCActivated())	{
			// cerr << "GC : " << GetModel()->GetGCIndex() << '\n';
			meangc = new MeanExpNormTree(GetModel()->GetFineGrainedTree(),true,printlog,printmean,printmed,printci,printstdev,withleaf,withinternal);
		}
		if (GetModel()->isGC3Activated())	{
			meangc1 = new MeanExpNormTree(GetModel()->GetFineGrainedTree(),true,printlog,printmean,printmed,printci,printstdev,withleaf,withinternal);
			meangc2 = new MeanExpNormTree(GetModel()->GetFineGrainedTree(),true,printlog,printmean,printmed,printci,printstdev,withleaf,withinternal);
			meangc3 = new MeanExpNormTree(GetModel()->GetFineGrainedTree(),true,printlog,printmean,printmed,printci,printstdev,withleaf,withinternal);
		}

		MeanExpNormTree** tree = new MeanExpNormTree*[Ncont];
		for (int k=0; k<Ncont; k++)	{
			tree[k] = new MeanExpNormTree(GetModel()->GetFineGrainedTree(),false,printlog,printmean,printmed,printci,printstdev,withleaf,withinternal);
		}

		int dim = GetModel()->GetCovMatrix()->GetDim();
		MeanCovMatrix*  mat = new MeanCovMatrix(dim);

		double logmixalpha = 0;
		ofstream los((GetName() + ".logmixalpha").c_str());

		for (int i=0; i<size; i++)	{
			cerr << '.';

			GetNextPoint();

			if (!GetModel()->Unconstrained() && ! GetModel()->SeparateSyn())	{
				GetModel()->GetSynRateTree()->specialUpdate();
			}
			if (! GetModel()->Unconstrained())	{
				GetModel()->GetChronogram()->specialUpdate();
			}
			if (!GetModel()->Unconstrained() && ! GetModel()->SeparateSyn())	{
				double r = GetModel()->GetMeanSynRate();
				meansyn += r;
				varsyn += r*r;
			}

			double t = 0;
			if (! GetModel()->Unconstrained())	{
				t = GetModel()->GetTotalTime();
			}
			meantime += t;
			vartime += t*t;

			double tmp = log(GetModel()->GetMixAlpha());
			logmixalpha += tmp;
			los << tmp << '\n';


			if (! GetModel()->Unconstrained())	{
				meanchrono->Add(GetModel()->GetChronogram());
			}

			if (GetModel()->Split())	{
				GetModel()->GetSplitLengthTree()->specialUpdate();
			}

			if (!GetModel()->Unconstrained() && ! GetModel()->SeparateSyn())	{
				meantree->Add(GetModel()->GetSynRateTree());
			}

			if (!GetModel()->Unconstrained() && ! GetModel()->SeparateSyn())	{
				meanbranchsynrate->Add(GetModel()->GetSynRateTree(), GetModel()->GetLengthTree(),true);
			}
			meansynrate->Add(GetModel()->GetMultiVariateProcess(), GetModel()->GetLengthTree(), 0);
			if (GetModel()->isGCActivated())	{
				meangc->Add(GetModel()->GetMultiVariateProcess(), GetModel()->GetLengthTree(), GetModel()->GetGCIndex());
			}
			if (GetModel()->isGC3Activated())	{
				meangc1->Add(GetModel()->GetMultiVariateProcess(), GetModel()->GetLengthTree(), GetModel()->GetGCIndex());
				meangc2->Add(GetModel()->GetMultiVariateProcess(), GetModel()->GetLengthTree(), GetModel()->GetGCIndex()+1);
				meangc3->Add(GetModel()->GetMultiVariateProcess(), GetModel()->GetLengthTree(), GetModel()->GetGCIndex()+2);
			}

			if (GetModel()->Has3Omega())	{
				meanomegats->Add(GetModel()->GetMultiVariateProcess(), GetModel()->GetLengthTree(), GetModel()->GetOmegaTsIndex());
				meanomegatv0->Add(GetModel()->GetMultiVariateProcess(), GetModel()->GetLengthTree(), GetModel()->GetOmegaTv0Index());
				meanomegatvgc->Add(GetModel()->GetMultiVariateProcess(), GetModel()->GetLengthTree(), GetModel()->GetOmegaTvGCIndex());
			}
			else if (GetModel()->Has2Omega())	{
				meanomegats->Add(GetModel()->GetMultiVariateProcess(), GetModel()->GetLengthTree(), GetModel()->GetOmegaTsIndex());
				meanomegatv0->Add(GetModel()->GetMultiVariateProcess(), GetModel()->GetLengthTree(), GetModel()->GetOmegaTv0Index());
			}
			else if (GetModel()->HasOmega())	{
				if (GetModel()->rawNstate == Naa)	{
					double radcons = GetModel()->GetRadConsNormFactor();
					meanomega->Add(GetModel()->GetMultiVariateProcess(), GetModel()->GetLengthTree(), GetModel()->GetOmegaIndex(),radcons);
				}
				else	{
					meanomega->Add(GetModel()->GetMultiVariateProcess(), GetModel()->GetLengthTree(), GetModel()->GetOmegaIndex());
				}
				if (! GetModel()->SeparateOmega())	{
					meanbranchomega->Add(GetModel()->GetNonSynRateOrOmegaTree(), GetModel()->GetLengthTree(),false);
				}
			}

			if (GetModel()->HasTsTv())	{
				meantstv->Add(GetModel()->GetMultiVariateProcess(), GetModel()->GetLengthTree(), GetModel()->GetTsTvIndex());
			}
			if (GetModel()->HasTvGC())	{
				meantvgc->Add(GetModel()->GetMultiVariateProcess(), GetModel()->GetLengthTree(), GetModel()->GetTvGCIndex());
			}

			for (int k=0; k<Ncont; k++)	{
				tree[k]->Add(GetModel()->GetMultiVariateProcess(), GetModel()->GetLengthTree(), GetModel()->GetL()+k);
			}

			CovMatrix& m = *(GetModel()->GetCovMatrix());
			mat->Add(&m);
		}
		cerr << '\n';
		cerr << '\n';

		if (!GetModel()->Unconstrained() && ! GetModel()->SeparateSyn())	{
			meansyn /= size;
			varsyn /= size;
			varsyn -= meansyn * meansyn;
			// cerr << '\n';
			// cerr << "mean total syn : " << meansyn << " +/- " << sqrt(varsyn) << '\n';
		}

		meantime /= size;
		vartime /= size;
		vartime -= meantime * meantime;
		// cerr << "mean total time: " << meantime << " +/- " << sqrt(vartime) << '\n';
		// cerr << "average syn rate per tree depth : " << meansyn / meantime << '\n';
		// cerr << '\n';

		logmixalpha /= size;
		if (logmixalpha)	{
			cerr << "mean log mix alpha: " << logmixalpha << '\n';
		}

		mat->Normalize();

		if (mulreg != "")	{
			ofstream cout((GetName() + ".marginalcov").c_str());
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
			PartialMeanCovMatrix* partmat = new PartialMeanCovMatrix(mat,indexarray,dim);

			ReducedMeanCovMatrix* redmat = new ReducedMeanCovMatrix(mat,indexarray,dim);
			cout << '\n';
			cout << "partial correl\n";
			cout << '\n';
			cout << *partmat;

			ofstream cout2((GetName() + ".controlcov").c_str());
			cout2 << "entries are in the following order:\n";
			GetModel()->PrintEntries(cout2, indexarray);
			cout2 << '\n';
			cout2 << "reduced correl\n";
			cout2 << '\n';
			redmat->PrintCovariances(cout2);
			redmat->PrintCorrel(cout2);
			redmat->PrintPosteriorProbs(cout2);
			// cout2 << *redmat;
			cout2 << '\n';
			delete partmat;
			delete redmat;
		}
		ofstream cout((GetName() + ".cov").c_str());
		cout << "entries are in the following order:\n";
		GetModel()->PrintEntries(cout);

		cout << '\n';
		mat->SetLatex(tex);
		cout << *mat;

		cerr << "covariance matrix in " << name << ".cov\n";
		cerr << '\n';

		if (postdist)	{

			if (!GetModel()->Unconstrained() && ! GetModel()->SeparateSyn())	{
				ofstream os((GetName() + ".postdistsynrate.tab").c_str());
				meansynrate->TabulateDistribution(os);
				os.close();
			}

			ofstream os((GetName() + ".postdistbranchsynrate.tab").c_str());
			meanbranchsynrate->TabulateDistribution(os);
			os.close();

			if (GetModel()->Has3Omega())	{
				ofstream otsoos((GetName() + ".postdistomegats.tab").c_str());
				meanomegats->TabulateDistribution(otsoos);
				otsoos.close();

				ofstream otv0oos((GetName() + ".postdistomegatv0.tab").c_str());
				meanomegatv0->TabulateDistribution(otv0oos);
				otv0oos.close();

				ofstream otvgcoos((GetName() + ".postdistomegatvgc.tab").c_str());
				meanomegatvgc->TabulateDistribution(otvgcoos);
				otvgcoos.close();
			}
			else if (GetModel()->Has2Omega())	{
				ofstream otsoos((GetName() + ".postdistomegats.tab").c_str());
				meanomegats->TabulateDistribution(otsoos);
				otsoos.close();

				ofstream otv0oos((GetName() + ".postdistomegatv0.tab").c_str());
				meanomegatv0->TabulateDistribution(otv0oos);
				otv0oos.close();
			}
			else if (GetModel()->Has1Omega())	{
				ofstream os((GetName() + ".postdistomega.tab").c_str());
				meanomega->TabulateDistribution(os);
				os.close();

				ofstream bos((GetName() + ".postdistbranchomega.tab").c_str());
				meanbranchomega->TabulateDistribution(bos);
				bos.close();
			}

			if (GetModel()->isGCActivated())	{
				ofstream gcos((GetName() + ".postdistgc.tab").c_str());
				meangc->TabulateDistribution(gcos);
				gcos.close();
			}

			for (int k=0; k<Ncont; k++)	{
				ostringstream s;
				s << GetName() << ".postdist" << k+1 << ".tab";
				ofstream os(s.str().c_str());
				tree[k]->TabulateDistribution(os);
			}
		}

		meantree->Normalise();
		ofstream tos((GetName() + ".postmean.tre").c_str());
		meantree->ToStream(tos);

		if (! GetModel()->Unconstrained())	{
			meanchrono->Normalise();
			ofstream chos((GetName() + ".postmeandates.tre").c_str());
			meanchrono->ToStream(chos);

			ofstream cchos((GetName() + ".postmeandates.tab").c_str());
			meanchrono->Tabulate(cchos);
		}

		if (GetModel()->Has3Omega())	{
			meanomegats->Normalise();
			ofstream otsos((GetName() + ".postmeanomegats.tre").c_str());
			meanomegats->ToStream(otsos);
			meanomegatv0->Normalise();
			ofstream otv0os((GetName() + ".postmeanomegatv0.tre").c_str());
			meanomegatv0->ToStream(otv0os);
			meanomegatvgc->Normalise();
			ofstream otvgcos((GetName() + ".postmeanomegatvgc.tre").c_str());
			meanomegatvgc->ToStream(otvgcos);
			cerr << "reconstructed variation in omega in " << name << ".postmeanomega.tre\n";
		}
		else if (GetModel()->Has2Omega())	{
			meanomegats->Normalise();
			ofstream otsos((GetName() + ".postmeanomegats.tre").c_str());
			meanomegats->ToStream(otsos);
			meanomegatv0->Normalise();
			ofstream otv0os((GetName() + ".postmeanomegatv0.tre").c_str());
			meanomegatv0->ToStream(otv0os);
			cerr << "reconstructed variation in omega in " << name << ".postmeanomega.tre\n";
		}
		else if (GetModel()->HasOmega())	{
			meanomega->Normalise();
			ofstream oos((GetName() + ".postmeanomega.tre").c_str());
			meanomega->ToStream(oos);
			cerr << "reconstructed variation in omega in " << name << ".postmeanomega.tre\n";
			// cerr << "pp of mean leaf values > root value : " << meanomega->GetPPLeafRoot() << '\n';
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
		else if (GetModel()->isGCActivated())	{
			meangc->Normalise();
			ofstream gcos((GetName() + ".postmeangc.tre").c_str());
			meangc->ToStream(gcos);
			cerr << "reconstructed variation in gc in " << name << ".postmeangc.tre\n";
			// cerr << "pp of mean leaf values > root value : " << meangc->GetPPLeafRoot() << '\n';
		}

		if (GetModel()->HasTvGC())	{
			meantstv->Normalise();
			meantvgc->Normalise();
			ofstream tstvos((GetName() + ".postmeants_tv0.tre").c_str());
			ofstream tvgcos((GetName() + ".postmeantvgc_tv0.tre").c_str());
			meantstv->ToStream(tstvos);
			meantvgc->ToStream(tvgcos);
			cerr << "reconstructed variation in Ts/Tv0 in " << name << ".postmeants_tv0.tre\n";
			cerr << "reconstructed variation in TvGC/Tv0 in " << name << ".postmeantvgc_tv0.tre\n";
		}
		else if (GetModel()->HasTsTv())	{
			meantstv->Normalise();
			ofstream tstvos((GetName() + ".postmeants_tv.tre").c_str());
			meantstv->ToStream(tstvos);
			cerr << "reconstructed variation in Ts/Tv in " << name << ".postmeants_tv.tre\n";
		}

		if (!GetModel()->Unconstrained() && ! GetModel()->SeparateSyn())	{
			meansynrate->Normalise();
			ofstream sos((GetName() + ".postmeansynrate.tre").c_str());
			meansynrate->ToStream(sos);
			cerr << "reconstructed variation in dS in " << name << ".postmeansynrate.tre\n";
			// cerr << "pp of mean leaf values > root value : " << meansynrate->GetPPLeafRoot() << '\n';
			if (tex)	{
				ostringstream s;
				s << GetName() << ".postmeansynrate.tre";
				MeanChronoBubbleTree* textree = new MeanChronoBubbleTree(meanchrono,meansynrate,xscale,yscale,nodescale,nodepower,barwidth,fontsize,bubbletext,withheader,leafnameshift);
				textree->Draw((s.str() + ".tex").c_str());
			}
		}

		meanbranchsynrate->Normalise();
		ofstream bsos((GetName() + ".postmeanbranchsynrate.tre").c_str());
		meanbranchsynrate->ToStream(bsos);
		cerr << "reconstructed mean dS per branch in " << name << ".postmeanbranchsynrate.tre\n";

		for (int k=0; k<Ncont; k++)	{
			tree[k]->Normalise();
			ostringstream s;
			s << GetName() << ".postmean" << k+1 << ".tre";
			ofstream os(s.str().c_str());
			tree[k]->ToStream(os);
			cerr << "reconstructed variation in " << GetModel()->GetContinuousData()->GetCharacterName(k) << " in "  << name << ".postmean" << k+1 << ".tre\n";
			// cerr << "pp of mean leaf values > root value : " << tree[k]->GetPPLeafRoot() << '\n';
			if (tex && (! GetModel()->Unconstrained()))	{
				MeanChronoBubbleTree* textree = new MeanChronoBubbleTree(meanchrono,tree[k],xscale,yscale,nodescale,nodepower,barwidth,fontsize,bubbletext,withheader,leafnameshift);
				textree->Draw((s.str() + ".tex").c_str());
			}
		}

		if (!GetModel()->Unconstrained() && ! GetModel()->SeparateSyn())	{
			ofstream ssos((GetName() + ".postmeansynrate.tab").c_str());
			meansynrate->Tabulate(ssos);
			ssos.close();
		}

		ofstream bssos((GetName() + ".postmeanbranchsynrate.tab").c_str());
		meanbranchsynrate->Tabulate(bssos);
		bssos.close();

		if (GetModel()->Has3Omega())	{
			ofstream otsoos((GetName() + ".postmeanomegats.tab").c_str());
			meanomegats->Tabulate(otsoos);
			otsoos.close();

			ofstream otv0oos((GetName() + ".postmeanomegatv0.tab").c_str());
			meanomegatv0->Tabulate(otv0oos);
			otv0oos.close();

			ofstream otvgcoos((GetName() + ".postmeanomegatvgc.tab").c_str());
			meanomegatvgc->Tabulate(otvgcoos);
			otvgcoos.close();

		}
		if (GetModel()->Has2Omega())	{
			ofstream otsoos((GetName() + ".postmeanomegats.tab").c_str());
			meanomegats->Tabulate(otsoos);
			otsoos.close();

			ofstream otv0oos((GetName() + ".postmeanomegatv0.tab").c_str());
			meanomegatv0->Tabulate(otv0oos);
			otv0oos.close();
		}
		else if (GetModel()->Has1Omega())	{
			ofstream ooos((GetName() + ".postmeanomega.tab").c_str());
			meanomega->Tabulate(ooos);
			ooos.close();

			meanbranchomega->Normalise();
			ofstream bsos((GetName() + ".postmeanbranchnonsynrate.tre").c_str());
			meanbranchomega->ToStream(bsos);
			ofstream bssos((GetName() + ".postmeanbranchnonsynrate.tab").c_str());
			meanbranchomega->Tabulate(bssos);
			bssos.close();

			
		}
		if (GetModel()->isGCActivated())	{
			ofstream gcos((GetName() + ".postmeangc.tab").c_str());
			meangc->Tabulate(gcos);
			gcos.close();
		}

		for (int k=0; k<Ncont; k++)	{
			ostringstream s;
			s << GetName() << ".postmean" << k+1 << ".tab";
			ofstream os(s.str().c_str());
			tree[k]->Tabulate(os);
		}
		cerr << '\n';
	}

	int rrindex(int i, int j, int nstate)	{
		return (i<j) ? (2 * nstate - i - 1) * i / 2 + j - i - 1 : (2 * nstate - j - 1) * j / 2 + i - j - 1 ;
	}

	void ReadRelRates()	{

		int nrr = Nnuc * (Nnuc-1) / 2;
		double* meanrr = new double[nrr];
		double* varrr = new double[nrr];
		for (int k=0; k<nrr; k++)	{
			meanrr[k] = 0;
			varrr[k] = 0;
		}

		for (int i=0; i<size; i++)	{
			GetNextPoint();
			const double* rr = GetModel()->GetRelativeRates();
			for (int k=0; k<nrr; k++)	{
				meanrr[k] += rr[k];
				varrr[k] += rr[k] * rr[k];
			}
		}

		for (int k=0; k<nrr; k++)	{
			meanrr[k] /= size;
			varrr[k] /= size;
			varrr[k] -= meanrr[k] * meanrr[k];
			cout << meanrr[k] << '\t' << sqrt(varrr[k]) << '\n';
		}
		cout << '\n';
		for (int k=0; k<Nnuc; k++)	{
			for (int l=k+1; l<Nnuc; l++)	{
				cout << k << '\t' << l << '\t' << meanrr[rrindex(k,l,Nnuc)] << '\t' << sqrt(varrr[rrindex(k,l,Nnuc)]) << '\n';
			}
		}

		delete[] meanrr;
		delete[] varrr;
	}

	void ReadNuc(string taxon)	{
		list<double> meanat2gc;
		list<double> meangc2at;
		list<double> meangc2ta;
		list<double> meanat2cg;
		list<double> meangc2cg;

		for (int i=0; i<size; i++)	{
			cerr << '.';

			GetNextPoint();

			double* rr = GetModel()->exprelrate->GetArray();
			double* stat = 0;
			if (GetModel()->isGCActivated())	{
				GetModel()->GetMeanLogitGCTree()->specialUpdate();
				GetModel()->GetGCStatTree()->specialUpdate();
				const Link* link = GetModel()->GetFineGrainedTree()->GetLCA(taxon,taxon);
				stat = GetModel()->GetGCStatTree()->GetBranchVal(link->GetBranch())->GetArray();
				// cerr << stat[0] << '\t' << stat[1] << '\t' << stat[2] << '\t' << stat[3] << '\t' << stat[0] + stat[1] + stat[2] + stat[3] << '\n';
			}
			else	{
				stat = GetModel()->stationary->GetArray();
			}

			double tmp = rr[2] * stat[3] + rr[2] * stat[0];
			meanat2gc.push_front((rr[1] * stat[2] + rr[4] * stat[1]) / tmp);
			meangc2at.push_front((rr[1] * stat[0] + rr[5] * stat[3]) / tmp);
			meangc2ta.push_front((rr[5] * stat[3] + rr[0] * stat[0]) / tmp);
			meanat2cg.push_front((rr[0] * stat[1] + rr[5] * stat[2]) / tmp);
			meangc2cg.push_front((rr[3] * stat[1] + rr[3] * stat[2]) / tmp);
		}

		cerr << '\n';

		ofstream cout((GetName() + "_table2.tex").c_str());
		cout << taxon << "&\t";
		cout << textabentry(meanat2gc,true,false,true,2);
		cout <<  "&\t";
		cout << textabentry(meangc2at,true,false,true,2);
		cout <<  "&\t";
		cout << textabentry(meangc2ta,true,false,true,2);
		cout <<  "&\t";
		cout << textabentry(meanat2cg,true,false,true,2);
		cout <<  "&\t";
		cout << textabentry(meangc2cg,true,false,true,2);
		cout << "\\\\\n";
		cout.close();

	}

	void PostPred()	{

		double obscorrel = GetModel()->GetObsCorrelStat();

		double mean = 0;
		double var = 0;
		double pp = 0;

		for (int i=0; i<size; i++)	{
			cerr << '.';
			// get next point -> will be stored into "model", and thus, will be accessible through GetModel()
			GetNextPoint();

			GetModel()->Update();

			double tmp = GetModel()->PostPredCompo();
			mean += tmp;
			var += tmp * tmp;
			if (tmp > obscorrel)	{
				pp ++;
			}

			// GetModel()->PostPredSimu(s.str());
		}

		mean /= size;
		var /= size;
		var -= mean * mean;
		pp /= size;

		cout << '\n';
		cout << "observed  : " << obscorrel << '\n';
		cout << "post pred : " << mean << '\n';
		cout << "z score   : " << (mean - obscorrel)/sqrt(var) << '\n';
		cout << "pp        : " << pp << '\n';
		cout << '\n';
	}

	void CheckCov(string truefile)	{

		int dim = GetModel()->GetCovMatrix()->GetDim();

		double trueval[dim][dim];
		double mean[dim][dim];
		double var[dim][dim];
		double totval[dim][dim][size];
		double meancorrel[dim][dim];
		double varcorrel[dim][dim];
		double truecorrel[dim][dim];
		double totcorrel[dim][dim][size];

		for (int j=0; j<dim; j++)	{
			for (int k=0; k<dim; k++)	{
				mean[j][k] = 0;
				var[j][k] = 0;
				meancorrel[j][k] = 0;
				varcorrel[j][k] = 0;
			}
		}

		cerr << "TRUE  " << truefile << '\n';
		ifstream is(truefile.c_str());
		// ifstream is((truefile + "cov").c_str());
		if (!is)	{
			cerr << "error: did not find file : " << truefile << '\n';
			exit(1);
		}

		/*
		int tmpdim;
		is >> tmpdim;
		if (tmpdim != dim)	{
			cerr << "error when reading true file\n";
			cerr << tmpdim << '\t' << dim << '\n';
			exit(1);
		}
		*/
		for (int j=0; j<dim; j++)	{
			for (int k=0; k<dim; k++)	{
				is >> trueval[j][k];
			}
		}
		for (int j=0; j<dim; j++)	{
			for (int k=j+1; k<dim; k++)	{
				truecorrel[j][k] = trueval[j][k] / sqrt(trueval[j][j] * trueval[k][k]);
			}
		}

		for (int i=0; i<size; i++)	{
			// get next point -> will be stored into "model", and thus, will be accessible through GetModel()
			GetNextPoint();
			for (int j=0; j<dim; j++)	{
				for (int k=j; k<dim; k++)	{
					double tmp = (*GetModel()->GetCovMatrix())[j][k];
					mean[j][k] += tmp;
					var[j][k] += tmp * tmp;
					totval[j][k][i] = tmp;
				}
			}
			for (int j=0; j<dim; j++)	{
				for (int k=j+1; k<dim; k++)	{
					double tmp = (*GetModel()->GetCovMatrix())[j][k] / sqrt( (*GetModel()->GetCovMatrix())[j][j] * (*GetModel()->GetCovMatrix())[k][k]);
					meancorrel[j][k] += tmp;
					totcorrel[j][k][i] = tmp;
				}
			}
		}

		ofstream vos((name + ".varcheck").c_str());
		for (int j=0; j<dim; j++)	{
			// for (int k=j; k<dim; k++)	{
				int k = j;
				mean[j][k] /= size;
				var[j][k] /= size;
				var[j][k] -= mean[j][k] * mean[j][k];
				int n = 0;
				for (int i=0; i<size; i++)	{
					if (trueval[j][k] > totval[j][k][i])	{
						n++;
					}
				}

				for (int i=0; i<size; i++)	{
					for (int l=size-1; l>i; l--)	{
						if (totval[j][k][i] > totval[j][k][l])	{
							double tmp = totval[j][k][i];
							totval[j][k][i] = totval[j][k][l];
							totval[j][k][l] = tmp;
						}
					}
				}
				int min = ((int) (((double) (size)) / 100 * 2.5));
				vos << trueval[j][k] << '\t' << mean[j][k] << '\t' << totval[j][k][min] << '\t' << totval[j][k][size-min-1] << '\t' << var[j][k] << '\t' << ((double) n)/size << '\n';
				// os << trueval[j][k] << '\t' << mean[j][k] << '\t' << totval[j][k][min] << '\t' << totval[j][k][size-min-1] << '\t';
			// }
		}
		// os << '\n';

		ofstream cos((name + ".covcheck").c_str());
		for (int j=0; j<dim; j++)	{
			for (int k=j+1; k<dim; k++)	{
				mean[j][k] /= size;
				var[j][k] /= size;
				var[j][k] -= mean[j][k] * mean[j][k];

				meancorrel[j][k] /= size;
				varcorrel[j][k] /= size;
				varcorrel[j][k] -= meancorrel[j][k] * meancorrel[j][k];

				int n = 0;
				for (int i=0; i<size; i++)	{
					if (trueval[j][k] > totval[j][k][i])	{
						n++;
					}
				}

				for (int i=0; i<size; i++)	{
					for (int l=size-1; l>i; l--)	{
						if (totval[j][k][i] > totval[j][k][l])	{
							double tmp = totval[j][k][i];
							totval[j][k][i] = totval[j][k][l];
							totval[j][k][l] = tmp;
						}
					}
				}

				for (int i=0; i<size; i++)	{
					for (int l=size-1; l>i; l--)	{
						if (totcorrel[j][k][i] > totcorrel[j][k][l])	{
							double tmp = totcorrel[j][k][i];
							totcorrel[j][k][i] = totcorrel[j][k][l];
							totcorrel[j][k][l] = tmp;
						}
					}
				}
				int min = ((int) (((double) (size)) / 100 * 2.5));
				cos << trueval[j][k] << '\t' << mean[j][k] << '\t' << totval[j][k][min] << '\t' << totval[j][k][size-min-1] << '\t' << sqrt(var[j][k]) << '\t' << ((double) n)/size << '\t' << truecorrel[j][k] << '\t' << meancorrel[j][k] << '\t' << totcorrel[j][k][min] << '\t' << totcorrel[j][k][size-min-1] << '\t' << sqrt(varcorrel[j][k]) << '\n';
				// os << trueval[j][k] << '\t' << mean[j][k] << '\t' << totval[j][k][min] << '\t' << totval[j][k][size-min-1] << '\t';
			}
		}
		// os << '\n';

	}


	void Check(string truefile)	{

		int Ncheck = GetModel()->GetNcheck();
		cerr << Ncheck << '\n';
		cerr << "size : " << size << '\n';
		double** tab = new double*[size];
		for (int i=0; i<size; i++)	{
			tab[i] = new double[Ncheck];
		}

		double** ttab = new double*[Ncheck];
		for (int i=0; i<size; i++)	{
			ttab[i] = new double[size];
		}

		double* trueval = new double[Ncheck];
		double* mean = new double[Ncheck];
		for (int k=0; k<Ncheck; k++)	{
			mean[k] = 0;
		}
		double* var = new double[Ncheck];
		for (int k=0; k<Ncheck; k++)	{
			var[k] = 0;
		}

		// read true value

		cerr << "TRUE  " << truefile << '\n';
		ifstream is((truefile + "param").c_str());
		if (!is)	{
			cerr << "error: did not find file : " << truefile << '\n';
			exit(1);
		}

		GetModel()->FromStream(is);
		ofstream cos((truefile  + "cov").c_str());
		cos << *(GetModel()->GetCovMatrix()) << '\n';
		cos.close();
		exit(1);
		GetModel()->GetCheck(trueval);

		for (int i=0; i<size; i++)	{
			// get next point -> will be stored into "model", and thus, will be accessible through GetModel()
			GetNextPoint();
			GetModel()->GetCheck(tab[i]);
			for (int k=0; k<Ncheck; k++)	{
				ttab[k][i] = tab[i][k];
			}
			for (int k=0; k<Ncheck; k++)	{
				mean[k] += tab[i][k];
				var[k] += tab[i][k] * tab[i][k];
			}
		}
		ofstream os((name + ".check").c_str());
		ofstream vos((name + ".vcheck").c_str());
		for (int k=0; k<Ncheck; k++)	{
			for (int i=0; i<size; i++)	{
				for (int j=size-1; j>i; j--)	{
					if (ttab[k][i] > ttab[k][j])	{
						double tmp = ttab[k][i];
						ttab[k][i] = ttab[k][j];
						ttab[k][j] = tmp;
					}
				}
			}
		}

		for (int k=0; k<Ncheck; k++)	{
			mean[k] /= size;
			var[k] /= size;
			var[k] -= mean[k] * mean[k];
			int n = 0;
			for (int i=0; i<size; i++)	{
				if (trueval[k] > tab[i][k])	{
					n++;
				}
			}
			os << ((double) n) / size << '\t';
			int min = ((int) (((double) (size)) / 100 * 2.5));
			vos << trueval[k] << '\t' << mean[k] << '\t' << ttab[k][min] << '\t' << ttab[k][size-min-1] << '\t';
			// vos << trueval[k] << '\t' << mean[k] << '\t' << sqrt(var[k]) << '\t' << ((double) n) / size << '\t';
		}
		os << '\n';
		vos << '\n';

		for (int k=0; k<Ncheck; k++)	{
			delete[] tab[k];
		}
		delete[] tab;
		delete[] trueval;
	}

	void ReadNe(bool printlog, bool printmean, bool printmed, bool printci, bool printstdev, bool withleaf, bool withinternal, double meanreg, double stdevreg)	{

		int Ncont = GetModel()->Ncont;
		int dim = GetModel()->GetCovMatrix()->GetDim();

		MeanChronogram* meanchrono = new MeanChronogram(GetModel()->GetTree());
		MeanExpNormTree* meansynrate = new MeanExpNormTree(GetModel()->GetTree(),false,printlog,printmean,printmed,printci,printstdev,withleaf,withinternal,meanreg,stdevreg);

        meansynrate->SetLogScale(10.0);

		MeanExpNormTree* meanomega = new MeanExpNormTree(GetModel()->GetTree(),false,printlog,printmean,printmed,printci,printstdev,withleaf,withinternal);

		MeanExpNormTree* meanNe = new MeanExpNormTree(GetModel()->GetTree(),false,printlog,printmean,printmed,printci,printstdev,withleaf,withinternal);
        meanNe->SetLogScale(10.0);

		MeanExpNormTree* meanu = new MeanExpNormTree(GetModel()->GetTree(),false,printlog,printmean,printmed,printci,printstdev,withleaf,withinternal);
        meanu->SetLogScale(10.0);

        // index of dS: 0
        int idxdS = 0;
        // index of dNdS: 1
        int idxdNdS = 1;
        // index of generation_time in continuous data
		int idxgentime(-1);
        // index of piS in continuous data
		int idxpiS(-1);
        // index of piNpiS in cont data
        int idxpiNpiS(-1);

        int idxu(dim);
        int idxNe(dim+1);
		
		for (int k=0; k<Ncont; k++)	{
			if (GetModel()->GetContinuousData()->GetCharacterName(k) == "generation_time") {
				idxgentime= k+dim-Ncont;
			}	
			else if (GetModel()->GetContinuousData()->GetCharacterName(k) == "piS") {
				idxpiS = k+dim-Ncont;
			}	
            else if (GetModel()->GetContinuousData()->GetCharacterName(k) == "piNpiS")  {
                idxpiNpiS = k+dim-Ncont;
            }
		}
		
		if (idxgentime == -1)  {
            cerr << "error: cannot find entry generation_time in continuous data matrix\n";
			exit(1);
		}
		
		if (idxpiS == -1)  {
            cerr << "error: cannot find entry piS in continuous data matrix\n";
			exit(1);
		}

        /*
		if (idxpiNpiS == -1)  {
            cerr << "error: cannot find entry piNpiS in continuous data matrix\n";
			exit(1);
		}
        */

        /*
        cerr << "dim       : " << dim << '\n';
        cerr << "idxdS     : " << idxdS << '\n';
        cerr << "idxpiS    : " << idxpiS << '\n';
        cerr << "idxpiNpiS : " << idxpiNpiS << '\n';
        cerr << "idxgentime: " << idxgentime << '\n';
        */
		
        // dS: mutation rate per tree depth
        // tau: generation time in days
        // rootage: tree depth in myr
        // rootage * 365.10^6: tree depth in days
        //
        // mutation rate per generation:
        // u = dS * tau / (rootage * 365.10^6)
        // log u = log dS + log tau - log(rootage * 365.10^6)
        //
        // effective population size:
        // Ne = pi_S / 4 / u = pi_s / 4 / dS / tau * (rootage * 365.10^6)
        // log Ne = log pi_S - log dS - log tau + log(rootage * 365.10^6 / 4)

        // alphau and alphaNe: slopes
        vector<double> alphau(dim, 0);
		alphau[idxdS] = 1;
		alphau[idxgentime] = 1;
		
        vector<double> alphaNe(dim, 0);
		alphaNe[idxdS] = -1;
		alphaNe[idxgentime] = -1;
		alphaNe[idxpiS] = 1;

        // betau and betaNe: offsets (but dependent on root age, and so, defined point by point below)
        
        // matrix giving all original variables + logu and logNe as a function of original variables
        vector<vector<double>> A(dim+2, vector<double>(dim+2, 0));
        // original variables
        for (int k=0; k<dim; k++)   {
            A[k][k] = 1.0;
        }
        // log u
        A[idxu][idxdS] = 1.0;
        A[idxu][idxgentime] = 1.0;
        // log Ne
        A[idxNe][idxdS] = -1.0;
        A[idxNe][idxgentime] = -1.0;
        A[idxNe][idxpiS] = 1.0;

		MeanExpNormTree** tree = new MeanExpNormTree*[Ncont];
		for (int k=0; k<Ncont; k++)	{
			tree[k] = new MeanExpNormTree(GetModel()->GetTree(),false,printlog,printmean,printmed,printci,printstdev,withleaf,withinternal);
		}

		MeanCovMatrix*  meancov = new MeanCovMatrix(dim);
        MeanCovMatrix* expandedmeancov = new MeanCovMatrix(dim+2, false);

		// cycle over the sample
		for (int i=0; i<size; i++)	{
			cerr << '.';

			GetNextPoint();
            // GetModel()->UpdateLengthTree();
			GetModel()->GetSynRateTree()->specialUpdate();

			double t0 = GetModel()->GetRootAge();
            /*
            double t0 = 0;
			if (!iscalspe) {
				t0 = GetModel()->GetRootAge();
			}
			else {
				t0 = rootage;
			}
            */

			meanchrono->Add(GetModel()->GetChronogram());

            // dS: mutation rate per tree depth
            // tree depth in years: rootage * 10^6
            // mutation rate per year: dS / rootage / 10^6
            double synrate_offset = -log(t0 * 1000000);
			meansynrate->Add(GetModel()->GetMultiVariateProcess(), GetModel()->GetLengthTree(), 0, synrate_offset);

			meanomega->Add(GetModel()->GetMultiVariateProcess(), GetModel()->GetLengthTree(), 1);

            // offsets for log-linear combinations for logNe and logu phylogenetic histories
            double betau = -log(365.0*t0*1000000.0);
            double betaNe = log(365.0*t0*1000000.0/4.0);
				
			meanNe->AddLogLinearCombination(GetModel()->GetMultiVariateProcess(), GetModel()->GetLengthTree(), alphaNe, betaNe);
			meanu->AddLogLinearCombination(GetModel()->GetMultiVariateProcess(), GetModel()->GetLengthTree(), alphau, betau);

			for (int k=0; k<Ncont; k++)	{
				tree[k]->Add(GetModel()->GetMultiVariateProcess(), GetModel()->GetLengthTree(), GetModel()->GetL()+k);
			}

            if (! clampdiag)    {
                // original covariance matrix (dS, dN/dS, traits, generation time, piS and piN/piS)
                CovMatrix& m = *(GetModel()->GetCovMatrix());
                meancov->Add(&m);
                
                CovMatrix expandedcov(dim+2);
                double midmat[dim+2][dim];
                for (int i=0; i<dim+2; i++) {
                    for (int j=0; j<dim; j++)    {
                        double tmp = 0;
                        for (int k=0; k<dim; k++)   {
                            tmp += A[i][k] * m[k][j];
                        }
                        midmat[i][j] = tmp;
                    }
                }
                for (int i=0; i<dim+2; i++) {
                    for (int j=0; j<dim+2; j++) {
                        double tmp = 0;
                        for (int k=0; k<dim; k++)   {
                            tmp += midmat[i][k] * A[j][k];
                        }
                        expandedcov[i][j] = tmp;
                    }
                }
                expandedmeancov->Add(&expandedcov);
            }
			
		}
		cerr << '\n';

        if (! clampdiag)    {

            meancov->Normalize();
            expandedmeancov->Normalize();

            ofstream cov_os((GetName() + ".cov").c_str());
            cov_os.precision(3);
            cov_os << "entries are in the following order:\n";
            cov_os << "dS\n";
            cov_os << "dN/dS\n";
            for (int k=0; k<GetModel()->Ncont; k++)  {
				cov_os << GetModel()->GetContinuousData()->GetCharacterName(k) << '\n';
            }
            cov_os << "u\n";
            cov_os << "Ne\n";
            cov_os << '\n';

            expandedmeancov->ToStream(cov_os);

            ofstream slope_os((GetName() + ".slopes").c_str());
            expandedmeancov->PrintSlopesNe(slope_os, idxdNdS, idxpiNpiS, idxpiS, idxNe);

            cerr << "covariance matrix in " << name << ".cov\n";
            cerr << "slopes of log dN/dS and log piN/piS ~ log Ne in " << name << ".slopes\n";
            cerr << '\n';
        }

		meanchrono->Normalise();
		ofstream chos((GetName() + ".postmeandates.tre").c_str());
		meanchrono->ToStream(chos);

		ofstream cchos((GetName() + ".postmeandates.tab").c_str());
		meanchrono->Tabulate(cchos);

		meanomega->Normalise();
		ofstream oos((GetName() + ".postmeanomega.tre").c_str());
		meanomega->ToStream(oos);
		cerr << "reconstructed variations of omega in " << name << ".postmeanomega.tre\n";

		meanNe->Normalise();
		ofstream Neos((GetName() + ".postmeanNe.tre").c_str());
		meanNe->ToStream(Neos);
		cerr << "reconstructed variations of Ne in " << name << ".postmeanNe.tre\n";

		meanu->Normalise();
		ofstream uos((GetName() + ".postmeanu.tre").c_str());
		meanu->ToStream(uos);
		cerr << "reconstructed variations of u in " << name << ".postmeanu.tre\n";

		meansynrate->Normalise();
		ofstream sos((GetName() + ".postmeansynrate.tre").c_str());
		meansynrate->ToStream(sos);
		cerr << "reconstructed variations of Ks in " << name << ".postmeansynrate.tre\n";

		for (int k=0; k<Ncont; k++)	{
			tree[k]->Normalise();
			ostringstream s;
			s << GetName() << ".postmean" << k+1 << ".tre";
			ofstream os(s.str().c_str());
			tree[k]->ToStream(os);
			cerr << "reconstructed variations of continuous character # " << k+1 << " in "  << name << ".postmean" << k+1 << ".tre\n";
		}

		ofstream ssos((GetName() + ".postmeansynrate.tab").c_str());
		meansynrate->Tabulate(ssos);
		ssos.close();

		ofstream ooos((GetName() + ".postmeanomega.tab").c_str());
		meanomega->Tabulate(ooos);
		ooos.close();

		ofstream NeNeos((GetName() + ".postmeanNe.tab").c_str());
		meanNe->Tabulate(NeNeos);
		NeNeos.close();

		ofstream uuos((GetName() + ".postmeanu.tab").c_str());
		meanu->Tabulate(uuos);
		uuos.close();

		for (int k=0; k<Ncont; k++)	{
			ostringstream s;
			s << GetName() << ".postmean" << k+1 << ".tab";
			ofstream os(s.str().c_str());
			tree[k]->Tabulate(os);
		}

		cerr << '\n';
	}	 
};


int main(int argc, char* argv[])	{

	int burnin = -1;
	int every = 1;
	int until = -1;
	string name;

	int fullsample = 0;

	bool check = false;
	string truefile = "";
	bool ic = false;

	int ppred = 0;

	bool printlog = false;
	bool printmean = false;
    bool printmed = false;
	bool printci = true;
	bool printstdev = false;
	bool withleaf = true;
	bool withinternal = true;

	bool withheader = true;
	double leafnameshift = 0.04;

	string mulreg = "";

	bool timeline = false;
	bool trend = false;

	bool strandsym = false;

	bool tex = false;
	double nodescale = 5;
	double nodepower = 1;
	double barwidth = 0.04;
	double x = 6;
	double y = 8;
	// bool withnodebubbles = true;
	// bool withtimeintervals = true;
	int fontsize = 4;

	bool bubbletext = false;

	double meanreg = 0;
	double stdevreg = 0;

	bool rr = false;

	bool nuc = false;
	string taxon = "";

	bool bf = false;
	bool kappa = false;
	bool normtest = false;

	bool densities = false;
	int cont = 0;
	string taxpairfile = "";

	bool postdist = false;
	bool Ne = false;

	try	{

		if (argc == 1)	{
			throw(0);
		}

		int i = 1;
		while (i < argc)	{
			string s = argv[i];

			if (s == "-pp")	{
				ppred = 1;
			}
			else if (s == "-nodata")	{
				fullsample = -1;
			}
			else if (s == "-strandsym")	{
				strandsym = true;
			}
			else if (s == "-dens")	{
				densities = true;
				i++;
				cont = atoi(argv[i]);
				i++;
				taxpairfile = argv[i];
			}
			else if (s == "-normtest")	{
				normtest = true;
			}
			else if (s == "-bf")	{
				bf = true;
			}
			else if (s == "-postdist")	{
				postdist = true;
			}
			else if (s == "-kappa")	{
				kappa = true;
			}
			else if (s == "-nuc")	{
				nuc = true;
			}
			else if (s == "-tax")	{
				i++;
				taxon = argv[i];
			}
			else if (s == "-tex")	{
				tex = true;
			}
			else if (s == "-header")	{
				withheader = false;
			}
			else if (s == "-shiftname")	{
				i++;
				leafnameshift = atof(argv[i]);
			}
			else if (s == "-ns")	{
				i++;
				nodescale = atof(argv[i]);
			}
			else if (s == "-np")	{
				i++;
				nodepower= atof(argv[i]);
			}
			else if (s == "-bw")	{
				i++;
				barwidth = atof(argv[i]);
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
			else if (s == "+bt")	{
				bubbletext = true;
			}
			else if (s == "-bt")	{
				bubbletext = false;
			}
			else if (s == "-pointwise")	{
				i++;
				meanreg = atof(argv[i]);
				i++;
				stdevreg = atof(argv[i]);
			}
			/*
			else if (s == "+nb")	{
				withnodebubbles = true;
			}
			else if (s == "-nb")	{
				withnodebubbles = false;
			}
			*/
			else if (s == "-ic")	{
				ic = true;
			}
			else if (s == "-timeline")	{
				timeline = true;
			}
			else if (s == "-trend")	{
				trend = true;
			}
			else if (s == "-ch")	{
				check = true;
				i++;
				truefile = argv[i];
			}
			else if (s == "+log")	{
				printlog = true;
			}
			else if (s == "-log")	{
				printlog = false;
			}
			else if (s == "+mean")	{
				printmean = true;
                printmed = false;
			}
			else if (s == "-mean")	{
				printmean = false;
			}
			else if ((s == "+med") || (s == "+median"))	{
				printmed = true;
                printmean = false;
			}
			else if ((s == "-med") || (s == "-median"))	{
				printmed = false;
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
			else if (s == "-rr")	{
				rr = true;
			}
			else if (s == "-Ne") {
				Ne = true;
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

	BranchOmegaMultivariateSample sample(name,burnin,every,until,fullsample);

	if (densities)	{
		sample.DrawDensities(taxpairfile,cont);
		exit(1);
	}
	if (normtest)	{
		sample.PostPredNormalityTest();
		exit(1);
	}
	if (nuc)	{
		sample.ReadNuc(taxon);
		exit(1);
	}
	if (Ne) {
		sample.ReadNe(printlog,printmean,printmed,printci,printstdev,withleaf,withinternal,meanreg,stdevreg);
		exit(1);
	}	
	if (rr)	{
		sample.ReadRelRates();
		exit(1);
	}

	if (trend)	{
		sample.PostPredTrend();
		exit(1);
	}

	if (timeline)	{
		if (sample.GetModel()->withtimeline)	{
			sample.ReadTimeLine();
		}
		else	{
			sample.MakeTimeLine();
		}
		exit(1);
	}
	if (check)	{
		cerr << "check cov\n";
		sample.CheckCov(truefile);
		exit(1);
	}
	if (ppred)	{
		sample.PostPred();
		exit(1);
	}
	if (ic)	{
		sample.ReadIC();
		exit(1);
	}

	sample.Read(printlog,printmean,printmed,printci,printstdev,withleaf,withinternal,mulreg,tex,x,y,nodescale,nodepower,barwidth,fontsize,bubbletext,withheader,leafnameshift,meanreg,stdevreg,postdist);

}



