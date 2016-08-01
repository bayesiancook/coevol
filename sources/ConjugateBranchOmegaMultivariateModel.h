
#include "BranchOmegaMultivariateModel.h"
#include "ConjugateMultiVariateTreeProcess.h"

class ConjugateBranchOmegaMultivariateModel : public BranchOmegaMultivariateModel {

	public:

	ConjugateBranchOmegaMultivariateModel(string datafile, string treefile, string contdatafile, string calibfile, double rootage, double rootstdev, int inchronoprior, double insofta, double inmeanchi, double inmeanchi2, double priorsigma, string priorsigmafile, int indf, int inmutmodel, int ingc, bool inautoregressive, int inconjpath, double inmappingfreq, int contdatatype, int inomegaratiotree, bool inclamproot, bool inclamptree, bool inmeanexp, bool innormalise, int innrep, int inncycle, string inbounds, string inmix, int inNinterpol, int inwithdrift, int inuniformprior, string rootfile, string insuffstatfile, bool intimeline, bool inseparatesyn, bool inseparateomega, int inkrkctype, int injitter, int inmyid, int innprocs, int insample, GeneticCodeType type)	{

		sample = insample;

		nprocs = innprocs;
		myid = inmyid;

		dsindex = -1;
		whitenoise = false;

		contjitter = 0;
		if (injitter > 2)	{
			injitter -= 2;
			contjitter = 1;
		}
		jitter = injitter;

		krkctype = inkrkctype;
		splitkrkctype = 0;
		if (krkctype >= 4)	{
			splitkrkctype = 1;
			krkctype -= 4;
			krkctype = inkrkctype;
		}
		codetype = type;
		codonstatespace = new CodonStateSpace(codetype);

		separatesyn = inseparatesyn;
		separateomega = inseparateomega;

		withtimeline = intimeline;

		suffstatfile = insuffstatfile;
		clampsuffstat = (suffstatfile != "None");

		maxphi = 10.0;
		maxdrift = 100.0;
		uniformprior = inuniformprior;

		if (inwithdrift == 4)	{
			withdrift = false;
			withexpdrift = false;
			withreldrift = true;
		}
		else if (inwithdrift == 3)	{
			withdrift = true;
			withexpdrift = true;
			withreldrift = true;
		}
		else if (inwithdrift == 2)	{
			withdrift = true;
			withexpdrift = true;
			withreldrift = false;
		}
		else if (inwithdrift == 1)	{
			withdrift = true;
			withexpdrift = false;
			withreldrift = false;
		}
		else {
			withdrift = false;
			withexpdrift = false;
			withreldrift = false;
		}

		Ninterpol = inNinterpol;
		bounds = inbounds;
		mix = inmix;

		df = indf;
		mutmodel = inmutmodel;

		normalise = innormalise;

		chronoprior = inchronoprior;
		softa = insofta;
		meanchi = inmeanchi;
		meanchi2 = inmeanchi2;

		autoregressive = inautoregressive;
		clampdiag = false;
		clamproot = inclamproot;
		clamptree = inclamptree;
		meanexp = inmeanexp;

		// A FROZEN ACCIDENT...
		if (inomegaratiotree == 3)	{
			omegaratiotree = 4;
		}
		else if (inomegaratiotree == 4)	{
			omegaratiotree = 3;
		}
		else	{
			omegaratiotree = inomegaratiotree;
		}

		if (withtimeline && ! clamptree)	{
			cerr << "error : should clamp the tree to use timeline\n";
			exit(1);
		}

		gc = ingc;
		if (gc == 2)	{
			gc = 1;
			whitenoise = true;
		}

		if (gc == 3)	{
			if (! omegaratiotree)	{
				cerr << "error : can apply gc3 formalism only under codon model\n";
				cerr << "to activate codon model, use -dSdN or -dSomega options\n";
				exit(1);
			}
		}
		gcstat = 0;
		if (gc == -1)	{
			gc = 0;
			gcstat = 1;
		}


		// 4 : 3 independent processes : ts tv0 tvgc
		// OLD CODE
		// 5 : 2 independent processes : mixture of the log
		// 6 : 2 independent processes : mixture of the exp

		// NEW CODE
		// 5 : ts and tv0 but only one omega
		if (omegaratiotree == 5) 	{
			L = 3 + gc;
		}
		else if (omegaratiotree == 4) 	{
			L = 6 + gc;
		}
		else if (omegaratiotree == 3)	{
			// transition transversion
			L = 4 + gc;
		}
		else if (omegaratiotree)	{
			L = 2 + gc;
		}
		else	{
			L = 1 + gc;
		}

		if (mutmodel == 5)	{
			L += 2;
		}
		else if (mutmodel == 4)	{
			L += 1;
		}
		else if (mutmodel == 3)	{
			L += 3;
		}
		else if (mutmodel)	{
			L++;
		}

		if (separatesyn)	{
			L--;
		}
		if (separateomega)	{
			L--;
		}

		priorsampling = false;

		if (inconjpath == -1)	{
			conjpath = 2;
		}
		else if (inconjpath == 3)	{
			conjpath = 0;
			priorsampling = true;
		}
		else	{
			conjpath = inconjpath;
		}

		nrep = innrep;
		if (nrep == 0)	{
			nrep = conjpath ? 30 : 1;
		}
		ncycle = inncycle;

		matrixtype = 0;
		mappingfreq = inmappingfreq;

		if (clampsuffstat && ! conjpath)	{
			cerr << "error : suffstat requires conjugate path sampling\n";
			exit(1);
		}

		// get data from file
		nucdata = new FileSequenceAlignment(datafile);

		rawNstate = nucdata->GetNstate();

		if (rawNstate == Naa)	{
			if ((! omegaratiotree) || (krkctype == -1))	{
				cerr << "error : with amino acid data, should specify a krkc model\n";
				exit(1);
			}
		}

		if (omegaratiotree && (rawNstate == Nnuc))	{
			codondata = new CodonSequenceAlignment(nucdata, true, type);
		}
		else	{
			codondata = 0;
		}

		Nsite = GetData()->GetNsite();

		if (GetNprocs() > 1)	{
			MakeMPI();
		}

		taxonset = nucdata->GetTaxonSet();

		if (GetMyid())	{
			nucdata->SubSelect(GetSiteMin(),GetSiteMax());
			if (codondata)	{
				codondata->SubSelect(GetSiteMin(),GetSiteMax());
			}
		}

		// get tree from file (newick format)
		tree = new Tree(treefile);
		// check whether tree and data fit together
		tree->RegisterWith(taxonset);

		if (Split())	{
			if (withexpdrift || withreldrift)	{
				cerr << "error : split and exp drift not yet compatible\n";
				exit(1);
			}
			// tree->Subdivide(tree->GetRoot(),Ninterpol);
			splittree = new SplitTree(tree,Ninterpol);
		}
		else	{
			splittree = tree;
		}

		// get continuous data from file
		if (contdatafile != "None")	{
			contdata = new FileContinuousData(contdatafile);
			Ncont = contdata->GetNsite();
		}
		else	{
			contdata = 0;
			Ncont = 0;
		}

		cerr << "tree and data ok\n";
		cerr << '\n';

		// ----------
		// construction of the graph
		// ----------

		Zero = new Const<Real>(0);
		One = new Const<PosReal>(1);

		MeanChi = 0;
		MeanChi2 = 0;
		Chi = 0;
		Chi2 = 0;

		syngammatree = 0;
		iscalib = false;

		if (calibfile == "Unconstrained")	{
			if (withexpdrift || withreldrift)	{
				cerr << "error : unconstrained and exp drift not compatible\n";
				exit(1);
			}
			L--;
			PriorMu = new Const<PosReal>(0.1);
			mu = new Gamma(One,PriorMu);
			syngammatree = new GammaTree(tree,One,mu);
			lengthtree = syngammatree;
			synratetree = syngammatree;
		}
		else	{
			PriorMu = new Const<PosReal>(1);
			mu = new Gamma(One,PriorMu);
			mu->ClampAt(1);

			if (calibfile != "None")	{
				iscalib = true;

				double a = rootage * rootage / rootstdev / rootstdev;
				double b = rootage / rootstdev / rootstdev;
				if (rootage == -1)	{
					a = b = -1;
				}
				CalibrationSet* calibset = new FileCalibrationSet(calibfile, tree);

				if (chronoprior == 0)	{
					chronogram = new CalibratedChronogram(tree,mu,a,b,calibset);
				}
				else	{
					if (meanchi != -1)	{
						MeanChi = new Const<PosReal>(meanchi);
						MeanChi2 = new Const<PosReal>(meanchi2);
						Chi = new Exponential(MeanChi,Exponential::MEAN);
						Chi2 = new Exponential(MeanChi2,Exponential::MEAN);
					}
					else	{
						double min = 1e-6;
						double max = 1e6;
						Chi = new Jeffreys(min,max,Zero);
						Chi2 = new Jeffreys(min,max,Zero);
					}
					chronogram = new BDCalibratedChronogram(tree,mu,Chi,Chi2,a,b,calibset,chronoprior,softa);
				}
			}
			else	{
				if (chronoprior == 0)	{
					chronogram = new Chronogram(tree,mu);
				}
				else {
					if (meanchi != -1)	{
						MeanChi = new Const<PosReal>(meanchi);
						MeanChi2 = new Const<PosReal>(meanchi2);
						Chi = new Exponential(MeanChi,Exponential::MEAN);
						Chi2 = new Exponential(MeanChi2,Exponential::MEAN);
					}
					else	{
						double min = 1e-6;
						double max = 1e6;
						Chi = new Jeffreys(min,max,Zero);
						Chi2 = new Jeffreys(min,max,Zero);
					}
					chronogram = new BDChronogram(tree,mu,Chi,Chi2);
				}
			}

			if (clamptree)	{
				chronogram->Clamp();
			}

			if (Split())	{
				lengthtree = new SplitLengthTree(chronogram,GetSplitTree());
			}
			else	{
				lengthtree = chronogram;
			}
		}

		if (mix == "branch")	{
			MixAlpha = new Jeffreys(0.1,10,Zero);
			MixAlpha->setval(1.0);
			gammatree = new GammaTree(lengthtree->GetTree(),MixAlpha,MixAlpha);
			// gammatree->GetBranchVal(gammatree->GetRoot()->Next()->Out()->GetBranch())->ClampAt(1.0);
			partition = new BranchPartition(tree);
			gammamixtree = 0;
		}
		else if (mix != "None")	{
			partition = new BranchPartition(tree,mix);
			if (Split())	{
				splitpartition = new SplitBranchPartition(partition,GetSplitTree());
			}
			gammamixtree = 0;
			// gammamixtree = new GammaMixTree(partition,One,One,true); // clamp first component to 1
			gammatree = 0;
		}
		else	{
			partition = new BranchPartition(tree);
			if (Split())	{
				splitpartition = new SplitBranchPartition(partition,GetSplitTree());
			}
			gammamixtree = 0;
			gammatree = 0;
		}

		Nmat = GetPartition()->GetNComponent();

		double mindiag = 0.001;
		double maxdiag = 1000;
		DiagArray = new JeffreysIIDArray(Ncont+L,mindiag,maxdiag,Zero);
		if (priorsigma == -2)	{
			ifstream is(priorsigmafile.c_str());
			for (int i=0; i<Ncont+L; i++)	{
				double tmp;
				is >> tmp;
				DiagArray->ClampAt(tmp,i);
			}
		}
		else if (priorsigma == -1)	{
			DiagArray->setval(1.0);
		}
		else	{
			DiagArray->ClampAt(priorsigma);
		}
		DiagArray0 = new JeffreysIIDArray(Ncont+L,mindiag,maxdiag,Zero);
		if (priorsigma == -2)	{
			ifstream is(priorsigmafile.c_str());
			for (int i=0; i<Ncont+L; i++)	{
				double tmp;
				is >> tmp;
				DiagArray0->ClampAt(tmp,i);
			}
		}
		else if (priorsigma == -1)	{
			DiagArray0->ClampAt(1.0);
		}
		else	{
			DiagArray0->ClampAt(priorsigma);
		}

		sigmaZero = new SigmaZero(DiagArray);
		sigmaZero0 = new SigmaZero(DiagArray0);

		sigmaarray = new ConjugateInverseWishartArray(Nmat,df,sigmaZero,sigmaZero0);
		sigma = sigmaarray->GetVal(0);
		// sigma->SetIdentity();

		// driftarray = new MultiVarArray(Nmat, Zero, ContDiagArray, ContDiagArray0);
		if (uniformprior)	{
			MultiUniArray* tmparray = new MultiUniArray(Nmat, Ncont+L, Zero, maxdrift);
			driftarray = tmparray;
			if (! withdrift)	{
				tmparray->ClampAtZero();
			}
		}
		else	{
			MultiVarArray* tmparray = new MultiVarArray(Nmat, Zero, DiagArray, DiagArray0);
			driftarray = tmparray;
			if (! withdrift)	{
				tmparray->ClampAtZero();
			}
		}
		drift = driftarray->GetVal(0);

		driftphiarray = 0;
		driftphi = 0;
		driftarray2 = 0;
		drift2 = 0;
		driftphiarray2 = 0;
		driftphi2 = 0;

		if (withexpdrift)	{
			if (uniformprior)	{
				driftphiarray = new PosUniIIDArray(Nmat,One,maxphi);
			}
			else	{
				driftphiarray = new GammaIIDArray(Nmat,One,One);
			}
			driftphi = driftphiarray->GetVal(0);
		}

		if (withreldrift)	{
			if (uniformprior)	{
				driftarray2 = new MultiUniArray(Nmat, Ncont+L, Zero, maxdrift);
				driftphiarray2 = new PosUniIIDArray(Nmat,One,maxphi);
			}
			else	{
				driftarray2 = new MultiVarArray(Nmat, Zero, DiagArray, DiagArray0);
				driftphiarray2 = new GammaIIDArray(Nmat,One,One);
			}
			drift2 = driftarray2->GetVal(0);
			driftphi2 = driftphiarray2->GetVal(0);
		}


		rootmean = 0;
		rootvar = 0;
		if (rootfile != "None")	{
			rootmean = new Const<RealVector>(RealVector(Ncont+L));
			rootvar = new Const<PosRealVector>(PosRealVector(Ncont+L));
			ifstream is(rootfile.c_str());
			if (! is)	{
				cerr << "error: cannot open root file " << rootfile << '\n';
				exit(1);
			}
			cerr << "root constraints \n";
			for (int i=0; i<Ncont + L; i++)	{
				(*rootmean)[i] = 0;
				(*rootvar)[i] = 0;
			}
			int  n;
			is >> n;
			if (n != Ncont)	{
				cerr << "error : number of entries in root calibration file " << rootfile << " does not match number of quantitative traits (" << Ncont << ")\n";
				exit(1);
			}
			for (int i=0; i<Ncont; i++)	{
				double mean, stdev;
				is >> mean >> stdev;
				(*rootmean)[i+L] = mean;
				(*rootvar)[i+L] = stdev * stdev;
				cerr << GetContinuousData()->GetCharacterName(i) << '\t' << mean << '\t' << stdev << '\n';
			}
		}

		if (withtimeline)	{
			if (Unconstrained())	{
				cerr << "error : timeline and unconstrained incompatible\n";
				exit(1);
			}
			TimeLineDiagArray = new JeffreysIIDArray(Ncont+L,mindiag,maxdiag,Zero);
			TimeLineDiagArray->ClampAt(1.0);
			TimeLineSigmaZero = new SigmaZero(TimeLineDiagArray);
			timelinesigma = new DiagonalCovMatrix(TimeLineSigmaZero, Ncont+L+df);
			timelinesigma->SetIdentity();

			timeintervals = new TimeIntervals(chronogram);
			timeline = new TimeLine(chronogram,timeintervals, timelinesigma);
		}

		ugam = 0;
		ugamtree = 0;
		wngamtree = 0;
		if (jitter == 1)	{
			ugam = new Gamma(One,One);
			ugam->setval(10.0);
			ugamtree = new GammaTree(tree,ugam,ugam);
		}
		if (jitter == 2)	{
			ugam = new Gamma(One,One);
			ugam->setval(10.0);
			wngamtree = new GammaWhiteNoiseProcess(lengthtree,ugam);
			wngamtree->Reset();
		}

		synsigma = 0;
		omegasigma = 0;
		lognormalsyntree = 0;
		lognormalomegatree = 0;
		if (separatesyn)	{
			synsigma = new Jeffreys(mindiag,maxdiag,Zero);
			synsigma->setval(1);
			lognormalsyntree = new LogNormalTreeProcess(lengthtree,synsigma,INTEGRAL);
			lognormalsyntree->Reset(-1.0);
		}
		if (separateomega)	{
			omegasigma = new Jeffreys(mindiag,maxdiag,Zero);
			omegasigma->setval(1);
			lognormalomegatree = new LogNormalTreeProcess(lengthtree,omegasigma,MEAN);
			lognormalomegatree->Reset(-2.0);
		}

		if (autoregressive)	{
			phi = new Gamma(One,One);
			mean = new IIDUniform(Zero,Ncont+L,100);
			if (gammamixtree)	{
				cerr << "error: cannot use partitions for autocorrelated processes\n";
				exit(1);
			}
			process = new ConjugateAutoRegressiveMultiVariateTreeProcess(GetConjugateInverseWishart(0),mean,phi,lengthtree);
		}
		else	{
			phi = 0;
			mean = 0;
			if (gammatree)	{
				// process = new ConjugatePartitionMultiVariateTreeProcess(sigmaarray,lengthtree,GetPartition(),gammatree,driftarray, driftphiarray, chronogram);
				if (withexpdrift)	{
					process = new ConjugatePartitionMultiVariateTreeProcess(sigmaarray,lengthtree,GetPartition(),gammatree,driftarray,driftphiarray,chronogram, rootmean, rootvar, driftarray2, driftphiarray2, GetScale(), 65);
				}
				else	{
				}
			}
			else if (withtimeline)	{
				process = new ConjugateTimeLineMultiVariateTreeProcess(GetConjugateInverseWishart(0),timeline,chronogram);
			}
			else	{
				// process = new ConjugatePartitionMultiVariateTreeProcess(sigmaarray,lengthtree,GetPartition(),0,driftarray, driftphiarray, chronogram);
				process = new ConjugatePartitionMultiVariateTreeProcess(sigmaarray,lengthtree,GetPartition(),gammamixtree,driftarray,driftphiarray,chronogram, rootmean, rootvar, driftarray2, driftphiarray2, GetScale(), 65);
			}
		}

		/*
		if (Mix2Omega())	{
			process->Reset();
		}
		*/
		cerr << "process RESET\n";
		process->Reset();

		if (clamproot)	{
			process->ClampRoot();
		}

		if (whitenoise)	{
			cerr << "white noise\n";
			wnvar = new Gamma(One,One);
			cerr << "white noise tree\n";
			wntree = new WhiteNoiseProcess(lengthtree,Zero,wnvar);
			cerr << "ok\n";
		}
		else	{
			wnvar = 0;
			wntree = 0;
		}

		leafstates = 0;
		leafvar = 0;
		if (contjitter && (! contdata))	{
			cerr << "error: should have continuous data with contjitter option\n";
			exit(1);
		}

		if (contdata)	{
			if (contjitter)	{
				if (contdatatype)	{
					cerr << "error: model with leaf polymorphism only with log transformation\n";
					exit(1);
				}
				leafvar = new Gamma(One,One);
				leafvar->setval(0.1);
				leafstates = new Normal**[GetNtaxa()];
				for (int i=0; i<GetNtaxa(); i++)	{
					leafstates[i] = new Normal*[Ncont];
					for (int j=0; j<Ncont; j++)	{
						leafstates[i][j] = 0;
					}
				}
				for (int i=0; i<Ncont; i++)	{
					process->SetLeafStates(leafstates,contdata,leafvar,L+i,i);
				}
			}
			else	{
				for (int i=0; i<Ncont; i++)	{
					process->SetAndClamp(contdata,L+i,i,contdatatype);
				}
			}
		}


		for (int l=0; l<L; l++)	{
			process->CutOff(1,l);
		}

		if (bounds != "None")	{
			FileBoundSet boundset(bounds,process->GetTree());
			process->SetBounds(boundset,L);
		}

		CreateSubstitutionProcess(sample);

		if (phyloprocess)	{
			phyloprocess->Unfold();
		}
		if (sample > 0)	{
			if (phyloprocess)	{
				phyloprocess->Sample();
			}
		}


		RootRegister(PriorMu);
		RootRegister(Zero);
		RootRegister(One);
		if (MeanChi)	{
			RootRegister(MeanChi);
			RootRegister(MeanChi2);
		}
		if (rawNstate == Naa)	{
			RootRegister(exprelrate);
		}
		else if ((mutmodel == 0) || (mutmodel >= 4))	{
			RootRegister(exprelrate);
		}
		else if (mutmodel == 2)	{
			RootRegister(tsrelrate);
			RootRegister(tvrelrate);
		}

		if (freestationary)	{
			RootRegister(freestationary);
		}
		if (rootmean)	{
			RootRegister(rootmean);
			RootRegister(rootvar);
		}
		Register();

		MakeScheduler();
		PreUpdate();
	}

	~ConjugateBranchOmegaMultivariateModel() {}

	ConjugateInverseWishart* GetConjugateInverseWishart(int mat) {
		ConjugateInverseWishart* tmp = dynamic_cast<ConjugateInverseWishart*>(sigmaarray->GetVal(mat));
		if (! tmp)	{
			cerr << "error : dynamic castt of conjugate inverse wishart : " << sigma << '\t' << tmp << '\n';
			exit(1);
		}
		return tmp;
	}

	ConjugatePartitionMultiVariateTreeProcess* GetConjugatePartitionMultiVariateTreeProcess() {
		ConjugatePartitionMultiVariateTreeProcess* tmp = dynamic_cast<ConjugatePartitionMultiVariateTreeProcess*>(process);
		if (! tmp)	{
			cerr << "error : dynamic cast of partition multivariate tree process : " << process << '\t' << tmp << '\n';
			exit(1);
		}
		return tmp;
	}

	ConjugateMultiVariateTreeProcess* GetConjugateMultiVariateTreeProcess() {
		ConjugateMultiVariateTreeProcess* tmp = dynamic_cast<ConjugateMultiVariateTreeProcess*>(process);
		if (! tmp)	{
			cerr << "error : dynamic cast of multivariate tree process : " << process << '\t' << tmp << '\n';
			exit(1);
		}
		return tmp;
	}

	void MakeScheduler()	{

		if (sample != -1)	{
		if (conjpath)	{
			if (! clampsuffstat)	{
				if ((GetNprocs() == 1) || (GetMyid() > 0))	{
					scheduler.Register(new DSemiConjugateMappingMove(phyloprocess,pathconjtree,mappingfreq),1,"mapping + sufficient stat");
				}
			}
		}
		else	{
			if (phyloprocess)	{
				scheduler.Register(new SimpleMove(phyloprocess,1),1,"mapping");
			}
		}
		}

		if (! GetMyid())	{

		vector <SplitMultiVariateNodeMove*> nodesplitarray;
		vector <SplitMultiVariateBranchMove*> branchsplitarray;
		if (Split())	{
			nodesplitarray.push_back(new SplitMultiVariateNodeMove(process,10));
			nodesplitarray.push_back(new SplitMultiVariateNodeMove(process,1));
			nodesplitarray.push_back(new SplitMultiVariateNodeMove(process,0.1));
			nodesplitarray.push_back(new SplitMultiVariateNodeMove(process,0.01));
			branchsplitarray.push_back(new SplitMultiVariateBranchMove(tree,process,10));
			branchsplitarray.push_back(new SplitMultiVariateBranchMove(tree,process,1));
			branchsplitarray.push_back(new SplitMultiVariateBranchMove(tree,process,0.1));
			branchsplitarray.push_back(new SplitMultiVariateBranchMove(tree,process,0.01));
		}

		for (int i=0; i<nrep; i++)	{
			if (Unconstrained())	{
				scheduler.Register(new SimpleMove(mu,1),10,"syngamtree hyper");
				scheduler.Register(new SimpleMove(mu,0.1),10,"syngamtree hyper");
				scheduler.Register(new SimpleMove(syngammatree,3),30,"syngamtree");
				scheduler.Register(new SimpleMove(syngammatree,0.3),30,"syngamtree");
				scheduler.Register(new SimpleMove(syngammatree,0.03),30,"syngamtree");
			}
			else if (! clamptree)	{
				if (chronoprior)	{
					scheduler.Register(new SimpleMove(Chi,1),10,"bd hyper");
					scheduler.Register(new SimpleMove(Chi,0.1),10,"bd hyper");
					scheduler.Register(new SimpleMove(Chi2,1),10,"bd hyper");
					scheduler.Register(new SimpleMove(Chi2,0.1),10,"bd hyper");
				}
				scheduler.Register(new SimpleMove(chronogram,1),50,"chrono");
				scheduler.Register(new SimpleMove(chronogram,0.1),50,"chrono");
				scheduler.Register(new SimpleMove(chronogram,0.01),50,"chrono");

				/*
				if (separatesyn)	{
					scheduler.Register(new LogNormalChronoCompMove(chronogram,lognormalsyntree,1),10,"chrono lognormalsynrate comp");
					scheduler.Register(new LogNormalChronoCompMove(chronogram,lognormalsyntree,0.1),10,"chrono lognormalsynrate comp");
					scheduler.Register(new LogNormalChronoCompMove(chronogram,lognormalsyntree,0.01),10,"chrono lognormalsynrate comp");
				}
				else	{
					scheduler.Register(new MultiVariateChronoCompMove(chronogram,process,dsindex,1),10,"chrono process comp");
					scheduler.Register(new MultiVariateChronoCompMove(chronogram,process,dsindex,0.1),10,"chrono process comp");
					scheduler.Register(new MultiVariateChronoCompMove(chronogram,process,dsindex,0.01),10,"chrono process comp");
				}
				*/
			}

			if (gammamixtree)	{
				scheduler.Register(new SimpleMove(gammamixtree,10),10,"gammamixtree");
				scheduler.Register(new SimpleMove(gammamixtree,1),10,"gammamixtree");
				scheduler.Register(new SimpleMove(gammamixtree,0.1),10,"gammamixtree");
			}
			if (gammatree)	{
				scheduler.Register(new SimpleMove(MixAlpha,10),10,"mixalpha");
				scheduler.Register(new SimpleMove(MixAlpha,1),10,"mixalpha");
				scheduler.Register(new SimpleMove(MixAlpha,0.1),10,"mixalpha");
				scheduler.Register(new SimpleMove(gammatree,10),10,"gammatree");
				scheduler.Register(new SimpleMove(gammatree,1),10,"gammatree");
				scheduler.Register(new SimpleMove(gammatree,0.1),10,"gammatree");
			}
			if (whitenoise)	{
				scheduler.Register(new SimpleMove(wntree,10),10,"gammatree");
				scheduler.Register(new SimpleMove(wntree,1),10,"gammatree");
				scheduler.Register(new SimpleMove(wntree,0.1),10,"gammatree");
				scheduler.Register(new SimpleMove(wnvar,10),10,"gammatree");
				scheduler.Register(new SimpleMove(wnvar,1),10,"gammatree");
				scheduler.Register(new SimpleMove(wnvar,0.1),10,"gammatree");
			}

			if (isCalibrated())	{
				/*
				scheduler.Register(new SimpleMove(GetCalibratedChronogram()->GetScale(),1),1,"root age");
				scheduler.Register(new SimpleMove(GetCalibratedChronogram()->GetScale(),0.1),1,"root age");
				scheduler.Register(new SimpleMove(GetCalibratedChronogram()->GetScale(),0.01),1,"root age");
				scheduler.Register(new SimpleMove(GetCalibratedChronogram()->GetScale(),0.001),1,"root age");
				*/

				scheduler.Register(new SimpleMove(GetCalibratedChronogram()->GetScale(),1),100,"root age");
				scheduler.Register(new SimpleMove(GetCalibratedChronogram()->GetScale(),0.1),100,"root age");
				scheduler.Register(new SimpleMove(GetCalibratedChronogram()->GetScale(),0.01),100,"root age");
				scheduler.Register(new SimpleMove(GetCalibratedChronogram()->GetScale(),0.001),100,"root age");

				/*
				scheduler.Register(new RootMove(GetCalibratedChronogram(),1),10,"root age additive");
				scheduler.Register(new RootMove(GetCalibratedChronogram(),0.1),10,"root age additive");
				scheduler.Register(new RootMove(GetCalibratedChronogram(),0.01),10,"root age additive");
				*/

				/*
				scheduler.Register(new MultiVariateScaleCompMove(GetCalibratedChronogram()->GetScale(),process,dsindex,1),10,"chrono scale process comp");
				scheduler.Register(new MultiVariateScaleCompMove(GetCalibratedChronogram()->GetScale(),process,dsindex,0.1),10,"chrono scale process comp");
				scheduler.Register(new MultiVariateScaleCompMove(GetCalibratedChronogram()->GetScale(),process,dsindex,0.1),10,"chrono scale process comp");


				scheduler.Register(new RootProcessMove(GetCalibratedChronogram(),process,dsindex,1),10,"root age additive");
				scheduler.Register(new RootProcessMove(GetCalibratedChronogram(),process,dsindex,0.1),10,"root age additive");
				scheduler.Register(new RootProcessMove(GetCalibratedChronogram(),process,dsindex,0.01),10,"root age additive");

				scheduler.Register(new RootSigmaMove(GetCalibratedChronogram(),GetConjugateInverseWishart(0),1),10,"root age sigma additive");
				scheduler.Register(new RootSigmaMove(GetCalibratedChronogram(),GetConjugateInverseWishart(0),0.1),10,"root age sigma additive");
				scheduler.Register(new RootSigmaMove(GetCalibratedChronogram(),GetConjugateInverseWishart(0),0.01),10,"root age sigma additive");

				scheduler.Register(new RootSigmaProcessMove(GetCalibratedChronogram(),GetConjugateInverseWishart(0),process,dsindex,1),10,"root age sigma additive");
				scheduler.Register(new RootSigmaProcessMove(GetCalibratedChronogram(),GetConjugateInverseWishart(0),process,dsindex,0.1),10,"root age sigma additive");
				scheduler.Register(new RootSigmaProcessMove(GetCalibratedChronogram(),GetConjugateInverseWishart(0),process,dsindex,0.01),10,"root age sigma additive");
				*/

				/*
				if (separatesyn)	{
					scheduler.Register(new LogNormalScaleCompMove(GetCalibratedChronogram()->GetScale(),lognormalsyntree,1),10,"chrono scale lognormalsynrate comp");
					scheduler.Register(new LogNormalScaleCompMove(GetCalibratedChronogram()->GetScale(),lognormalsyntree,0.1),10,"chrono scale lognormalsynrate comp");
					scheduler.Register(new LogNormalScaleCompMove(GetCalibratedChronogram()->GetScale(),lognormalsyntree,0.01),10,"chrono scale lognormalsynrate comp");
				}
				else	{
					scheduler.Register(new MultiVariateScaleCompMove(GetCalibratedChronogram()->GetScale(),process,dsindex,1),10,"chrono scale process comp");
					scheduler.Register(new MultiVariateScaleCompMove(GetCalibratedChronogram()->GetScale(),process,dsindex,0.1),10,"chrono scale process comp");
					scheduler.Register(new MultiVariateScaleCompMove(GetCalibratedChronogram()->GetScale(),process,dsindex,0.1),10,"chrono scale process comp");
				}
				*/
			}

			if (jitter == 1)	{
				scheduler.Register(new SimpleMove(ugamtree,10),10,"ugamtree");
				scheduler.Register(new SimpleMove(ugamtree,1),10,"ugamtree");
				scheduler.Register(new SimpleMove(ugamtree,0.1),10,"ugamtree");
				scheduler.Register(new SimpleMove(ugam,10),10,"ugam");
				scheduler.Register(new SimpleMove(ugam,1),10,"ugam");
				scheduler.Register(new SimpleMove(ugam,0.1),10,"ugam");
				scheduler.Register(new MultiVariateUgamCompMove(ugamtree,process,omegaindex,1),10,"ugam process comp");
				scheduler.Register(new MultiVariateUgamCompMove(ugamtree,process,omegaindex,0.1),10,"ugam process comp");
				scheduler.Register(new MultiVariateUgamCompMove(ugamtree,process,omegaindex,0.01),10,"ugam process comp");
			}

			if (jitter == 2)	{
				scheduler.Register(new SimpleMove(wngamtree,10),10,"wntree");
				scheduler.Register(new SimpleMove(wngamtree,1),10,"wntree");
				scheduler.Register(new SimpleMove(wngamtree,0.1),10,"wntree");
				scheduler.Register(new SimpleMove(ugam,10),10,"wn");
				scheduler.Register(new SimpleMove(ugam,1),10,"wn");
				scheduler.Register(new SimpleMove(ugam,0.1),10,"wn");
				scheduler.Register(new MultiVariateUgamCompMove(wngamtree,process,omegaindex,1),10,"wnugam process comp");
				scheduler.Register(new MultiVariateUgamCompMove(wngamtree,process,omegaindex,0.1),10,"wnugam process comp");
				scheduler.Register(new MultiVariateUgamCompMove(wngamtree,process,omegaindex,0.01),10,"wnugam process comp");
			}

			if (contjitter)	{
				scheduler.Register(new SimpleMove(leafvar,1),10,"leafvar");
				scheduler.Register(new SimpleMove(leafvar,0.1),10,"leafvar");
				scheduler.Register(new SimpleMove(leafvar,0.01),10,"leafvar");
			}

			if (bounds == "None")	{
				int n = taxonset->GetNtaxa() * 100;
				scheduler.Register(new MultiVariatePropagateMove(process,1,0.1,0.1),n,"propmove");
				scheduler.Register(new MultiVariatePropagateMove(process,0.1,0.5,0.5),n,"propmove");
				scheduler.Register(new MultiVariatePropagateMove(process,0.1,0.9,0.9),n,"propmove");
				scheduler.Register(new MultiVariatePropagateMove(process,0.01,0.99,0.99),n,"propmove");
			}

			if (Split())	{

				scheduler.Register(nodesplitarray[0],10,"split node process");
				scheduler.Register(nodesplitarray[1],10,"split node process");
				scheduler.Register(nodesplitarray[2],10,"split node process");
				scheduler.Register(nodesplitarray[3],10,"split node process");
				scheduler.Register(branchsplitarray[0],10,"split node process");
				scheduler.Register(branchsplitarray[1],10,"split node process");
				scheduler.Register(branchsplitarray[2],10,"split node process");
				scheduler.Register(branchsplitarray[3],10,"split node process");
			}

			// scheduler.Register(new SimpleMove(process,10),30,"multinormal");
			scheduler.Register(new SimpleMove(process,1),150,"multinormal");
			scheduler.Register(new SimpleMove(process,0.1),150,"multinormal");
			scheduler.Register(new SimpleMove(process,0.01),150,"multinormal");

			// scheduler.Register(new ConjugateMultiVariateMove(GetConjugateInverseWishart(),GetConjugateMultiVariateTreeProcess(),1,0),1,"conjugate sigma - process");

			/*
			scheduler.Register(new ConjugatePartitionMultiVariateMove(GetConjugatePartitionMultiVariateTreeProcess(),10,10),1,"conjugate sigma - process");
			scheduler.Register(new ConjugatePartitionMultiVariateMove(GetConjugatePartitionMultiVariateTreeProcess(),1,10),1,"conjugate sigma - process");
			scheduler.Register(new ConjugatePartitionMultiVariateMove(GetConjugatePartitionMultiVariateTreeProcess(),0.1,10),1,"conjugate sigma - process");
			scheduler.Register(new ConjugatePartitionMultiVariateMove(GetConjugatePartitionMultiVariateTreeProcess(),0.01,10),1,"conjugate sigma - process");
			scheduler.Register(new ConjugatePartitionMultiVariateMove(GetConjugatePartitionMultiVariateTreeProcess(),0.001,10),1,"conjugate sigma - process");
			*/

			/*
			if (!withtimeline)	{
			// if (Split())	{
				scheduler.Register(new ConjugatePartitionMultiVariateMove(GetConjugatePartitionMultiVariateTreeProcess(),10,10),1,"conjugate sigma - process");
				scheduler.Register(new ConjugatePartitionMultiVariateMove(GetConjugatePartitionMultiVariateTreeProcess(),1,10),1,"conjugate sigma - process");
				scheduler.Register(new ConjugatePartitionMultiVariateMove(GetConjugatePartitionMultiVariateTreeProcess(),0.1,10),1,"conjugate sigma - process");
				scheduler.Register(new ConjugatePartitionMultiVariateMove(GetConjugatePartitionMultiVariateTreeProcess(),0.01,10),1,"conjugate sigma - process");
				scheduler.Register(new ConjugatePartitionMultiVariateMove(GetConjugatePartitionMultiVariateTreeProcess(),0.001,10),1,"conjugate sigma - process");
			}
			*/
			// else	{
				/*
				scheduler.Register(new ConjugateMultiVariateExternalMove(GetConjugateInverseWishart(0),GetConjugateMultiVariateTreeProcess(),10,10),1,"external conjugate sigma - process");
				scheduler.Register(new ConjugateMultiVariateExternalMove(GetConjugateInverseWishart(0),GetConjugateMultiVariateTreeProcess(),1,10),1,"external conjugate sigma - process");
				scheduler.Register(new ConjugateMultiVariateExternalMove(GetConjugateInverseWishart(0),GetConjugateMultiVariateTreeProcess(),0.1,10),1,"external conjugate sigma - process");
				scheduler.Register(new ConjugateMultiVariateExternalMove(GetConjugateInverseWishart(0),GetConjugateMultiVariateTreeProcess(),0.01,10),1,"external conjugate sigma - process");
				scheduler.Register(new ConjugateMultiVariateExternalMove(GetConjugateInverseWishart(0),GetConjugateMultiVariateTreeProcess(),0.001,10),1,"external conjugate sigma - process");

				*/
				scheduler.Register(new ConjugateMultiVariateMove(GetConjugateInverseWishart(0),GetConjugateMultiVariateTreeProcess(),10,10),1,"conjugate sigma - process");
				scheduler.Register(new ConjugateMultiVariateMove(GetConjugateInverseWishart(0),GetConjugateMultiVariateTreeProcess(),1,10),1,"conjugate sigma - process");
				scheduler.Register(new ConjugateMultiVariateMove(GetConjugateInverseWishart(0),GetConjugateMultiVariateTreeProcess(),0.1,10),1,"conjugate sigma - process");
				scheduler.Register(new ConjugateMultiVariateMove(GetConjugateInverseWishart(0),GetConjugateMultiVariateTreeProcess(),0.01,10),1,"conjugate sigma - process");
			// }

			/*
			scheduler.Register(new ConjugateInverseWishartPriorMove(GetConjugateInverseWishart(),priorOnSigmaZero,10,10),1,"conjugate sigma - prior");
			scheduler.Register(new ConjugateInverseWishartPriorMove(GetConjugateInverseWishart(),priorOnSigmaZero,1,10),1,"conjugate sigma - prior");
			scheduler.Register(new ConjugateInverseWishartPriorMove(GetConjugateInverseWishart(),priorOnSigmaZero,0.1,10),1,"conjugate sigma - prior");
			scheduler.Register(new ConjugateInverseWishartPriorMove(GetConjugateInverseWishart(),priorOnSigmaZero,0.01,10),1,"conjugate sigma - prior");
			scheduler.Register(new ConjugateInverseWishartPriorMove(GetConjugateInverseWishart(),priorOnSigmaZero,0.001,10),1,"conjugate sigma - prior");
			*/

			if (withdrift)	{
				scheduler.Register(new SimpleMove(driftarray,10),10,"drift");
				scheduler.Register(new SimpleMove(driftarray,1),10,"drift");
				scheduler.Register(new SimpleMove(driftarray,0.1),10,"drift");
				scheduler.Register(new SimpleMove(driftarray,0.01),10,"drift");
			}

			if (withexpdrift)	{
				scheduler.Register(new SimpleMove(driftphiarray,10),10,"drift exp relax");
				scheduler.Register(new SimpleMove(driftphiarray,1),10,"drift exp relax");
				scheduler.Register(new SimpleMove(driftphiarray,0.1),10,"drift exp relax");
				scheduler.Register(new SimpleMove(driftphiarray,0.01),10,"drift exp relax");
			}
			if (withreldrift)	{
				scheduler.Register(new SimpleMove(driftphiarray2,10),10,"drift exp relax");
				scheduler.Register(new SimpleMove(driftphiarray2,1),10,"drift exp relax");
				scheduler.Register(new SimpleMove(driftphiarray2,0.1),10,"drift exp relax");
				scheduler.Register(new SimpleMove(driftphiarray2,0.01),10,"drift exp relax");

				scheduler.Register(new SimpleMove(driftarray2,10),10,"drift");
				scheduler.Register(new SimpleMove(driftarray2,1),10,"drift");
				scheduler.Register(new SimpleMove(driftarray2,0.1),10,"drift");
				scheduler.Register(new SimpleMove(driftarray2,0.01),10,"drift");
			}

			scheduler.Register(new SimpleMove(DiagArray,10),100,"theta");
			scheduler.Register(new SimpleMove(DiagArray,1),100,"theta");
			scheduler.Register(new SimpleMove(DiagArray,0.1),100,"theta");

			if (separatesyn)	{
				scheduler.Register(new SimpleMove(lognormalsyntree,10),100,"log normal syn");
				scheduler.Register(new SimpleMove(lognormalsyntree,1),100,"log normal syn");
				scheduler.Register(new SimpleMove(lognormalsyntree,0.1),100,"log normal syn");
				scheduler.Register(new SimpleMove(lognormalsyntree,0.01),100,"log normal syn");
				scheduler.Register(new SimpleMove(synsigma,10),100,"syn sigma");
				scheduler.Register(new SimpleMove(synsigma,1),100,"syn sigma");
				scheduler.Register(new SimpleMove(synsigma,0.1),100,"syn sigma");
				scheduler.Register(new SimpleMove(synsigma,0.01),100,"syn sigma");
			}

			if (separateomega)	{
				scheduler.Register(new SimpleMove(lognormalomegatree,10),100,"log normal omega");
				scheduler.Register(new SimpleMove(lognormalomegatree,1),100,"log normal omega");
				scheduler.Register(new SimpleMove(lognormalomegatree,0.1),100,"log normal omega");
				scheduler.Register(new SimpleMove(lognormalomegatree,0.01),100,"log normal omega");
				scheduler.Register(new SimpleMove(omegasigma,10),100,"omega sigma");
				scheduler.Register(new SimpleMove(omegasigma,1),100,"omega sigma");
				scheduler.Register(new SimpleMove(omegasigma,0.1),100,"omega sigma");
				scheduler.Register(new SimpleMove(omegasigma,0.01),100,"omega sigma");
			}

			if (withtimeline)	{
				scheduler.Register(new SimpleMove(timelinesigma,10),100,"timelinesigma");
				scheduler.Register(new SimpleMove(timelinesigma,1),100,"timelinesigma");
				scheduler.Register(new SimpleMove(timelinesigma,0.1),100,"timelinesigma");
				scheduler.Register(new SimpleMove(timelinesigma,0.01),100,"timelinesigma");

				scheduler.Register(new SimpleMove(timeline,10),100,"timeline");
				scheduler.Register(new SimpleMove(timeline,1),100,"timeline");
				scheduler.Register(new SimpleMove(timeline,0.1),100,"timeline");
				scheduler.Register(new SimpleMove(timeline,0.01),100,"timeline");
			}

			if (autoregressive)	{
				scheduler.Register(new SimpleMove(phi,10),100,"phi");
				scheduler.Register(new SimpleMove(phi,1),100,"phi");
				scheduler.Register(new SimpleMove(phi,0.1),100,"phi");
				scheduler.Register(new SimpleMove(phi,0.01),100,"phi");

				scheduler.Register(new SimpleMove(mean,1),100,"mean");
				scheduler.Register(new SimpleMove(mean,0.1),100,"mean");
				scheduler.Register(new SimpleMove(mean,0.01),100,"mean");
			}

			if (rawNstate == Naa)	{
				scheduler.Register(new PosRealVectorMove(exprelrate,1,1),10,"relrates");
				scheduler.Register(new PosRealVectorMove(exprelrate,0.1,2),10,"relrates");
				scheduler.Register(new PosRealVectorMove(exprelrate,0.03,4),10,"relrates");
				scheduler.Register(new SimpleMove(exprelrate,0.01),10,"relrates");
				scheduler.Register(new PosRealVectorTranslationMove(exprelrate,1),3,"relrates");
				scheduler.Register(new PosRealVectorTranslationMove(exprelrate,0.1),3,"relrates");
				scheduler.Register(new PosRealVectorTranslationMove(exprelrate,0.01),3,"relrates");

				if (omegaratiotree == 3)	{
					scheduler.Register(new SplitMultiVariateRelRateCompensatoryMove(process,exprelrate,aasplitMatrix,splitkrkctype,aasimilarityMatrix,1,tstvindex,omegatsindex,omegatv0index),10,"process relrate comp move");
					scheduler.Register(new SplitMultiVariateRelRateCompensatoryMove(process,exprelrate,aasplitMatrix,splitkrkctype,aasimilarityMatrix,0.1,tstvindex,omegatsindex,omegatv0index),10,"process relrate comp move");
					scheduler.Register(new SplitMultiVariateRelRateCompensatoryMove(process,exprelrate,aasplitMatrix,splitkrkctype,aasimilarityMatrix,0.01,tstvindex,omegatsindex,omegatv0index),10,"process relrate comp move");
				}
				else	{
					scheduler.Register(new MultiVariateRelRateCompensatoryMove(process,exprelrate,aasimilarityMatrix,1,omegaindex),10,"process relrate comp move");
					scheduler.Register(new MultiVariateRelRateCompensatoryMove(process,exprelrate,aasimilarityMatrix,0.1,omegaindex),10,"process relrate comp move");
					scheduler.Register(new MultiVariateRelRateCompensatoryMove(process,exprelrate,aasimilarityMatrix,0.01,omegaindex),10,"process relrate comp move");
				}
			}
			else	{
				if ((mutmodel == 0) || (mutmodel >= 4))	{
				// if (! mutmodel)	{
					/*
					scheduler.Register(new ProfileMove(relrate,0.1,1),10,"relrates");
					scheduler.Register(new ProfileMove(relrate,0.01,1),10,"relrates");
					scheduler.Register(new ProfileMove(relrate,0.001,1),10,"relrates");
					scheduler.Register(new ProfileMove(relrate,0.03,2),10,"relrates");
					scheduler.Register(new ProfileMove(relrate,0.003,2),10,"relrates");
					scheduler.Register(new SimpleMove(relrate,0.01),10,"relrates");
					*/

					scheduler.Register(new PosRealVectorMove(exprelrate,1,1),10,"relrates");
					scheduler.Register(new PosRealVectorMove(exprelrate,0.1,1),10,"relrates");
					scheduler.Register(new PosRealVectorMove(exprelrate,0.03,2),10,"relrates");
					scheduler.Register(new PosRealVectorMove(exprelrate,0.01,4),10,"relrates");
					scheduler.Register(new SimpleMove(exprelrate,0.01),10,"relrates");
					scheduler.Register(new PosRealVectorTranslationMove(exprelrate,1),3,"relrates");
					scheduler.Register(new PosRealVectorTranslationMove(exprelrate,0.1),3,"relrates");
					scheduler.Register(new PosRealVectorTranslationMove(exprelrate,0.01),3,"relrates");

				}
				else if (mutmodel == 2)	{
					scheduler.Register(new ProfileMove(tsrelrate,0.1,1),10,"relrates");
					scheduler.Register(new ProfileMove(tvrelrate,0.1,1),10,"relrates");
					scheduler.Register(new ProfileMove(tvrelrate,0.03,2),10,"relrates");
				}

				if ((mutmodel >= 4) || (omegaratiotree >= 4))	{
					scheduler.Register(new MultiVariateRelRateMutCompensatoryMove(process,exprelrate,1,tstvindex,tvgcindex),10,"process relrate comp move");
					scheduler.Register(new MultiVariateRelRateMutCompensatoryMove(process,exprelrate,0.1,tstvindex,tvgcindex),10,"process relrate comp move");
					scheduler.Register(new MultiVariateRelRateMutCompensatoryMove(process,exprelrate,0.1,tstvindex,tvgcindex),10,"process relrate comp move");
				}
			}

			if (gc == 3)	{
				scheduler.Register(new SimpleMove(rootgc1,1),10,"root gc");
				scheduler.Register(new SimpleMove(rootgc1,0.1),10,"root gc");
				scheduler.Register(new SimpleMove(rootgc1,0.01),10,"root gc");
				scheduler.Register(new SimpleMove(rootgc1,0.001),10,"root gc");
				scheduler.Register(new SimpleMove(rootgc2,1),10,"root gc");
				scheduler.Register(new SimpleMove(rootgc2,0.1),10,"root gc");
				scheduler.Register(new SimpleMove(rootgc2,0.01),10,"root gc");
				scheduler.Register(new SimpleMove(rootgc2,0.001),10,"root gc");
				scheduler.Register(new SimpleMove(rootgc3,1),10,"root gc");
				scheduler.Register(new SimpleMove(rootgc3,0.1),10,"root gc");
				scheduler.Register(new SimpleMove(rootgc3,0.01),10,"root gc");
				scheduler.Register(new SimpleMove(rootgc3,0.001),10,"root gc");
			}
			if (rootgc)	{
				scheduler.Register(new SimpleMove(rootgc,1),10,"root gc");
				scheduler.Register(new SimpleMove(rootgc,0.1),10,"root gc");
				scheduler.Register(new SimpleMove(rootgc,0.01),10,"root gc");
				scheduler.Register(new SimpleMove(rootgc,0.001),10,"root gc");
			}
			if (freestationary)	{
				scheduler.Register(new ProfileMove(freestationary,0.01,2),10,"stat4");
				scheduler.Register(new ProfileMove(freestationary,0.03,2),10,"stat4");
				scheduler.Register(new ProfileMove(freestationary,0.01,5),10,"stat10");
				scheduler.Register(new SimpleMove(freestationary,0.001),10,"stat");
			}
		}

		}
	}

};

