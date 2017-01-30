#include "BrownianModel.h"

BrownianModel::~BrownianModel()
{
	//dtor
}

double BrownianModel::GetLogPrior()	{
	double total = 0;
	total += chronogram->GetLogProb();
	total += sigma->GetLogProb();
	total += brownianProcess->GetLogProb();
	total += relrate->GetLogProb();
	total += stationary->GetLogProb();
	return total;
}

double BrownianModel::GetLogLikelihood()	{
	double ret = phyloprocess->GetLogProb();
	return ret;
}

void BrownianModel::drawSample()	{
	chronogram->Sample();
	sigma->Sample();
	brownianProcess->Sample();
	stationary->Sample();
	relrate->Sample();
	phyloprocess->Sample();
	cerr << "ok\n";
}

void BrownianModel::MakeScheduler()	{

	if (pathconjugate)	{
		if(commut)
			scheduler.Register(new DSemiConjugateMappingMove(phyloprocess,commutPathconjtree),1,"mapping + sufficient stat");
		else if(mapSegment)
			scheduler.Register(new DSemiConjugateMappingMove(phyloprocess,segmentedPathconjtree),1,"mapping + sufficient stat");
		else
			scheduler.Register(new DSemiConjugateMappingMove(phyloprocess,pathconjtree),1,"mapping + sufficient stat");
	}
	else	{
		scheduler.Register(new SimpleMove(phyloprocess,1),1,"mapping");
	}

	scheduler.OpenLoop(100);

	if(!fixtimes) {
		if(segm == SEGM_ABSOLUTE) {
			scheduler.Register(new BrownianHorizontalMove(chronogram, brownianProcess, 1), 10, "horizontal move");
			scheduler.Register(new BrownianHorizontalMove(chronogram, brownianProcess, 0.1), 10, "horizontal move");
			scheduler.Register(new BrownianHorizontalMove(chronogram, brownianProcess, 0.01), 10, "horizontal move");
		}
		else {
			scheduler.Register(new SimpleMove(chronogram, 1), 10, "chronogram");
			scheduler.Register(new SimpleMove(chronogram, 0.1), 10, "chronogram");
			scheduler.Register(new SimpleMove(chronogram, 0.01), 10, "chronogram");
		}
	}

 
	if (conjsigma == 2)	{

		scheduler.Register(new ConjugateBrownianExternalMove(conjugatesigma, conjugateBrownianProcess, 0.1,10) ,1,"External Conjugate Brownian sigma");
		scheduler.Register(new ConjugateBrownianExternalMove(conjugatesigma, conjugateBrownianProcess, 0.01,10) ,1,"External Conjugate Brownian sigma");
		scheduler.Register(new ConjugateBrownianExternalAllBranchMove(conjugatesigma, conjugateBrownianProcess, 0.1,10) ,1,"External Conjugate Brownian all branch sigma");
		scheduler.Register(new ConjugateBrownianExternalAllBranchMove(conjugatesigma, conjugateBrownianProcess, 0.01,10) ,1,"External Conjugate Brownian all branch sigma");
		scheduler.Register(new ConjugateBrownianExternalAllBranchMove(conjugatesigma, conjugateBrownianProcess, 0.001,10) ,1,"External Conjugate Brownian all branch sigma");
		/*
		scheduler.Register(new FullConjugateBrownianSigmaMove(conjugatesigma, conjugateBrownianProcess, 0.1,10) ,1,"Conjugate PureBrownian sigma");
		scheduler.Register(new FullConjugateBrownianSigmaMove(conjugatesigma, conjugateBrownianProcess, 0.01,10) ,1,"Conjugate PureBrownian sigma");
		*/
	}
	else	{

		scheduler.Register(new SimpleMove(brownianProcess->GetPureBrownianProcess(),0.1),10,"PureBrownian");
		scheduler.Register(new SimpleMove(brownianProcess->GetPureBrownianProcess(),0.01),10,"PureBrownian");
		scheduler.Register(new SimpleMove(brownianProcess->GetInstantProcess(),1),10,"InstantValues");
		scheduler.Register(new SimpleMove(brownianProcess->GetInstantProcess(),0.1),10,"InstantValues");
	 
		scheduler.Register(new BrownianSigmaMove(sigma, brownianProcess, 0.1), 10, "linear sigma brownian");
		scheduler.Register(new BrownianSigmaMove(sigma, brownianProcess, 0.01), 10, "linear sigma brownian");
		scheduler.Register(new BrownianSigmaMove(sigma, brownianProcess, 0.001), 10, "linear sigma brownian");

		if(conjsigma) {
			scheduler.Register(new ConjugateBrownianSigmaMove(dynamic_cast<ConjugateInverseWishart*>(sigma), brownianProcess), 1, "sigma(Gibbs)");
		}
		else {
			scheduler.Register(new SimpleMove(sigma,1),100,"sigma");
			scheduler.Register(new SimpleMove(sigma,0.1),100,"sigma");
			scheduler.Register(new SimpleMove(sigma,0.01),100,"sigma");
		}
	}

	scheduler.Register(new ProfileMove(relrate,0.1,1),10,"relrates");
	scheduler.Register(new ProfileMove(relrate,0.03,2),10,"relrates");
	scheduler.Register(new SimpleMove(relrate,0.01),10,"relrates");

	if(gc) {
		scheduler.Register(new SimpleMove(rootGCrate,0.1),10,"rootstat");
		scheduler.Register(new SimpleMove(rootGCrate,0.01),10,"rootstat");
		scheduler.Register(new SimpleMove(rootGCrate,0.001),10,"rootstat");
	}

	scheduler.Register(new ProfileMove(stationary,0.01,2),10,"stat4");
	scheduler.Register(new ProfileMove(stationary,0.03,2),10,"stat4");
	scheduler.Register(new ProfileMove(stationary,0.01,5),10,"stat10");
	scheduler.Register(new SimpleMove(stationary,0.001),10,"stat");

	scheduler.CloseLoop();

}

double BrownianModel::GetMeanRho()	{
	return brownianProcess->GetMeanRate(0);
}

double BrownianModel::GetIntegralRho()	{
	return brownianProcess->GetIntegralRate(0);
}


double BrownianModel::GetVarRho()	{
	return brownianProcess->GetVarRate(0);
}


double BrownianModel::GetMeanOmega()	{
	return brownianProcess->GetMeanRate(1);
}

double BrownianModel::GetVarOmega()	{
	return brownianProcess->GetVarRate(1);
}

double BrownianModel::GetMeanGC()	{
	return brownianProcess->GetMeanGC(1);
}

double BrownianModel::GetVarGC()	{
	return brownianProcess->GetVarGC(1);
}

double BrownianModel::GetRootGC()	{
	return rootGCrate->val();
}


void BrownianModel::TraceHeader(ostream& os)	{
	os << "#logprior\tlnL\ttime\tdS";
	for (int k=0; k<Ncont+L; k++)	{
		for (int l=k+1; l<Ncont+L; l++)	{
			os << '\t' << "sigma_" << k << '_' << l;
		}
	}
	for (int k=0; k<Ncont+L; k++)	{
		os << '\t' << "sigma_" << k << '_' << k;
	}
	os << "\tmeannsub\tminnsub\tmaxnsub\tnonesub\tstatent\trrent\tmeanrho\tvarrho";

	if(gc)
		os <<"\tmeangc\tvargc\trootgamma";
	if(omega)
		os <<"\tmeanomega\tvaromega";

	if (gc || omega)	{
		os << "\tmeanz\tmaxz";
	}

	os << '\n';

}

void BrownianModel::Trace(ostream& os)	{
	os << GetLogPrior() << '\t' << GetLogLikelihood() << '\t' << chronogram->GetTotalTime() << '\t' << GetIntegralRho();
	for (int k=0; k<Ncont+L; k++)	{
		for (int l=k+1; l<Ncont+L; l++)	{
			os << '\t' << (*sigma)[k][l];
		}
	}
	for (int k=0; k<Ncont+L; k++)	{
		os << '\t' << (*sigma)[k][k];
	}
	os << '\t' << brownianProcess->GetMeanNSubBranch() << '\t' << brownianProcess->GetMinNSubBranch() << '\t' << brownianProcess->GetMaxNSubBranch() << '\t' << brownianProcess->GetNOneSubBranch() << '\t' << stationary->val().GetEntropy() << '\t' << relrate->val().GetEntropy() << '\t' << GetMeanRho() << '\t' << GetVarRho() ;
	if(gc)
		os << '\t' << GetMeanGC() << '\t' << GetVarGC() << '\t' << GetRootGC();
	if(omega)
		 os << '\t' << GetMeanOmega() << '\t' << GetVarOmega() ;

	if (gc || omega)	{
		os << '\t' << SubMatrix::meanz / SubMatrix::nz << '\t' << SubMatrix::maxz;
	}
	SubMatrix::nz = 0;
	SubMatrix::meanz = 0;
	SubMatrix::maxz = 0;

	os << '\n';
	os.flush();

}

void BrownianModel::Details(ostream& os)	{

	if(mapSegment) {
		int* x = evolutionSegmentedProcess->GetNbOrders();
		double tot = 0;
		for(int i=0; i< 11; i++)
			tot += x[i];
		for(int i=0; i<11; i++) {
			os << "Order " << (i+1) << " : " << x[i] << " (" << 100*x[i]/tot << "%)" << endl;
		}
		os << "Total : " << tot << endl << endl;
	}
}

void BrownianModel::ToStream(ostream& os)	{
	os << *chronogram << '\n';
	os << *DiagArray << '\n';
	os << *sigma << '\n';
	if (gc)	{
		os << *rootGCrate << '\n';
	}
	os << *brownianProcess->GetPureBrownianProcess() << '\n';
	os << *brownianProcess->GetInstantProcess() << '\n';
	os << *relrate << '\n';
	os << *stationary << '\n';
}

void BrownianModel::FromStream(istream& is)	{
	is >> *chronogram;
	is >> *DiagArray;
	is >> *sigma;
	if (gc)	{
		is >> *rootGCrate;
	}
	is >> *brownianProcess->GetPureBrownianProcess();
	is >> *brownianProcess->GetInstantProcess();
	is >> *relrate;
	is >> *stationary;
}

void BrownianModel::ToStreamShort(ostream& os)	{
	os << *chronogram << '\n';
	os << *DiagArray << '\n';
	os << *sigma << '\n';
	if (gc)	{
		os << *rootGCrate << '\n';
	}
	os << *relrate << '\n';
	os << *stationary << '\n';
}

void BrownianModel::FromStreamShort(istream& is)	{
	is >> *chronogram;
	is >> *DiagArray;
	is >> *sigma;
	if (gc)	{
		is >> *rootGCrate;
	}
	is >> *relrate;
	is >> *stationary;
}


void BrownianModel::PrintEntries(ostream& os)	{

	os << "dS\n";

	for (int k=0; k<Ncont; k++)	{
		os << contdata->GetCharacterName(k) << '\n';


	}
}

void BrownianModel::Load(string simufile)	{

	ifstream is(simufile.c_str());
	FromStreamShort(is);
	Update();
}

void BrownianModel::PostPredSimu(string name, bool resamplecov, bool resampleprocess, int nsegment)	{

	int cont = 1;
	while (cont)	{
		if (resamplecov)	{
			sigma->Sample();
		}

		if (resampleprocess)	{
			MultiVariateTreeProcess* process = brownianProcess->GetInstantProcess();
			for (int i=0; i<sigma->GetDim(); i++)	{
				process->Unclamp(i);
			}
			process->GetNodeVal(process->GetRoot()->GetNode())->Clamp();

			/*
			if (nsegment)	{
				BrownianBridge::setNTreeSegments(nSegments);
			}
			*/
			brownianProcess->Sample();

			cont = 0;
			for (int i=0; i<L; i++)	{
				if ((process->GetMin(i) < -5) || (process->GetMax(i) > 5))	{
					cont = 1;
				}
			}
			cerr << cont << '\t';
		}
	}
	cerr << '\n';

	cerr << "update\n";

	Update();

	cerr << "grep alignment\n";
	SequenceAlignment* datacopy = 0;
	if (codondata)	{
		datacopy = new CodonSequenceAlignment(codondata);
	}
	else	{
		datacopy = new SequenceAlignment(data);
	}

	cerr << "phyloprocess\n";
	phyloprocess->PostPredSample();
	phyloprocess->GetLeafData(datacopy);

	cerr << "print out\n";
	ostringstream s;
	ofstream os((name + ".ali").c_str());
	datacopy->ToStream(os);

	ofstream cos((name + ".cov").c_str());
	cos << *sigma << '\n';
	/*
	cerr << "MAX dS : " << GetMaxdS() << '\n';
	// cerr << "MAX dN : " << GetMaxdN() << '\n';
	*/

	if (contdata)	{
		cerr << "get cont data\n";
		ContinuousData* contdatacopy = new ContinuousData(contdata);
		cerr << L << '\n';
		brownianProcess->GetInstantProcess()->GetLeafData(contdatacopy,L);
		ofstream cos((name + ".cont").c_str());
		cerr << "to stream\n";
		contdatacopy->ToStream(cos);
		cerr << "delete cont data copy\n";
		delete contdatacopy;
		cerr << "ok\n";
	}

	// times 
	ofstream tos((name + ".dates").c_str());
	GetChronogram()->PrintAges(tos);

	// log rates amd lengths
	ofstream los((name + ".lengths").c_str());
	brownianProcess->PrintLengths(los);

	// node vals
	ofstream vos((name + ".nodevals").c_str());
	brownianProcess->PrintNodeVals(vos);

	ofstream pros((name + ".param").c_str());
	ToStream(pros);
	cerr << "delete\n";
	delete datacopy;
	cerr << "ok\n";
}

/*
void BrownianModel::SpecialUpdate()	{

}
*/
