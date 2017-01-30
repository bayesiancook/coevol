
#include "Move.h"
#include "Chrono.h"
#include "Random.h"
#include "ProbModel.h"

#include <iostream>



void MCScheduler::Reset()	{

	size = update.size();
	time = vector<double>(size);
	success = vector<double>(size);
	ncall = vector<int>(size);

	totalweight = 0;
	for (int i=0; i<size; i++)	{
		totalweight += weight[i];
		time[i] = 0;
		success[i] = 0;
	}
	totaltime = 0;
	ncycle = 0;
	closed = true;
 
}

void MCScheduler::Register(MCUpdate* inupdate, int inweight, string inname)	{

	if (closed)	{
		cerr << "error in MCScheduler::Register: registration of new updates is closed\n";
		throw;
	}
	update.push_back(inupdate);
	weight.push_back(inweight);
	name.push_back(inname);

		ostringstream oss;
		oss << (update.size()-1) << ",";
		command += oss.str();

}

void MCScheduler::Cycle(double tuning_modulator, int nrep, bool verbose, bool check)	{

	if (! closed) Reset();

	for (int rep = 0; rep<nrep; rep++)	{
			unsigned int n=0;
			vector<int> moves = ReadCommand(n);

			if(n!=command.size()) {
				cerr << "Error in scheduler command line : " << endl << command << endl;
				exit(1);
			}

			for (unsigned int i=0; i<moves.size(); i++)	{
					Move(tuning_modulator,moves[i],verbose,check,weight[moves[i]]);
			}
		/*for (int i=0; i<size; i++)	{
			Move(tuning_modulator,i,verbose,check,weight[i]);
		}*/
	}
	ncycle+=nrep;


}

vector<int> MCScheduler::ReadCommand(unsigned int &n)	{
	vector<int> v;
	while(n<command.size() && command[n]!=')') {
		if(command[n]=='(') {
			unsigned int n2 = command.find(":", n+1);
			/*
			if(n2 == string::npos) {
				cerr << "Error in scheduler command line : " << endl << command << endl;
				exit(1);
			}
			*/
			int x= atoi(command.substr(n+1, n2).c_str());
			n = n2+1;
			vector<int> vv = ReadCommand(n);
			if(command[n] != ')') {
				cerr << "Error in scheduler command line : " << endl << command << endl;
				exit(1);
			}

			for(int i = 0; i<x; i++)
				for(unsigned int j=0; j<vv.size(); j++)
					v.push_back(vv[j]);
			n++;
		}
		else {
			unsigned int n2 = command.find(",", n+1);
			/*
			if(n2 == string::npos) {
				cerr << "Error in scheduler command line : " << endl << command << endl;
				exit(1);
			}
			*/
			int x= atoi(command.substr(n, n2).c_str());
			v.push_back(x);
			n = n2+1;
		}
	}
	return v;
}


void MCScheduler::Move(double tuning_modulator, int i, bool verbose, bool check, int nrep)	{

	Chrono chrono;
	chrono.Reset();
	chrono.Start();
	double meansuccess = 0;
	if (verbose)	{
		cerr << i << ' ' << name[i] << '\n';
	}
	for (int k=0; k<nrep; k++)	{
		meansuccess += update[i]->Move(tuning_modulator);
	}
	chrono.Stop();
	totaltime += chrono.GetTime() / 1000;
	time[i] += chrono.GetTime() / 1000;
	ncall[i] += nrep;
	success[i] += meansuccess;

	if (verbose)	{
		cerr << "success "  << meansuccess << '\n';
	}
	if (check)	{
		try	{
			if (verbose)	{
				cerr << "check\n";
			}
			/*
			if (! model->CheckUpdateFlags())	{
				cerr << "error : flags are not all up\n";
				// exit(1);
			}
			*/

			double tmp = model->Update(true);

			if (fabs(tmp) > MCMC::MAXDIFF)	{
				cerr << "total is wrong\n";
				throw CheckSumException(tmp);
			}

			if (verbose)	{
				cerr << "check ok\n";
			}
		}
		catch(CheckSumException e)	{
			cerr << "NON ZERO CHECKSUM. After " << name[i] << " : " << e.GetCheckSum() << '\n';

			throw e;
		}
	}
}

void MCScheduler::RandomCycle(double tuning_modulator, int nrep, bool verbose, bool check)	{

	if (! closed) Reset();

	int N = (int) (nrep * totalweight);
	for (int rep = 0; rep<N; rep++)	{

		double q = totalweight * Random::Uniform();
		double total = weight[0];
		int choose = 0;
		while ((choose < size) && (q > total))	{
			choose++;
			total += weight[choose];
		}
		if (choose == size)	{
			cerr << "error in MCScheduler::RandomCycle: overflow\n";
			exit(1);
		}
		Move(tuning_modulator, choose,verbose,check,1);
	}
	ncycle+=nrep;
}

void MCScheduler::ToStream(ostream& os, ostream& osdetail)	{

	os << '\n';
	os << "total number of cycles : " << ncycle << '\n';
	os << "total time (s)		 : " << totaltime << '\n';
	os << "time per cycle (s)	 : " << totaltime / ncycle << '\n';
	os << '\n';
	os << "#\t\%time\tsuccess\tname\n";
	os << '\n';
	for (int i=0; i<size; i++)	{
		// os << i+1 << '\t' << (int) (time[i]/ totaltime * 100)  << '\t' << ((int) (100 * success[i]/ncall[i])) << '\t' << name[i]  << '\n';
		os << i+1 << '\t' << ((double) (int) (time[i]/ totaltime * 1000)) / 10  << '\t' << ((int) (100 * success[i]/ncall[i])) << '\t' << name[i]  << '\n';
	}
	os << '\n';

		for(int i=0; i<size; i++)
			update[i]->ToStream(osdetail);

}

RealVectorComponentwiseCompensatoryMove::RealVectorComponentwiseCompensatoryMove(Rvar<RealVector>* ina1, Rvar<RealVector>* ina2, double intuning) : a1(ina1), a2(ina2), tuning(intuning) {

	cerr << "error in comp move\n";
	cerr << "register should be the other way around\n";
	exit(1);

	Register(a1);
	Register(a2);
}

double RealVectorComponentwiseCompensatoryMove::Move(double tuning_modulator)	{

	Corrupt(true);

	int dim = a1->GetDim();
	for (int i=0; i<dim; i++)	{
		double h = tuning_modulator * tuning * (Random::Uniform() - 0.5);
		(*a1)[i] += h;
		(*a2)[i] -= h;
	}

	double logratio = Update();
	bool accepted = (log(Random::Uniform()) < logratio);

	if (! accepted)	{
		Corrupt(false);
		Restore();
	}

	return (double) accepted;
}

/*
OneToManyRealVectorComponentwiseCompensatoryMove::OneToManyRealVectorComponentwiseCompensatoryMove(Rvar<RealVector>* ina1, IIDNormal** ina2, int inK, double intuning) : a1(ina1), a2(ina2), K(inK), tuning(intuning) {

	a1->Register(this);
	for (int i=0; i<K; i++)	{
		a2[i]->Register(this);
	}
}

double OneToManyRealVectorComponentwiseCompensatoryMove::Move(double tuning_modulator)	{

	Corrupt(true);

	int dim = a1->GetDim();
	for (int i=0; i<dim; i++)	{
		double h = tuning_modulator * tuning * (Random::Uniform() - 0.5);
		(*a1)[i] += h;
		for (int k=0; k<K; k++)	{
			(*(a2[k]))[i] -= h;
		}
	}

	double logratio = Update();
	bool accepted = (log(Random::Uniform()) < logratio);

	if (! accepted)	{
		Corrupt(false);
		Restore();
	}

	return (double) accepted;
}
*/

RealVectorMove::RealVectorMove(Rvar<RealVector>* invar, double intuning, int inm) : var(invar), tuning(intuning), m(inm) {}

double RealVectorMove::Move(double tuning_modulator)	{
	if (! var->isClamped())	{
		var->Corrupt(true);
		double logHastings = var->RealVector::ProposeMove(tuning * tuning_modulator,m);
		double deltaLogProb = var->Update();
		double logRatio = deltaLogProb + logHastings;
		bool accepted = (log(Random::Uniform()) < logRatio);
		if (! accepted)	{
			var->Corrupt(false);
			var->Restore();
		}
		return (double) accepted;
	}
	return 1;
}

RealVectorTranslationMove::RealVectorTranslationMove(Rvar<RealVector>* invar, double intuning) : var(invar), tuning(intuning) {}

double RealVectorTranslationMove::Move(double tuning_modulator)	{
	if (! var->isClamped())	{
		var->Corrupt(true);
		double logHastings = var->RealVector::ScalarAddition(tuning * tuning_modulator);
		double deltaLogProb = var->Update();
		double logRatio = deltaLogProb + logHastings;
		bool accepted = (log(Random::Uniform()) < logRatio);
		if (! accepted)	{
			var->Corrupt(false);
			var->Restore();
		}
		return (double) accepted;
	}
	return 1;
}


ProfileMove::ProfileMove(Rvar<Profile>* invar, double intuning, int inn) : var(invar), tuning(intuning), n(inn) {}

double ProfileMove::Move(double tuning_modulator)	{

	if (! var->isClamped())	{
		var->Corrupt(true);
		double logHastings = var->Profile::ProposeMove(tuning * tuning_modulator, n);
		double deltaLogProb = var->Update();
		double logRatio = deltaLogProb + logHastings;
		bool accepted = (log(Random::Uniform()) < logRatio);
		if (! accepted)	{
			var->Corrupt(false);
			var->Restore();
		}
		return (double) accepted;
	}
	return 1;
}

MultiplicativeCompensatoryMove::MultiplicativeCompensatoryMove(Multiplicative* inm1, Multiplicative* inm2, double intuning) : m1(inm1), m2(inm2), tuning(intuning) {

	if (m1)	{
		m1->Register(this);
	}
	if (m2)	{
		m2->Register(this);
	}
}

double MultiplicativeCompensatoryMove::Move(double tuning_modulator)	{

	// if ((! m1->isClamped()) && (! m2->isClamped()))	{
		Corrupt(true);
		double h = tuning_modulator * tuning * (Random::Uniform() - 0.5);
		double e = exp(h);
		double loghastings = 0;
		if (m1)	{
			loghastings += h * m1->ScalarMultiplication(e);
		}
		if (m2)	{
			loghastings -= h * m2->ScalarMultiplication(1.0/e);
		}

		double logratio = Update();
		logratio += loghastings;
		bool accepted = (log(Random::Uniform()) < logratio);
		if (! accepted)	{
			Corrupt(false);
			Restore();
		}
		return (double) accepted;
	// }
	// return 1;
}

JointSimpleMove::JointSimpleMove(Rnode* ina1, Rnode* ina2, double intuning) : a1(ina1), a2(ina2), tuning(intuning) {

	a1->Register(this);
	a2->Register(this);
}

double JointSimpleMove::Move(double tuning_modulator)	{

	if ((! a1->isClamped()) && (! a2->isClamped()))	{

		Corrupt(true);

		double logHastings = 0;
		logHastings += a1->ProposeMove(tuning*tuning_modulator);
		logHastings += a2->ProposeMove(tuning*tuning_modulator);

		double logratio = Update() + logHastings;
		bool accepted = (log(Random::Uniform()) < logratio);

		if (! accepted)	{
			Corrupt(false);
			Restore();
		}

		return (double) accepted;
	}
	return 1;
}

AdditiveCompensatoryMove::AdditiveCompensatoryMove(Additive* ina1, Additive* ina2, double intuning) : a1(ina1), a2(ina2), tuning(intuning) {

	a1->Register(this);
	a2->Register(this);
}

double AdditiveCompensatoryMove::Move(double tuning_modulator)	{

	// if ((! a1->isClamped()) && (! a2->isClamped()))	{

		Corrupt(true);

		double h = tuning_modulator * tuning * (Random::Uniform() - 0.5);
		a1->ScalarAddition(h);
		a2->ScalarAddition(-h);

		double logratio = Update();
		bool accepted = (log(Random::Uniform()) < logratio);

		if (! accepted)	{
			Corrupt(false);
			Restore();
		}

		return (double) accepted;
	// }
	// return 1;
}

AdditiveAntiCompensatoryMove::AdditiveAntiCompensatoryMove(Additive* ina1, Additive* ina2, double intuning) : a1(ina1), a2(ina2), tuning(intuning) {

	a1->Register(this);
	a2->Register(this);
}

double AdditiveAntiCompensatoryMove::Move(double tuning_modulator)	{

	// if ((! a1->isClamped()) && (! a2->isClamped()))	{

		Corrupt(true);

		double h = tuning_modulator * tuning * (Random::Uniform() - 0.5);

		a1->ScalarAddition(h);
		a2->ScalarAddition(h);

		double logratio = Update();
		bool accepted = (log(Random::Uniform()) < logratio);

		if (! accepted)	{
			Corrupt(false);
			Restore();
		}
		else	{
		}

		return (double) accepted;
	// }
	// return 1;
}

