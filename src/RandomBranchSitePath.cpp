#include <vector>
#include <algorithm>
#include <utility>

#include "RandomBranchSitePath.h"
#include "PhyloProcess.h"

//-------------------------------------------------------------------------
//	* RandomBranchSitePath
//-------------------------------------------------------------------------

StateSpace* RandomBranchSitePath::GetStateSpace() {return GetPhyloProcess()->GetStateSpace();}

bool RandomBranchSitePath::SampleBranchMapping()	{
	return myprocess->SampleBranchMapping();
}

void RandomBranchSitePath::SetUp(RandomBranchSitePath* inup)	{
	// pathup = inup;
	// Register(up);
}

bool RandomBranchSitePath::ResampleAcceptReject(int maxtrial)	{

	int ntrial = 0;
	do	{
		ntrial++;
		Reset(stateup);
		double t = 0;
		double totaltime = GetTime();
		int state = stateup;

		if (state != statedown)	{
			double u = DrawWaitingTimeGivenAtLeastOne(state,totaltime);
			t += u;
			int newstate = DrawOneStep(state);
			Append(newstate,u/totaltime);
			state = newstate;
		}
		while (t < totaltime)	{
			double u = DrawWaitingTime(state);
			t += u;
			if (t < totaltime)	{
				int newstate = DrawOneStep(state);
				Append(newstate,u/totaltime);
				state = newstate;
			}
			else	{
				t -= u;
				u = totaltime - t;
				last->SetRelativeTime(u/totaltime);
				t = totaltime;
			}
		}
	} while ((ntrial < maxtrial) && (last->GetState() != statedown));
	return (last->GetState() == statedown);
}

void RandomBranchSitePath::ResampleUniformized()	{

	// CheckUniformizedSubstitutionNumber(stateup, statedown);
	int m = DrawUniformizedSubstitutionNumber(stateup, statedown);

	vector<double> y(m+1);
	for (int r=0; r<m; r++)	{
		y[r] = Random::Uniform();
	}
	y[m] = 1;
	sort(y.begin(),y.end());

	int state = stateup;

	Reset(stateup);

	double t = y[0];
	for (int r=0; r<m; r++)	{
		int k = (r== m-1) ? statedown : DrawUniformizedTransition(state,statedown,m-r-1);
		if (k != state)	{
			Append(k,t);
			t = 0;
		}
		state = k;
		t += y[r+1] - y[r];
	}
	last->SetRelativeTime(t);
}

/*
void RandomBranchSitePath::ResampleUniformized()	{

	uninsub = DrawUniformizedSubstitutionNumber(stateup, statedown);

	unitime.resize(uninsub+1);
	unistate.resize(uninsub+1);

	for (int r=0; r<uninsub; r++)	{
		unitime[r] = Random::Uniform();
	}
	unitime[uninsub] = 1;
	sort(unitime.begin(),unitime.end());

	int state = stateup;
	// unistate[0] = stateup;

	Reset(stateup);

	double t = unitime[0];
	for (int r=0; r<uninsub; r++)	{
		int k = (r== uninsub-1) ? statedown : DrawUniformizedTransition(state,statedown,uninsub-r-1);
		unistate[r] = k;
		if (k != state)	{
			Append(k,t);
			t = 0;
		}
		state = k;
		t += unitime[r+1] - unitime[r];
	}
	last->SetRelativeTime(t);
}

double RandomBranchSitePath::RecordResampleUniformizedLogProb()	{

	double total = 0;

	total += RecordDrawUniformizedSubstitutionNumberLogProb(stateup, statedown, uninsub);

	int state = stateup;
	for (int r=0; r<uninsub; r++)	{
		int k = unistate[r];
		total += RecordDrawUniformizedTransitionLogProb(state,statedown,uninsub-r-1,k);
		state = k;
	}
	return total;
}
*/

double RandomBranchSitePath::logProb()	{

	if (SampleBranchMapping())	{
		return PathLogProb();
	}
	return NoPathLogProb();
}

double RandomBranchSitePath::NoPathLogProb()	{
	if (isRoot())	{
		return log((GetStationary())[stateup]);
	}
	return log((*GetTransitionMatrix())(stateup,statedown));
}

double RandomBranchSitePath::PathLogProb()	{
	if (! isActivated())	{
		return 0;
	}
	if (isRoot())	{
		return StationaryLogProb(init->GetState());
	}
	double total = 0;
	Plink* link = init;
	while (link)	{
		if (link != last)	{
			total += ReducedWaitingTimeLogProb(link->GetState(), GetAbsoluteTime(link));
			total += ReducedOneStepLogProb(link->GetState(),link->next->GetState());
			/*
			if (isinf(ReducedOneStepLogProb(link->GetState(),link->next->GetState())))	{
				cerr << "forbidden transition\n";
				cerr << GetStateSpace()->GetState(link->GetState()) << '\t' << GetStateSpace()->GetState(link->next->GetState()) << '\n';
				exit(1);
			}
			*/
		}
		else	{
			total += ReducedWaitingTimeLogProb(link->GetState(), GetAbsoluteTime(link));
		}
		link = link->next;
	}
	if (GetTime() > 0)	{
		total += GetNsub() * log(GetRate() * GetTime());
	}
	else	{
		cerr << "error in RandomBranchSitePath::logProb : null branch length\n";
		cerr << length << '\n';
		cerr << length->val() << '\n';
		exit(1);
	}
	return total;
}


int RandomBranchSitePath::GetMaxTrial()	{
	return myprocess->GetMaxTrial();
}
