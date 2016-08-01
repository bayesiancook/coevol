
#include <sstream>
#include "BranchSitePath.h"

#include <cmath>

//-------------------------------------------------------------------------
//	* BranchSitePath
//-------------------------------------------------------------------------

BranchSitePath::BranchSitePath()	{
	init = last = new Plink;
	bkinit = bklast = new Plink;
	nsub = 0;
}

BranchSitePath::BranchSitePath(StateSpace* instatespace)	{
	statespace = instatespace;
	init = last = new Plink;
	nsub = 0;
}

BranchSitePath::~BranchSitePath()	{
	Reset(0);
	delete init;
}

void BranchSitePath::Reset(int state)	{
	Plink* link=last;
	while (link!=init)	{
		Plink* prev = link->Prev();
		delete link;
		link = prev;
	}
	nsub = 0;
	init->SetState(state);
	init->SetRelativeTime(0);
	last = init;
}

void BranchSitePath::BKReset(int state)	{
	Plink* link=bklast;
	while (link!=bkinit)	{
		Plink* prev = link->Prev();
		delete link;
		link = prev;
	}
	bknsub = 0;
	bkinit->SetState(state);
	bkinit->SetRelativeTime(0);
	bklast = bkinit;
}

void BranchSitePath::BackupPath()	{

	BKReset(init->GetState());
	Plink* link = init;
	while (link != last)	{
		BKAppend(link->Next()->GetState(),link->GetRelativeTime());
		link = link->Next();
	}
	bklast->SetRelativeTime(last->GetRelativeTime());
	bknsub = nsub;
}

void BranchSitePath::RestorePath()	{

	Reset(bkinit->GetState());
	Plink* link = bkinit;
	while (link != bklast)	{
		Append(link->Next()->GetState(),link->GetRelativeTime());
		link = link->Next();
	}
	last->SetRelativeTime(bklast->GetRelativeTime());
	nsub = bknsub;
}

string BranchSitePath::ToString(bool redundant)	{

	ostringstream s;
	if (redundant)	{
		s << GetCharFinalState() << ':';
	}
	Plink* link = last;
	while (link)	{
		if (link != last)	{
			s << ':' << GetState(link->Next()) << ':' << GetAbsoluteTime(link);
		}
		else	{
			s << GetAbsoluteTime(link);
		}
		/*
		if (link != last)	{
			s << ':';
		}
		s << GetState(link) << ':' << GetAbsoluteTime(link);
		*/
		link = link->prev;
	}
	if (redundant)	{
		s << ':' << GetCharInitState();
	}
	return s.str();
}

void BranchSitePath::AddCounts(int** paircounts, int* statecounts)	{

	if (last->GetRelativeTime())	{
		Plink* link = init;
		while (link != last)	{
			paircounts[link->GetState(),link->Next()->GetState()]++;
			statecounts[link->Next()->GetState()]++;
			link = link->Next();
		}
	}
	else	{
		statecounts[init->GetState()]++;
	}
}


void BranchSitePath::Prefix(BranchSitePath* p, BranchSitePath* root, double abstime)	{

	if ((init->GetState() != root->init->GetState()) || (p->init->GetState() != root->init->GetState()))	{
		cerr << "error in prefix: no matching states at root\n";
	}

	if (abstime > p->GetTotalTime())	{
		cerr << "error : time to be removed is too large\n";
		cerr << abstime << '\t' << p->GetTotalTime() << '\n';
		exit(1);
	}

	SetTimesRelativeToAbsolute();
	p->SetTimesRelativeToAbsolute();

	double t = abstime;
	while (t>0)	{
		double dt = p->init->GetRelativeTime();
		if (dt > t)	{
			// shift only t
			p->init->SetRelativeTime(p->init->GetRelativeTime() - t);
			init->SetRelativeTime(init->GetRelativeTime() + t);
			t = 0;
		}
		else	{
			if (p->init->IsLast())	{
				cerr << "error :  is last\n";
				exit(1);
			}
			// shift dt
			init->SetRelativeTime(init->GetRelativeTime() + dt);

			// splice first link of p
			Plink* tmp = p->init->next;
			if (! tmp)	{
				cerr << "error : next is 0\n";
				exit(1);
			}
			tmp->prev = 0;

			// put it at the base of this
			init->prev = p->init;
			p->init->next = init;
			init = p->init;
			p->init = tmp;

			// set time and state
			init->SetRelativeTime(0);
			init->SetState(p->init->GetState());
			t -= dt;
		}
	}
	if (p->init->GetState() != init->GetState())	{
		cerr << "error in Prefix\n";
		exit(1);
	}
	root->init->SetState(init->GetState());

	SetTimesAbsoluteToRelative();
	p->SetTimesAbsoluteToRelative();
}

void BranchSitePath::SetTimesRelativeToAbsolute()	{

	if (fabs(CheckTotalTime() - 1) > 1e-6)	{
		cerr << "error in BranchSitePath: relative time does not sum to 1 : " << CheckTotalTime()  << '\n';
		cerr << GetTotalTime() << '\n';
		exit(1);
	}
	Plink* link = init;
	double totaltime = GetTotalTime();
	while (link != last)	{
		link->SetRelativeTime(link->GetRelativeTime() * totaltime);
		link = link->next;
	}
	link->SetRelativeTime(link->GetRelativeTime() * totaltime);
}


void BranchSitePath::SetTimesAbsoluteToRelative()	{

	/*
	if (fabs(CheckTotalTime() - GetTotalTime()) > 1e-6)	{
		cerr << "error in BranchSitePath: relative time does not sum to what it should : " << CheckTotalTime() << '\t' << GetTotalTime() << '\n';
		exit(1);
	}
	*/
	Plink* link = init;
	double totaltime = CheckTotalTime();
	while (link != last)	{
		link->SetRelativeTime(link->GetRelativeTime() / totaltime);
		link = link->next;
	}
	link->SetRelativeTime(link->GetRelativeTime() / totaltime);
}

double BranchSitePath::CheckTotalTime()	{

	double tot = 0;
	Plink* link = init;
	while (link != last)	{
		if (! link)	{
			cerr << "error : null link\n";
			exit(1);
		}
		tot += link->GetRelativeTime();
		link = link->next;
	}
	if (! link)	{
		cerr << "error : null link\n";
		exit(1);
	}
	tot += link->GetRelativeTime();
	return tot;
}

