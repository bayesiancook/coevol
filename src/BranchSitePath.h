
#ifndef SITEPATH_H
#define SITEPATH_H

#include <string>
#include "Tree.h"
#include "StateSpace.h"

class Plink	{

	friend class BranchSitePath;
	friend class RandomBranchSitePath;

	public:

				Plink();
				Plink(int instate, double inrel_time);
				~Plink();
	Plink* 			Prev();
	Plink* 			Next();

	bool 			IsFirst();
	bool 			IsLast();

	void 			Splice();
	void 			Insert(Plink* link);

	void 			SetState(int instate);
	int 			GetState();

	void 			SetRelativeTime(double inrel_time);
	double 			GetRelativeTime();

	private:

	Plink* next;
	Plink* prev;

	int state;
	double rel_time;

};

class BranchSitePath {

	friend class RandomBranchSitePath;

	public:

				BranchSitePath();
				BranchSitePath(StateSpace* instatespace);
	virtual			~BranchSitePath();

	virtual double 		GetTotalTime() = 0;
	virtual void		SetTotalTime(double intime) = 0;

	Plink* 			Init();
	Plink*  		Last();
	StateSpace* 		GetStateSpace();
	void			SetStateSpace(StateSpace* instatespace) {statespace = instatespace;}
	int 			GetNsub();
	string			GetState(Plink* link);
	string 			GetCharInitState();
	string 			GetCharFinalState();
	int 			GetInitState();
	int 			GetFinalState();
	double 			GetAbsoluteTime(Plink* link);
	double 			GetRelativeTime(Plink* link) {return link->GetRelativeTime();}

	void			AddCounts(int** paircount, int* statecount);

	void			SetTimesRelativeToAbsolute();
	void			SetTimesAbsoluteToRelative();
	double 			CheckTotalTime();

	void 			Reset(int state);
	void 			Append(int instate, double reltimelength);

	void 			BKReset(int state);
	void 			BKAppend(int instate, double reltimelength);
	void			BackupPath();
	void			RestorePath();

	// pulley: the beginning of path p, up to absolute time point t,
	// is inverted and transferred at the base of this path
	// path p is modified accordingly
	void 			Prefix(BranchSitePath* p, BranchSitePath* root, double abstime);

	string 			ToString(bool redundant = false);

	protected :

	Plink* init;
	Plink* last;
	int nsub;
	StateSpace* statespace;

	Plink* bkinit;
	Plink* bklast;
	int bknsub;
};




//-------------------------------------------------------------------------
//-------------------------------------------------------------------------
//	* Inline definitions
//-------------------------------------------------------------------------
//-------------------------------------------------------------------------


//-------------------------------------------------------------------------
//	* Plink
//-------------------------------------------------------------------------

inline Plink::Plink() : next(0), prev(0), state(0), rel_time(0)  {}
inline Plink::Plink(int instate, double inrel_time) : next(0), prev(0), state(instate), rel_time(inrel_time){}
inline Plink::~Plink()	{Splice();}

inline Plink* Plink::Prev() {return prev;}
inline Plink* Plink::Next() {return next;}

inline bool Plink::IsFirst()	{return prev == 0;}
inline bool Plink::IsLast() {return next == 0;}

inline void Plink::Insert(Plink* link)	{
	link->next = next;
	if (next)	{
		next->prev = link;
	}
	link->prev = this;
	next = link;
}

inline void Plink::Splice()	{
	if (prev)	{
		prev->next = next;
	}
	if (next)	{
		next->prev = prev;
	}
	prev = next = 0;
}

inline void Plink::SetState(int instate) {state = instate;}
inline void Plink::SetRelativeTime(double inrel_time) {rel_time = inrel_time;}
inline double Plink::GetRelativeTime() {return rel_time;}
inline int Plink::GetState() {return state;}

inline Plink* BranchSitePath::Init() {return init;}
inline Plink*  BranchSitePath::Last() {return last;}

//-------------------------------------------------------------------------
//	* BranchSitePath
//-------------------------------------------------------------------------

inline StateSpace* BranchSitePath::GetStateSpace() {if (!statespace) {cerr << "null pointer : BranchSitePath::statespace\n"; throw(0); } return statespace;}

inline void BranchSitePath::Append(int instate, double reltimelength)	{
	last->SetRelativeTime(reltimelength);
	Plink* link = new Plink(instate,0);
	last->Insert(link);
	last = link;
	nsub++;
}

inline void BranchSitePath::BKAppend(int instate, double reltimelength)	{
	bklast->SetRelativeTime(reltimelength);
	Plink* link = new Plink(instate,0);
	bklast->Insert(link);
	bklast = link;
	bknsub++;
}

inline int BranchSitePath::GetNsub()	{
	return nsub;
}

inline string BranchSitePath::GetState(Plink* link) {return GetStateSpace()->GetState(link->GetState());}
inline string BranchSitePath::GetCharInitState() {return GetState(init);}
inline string BranchSitePath::GetCharFinalState() {return GetState(last);}
inline int BranchSitePath::GetInitState() {return init->GetState();}
inline int BranchSitePath::GetFinalState() {return last->GetState();}

inline double BranchSitePath::GetAbsoluteTime(Plink* link) {return link->GetRelativeTime() * GetTotalTime();}


#endif // SITEPATH_H
