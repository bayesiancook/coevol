
#include <algorithm>

#include "Random.h"
#include "ProbModel.h"
#include "RandomSubMatrix.h"

bool DAGnode::initmode = true;


//-------------------------------------------------------------------------
//-------------------------------------------------------------------------
//	* DAGnode
//-------------------------------------------------------------------------
//-------------------------------------------------------------------------
const double MCMC::MAXDIFF = 1e-4;

DAGnode::~DAGnode()	{
	Detach();
}

int DAGnode::GetChildrenNumber()	{
	int tot = 0;
	for (clit i=down.begin(); i!=down.end(); i++)	{
		tot++;
	}
	return tot;
}

void DAGnode::Detach()	{
	while (up.size())	{
		DeregisterFrom(*up.begin());
	}
	while(down.size())	{
		(*down.begin())->DeregisterFrom(this);
	}
	/*
	for (clit i=up.begin(); i!=up.end(); i++)	{
		DeregisterFrom(*i);
	}
	*/
}

void DAGnode::DeregisterFrom(DAGnode* parent)	{
	if (parent)	{
		parent->down.erase(this);
		up.erase(parent);
	}
}

void DAGnode::Register(DAGnode* parent)	{
	if (parent)	{
		parent->down.insert(this);
		up.insert(parent);
	}
}

void DAGnode::RecursiveRegister(ProbModel* model)	{
	clit i=up.begin();
	while ((i!=up.end()) && (*i)->flag)  i++;
	bool up_ok = (i == up.end());

	for (clit i=up.begin(); i!=up.end(); i++)	{
		up_ok &= (*i)->flag;
	}
	if (up_ok)	{
		model->Register(this);
		flag = true;
		for (clit i=down.begin(); i!=down.end(); i++)	{
			(*i)->RecursiveRegister(model);
		}
	}
}

bool DAGnode::CheckUpdateFlags()	{
	bool ret = flag;
	if (! flag)	{
		cerr << "flag error : " << GetName() << '\n';
	}
	for (clit i=down.begin(); i!=down.end(); i++)	{
		ret &= (*i)->CheckUpdateFlags();
	}
	return ret;
}

//-------------------------------------------------------------------------
//-------------------------------------------------------------------------
//	* Rnode
//-------------------------------------------------------------------------
//-------------------------------------------------------------------------

//-------------------------------------------------------------------------
//	* corrupt
//-------------------------------------------------------------------------

void Rnode::Corrupt(bool bk)	{
	value_updated = false;
	localCorrupt(bk);
	for (clit i=down.begin(); i!=down.end(); i++)	{
		(*i)->NotifyCorrupt(bk);
	}
}

void Rnode:: NotifyCorrupt(bool bk)	{
	localCorrupt(bk);
}

void Rnode::localCorrupt(bool bk)	{
	if (bk)	{
		bklogprob = logprob;
	}
	flag = false;
}

void Rnode::FullCorrupt(map<DAGnode*,int>& m)	{
	localCorrupt(true);
	if (m.find(this) == m.end())	{
		m[this] = 1;
		for (clit i=down.begin(); i!=down.end(); i++)	{
			(*i)->FullCorrupt(m);
		}
	}
}

//-------------------------------------------------------------------------
//	* update
//-------------------------------------------------------------------------

double Rnode::Update()	{
	double ret = 0;
	if (! flag)	{
		clit i=up.begin();
		while ((i!=up.end()) && ((*i)->isValueUpdated())) i++;
		bool up_ok = (i == up.end());
		/*
		bool up_ok = true;
		for (clit i=up.begin(); i!=up.end(); i++)	{
			up_ok &= (*i)->isValueUpdated();
			// up_ok &= (*i)->flag;
		}
		*/
		if (up_ok)	{
			ret = localUpdate();
			value_updated = true;
			for (clit i=down.begin(); i!=down.end(); i++)	{
				ret += (*i)->NotifyUpdate();
			}
		}
	}
	return ret;
}

double Rnode::NotifyUpdate()	{
	double ret = 0;
	if (! flag)	{
		clit i=up.begin();
		while ((i!=up.end()) && ((*i)->isValueUpdated())) i++;
		bool up_ok = (i == up.end());
		/*
		bool up_ok = true;
		for (clit i=up.begin(); i!=up.end(); i++)	{
			up_ok &= (*i)->isValueUpdated();
			// up_ok &= (*i)->flag;
		}
		if (GetName() == "BD Chrono")	{
			cerr << "BD::NotifyUpdate : " << up_ok << '\n';
		}
		*/
		if (up_ok)	{
			ret = localUpdate();
			if (! value_updated)	{
				value_updated = true;
				for (clit i=down.begin(); i!=down.end(); i++)	{
					ret += (*i)->NotifyUpdate();
				}
			}
		}
	}
	return ret;
}

double Rnode::localUpdate()	{
	logprob =  logProb();
	flag = true;
	return logprob-bklogprob;
}

double Rnode::FullUpdate(bool check)	{
	double ret = 0;
	if (! flag)	{
		clit i=up.begin();
		while ((i!=up.end()) && ((*i)->isValueUpdated())) i++;
		bool up_ok = (i == up.end());
		/*
		bool up_ok = true;
		for (clit i=up.begin(); i!=up.end(); i++)	{
			up_ok &= (*i)->isValueUpdated();
			// up_ok &= (*i)->flag;
		}
		*/
		if (up_ok)	{
			ret = localUpdate();
			value_updated = true;
			if ((fabs(ret)>MCMC::MAXDIFF) && check)	{
				cerr << "NON ZERO CHECKSUM : " << GetName() << '\n';
				cerr << "number of parents : " << up.size() << '\n';
				throw CheckSumException(ret);
			}
			for (clit i=down.begin(); i!=down.end(); i++)	{
				ret += (*i)->FullUpdate(check);
			}
		}
	}
	return ret;

}

/*
void Rnode::Initialise()	{

	if (! flag)	{
		bool up_ok = true;
		for (clit i=up.begin(); i!=up.end(); i++)	{
			up_ok &= (*i)->flag;
		}
		if (up_ok)	{
			Sample();
			localUpdate();
			for (clit i=down.begin(); i!=down.end(); i++)	{
				(*i)->Initialise();
			}
		}
	}
}
*/

//-------------------------------------------------------------------------
//	* restore
//-------------------------------------------------------------------------

void Rnode::Restore()	{
	if (! flag)	{
		clit i=up.begin();
		while ((i!=up.end()) && ((*i)->isValueUpdated())) i++;
		bool up_ok = (i == up.end());
		/*
		bool up_ok = true;
		for (clit i=up.begin(); i!=up.end(); i++)	{
			up_ok &= (*i)->isValueUpdated();
			// up_ok &= (*i)->flag;
		}
		*/
		if (up_ok)	{
			localRestore();
			value_updated = true;
			for (clit i=down.begin(); i!=down.end(); i++)	{
				(*i)->NotifyRestore();
			}
		}
	}
}

void Rnode::NotifyRestore()	{
	if (! flag)	{
		clit i=up.begin();
		while ((i!=up.end()) && ((*i)->isValueUpdated())) i++;
		bool up_ok = (i == up.end());
		/*
		bool up_ok = true;
		for (clit i=up.begin(); i!=up.end(); i++)	{
			up_ok &= (*i)->isValueUpdated();
			// up_ok &= (*i)->flag;
		}
		*/
		if (up_ok)	{
			localRestore();
			if (! value_updated)	{
				value_updated = true;
				for (clit i=down.begin(); i!=down.end(); i++)	{
					(*i)->NotifyRestore();
				}
			}
		}
	}
}

void Rnode::localRestore()	{
	logprob = bklogprob;
	/*
	if (logprob != logProb())	{
		cerr << "error in local restore\n";
		exit(1);
	}
	*/
	flag = true;
}


//-------------------------------------------------------------------------
//-------------------------------------------------------------------------
//	* Dnode
//-------------------------------------------------------------------------
//-------------------------------------------------------------------------

//-------------------------------------------------------------------------
//	* corrupt
//-------------------------------------------------------------------------

void Dnode::Corrupt(bool bk)	{
	localCorrupt(bk);
	for (clit i=down.begin(); i!=down.end(); i++)	{
		(*i)->NotifyCorrupt(bk);
	}
}

void Dnode:: NotifyCorrupt(bool bk)	{
	Corrupt(bk);
}

void Dnode::localCorrupt(bool bk)	{
	flag = false;
}

void Dnode::FullCorrupt(map<DAGnode*,int>& m)	{
	localCorrupt(true);
	if (m.find(this) == m.end())	{
		m[this] = 1;
		for (clit i=down.begin(); i!=down.end(); i++)	{
			(*i)->FullCorrupt(m);
		}
	}
}

//-------------------------------------------------------------------------
//	* update
//-------------------------------------------------------------------------

double Dnode::Update()	{
	double ret = 0;
	if (! flag)	{
		clit i=up.begin();
		while ((i!=up.end()) && ((*i)->isValueUpdated())) i++;
		bool up_ok = (i == up.end());
		/*
		bool up_ok = true;
		for (clit i=up.begin(); i!=up.end(); i++)	{
			up_ok &= (*i)->isValueUpdated();
			// up_ok &= (*i)->flag;
		}
		*/
		if (up_ok)	{
			ret = localUpdate();
			for (clit i=down.begin(); i!=down.end(); i++)	{
				ret += (*i)->NotifyUpdate();
			}
		}
	}
	return ret;
}

double Dnode::NotifyUpdate()	{
	return Update();
}

double Dnode::localUpdate()	{
	specialUpdate();
	flag = true;
	return 0;
}

double Dnode::FullUpdate(bool check)	{
	double ret = 0;
	if (! flag)	{
		clit i=up.begin();
		while ((i!=up.end()) && ((*i)->isValueUpdated())) i++;
		bool up_ok = (i == up.end());
		/*
		bool up_ok = true;
		for (clit i=up.begin(); i!=up.end(); i++)	{
			up_ok &= (*i)->isValueUpdated();
			// up_ok &= (*i)->flag;
		}
		*/
		if (up_ok)	{
			ret = localUpdate();
			for (clit i=down.begin(); i!=down.end(); i++)	{
				ret += (*i)->FullUpdate(check);
			}
		}
	}
	return ret;
}

/*
void Dnode::Initialise()	{
	if (! flag)	{
		bool up_ok = true;
		for (clit i=up.begin(); i!=up.end(); i++)	{
			up_ok &= (*i)->flag;
		}
		if (up_ok)	{
			localUpdate();
			for (clit i=down.begin(); i!=down.end(); i++)	{
				(*i)->Initialise();
			}
		}
	}
}
*/

//-------------------------------------------------------------------------
//	* restore
//-------------------------------------------------------------------------

void Dnode::Restore()	{
	if (! flag)	{
		clit i=up.begin();
		while ((i!=up.end()) && ((*i)->isValueUpdated())) i++;
		bool up_ok = (i == up.end());
		/*
		bool up_ok = true;
		for (clit i=up.begin(); i!=up.end(); i++)	{
			up_ok &= (*i)->isValueUpdated();
			// up_ok &= (*i)->flag;
		}
		*/
		if (up_ok)	{
			localRestore();
			for (clit i=down.begin(); i!=down.end(); i++)	{
				(*i)->NotifyRestore();
			}
		}
	}
}

void Dnode::NotifyRestore()	{
	Restore();
}

void Dnode::localRestore()	{
	flag = true;
}


//-------------------------------------------------------------------------
//-------------------------------------------------------------------------
//	* Mnode
//-------------------------------------------------------------------------
//-------------------------------------------------------------------------

//-------------------------------------------------------------------------
//	* corrupt
//-------------------------------------------------------------------------

void Mnode::Corrupt(bool bk)	{
	flag = false;
	for (clit i=down.begin(); i!=down.end(); i++)	{
		(*i)->Corrupt(bk);
	}
}

//-------------------------------------------------------------------------
//	* update
//-------------------------------------------------------------------------

double Mnode::Update()	{
	flag = true;
	double ret = 0;
	for (clit i=down.begin(); i!=down.end(); i++)	{
		ret += (*i)->Update();
	}
	return ret;
}

//-------------------------------------------------------------------------
//	* restore
//-------------------------------------------------------------------------

void Mnode::Restore()	{
	flag = true;
	for (clit i=down.begin(); i!=down.end(); i++)	{
		(*i)->Restore();
	}
}

