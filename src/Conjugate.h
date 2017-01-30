
#ifndef CONJUGATE_H
#define CONJUGATE_H

using namespace std;
#include "MCMC.h"
#include <set>

class Rnode;

typedef set<DAGnode*>::const_iterator clit;

class SemiConjPrior {

	public:

	// to be implemented in non-abstract subclasses

	//
	virtual void ResetSufficientStatistic() = 0;
	virtual void SaveSufficientStatistic() = 0;
	virtual void RestoreSufficientStatistic() = 0;

	// this function
	virtual double SuffStatLogProb() = 0;


	SemiConjPrior()	{
		suffstat_flag = false;
		active_flag = false;
		corrupt_counter = 0;
	}

	virtual ~SemiConjPrior() {}

	// are already implemented in SemiConjugatePrior template
	virtual void InactivateSufficientStatistic() = 0;
	virtual void ConjugateNotifyCorrupt(bool bk) = 0;
	virtual double ConjugateNotifyUpdate() = 0;
	virtual void ConjugateNotifyRestore() = 0;


	virtual void CorruptSufficientStatistic()	{
		suffstat_flag = false;
	}

	bool isActive()	{
		return active_flag;
	}

	protected:

	bool suffstat_flag;
	bool active_flag;
	int  corrupt_counter;
};

class ConjSampling	{

	public:

	typedef set<SemiConjPrior*>::const_iterator cplit;

	virtual ~ConjSampling() {}

	// to be implemented in non-abstract subclasses
	virtual void AddSufficientStatistic(SemiConjPrior* parent) = 0;

	//
	// these 3 methods are called during a corrupt/update/restore cycle
	// after the local work (i.e. localCorrupt/localUpdate/localRestore) has been done
	//
	// with these methods,the conjugate sampling node notifies its conjugate prior parent
	// that the sufficient statistic and the associated log prob should be backed-up/recomputed/restored.
	//
	virtual void ConjugateCorrupt(bool bk)	{
		for (cplit i=conjugate_up.begin(); i!=conjugate_up.end(); i++)	{
			if ((*i)->isActive())	{
				(*i)->ConjugateNotifyCorrupt(bk);
			}
		}
	}

	virtual double ConjugateUpdate(bool& active)	{
		double ret = 0;
		for (cplit i=conjugate_up.begin(); i!=conjugate_up.end(); i++)	{
			if ((*i)->isActive())	{
				active = true;
				ret += (*i)->ConjugateNotifyUpdate();
			}
		}
		return ret;
	}

	virtual void ConjugateRestore(bool& active)	{
		for (cplit i=conjugate_up.begin(); i!=conjugate_up.end(); i++)	{
			if ((*i)->isActive())	{
				active = true;
				(*i)->ConjugateNotifyRestore();
			}
		}
	}

	bool isUpActive()	{
		cplit i=conjugate_up.begin();
		while ((i != conjugate_up.end()) && (!(*i)->isActive())) i++;
		return (i != conjugate_up.end());
	}

	//
	// when a conjugate prior nodes activates its conjugate mode
	// it notifies its conjugate sampling child nodes by calling this method
	// the child nodes then check that other "conjugate canals"
	// (other conjugate pairs they are involved in as the conj sampling partner)
	// are not also activated (otherwise, they shut them down)
	//
	void NotifyActivateSufficientStatistic(SemiConjPrior* from)	{
		// should also check that the type matches
		for (cplit i=conjugate_up.begin(); i!=conjugate_up.end(); i++)	{
			if (*(i) != from)	{
				(*i)->InactivateSufficientStatistic();
			}
		}
	}

	protected:

	set<SemiConjPrior*> conjugate_up;

};


template<class T> class DSemiConjugatePrior : public virtual Dvar<T>, public SemiConjPrior {

	public:


	// subclasses will diamond-inherit from Rvar<T>
	// they will HAVE to reprogram logProb()
	virtual double logProb() = 0;

	// Corrupt/Update/Restore methods

	//
	// the modification, compared to a non-conjugate Dvar<T> node
	// is that the corrupt/update/restore wave should not be propagated to downstream nodes
	// (all of which are already accounted for through the sufficient statistic)
	//

	virtual void FullCorrupt(map<DAGnode*,int>& m)	{
		if (! isActive())	{
			Dvar<T>::Corrupt(true);
		}
		else	{
			/*
			if (bk)	{
				this->bkvalue = *this;
			}
			*/
			this->localCorrupt(true);
			bklogprob = logprob;
			// ???
			// value_updated = false;
		}
	}

	virtual void Corrupt(bool bk)	{
		if (! isActive())	{
			Dvar<T>::Corrupt(bk);
		}
		else	{
			/*
			if (bk)	{
				this->bkvalue = *this;
			}
			*/
			this->localCorrupt(bk);
			if (bk)	{
				bklogprob = logprob;
			}
			// ???
			// value_updated = false;
		}
	}

	virtual double FullUpdate(bool check)	{
		return Update();
	}

	virtual double Update()	{
		double ret = 0;
		if (! isActive())	{
			ret = Dvar<T>::Update();
		}
		else	{
			if (! this->flag)	{
				bool up_ok = true;
				for (clit i=this->up.begin(); i!=this->up.end(); i++)	{
					up_ok &= (*i)->isUpdated();
				}
				if (up_ok)	{
					this->localUpdate();
					logprob = this->logProb();
					ret += logprob - bklogprob;
				}
			}
		}
		return ret;
	}

	virtual void Restore()	{
		if (! isActive())	{
			Dvar<T>::Restore();
		}
		else	{
			// setval(this->bkvalue);
			if (! this->isUpdated())	{
				bool up_ok = true;
				for (clit i=this->up.begin(); i!=this->up.end(); i++)	{
					up_ok &= (*i)->isUpdated();
				}
				if (up_ok)	{
					this->localRestore();
					logprob = bklogprob;
					/*
					for (clit i=this->down.begin(); i!=this->down.end(); i++)	{
						(*i)->NotifyRestore();
					}
					*/
				}
			}
		}
	}

	virtual void NotifyCorrupt(bool bk)	{
		if (! isActive())	{
			Dvar<T>::NotifyCorrupt(bk);
		}
		else	{
			/*
			if (bk)	{
				this->bkvalue = *this;
			}
			*/
			this->localCorrupt(bk);
			if (bk)	{
				bklogprob = logprob;
			}
			// ???
			// value_updated = false;
		}
	}

	virtual double NotifyUpdate()	{
		double ret = 0;
		if (! isActive())	{
			ret = Dvar<T>::NotifyUpdate();
		}
		else	{
			if (! this->flag)	{
				bool up_ok = true;
				for (clit i=this->up.begin(); i!=this->up.end(); i++)	{
					up_ok &= (*i)->isUpdated();
				}
				if (up_ok)	{
					this->localUpdate();
					logprob = this->logProb();
					ret += logprob - bklogprob;
				}
			}
		}
		return ret;
	}

	virtual void NotifyRestore()	{
		if (! isActive())	{
			Dvar<T>::NotifyRestore();
		}
		else	{
			// setval(this->bkvalue);
			if (! this->isUpdated())	{
				bool up_ok = true;
				for (clit i=this->up.begin(); i!=this->up.end(); i++)	{
					up_ok &= (*i)->isUpdated();
				}
				if (up_ok)	{
					this->localRestore();
					logprob = bklogprob;
					/*
					for (clit i=this->down.begin(); i!=this->down.end(); i++)	{
						(*i)->NotifyRestore();
					}
					*/
				}
			}
		}
	}

	// Conjugate Notify methods
	//
	// when a child (conjugate sampling) node makes a corrupt/update/restore cycle
	// it notifies its (conjugate prior) parent using these 3 methods
	// the parent then knows that the sufficient statistics are corrupted,
	// it backups/recomputes/restores them accordingly
	// and backups/updates/restores its own probability
	//
	virtual double ConjugateNotifyUpdate()	{
		cerr << "conjugate notify\n";
		exit(1);
		corrupt_counter--;
		if (! corrupt_counter)	{
			this->flag = false;
			ComputeSufficientStatistic();
			// return  this->localUpdate();
			this->localUpdate();
			logprob = this->logProb();
			return logprob - bklogprob;
		}
		return 0;
	}

	virtual void ConjugateNotifyCorrupt(bool bk)	{
		cerr << "conjugate notify corrupt\n";
		exit(1);
		if (! isActive())	{
			cerr << "error in Conjugate corrupt\n";
			exit(1);
		}
		if (! corrupt_counter)	{
			if (bk)	{
				SaveSufficientStatistic();
			}
			CorruptSufficientStatistic();
			this->localCorrupt(bk);
			if (bk)	{
				bklogprob = logprob;
			}
			this->flag = true;
		}
		corrupt_counter++;
	}

	virtual void ConjugateNotifyRestore()	{
		cerr << "conjugate notify\n";
		exit(1);
		corrupt_counter--;
		if (! corrupt_counter)	{
			RestoreSufficientStatistic();
			this->localRestore();
			logprob = bklogprob;
		}
	}

	//
	// sufficient statistics methods
	//
	// ActivateSufficientStatistic, well, activates the sufficient statistics,
	// but also, notifies the children (so that they can inactivate other potential conjugate prior nodes)
	// this is useful for nodes that are involved in several conjugate pairs (each time as the conjugate sampling pole),
	// with several of their parents: in such cases, one should avoid potential interferences
	//
	// InactivateSufficientStatistics restores the plain MCMC mode
	//
	// ComputeSufficientStatistic gathers the sufficient statistics of all its (conjugate) child nodes,
	// by requesting them one by one
	//
	virtual void ActivateSufficientStatistic()	{
			if (! active_flag)	{

				for (clit i=this->down.begin(); i!=this->down.end(); i++)	{
					ConjSampling* p = dynamic_cast<ConjSampling*> (*i);
					if (!p)	{
						Dnode* q = dynamic_cast<Dnode*> (*i);
						if (! q)	{
							cerr << "error : non conjugate child nodes, cannot activate sufficient statistic\n";
							cerr << (*i) << '\n';
							exit(1);
						}
					}
					if (p)	{
						p->NotifyActivateSufficientStatistic(this);
					}
				}
				suffstat_flag = false;
				active_flag = true;
				ComputeSufficientStatistic();
				this->localUpdate();
				logprob = logProb();
			}
	}

	virtual void InactivateSufficientStatistic()	{
		if (active_flag)	{
			active_flag = false;
			suffstat_flag = false;
			Corrupt(false);
			Update();
		}
	}

	virtual void ComputeSufficientStatistic()	{
		if (active_flag)	{
			if (! suffstat_flag)	{
				ResetSufficientStatistic();
				for (clit i=this->down.begin(); i!=this->down.end(); i++)	{
					ConjSampling* p = dynamic_cast<ConjSampling*> (*i);
					if (p)	{
						p->AddSufficientStatistic(this);
					}
				}
				suffstat_flag = true;
			}
		}
	}

	private:

	double bklogprob;
	double logprob;

};

template<class T> class SemiConjugatePrior : public virtual Rvar<T>, public SemiConjPrior {

	public:


	// subclasses will diamond-inherit from Rvar<T>
	// they will HAVE to reprogram logProb()
	virtual double logProb() = 0;

	// Corrupt/Update/Restore methods

	//
	// the modification, compared to a non-conjugate Rvar<T> node
	// is that the corrupt/update/restore wave should not be propagated to downstream nodes
	// (all of which are already accounted for through the sufficient statistic)
	//

	virtual void Corrupt(bool bk)	{
		if (! isActive())	{
			Rvar<T>::Corrupt(bk);
		}
		else	{
			if (bk)	{
				this->bkvalue = *this;
			}
			this->localCorrupt(bk);
			// ???
			// value_updated = false;
		}
	}

	virtual double Update()	{
		double ret = 0;
		if (! isActive())	{
			ret = Rvar<T>::Update();
		}
		else	{
			if (! this->flag)	{
				bool up_ok = true;
				for (clit i=this->up.begin(); i!=this->up.end(); i++)	{
					up_ok &= (*i)->isUpdated();
				}
				if (up_ok)	{
					ret = this->localUpdate();
				}
			}
		}
		return ret;
	}

	virtual void Restore()	{
		if (! isActive())	{
			Rvar<T>::Restore();
		}
		else	{
			this->setval(this->bkvalue);
			if (! this->isUpdated())	{
				bool up_ok = true;
				for (clit i=this->up.begin(); i!=this->up.end(); i++)	{
					up_ok &= (*i)->isUpdated();
				}
				if (up_ok)	{
					this->localRestore();
					for (clit i=this->down.begin(); i!=this->down.end(); i++)	{
						(*i)->NotifyRestore();
					}
				}
			}
		}
	}

	// Conjugate Notify methods
	//
	// when a child (conjugate sampling) node makes a corrupt/update/restore cycle
	// it notifies its (conjugate prior) parent using these 3 methods
	// the parent then knows that the sufficient statistics are corrupted,
	// it backups/recomputes/restores them accordingly
	// and backups/updates/restores its own probability
	//
	virtual double ConjugateNotifyUpdate()	{
		corrupt_counter--;
		if (! corrupt_counter)	{
			this->flag = false;
			ComputeSufficientStatistic();
			return  this->localUpdate();
		}
		return 0;
	}

	virtual void ConjugateNotifyCorrupt(bool bk)	{
		if (! isActive())	{
			cerr << "error in Conjugate corrupt\n";
			exit(1);
		}
		if (! corrupt_counter)	{
			if (bk)	{
				SaveSufficientStatistic();
			}
			CorruptSufficientStatistic();
			this->localCorrupt(bk);
			this->flag = true;
		}
		corrupt_counter++;
	}

	virtual void ConjugateNotifyRestore()	{
		corrupt_counter--;
		if (! corrupt_counter)	{
			RestoreSufficientStatistic();
			this->localRestore();
		}
	}

	//
	// sufficient statistics methods
	//
	// ActivateSufficientStatistic, well, activates the sufficient statistics,
	// but also, notifies the children (so that they can inactivate other potential conjugate prior nodes)
	// this is useful for nodes that are involved in several conjugate pairs (each time as the conjugate sampling pole),
	// with several of their parents: in such cases, one should avoid potential interferences
	//
	// InactivateSufficientStatistics restores the plain MCMC mode
	//
	// ComputeSufficientStatistic gathers the sufficient statistics of all its (conjugate) child nodes,
	// by requesting them one by one
	//
	virtual void ActivateSufficientStatistic()	{
		if (! this->isClamped())	{
			if (! active_flag)	{

				for (clit i=this->down.begin(); i!=this->down.end(); i++)	{
					ConjSampling* p = dynamic_cast<ConjSampling*> (*i);
					/*if (!p)	{
						Dnode* q = dynamic_cast<Dnode*> (*i);
						if (! q)	{
							cerr << "error : non conjugate child nodes, cannot activate sufficient statistic\n";
							cerr << (*i) << '\n';
							exit(1);
						}
					}*/
					if (p)	{
						p->NotifyActivateSufficientStatistic(this);
					}
				}
				suffstat_flag = false;
				active_flag = true;
				ComputeSufficientStatistic();
				this->localUpdate();
			}
		}
	}

	virtual void InactivateSufficientStatistic()	{
		if (active_flag)	{
			active_flag = false;
			suffstat_flag = false;
			Corrupt(false);
			Update();
		}
	}

	virtual void ComputeSufficientStatistic()	{
		if (active_flag)	{
			if (! suffstat_flag)	{
				ResetSufficientStatistic();
				for (clit i=this->down.begin(); i!=this->down.end(); i++)	{
					ConjSampling* p = dynamic_cast<ConjSampling*> (*i);
					if (p)	{
						p->AddSufficientStatistic(this);
					}
				}
				suffstat_flag = true;
			}
		}
	}

};

template <class R> class ConjugatePrior : public virtual SemiConjugatePrior<R>	{

	public:

	ConjugatePrior()	{
		integrated_flag = false;
	}

	virtual ~ConjugatePrior() {}

	// to be implemented in non-abstract subclasses
	virtual void GibbsResample() = 0;

	void InactivateSufficientStatistic()	{
		if (this->active_flag)	{
			if (isIntegrated())	{
				GibbsResample();
				/*
				this->CorruptSufficientStatistic();
				this->ComputeSufficientStatistic();
				GibbsResample();
				*/
				integrated_flag = false;
			}
			this->active_flag = false;
			this->suffstat_flag = false;
			// this->localUpdate();
			this->Corrupt(false);
			this->Update();
		}
	}

	bool isIntegrated()	{
		return integrated_flag;
	}

	void Integrate()	{
		if (! this->isClamped())	{
			integrated_flag = true;
			this->ActivateSufficientStatistic();
		}
	}

	void Resample()	{
		InactivateSufficientStatistic();
	}

	private:

	bool integrated_flag;
};


template <class T> class ConjugateSampling : public virtual Rvar<T>, public ConjSampling	{

	public:

	virtual void Corrupt(bool bk)	{
		Rvar<T>::Corrupt(bk);
		ConjugateCorrupt(bk);
	}

	virtual void NotifyCorrupt(bool bk)	{
		Rvar<T>::NotifyCorrupt(bk);
		ConjugateCorrupt(bk);
	}

	virtual double Update()	{
		double ret = 0;
		if (! this->isUpdated())	{
			bool up_ok = true;
			for (clit i=this->up.begin(); i!=this->up.end(); i++)	{
				up_ok &= (*i)->isValueUpdated();
				// up_ok &= (*i)->isUpdated();
			}
			if (up_ok)	{
				bool active = false;
				ret += ConjugateUpdate(active);
				if (! active)	{
					ret += this->localUpdate();
				}
				else	{
					this->flag = true;
				}
				this->value_updated = true;
				for (clit i=this->down.begin(); i!=this->down.end(); i++)	{
					ret += (*i)->NotifyUpdate();
				}
			}
		}
		return ret;
	}

	virtual void Restore()	{
		if (! this->isUpdated())	{
			Rvar<T>::RestoreBackup();
			bool up_ok = true;
			for (clit i=this->up.begin(); i!=this->up.end(); i++)	{
				up_ok &= (*i)->isValueUpdated();
				// up_ok &= (*i)->isUpdated();
			}
			if (up_ok)	{
				Rvar<T>::RestoreBackup();
				bool active = false;
				ConjugateRestore(active);
				if (! active)	{
					this->localRestore();
				}
				else	{
					this->flag = true;
				}
				this->value_updated = true;
				for (clit i=this->down.begin(); i!=this->down.end(); i++)	{
					(*i)->NotifyRestore();
				}
			}
		}
	}

	virtual double NotifyUpdate()	{
		double ret = 0;
		if (! this->isUpdated())	{
			bool up_ok = true;
			for (clit i=this->up.begin(); i!=this->up.end(); i++)	{
				up_ok &= (*i)->isValueUpdated();
				// up_ok &= (*i)->isUpdated();
			}
			if (up_ok)	{
				bool active = false;
				ret += ConjugateUpdate(active);
				if (! active)	{
					ret += this->localUpdate();
					if (! this->value_updated)	{
						this->value_updated = true;
						for (clit i=this->down.begin(); i!=this->down.end(); i++)	{
							ret += (*i)->NotifyUpdate();
						}
					}
				}
				else	{
					this->flag = true;
					this->value_updated = true;
				}
			}
		}
		return ret;
	}

	virtual void NotifyRestore()	{
		if (! this->isUpdated())	{
			bool up_ok = true;
			for (clit i=this->up.begin(); i!=this->up.end(); i++)	{
				up_ok &= (*i)->isValueUpdated();
				// up_ok &= (*i)->isUpdated();
			}
			if (up_ok)	{
				bool active = false;
				ConjugateRestore(active);
				if (! active)	{
					if (! this->value_updated)	{
						Rvar<T>::RestoreBackup();
					}
					this->localRestore();
					if (! this->value_updated)	{
						this->value_updated = true;
						for (clit i=this->down.begin(); i!=this->down.end(); i++)	{
							(*i)->NotifyRestore();
						}
					}
				}
				else	{
					this->flag = true;
					this->value_updated = true;
				}
			}
		}
	}
};


#endif // CONJUGATE_H

