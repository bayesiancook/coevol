
#ifndef VAR_H
#define VAR_H

#include "DAGnode.h"

template<class T> class Var : public virtual DAGnode , public T{

	public:

	Var() {}

	Var(const T& from) : T(from) , bkvalue(from) {}

	virtual ~Var() {}

	Var& operator=(const T& from)	{
		return T::operator=(from);
	}

	void  localcorrupt(bool bk)	{
		if (bk)	{
			this->bkvalue = *this;
		}
	}

	void  localrestore()	{
		T::operator=(bkvalue);
	}

	// returns the current value
	const T& val() {return *this;}

	// const T& val2() {return *this;}

	// sets the new value
	void setval(const T& inval) {T::operator=(inval);}

	virtual void Register(DAGnode* in) {
		DAGnode::Register(in);
	}

	protected:

	T bkvalue;
};

template<class T> class Rvar : public Var<T>, public Rnode	{

	public:

	Rvar() {}

	Rvar(const T& from) : Var<T>(from) {}

	virtual double ProposeMove(double tuning)	{
		return T::ProposeMove(tuning);
	}

	void ClampAt(const T& inval)	{
		T::operator=(inval);
		Clamp();
	}

	virtual void	Corrupt(bool bk)	{
		Var<T>::localcorrupt(bk);
		Rnode::Corrupt(bk);
	}

	virtual void 	Restore()	{
		Var<T>::localrestore();
		Rnode::Restore();
	}

	void RestoreBackup()	{
		Var<T>::localrestore();
		// value_updated = true;
	}

};

template<class T> class Dvar : public Var<T>, public Dnode {

	public:

	Dvar() {}

	Dvar(const T& from) : Var<T>(from) {}

	void	localCorrupt(bool bk)	{
		Var<T>::localcorrupt(bk);
		Dnode::localCorrupt(bk);
	}

	void 	localRestore()	{
		Var<T>::localrestore();
		Dnode::localRestore();
	}

	// returns the current value
	const T& val() {
		if (! isUpdated())	{
			if (! DAGnode::initmode)	{
				cerr << "error : corrupting Dvar\n";
			}
			localUpdate();
		}
		return (*this);
	}
};


template <class T> class Const : public Dvar<T>	{

	public:

	Const(const T& from) : Dvar<T>(from) {}

	void specialUpdate() {}

};

template<> class Rvar<void> : public virtual Rnode {};
template<> class Dvar<void> : public virtual Dnode {};

#endif

