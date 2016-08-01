
#ifndef VALARRAY_H
#define VALARRAY_H

#include "BaseType.h"

class AbstractArray	{

	public:

	virtual ~AbstractArray() {}

	virtual int GetSize() = 0;
};

template <class V> class VarArray : public virtual AbstractArray {

	public:

	virtual ~VarArray() {}

	virtual Var<V>* GetVal(int site) = 0;

	void SetName(string inname)	{
		for (int i=0; i<GetSize(); i++)	{
			GetVal(i)->SetName(inname);
		}
	}

	void ToStream(ostream& os)	{
		for (int i=0; i<GetSize(); i++)	{
			os << *(GetVal(i)) << '\t';
		}
		os << '\n';
	}

	void FromStream(istream& is)	{
		for (int i=0; i<GetSize(); i++)	{
			is >> *(GetVal(i));
		}
	}

	friend ostream& operator<<(ostream& os, VarArray<V>& a)	{
		a.ToStream(os);
		return os;
	}

	friend istream& operator>>(istream& is, VarArray<V>& a)  {
		a.FromStream(is);
		return is;
	}

	void RegisterChild(DAGnode* node)	{
		for (int i=0; i<this->GetSize(); i++)	{
			node->Register(this->GetVal(i));
		}
	}

	virtual void RegisterArray(DAGnode* in)	{
		for (int i=0; i<GetSize(); i++)	{
			GetVal(i)->Register(in);
		}
	}

};

template<class V> class RandomVarArray : public VarArray<V>, public MCMC {

	public:
};

template<> class RandomVarArray<PosReal> : public VarArray<PosReal>, public MCMC, public Multiplicative {

	public:

	virtual void Register(DAGnode* in)	{
		RegisterArray(in);
	}

	int ScalarMultiplication(double d)	{
		for (int i=0; i<GetSize(); i++)	{
			GetVal(i)->ScalarMultiplication(d);
		}
		return GetSize();
	}
};


template <class T> class _ValPtrArray : public virtual AbstractArray {

	public:

	_ValPtrArray(int insize)	{
		size = insize;
		array = new T*[size];
		for (int i=0; i<GetSize(); i++)	{
			array[i] = 0;
		}
	}

	virtual ~_ValPtrArray()	{
		for (int i=0; i<GetSize(); i++)	{
			delete array[i];
			array[i] = 0;
		}
		delete[] array;
	}

	int GetSize()	{
		return size;
	}

	T* GetVal(int site)	{
		return array[site];
	}

	T* operator[](int site)	{
		return GetVal(site);
	}
	void SetVal(int site, T* in)	{
		array[site] = in;
	}

	void Create()	{
		for (int i=0; i<GetSize(); i++)	{
			array[i] = CreateVal(i);
		}
	}

	void Delete()	{
		for (int i=0; i<GetSize(); i++)	{
			delete array[i];
			array[i] = 0;
		}
	}

	/*
	void ToStream(ostream& os)	{
		for (int i=0; i<GetSize(); i++)	{
			os << *(array[i]) << '\t';
		}
		os << '\n';
	}

	void FromStream(istream& is)	{
		for (int i=0; i<GetSize(); i++)	{
			is >> *(array[i]);
		}
	}

	friend ostream& operator<<(ostream& os, _ValPtrArray<T>& a)	{
		a.ToStream(os);
		return os;
	}

	friend istream& operator>>(istream& is, _ValPtrArray<T>& a)  {
		a.FromStream(is);
		return is;
	}
	*/

	protected:

	virtual T* CreateVal(int site) = 0;

	int size;
	T** array;

};


template<class T> class ValPtrArray: public _ValPtrArray<T> {

	ValPtrArray(int insize) : _ValPtrArray<T>(insize) {}
};

template<class V> class ValPtrArray<Rvar<V> > : public _ValPtrArray<Rvar<V> >, public virtual RandomVarArray<V> {

	public:

	ValPtrArray<Rvar<V> >(int insize) : _ValPtrArray<Rvar<V> >(insize) {}

	virtual Rvar<V>* GetVal(int site)	{
		return _ValPtrArray<Rvar<V> >::GetVal(site);
	}
};


template<class V> class ValPtrArray<Dvar<V> > : public _ValPtrArray<Dvar<V> >, public virtual VarArray<V>	{

	public:

	ValPtrArray<Dvar<V> >(int insize) : _ValPtrArray<Dvar<V> >(insize) {}

	virtual Dvar<V>* GetVal(int site)	{
		return _ValPtrArray<Dvar<V> >::GetVal(site);
	}
};


#endif
