#ifndef VAR_H
#define VAR_H

#include "DAGnode.hpp"
#include <iostream>// FIXME could be removed if implementation is moved into a file
#include <sstream>


template<class T> class Var : public virtual DAGnode , public T{
public:
  Var() {}
  Var(const T& from) : T(from) , bkvalue(from) {}
  virtual ~Var() {}

  inline Var& operator=(const T& from) { return T::operator=(from); }

  inline void localcorrupt(bool bk)	{ if (bk) this->bkvalue = *this; }
  inline void localrestore()	{ T::operator=(bkvalue); }

  // returns the current value
  inline const T& val() { return *this; }

  // sets the new value
  inline void setval(const T& inval) { T::operator=(inval); }

  inline virtual void Register(DAGnode* in) { DAGnode::Register(in); }

protected:
  T bkvalue;
};


template<class T> class Rvar : public Var<T>, public Rnode	{
public:
  Rvar() {}
  Rvar(const T& from) : Var<T>(from) {}

  inline virtual double ProposeMove(double tuning) { return T::ProposeMove(tuning); }

  inline void ClampAt(const T& inval) { T::operator=(inval); Clamp(); }

  inline virtual void	Corrupt(bool bk) { Var<T>::localcorrupt(bk); Rnode::Corrupt(bk); }
  inline virtual void  Restore()	{ Var<T>::localrestore(); Rnode::Restore(); }
  inline void RestoreBackup(){ Var<T>::localrestore(); /*value_updated = true;*/ }
};


template<class T> class Dvar : public Var<T>, public Dnode {
public:
  Dvar() {}

  Dvar(const T& from) : Var<T>(from) {}

  inline void	localCorrupt(bool bk)	{ Var<T>::localcorrupt(bk); Dnode::localCorrupt(bk); }
  inline void localRestore() { Var<T>::localrestore(); Dnode::localRestore(); }

  // returns the current value (VL: won't move into .cpp for just one function)
  const T& val() {
    if (! isUpdated())	{
      if (! DAGnode::initmode)	{
        std::cerr << "error : corrupting Dvar\n";
      }
      localUpdate();
    }
    return (*this);
  }
};


template <class T> class Const : public Dvar<T>	{
public:
  Const(const T& from) : Dvar<T>(from) {
    std::ostringstream tmp;
    tmp << "const = " << from ;
    Dvar<T>::SetName(tmp.str());
  }

  void specialUpdate() {}
};


template<> class Rvar<void> : public virtual Rnode {};
template<> class Dvar<void> : public virtual Dnode {};


#endif
