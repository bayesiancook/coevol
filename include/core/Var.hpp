#ifndef VAR_H
#define VAR_H

#include <iostream>  // FIXME could be removed if implementation is moved into a file
#include <sstream>
#include "DAGnode.hpp"


template <class T>
class Var : public virtual DAGnode, public T {
  public:
    Var() = default;
    Var(const T& from) : T(from), bkvalue(from) {}
    ~Var() override = default;

    inline Var& operator=(const T& from) { return T::operator=(from); }

    inline void localcorrupt(bool bk) {
        if (bk) {
            this->bkvalue = *this;
        }
    }
    inline void localrestore() { T::operator=(bkvalue); }

    // returns the current value
    inline T& val() { return *this; }

    // sets the new value
    inline void setval(const T& inval) { T::operator=(inval); }

    inline void Register(DAGnode* in) override { DAGnode::Register(in); }

  protected:
    T bkvalue;
};


template <class T>
class Rvar : public Var<T>, public Rnode {
  public:
    Rvar() = default;
    Rvar(const T& from) : Var<T>(from) {}

    inline double ProposeMove(double tuning) override { return T::ProposeMove(tuning); }

    inline void ClampAt(const T& inval) {
        T::operator=(inval);
        Clamp();
    }

    inline void Corrupt(bool bk) override {
        Var<T>::localcorrupt(bk);
        Rnode::Corrupt(bk);
    }
    inline void Restore() override {
        Var<T>::localrestore();
        Rnode::Restore();
    }
    inline void RestoreBackup() final { Var<T>::localrestore(); /*value_updated = true;*/ }
};


template <class T>
class Dvar : public Var<T>, public Dnode {
  public:
    Dvar() = default;

    Dvar(const T& from) : Var<T>(from) {}

    inline void localCorrupt(bool bk) override {
        Var<T>::localcorrupt(bk);
        Dnode::localCorrupt(bk);
    }
    inline void localRestore() override {
        Var<T>::localrestore();
        Dnode::localRestore();
    }

    // returns the current value (VL: won't move into .cpp for just one function)
    const T& val() {
        if (!isUpdated()) {
            if (!DAGnode::initmode) {
                std::cerr << "error : corrupting Dvar\n";
            }
            localUpdate();
        }
        return (*this);
    }
};


template <class T>
class Const : public Dvar<T> {
  public:
    Const(const T& from) : Dvar<T>(from) {
        std::ostringstream tmp;
        tmp << "const = " << from;
        Dvar<T>::SetName(tmp.str());
    }

    void specialUpdate() {}
};


template <>
class Rvar<void> : public virtual Rnode {};
template <>
class Dvar<void> : public virtual Dnode {};


#endif
