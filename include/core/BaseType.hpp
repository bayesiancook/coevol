#ifndef BASETYPE_H
#define BASETYPE_H

#include <iostream>
#include <vector>
#include "Random.hpp"

class DAGnode; // forward declaration

/// An interface for all base types that will be used to make Random variables
// e.g. Real, PosReal, Int, Profile
// base types:
// - know their domain of definition (which they should check in Check())
// - propose default kernels (in ProposeMove()) for Metropolis Hastings resampling
class BaseType {
public:
  virtual ~BaseType() {}

  // default kernel for Metropolis Hastings updates
  // smaller tuning: smaller moves
  // returns the log of the Hastings ratio
  virtual double ProposeMove(double tuning) = 0;
};


class Additive {
public:
  virtual ~Additive() {}

  // returns the number of components that have been added d
  virtual int ScalarAddition(double d) = 0;

  virtual void Register(DAGnode*) ;
};


class Multiplicative {
public:
  virtual ~Multiplicative() {}

  // returns the number of components that have been multiplied
  virtual int ScalarMultiplication(double d) = 0;

  virtual void Register(DAGnode*) ;
};


/// A wrap-up class for real numbers
// implements a simple random translational (additive) move
class Real : public BaseType , public Additive {
public:
  Real(double d=0) : value(d) {}
  Real(const Real& from) : value(from.value) {}
  ~Real() {}

  /* (VL) Operators */
  inline Real& operator=(const Real& from) { value = from.value; return *this; }
  inline Real& operator=(double from) { value = from; return *this; }
  inline Real& operator+=(const Real& from) { value += from.value; return *this; }
  inline Real& operator/=(double from) { value /= from; return *this; }
  inline operator double() { return value; }
  inline operator double() const { return value; }

  /* (VL) Stream operator friend function */
  friend std::istream& operator>>(std::istream& is, Real& r) ;

  /* (VL) Overriding virtual methods from parents */
  int ScalarAddition(double d) override;
  double ProposeMove(double tuning) override ;

  /* (VL) what exactly are check functions? */
  int check() ;

protected:
  double value;
};


class UnitReal : public BaseType {
public:
  UnitReal(double d=0) : value(d) {}
  UnitReal(const UnitReal& from) : value(from.value) {}
  ~UnitReal() {}

  /* (VL) Operators */
  inline UnitReal& operator=(const UnitReal& from) { value = from.value; return *this; }
  inline UnitReal& operator=(double from) { value = from; return *this; }
  inline UnitReal& operator+=(const UnitReal from) { value += from.value; return *this; }
  inline operator double() { return value; }
  inline operator double() const { return value; }

  /* (VL) Stream operator friend function */
  friend std::istream& operator>>(std::istream& is, UnitReal& r) ;

  /* (VL) Overriding ProposeMove from BaseType class */
  double ProposeMove(double tuning) override ;

  int check() ;


protected:
  double value;
};


/// A wrap-up class for positive real numbers
// implements a simple random multiplicative move
class PosReal : public BaseType , public Multiplicative {
public:
  PosReal(double d=0) : value(d) {}
  PosReal(const PosReal& from) : value(from.value) {}  virtual ~PosReal() {}

  /* (VL) Operators */
  inline PosReal& operator=(const PosReal& from) { value = from.value; return *this; }
  inline PosReal& operator=(const double& from) { value = from; return *this; }
  inline PosReal& operator+=(const PosReal from) { value += from.value; return *this; }
  inline PosReal& operator/=(double from) { value /= from; return *this; }
  operator double() { return value; }
  operator double() const { return value; }
  operator Real() { return Real(value); }
  operator Real() const { return Real(value); }

  /* (VL) Stream operator friend function */
  friend std::istream& operator>>(std::istream& is, PosReal& r) ;

  /* (VL) Overriding virtual methods from base classes */
  int ScalarMultiplication(double d) override ;
  double ProposeMove(double tuning) override ;

  int check() ;


protected:
  double value;
};


/// A wrap-up class for integers
// discretized additive move
class Int : public BaseType {
public:
  Int(int d=0) : value(d) {}
  Int(const Int& from) : value(from.value) {}
  ~Int() {}

  /* (VL) Operators */
  inline Int& operator=(const Int& from) { value = from.value; return *this; }
  inline Int& operator=(const int& from) { value = from; return *this; }
  inline operator int() { return value; }
  inline operator int() const { return value; }
  inline operator Real() { return Real(double(value)); }
  inline operator Real() const { return Real(double(value)); }

  /* (VL) Stream operator friend function */
  friend std::istream& operator>>(std::istream& is, Int& r) ;

  /* (VL) Overriding virtual methods from base classes */
  virtual double ProposeMove(double) override ;

  inline int check() { return 1; }

protected:
  int value;
};


/// A probability profile
class Profile : public BaseType {
protected:
  static const double MIN;
  std::vector<double> profile;

public:
  Profile() : profile(0) {}
  Profile(int indim, double* v=0) ; //TODO should be deleted when vectorification is done
  Profile(const Profile& from) : profile(from.profile) {}
  virtual ~Profile() {}

  Profile& operator=(const Profile& from) ;
  inline double& operator[](int i) { return profile[i]; }
  inline const double& operator[](int i) const { return profile[i]; }

  // Getters FIXME (these and the setters below should probably be inlined)
  inline const double* GetArray() const { return &profile[0]; }
  inline double* GetArray() { return &profile[0]; }
  inline int GetDim() const { return profile.size(); }
  double GetEntropy() const ;

  // Setters
  inline void setAtZero() { for (auto& i : profile) i = 0; }
  void setuniform() ;
  void setarray(double* in) ;

  inline void scalarMultiplication(double d) { for (auto& i : profile) i *= d; }
  inline void add(const Profile& in) { for (unsigned int i=0; i<profile.size(); i++) profile[i] += in[i]; }

  inline int check() { return 1; }

  double ProposeMove(double tuning, int dim);
  double ProposeMove(double tuning) override { return ProposeMove(tuning, profile.size()); }

  friend std::ostream& operator<<(std::ostream& os, const Profile& r) ;
  friend std::istream& operator>>(std::istream& is, Profile& r) ;
};


class RealVector : public BaseType, public Additive {
protected:
  std::vector<double> vec;

public:
  RealVector(int indim = 0) : vec(indim) {}
  RealVector(const RealVector& from) : vec(from.vec) {}
  RealVector(const double* from, int indim) ;
  virtual ~RealVector() {}

  RealVector& operator=(const RealVector& from) ;

  inline double& operator[](int i) { return vec[i]; }
  inline const double& operator[](int i) const { return vec[i]; }

  inline double* GetArray() { return &vec[0]; }
  inline int GetDim() const { return vec.size(); }
  double GetMean() const ;
  double GetVar() const ;

  inline int check() { return 1; }

  int ScalarAddition(double d) ;
  void ScalarMultiplication(double d) ;

  void add(const RealVector& in) ;
  void add(const double* in, double f = 1) ;

  double ProposeMove(double tuning, int n) ;
  inline double ProposeMove(double tuning) { return ProposeMove(tuning, GetDim()); }

  inline void setAtZero() { for (auto& i : vec) i = 0; }

  friend std::ostream& operator<<(std::ostream& os, const RealVector& r) ;
  friend std::istream& operator>>(std::istream& is, RealVector& r) ;
};


class PosRealVector : public RealVector, public Multiplicative {
public:
  PosRealVector() : RealVector() {}
  PosRealVector(int indim) : RealVector(indim) {}
  PosRealVector(const PosRealVector& from) : RealVector(from) {};
  PosRealVector(const double* from, int indim) ;
  virtual ~PosRealVector() {}

  PosRealVector& operator=(const PosRealVector& from) ;

  inline void SetAtOne() { for (auto& i : vec) i = 1; }

  double GetMean() const ;
  double GetVar() const ;
  double GetEntropy() const ;

  double ProposeMove(double tuning, int n) ;
  inline double ProposeMove(double tuning) { return ProposeMove(tuning, GetDim()); }

  inline int ScalarMultiplication(double d) { for (auto& i : vec) i *= d; return GetDim(); }
};


class IntVector : public BaseType {
protected:
  int dim;
  int* vec;

public:
  IntVector() : dim(0), vec(0) {}
  IntVector(int indim) { dim = indim; vec = new int[dim]; }
  IntVector(const IntVector& from) ;
  IntVector(const int* from, int indim) ;
  virtual   ~IntVector() { delete[] vec; }

  IntVector& operator=(const IntVector& from) ;
  IntVector& operator=(const int* from) ;

  int& operator[](int i) { return vec[i]; }
  int& operator[](int i) const { return vec[i]; }

  inline const int* GetArray() const { return vec; }
  inline int GetDim() { return dim; }
  double GetMean() const ;
  double GetVar() const ;

  inline int check() { return 1; }

  int ProposeMove(double tuning, int n) ;
  inline double ProposeMove(double tuning) { return ProposeMove(tuning,dim); }

  friend std::ostream& operator<<(std::ostream& os, const IntVector& r) ;
  friend std::istream& operator>>(std::istream& is, IntVector& r) ;
};


#endif // BASETYPE_H
