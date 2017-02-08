#ifndef BASETYPE_H
#define BASETYPE_H

#include <iostream>
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
  virtual ~Real() {}

  inline Real& operator=(const Real& from) { value = from.value; return *this; }
  inline Real& operator=(double from) { value = from; return *this; }
  inline Real& operator+=(const Real& from) { value += from.value; return *this; }
  inline Real& operator/=(double from) { value /= from; return *this; }
  inline operator double() { return value; }
  inline operator double() const { return value; }

  friend std::istream& operator>>(std::istream& is, Real& r) ;

  int ScalarAddition(double d) ;

  virtual double ProposeMove(double tuning) ;

  int check() ;

protected:
  double value;

};


class UnitReal : public BaseType {
public:
  UnitReal(double d=0) : value(d) {}
  UnitReal(const UnitReal& from) : value(from.value) {}
  virtual   ~UnitReal() {}

  inline UnitReal& operator=(const UnitReal& from) { value = from.value; return *this; }
  inline UnitReal& operator=(double from) { value = from; return *this; }
  inline UnitReal& operator+=(const UnitReal from) { value += from.value; return *this; }
  operator double() { return value; }
  operator double() const { return value; }

  virtual double ProposeMove(double tuning) ;

  int check() ;

  friend std::istream& operator>>(std::istream& is, UnitReal& r) ;

protected:
  double value;

};


/// A wrap-up class for positive real numbers
// implements a simple random multiplicative move
class PosReal : public BaseType , public Multiplicative {
public:
  PosReal(double d=0) : value(d) {}
  PosReal(const PosReal& from) : value(from.value) {}  virtual ~PosReal() {}

  inline PosReal& operator=(const PosReal& from) { value = from.value; return *this; }
  inline PosReal& operator=(const double& from) { value = from; return *this; }
  inline PosReal& operator+=(const PosReal from) { value += from.value; return *this; }
  inline PosReal& operator/=(double from) { value /= from; return *this; }
  operator double() { return value; }
  operator double() const { return value; }
  operator Real() { return Real(value); }
  operator Real() const { return Real(value); }

  int ScalarMultiplication(double d) ;

  virtual double ProposeMove(double tuning) ;

  int check() ;

  friend std::istream& operator>>(std::istream& is, PosReal& r) ;

protected:
  double value;

};


/// A wrap-up class for integers
// discretized additive move
class Int : public BaseType {
public:
  Int(int d=0) : value(d) {}
  Int(const Int& from) : value(from.value) {}
  virtual ~Int() {}

  inline Int& operator=(const Int& from) { value = from.value; return *this; }
  inline Int& operator=(const int& from) { value = from; return *this; }
  inline operator int() { return value; }
  inline operator int() const { return value; }
  inline operator Real() { return Real(double(value)); }
  inline operator Real() const { return Real(double(value)); }

  virtual double ProposeMove(double) ;

  inline int check() { return 1; }

  friend std::istream& operator>>(std::istream& is, Int& r) ;

protected:
  int value;

};


/// A probability profile
class Profile : public BaseType {
protected:
  static const double MIN;
  int dim;
  double* profile;

public:
  Profile() : dim(0) , profile(0) {}
  Profile(int indim, double* v=0) ;
  Profile(const Profile& from) ;
  virtual ~Profile() ;

  Profile& operator=(const Profile& from) ;
  double& operator[](int i) ;
  double& operator[](int i) const  ;

  // Getters FIXME (these and the setters below should probably be inlined)
  const double* GetArray() const ;
  double* GetArray() ;
  int GetDim() const ;
  double GetEntropy() const ;

  // Setters
  void setAtZero() ;
  void setuniform() ;
  void setarray(double* in) ;

  void scalarMultiplication(double d) ;
  void add(const Profile& in) ;

  int check() ;

  double ProposeMove(double tuning, int dim);
  double ProposeMove(double tuning) ;

  friend std::ostream& operator<<(std::ostream& os, const Profile& r) ;
  friend std::istream& operator>>(std::istream& is, Profile& r) ;

};

class RealVector : public BaseType, public Additive {
protected:
  int dim;
  double* vec;

public:
  RealVector() : dim(0), vec(0) {}

  RealVector(int indim) {
    dim = indim;
    vec = new double[dim];
  }

  RealVector(const RealVector& from) {
    dim = from.dim;
    vec = new double[dim];
    for (int i=0; i<dim; i++) {
      vec[i] = from.vec[i];
    }
  }

  RealVector(const double* from, int indim) {
    dim = indim;
    vec = new double[dim];
    for (int i=0; i<dim; i++) {
      vec[i] = from[i];
    }
  }

  virtual   ~RealVector() {
    delete[] vec;
  }

  RealVector& operator=(const RealVector& from) {
    if (!dim) {
      dim = from.dim;
      vec = new double[dim];
    }
    if (dim != from.dim) {
      std::cerr << "error : non matching dimenstion for vectors\n";
      std::cerr << dim << '\t' << from.dim << '\n';
      exit(1);
      delete[] vec;
      dim = from.dim;
      vec = new double[dim];
    }
    for (int i=0; i<dim; i++) {
      vec[i] = from.vec[i];
    }
    return *this;
  }

  double* GetArray() const { return vec; }

  double& operator[](int i) {
    return vec[i];
  }

  double& operator[](int i) const {
    return vec[i];
  }

  int GetDim() { return dim; }
  int check() { return 1; }

  double GetMean() const {
    double total = 0;
    for (int i=0; i<dim; i++) {
      total += vec[i];
    }
    return total / dim;
  }

  double GetVar() const {
    double mean = 0;
    double var = 0;
    for (int i=0; i<dim; i++) {
      var += vec[i] * vec[i];
      mean += vec[i];
    }
    mean /= dim;
    var /= dim;
    var -= mean * mean;
    return var;
  }

  int ScalarAddition(double d) {
    for (int i=0; i<dim; i++) {
      vec[i] += d;
    }
    return dim;
  }

  void ScalarMultiplication(double d) {
    for (int i=0; i<dim; i++) {
      vec[i] *= d;
    }
  }

  void add(const RealVector& in) {
    for (int i=0; i<dim; i++) {
      vec[i] += in[i];
    }
  }

  void add(const double* in, double f = 1) {
    for (int i=0; i<dim; i++) {
      vec[i] += f * in[i];
    }
  }




  double ProposeMove(double tuning, int n) {
    if ((n<=0) || (n > dim)) {
      n = dim;
    }
    int* indices = new int[n];
    Random::DrawFromUrn(indices,n,dim);
    for (int i=0; i<n; i++) {
      vec[indices[i]] += tuning * (Random::Uniform() - 0.5);
    }
    delete[] indices;
    return 0;
  }

  double ProposeMove(double tuning) {
    return ProposeMove(tuning,dim);
  }

  void setAtZero() {
    for (int i=0; i<dim; i++) {
      vec[i] = 0;
    }
  }

  friend std::ostream& operator<<(std::ostream& os, const RealVector& r) {
    os << r.dim;
    for (int i=0; i<r.dim; i++) {
      os << '\t' << r.vec[i];
    }
    return os;
  }

  friend std::istream& operator>>(std::istream& is, RealVector& r) {
    int indim;
    is >> indim;
    if (r.dim != indim) {
      r.dim = indim;
      delete[] r.vec;
      r.vec = new double[r.dim];
    }
    for (int i=0; i<r.dim; i++) {
      is >> r.vec[i];
    }
    return is;
  }

};

class PosRealVector : public RealVector, public Multiplicative {
public:
  PosRealVector() : RealVector() {}

  PosRealVector(int indim) {
    dim = indim;
    vec = new double[dim];
  }

  PosRealVector(const PosRealVector& from): RealVector() {
    dim = from.dim;
    vec = new double[dim];
    for (int i=0; i<dim; i++) {
      vec[i] = from.vec[i];
    }
  }

  PosRealVector(const double* from, int indim) {
    dim = indim;
    vec = new double[dim];
    for (int i=0; i<dim; i++) {
      vec[i] = from[i];
    }
  }

  virtual   ~PosRealVector() {}

  PosRealVector& operator=(const PosRealVector& from) {
    if (!dim) {
      dim = from.dim;
      vec = new double[dim];
    }
    if (dim != from.dim) {
      std::cerr << "error : non matching dimenstion for pos vectors\n";
      std::cerr << dim << '\t' << from.dim << '\n';
      exit(1);
      delete[] vec;
      dim = from.dim;
      vec = new double[dim];
    }
    for (int i=0; i<dim; i++) {
      vec[i] = from.vec[i];
    }
    return *this;
  }

  double GetMean() const {
    double total = 0;
    for (int i=0; i<dim; i++) {
      total += vec[i];
    }
    return total / dim;
  }

  void SetAtOne() {
    for (int i=0; i<dim; i++) {
      vec[i] = 1;
    }
  }

  double GetVar() const {
    double mean = 0;
    double var = 0;
    for (int i=0; i<dim; i++) {
      var += vec[i] * vec[i];
      mean += vec[i];
    }
    mean /= dim;
    var /= dim;
    var -= mean * mean;
    return var;
  }

  double GetEntropy() const ;


  double ProposeMove(double tuning, int n) ;

  double ProposeMove(double tuning) {
    return ProposeMove(tuning,dim);
  }

  int ScalarMultiplication(double d) {
    for (int i=0; i<dim; i++) {
      vec[i] *= d;
    }
    return dim;
  }

};


class IntVector : public BaseType {
protected:
  int dim;
  int* vec;

public:
  IntVector() : dim(0), vec(0) {}

  IntVector(int indim) {
    dim = indim;
    vec = new int[dim];
  }

  IntVector(const IntVector& from) {
    dim = from.dim;
    vec = new int[dim];
    for (int i=0; i<dim; i++) {
      vec[i] = from.vec[i];
    }
  }

  IntVector(const int* from, int indim) {
    dim = indim;
    vec = new int[dim];
    for (int i=0; i<dim; i++) {
      vec[i] = from[i];
    }
  }

  virtual   ~IntVector() {
    delete[] vec;
  }

  IntVector& operator=(const IntVector& from) {
    if (!dim) {
      dim = from.dim;
      vec = new int[dim];
    }
    if (dim != from.dim) {
      std::cerr << "error : non matching dimenstion for vectors\n";
      std::cerr << dim << '\t' << from.dim << '\n';
      exit(1);
      delete[] vec;
      dim = from.dim;
      vec = new int[dim];
    }
    for (int i=0; i<dim; i++) {
      vec[i] = from.vec[i];
    }
    return *this;
  }

  IntVector& operator=(const int* from) {
    if (!dim) {
      std::cerr << "error in IntVector::operator=(const int*)\n";
      exit(1);
    }
    for (int i=0; i<dim; i++) {
      vec[i] = from[i];
    }
    return *this;
  }

  const int* GetArray() const { return vec; }

  int& operator[](int i) {
    return vec[i];
  }

  int& operator[](int i) const {
    return vec[i];
  }

  int GetDim() { return dim; }
  int check() { return 1; }

  double    GetMean() const {
    double total = 0;
    for (int i=0; i<dim; i++) {
      total += vec[i];
    }
    return total / dim;
  }

  double GetVar() const {
    double mean = 0;
    double var = 0;
    for (int i=0; i<dim; i++) {
      var += vec[i] * vec[i];
      mean += vec[i];
    }
    mean /= dim;
    var /= dim;
    var -= mean * mean;
    return var;
  }

  int ProposeMove(double tuning, int n) {
    if ((n<=0) || (n > dim)) {
      n = dim;
    }
    int* indices = new int[n];
    Random::DrawFromUrn(indices,n,dim);
    for (int i=0; i<n; i++) {
      vec[indices[i]] += (int) (tuning * (Random::Uniform() - 0.5));
    }
    delete[] indices;
    return 0;
  }

  double ProposeMove(double tuning) {
    return ProposeMove(tuning,dim);
  }

  friend std::ostream& operator<<(std::ostream& os, const IntVector& r) {
    os << r.dim;
    for (int i=0; i<r.dim; i++) {
      os << '\t' << r.vec[i];
    }
    return os;
  }

  friend std::istream& operator>>(std::istream& is, IntVector& r) {
    int indim;
    is >> indim;
    if (r.dim != indim) {
      r.dim = indim;
      delete[] r.vec;
      r.vec = new int[r.dim];
    }
    for (int i=0; i<r.dim; i++) {
      is >> r.vec[i];
    }
    return is;
  }

};

#endif // BASETYPE_H
