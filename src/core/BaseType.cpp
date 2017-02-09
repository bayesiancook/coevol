#include "core/BaseType.hpp"
#include <cmath>
using namespace std;


//-------------------------------------------------------------------------
//-------------------------------------------------------------------------
//	* Additive
//-------------------------------------------------------------------------
//-------------------------------------------------------------------------
void Additive::Register(DAGnode*) {std::cerr << "error in Additive::Register\n"; throw;}


//-------------------------------------------------------------------------
//-------------------------------------------------------------------------
//	* Multiplicative
//-------------------------------------------------------------------------
//-------------------------------------------------------------------------
void Multiplicative::Register(DAGnode*)  {std::cerr << "error in Multiplicative::Register\n"; throw;}


//-------------------------------------------------------------------------
//-------------------------------------------------------------------------
//	* Real
//-------------------------------------------------------------------------
//-------------------------------------------------------------------------
int Real::ScalarAddition(double d) {
  value += d;
  return 1;
}

double Real::ProposeMove(double tuning) {
  // simple additive move
  double m = tuning*(Random::Uniform() - 0.5);
  value += m;
  return 0;
}

int Real::check() {return 1;}

std::istream& operator>>(std::istream& is, Real& r) {
  is >> r.value;
  return is;
}


//-------------------------------------------------------------------------
//-------------------------------------------------------------------------
//	* UnitReal
//-------------------------------------------------------------------------
//-------------------------------------------------------------------------
double UnitReal::ProposeMove(double tuning) {
  // simple additive move
  double m = tuning*(Random::Uniform() - 0.5);
  value += m;
  while ((value<0) || (value>1)) {
    if (value<0) {
      value = -value;
    }
    if (value>1) {
      value = 2 - value;
    }
  }
  return 0.0;
}

int UnitReal::check() {return 1;}

std::istream& operator>>(std::istream& is, UnitReal& r) {
  is >> r.value;
  return is;
}


//-------------------------------------------------------------------------
//-------------------------------------------------------------------------
//	* PosReal
//-------------------------------------------------------------------------
//-------------------------------------------------------------------------
int PosReal::ScalarMultiplication(double d) {
  value *= d;
  return 1;
}

double PosReal::ProposeMove(double tuning) {
  // simple multiplicative move
  double m = tuning*(Random::Uniform() - 0.5);
  value *= exp(m);
  return m;
}

int PosReal::check() {
  if (value<=0) {
    std::cerr << "error : positive double is not positive : " << value << '\n';
    return 0;
  }
  return 1;
}

std::istream& operator>>(std::istream& is, PosReal& r) {
  is >> r.value;
  return is;
}


//-------------------------------------------------------------------------
//-------------------------------------------------------------------------
//	* Int
//-------------------------------------------------------------------------
//-------------------------------------------------------------------------
double Int::ProposeMove(double) {
  if (Random::Uniform() < 0.5) {
    value++;
  }
  else {
    value --;
  }
  return 0;
}

std::istream& operator>>(std::istream& is, Int& r) {
  is >> r.value;
  return is;
}


//-------------------------------------------------------------------------
//-------------------------------------------------------------------------
//	* Profile
//-------------------------------------------------------------------------
//-------------------------------------------------------------------------
const double Profile::MIN = 1e-20;

Profile::Profile(int indim, double* v) {
  profile.assign(indim, 0);
  if (v) {
    double total = 0;
    for (int k=0; k<indim; k++) { // (VL) copy v into profile + checking everything is positive
      if (!(v[k]>0)) { // (VL) + computing total to normalize afterwards
        std::cerr << "error : profiles should be strictly positive\n";
        exit(1);
      }
      profile[k] = v[k];
      total += profile[k];
    }
    for (auto k : profile)
      profile[k] /= total; // (VL) Normalizing happens here
  }
}

double	Profile::ProposeMove(double tuning, int n)	{ // n==0dirichlet resampling, otherwise, vase communiquants
  double ret = 0;
  int dim = profile.size();
  if (!n)	{ // dirichlet
    vector<double> oldprofile = profile;
    double total = 0;
    for (int i=0; i<dim; i++)	{
      profile[i] = Random::sGamma(tuning*oldprofile[i]);
      if (profile[i] == 0)	{
        std::cerr << "error in dirichlet resampling : 0 \n";
        exit(1);
      }
      total += profile[i];
    }
    double logHastings = 0;
    for (int i=0; i<dim; i++)	{
      profile[i] /= total;
      logHastings += - Random::logGamma(tuning*oldprofile[i]) + Random::logGamma(tuning*profile[i])
        -  (tuning*profile[i] -1.0) * log(oldprofile[i]) + (tuning * oldprofile[i] -1.0) * log(profile[i]);
    }
    return logHastings;
  }
  else	{
    if (2*n > dim)	{
      n = dim / 2;
    }
    int* indices = new int[2*n];
    Random::DrawFromUrn(indices,2*n,dim);
    for (int i=0; i<n; i++)	{
      int i1 = indices[2*i];
      int i2 = indices[2*i+1];
      double tot = profile[i1] + profile[i2];
      double x = profile[i1];

      // double h = tuning * (Random::Uniform() - 0.5);
      double h = tot * tuning * (Random::Uniform() - 0.5);
      /*
        int c = (int) (h / (2 * tot));
        h -= c*2*tot;
      */
      x += h;
      while ((x<0) || (x>tot))	{
        if (x<0)	{
          x = -x;
        }
        if (x>tot)	{
          x = 2*tot - x;
        }
      }
      profile[i1] = x;
      profile[i2] = tot - x;
    }
    delete[] indices;
  }
  return ret;
}

Profile& Profile::operator=(const Profile& from) {
  profile = from.profile;
  return *this;
}

void Profile::setuniform() {
  double value = 1.0 /profile.size();
  for (auto& i: profile)
    i = value;
}

void Profile::setarray(double* in) {
  for (unsigned int i=0; i<profile.size(); i++) {
    profile[i] = in[i];
  }
}

/* (VL) TODO inline all thos little methods TODO */
const double* Profile::GetArray() const { return &profile[0]; }

double* Profile::GetArray() { return &profile[0]; }

double& Profile::operator[](int i) { return profile[i]; }

const double& Profile::operator[](int i) const { return profile[i]; }

int Profile::GetDim() const { return profile.size(); }

void Profile::setAtZero() { for (auto& i : profile) i = 0; }

void Profile::scalarMultiplication(double d) { for (auto& i : profile) i *= d; }

void Profile::add(const Profile& in) {
  for (unsigned int i=0; i<profile.size(); i++)
    profile[i] += in[i];
}

int Profile::check() { return 1; }

double Profile::GetEntropy() const {
  double total = 0;
  for (auto i : profile)
    total += (i>1e-8) ? -i*log(i) : 0;
  return total;
}

double Profile::ProposeMove(double tuning) {
  return ProposeMove(tuning, profile.size());
}

std::ostream& operator<<(std::ostream& os, const Profile& r) {
  os << r.GetDim();
  for (int i=0; i<r.GetDim(); i++) {
    os << '\t' << r.profile[i];
  }
  return os;
}

std::istream& operator>>(std::istream& is, Profile& r) {
  int indim;
  is >> indim;
  for (auto& i : r.profile)
    is >> i;
  return is;
}


//-------------------------------------------------------------------------
//-------------------------------------------------------------------------
//	* RealVector
//-------------------------------------------------------------------------
//-------------------------------------------------------------------------
RealVector::RealVector(const RealVector& from) {
  dim = from.dim;
  vec = new double[dim];
  for (int i=0; i<dim; i++) {
    vec[i] = from.vec[i];
  }
}

RealVector::RealVector(const double* from, int indim) {
  dim = indim;
  vec = new double[dim];
  for (int i=0; i<dim; i++) {
    vec[i] = from[i];
  }
}

RealVector& RealVector::operator=(const RealVector& from) {
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

double RealVector::GetMean() const {
  double total = 0;
  for (int i=0; i<dim; i++) {
    total += vec[i];
  }
  return total / dim;
}

double RealVector::GetVar() const {
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

int RealVector::ScalarAddition(double d) {
  for (int i=0; i<dim; i++) {
    vec[i] += d;
  }
  return dim;
}

void RealVector::ScalarMultiplication(double d) {
  for (int i=0; i<dim; i++) {
    vec[i] *= d;
  }
}

void RealVector::add(const RealVector& in) {
  for (int i=0; i<dim; i++) {
    vec[i] += in[i];
  }
}

void RealVector::add(const double* in, double f) {
  for (int i=0; i<dim; i++) {
    vec[i] += f * in[i];
  }
}

double RealVector::ProposeMove(double tuning, int n) {
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

std::ostream& operator<<(std::ostream& os, const RealVector& r) {
  os << r.dim;
  for (int i=0; i<r.dim; i++) {
    os << '\t' << r.vec[i];
  }
  return os;
}

std::istream& operator>>(std::istream& is, RealVector& r) {
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


//-------------------------------------------------------------------------
//-------------------------------------------------------------------------
//	* PosRealVector
//-------------------------------------------------------------------------
//-------------------------------------------------------------------------
PosRealVector::PosRealVector(const PosRealVector& from): RealVector() {
  dim = from.dim;
  vec = new double[dim];
  for (int i=0; i<dim; i++) {
    vec[i] = from.vec[i];
  }
}

PosRealVector::PosRealVector(const double* from, int indim) {
  dim = indim;
  vec = new double[dim];
  for (int i=0; i<dim; i++) {
    vec[i] = from[i];
  }
}

PosRealVector& PosRealVector::operator=(const PosRealVector& from) {
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

double PosRealVector::GetMean() const {
  double total = 0;
  for (int i=0; i<dim; i++) {
    total += vec[i];
  }
  return total / dim;
}

double PosRealVector::GetVar() const {
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

double PosRealVector::ProposeMove(double tuning, int n) {
  if ((n<=0) || (n > dim)) {
    n = dim;
  }
  int* indices = new int[n];
  Random::DrawFromUrn(indices,n,dim);
  double ret = 0;
  for (int i=0; i<n; i++) {
    double m = tuning * (Random::Uniform() - 0.5);
    vec[indices[i]] *= exp(m);
    ret += m;
  }
  delete[] indices;
  return ret;
}

double PosRealVector::GetEntropy() const {
  double total = 0;
  for (int i=0; i<dim; i++) {
    total += vec[i];
  }
  double ent = 0;
  for (int i=0; i<dim; i++) {
    double tmp = vec[i]/total;
    ent += (tmp>1e-8) ? -tmp*log(tmp) : 0;
  }
  return ent;
}


//-------------------------------------------------------------------------
//-------------------------------------------------------------------------
//	* IntVector
//-------------------------------------------------------------------------
//-------------------------------------------------------------------------
IntVector::IntVector(const IntVector& from) {
  dim = from.dim;
  vec = new int[dim];
  for (int i=0; i<dim; i++) {
    vec[i] = from.vec[i];
  }
}

IntVector::IntVector(const int* from, int indim) {
  dim = indim;
  vec = new int[dim];
  for (int i=0; i<dim; i++) {
    vec[i] = from[i];
  }
}

IntVector& IntVector::operator=(const IntVector& from) {
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

IntVector& IntVector::operator=(const int* from) {
  if (!dim) {
    std::cerr << "error in IntVector::operator=(const int*)\n";
    exit(1);
  }
  for (int i=0; i<dim; i++) {
    vec[i] = from[i];
  }
  return *this;
}

double IntVector::GetMean() const {
  double total = 0;
  for (int i=0; i<dim; i++) {
    total += vec[i];
  }
  return total / dim;
}

double IntVector::GetVar() const {
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

int IntVector::ProposeMove(double tuning, int n) {
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

std::ostream& operator<<(std::ostream& os, const IntVector& r) {
  os << r.dim;
  for (int i=0; i<r.dim; i++) {
    os << '\t' << r.vec[i];
  }
  return os;
}

std::istream& operator>>(std::istream& is, IntVector& r) {
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
