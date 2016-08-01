
#ifndef BASETYPE_H
#define BASETYPE_H


#include "Random.h"

class DAGnode;

/// An interface for all base types that will be used to make Random variables
// e.g. Real, PosReal, Int, Profile
// base types:
// - know their domain of definition (which they should check in Check())
// - propose default kernels (in ProposeMove()) for Metropolis Hastings resampling

class BaseType {

	public:

	virtual ~BaseType() {};

	// default kernel for Metropolis Hastings updates
	// smaller tuning: smaller moves
	// returns the log of the Hastings ratio
	virtual double ProposeMove(double tuning) = 0;


};

class Additive	{

	public:

	virtual ~Additive() {}

	// returns the number of components that have been added d
	virtual int ScalarAddition(double d) = 0;

	virtual void Register(DAGnode* in)  {cerr << "error in Additive::Register\n"; throw;}
};

class Multiplicative	{

	public:

	virtual ~Multiplicative() {}

	// returns the number of components that have been multiplied
	virtual int ScalarMultiplication(double d) = 0;

	virtual void Register(DAGnode* in)  {cerr << "error in Multiplicative::Register\n"; throw;}
};

/// A wrap-up class for real numbers
// implements a simple random translational (additive) move

class Real : public BaseType , public Additive {


	public:
			Real(double d=0) : value(d) {}
			Real(const Real& from) : value(from.value) {}

	virtual 	~Real() {}

	Real&		operator=(const Real& from)	{
				value = from.value;
				return *this;
			}

	Real&		operator=(double from)	{
				value = from;
				return *this;
			}

	operator 	double() {return value;}
	operator 	double() const {return value;}

	Real&		operator+=(const Real& from)	{
				value += from.value;
				return *this;
			}

	Real&		operator/=(double from)	{
				value /= from;
				return *this;
			}

	int		ScalarAddition(double d)	{
				value += d;
				return 1;
			}

	virtual double	ProposeMove(double tuning)	{
				// simple additive move
				double m = tuning*(Random::Uniform() - 0.5);
				value += m;
				return 0;
			}

	int		Check() {return 1;}

	friend istream& operator>>(istream& is, Real& r)  {
		is >> r.value;
		return is;
	}

	protected:
	double value;
};

class UnitReal : public BaseType {


	public:
			UnitReal(double d=0) : value(d) {}
			UnitReal(const UnitReal& from) : value(from.value) {}

	virtual 	~UnitReal() {}

	UnitReal&		operator=(const UnitReal& from)	{
				value = from.value;
				return *this;
			}

	UnitReal&		operator=(double from)	{
				value = from;
				return *this;
			}

	operator 	double() {return value;}
	operator 	double() const {return value;}

	UnitReal&		operator+=(const UnitReal from)	{
				value += from.value;
				return *this;
			}


	virtual double	ProposeMove(double tuning)	{
				// simple additive move
				double m = tuning*(Random::Uniform() - 0.5);
				value += m;
				while ((value<0) || (value>1))	{
					if (value<0)	{
						value = -value;
					}
					if (value>1)	{
						value = 2 - value;
					}
				}
				return 0.0;
			}

	int		Check() {return 1;}

	friend istream& operator>>(istream& is, UnitReal& r)  {
		is >> r.value;
		return is;
	}

	protected:
	double value;
};


/// A wrap-up class for positive real numbers
// implements a simple random multiplicative move

class PosReal : public BaseType	, public Multiplicative {


	public:
			PosReal(double d=0) : value(d) {}
			PosReal(const PosReal& from) : value(from.value) {}

	virtual 	~PosReal() {}

	PosReal&	operator=(const PosReal& from)	{
				value = from.value;
				return *this;
			}

	PosReal&	operator=(const double& from)	{
				value = from;
				return *this;
			}

	operator 	double() {return value;}
	operator 	double() const {return value;}

	PosReal&	operator+=(const PosReal from)	{
				value += from.value;
				return *this;
			}

	PosReal&	operator/=(double from)	{
				value /= from;
				return *this;
			}

	int 		ScalarMultiplication(double d)	{
				value *= d;
				return 1;
			}


 	operator	Real() {return Real(value);}
 	operator	Real() const {return Real(value);}

	virtual double	ProposeMove(double tuning)	{
				// simple multiplicative move
				double m = tuning*(Random::Uniform() - 0.5);
				value *= exp(m);
				return m;
			}

	int		Check()	{
				if (value<=0)	{
					cerr << "error : positive double is not positive : " << value << '\n';
					return 0;
				}
				return 1;
			}

	friend istream& operator>>(istream& is, PosReal& r)  {
		is >> r.value;
		return is;
	}

	/*
	friend ostream& operator<<(ostream& os, PosReal& r)  {
		os << r.value;
		return os;
	}
	*/

	protected:
	double value;
};


/// A wrap-up class for integers
// discretized additive move

class Int : public BaseType	{


	public:
			Int(int d=0) : value(d) {}
			Int(const Int& from) : value(from.value) {}

	virtual 	~Int() {}

	Int&		operator=(const Int& from)	{
				value = from.value;
				return *this;
			}

	Int&		operator=(const int& from)	{
				value = from;
				return *this;
			}

	operator 	int() {return value;}
	operator 	int() const {return value;}

	operator 	Real() {return Real(double(value));}
	operator 	Real() const {return Real(double(value));}

	virtual double	ProposeMove(double tuning)	{
				/*
				int m = (int) (tuning*(Random::Uniform() - 0.5));
				value += m;
				if (value < 0)	{
					value = -value;
				}
				*/
				if (Random::Uniform() < 0.5)	{
					value++;
				}
				else	{
					value --;
				}
				return 0;
			}

	int		Check() {return 1;}

	friend istream& operator>>(istream& is, Int& r)  {
		is >> r.value;
		return is;
	}

	protected:
	int value;
};

/*
class FinitePosInt : public Int	{


	public:
			FinitePosInt(int inmax) : value(0), max(inmax) {}
			FinitePosInt(const FinietPosInt& from) : value(from.value), max(from.max) {}

	virtual 	~FinitePosInt() {}

	FinitePosInt&	operator=(const FinitePosInt& from)	{
				value = from.value;
				return *this;
			}

	FinitePosInt&	operator=(const int& from)	{
				value = from;
				return *this;
			}

	operator 	int() {return value;}
	operator 	int() const {return value;}

	operator 	Real() {return Real(double(value));}
	operator 	Real() const {return Real(double(value));}

	virtual double	ProposeMove(double tuning)	{
				int m = (int) (tuning*(Random::Uniform() - 0.5));
				value += m;
				while ((value < 0) || (value > max))	{
					if (value < 0)
						value = -value;
					}
					if (value > max)	{
						value = 2*max - value;
					}
				}
				return 0;
			}

	int		Check() {return 1;}

	friend istream& operator>>(istream& is, FinitePosInt& r)  {
		is >> r.value;
		return is;
	}

	protected:
	int value;
};
*/

/// A probability profile

class Profile : public BaseType	{

	public:
	static const double MIN;

	protected:
	int		dim;
	double*		profile;

	public:

			Profile() : dim(0) , profile(0) {}

			Profile(int indim, double* v=0)	{
				dim = indim;
				profile = new double[dim];
				if (v)	{
					double total = 0;
					for (int k=0; k<dim; k++)	{
						if (!v[k]>0)	{
							cerr << "error : profiles should be strictly positive\n";
							exit(1);
						}
						profile[k] = v[k];
						total += profile[k];
					}
					for (int k=0; k<dim; k++)	{
						profile[k] /= total;
					}
				}
			}

			Profile(const Profile& from)	{
				dim = from.dim;
				profile = new double[dim];
				for (int k=0; k<dim; k++)	{
					profile[k] = from.profile[k];
				}
			}

	virtual		~Profile()	{
				delete[] profile;
			}

	Profile&	operator=(const Profile& from)	{

				if (!dim)	{
					dim = from.dim;
					profile = new double[dim];
				}
				if (dim != from.dim)	{
					cerr << "error : non matching dimenstion for profiles\n";
					cerr << dim << '\t' << from.dim << '\n';
					exit(1);
					dim = from.dim;
					delete[] profile;
					profile = new double[dim];
				}
				for (int k=0; k<dim; k++)	{
					profile[k] = from.profile[k];
				}
				return *this;
			}

	void		setuniform()	{
				for (int i=0; i<dim; i++)	{
					profile[i] = 1.0 / dim;
				}
			}

	void		setarray(double* in)	{
				for (int i=0; i<dim; i++)	{
					profile[i] = in[i];
				}
			}


	const double*		GetArray() const {return profile;}
	double*		GetArray() {return profile;}

	double&		operator[](int i)	{
				return profile[i];
			}

	double&		operator[](int i) const 	{
				return profile[i];
			}

	int		GetDim() const {return dim;}

	void SetAtZero()	{
		for (int i=0; i<dim; i++)	{
			profile[i] = 0;
		}
	}

	void		ScalarMultiplication(double d)	{
				for (int i=0; i<dim; i++)	{
					profile[i] *= d;
				}
			}

	void		Add(const Profile& in)	{
				for (int i=0; i<dim; i++)	{
					profile[i] += in[i];
				}
			}

	int		Check() {return 1;}

	double		GetEntropy()	const {
		double total = 0;
		for (int i=0; i<dim; i++)	{
			total += (profile[i]>1e-8) ? -profile[i]*log(profile[i]) : 0;
		}
		return total;
	}

	double ProposeMove(double tuning, int dim);

	double ProposeMove(double tuning)	{
		return ProposeMove(tuning,dim);
	}

	friend ostream& operator<<(ostream& os, const Profile& r)  {
		os << r.dim;
		for (int i=0; i<r.dim; i++)	{
			os << '\t' << r.profile[i];
		}
		return os;
	}

	friend istream& operator>>(istream& is, Profile& r)  {
		int indim;
		is >> indim;
		if (r.dim != indim)	{
			r.dim = indim;
			delete[] r.profile;
			r.profile = new double[r.dim];
		}
		for (int i=0; i<r.dim; i++)	{
			is >> r.profile[i];
		}
		return is;
	}
};

class RealVector : public BaseType, public Additive {

	protected:

	int dim;
	double* vec;

	public:
			RealVector() : dim(0), vec(0) {}

			RealVector(int indim)	{
				dim = indim;
				vec = new double[dim];
			}

			RealVector(const RealVector& from)	{
				dim = from.dim;
				vec = new double[dim];
				for (int i=0; i<dim; i++)	{
					vec[i] = from.vec[i];
				}
			}

			RealVector(const double* from, int indim)	{
				dim = indim;
				vec = new double[dim];
				for (int i=0; i<dim; i++)	{
					vec[i] = from[i];
				}
			}

	virtual 	~RealVector() {
				delete[] vec;
			}

	RealVector&		operator=(const RealVector& from)	{
				if (!dim)	{
					dim = from.dim;
					vec = new double[dim];
				}
				if (dim != from.dim)	{
					cerr << "error : non matching dimenstion for vectors\n";
					cerr << dim << '\t' << from.dim << '\n';
					exit(1);
					delete[] vec;
					dim = from.dim;
					vec = new double[dim];
				}
				for (int i=0; i<dim; i++)	{
					vec[i] = from.vec[i];
				}
				return *this;
			}

	double*		GetArray() const {return vec;}

	double&		operator[](int i)	{
				return vec[i];
			}

	double&		operator[](int i) const	{
				return vec[i];
			}

	int		GetDim() {return dim;}
	int		Check() {return 1;}

	double		GetMean() const	{
		double total = 0;
		for (int i=0; i<dim; i++)	{
			total += vec[i];
		}
		return total / dim;
	}

	double		GetVar() const 	{
		double mean = 0;
		double var = 0;
		for (int i=0; i<dim; i++)	{
			var += vec[i] * vec[i];
			mean += vec[i];
		}
		mean /= dim;
		var /= dim;
		var -= mean * mean;
		return var;
	}

	int		ScalarAddition(double d)	{
				for (int i=0; i<dim; i++)	{
					vec[i] += d;
				}
				return dim;
			}

	void		ScalarMultiplication(double d)	{
				for (int i=0; i<dim; i++)	{
					vec[i] *= d;
				}
			}

	void		Add(const RealVector& in)	{
				for (int i=0; i<dim; i++)	{
					vec[i] += in[i];
				}
			}

	void		Add(const double* in, double f = 1)	{
				for (int i=0; i<dim; i++)	{
					vec[i] += f * in[i];
				}
			}




	double	ProposeMove(double tuning, int n)	{
		if ((n<=0) || (n > dim))	{
			n = dim;
		}
		int* indices = new int[n];
		Random::DrawFromUrn(indices,n,dim);
		for (int i=0; i<n; i++)	{
			vec[indices[i]] += tuning * (Random::Uniform() - 0.5);
		}
		delete[] indices;
		return 0;
	}

	double ProposeMove(double tuning)	{
		return ProposeMove(tuning,dim);
	}

	void SetAtZero()	{
		for (int i=0; i<dim; i++)	{
			vec[i] = 0;
		}
	}

	friend ostream& operator<<(ostream& os, const RealVector& r)  {
		os << r.dim;
		for (int i=0; i<r.dim; i++)	{
			os << '\t' << r.vec[i];
		}
		return os;
	}

	friend istream& operator>>(istream& is, RealVector& r)  {
		int indim;
		is >> indim;
		if (r.dim != indim)	{
			r.dim = indim;
			delete[] r.vec;
			r.vec = new double[r.dim];
		}
		for (int i=0; i<r.dim; i++)	{
			is >> r.vec[i];
		}
		return is;
	}
};

class PosRealVector : public RealVector, public Multiplicative	{

	public:
			PosRealVector() : RealVector()	{}

			PosRealVector(int indim) 	{
				dim = indim;
				vec = new double[dim];
			}

			PosRealVector(const PosRealVector& from)	{
				dim = from.dim;
				vec = new double[dim];
				for (int i=0; i<dim; i++)	{
					vec[i] = from.vec[i];
				}
			}

			PosRealVector(const double* from, int indim)	{
				dim = indim;
				vec = new double[dim];
				for (int i=0; i<dim; i++)	{
					vec[i] = from[i];
				}
			}

	virtual 	~PosRealVector() {}

	PosRealVector&		operator=(const PosRealVector& from)	{
				if (!dim)	{
					dim = from.dim;
					vec = new double[dim];
				}
				if (dim != from.dim)	{
					cerr << "error : non matching dimenstion for pos vectors\n";
					cerr << dim << '\t' << from.dim << '\n';
					exit(1);
					delete[] vec;
					dim = from.dim;
					vec = new double[dim];
				}
				for (int i=0; i<dim; i++)	{
					vec[i] = from.vec[i];
				}
				return *this;
			}

	double		GetMean()	const {
		double total = 0;
		for (int i=0; i<dim; i++)	{
			total += vec[i];
		}
		return total / dim;
	}

	void SetAtOne()	{
		for (int i=0; i<dim; i++)	{
			vec[i] = 1;
		}
	}

	double		GetVar()	const {
		double mean = 0;
		double var = 0;
		for (int i=0; i<dim; i++)	{
			var += vec[i] * vec[i];
			mean += vec[i];
		}
		mean /= dim;
		var /= dim;
		var -= mean * mean;
		return var;
	}

	double		GetEntropy()	const {
		double total = 0;
		for (int i=0; i<dim; i++)	{
			total += vec[i];
		}
		double ent = 0;
		for (int i=0; i<dim; i++)	{
			double tmp = vec[i]/total;
			ent += (tmp>1e-8) ? -tmp*log(tmp) : 0;
		}
		return ent;
	}


	double	ProposeMove(double tuning, int n)	{
		if ((n<=0) || (n > dim))	{
			n = dim;
		}
		int* indices = new int[n];
		Random::DrawFromUrn(indices,n,dim);
		double ret = 0;
		for (int i=0; i<n; i++)	{
			double m = tuning * (Random::Uniform() - 0.5);
			vec[indices[i]] *= exp(m);
			ret += m;
		}
		delete[] indices;
		return ret;
	}

	double ProposeMove(double tuning)	{
		return ProposeMove(tuning,dim);
	}

	int		ScalarMultiplication(double d)	{
		for (int i=0; i<dim; i++)	{
			vec[i] *= d;
		}
		return dim;
	}
};


class IntVector : public BaseType	{

	protected:

	int dim;
	int* vec;

	public:
			IntVector() : dim(0), vec(0) {}

			IntVector(int indim)	{
				dim = indim;
				vec = new int[dim];
			}

			IntVector(const IntVector& from)	{
				dim = from.dim;
				vec = new int[dim];
				for (int i=0; i<dim; i++)	{
					vec[i] = from.vec[i];
				}
			}

			IntVector(const int* from, int indim)	{
				dim = indim;
				vec = new int[dim];
				for (int i=0; i<dim; i++)	{
					vec[i] = from[i];
				}
			}

	virtual 	~IntVector() {
				delete[] vec;
			}

	IntVector&		operator=(const IntVector& from)	{
				if (!dim)	{
					dim = from.dim;
					vec = new int[dim];
				}
				if (dim != from.dim)	{
					cerr << "error : non matching dimenstion for vectors\n";
					cerr << dim << '\t' << from.dim << '\n';
					exit(1);
					delete[] vec;
					dim = from.dim;
					vec = new int[dim];
				}
				for (int i=0; i<dim; i++)	{
					vec[i] = from.vec[i];
				}
				return *this;
			}

	IntVector&		operator=(const int* from)	{
				if (!dim)	{
					cerr << "error in IntVector::operator=(const int*)\n";
					exit(1);
				}
				for (int i=0; i<dim; i++)	{
					vec[i] = from[i];
				}
				return *this;
			}

	const int*		GetArray() const {return vec;}

	int&		operator[](int i)	{
				return vec[i];
			}

	int&		operator[](int i) const {
				return vec[i];
			}

	int		GetDim() {return dim;}
	int		Check() {return 1;}

	double 		GetMean() const	{
		double total = 0;
		for (int i=0; i<dim; i++)	{
			total += vec[i];
		}
		return total / dim;
	}

	double		GetVar() const 	{
		double mean = 0;
		double var = 0;
		for (int i=0; i<dim; i++)	{
			var += vec[i] * vec[i];
			mean += vec[i];
		}
		mean /= dim;
		var /= dim;
		var -= mean * mean;
		return var;
	}

	int	ProposeMove(double tuning, int n)	{
		if ((n<=0) || (n > dim))	{
			n = dim;
		}
		int* indices = new int[n];
		Random::DrawFromUrn(indices,n,dim);
		for (int i=0; i<n; i++)	{
			vec[indices[i]] += (int) (tuning * (Random::Uniform() - 0.5));
		}
		delete[] indices;
		return 0;
	}

	double ProposeMove(double tuning)	{
		return ProposeMove(tuning,dim);
	}

	friend ostream& operator<<(ostream& os, const IntVector& r)  {
		os << r.dim;
		for (int i=0; i<r.dim; i++)	{
			os << '\t' << r.vec[i];
		}
		return os;
	}

	friend istream& operator>>(istream& is, IntVector& r)  {
		int indim;
		is >> indim;
		if (r.dim != indim)	{
			r.dim = indim;
			delete[] r.vec;
			r.vec = new int[r.dim];
		}
		for (int i=0; i<r.dim; i++)	{
			is >> r.vec[i];
		}
		return is;
	}
};


#endif // BASETYPE_H


