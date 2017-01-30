
#ifndef CONJUGATEPOISSON_H
#define CONJUGATEPOISSON_H

#include "Conjugate.h"
#include "RandomTypes.h"
#include "IID.h"

// Imagine a collection of N integer variables, x_i, i=1..N, all Poisson distributed with the same mean mu
// and imagine that this mean mu is itself distributed according to a gamma, of shape parameter alpha and scale parameter beta.
//
// Each time the mean mu makes a Metropolis Hastings move,
// it has to re-calculate its own log probability p(mu) (the gamma density),
// plus the sum of the log probability of all its children, p(x|mu) = \prod_i p(x_i | mu) (a product of N Poisson densities)
// which means N+1 calls to a logProb() method.
//
// Yet, computations can be simplified here
// specifically, if we concentrate on the product of the probabilities of all N Poisson distributed children:
//
// p(x | mu)  = \prod_i p(x_i | mu) \propto \prod_i e^(-mu) mu^(x_i) = mu^(\sum_i x_i) e^(-N mu) = mu^A e^(-B mu)
//
// let us call A = \sum_i x_i the shape statistic, and B = N the scale statistic.
// they are called that way because if mu is itself a gamma, of shape parameter alpha and scale parameter beta:
//
// p(mu) \propto mu^(alpha-1) e^(-beta mu)
//
// then the posterior is proportional to:
//
// p(mu | x) \propto p(x|mu) p(mu) \propto mu^(alpha + A- 1) e^[(-beta + B) mu]
//
// i.e. the posterior is a Gamma of shape alpha + A and scale beta + B
// in particular, this means that mu can be integrated out analytically, and directly resampled conditional on alpha, beta, A and B
//
// If, on the other hand, mu is not itself distributed as a gamma (e.g. mu is log normally distributed)
// we cannot anymore make this analytical integration and this Gibbs resampling.
// But we can still compute p(x|mu) based on sufficient statistics A and B.
//
// The following classes implement the complete series of behaviors implied by this short mathematical development.
// Let us take each class in turn:
//
//
// class PoissonSemiConjugate: public SemiConjugatePrior<PosReal>
//
// A PoissonSemiConjugate is a random variable representing the mean mu (it is a Rvar<PosReal>),
// which allows for fast computation of log p(x | mu) = A log(mu) - B mu (forgetting about uninteresting additive constants)
// This log p(x | mu) is returned by the function SuffStatLogProb() declared pure virtual in SemiConjugatePrior<PosReal>.
//
//
// class ConjugatePoisson: public ConjugateSampling<Int>, public Poisson
//
// In order to collect the sufficient statistics, a PoissonSemiConjugate has to ask them to its children.
// The children should therefore be Poisson random variables, but with an additional method,
// allowing them to 'communicate' their personal contribution to the shape and scale statistics collected by their parent.
// This is what the ConjugatePoisson class does.
// It derives both from ConjugateSampling<Int> (providing the interface for communicating ones contribution to the sufficient statistics)
// and from Poisson.
//
//
// class SemiConjugateGamma : public PoissonSemiConjugate, public Gamma
//
// On the other hand, a PoissonSemiConjugate does not yet know its own probability p(mu)
// this is implemented in derived classes
// Assuming p(mu) is a Gamma, of shape parameter alpha and scale parameter beta.
// we can represent mu as a SemiConjugateGamma, deriving both from Gamma and from PoissonSemiConjugate
// When asked for its log prob, this class returns log p(mu) + log p(x|mu):
// SemiConjugateGamma::logProb() returns SuffStatLogProb() + Gamma::logProb()
// Likewise, we could define SemiConjugateLogNormal(), whose SuffStatLogProb() would return SuffStatLogProb() + LogNormal::logProb(), etc.
//
//
// class ConjugateGamma : public SemiConjugateGamma, public ConjugatePrior<Real>
//
// Going one step further: ConjugateGamma is a SemiConjugateGamma that can integrate itself out:
//
// thus, for instance, when a PARENT of ConjugateGamma makes a Metropolis Hastings Move,
// ConjugateGamma returns the log of p(x | alpha, beta).
// (the integral of the probabilty of all its children, averaged over all possible values of mu):
//
// p(x | alpha, beta) = \int \prod_i p(x_i | mu) p(mu | alpha, beta) dmu
//
// Conversely, a ConjugateGamma can GibbsResample itself back in:
//
// mu ~ p(mu | x, alpha, beta)
//
// which, given what we explained above, simply amounts to drawing a Gamma random number of shape alpha + A and scale beta + B.
// This latter behavior is implemented by GibbsResample() method (declared pure virtual in ConjugatePrior<PosReal>).
//
//
// ConjugateGamma is the only possible ConjugatePrior<Real> than can be associated to a ConjugatePoisson.
// For instance, we cannot define a ConjugateLogNormal.
// This is because the integral over mu would not be analytical, and the Gibbs resampling method not possible.
// On the other hand, as we have seen above (e.g. the SemiConjugateLogNormal case),
// one can derive several SemiConjugatePrior<PosReal>, associated to a Poisson Distribution
// More generally, to a particular ConjugateSampling<U>, one can associate many SemiConjugatePrior<T>,
// but at most one ConjugatePrior<T> (and still, there are very few cases where this is possible)
//
//
// Finally, note that in the present case, we have defined
// (1) class ConjugateGamma : public ConjugatePrior<PosReal>, public SemiConjugateGamma;
//
// but we could have defined ConjugateGamma directly, as
// (2) class ConjugateGamma : public ConjugatePrior<PosReal>, public PoissonSemiConjugate, public Gamma;
//
// or even
// (3) class ConjygateGamma : public ConjugatePrior<PosReal>, public Gamma;
// (then, defining the PoissonSemiConjugate behavior directly in ConjugateGamma)
//
// case (1) is conceptually the most elegant one, because it discriminates all successive steps of the procedure,
// but is not really useful (because we do not really need a semi conjugate gamma)
// case (3) is not optimal either, because we will often need a PoissonSemiConjugate, whose ConjugateGamma is a particular case
// therefore, case (2) is probably the best compromise
//


// the mean parameter of a poisson distribution
class PoissonSemiConjugate : public virtual SemiConjugatePrior<PosReal> {

	public:

	virtual ~PoissonSemiConjugate() {}

	void ResetSufficientStatistic()	{
		shapestat = 0;
		scalestat = 0;
		logfactstat = 0;
	}

	void SaveSufficientStatistic()	{
		bkshapestat = shapestat;
		bkscalestat = scalestat;
		bklogfactstat = logfactstat;
	}

	void RestoreSufficientStatistic()	{
		shapestat = bkshapestat;
		scalestat = bkscalestat;
		logfactstat = bklogfactstat;
	}

	void AddToShape(double in)	{
		shapestat += in;
	}

	void AddToScale(double in)	{
		scalestat += in;
	}

	void AddToLogFact(double in)	{
		logfactstat += in;
	}

	double SuffStatLogProb()	{
		return shapestat * log(val()) - scalestat * val() - logfactstat;
	}

	protected:

	double shapestat;
	double scalestat;
	double logfactstat;
	double bkshapestat;
	double bkscalestat;
	double bklogfactstat;

};

// this mean parameter is (semi-)conjugate to a Poisson sampling distribution
class ConjugatePoisson : public ConjugateSampling<Int>, public Poisson	{

	public:

	ConjugatePoisson(PoissonSemiConjugate* inmu) : Poisson(inmu) {
		conjugate_up.insert(inmu);
	}

	~ConjugatePoisson() {}

	void AddSufficientStatistic(SemiConjPrior* parent)	{
		PoissonSemiConjugate* prior = dynamic_cast<PoissonSemiConjugate*>(parent);
		if (! prior)	{
			cerr << "cast error in ConjugatePoisson::AddSuffStat\n";
			exit(1);
		}
		prior->AddToShape(val());
		prior->AddToScale(1.0);
		prior->AddToLogFact(Random::logGamma(1 + val()));
	}

};

// as an example, we can give this mean parameter a gamma distribution
class SemiConjugateGamma : public PoissonSemiConjugate, public Gamma	{

	public:

	SemiConjugateGamma(Var<PosReal>* inshape, Var<PosReal>* inscale) : Gamma(inshape,inscale) {}

	double logProb()	{
		if (isActive())	{
			return Gamma::logProb() + SuffStatLogProb();
		}
		return Gamma::logProb();
	}
};

// but in fact, we have analytical conjugacy here, so let us implement it
class ConjugateGamma : public ConjugatePrior<PosReal>, public SemiConjugateGamma {

	public:

	ConjugateGamma(Var<PosReal>* inshape, Var<PosReal>* inscale) : SemiConjugateGamma(inshape,inscale) {}

	void GibbsResample()	{
		// draw a gamma of shape alpha + A and scale beta + B
		if (shapestat < 0)	{
			cerr << "error : negative shapestat in conjugate gamma\n";
			exit(1);
		}
		if (scalestat < 0)	{
			cerr << "error : negative scalestat in conjugate gamma\n";
			exit(1);
		}

		double value = Random::Gamma(shape->val() + shapestat, scale->val() + scalestat);

		// just be careful about mathematical errors
		if (value < GAMMAMIN)	{
			value = GAMMAMIN;
		}

		this->setval(value);
	}

	double logProb()	{
		if (isActive())	{

			// if in integrated (conjugate) mode, return the integrated log probability
			if (isIntegrated())	{
				double postshape = shape->val() + shapestat;
				double postscale = scale->val() + scalestat;
				return -Random::logGamma(shape->val()) + shape->val() * log(scale->val()) + Random::logGamma(postshape) - postshape * log(postscale) - logfactstat;
			}

			// if in active but not integrated mode (i.e. semi-conjugate mode), return the semi-conjugate log probabilty
			return SemiConjugateGamma::logProb();
		}
		// if in inactive mode, return the plain log probability, as a Gamma
		return Gamma::logProb();
	}
};

class ConjugateGammaIIDArray : public IIDArray<PosReal>	{

	public:
	ConjugateGammaIIDArray(int insize, Var<PosReal>* inalpha, Var<PosReal>* inbeta) : IIDArray<PosReal>(insize)	{
		alpha = inalpha;
		beta = inbeta;
		Create();
	}

	ConjugateGamma* operator[](int site)	{
		return dynamic_cast<ConjugateGamma*>(array[site]);
	}

	protected:

	Rvar<PosReal>* CreateVal(int site)	{
		return new ConjugateGamma(alpha, beta);
	}

	Var<PosReal>* alpha;
	Var<PosReal>* beta;
};

#endif

