#include "core/RandomTypes.hpp"
#include <cmath>
using namespace std;


//-------------------------------------------------------------------------
//-------------------------------------------------------------------------
//	* Normal
//-------------------------------------------------------------------------
//-------------------------------------------------------------------------
Normal::Normal(Var<Real>* inmean, Var<PosReal>* invariance) {
    meanvec = nullptr;
    mean = inmean;
    variance = invariance;
    Register(mean);
    Register(variance);
    Sample();
}

Normal::Normal(Var<RealVector>* inmeanvec, Var<PosReal>* invariance, int inindex) {
    mean = nullptr;
    meanvec = inmeanvec;
    variance = invariance;
    index = inindex;
    Register(meanvec);
    Register(variance);
    Sample();
}

double Normal::logProb() {
    if ((!variance) || (!variance->val())) {
        return 0;
    }
    if (meanvec) {
        return -0.5 * log(2 * M_PI * variance->val()) -
               0.5 * (this->val() - (*meanvec)[index]) * (this->val() - (*meanvec)[index]) /
                   variance->val();
    }
    return -0.5 * log(2 * M_PI * variance->val()) -
           0.5 * (this->val() - mean->val()) * (this->val() - mean->val()) / variance->val();
}

void Normal::drawSample() {
    if (meanvec) {
        if ((!variance) || (!variance->val())) {
            // this is just for starting the model
            setval((*meanvec)[index] + Random::sNormal() * sqrt(10));
        } else {
            setval((*meanvec)[index] + Random::sNormal() * sqrt(variance->val()));
        }
    } else {
        if ((!variance) || (!variance->val())) {
            // this is just for starting the model
            setval(mean->val() + Random::sNormal() * sqrt(10));
        } else {
            setval(mean->val() + Random::sNormal() * sqrt(variance->val()));
        }
    }
}


//-------------------------------------------------------------------------
//-------------------------------------------------------------------------
//	* Gamma
//-------------------------------------------------------------------------
//-------------------------------------------------------------------------
const double Gamma::GAMMAMIN = 1e-20;

Gamma::Gamma(Var<PosReal>* inshape, Var<PosReal>* inscale, bool inmeanvar) {
    SetName("gamma");
    meanvar = inmeanvar;
    scale = inscale;
    shape = inshape;
    Register(shape);
    Register(scale);
    Sample();
}

Gamma::Gamma(const Gamma& from) : Rvar<PosReal>() {
    SetName("gamma");
    meanvar = from.meanvar;
    scale = from.scale;
    shape = from.shape;
    Register(shape);
    Register(scale);
    Sample();
}

inline void Gamma::drawSample() {
    double v = 0;
    if (meanvar) {
        double mean = shape->val();
        double var = scale->val();
        double alpha = mean * mean / var;
        double beta = mean / var;
        v = Random::Gamma(alpha, beta);
    } else {
        v = Random::Gamma(shape->val(), scale->val());
    }
    if (v < GAMMAMIN) {
        v = GAMMAMIN;
    }
    setval(v);
}

double Gamma::logProb() {
    if (meanvar) {
        double mean = shape->val();
        double var = scale->val();
        double alpha = mean * mean / var;
        double beta = mean / var;
        return -Random::logGamma(alpha) + (alpha)*log(beta) - beta * *this -
               (1 - alpha) * log(*this);
    }
    return -Random::logGamma(shape->val()) + (shape->val()) * log(scale->val()) -
           *this * (scale->val()) - (1 - (shape->val())) * log(*this);
}


//-------------------------------------------------------------------------
//-------------------------------------------------------------------------
//	* Beta
//-------------------------------------------------------------------------
//-------------------------------------------------------------------------
Beta::Beta(Var<PosReal>* inalpha, Var<PosReal>* inbeta) {
    alpha = inalpha;
    beta = inbeta;
    Register(alpha);
    Register(beta);
    Sample();
}

Beta::Beta(const Beta& from) : Rvar<UnitReal>() {
    alpha = from.alpha;
    beta = from.beta;
    Register(alpha);
    Register(beta);
    Sample();
}

inline void Beta::drawSample() {
    double value = 0;
    int count = 0;
    while ((!value) && (count < 10000)) {
        double x = Random::sGamma(alpha->val());
        double y = Random::sGamma(beta->val());
        value = x / (x + y);
        count++;
    }
    if (!value) {
        cerr << "error : zero beta random variable\n";
        exit(1);
    }
    setval(value);
}

double Beta::logProb() {
    return Random::logGamma(alpha->val() + beta->val()) - Random::logGamma(alpha->val()) -
           Random::logGamma(beta->val()) + (alpha->val() - 1) * log(*this) +
           (beta->val() - 1) * log(1 - *this);
}


//-------------------------------------------------------------------------
//-------------------------------------------------------------------------
//	* Exponential
//-------------------------------------------------------------------------
//-------------------------------------------------------------------------
Exponential::Exponential(Var<PosReal>* inscale, ParentType intype) {
    SetName("exponential");
    scale = inscale;
    type = intype;
    Register(scale);
    Sample();
}

Exponential::Exponential(const Exponential& from) : Rvar<PosReal>() {
    SetName("exponential");
    scale = from.scale;
    Register(scale);
    Sample();
}

void Exponential::drawSample() {
    if (type == MEAN) {
        setval(Random::sExpo() * scale->val());
    } else {
        setval(Random::sExpo() / scale->val());
    }
}

double Exponential::logProb() {
    if (type == MEAN) {
        return -*this / scale->val() - log(scale->val());
    }
    return -(scale->val()) * val() + log(scale->val());
}


//-------------------------------------------------------------------------
//-------------------------------------------------------------------------
//	* PosUniform
//-------------------------------------------------------------------------
//-------------------------------------------------------------------------
PosUniform::PosUniform(Var<PosReal>* inroot, double inmax) {
    root = inroot;
    max = inmax;
    Register(root);
    Sample();
}

PosUniform::PosUniform(const PosUniform& from) : Rvar<PosReal>() {
    root = from.root;
    max = from.max;
    Register(root);
    Sample();
}

void PosUniform::drawSample() {
    if (max > 0) {
        setval(Random::Uniform() * max);
    } else {
        setval(1.0);
    }
}

double PosUniform::logProb() {
    if (val() > max) {
        return log(0);
    }
    return 0;
}


//-------------------------------------------------------------------------
//-------------------------------------------------------------------------
//	* Binomial
//-------------------------------------------------------------------------
//-------------------------------------------------------------------------

Binomial::Binomial(int inN, Var<UnitReal>* intheta) {
    N = inN;
    theta = intheta;
    Register(theta);
    Sample();
}

void Binomial::drawSample() {
    int val = 0;
    for (int k = 0; k < N; k++) {
        if (Random::Uniform() < theta->val()) {
            val++;
        }
    }
    setval(val);
}

double Binomial::logProb() {
    double ret = *this * log(theta->val()) + (N - *this) * log(1 - theta->val());
    for (int i = 1; i <= *this; i++) {
        ret -= log((double)i);
        ret += log((double)(N + 1 - i));
    }
    return ret;
}

double Binomial::ProposeMove(double) {  // FIXME unused parameter!
    bkvalue = *this;
    flag = false;
    drawSample();
    return 0;
}


//-------------------------------------------------------------------------
//-------------------------------------------------------------------------
//	* Poisson
//-------------------------------------------------------------------------
//-------------------------------------------------------------------------
Poisson::Poisson(Var<PosReal>* inmu) {
    SetName("poisson");
    mu = inmu;
    Register(mu);
    Sample();
}

void Poisson::drawSample() {
    double p = exp(mu->val()) * Random::Uniform();
    int n = 0;
    double ratio = 1;
    double total = 1;
    while (total < p) {
        n++;
        ratio *= mu->val() / n;
        total += ratio;
    }
    setval(n);
}

double Poisson::logProb() {
    if (this->val() < 0) {
        return -200;
    }
    double ret = -mu->val() + *this * log(mu->val());
    for (int i = 2; i <= *this; i++) {
        ret -= log((double)i);
    }
    return ret;
}
/*double Poisson::ProposeMove(double tuning)	{
        bkvalue = *this;
        flag = false;
        drawSample();
        return 0;
}*/


//-------------------------------------------------------------------------
//-------------------------------------------------------------------------
//	* Dirichlet
//-------------------------------------------------------------------------
//-------------------------------------------------------------------------
Dirichlet::Dirichlet(int dimension) {
    setval(Profile(dimension));
    bkvalue = Profile(dimension);
    center = nullptr;
    concentration = nullptr;
    if (dimension == 1) {
        (*this)[0] = 1;
        Clamp();
    }
    Sample();
}

Dirichlet::Dirichlet(Var<Profile>* incenter, Var<PosReal>* inconcentration) {
    setval(Profile(incenter->val()));
    bkvalue = Profile(incenter->val());
    center = incenter;
    concentration = inconcentration;
    Register(center);
    Register(concentration);
    Sample();
}

void Dirichlet::drawSample() {
    Profile& profile = *this;
    double total = 0;
    for (int k = 0; k < GetDim(); k++) {
        if (concentration) {
            profile[k] = Random::sGamma(concentration->val() * center->val()[k]);
        } else {
            profile[k] = Random::sGamma(1);
        }
        if (profile[k] < Profile::MIN) {
            profile[k] = Profile::MIN;
        }
        total += profile[k];
    }
    for (int k = 0; k < GetDim(); k++) {
        profile[k] /= total;
    }
}

double Dirichlet::Move(double tuning, int n) {
    if (!isClamped()) {
        // Metropolis Hastings here
        Corrupt(true);
        double logHastings = Profile::ProposeMove(tuning, n);
        double deltaLogProb = Update();
        double logRatio = deltaLogProb + logHastings;
        bool accepted = (log(Random::Uniform()) < logRatio);
        if (!accepted) {
            Corrupt(false);
            Restore();
        }
        return (double)accepted;
    }
    return 1;
}

double Dirichlet::Move(double tuning) { return Move(tuning, GetDim()); }

double Dirichlet::ProposeMove(double tuning) { return Profile::ProposeMove(tuning, GetDim()); }

double Dirichlet::logProb() {
    double total = 0;
    if (concentration) {
        total = Random::logGamma(concentration->val());
        for (int i = 0; i < GetDim(); i++) {
            total += -Random::logGamma(concentration->val() * center->val()[i]) +
                     (concentration->val() * center->val()[i] - 1) * log((*this)[i]);
        }
    } else {
        total = Random::logGamma(GetDim());
    }
    return total;
}

Multinomial::Multinomial(Var<Profile>* inprobarray, int inN) {
    probarray = inprobarray;
    N = inN;
    setval(IntVector(inprobarray->val().GetDim()));
    bkvalue = IntVector(inprobarray->val().GetDim());
    Register(probarray);
    Sample();
}

double Multinomial::logProb() {
    double ret = Random::logGamma(N + 1);
    for (int k = 0; k < GetDim(); k++) {
        ret -= Random::logGamma((*this)[k] + 1);
        ret += (*this)[k] * log(probarray->val()[k]);
    }
    return ret;
}

void Multinomial::drawSample() {
    for (int k = 0; k < GetDim(); k++) {
        (*this)[k] = 0;
    }
    for (int i = 0; i < N; i++) {
        (*this)[Random::FiniteDiscrete(N, probarray->val().GetArray())]++;
    }
}


//-------------------------------------------------------------------------
//-------------------------------------------------------------------------
//	* FiniteDiscrete
//-------------------------------------------------------------------------
//-------------------------------------------------------------------------
FiniteDiscrete::FiniteDiscrete(Var<Profile>* inprobarray) {
    probarray = inprobarray;
    Register(probarray);
    Sample();
}

double FiniteDiscrete::logProb() { return log((*probarray)[*this]); }

void FiniteDiscrete::drawSample() {
    setval(Random::FiniteDiscrete(probarray->GetDim(), probarray->val().GetArray()));
}


//-------------------------------------------------------------------------
//-------------------------------------------------------------------------
//	* IIDExp
//-------------------------------------------------------------------------
//-------------------------------------------------------------------------
IIDExp::IIDExp(int dimension) {
    dim = dimension;
    setval(PosRealVector(dimension));
    bkvalue = PosRealVector(dimension);
    mean = nullptr;
    Sample();
}

IIDExp::IIDExp(int dimension, Var<PosReal>* inmean) {
    dim = dimension;
    setval(PosRealVector(dimension));
    bkvalue = PosRealVector(dimension);
    mean = inmean;
    Register(mean);
    Sample();
}

void IIDExp::setall(double in) {
    for (int i = 0; i < GetDim(); i++) {
        (*this)[i] = in;
    }
}

void IIDExp::drawSample() {
    double m = mean ? (double)mean->val() : 1.0;
    for (int i = 0; i < GetDim(); i++) {
        (*this)[i] = Random::sExpo() * m;
        if (!(*this)[i]) {
            cerr << "error in IIDexp: 0 exp random variable\n";
            exit(1);
        }
    }
}

double IIDExp::logProb() {
    double total = 0;
    double m = mean ? (double)mean->val() : 1;
    for (int i = 0; i < GetDim(); i++) {
        total += -(*this)[i] / m - log(m);
    }
    return total;
}


//-------------------------------------------------------------------------
//-------------------------------------------------------------------------
//	* IIDGamma
//-------------------------------------------------------------------------
//-------------------------------------------------------------------------
IIDGamma::IIDGamma(int dimension) {
    dim = dimension;
    setval(PosRealVector(dimension));
    bkvalue = PosRealVector(dimension);
    alpha = nullptr;
    beta = nullptr;
    Sample();
}

IIDGamma::IIDGamma(int dimension, Var<PosReal>* inalpha, Var<PosReal>* inbeta) {
    dim = dimension;
    setval(PosRealVector(dimension));
    bkvalue = PosRealVector(dimension);
    alpha = inalpha;
    beta = inbeta;
    Register(alpha);
    Register(beta);
    Sample();
}

void IIDGamma::setall(double in) {
    for (int i = 0; i < GetDim(); i++) {
        (*this)[i] = in;
    }
}

void IIDGamma::drawSample() {
    double a = alpha ? (double)alpha->val() : 1.0;
    double b = beta ? (double)beta->val() : 1.0;
    for (int i = 0; i < GetDim(); i++) {
        (*this)[i] = Random::Gamma(a, b);
        if ((*this)[i] < Gamma::GAMMAMIN) {
            (*this)[i] = Gamma::GAMMAMIN;
        }
    }
}

double IIDGamma::logProb() {
    double a = alpha ? (double)alpha->val() : 1.0;
    double b = beta ? (double)beta->val() : 1.0;
    double loggammaa = Random::logGamma(a);
    double logb = log(b);
    double total = 0;
    for (int i = 0; i < GetDim(); i++) {
        total += -loggammaa + a * logb - (*this)[i] * b + (a - 1) * log((*this)[i]);
    }
    return total;
}
