#include "utils/Random.hpp"
#include <sys/time.h>
#include <cmath>
#include <iostream>


/* =======================================================
   (VL) Magical constants, to be used only in this file.
   ======================================================= */
#define MT_IA 397
#define MT_IB (MT_LEN - MT_IA)
#define UPPER_MASK 0x80000000
#define LOWER_MASK 0x7FFFFFFF
#define MATRIX_A 0x9908B0DF
#define TWIST(b, i, j) ((b)[i] & UPPER_MASK) | ((b)[j] & LOWER_MASK)
#define MAGIC(s) (((s)&1) * MATRIX_A)

// #define SAFE_EXP(x) ((x)<-200.0 ? 0.0 : exp(x))
#define SAFE_EXP(x) exp(x)

const double gammacoefs[] = {0.9999999999995183,  676.5203681218835,      -1259.139216722289,
                             771.3234287757674,   -176.6150291498386,     12.50734324009056,
                             -0.1385710331296526, 0.9934937113930748e-05, 0.1659470187408462e-06};
const double Pi = 3.1415926535897932384626;
// const double Logroot2pi = 0.918938533204673;


// -------------------------------------------------
// just a trick for random number initialisation
// function to be called before entering main()
class random_init {
  public:
    random_init() {
        std::cerr << '\n';
        Random::InitRandom();
        // Random::InitRandom(625878);
        // Random::InitRandom(645823);
        // Random::InitRandom(778389);
        // Random::InitRandom(818970);
        std::cerr << "random seed : " << Random::GetSeed() << '\n';
        std::cerr << '\n';
    }
};

static random_init init;

int Random::Seed = 0;
int Random::mt_index = 0;
unsigned long long Random::mt_buffer[MT_LEN];

const double Random::INFPROB = 250;


// ---------------------------------------------------------------------------------
//		� Random()
// ---------------------------------------------------------------------------------
void Random::InitRandom(int seed) {
    if (seed == -1) {
        struct timeval tod;
        gettimeofday(&tod, nullptr);
        seed = tod.tv_usec;
    }
    Seed = seed;
    srand(seed);
    int i;

    if (RAND_MAX == 32767) {
        for (i = 0; i < MT_LEN; i++) {
            mt_buffer[i] = rand() * 32768 + rand();
        }
    } else {
        for (i = 0; i < MT_LEN; i++) {
            mt_buffer[i] = rand();
        }
    }
    mt_index = 0;
}


Random::Random(int seed) { InitRandom(seed); }

int Random::GetSeed() { return Seed; }


// ---------------------------------------------------------------------------------
//		� Uniform()
// ---------------------------------------------------------------------------------
double Random::Uniform() {
    // Mersenne twister
    // Matsumora and Nishimora 1996
    // 32-bit generator

    // implementation adapted from
    // Michael Brundage
    // copyright 1995-2005
    // creative commons

    // check that number belongs to (0,1), boundaries excluded
    double ret = 0;
    while ((ret == 0) || (ret == 1)) {
        unsigned long long* b = mt_buffer;
        int idx = mt_index;
        unsigned long long s;
        int i;

        if (idx == MT_LEN * sizeof(unsigned long long)) {
            idx = 0;
            i = 0;
            for (; i < MT_IB; i++) {
                s = TWIST(b, i, i + 1);
                b[i] = b[i + MT_IA] ^ (s >> 1) ^ MAGIC(s);
            }
            for (; i < MT_LEN - 1; i++) {
                s = TWIST(b, i, i + 1);
                b[i] = b[i - MT_IB] ^ (s >> 1) ^ MAGIC(s);
            }

            s = TWIST(b, MT_LEN - 1, 0);
            b[MT_LEN - 1] = b[MT_IA - 1] ^ (s >> 1) ^ MAGIC(s);
        }
        mt_index = idx + sizeof(unsigned long long);
        unsigned long long r = *(unsigned long long*)((unsigned char*)b + idx);

        // counfounding the bits returned to the caller is done here
        // as  in Matsumoto and Nishimura, and
        // unlike in the initial implementation of Brundage (see below)
        r ^= (r >> 11);
        r ^= (r << 7) & 0x9D2C5680;
        r ^= (r << 15) & 0xEFC60000;
        r ^= (r >> 18);
        // ret =  ((double) r) / 4294967296;
        ret = ((double)r) / 65536;
        ret /= 65536;
    }
    return ret;
    // Matsumoto and Nishimura additionally confound the bits returned to the caller
    // but this doesn't increase the randomness, and slows down the generator by
    // as much as 25%.  So I omit these operations here.

    // r ^= (r >> 11);
    // r ^= (r << 7) & 0x9D2C5680;
    // r ^= (r << 15) & 0xEFC60000;
    // r ^= (r >> 18);
}


// ---------------------------------------------------------------------------------
//		� GPoisson()
// ---------------------------------------------------------------------------------
int Random::Poisson(double mu) {
    int n = 0;
    double tottime = Random::sExpo();
    while (tottime < mu) {
        n++;
        tottime += Random::sExpo();
    }
    return n;
}


// ---------------------------------------------------------------------------------
//		� Gamma()
// ---------------------------------------------------------------------------------
int Random::ApproxBinomial(int N, double p) { return Poisson(N * p); }


// ---------------------------------------------------------------------------------
//		� Gamma()
// ---------------------------------------------------------------------------------
double Random::Gamma(double alpha, double beta) { return sGamma(alpha) / beta; }


// ---------------------------------------------------------------------------------
//		� DrawFromDiscreteDistribution()
// ---------------------------------------------------------------------------------
int Random::DrawFromDiscreteDistribution(const double* prob, int nstate) {
    try {
        double total = 0;
        for (int k = 0; k < nstate; k++) {
            total += prob[k];
        }
        double p = total * Random::Uniform();
        double tot = 0;
        int k = -1;
        do {
            k++;
            tot += prob[k];
        } while ((k < nstate) && (tot < p));
        if (k == nstate) {
            std::cerr << "finite discrete overflow\n";
            for (int k = 0; k < nstate; k++) {
                std::cerr << prob[k] << '\n';
            }
            throw;
        }
        return k;
    } catch (...) {
        std::cerr << "error in draw from \n";
        throw;
    }
}


// ---------------------------------------------------------------------------------
//		� DrawFromUrn()
// ---------------------------------------------------------------------------------
void Random::DrawFromUrn(int* tab, int n, int N) {  // draw n out of N
    // assumes that tab is an Int16[n]
    for (int i = 0; i < n; i++) {
        tab[i] = 0;
    }
    auto index = new int[N];
    for (int i = 0; i < N; i++) {
        index[i] = 0;
    }
    for (int i = 0; i < n; i++) {
        int trial = (int)(Random::Uniform() * (N - i));
        for (int k = 0; k < N; k++) {
            if (index[k] != 0) {
                if (trial >= k) {
                    trial++;
                }
            }
        }
        if (trial == N) {
            std::cerr << "error in draw from urn: overflow\n";
            exit(1);
        }
        tab[i] = trial;
        if (index[trial] != 0) {
            std::cerr << "error in draw from urn: chose twice the same\n";
            exit(1);
        }
        index[trial] = 1;
    }
    delete[] index;
}


// ---------------------------------------------------------------------------------
//		� Choose()
// ---------------------------------------------------------------------------------
int Random::Choose(int scale) { return (int)(Random::Uniform() * scale); }


int Random::FiniteDiscrete(int n, const double* probarray) {
    double total = 0;
    auto cumul = new double[n];
    for (int k = 0; k < n; k++) {
        total += probarray[k];
        cumul[k] = total;
    }
    double u = total * Random::Uniform();
    int k = 0;
    while ((k < n) && (u > cumul[k])) {
        k++;
    }
    if (k == n) {
        std::cerr << "error in Random::FiniteDiscrete\n";
        exit(1);
    }
    delete[] cumul;
    return k;
}


// ---------------------------------------------------------------------------------
//		� sNormal()
// ---------------------------------------------------------------------------------
double Random::sNormal() {
    double u = Random::Uniform();
    if (u <= 0.8638) {
        double v = 2 * Random::Uniform() - 1;
        double w = 2 * Random::Uniform() - 1;
        return 2.3153508 * u - 1 + v + w;
    }
    if ((0.8638 < u) && (u <= 0.9745)) {
        double v = Random::Uniform();
        return 1.5 * (v - 1 + 9.0334237 * (u - 0.8638));
    } else if ((0.9973002 < u) && (u <= 1)) {
        double x, v;
        do {
            v = Random::Uniform();
            double w = Random::Uniform();
            x = 4.5 - log(w);
        } while (x * v * v > 4.5);
        // double ret = (u - 0.9986501) > 0 ? sqrt(2 * x) : - sqrt(2 * x);
        double ret = sqrt(2 * x);
        if (u - 0.9986501 > 0) {
            ret = -ret;
        }
        return ret;
    }
    double x, v, w, tot;
    do {
        x = 6 * Random::Uniform() - 3;
        u = Random::Uniform();
        v = (x > 0) ? x : -x;
        w = 6.6313339 * (3 - v) * (3 - v);
        tot = 0;
        if (v < 1.5) {
            tot += 6.0432809 * (1.5 - v);
        }
        if (v < 1) {
            tot += 13.2626678 * (3 - v * v) - w;
        }
    } while (u > 49.0024445 * SAFE_EXP(-0.5 * v * v) - tot - w);
    return x;
}


// ---------------------------------------------------------------------------------
//		� sExpo()
// ---------------------------------------------------------------------------------
double Random::sExpo() { return -log(Random::Uniform()); }


double fsign(double num, double sign)
/* Transfers sign of argument sign to argument num */
{
    if ((sign > 0.0f && num < 0.0f) || (sign < 0.0f && num > 0.0f)) {
        return -num;
    } else {
        return num;
    }
}


// ---------------------------------------------------------------------------------
//		� sGamma()
// ---------------------------------------------------------------------------------
double Random::sGamma(double a) {
    if (a > 1) {
        static double a1 = 0;
        static double a2 = 0;

        static double s2, s, d, t, x, u, q0, b, sigma, c, v, q, e;

        // step 1
        if (a != a1) {
            a1 = a;
            s2 = a - 0.5;
            s = sqrt(s2);
            d = 4 * sqrt(2.0) - 12 * s;
        }

        // step 2
        t = sNormal();
        x = s + 0.5 * t;
        if (t > 0) {
            if (x == 0.0) {
                std::cerr << "1\n";
            }
            return x * x;
        }

        // step 3
        u = Uniform();
        if (d * u < t * t * t) {
            if (x == 0.0) {
                std::cerr << "2\n";
            }
            return x * x;
        }

        // step 4
        if (a != a2) {
            a2 = a;
            q0 = log(sqrt(2 * Pi)) - logGamma(a) - s2 + s2 * log(s2);
            if (a < 3.686) {
                b = 0.463 + s + 0.178 * s2;
                sigma = 1.235;
                c = 0.195 / s - 0.079 + 0.16 * s;
            } else if (a < 13.022) {
                b = 1.654 + 0.0076 * s2;
                sigma = 1.68 / s + 0.275;
                c = 0.062 / s + 0.024;
            } else {
                b = 1.77;
                sigma = 0.75;
                c = 0.1515 / s;
            }
        }

        // step 5-7
        if (x > 0) {
            v = 0.5 * t / s;
            q = q0 - s * t + 0.25 * t * t + 2 * s2 * log(1.0 + v);
            if (log(1 - u) < q) {
                if (x == 0.0) {
                    std::cerr << "3\n";
                }
                return x * x;
            }
        }

        do {
            do {
                e = sExpo();
                u = Uniform();
                u = u + u - 1;
                t = fabs(e * sigma);
                if (u < 0) {
                    t = -t;
                }
                t += b;
            } while (t <= -0.71874483771719);

            v = 0.5 * t / s;
            q = q0 - s * t + 0.25 * t * t + 2 * s2 * log(1.0 + v);
        }

        while ((q < 0) || (c * fabs(u) > (exp(q) - 1) * exp(e - 0.5 * t * t)));

        x = s + 0.5 * t;
        /*
          if (!x)	{
          std::cerr << "4\n";
          std::cerr << s << '\t' << t << '\n';
          std::cerr << e << '\t' << sigma << '\t' << u << '\n';
          std::cerr << b << '\n';
          }
        */
        return x * x;
    }


    double x, y;
    do {
        double u = Random::Uniform();
        double v = Random::Uniform();
        x = SAFE_EXP(log(u) / a);
        y = SAFE_EXP(log(v) / (1 - a));
    } while (x + y > 1);

    double e = -log(Random::Uniform());

    return e * x / (x + y);
}


// ---------------------------------------------------------------------------------
//		* logGamma()
// ---------------------------------------------------------------------------------
double Random::logGamma(double alpha) {
    // adapted from statlib
    if (alpha < 0) {
        std::cerr << "error in loggamma: only positive argument\n";
        exit(1);
    }

    double tot = gammacoefs[0];
    double f = alpha;
    for (int i = 1; i < 8; i++) {
        tot += gammacoefs[i] / (f++);
    }
    return log(tot * sqrt(2 * Pi)) - alpha - 6.5 + (alpha - 0.5) * log(alpha + 6.5);
}

double Random::logMultivariateGamma(double a, int p) {
    if (2 * a <= (p - 1)) {
        std::cerr << "error in Random::logMultivariateGamma\n";
        std::cerr << "parameter out of bound: " << a << '\t' << p << '\n';
        exit(1);
    }
    double ret = p * (p - 1) / 4 * log(Pi);
    for (int j = 1; j <= p; j++) {
        ret += logGamma(a + (1 - j) / 2);
    }
    return ret;
}
