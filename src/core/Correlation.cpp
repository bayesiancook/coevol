
/********************

Coevol. Copyright 2010-2013 Nicolas Lartillot, Raphael Poujol

Coevol is free software: you can redistribute it and/or modify it under the terms of the GNU General
Public License
as published by the Free Software Foundation, either version 3 of the License, or (at your option)
any later version.
Coevol is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

See the GNU General Public License for more details. You should have received a copy of the GNU
General Public License
along with Coevol. If not, see <http://www.gnu.org/licenses/>.

**********************/


#include "Correlation.hpp"
#include <getopt.h>
#include <stdio.h>
#include <stdlib.h>
#include <string>
using namespace std;

Correlation::Correlation(double ci) {
    burnin = 0;
    nbsample = 0;
    nbparameter = 0;
    if (ci == -1) {
        CI = defaultCI;
    } else {
        CI = ci;
    }
    isSort = No;
    parameterName = NULL;
    parameters = NULL;
    sortparameters = NULL;
    covparam = NULL;
    covnorm = NULL;
    meanparam = NULL;
    isConstant = NULL;
    effectiveSize = NULL;
    variance = NULL;
    CIbuffer = NULL;
}

void Correlation::init() {
    parameterName = new string[nbparameter];
    parameters = new double*[nbparameter];
    sortparameters = new double*[nbparameter];
    covparam = new double*[nbparameter];
    covnorm = new double*[nbparameter];
    meanparam = new double[nbparameter];
    isConstant = new Switch[nbparameter];
    effectiveSize = new double[nbparameter];
    variance = new double[nbparameter];
    CIbuffer = new double*[nbparameter];
    for (int i = 0; i < nbparameter; i++) {
        CIbuffer[i] = new double[2];
        parameters[i] = NULL;
        sortparameters[i] = NULL;
        covparam[i] = NULL;
        covnorm[i] = NULL;
    }
    weight = NULL;
}

Correlation::~Correlation() {
    for (int i = 0; i < nbparameter; i++) {
        if (parameters[i] != NULL) delete[] parameters[i];
        if (sortparameters[i] != NULL) delete[] sortparameters[i];
        if (CIbuffer[i] != NULL) delete[] CIbuffer[i];
        if (covparam[i] != NULL) delete[] covparam[i];
        if (covnorm[i] != NULL) delete[] covnorm[i];
    }
    delete[] CIbuffer;
    delete[] covparam;
    delete[] parameters;
    delete[] sortparameters;
    delete[] meanparam;
    if (weight != NULL) delete[] weight;
    delete[] effectiveSize;
    delete[] isConstant;
    delete[] variance;
    delete[] covnorm;
    delete[] parameterName;
}

void Correlation::reccurquicksort(double* to, double* buf, int deb, int fin) {
    if (fin - deb > 1) {
        int last = fin;
        int first = deb;
        double pivot = to[deb];
        for (int i = deb + 1; i <= fin; i++) {
            if (to[i] > pivot) {
                buf[last] = to[i];
                last--;
            } else {
                buf[first] = to[i];
                first++;
            }
        }
        if (first == last)
            buf[first] = pivot;
        else {
            cerr << "Quiksort: erreur, ou mettre le pivot ???\n";
            cerr.flush();
            exit(0);
        }
        for (int i = deb; i <= fin; i++) {
            to[i] = buf[i];
        }
        if (deb != first) reccurquicksort(to, buf, deb, first - 1);
        if (fin != first) reccurquicksort(to, buf, first + 1, fin);
    } else {
        if (fin != deb) {
            double pivot;
            if (to[fin] < to[deb]) {
                pivot = to[fin];
                to[fin] = to[deb];
                to[deb] = pivot;
            }
        }
    }
}

void Correlation::quicksort(double* from, double* to, int s) {
    for (int i = 0; i < s; i++) to[i] = from[i];
    double* tmp = new double[s];
    reccurquicksort(to, tmp, 0, s - 1);
    delete[] tmp;
}

void Correlation::createWeight() { weight = new double[nbsample / 2]; }

void Correlation::createCovarianceBuffer() {
    if (nbsample != 0) {
        for (int i = 0; i < nbparameter; i++) {
            covnorm[i] = new double[nbsample / 2];
            for (int j = 0; j < nbsample / 2; j++) covnorm[i][j] = 0;
            covparam[i] = new double[nbsample / 2];
            for (int j = 0; j < nbsample / 2; j++) covparam[i][j] = 0;
        }
    } else {
        cerr << "ERROR: in Correlation::createCovarianceBuffer, chain size is 0, exit\n";
        cerr.flush();
        exit(0);
    }
}

void Correlation::createParameterBuffer() {
    if (nbsample != 0) {
        for (int i = 0; i < nbparameter; i++) {
            parameters[i] = new double[nbsample];
            sortparameters[i] = new double[nbsample];
        }
    } else {
        cerr << "ERROR: in Correlation::createParameterBuffer, chain size is 0, exit\n";
        cerr.flush();
        exit(0);
    }
}


void Correlation::getParameters(string filename, int start, int stop) {
    chainName = filename;
    double tmp;
    string strtmp;
    try {
        if (start == -1) {
            burnin = stop / 5;
        } else {
            burnin = start;
        }
        nbsample = stop - start;
        if (nbsample <= 0) {
            cerr << "ERROR: in Correlation::getParameters, asking sampling from point " << start
                 << " to point " << stop << ", exiting\n";
            cerr.flush();
            throw(0);
        }
        cerr << filename << '\t' << " burnin : " << burnin << '\t' << "sample size : " << nbsample
             << '\n';
        ifstream* is = new ifstream;
        is->open(filename.c_str(), ifstream::in);

        char line[100000];
        string strline;
        is->getline(line, 100000);
        strline = line;
        istringstream iss1(strline);

        nbparameter = 0;

        /*ALL*/
        /*
    iss1 >> strtmp;
    iss1 >> strtmp;
    iss1 >> strtmp;
        */
        /*ALL*/

        strtmp = "null";

        while (iss1.good()) {
            iss1 >> strtmp;
            nbparameter++;
        }
        init();
        istringstream iss2(strline);

        /*ALL*/
        /*
    iss2 >> strtmp;
    iss2 >> strtmp;
    iss2 >> strtmp;
        */
        /*ALL*/

        for (int j = 0; j < nbparameter; j++) iss2 >> parameterName[j];

        createParameterBuffer();
        createWeight();

        for (int i = 0; i < stop; i++) {
            /*ALL*/
            /*
            *is >> tmp;
            *is >> tmp;
            *is >> tmp;
            */
            /*ALL*/

            for (int j = 0; j < nbparameter; j++) {
                *is >> tmp;

                if (i >= burnin) {
                    parameters[j][i - burnin] = tmp;
                }
            }
        }
        delete is;
    } catch (...) {
        cerr << "ERROR while reading " << filename << "\n";
        cerr.flush();
        exit(0);
    }
}

void Correlation::sortParameters() {
    for (int i = 0; i < nbparameter; i++) {
        quicksort(parameters[i], sortparameters[i], nbsample);
    }
    isSort = Yes;
}

void Correlation::getCI(double* distrib, int s, double ci, double* inf, double* sup) {
    int iinf, isup;
    if (ci == 100) {
        iinf = 0;
        isup = s - 1;
    } else {
        iinf = (int)(0.5 * (1 - 0.01 * ci) * s);
        isup = (int)(s - (0.5 * (1 - 0.01 * ci) * s));
        isup--;
    }
    *inf = distrib[iinf];
    *sup = distrib[isup];
}

void Correlation::computeCI(double ci) {
    if (isSort == No) sortParameters();
    if (ci == -1) {
        ci = defaultCI;
    }
    for (int i = 0; i < nbparameter; i++) {
        getCI(sortparameters[i], nbsample, ci, &CIbuffer[i][0], &CIbuffer[i][1]);
    }
}

void Correlation::computeMean() {
    for (int j = 0; j < nbparameter; j++) meanparam[j] = 0;
    for (int j = 0; j < nbparameter; j++) {
        isConstant[j] = Yes;
        for (int i = 0; i < nbsample; i++) {
            if (parameters[j][i] != parameters[j][0]) isConstant[j] = No;
            meanparam[j] += parameters[j][i];
        }
    }
    for (int j = 0; j < nbparameter; j++) {
        if (isConstant[j] == No)
            meanparam[j] /= (double)nbsample;
        else
            meanparam[j] = parameters[j][0];
    }
}

void Correlation::computeCovariance() {
    if (nbsample == 0) {
        cerr << "ERROR: in Correlation::computeCovariance, chain size is 0, exit\n";
        cerr.flush();
        exit(0);
    } else {
        createCovarianceBuffer();
        computeMean();

        for (int t = 0; t < nbsample / 2; t++) {
            for (int j = 0; j < nbparameter; j++) {
                for (int i = 0; i < nbsample - t; i++) {
                    covparam[j][t] +=
                        (parameters[j][i] - meanparam[j]) * (parameters[j][i + t] - meanparam[j]);
                }
                if (isConstant[j] == Yes) covparam[j][t] = 0;
            }
            if (t == 0) {
                for (int j = 0; j < nbparameter; j++) {
                    variance[j] = covparam[j][0] / (double)(nbsample - 1);
                }
            }
            for (int j = 0; j < nbparameter; j++) {
                covparam[j][t] /= (double)(nbsample);
                // normalisation
                if (covparam[j][0] != 0)
                    covnorm[j][t] = covparam[j][t] / covparam[j][0];
                else
                    covnorm[j][t] = 0;
            }
        }  // fin for t
    }
}


void Correlation::computeWeight() {
    if (nbsample == 0) {
        cerr << "ERROR: in Correlation::computeWeight, chain size is 0, exit\n";
        cerr.flush();
        exit(0);
    } else {
        for (int i = 0; i < nbsample / 2; i++) {
            weight[i] = 0.5 * (1 + cos((2 * Pi * i) / (double)(nbsample)));
        }
    }
}

void Correlation::computeEffectiveSize() {
    if (nbsample == 0) {
        cerr << "ERROR: in Correlation::computeEffectiveSize, chain size is 0, exit\n";
        cerr.flush();
        exit(0);
    } else {
        computeWeight();
        double sum;
        for (int i = 0; i < nbparameter; i++) {
            sum = 1;  // covparam[i][0];*covparam[i][0];
            for (int j = 1; j < nbsample / 2; j++) {
                sum += 2 * weight[j] * covnorm[i][j];
            }
            effectiveSize[i] = nbsample / sum;
            if (sum < 1) effectiveSize[i] = nbsample;
        }
    }
}

void Correlation::outputResult() {
    string output = chainName;
    ofstream os(output.c_str());
    os << "\nparam\t\t\tmean\tSD\tSE \teff size\t " << CI << "%CI\n";
    os << std::setprecision(4) << std::setiosflags(std::ios::fixed | std::ios::showpoint)
       << /*0.123456789 <<*/ "\n";
    string space;
    for (int i = 0; i < nbparameter; i++) {
        space = "  ";
        os << i << " " << parameterName[i] << space << "\t" << meanparam[i] << "\t"
           << sqrt(variance[i]) << "\t" << sqrt(variance[i] / effectiveSize[i]) << "\t"
           << effectiveSize[i] << "\t" << CIbuffer[i][0] << "\t" << CIbuffer[i][1] << "\n";
    }
    os.close();
}

int Correlation::getNbParameter() { return nbparameter; }

string Correlation::getParameterName(int n) {
    if (n < nbparameter)
        return parameterName[n];
    else
        return "NULL";
}

double Correlation::getVariance(int n) { return variance[n]; }

double Correlation::getEffectiveSize(int n) { return effectiveSize[n]; }

double Correlation::getMean(int n) { return meanparam[n]; }

double Correlation::getInfCI(int n) { return CIbuffer[n][0]; }

double Correlation::getSupCI(int n) { return CIbuffer[n][1]; }

double Correlation::getNbSample() { return nbsample; }

// below: old Correl.cpp file

void compareChainsConvergence(string outname, Correlation** corr, int nchain, double& disc,
                              double& overlap, double& effsize) {
    int nbp = 0;
    for (int chain = 0; chain < nchain; chain++) {
        int n = corr[chain]->getNbParameter();
        if (!chain) {
            nbp = n;
        } else {
            if (nbp != n) {
                cerr << "ERROR in compChain compareChainsConvergence(), the two chain have not the "
                        "same dimensionality\n";
                cerr.flush();
                exit(0);
            }
        }
    }

    // int nbyes=0;
    // double tmp;
    double* meaneffsize = new double[nbp];
    double* absDiff = new double[nbp];
    double* meanse = new double[nbp];
    double* discrepancy = new double[nbp];

    double* maxCIoverlap = new double[nbp];

    //  int nbsesd[10];
    int nbrecci[20];

    for (int j = 0; j < 20; j++) nbrecci[j] = 0;
    // for(int i=0;i<10;i++)
    // nbsesd[i]=0;
    double* mean = new double[nchain];
    double* se = new double[nchain];

    ofstream os_conv((outname + ".contdiff").c_str());
    os_conv << "name                effsize\trel_diff\n";
    os_conv << '\n';

    disc = 0;
    overlap = 100;
    effsize = corr[0]->getNbSample();

    for (int i = 0; i < nbp; i++) {
        for (int chain = 0; chain < nchain; chain++) {
            mean[chain] = corr[chain]->getMean(i);
            se[chain] = corr[chain]->getVariance(i);
        }
        absDiff[i] = 0;
        for (int chain1 = 0; chain1 < nchain; chain1++) {
            for (int chain2 = chain1 + 1; chain2 < nchain; chain2++) {
                double tmp = fabs(mean[chain2] - mean[chain1]);
                if (absDiff[i] < tmp) {
                    absDiff[i] = tmp;
                }
            }
        }

        double max = 0;
        for (int chain = 0; chain < nchain; chain++) {
            if (max < mean[chain]) {
                max = mean[chain];
            }
        }

        meanse[i] = 0;
        for (int chain = 0; chain < nchain; chain++) {
            if (se[chain] > 1e-12) {
                meanse[i] += sqrt(se[chain]);
            }
        }
        meanse[i] /= nchain;
        if (meanse[i] > 1e-12) {
            discrepancy[i] = absDiff[i] / meanse[i];
        } else {
            discrepancy[i] = 0;
        }

        meaneffsize[i] = 0;
        for (int chain = 0; chain < nchain; chain++) {
            meaneffsize[i] += corr[chain]->getEffectiveSize(i);
        }
        meaneffsize[i] /= nchain;

        double* ciinf = new double[nchain];
        double* cisup = new double[nchain];
        for (int chain = 0; chain < nchain; chain++) {
            ciinf[chain] = corr[chain]->getInfCI(i);
            cisup[chain] = corr[chain]->getSupCI(i);
        }

        maxCIoverlap[i] = 0;
        for (int chain1 = 0; chain1 < nchain; chain1++) {
            for (int chain2 = chain1 + 1; chain2 < nchain; chain2++) {
                double CIoverlap, inf, sup, medi, meds;

                double ciinf1 = ciinf[chain1];
                double ciinf2 = ciinf[chain2];
                double cisup1 = cisup[chain1];
                double cisup2 = cisup[chain2];

                if (ciinf1 < ciinf2) {
                    inf = ciinf1;
                    medi = ciinf2;
                } else {
                    inf = ciinf2;
                    medi = ciinf1;
                }
                if (cisup1 < cisup2) {
                    sup = cisup2;
                    meds = cisup1;
                } else {
                    sup = cisup1;
                    meds = cisup2;
                }
                if (medi > meds)
                    CIoverlap = 0;  // aucun recouvrement
                else {
                    if (sup - inf != 0)
                        CIoverlap = (meds - medi) * 100 / (sup - inf);
                    else
                        CIoverlap = 100;
                }
                for (int j = 0; j < 20; j++) {
                    if (CIoverlap <= 5 * (j + 1)) {
                        nbrecci[j]++;
                        break;
                    }
                }

                if (maxCIoverlap[i] < CIoverlap) {
                    maxCIoverlap[i] = CIoverlap;
                }
            }
        }
        os_conv << corr[0]->getParameterName(i);
        for (unsigned int k = corr[0]->getParameterName(i).length(); k < 20; k++) os_conv << ' ';
        os_conv << (int)meaneffsize[i] << '\t' << '\t' << discrepancy[i] << "\n";

        if (disc < discrepancy[i]) {
            disc = discrepancy[i];
        }
        if (effsize > meaneffsize[i]) {
            effsize = meaneffsize[i];
        }
        if (overlap > maxCIoverlap[i]) {
            overlap = maxCIoverlap[i];
        }
    }
    os_conv.close();
    string cat = "cat " + outname + ".contdiff";
    system(cat.c_str());
}


int SamCompare(int nchain, int burnin, int stop, string* ChainName, double& disc, double& overlap,
               double& effsize, string outname) {
    // cerr << "burnin : " << burnin << '\n';
    // cerr << "stop : " << stop << '\n';
    Correlation** corr = new Correlation*[nchain];
    for (int chain = 0; chain < nchain; chain++) {
        corr[chain] = new Correlation();
        corr[chain]->getParameters(ChainName[chain], burnin, stop);
        corr[chain]->computeCovariance();
        corr[chain]->computeEffectiveSize();
        corr[chain]->computeCI(95);
        // corr[chain]->outputResult();
    }

    compareChainsConvergence(outname, corr, nchain, disc, overlap, effsize);

    for (int chain = 0; chain < nchain; chain++) {
        delete corr[chain];
    }
    delete[] corr;

    return 1;
}
