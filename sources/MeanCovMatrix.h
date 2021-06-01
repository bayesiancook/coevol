
#ifndef MEANCOVMATRIX_H
#define MEANCOVMATRIX_H

#include "CovMatrix.h"
#include "SumConstrained.h"


class MeanCovMatrix {

	public:

	MeanCovMatrix(int indim, bool ininverse = true) : tex(false), dim(indim), size(0), inverse(ininverse) {

		val = new vector<double>*[GetDim()];
		for (int i=0; i<GetDim(); i++)	{
			val[i] = new vector<double>[GetDim()];
		}

		invval = new vector<double>*[GetDim()];
		for (int i=0; i<GetDim(); i++)	{
			invval[i] = new vector<double>[GetDim()];
		}

		slope = new vector<double>*[GetDim()];
		for (int i=0; i<GetDim(); i++)	{
			slope[i] = new vector<double>[GetDim()];
		}

		mean = new double*[GetDim()];
		for (int i=0; i<GetDim(); i++)	{
			mean[i] = new double[GetDim()];
		}
		meaninv = new double*[GetDim()];
		for (int i=0; i<GetDim(); i++)	{
			meaninv[i] = new double[GetDim()];
		}
		meanslope= new double*[GetDim()];
		for (int i=0; i<GetDim(); i++)	{
			meanslope[i] = new double[GetDim()];
		}
		varslope= new double*[GetDim()];
		for (int i=0; i<GetDim(); i++)	{
			varslope[i] = new double[GetDim()];
		}
		var = new double*[GetDim()];
		for (int i=0; i<GetDim(); i++)	{
			var[i] = new double[GetDim()];
		}
		varinv = new double*[GetDim()];
		for (int i=0; i<GetDim(); i++)	{
			varinv[i] = new double[GetDim()];
		}
		pp = new double*[GetDim()];
		for (int i=0; i<GetDim(); i++)	{
			pp[i] = new double[GetDim()];
		}
		ppinv = new double*[GetDim()];
		for (int i=0; i<GetDim(); i++)	{
			ppinv[i] = new double[GetDim()];
		}
		correl = new double*[GetDim()];
		for (int i=0; i<GetDim(); i++)	{
			correl[i] = new double[GetDim()];
		}
		correl2 = new double*[GetDim()];
		for (int i=0; i<GetDim(); i++)	{
			correl2[i] = new double[GetDim()];
		}
		partialcorrel = new double*[GetDim()];
		for (int i=0; i<GetDim(); i++)	{
			partialcorrel[i] = new double[GetDim()];
		}
		partialcorrel2 = new double*[GetDim()];
		for (int i=0; i<GetDim(); i++)	{
			partialcorrel2[i] = new double[GetDim()];
		}

		meanpropvar = new double[GetDim()];
	}

	virtual ~MeanCovMatrix()	{
		for (int i=0; i<GetDim(); i++)	{
			delete[] mean[i];
			delete[] meaninv[i];
			delete[] meanslope[i];
			delete[] varslope[i];
			delete[] var[i];
			delete[] varinv[i];
			delete[] correl[i];
			delete[] correl2[i];
			delete[] partialcorrel[i];
			delete[] partialcorrel2[i];
			delete[] pp[i];
			delete[] ppinv[i];
			delete[] val[i];
			delete[] invval[i];
			delete[] slope[i];
		}
		delete[] mean;
		delete[] meaninv;
		delete[] meanslope;
		delete[] varslope;
		delete[] var;
		delete[] varinv;
		delete[] correl;
		delete[] correl2;
		delete[] partialcorrel;
		delete[] partialcorrel2;
		delete[] pp;
		delete[] ppinv;
		delete[] val;
		delete[] invval;
		delete[] slope;

		delete[] meanpropvar;
	}

	void SetLatex(bool b)	{
		tex = b;
	}

	int GetDim() const {
		return dim;
	}

	unsigned int GetSize() const	{
		return size;
	}

	double GetVal(int i, int j, int k) const	{
		if (val[i][j].size() != GetSize())	{
			cerr << "error : " << i << '\t' << j << '\t' << k << '\n';
			exit(1);
		}
		return val[i][j][k];
	}

	double GetInvVal(int i, int j, int k) const	{
        if (! inverse)  {
            cerr << "error in MeanCovMatrix::GetInvVal: inverse not activated\n";
            exit(1);
        }
		if (invval[i][j].size() != GetSize())	{
			cerr << "error : " << i << '\t' << j << '\t' << k << '\n';
			exit(1);
		}
		return invval[i][j][k];
	}

	virtual void Add(CovMatrix* sample)	{
        if (inverse)    {
            sample->Diagonalise();
        }
		for (int i=0; i<GetDim(); i++)	{
			for (int j=0; j<GetDim(); j++)	{
				val[i][j].push_back((*sample)[i][j]);
                if (inverse)    {
                    invval[i][j].push_back(sample->GetInvMatrix()[i][j]);
                }
			}
		}
		size++;
	}

	void ComputeSlopes()	{
		for (int i=0; i<GetDim(); i++)	{
			for (int j=0; j<GetDim(); j++)	{
				if (i != j)	{
					meanslope[i][j] = 0;
					varslope[i][j] = 0;
					for (unsigned int k=0; k<size; k++)	{
						double sxy = val[i][j][k];
						double sx = val[i][i][k];
						double sy = val[j][j][k];
						// double tmp = (sy + sqrt((sy - sx)*(sy - sx) - 4 *sxy *sxy - 2 * sx * sy)) / 2 / sxy;
						double tmp = (sy - sx + sqrt((sy - sx)*(sy - sx) + 4 *sxy *sxy)) / 2 / sxy;
						slope[i][j].push_back(tmp);
						meanslope[i][j] += tmp;
						varslope[i][j] += tmp * tmp;
					}
					meanslope[i][j] /= size;
					varslope[i][j] /= size;
					varslope[i][j] -= meanslope[i][j] * meanslope[i][j];
					sort(slope[i][j].begin(), slope[i][j].end());
				}
			}
		}
	}

	void Normalize()	{
		for (int i=0; i<GetDim(); i++)	{
			meanpropvar[i] = 0;
			for (int j=0; j<GetDim(); j++)	{
				mean[i][j] = 0;
				meaninv[i][j] = 0;
				var[i][j] = 0;
				varinv[i][j] = 0;
				pp[i][j] = 0;
				ppinv[i][j] = 0;
				correl[i][j] = 0;
				correl2[i][j] = 0;
				partialcorrel[i][j] = 0;
				partialcorrel2[i][j] = 0;
				if (val[i][j].size() != size)	{
					cerr << "error in meancovmatrix: non matching length\n";
					exit(1);
				}
				for (unsigned int k=0; k<size; k++)	{
					double& tmp = val[i][j][k];
					mean[i][j] += tmp;
					var[i][j] += tmp * tmp;
					if (tmp > 0)	{
						pp[i][j]++;
					}
					double r = tmp / sqrt(val[i][i][k] * val[j][j][k]);
					correl[i][j] += r;
					correl2[i][j] += r*r;
                    if (inverse)    {
                        double& inv = invval[i][j][k];
                        meaninv[i][j] += inv;
                        varinv[i][j] += inv * inv;
                        if (inv < 0)	{
                            ppinv[i][j] ++;
                        }
                        double pr = - inv / sqrt(invval[i][i][k] * invval[j][j][k]);
                        partialcorrel[i][j] += pr;
                        partialcorrel2[i][j] += pr * pr;

                        meanpropvar[i] += 1.0 - 1.0 /(val[i][i][k] * invval[i][i][k]);
                    }
				}
				mean[i][j] /= size;
				var[i][j] /= size;
				var[i][j] -= mean[i][j] * mean[i][j];
				pp[i][j] /= size;
                if (inverse)    {
                    meaninv[i][j] /= size;
                    varinv[i][j] /= size;
                    varinv[i][j] -= meaninv[i][j] * meaninv[i][j];
                    ppinv[i][j] /= size;
                }
				correl[i][j] /= size;
				correl2[i][j] /= size;
                if (inverse)    {
                    partialcorrel[i][j] /= size;
                    partialcorrel2[i][j] /= size;
                    // sort(val[i][j].begin(), val[i][j].end());
                    meanpropvar[i] /= size;
                }
			}
		}
		ComputeSlopes();
	}

	double GetPropVariance(int i)	const {
		return meanpropvar[i];
		// return 1.0 - 1.0 / (mean[i][i] * meaninv[i][i]);
	}

	// MeanCovMatrix Project(bool* array);

	void PrintPropVariances(ostream& os) const {

		os << "proportion of variance of each trait explained by all other traits:\n";
		for (int i=0; i<GetDim(); i++)	{
			os << "trait " << i << " : " << GetPropVariance(i) << '\n';
		}
		os << '\n';
	}

	void PrintCovariances(ostream& os) const	{
		if (tex)	{
			os << "{\\sc covariance}";
			for (int i=0; i<GetDim(); i++)	{
				os << "&" << i;
			}
			os << "\\\\\n";
			os << "\\hline\n";
			os << '\n';
			for (int i=0; i<GetDim(); i++)	{
				os << i;
				for (int j=0; j<i; j++)	{
					os << "&-";
				}
				for (int j=i; j<GetDim(); j++)	{
					os << " & ";
					os << "$";
					os << ((double) ((int) (100 * mean[i][j]))) / 100;
					if ((pp[i][j] > 0.975) || (pp[i][j] < 0.025))	{
						os << "^*";
					}
					os << "$";
				}
				os << "\\\\\n";
			}
			os << "\\\\\n";
		}
		else	{
			os.precision(3);
			// output matrix of covariance
			os << "covariances\n";
			os << '\n';
			for (int i=0; i<GetDim(); i++)	{
				for (int j=0; j<GetDim(); j++)	{
					os << setw(7) << mean[i][j] << '\t';
				}
				os << '\n';
			}
			os << '\n';
		}
	}

	void PrintPrecisions(ostream& os) const	{
		os.precision(3);
		// output matrix of covariance
		os << "precisions\n";
		os << '\n';
		for (int i=0; i<GetDim(); i++)	{
			for (int j=0; j<GetDim(); j++)	{
				os << setw(7) << meaninv[i][j] << '\t';
			}
			os << '\n';
		}
		os << '\n';

	}

	void PrintR2(ostream& os) const	{
		if (tex)	{
			os << "{\\sc correlation ($R$)}";
			for (int i=0; i<GetDim(); i++)	{
				os << "&" << i;
			}
			os << "\\\\\n";
			os << "\\hline\n";
			os << '\n';
			for (int i=0; i<GetDim(); i++)	{
				os << i;
				for (int j=0; j<=i; j++)	{
					os << "&-";
				}
				for (int j=i+1; j<GetDim(); j++)	{
					os << " & ";
					os << "$";
					os << ((double) ((int) (100 * correl[i][j]))) / 100;
					// os << ((double) ((int) (100 * correl[i][j] * correl[i][j]))) / 100;
					if ((pp[i][j] > 0.975) || (pp[i][j] < 0.025))	{
						os << "^{**}";
					}
					else if ((pp[i][j] > 0.95) || (pp[i][j] < 0.05))	{
						os << "^*";
					}
					os << "$";
				}
				os << "\\\\\n";
			}
			os << "\\\\\n";
		}
		else	{

			os.precision(3);

			// correlation coefficients
			os << "r2\n";
			os << '\n';
			for (int i=0; i<GetDim(); i++)	{
				for (int j=0; j<GetDim(); j++)	{
					os << setw(7) << correl2[i][j] << '\t';
					// os << setw(7) << correl[i][j] * correl[i][j] << '\t';
				}
				os << '\n';
			}
			os << '\n';
		}

	}

	void PrintCorrel(ostream& os) const	{

		os.precision(3);

		// correlation coefficients
		os << "correlation coefficients\n";
		os << '\n';
		for (int i=0; i<GetDim(); i++)	{
			for (int j=0; j<GetDim(); j++)	{
				os << setw(7) << correl[i][j] << '\t';
			}
			os << '\n';
		}
		os << '\n';

	}

	void PrintPartialR2(ostream& os) const	{

        if (! inverse)  {
            cerr << "error in MeanCovMatrix::PrintPartialR2: inverse not activated\n";
            exit(1);
        }
		os.precision(3);

		// correlation coefficients
		os << "partial correlation coefficients\n";
		os << '\n';
		for (int i=0; i<GetDim(); i++)	{
			for (int j=0; j<GetDim(); j++)	{
				os << setw(7) << partialcorrel2[i][j] << '\t';
				// os << setw(7) << partialcorrel[i][j] * partialcorrel[i][j] << '\t';
			}
			os << '\n';
		}
		os << '\n';

	}

	void PrintPartialCorrel(ostream& os) const	{

        if (! inverse)  {
            cerr << "error in MeanCovMatrix::PrintPartialCorrel: inverse not activated\n";
            exit(1);
        }
		os.precision(3);

		// correlation coefficients
		os << "partial correlation coefficients\n";
		os << '\n';
		for (int i=0; i<GetDim(); i++)	{
			for (int j=0; j<GetDim(); j++)	{
				os << setw(7) << partialcorrel[i][j] << '\t';
			}
			os << '\n';
		}
		os << '\n';

	}

	void PrintPosteriorProbs(ostream& os) const	{

		if (tex)	{
			os << "{\\sc posterior prob.}";
			for (int i=0; i<GetDim(); i++)	{
				os << "&" << i;
			}
			os << "\\\\\n";
			os << "\\hline\n";
			os << '\n';
			for (int i=0; i<GetDim(); i++)	{
				os << i;
				for (int j=0; j<i; j++)	{
					os << "&-";
				}
				for (int j=i; j<GetDim(); j++)	{
					os << " & ";
					os << "$";
					os << ((double) ((int) (100 * pp[i][j]))) / 100;
					if ((pp[i][j] > 0.975) || (pp[i][j] < 0.025))	{
						os << "^*";
					}
					os << "$";
				}
				os << "\\\\\n";
			}
			os << "\\\\\n";
		}
		else	{
			os.precision(2);
			// pp
			os << "posterior probs\n";
			os << '\n';
			for (int i=0; i<GetDim(); i++)	{
				for (int j=0; j<GetDim(); j++)	{
					if (i != j)	{
						os << setw(3) << pp[i][j] << '\t';
					}
					else	{
						os << setw(3) << '-'  << '\t';
					}
				}
				os << '\n';
			}
			os << '\n';
			}
	}

	void PrintPrecisionsPosteriorProbs(ostream& os) const	{

        if (! inverse)  {
            cerr << "error in MeanCovMatrix::PrintPrecisionPosterioProbs: inverse not activated\n";
            exit(1);
        }
		os.precision(2);

		// pp
		os << "posterior probs\n";
		os << '\n';
		for (int i=0; i<GetDim(); i++)	{
			for (int j=0; j<GetDim(); j++)	{
				if (i != j)	{
					os << setw(3) << ppinv[i][j] << '\t';
				}
				else	{
					os << setw(3) << '-'  << '\t';
				}
			}
			os << '\n';
		}
		os << '\n';
	}

    void PrintSlopesNe(ostream& os, int idxdNdS, int idxpiNpiS, int idxpiS, int idxNe) const   {

        int l = (int) (0.025 * size);
        int j = idxNe; 
        os << "relation\tmedian\tCI95min\tCI95max\n";
        os << "logdN/dS~logNe" << '\t' << meanslope[j][idxdNdS] << '\t' << slope[j][idxdNdS][l] << '\t' << slope[j][idxdNdS][size-1-l] << '\n';
        if (idxpiNpiS != -1)    {
            os << "logpiN/piS~logNe" << '\t' << meanslope[j][idxpiNpiS] << '\t' << slope[j][idxpiNpiS][l] << '\t' << slope[j][idxpiNpiS][size-1-l] << '\n';
        }
        os << '\n';
        j = idxpiS;
        os << "logdN/dS~logpiS" << '\t' << meanslope[j][idxdNdS] << '\t' << slope[j][idxdNdS][l] << '\t' << slope[j][idxdNdS][size-1-l] << '\n';
        if (idxpiNpiS != -1)    {
            os << "logpiN/piS~logpiS" << '\t' << meanslope[j][idxpiNpiS] << '\t' << slope[j][idxpiNpiS][l] << '\t' << slope[j][idxpiNpiS][size-1-l] << '\n';
        }
    }

	void PrintSlopes(ostream& os) const	{

		os.precision(3);

		// slopes
		os << "slopes \n";
		os << '\n';
		for (int i=0; i<GetDim(); i++)	{
			for (int j=0; j<GetDim(); j++)	{
				if (i != j)	{
					int l = (int) (0.025 * size);
					// l = 0;
					// os << setw(7) << meanslope[j][i] << " ( " << setw(7) << slope[j][i][l] << " , " << setw(7) << slope[j][i][size-1-l] << " ) " << '\t';
					os << setw(7) << meanslope[j][i] << " : " << sqrt(varslope[j][i]) << " ( " << setw(7) << slope[j][i][l] << " , " << setw(7) << slope[j][i][size-1-l] << " ) " << '\t';
				}
				else	{
					os << setw(7) << '-' << '(' << setw(7) << '-' << " , " << setw(7) << '-' << " ) " << '\t';
				}

			}
			os << '\n';
		}
		os << '\n';
	}

	void ToStream(ostream& os) const {
		if (tex)	{
			os << "\\begin{tabular}{";
			for (int i=0; i<GetDim() + 1; i++)	{
				os << "r";
			}
			os << "}\n";
			os << '\n';

			// PrintCovariances(os);
			PrintR2(os);
			// PrintPosteriorProbs(os);

			os << "\\end{tabular}\n";
		}
		else	{
			PrintCovariances(os);
			PrintCorrel(os);
			PrintPosteriorProbs(os);
            if (inverse)    {
                PrintPrecisions(os);
                PrintPartialCorrel(os);
                PrintPrecisionsPosteriorProbs(os);
            }
			// PrintSlopes(os);
		}
	}

	friend ostream& operator<<(ostream& os, const MeanCovMatrix& m)	{
		m.ToStream(os);
		return os;
	}

	protected:

	bool tex;
	int dim;
	unsigned int size;
    bool inverse;

	// arrays of values
	vector<double>** val;
	vector<double>** invval;
	vector<double>** slope;

	double** mean;
	double** meanslope;
	double** varslope;
	double** var;
	double** pp;
	double** correl;
	double** partialcorrel;
	double** correl2;
	double** partialcorrel2;
	double** meaninv;
	double** varinv;
	double** ppinv;

	double* meanpropvar;
};

class ReducedMeanCovMatrix : public MeanCovMatrix	{

	public:

	ReducedMeanCovMatrix(MeanCovMatrix* infrom, int* inmultindexarray, int dim) : MeanCovMatrix(dim)	{

		from = infrom;
		int fulldim = from->GetDim();
		int checkdim = 0;
		multindexarray = new int[fulldim];
		multindextable = new int[fulldim];
		for (int i=0; i<fulldim; i++)	{
			multindexarray[i] = inmultindexarray[i];
			if (multindexarray[i])	{
				multindextable[checkdim] = i;
				checkdim++;
			}
		}
		if (dim != checkdim)	{
			cerr << "error in reduced cov matrix : non matching reduced dimension\n";
			exit(1);
		}

		double** tmp = new double*[fulldim];
		for (int i=0; i<fulldim; i++)	{
			tmp[i] = new double[fulldim];
		}
		double** u = new double*[dim];
		double** invu = new double*[dim];
		for (int i=0; i<dim; i++)	{
			u[i] = new double[dim];
			invu[i] = new double[dim];
		}
		for (unsigned int k=0; k<from->GetSize(); k++)	{
			for (int i=0; i<fulldim; i++)	{
				for (int j=0; j<fulldim;j++)	{
					tmp[i][j] = from->GetVal(i,j,k);
				}
			}
			for (int l=0; l<fulldim; l++)	{
				if (!multindexarray[l])	{
					for (int i=0; i<fulldim; i++)	{
						if (i != l)	{
							for (int j=0; j<fulldim;j++)	{
								if (j != l)	{
									tmp[i][j] -= tmp[i][l] * tmp[l][j] / tmp[l][l];
								}
							}
						}
					}
				}
			}
			for (int i=0; i<dim; i++)	{
				int ii = multindextable[i];
				for (int j=0; j<dim; j++)	{
					int jj = multindextable[j];
					val[i][j].push_back(tmp[ii][jj]);
					u[i][j] = tmp[ii][jj];
				}
			}
			LinAlg::Gauss(u,dim,invu);
			for (int i=0; i<dim; i++)	{
				for (int j=0; j<dim; j++)	{
					invval[i][j].push_back(invu[i][j]);
				}
			}
		}

		for (int i=0; i<dim; i++)	{
			delete[] invu[i];
			delete[] u[i];
		}
		delete[] invu;
		delete[] u;

		for (int i=0; i<fulldim; i++)	{
			delete[] tmp[i];
		}
		delete[] tmp;

		size = from->GetSize();
		Normalize();
	}

	~ReducedMeanCovMatrix()	{
		delete[] multindexarray;
		delete[] multindextable;
	}


	private:

	int* multindexarray;
	int* multindextable;
	MeanCovMatrix* from;

};

class PartialMeanCovMatrix : public MeanCovMatrix	{

	public:

	PartialMeanCovMatrix(MeanCovMatrix* infrom, int* inmultindexarray, int dim) : MeanCovMatrix(dim)	{

		from = infrom;
		int fulldim = from->GetDim();
		int checkdim = 0;
		cerr << "in constructor\n";
		cerr << dim << '\t' << fulldim << '\n';
		multindexarray = new int[fulldim];
		multindextable = new int[fulldim];
		for (int i=0; i<fulldim; i++)	{
			multindexarray[i] = inmultindexarray[i];
			if (multindexarray[i])	{
				multindextable[checkdim] = i;
				checkdim++;
			}
		}
		if (dim != checkdim)	{
			cerr << "error in reduced cov matrix : non matching reduced dimension\n";
			cerr << dim << '\t' << checkdim << '\n';
			exit(1);
		}

		double** u = new double*[dim];
		double** invu = new double*[dim];
		for (int i=0; i<dim; i++)	{
			u[i] = new double[dim];
			invu[i] = new double[dim];
		}
		for (unsigned int k=0; k<from->GetSize(); k++)	{
			for (int i=0; i<dim; i++)	{
				int ii = multindextable[i];
				for (int j=0; j<dim; j++)	{
					int jj = multindextable[j];
					double temp = from->GetVal(ii,jj,k);
					val[i][j].push_back(temp);
					u[i][j] = temp;
					// invval[i][j].push_back(from->GetInvVal(ii,jj,k));
				}
			}
			LinAlg::Gauss(u,dim,invu);
			for (int i=0; i<dim; i++)	{
				for (int j=0; j<dim; j++)	{
					invval[i][j].push_back(invu[i][j]);
				}
			}
		}

		for (int i=0; i<dim; i++)	{
			delete[] invu[i];
			delete[] u[i];
		}
		delete[] invu;
		delete[] u;

		size = from->GetSize();
		Normalize();
	}

	~PartialMeanCovMatrix()	{
		delete[] multindexarray;
		delete[] multindextable;
	}


	private:

	int* multindexarray;
	int* multindextable;
	MeanCovMatrix* from;

};

class SumConstrainedMeanCovMatrix : public MeanCovMatrix	{

	protected:

	int reduceddim;
	int Ncomp;
	int* Kcomp;
	SumConstrainedMapping** mapping;
	double** extended;
	double** invextended;
	double** aux;
	double** basis;
	int * indexmap;

	public:

	SumConstrainedMeanCovMatrix(int inreduceddim, int inNcomp, int* inKcomp, SumConstrainedMapping** inmapping) : MeanCovMatrix(inreduceddim + inNcomp)	{

		reduceddim = inreduceddim;
		Ncomp = inNcomp;
		Kcomp = inKcomp;
		mapping = inmapping;

		// make global change of variables
		indexmap = new int[reduceddim];

		basis = new double*[GetDim()];
		extended = new double*[GetDim()];
		invextended = new double*[GetDim()];
		aux = new double*[GetDim()];
		for (int i=0; i<GetDim(); i++)	{
			basis[i] = new double[GetDim()];
			extended[i] = new double[GetDim()];
			invextended[i] = new double[GetDim()];
			aux[i] = new double[GetDim()];
		}
		MakeBasis();
	}

	~SumConstrainedMeanCovMatrix()	{
		for (int i=0; i<GetDim(); i++)	{
			delete[] basis[i];
			delete[] extended[i];
			delete[] invextended[i];
			delete[] aux[i];
		}
		delete[] basis;
		delete[] extended;
		delete[] invextended;
		delete[] aux;

		delete[] indexmap;
	}

	void MakeBasis()	{

		int i = 0;
		int k = 0;
		int j = 0;
		while (i < reduceddim)	{
			if ((k<Ncomp) && (i == Kcomp[k]))	{
				j++;
				k++;
			}
			indexmap[i] = j;
			i++;
			j++;
		}
		if (k != Ncomp)	{
			cerr << "error in sum constrained mean cov matrix; make basis: non matching k\n";
			exit(1);
		}
		if (j != GetDim())	{
			cerr << "error in sum constrained mean cov matrix; make basis: non matching j\n";
			exit(1);
		}

		for (int i=0; i<GetDim(); i++)	{
			for (int j=0; j<GetDim(); j++)	{
				basis[i][j] = 0;
				extended[i][j] = 0;
				invextended[i][j] = 0;
			}
		}
		for (int i=0; i<GetDim(); i++)	{
			basis[i][i] = 1;
		}

		for (int k=0; k<Ncomp; k++)	{
			int a = indexmap[Kcomp[k]]-1;
			for (int i=0; i<mapping[k]->GetDim(); i++)	{
				for (int j=0; j<mapping[k]->GetDim(); j++)	{
					basis[i+a][j+a] = mapping[k]->base[j][i];
				}
			}
		}
	}

	void Add(CovMatrix* sample)	{
		sample->Diagonalise();
		double** m = sample->GetMatrix();
		double** invm = sample->GetInvMatrix();
		if (sample->GetDim() != reduceddim)	{
			cerr << "error in sum constrained cov matrix : non matching dim\n";
			exit(1);
		}
		for (int i=0; i<reduceddim; i++)	{
			for (int j=0; j<reduceddim; j++)	{
				extended[indexmap[i]][indexmap[j]] = m[i][j];
				invextended[indexmap[i]][indexmap[j]] = invm[i][j];
			}
		}

		for (int i=0; i<GetDim(); i++)	{
			for (int j=0; j<GetDim(); j++)	{
				double tmp = 0;
				for (int k=0; k<GetDim(); k++)	{
					tmp += extended[i][k] * basis[j][k];
				}
				aux[i][j] = tmp;
			}
		}
		for (int i=0; i<GetDim(); i++)	{
			for (int j=0; j<GetDim(); j++)	{
				double tmp = 0;
				for (int k=0; k<GetDim(); k++)	{
					tmp += basis[i][k] * aux[k][j];
				}
				val[i][j].push_back(tmp);
			}
		}

		for (int i=0; i<GetDim(); i++)	{
			for (int j=0; j<GetDim(); j++)	{
				double tmp = 0;
				for (int k=0; k<GetDim(); k++)	{
					tmp += invextended[i][k] * basis[j][k];
				}
				aux[i][j] = tmp;
			}
		}
		for (int i=0; i<GetDim(); i++)	{
			for (int j=0; j<GetDim(); j++)	{
				double tmp = 0;
				for (int k=0; k<GetDim(); k++)	{
					tmp += basis[i][k] * aux[k][j];
				}
				invval[i][j].push_back(tmp);
			}
		}
		size++;
	}
};


#endif

