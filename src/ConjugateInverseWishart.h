
#ifndef CONJUGATEWISHART_H
#define CONJUGATEWISHART_H

#include "Conjugate.h"
#include "InverseWishartMatrix.h"




// the mean parameter of an invert wishart distribution
class MultiNormalSemiConjugate : virtual public SemiConjugatePrior<CovMatrix> {

	public:

	MultiNormalSemiConjugate(int indim) :
		Rvar<CovMatrix>(CovMatrix(indim)),
		scalestat(indim), bkscalestat(indim)
	{}

	~MultiNormalSemiConjugate(){}

	// Must be initialize with
	void ResetSufficientStatistic()	{
		shapestat = 0;
		for( int i=0; i<dim; i++){
			for( int j=0; j<dim; j++){
				scalestat[i][j] = 0;
			}
		}
	}


	void SaveSufficientStatistic()	{
		bkshapestat = shapestat;
		for( int i=0; i<dim; i++){
			bkscalestat[i][i] = scalestat[i][i];
			for( int j=0; j<i; j++){
				bkscalestat[i][j] = scalestat[i][j];
			}
		}
	}

	void RestoreSufficientStatistic()	{
		bkshapestat = shapestat;
		for( int i=0; i<dim; i++){
			scalestat[i][i] = bkscalestat[i][i];
			for( int j=0; j<i; j++){
				scalestat[i][j] = bkscalestat[i][j];
				scalestat[j][i] = bkscalestat[i][j];
			}
		}
	}

	void AddToShape()	{
		shapestat ++;
	}

	void AddToScale(const double* in)	{
		for( int i=0; i<dim; i++){
			scalestat[i][i] += in[i] * in[i];
			for( int j=0; j<i; j++){
				scalestat[i][j] += in[i] * in[j];
				scalestat[j][i] = scalestat[i][j];
			}
		}
	}

	void AddToShape(int N)	{
		shapestat +=N;
	}

	void AddToScale(double** in)	{
		for( int i=0; i<dim; i++){
			for(int j=0; j<dim; j++) {
				 scalestat[i][j] += in[i][j];
			}
		}
	}

	double SuffStatLogProb()	{
		double t, trace = 0;
		for( int i=0; i<dim; i++){
			t = 0;
			for( int j=0; j<dim; j++){
				t += scalestat[i][j] * GetInvMatrix()[j][i];
			}
			trace += t;
		}
		return -0.5 * (shapestat * GetLogDeterminant() + trace);
	}

	protected:

	int shapestat;
	CovMatrix scalestat;
	int bkshapestat;
	CovMatrix bkscalestat;
};


// as an example, we can give this mean parameter an inverse wishart  distribution
class SemiConjugateInverseWishart : virtual public InverseWishartMatrix, virtual public MultiNormalSemiConjugate{

	public:

	SemiConjugateInverseWishart(SigmaZero* inDiagonalMatrix, int inP) :
		Rvar<CovMatrix>(inDiagonalMatrix->GetDim()),
		InverseWishartMatrix(inDiagonalMatrix, inP),
		MultiNormalSemiConjugate(inDiagonalMatrix->GetDim())
	{}

	double logProb()	{
		if (isActive())	{
			return InverseWishartMatrix::logProb() + SuffStatLogProb();
		}
		else { return InverseWishartMatrix::logProb(); }
	}

	void localRestore()	{
		RestoreSufficientStatistic();
		InverseWishartMatrix::localRestore();
	}

	void localCorrupt(bool bk)	{
		if (bk)	{
			SaveSufficientStatistic();
		}
		InverseWishartMatrix::localCorrupt(bk);
	}
};



// but in fact, we have analytical conjugacy here, so let us implement it
 class ConjugateInverseWishart : public ConjugatePrior<CovMatrix>, public SemiConjugateInverseWishart {

	public:

	ConjugateInverseWishart(SigmaZero* inDiagonalMatrix, int inP) :
		Rvar<CovMatrix>(inDiagonalMatrix->GetDim()),
		SemiConjugatePrior<CovMatrix>(),
		InverseWishartMatrix(inDiagonalMatrix, inP),
		MultiNormalSemiConjugate(inDiagonalMatrix->GetDim()),
		ConjugatePrior<CovMatrix>(),
		SemiConjugateInverseWishart(inDiagonalMatrix, inP)
	{}

	void GibbsResample()	{

		for(int i=0; i<GetDim(); i++){
			scalestat[i][i] += diagonalMatrix->val(i);
		}
		P += shapestat;
		scalestat.Diagonalise();
		drawSample(&scalestat);
		for(int i=0; i<GetDim(); i++){
			scalestat[i][i] -= diagonalMatrix->val(i);
		}
		P -= shapestat;
		// just be careful about mathematical errors
		if (!isPositiveDefine())	{
			cout << " drawForWishart de SigmaZero renvoie un sigma non défini positif\n";
			exit(1);
		}

	}

	double logProb()	{
		if (isActive())	{

			// if in integrated (conjugate) mode, return the integrated log probability
			if (isIntegrated())	{

				for(int i=0; i<GetDim(); i++){
					scalestat[i][i] += diagonalMatrix->val(i);
				}
				CovMatrix* postscale = new CovMatrix(scalestat);
				int postshape = P + shapestat;
				double l = diagonalMatrix->GetLogDeterminant() * P * 0.5 - postscale->GetLogDeterminant() * postshape * 0.5;
				delete postscale;
				for(int i=0; i<GetDim(); i++){
					scalestat[i][i] -= diagonalMatrix->val(i);
				}
				return l;
			}

			// if in active but not integrated mode (i.e. semi-conjugate mode), return the semi-conjugate log probabilty
			return SemiConjugateInverseWishart::logProb();
		}
		// if in inactive mode, return the plain log probability, as a Gamma
		return InverseWishartMatrix::logProb();
	}

	double DiagonalLogPrior() 	{ // double min, double max)	{
		double ret = GetDim() * Random::logGamma(0.5 * (P + GetDim() - 1)) + 0.5 * GetDim() * (P + GetDim() -1) * log(2.0) - 0.5 * P*GetDim() * log(2.0) - Random::logMultivariateGamma(0.5 * P,GetDim());
		for (int i=0; i<GetDim(); i++)	{
			ret -= 0.5 * (GetDim() - 1) * log(diagonalMatrix->val(i));
		}
		/*
		if (GetDim() != 1)	{
			double tmp = 2.0 / (GetDim() - 1) * (exp(-0.5*(GetDim()-1) * log(min)) - exp(-0.5*(GetDim()-1) * log(max))) / (log(max) - log(min));
			ret += GetDim() * log(tmp);
		}
		*/
		return ret;
	}

	double DiagonalLogPrior(double min, double max)	{
		double ret = GetDim() * Random::logGamma(0.5 * (P + GetDim() - 1)) + 0.5 * GetDim() * (P + GetDim() -1) * log(2.0) - 0.5 * P*GetDim() * log(2.0) - Random::logMultivariateGamma(0.5 * P,GetDim());
		if (GetDim() != 1)	{
			double tmp = 2.0 / (GetDim() - 1) * (exp(-0.5*(GetDim()-1) * log(min)) - exp(-0.5*(GetDim()-1) * log(max))) / (log(max) - log(min));
			ret += GetDim() * log(tmp);
		}
		return ret;
	}

	double DiagonalLogPosterior()	{
		double ret = GetDim() * Random::logGamma(0.5 * (P + shapestat + GetDim() - 1)) + 0.5 * GetDim() * (P + +shapestat + GetDim() -1) * log(2.0) - 0.5 * (P + shapestat) *GetDim() * log(2.0) - Random::logMultivariateGamma(0.5 * (P + shapestat),GetDim());
		for (int i=0; i<GetDim(); i++)	{
			ret -= 0.5 * (GetDim() - 1) * log(diagonalMatrix->val(i) + scalestat[i][i]);
			// ret += 0.5 * (P + shapestat - GetDim() + 1) * log(1.0 + scalestat[i][i]);
			// ret += 0.5 * (P + shapestat - GetDim() + 1) * log(diagonalMatrix->val(i) + scalestat[i][i]);
		}
		return ret;
	}

};

// this mean parameter is (semi-)conjugate to a Wishart sampling distribution
class ConjugateMultiNormal : public ConjugateSampling<RealVector>, public MultiNormal	{

	public:

	ConjugateMultiNormal(MultiNormalSemiConjugate* insigma, Var<RealVector>* inup = 0, Var<PosReal>* intime = 0, Var<PosReal>* inscale = 0, Var<RealVector>* indrift = 0, Var<PosReal>* indriftphi = 0, Var<PosReal>* indate=0, Var<RealVector>* indrift2 = 0, Var<PosReal>* indriftphi2 = 0, Var<PosReal>* inagescale = 0, double inkt = 0, GenericTimeLine* timeline = 0, Var<Real>* alpha = 0) :

		Rvar<RealVector>(),
		ConjugateSampling<RealVector>(),
		MultiNormal(insigma, inup, intime, inscale, indrift, indriftphi, indate, indrift2, indriftphi2, inagescale, inkt, timeline, alpha)
	{
		conjugate_up.insert(insigma);
		contrast = 0;
	}

	ConjugateMultiNormal(int unused, MultiNormalSemiConjugate* insigma, Var<RealVector>* inrootmean, Var<PosRealVector>* inrootvar)	 : MultiNormal(unused,insigma,inrootmean,inrootvar)	{

		conjugate_up.insert(insigma);
		contrast = 0;
	}

	ConjugateMultiNormal(MultiNormalSemiConjugate* insigma, Var<RealVector>* inup, Var<PosReal>* intime, Var<PosReal>* indate, Var<PosReal>* inagescale, GlobalScalingFunction* inscalefunction) : MultiNormal(insigma,inup,intime,indate,inagescale,inscalefunction) {

		conjugate_up.insert(insigma);
		contrast = 0;
	}

	~ConjugateMultiNormal() {}

	void AddSufficientStatistic(SemiConjPrior* parent)	{
		if (up)	{
			MultiNormalSemiConjugate* prior = dynamic_cast<MultiNormalSemiConjugate*>(parent);
			if (! prior)	{
				cout << "cast error in ConjugateMultiNormal::AddSuffStat\n";
				exit(1);
			}
			prior->AddToShape();
			ComputeContrast();
			prior->AddToScale(contrast);
		}
	}

	// overriding simple behavior, because this simple behavior assumes a specified covmatrix
	double ProposeMove(double tuning)	{
		// return SimpleProposeMove(tuning);
		if (isUpActive())	{
			return SimpleProposeMove(tuning);
		}
		return MultiNormal::ProposeMove(tuning);
	}

	void ComputeContrast()	{
		if (! contrast)	{
			contrast = new double[GetDim()];
		}
		double tt = time->val();
		if (scalefunction)	{
			double t2 = date->val() * agescale->val();
			double t1 = (date->val() + time->val()) * agescale->val();
			double scalefactor = scalefunction->GetScalingFactor(t1,t2);
			tt *= scalefactor;
		}
		else if (scale)	{
			tt *= scale->val();
		}
		double roott = sqrt(tt);
		if (drift)	{
			if (driftphi)	{
				double u = time->val();
				if (driftphi->val() > 1e-10)	{
					u = exp(-driftphi->val() * date->val()) * (1 - exp(-driftphi->val() * time->val())) / driftphi->val();
				}
				for (int i=0; i<GetDim(); i++)	{
					contrast[i] = ((*this)[i] - u * (*drift)[i] - (*up)[i]) / roott;
				}
			}
			else	{
				for (int i=0; i<GetDim(); i++)	{
					contrast[i] = ((*this)[i] - tt * (*drift)[i] - (*up)[i]) / roott;
				}
			}
		}
		else	{
			for (int i=0; i<GetDim(); i++)	{
				contrast[i] = ((*this)[i] - (*up)[i]) / roott;
			}
		}
	}

	double* GetContrast()	{
		return contrast;
	}

	private:
	double* contrast;

};


// this mean parameter is (semi-)conjugate to a Wishart sampling distribution
class ConjugateMultivariateNormal : public ConjugateSampling<RealVector>, public MultivariateNormal	{

	public:

	ConjugateMultivariateNormal(MultiNormalSemiConjugate* insigma) :
		Rvar<RealVector>(),
		ConjugateSampling<RealVector>(),
		MultivariateNormal(insigma)
	{
		conjugate_up.insert(insigma);
	}

	~ConjugateMultivariateNormal() {}

	void AddSufficientStatistic(SemiConjPrior* parent)	{
		MultiNormalSemiConjugate* prior = dynamic_cast<MultiNormalSemiConjugate*>(parent);
		if (! prior)	{
			cout << "cast error in ConjugateMultiNormal::AddSuffStat\n";
			exit(1);
		}
		prior->AddToShape();
		prior->AddToScale(GetArray());
	}
};


#endif

