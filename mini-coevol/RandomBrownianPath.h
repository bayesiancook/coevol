
#ifndef RANDOMBROWNIANPATH_H
#define	RANDOMBROWNIANPATH_H


#include "BrownianBridge.h"
#include "Var.h"

#include "ConjugateInverseWishart.h"

class RandomBrownianPath : virtual public Rvar<BrownianBridge> {

	 public:
		RandomBrownianPath(Var<PosReal> *inup, Var<PosReal> *indown, Var<CovMatrix> *insigma) :
		 Up(inup), Down(indown), Sigma(insigma) {

			setAges(inup->val(), indown->val(), true);

			Register(Up);
			Register(Down);
			Register(Sigma);
			SetName("Random multi-variate brownian path");

			Sample();
		}

		virtual void drawSample() {//Implements MCMC::DrawSample
			updateBrownianBridgeParameters();
			generateBridge();
		} 


		virtual double logProb() {  //Implements DAGNode::LogProb
			updateBrownianBridgeParameters();
			return BrownianBridge::getLogProb();
		}

		double ProposeMoveCov(CovMatrix& cov, double tuning=1)	{
			return BrownianBridge::ProposeMove(tuning,cov);
		}

		void updateBrownianBridgeParameters() { 
			setUnitVariance(Sigma->val());

			if(BrownianBridge::getSegmentation() == SEGM_REGULAR) {
				setAges(Up->val(), Down->val(), false);
			}
			/*
			// NICO
			else	{
				setAges(Up->val(), Down->val(), true);
			}
			*/
		}


		//Accessors
		Var<PosReal> *getUp() {return Up;}
		Var<PosReal> *getDown() {return Down;}
		virtual double getIntegral(int param, double initValue, double finalValue) {
			updateBrownianBridgeParameters();
			return BrownianBridge::getIntegral(param, initValue, finalValue);
		}

		virtual ~RandomBrownianPath() {
		}

		protected:

		Var<PosReal> *Up;
		Var<PosReal> *Down;
		Var<CovMatrix> *Sigma;


};

class ConjugateRandomBrownianPath : public virtual ConjugateSampling<BrownianBridge>, public virtual RandomBrownianPath	{

	public:

	ConjugateRandomBrownianPath(Var<PosReal> *inup, Var<PosReal> *indown, ConjugateInverseWishart *insigma)	:
		RandomBrownianPath(inup,indown,insigma)	{

		conjugate_up.insert(insigma);

		contrast = new double*[dim];
		for (int i=0; i<dim; i++)	{
			contrast[i] =  new double[dim];
		}
	}

	~ConjugateRandomBrownianPath() {

		for (int i=0; i<dim; i++)	{
			delete[] contrast[i];
		}
		delete[] contrast;
	}

	void AddSufficientStatistic(SemiConjPrior* parent)	{

		MultiNormalSemiConjugate* prior = dynamic_cast<MultiNormalSemiConjugate*>(parent);
		if (! prior)	{
			cout << "cast error in ConjugateMultiNormal::AddSuffStat\n";
			exit(1);
		}
		prior->AddToShape(nSegments-1);
		ComputeContrast();
		prior->AddToScale(contrast);
	}

	// overriding default behavior, which assumes a specified covmatrix
	double ProposeMove(double tuning)	{

		if (isUpActive())	{
			return SimpleProposeMove(tuning);
		}
		return BrownianBridge::ProposeMove(tuning);
	}

	void ComputeContrast()	{

		for(int i=0; i<dim; i++)	{
			for(int j=0; j<dim; j++)	{
				contrast[i][j] = 0;
			}
		}

		int N = getNSegments();

		double* y = new double[dim];

		for(int k=0; k<N; k++) {
			double segmentLength = getSegmentLength(k);
			for(int i=0; i<dim; i++)	{
				y[i] = getValue(k+1, i) - getValue(k, i);
			}
			for(int i=0; i<dim; i++)	{
				for(int j=0; j<dim; j++)	{
					contrast[i][j] += y[i]*y[j]/segmentLength;
				}
			}
		}

		/*
		for(int i=0; i<dim; i++)	{
			for(int j=0; j<dim; j++)	{
				if (isinf(contrast[i][j]))	{
					cerr << "inf\n";
					exit(1);
				}
				if (isnan(contrast[i][j]))	{
					cerr << "nan\n";
					exit(1);
				}
			}
		}
		*/

		delete[] y;
	}

	double** GetContrast()	{
		return contrast;
	}

	private:
	double** contrast;

};

#endif	/* RANDOMBROWNIANPATH_H */

