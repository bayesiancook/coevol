
#ifndef BROWNIANBRIDGE_H
#define	BROWNIANBRIDGE_H

#include "Random.h"
#include "CovMatrix.h"

enum Segmentation {SEGM_REGULAR, SEGM_ABSOLUTE};

class BrownianBridge {

public:
	BrownianBridge();
	BrownianBridge(double ininitAge, double infinalAge, CovMatrix variance, int innSegments, double* insegmentLength);
	BrownianBridge(const BrownianBridge& from);
	virtual ~BrownianBridge();

	void operator=(const BrownianBridge& from);
	double* operator[](int i){return bridge[i];}

	virtual void buildSegments();
	virtual void buildRegularSegments();
	virtual void buildAbsoluteSegments();

	void setAges(double ininitAge, double infinalAge, bool recomputeNSegments);
	double getLength() {return length;}
	double getSegmentLength(int i) {return segmentLength[i];}
	double* getSegmentLength() {return segmentLength;}

	void setUnitVariance(CovMatrix unitVariance);
	CovMatrix& getUnitVariance() {return unitVariance;};
	double** getBridge() {return bridge;}
	int getNSegments() {return nSegments;}
	int getDim() {return dim;}
	double getInitAge() {return initAge;}
	double getFinalAge() {return finalAge;}

	double getValue(int node, int param) {return bridge[node][param];}
	double operator()(int node, int param) {return getValue(node,param);}

	double getMaxValue();
	double logProbNormal(double* x, double* mean, double slength, double** prec, double logPrecDet);
	double getLogProb();
	double getLogProb(CovMatrix& covar);
	double getLogProbSubBridge(int beg, int end);

	void generateBridge();
	void generateSubBridge(int beg, int end);
	void setBridge(double** inBridge) {bridge = inBridge;}
	void showBridge();

	virtual double ProposeMove(double tuning);
	virtual double ProposeMove(double tuning, CovMatrix& othervar);
	virtual double SimpleProposeMove(double tuning);
	double initAgeMove(double newAge, double* initVal);
	double finalAgeMove(double newAge, double* finalVal);

	virtual double getIntegral(int param, double initValue, double finalValue);

	void Shift(BrownianBridge* bridge, double tuning);

	void LeftMultiply(double** P);

	friend ostream& operator<<(ostream& os, const BrownianBridge& r)  {



		os << r.nSegments << '\t' << r.initAge << '\t' << r.finalAge << '\t' << r.dim;

		for (int i=0; i<r.nSegments; i++) {
			os  << '\t' << r.segmentLength[i];
		}

		for (int i=0; i<r.nSegments+1; i++)
			for (int j=0; j<r.dim; j++)
				os << '\t' << r.bridge[i][j];

		return os;
	}
	friend istream& operator>>(istream& is, BrownianBridge& r)  {


		if (r.bridge) {		 
			r.removeBridge(r.bridge);
		}
		if(r.segmentLength)
			delete[] r.segmentLength;

		is >> r.nSegments;
		is >> r.initAge;
		is >> r.finalAge;
		is >> r.dim;

		r.length = r.initAge-r.finalAge;

		r.segmentLength = new double[r.nSegments];

		double s=0;

		for (int i=0; i<r.nSegments; i++) {
			is >> r.segmentLength[i];
			s+=r.segmentLength[i];
		}

		r.bridge = new double*[r.nSegments+1];

		for (int i=0; i<r.nSegments+1; i++)	 {
			r.bridge[i] = new double[r.dim];
			for(int j=0; j<r.dim; j++)
				is >> r.bridge[i][j];
		}

		return is;
	}


	static void setNTreeSegments(int innTreeSegments) {nTreeSegments = innTreeSegments;}
	static void setSegmentation(Segmentation insegmentation) {segmentation = insegmentation;}
	static Segmentation getSegmentation() {return segmentation;}


	void findPointBefore(double time, double &timeBefore, double* &valueBefore);
	void findPointAfter(double time, double &timeAfter, double* &valueAfter);
	int findPointBefore(double time);
	int findPointAfter(double time);

	void addSlope(double* initValue, double* finalValue);
	void removeSlope();

	void checkLength(const double _initAge, const double _finalAge, const double _length) {

		if(abs(initAge-_initAge)>1e-5) {
			cerr << "Error in brownian bridge : initAge different from the chronogram value" << endl;
			cerr << initAge << " " << _initAge << endl;

		}

		if(abs(finalAge-_finalAge)>1e-5) {
			cerr << "Error in brownian bridge : finalAge different from the chronogram value" << endl;
			cerr << finalAge << " " << _finalAge << endl;

		}

		if(abs(initAge-_initAge)>1e-5) {
			cerr << "Error in brownian bridge : length different from the difference between the extremities ages" << endl;
			exit(0);
		}

		if(abs(length-(initAge-finalAge))>1e-5) {
			cerr << "Error in brownian bridge : length different from the difference between the extremities ages" << endl;
			exit(0);
		}
		if(abs(length - _length) > 1e-5) {
			cerr << "Error in brownian bridge : length different from the chronogram value" << endl;
			cerr << length << " " << _length << endl;
			exit(0);
		}
		double s = 0;
		for(int i=0; i<nSegments; i++) {
			s+=segmentLength[i];
		}
		if(abs(s-length)>1e-5) {
			cerr << "Error in brownian bridge : branch length different from the segment length sum" << endl;
			cerr << s << " " << length;
			exit(0);
		}
	}

protected:

	double** createBridge(int N = -1);
	void removeBridge(double**b, int N = -1);


	static int nTreeSegments;  //number of segments in the whole tree depth
	static Segmentation segmentation;   //segmentation type
	int nSegments;
	int dim;
	double length, initAge, finalAge;
	CovMatrix unitVariance;
	double** bridge;
	double* segmentLength;


};

#endif	/* BROWNIANBRIDGE_H */

