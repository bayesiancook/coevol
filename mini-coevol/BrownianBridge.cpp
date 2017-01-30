#include "BrownianBridge.h"


int BrownianBridge::nTreeSegments = 0;
Segmentation BrownianBridge::segmentation = SEGM_ABSOLUTE;

BrownianBridge::BrownianBridge() {
	bridge = 0;
	segmentLength = 0;
}


BrownianBridge::BrownianBridge(double ininitAge, double infinalAge, CovMatrix inunitVariance, int innSegments, double* insegmentLength) {
	 bridge = 0;
	segmentLength = 0;

	if(nTreeSegments==0) {
		cerr << "Error : BrownianBridge::nTreeSegments has not been initialized" << endl;
		exit(0);
	}
  
	initAge = ininitAge;
	finalAge = infinalAge;
	length = initAge-finalAge;

	nSegments = innSegments;
	segmentLength = new double[nSegments];
	for(int i=0; i<nSegments; i++)
		segmentLength[i] = insegmentLength[i];
	unitVariance = inunitVariance;
	dim = unitVariance.GetDim();
	bridge = 0;

}

BrownianBridge::BrownianBridge(const BrownianBridge& from) {
	nSegments = from.nSegments;
	initAge = from.initAge;
	finalAge = from.finalAge;
	length = from.length;
	unitVariance = from.unitVariance;
	dim = from.dim;
	bridge = createBridge();
	for(int i=0; i<nSegments+1; i++) {
		for(int j=0; j<dim; j++)
			bridge[i][j] = from.bridge[i][j];
	}
	segmentLength = new double[nSegments];
	for(int i=0; i<nSegments; i++)
		segmentLength[i] = from.segmentLength[i];

}

void BrownianBridge::operator=(const BrownianBridge& from) {

	removeBridge(bridge);
	if(segmentLength)
		delete[] segmentLength;

	nSegments = from.nSegments;
	initAge = from.initAge;
	finalAge = from.finalAge;
	length = from.length;
	unitVariance = from.unitVariance;
	dim = from.dim;
	bridge = createBridge();
	for(int i=0; i<nSegments+1; i++) {
		for(int j=0; j<dim; j++)
			bridge[i][j] = from.bridge[i][j];
	}
	segmentLength = new double[nSegments];
	for(int i=0; i<nSegments; i++)
		segmentLength[i] = from.segmentLength[i];


}


BrownianBridge::~BrownianBridge()
{
	removeBridge(bridge);
	if(segmentLength) {
		delete[] segmentLength;
	}
}


void BrownianBridge::buildSegments() {

	length = initAge-finalAge;

	if(segmentation == SEGM_REGULAR)
		buildRegularSegments();
	else
		buildAbsoluteSegments();

}

void BrownianBridge::buildRegularSegments() {

	nSegments = (int) (nTreeSegments*length);

	if(segmentLength)
		delete[] segmentLength;
	if(nSegments<3) {
		nSegments=3;
	}
	segmentLength = new double[nSegments];
	for(int i=0; i<nSegments; i++) {
		segmentLength[i] = length/nSegments;
	}

}
void BrownianBridge::buildAbsoluteSegments() {

	double firstPoint = ((double) floor(nTreeSegments * initAge)) / nTreeSegments; //The first point after the initAge
	double lastPoint = ((double) ceil(nTreeSegments * finalAge)) / nTreeSegments;   //The last point before the finalAge

  
	double d = (firstPoint-lastPoint)*nTreeSegments;
	nSegments = (int) (d-floor(d)<0.1 ? floor(d) : ceil(d));	//To avoid numerical problems

	if(initAge-firstPoint > 10e-10)
		nSegments++;	//Adding a first, shorter segment between the initAge and the first point
	if(lastPoint-finalAge > 10e-10)
		nSegments++;	//Adding a last, shorter segment between the last point and the final age

	if (! nSegments)	{
		nSegments++;
	}

	if(segmentLength)
		delete[] segmentLength;

	/*
	if(nSegments==0) {
		cerr << "Error : branch with 0 segments" << endl;
		cerr << "Init age : " << initAge << " " << finalAge;
		exit(0);
	}
	*/

	segmentLength = new double[nSegments];

	if(nSegments>1) {

		for(int i=0; i<nSegments; i++) {
			segmentLength[i] = 1.0/nTreeSegments;
		}
		if(initAge-firstPoint > 10e-10)
			segmentLength[0] = initAge - firstPoint;
		if(lastPoint-finalAge > 10e-10)
			segmentLength[nSegments-1] = lastPoint - finalAge;
	}
	else {
		segmentLength[0] = length;
	}


}

void BrownianBridge::setUnitVariance(CovMatrix inunitVariance) {
	unitVariance = inunitVariance;
	dim = unitVariance.GetDim();
}


void BrownianBridge::setAges(double ininitAge, double infinalAge, bool recomputeNSegments) {
	initAge = ininitAge;
	finalAge = infinalAge;
	if(recomputeNSegments) {
		if(nTreeSegments==0) {
			cerr << "Error : BrownianBridge::nTreeSegments has not been initialized" << endl;
			exit(0);
		}
		removeBridge(bridge);
		buildSegments();
	}
	else {			  //Just a dilatation of the bridge
		double prop = (initAge-finalAge)/length;
		length*=prop;
		for(int i=0; i<nSegments; i++) {
			segmentLength[i]*=prop;
		}
	  
	}

	if(length <1e-12) {
		length = 1e-8;
		/*
		cerr << "Error in brownian bridge : null length" << endl;
		exit(0);
		*/
	}

}

double BrownianBridge::getMaxValue() {
	double value = 0;
	for(int i=0; i<nSegments; i++) {
		double mean = 0;
		for(int j=0; j<dim; j++)
			mean+=bridge[i][j];
		mean/=dim;
		value = value < mean ? mean : value;
	}
	return value;
}

double BrownianBridge::getLogProb() {

	double** prec = unitVariance.GetInvMatrix();
	double logDet = -unitVariance.GetLogDeterminant();

	double logprob = 0;

	for(int i=0; i<nSegments; i++) {
		logprob += logProbNormal(bridge[i+1], bridge[i], segmentLength[i], prec, logDet);
	} 

	logprob -= logProbNormal(bridge[nSegments], bridge[0], length, prec, logDet);
	return logprob;
 
}

double BrownianBridge::getLogProb(CovMatrix& covar) {

	double** prec = covar.GetInvMatrix();
	double logDet = -covar.GetLogDeterminant();

	double logprob = 0;

	for(int i=0; i<nSegments; i++) {
		logprob += logProbNormal(bridge[i+1], bridge[i], segmentLength[i], prec, logDet);
	} 

	logprob -= logProbNormal(bridge[nSegments], bridge[0], length, prec, logDet);
	return logprob;
 
}

void BrownianBridge::Shift(BrownianBridge* other, double tuning)	{

	for(int i=0; i<nSegments+1; i++) {
		for(int j=0; j<dim; j++)
			bridge[i][j] += tuning * other->bridge[i][j];
	}
}

double BrownianBridge::getLogProbSubBridge(int beg, int end) {
	double** prec = unitVariance.GetInvMatrix();
	double logDet = -unitVariance.GetLogDeterminant();

	double logprob = 0;
	double sublength = 0;
	for(int i=beg; i<end; i++) {
		logprob += logProbNormal(bridge[i+1], bridge[i], segmentLength[i], prec, logDet);
		sublength+=segmentLength[i];
	} 


	logprob -= logProbNormal(bridge[end], bridge[beg], sublength, prec, logDet);
	return logprob;
}


double BrownianBridge::logProbNormal(double* x, double* mean, double slength, double** prec, double logPrecDet) {

	if(slength==0) {
		cerr << "Error : null segment or bridge length" << endl;
		throw 0;
	}

	double value = -0.5*dim*log(2*Pi);
	value+= 0.5*(logPrecDet-dim*log(slength));
	for(int i=0; i<dim; i++)
	for(int j=0; j<dim; j++) {
		value-= 0.5*prec[i][j]*(x[i]-mean[i])*(x[j]-mean[j])/slength;
	}

	return value;
}



void BrownianBridge::generateBridge() {


	removeBridge(bridge);

	bridge = createBridge();

	for(int j=0; j<dim; j++) {
		bridge[0][j]= 0;
		bridge[nSegments][j]= 0;
	}
	double lengthRight = length;
	if(length <1e-10) {
		lengthRight = 1e-8;
		/*
		cerr << "Error in brownian bridge::generate bridge : null length" << endl;
		exit(0);
		*/
	}
	for(int i=0; i<nSegments-1; i++) {

		RealVector vectorMean(bridge[i], dim);
		vectorMean.ScalarMultiplication((lengthRight-segmentLength[i])/lengthRight);

		double pointVariance = segmentLength[i]*(lengthRight-segmentLength[i])/lengthRight;

		unitVariance.drawVal(bridge[i+1]);

		for(int j=0; j<dim; j++)
			bridge[i+1][j] = bridge[i+1][j]*sqrt(pointVariance)+vectorMean[j];

		lengthRight-=segmentLength[i];

	}
  
 
}

void BrownianBridge::generateSubBridge(int beg, int end) {

	if(beg <0 || beg > end || end>nSegments) {
		cerr << "Error in generateSubBridge : " << endl;
		cerr << "beg = " << beg << " ; end = " << end << " ; nSegments = " << nSegments << endl;
		exit(1);
	}

	double lengthRight = 0;
	for(int i=beg; i<end; i++)
		lengthRight += segmentLength[i];

	for(int i=beg; i<end-1; i++) {

		double* vectorMean = new double[dim];
		for(int j=0; j<dim; j++) {
			vectorMean[j] = ((lengthRight-segmentLength[i])*bridge[i][j] + segmentLength[i]*bridge[end][j]) / lengthRight;
		}

		double pointVariance = segmentLength[i]*(lengthRight-segmentLength[i])/lengthRight;
  
		unitVariance.drawVal(bridge[i+1]);

		for(int j=0; j<dim; j++)
			bridge[i+1][j] = bridge[i+1][j]*sqrt(pointVariance)+vectorMean[j];

		lengthRight-=segmentLength[i];
		delete[] vectorMean;
	}
  
		  
}

void BrownianBridge::showBridge() {
	for(int i=0; i<nSegments+1; i++) {
		for(int j=0; j<dim; j++) {
			cout << bridge[i][j] << ";";
		}
		cout << endl;
	}
}


double BrownianBridge::ProposeMove(double tuning) {

	CovMatrix otherVar = unitVariance;
	// otherVar.SetIdentity();
	otherVar*=tuning;
	BrownianBridge *other = new BrownianBridge(initAge, finalAge, otherVar, nSegments, segmentLength);
	other->generateBridge();
	for(int i=0; i<nSegments+1; i++) {
		for(int j=0; j<dim; j++)
			bridge[i][j]+=other->bridge[i][j];
	}
	delete other;
	return 0;
}

double BrownianBridge::ProposeMove(double tuning, CovMatrix& otherVar) {

	otherVar*=tuning;
	BrownianBridge *other = new BrownianBridge(initAge, finalAge, otherVar, nSegments, segmentLength);
	other->generateBridge();
	for(int i=0; i<nSegments+1; i++) {
		for(int j=0; j<dim; j++)
			bridge[i][j]+=other->bridge[i][j];
	}
	delete other;
	otherVar*= 1.0 / tuning;
	return 0;
}

double BrownianBridge::SimpleProposeMove(double tuning) {

	CovMatrix otherVar = unitVariance;
	otherVar.SetIdentity();
	otherVar*=tuning;
	BrownianBridge *other = new BrownianBridge(initAge, finalAge, otherVar, nSegments, segmentLength);
	other->generateBridge();
	for(int i=0; i<nSegments+1; i++) {
		for(int j=0; j<dim; j++)
			bridge[i][j]+=other->bridge[i][j];
	}
	delete other;
	return 0;
}

double BrownianBridge::getIntegral(int param, double initValue, double finalValue) {

	double deltaValue = finalValue-initValue;
	double integral = 0.5*(segmentLength[0] + segmentLength[nSegments-1]*exp(deltaValue));
	double cumulateLength = segmentLength[0];
	for(int i=1; i<nSegments; i++) {
		integral += 0.5*(segmentLength[i-1]+segmentLength[i]) * exp(bridge[i][param]) * exp(deltaValue*cumulateLength/length);
		cumulateLength+=segmentLength[i];
	}
	integral *= exp(initValue);

	return integral;
}

void BrownianBridge::LeftMultiply(double** P)	{

	double* aux = new double[dim];
	for (int t=0; t<=nSegments; t++)	{
		for (int i=0; i<dim; i++)	{
			double tmp = 0;
			for (int j=0; j<dim; j++)	{
				tmp += P[i][j] * bridge[t][j];
			}
			aux[i] = tmp;
		}
		for (int i=0; i<dim; i++)	{
			bridge[t][i] = aux[i];
		}
	}
	delete[] aux;
}

double** BrownianBridge::createBridge(int N) {
	N = (N==-1 ? nSegments : N);
 
	double** b = new double*[N+1];
	for(int i=0; i<N+1; i++)
		b[i] = new double[dim];
	return b;
}

void BrownianBridge::removeBridge(double** b, int N) {
	N = (N==-1 ? nSegments : N);
	if(b) {
		for(int i=0; i<N+1; i++)
			delete[] b[i];
		delete[] b;
		b = 0;
	}
}


void BrownianBridge::findPointBefore(double time, double &timeBefore, double* &valueBefore) {
	if(time > initAge) {
		cerr << "Error : trying to find a time before the beginning of a branch" << endl;
		exit(0);
	}

	int n=0;
	double t = initAge;

	while(n<nSegments && t-segmentLength[n]-10e-8 > time) {
		t-=segmentLength[n];
		n++;
	}


	timeBefore = t;
	valueBefore = bridge[n];



}


void BrownianBridge::findPointAfter(double time, double &timeAfter, double* &valueAfter) {
	if(time < finalAge) {
		cerr << "Error : trying to find a time after the end of a branch" << endl;
		exit(0);
	}


	int n=nSegments;
	double t = finalAge;
	while(n>0 && t+segmentLength[n-1] < time-10e-8) {
		t+=segmentLength[n-1];
		n--;
	}

	timeAfter = t;
	valueAfter = bridge[n];
}


int BrownianBridge::findPointBefore(double time) {
	if(time > initAge) {
		cerr << "Error : trying to find a time before the beginning of a branch" << endl;
		exit(0);
	}

	int n=0;
	double t = initAge;

	while(n<nSegments && t-segmentLength[n]-10e-8 > time) {
		t-=segmentLength[n];
		n++;
	}



	return n;
}


int BrownianBridge::findPointAfter(double time) {
	if(time < finalAge) {
		cerr << "Error : trying to find a time after the end of a branch" << endl;
		exit(0);
	}


	int n=nSegments;
	double t = finalAge;
	while(n>0 && t+segmentLength[n-1] < time-10e-8) {
		t+=segmentLength[n-1];
		n--;
	}


	return n;
}


void BrownianBridge::addSlope(double* initValue, double* finalValue) {

	if(abs(bridge[0][0])>10e-10 || abs(bridge[nSegments][0])>10e-10) {
		cerr << "Bridge slope error (adding)" << endl;
		cerr << bridge[0][0] << " " << bridge[nSegments][0] << endl;
		exit(0);
	}

	double lengthLeft = 0;
	for(int i=0; i<nSegments; i++) {
		for(int j=0; j<dim; j++)
			bridge[i][j] += initValue[j] + lengthLeft*(finalValue[j]-initValue[j])/length;
		lengthLeft+=segmentLength[i];
	}
	for(int j=0; j<dim; j++)
		bridge[nSegments][j] = finalValue[j];
}


void BrownianBridge::removeSlope() {

	double lengthLeft = 0;
	double *initValue = new double[dim];
	for(int j=0; j<dim; j++)
		initValue[j] = bridge[0][j];

	for(int i=0; i<nSegments; i++) {
		for(int j=0; j<dim; j++)
			bridge[i][j] = bridge[i][j] - initValue[j] - lengthLeft*(bridge[nSegments][j]-initValue[j])/length;
		lengthLeft+=segmentLength[i];
	}
	for(int j=0; j<dim; j++)
		bridge[nSegments][j] = 0;

	if(abs(bridge[0][0])>10e-10 || abs(bridge[nSegments][0])>10e-10) {
		cerr << "Bridge slope error (removing)" << endl;
		cerr << bridge[0][0] << " " << bridge[nSegments][0] << " " << nSegments << endl;
		exit(0);
	}

	delete[] initValue;
}




double BrownianBridge::initAgeMove(double newAge, double *initVal) {
	if(segmentation!=SEGM_ABSOLUTE) {
		cerr << "Error : horizontal move with non absolute segmentation" << endl;
		exit(0);
	}
	double hastings = 0;

	try {
		hastings+=getLogProbSubBridge(0, findPointAfter(initAge<newAge ? initAge : newAge));
	}
	catch(int i) {
		cerr << "init age move - first" << endl;
		cerr << 0 << " " << findPointAfter(initAge<newAge ? initAge : newAge) << endl;
		exit(0);
	}
	int oldNSegments = nSegments;
	double oldAge = initAge;
	initAge = newAge;
	buildSegments();

	double** oldBridge = bridge;
	bridge = createBridge();

	for(int j=0; j<dim; j++)
		bridge[0][j] = initVal[j];			   //New extremities

	if(newAge > oldAge) {	  //Adding points at the beginning of the bridge
		for(int i=1; i<oldNSegments+1; i++) {
			for(int j=0; j<dim; j++)
				bridge[i+nSegments-oldNSegments][j] = oldBridge[i][j];
		}
	  
		generateSubBridge(0, nSegments-oldNSegments+1);
		try {
			hastings-=getLogProbSubBridge(0, nSegments-oldNSegments+1);
		}
		catch(int i) {
			cerr << "init age move - second" << endl;
			cerr << 0 << " " << nSegments-oldNSegments+1 << endl;
			exit(0);
		}

	}

	else {	   //Removing points at the beginning of the bridge

		for(int i=1; i<nSegments+1; i++) {
			for(int j=0; j<dim; j++)
				bridge[i][j] = oldBridge[i+oldNSegments-nSegments][j];
		}

	}
	removeBridge(oldBridge, oldNSegments); 

	return hastings;

}
double BrownianBridge::finalAgeMove(double newAge, double *finalVal) {
	if(segmentation!=SEGM_ABSOLUTE) {
		cerr << "Error : horizontal move with non absolute segmentation" << endl;
		exit(0);
	}
	double hastings = 0;
	try {
		hastings+=getLogProbSubBridge(findPointBefore(finalAge>newAge ? finalAge : newAge), nSegments);
	}
	catch(int i) {
		cerr << "final age move - first" << endl;
		cerr << findPointBefore(finalAge>newAge ? finalAge : newAge) << " " << nSegments << endl;
		exit(0);
	}

	int oldNSegments = nSegments;
	double oldAge = finalAge;
	finalAge = newAge;
	buildSegments();

	double** oldBridge = bridge;
	bridge = createBridge();
	  
	for(int j=0; j<dim; j++)
		bridge[nSegments][j] = finalVal[j];			   //New extremities
	if(newAge < oldAge) {	  //Adding points at the end of the bridge

		//The first oldNSegments points are copied from the old bridge
		for(int i=0; i<oldNSegments; i++) {
			for(int j=0; j<dim; j++)
				 bridge[i][j] = oldBridge[i][j];
		}
	  
		//The next NSegments - oldNSegments are generated iteratively
		generateSubBridge(oldNSegments-1, nSegments);
		try {
			hastings-=getLogProbSubBridge(oldNSegments-1, nSegments);
		}
		catch(int i) {
			cerr << "final age move - second" << endl;
			cerr << oldNSegments-1 << " " << nSegments << endl;
			exit(0);
		}

	}
	else {	   //Removing points at the end of the bridge

		for(int i=0; i<nSegments; i++) {
			for(int j=0; j<dim; j++)
				bridge[i][j] = oldBridge[i][j];
		}

	}
	removeBridge(oldBridge, oldNSegments);  

	return hastings;
}


