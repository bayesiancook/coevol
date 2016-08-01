
#include "BaseType.h"
#include "RandomTypes.h"
#include "CalibratedChronogram.h"
#include "ProbModel.h"

class PriorChrono : public ProbModel	{

	Tree* tree;
	Const<PosReal>* One;
	Const<PosReal>* RootAlpha;
	Const<PosReal>* RootBeta;
	CalibrationSet* calibset;
	CalibratedChronogram* chronogram;

	public:

	PriorChrono(string treefile,string calibfile, double rootage, double rootstdev)	{

		tree = new Tree(treefile);
		One = new Const<PosReal>(1);

		double a = rootage * rootage / rootstdev / rootstdev;
		double b = rootage / rootstdev / rootstdev;

		RootAlpha = new Const<PosReal>(a);
		RootBeta = new Const<PosReal>(b);
		calibset = new FileCalibrationSet(calibfile, tree);
		chronogram = new CalibratedChronogram(tree,One,RootAlpha,RootBeta,calibset);
		/*
		cerr << "BD\n";
		MeanChi = new Const<PosReal>(meanchi);
		MeanChi2 = new Const<PosReal>(meanchi2);
		Chi = new Exponential(MeanChi,Exponential::MEAN);
		Chi2 = new Exponential(MeanChi2,Exponential::MEAN);
		chronogram = new BDCalibratedChronogram(tree,mu,Chi,Chi2,RootAlpha,RootBeta,calibset);
		*/
		RootRegister(One);
		RootRegister(RootAlpha);
		RootRegister(RootBeta);

		Register();
		MakeScheduler();
		Update();
		TraceHeader(cerr);
		Trace(cerr);

		cerr << "model created\n";
	}

	void MakeScheduler()	{
		scheduler.Register(new SimpleMove(chronogram,10),10,"chrono");
		scheduler.Register(new SimpleMove(chronogram,1),10,"chrono");
		scheduler.Register(new SimpleMove(chronogram,0.1),10,"chrono");
		scheduler.Register(new SimpleMove(chronogram->GetScale(),10),10,"scale");
		scheduler.Register(new SimpleMove(chronogram->GetScale(),1),10,"scale");
		scheduler.Register(new SimpleMove(chronogram->GetScale(),0.1),10,"scale");
	}

	double GetLogProb()	{
		return chronogram->GetLogProb();
	}

	void drawSample()	{
	}

	void TraceHeader(ostream& os)	{
		os << "logprob\trootage\n";
	}

	void Trace(ostream& os)	{
		os << GetLogProb() << '\t' << chronogram->GetScale()->val() << '\n';
	}

	void ToStream(ostream& os)	{
	}

	void FromStream(istream& is)	{
	}

};


int main(int argc, char* argv[])	{

	string treefile = argv[1];
	string calibfile = argv[2];
	double rootage = atof(argv[3]);
	double rootstdev = atof(argv[4]);
	int every = atoi(argv[5]);
	int nrep = atoi(argv[6]);
	string name = argv[7];

	PriorChrono* chrono = new PriorChrono(treefile,calibfile,rootage,rootstdev);
	ofstream os((name + ".trace").c_str());

	for (int rep=0; rep<nrep; rep++)	{
		for (int i=0; i<every; i++)	{
			chrono->Move(1);
		}
		chrono->Trace(os);
		os.flush();
	}
}
