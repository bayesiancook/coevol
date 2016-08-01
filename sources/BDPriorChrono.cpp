
#include "BaseType.h"
#include "RandomTypes.h"
#include "BDCalibratedChronogram.h"
#include "ProbModel.h"

class PriorChrono : public ProbModel	{

	Tree* tree;
	Const<PosReal>* One;
	Const<PosReal>* RootAlpha;
	Const<PosReal>* RootBeta;
	Const<PosReal>* Chi;
	Const<PosReal>* Chi2;
	CalibrationSet* calibset;
	CalibratedChronogram* chronogram;

	public:

	PriorChrono(string treefile,string calibfile, double rootage, double rootstdev, double chi, double chi2)	{

		tree = new Tree(treefile);
		One = new Const<PosReal>(1);

		double a = rootage * rootage / rootstdev / rootstdev;
		double b = rootage / rootstdev / rootstdev;

		RootAlpha = new Const<PosReal>(a);
		RootBeta = new Const<PosReal>(b);
		Chi = new Const<PosReal>(chi);
		Chi2 = new Const<PosReal>(chi2);

		calibset = new FileCalibrationSet(calibfile, tree);
		chronogram = new BDCalibratedChronogram(tree,One,Chi,Chi2,RootAlpha,RootBeta,calibset);

		RootRegister(One);
		RootRegister(Chi);
		RootRegister(Chi2);
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
	double chi = atof(argv[5]);
	double chi2 = atof(argv[6]);
	int every = atoi(argv[7]);
	int nrep = atoi(argv[8]);
	string name = argv[9];

	PriorChrono* chrono = new PriorChrono(treefile,calibfile,rootage,rootstdev,chi,chi2);
	ofstream os((name + ".trace").c_str());

	for (int rep=0; rep<nrep; rep++)	{
		for (int i=0; i<every; i++)	{
			chrono->Move(1);
		}
		chrono->Trace(os);
		os.flush();
	}
}
