
#include "BaseType.h"
#include "RandomTypes.h"
#include "ProbModel.h"
#include "ConjugateOneMatrixPhyloProcess.h"
#include "ConjugateGammaTree.h"

#include "Move.h"

#include "CodonSubMatrix.h"
#include "CodonSequenceAlignment.h"


#include "SiteMG.h"

// a small program,
// bypassing the use of the Chain class
//

int main(int argc, char* argv[])	{

	string datafile = argv[1];
	string treefile = argv[2];
	int burnin = atoi(argv[3]);
	string name = argv[4];

	SiteMGModel* model = new SiteMGModel(datafile,treefile,burnin);

	cerr << "start\n";

	ofstream os((name + ".trace").c_str());
	model->TraceHeader(os);

	ofstream ros((name + ".omega").c_str());

	while (1)	{
		model->Move(1);
		model->Move(0.1);
		model->Move(0.01);
		model->Trace(os);
		model->omega->ToStream(ros);
		ros << '\n';
		ros.flush();
		ofstream mon_os((name + ".monitor").c_str());
		ofstream dmon_os((name + ".detailedmonitor").c_str());
		model->Monitor(mon_os, dmon_os);
		model->PushOmega();
		model->WriteOmega(ros);
		ros.flush();
	}
}
