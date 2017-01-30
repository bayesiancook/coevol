#include "Chain.h"
#include "Chrono.h"
#include <cmath>

Chain::Chain()	{
	every = 1;
	until = -1;
	size = 0;
	model = 0;
	name = "";
}

void Chain::MakeFiles(int force)	{
	if (ifstream((name + ".param").c_str()) && (! force))	{
		cerr << "already existing chain, cannot override (unless in forcing mode)\n";
		exit(1);
	}
	ofstream param_os((name + ".param").c_str());
	ofstream chain_os((name + ".chain").c_str());
	ofstream mon_os((name + ".monitor").c_str());
	ofstream trace_os((name + ".trace").c_str());
	model->TraceHeader(trace_os);
}

void Chain::Monitor()	{
	ofstream trace_os((name + ".trace").c_str(), ios_base::app);
	model->Trace(trace_os);
	ofstream mon_os((name + ".monitor").c_str());
	ofstream mon_det_os((name + ".details").c_str());
	model->Monitor(mon_os, mon_det_os);
	// model->Monitor(mon_os);
}

void Chain::SavePoint()	{
	ofstream chain_os((name + ".chain").c_str(), ios_base::app);
	model->ToStream(chain_os);
	size++;
}

void Chain::Reset(int force)	{
	size = 0;
	MakeFiles(force);
	Save();
}

void Chain::Move()	{
	for (int i=0; i<every; i++)	{
		model->Move();
	}
	/*
	double delta = model->Update();
	if (fabs(delta) > 1e-4)	{
		cerr << "error : model corrupted: " << delta << ". Are you sure about the updates during Monte Carlo?\n";
		exit(1);
	}
	*/
	SavePoint();
	Save();
	Monitor();
}

void Chain::Start()	{
	ofstream run_os((name + ".run").c_str());
	run_os << 1 << '\n';
	run_os.close();
	Run();
}

int Chain::GetRunningStatus()	{
	ifstream run_is((name + ".run").c_str());
	int run;
	run_is >> run;
	return run;
}

void Chain::Run()	{
	while (GetRunningStatus() && ((until == -1) || (size <= until)))	{
		Chrono chrono;
		chrono.Reset();
		chrono.Start();
		Move();
		chrono.Stop();
		/*
		ofstream check_os((name + ".time").c_str());
		check_os << chrono.GetTime() / 1000 << '\n';
		*/
	}
	ofstream run_os((name + ".run").c_str());
	run_os << 0 << '\n';
}

