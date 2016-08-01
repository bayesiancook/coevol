#ifndef MEANTIMELINE
#define MEANTIMELINE

class MeanTimeLine {

	public:

	MeanTimeLine(int insize)	{
		size = insize;
		for (int i=0; i<size; i++)	{
			mean.push_back(0);
			var.push_back(0);
			count.push_back(0);
		}
		scale = 0;
		scalecount = 0;
	}

	void Normalize()	{
		for (int i=0; i<size; i++)	{
			mean[i] /= count[i];
			var[i] /= count[i];
			var[i] -= mean[i] * mean[i];
		}
		scale /= scalecount;
	}

	void ToStream(ostream& os)	{
		for (int i=0; i<size; i++)	{
			os <<  (1 - ((double) i) / size) * scale << '\t' << mean[i]/log(10.0) << '\t' << sqrt(var[i])/log(10.0) << '\n';
		}
	}

	void Add(double val, double mintime, double maxtime)	{
		int minindex = (int) (size * (1 - maxtime));
		int maxindex = (int) (size * (1 - mintime));
		for (int i=minindex; i<maxindex; i++)	{
			mean[i] += val;
			var[i] += val*val;
			count[i] ++;
		}
	}

	void Add(double up, double down, double mintime, double maxtime)	{
		int minindex = (int) (size * (1 - maxtime));
		int maxindex = (int) (size * (1 - mintime));
		for (int i=minindex; i<maxindex; i++)	{
			double val = up + (down-up) * ((double) (i - minindex)) / (maxindex - minindex);
			mean[i] += val;
			var[i] += val*val;
			count[i] ++;
		}
	}

	void Add(BranchVarTree<PosReal>* intree, Chronogram* inchrono, bool withmean)	{
		// RecursiveAdd(intree, inchrono, intree->GetTree()->GetLCA("Procavia","Pan"), withmean);
		RecursiveAdd(intree, inchrono, intree->GetRoot(), withmean);
		CalibratedChronogram* tmp = dynamic_cast<CalibratedChronogram*>(inchrono);
		if (tmp)	{
			scale += tmp->GetScale()->val();
		}
		else	{
			scale += 1;
		}
		scalecount ++;
	}

	void Add(BranchVarTree<UnitReal>* intree, Chronogram* inchrono, bool withmean)	{
		// RecursiveAdd(intree, inchrono, intree->GetTree()->GetLCA("Procavia","Pan"), withmean);
		RecursiveAdd(intree, inchrono, intree->GetRoot(), withmean);
		CalibratedChronogram* tmp = dynamic_cast<CalibratedChronogram*>(inchrono);
		if (tmp)	{
			scale += tmp->GetScale()->val();
		}
		else	{
			scale += 1;
		}
		scalecount ++;
	}

	void Add(MultiVariateTreeProcess* intree, Chronogram* inchrono, int index)	{
		// RecursiveAdd(intree, inchrono, intree->GetTree()->GetLCA("Procavia","Pan"), withmean, index);
		RecursiveAdd(intree, inchrono, intree->GetRoot(), index);
		CalibratedChronogram* tmp = dynamic_cast<CalibratedChronogram*>(inchrono);
		if (tmp)	{
			scale += tmp->GetScale()->val();
		}
		else	{
			scale += 1;
		}
		scalecount ++;
	}

	private:

	void RecursiveAdd(MultiVariateTreeProcess* intree, Chronogram* inchrono, const Link* from, int index)	{
		// if (from != intree->GetTree()->GetLCA("Procavia","Pan"))	{
		if (! from->isRoot())	{
			double tmp1 = (*intree->GetMultiNormal(from))[index];
			double tmp2 = (*intree->GetMultiNormal(from->Out()))[index];
			// double val = (tmp1 + tmp2) / 2;
			double mintime = inchrono->GetNodeVal(from->GetNode())->val();
			double maxtime = inchrono->GetNodeVal(from->Out()->GetNode())->val();
			// Add(val,mintime,maxtime);
			Add(tmp2,tmp1,mintime,maxtime);
		}
		for (const Link* link=from->Next(); link!=from; link=link->Next())	{
			RecursiveAdd(intree,inchrono,link->Out(), index);
		}
	}

	void RecursiveAdd(BranchVarTree<PosReal>* intree, Chronogram* inchrono, const Link* from, bool withmean)	{
		// if (from != intree->GetTree()->GetLCA("Procavia","Pan"))	{
		if (! from->isRoot())	{
			double tmp = intree->GetBranchVal(from->GetBranch())->val();
			double time = inchrono->GetBranchVal(from->GetBranch())->val();
			if (withmean)	{
				tmp /= time;
			}
			// double val = tmp;
			double val = log(tmp);
			double mintime = inchrono->GetNodeVal(from->GetNode())->val();
			double maxtime = inchrono->GetNodeVal(from->Out()->GetNode())->val();
			if (mintime > maxtime)	{
				cerr << mintime << '\t' << maxtime << '\n';
				exit(1);
			}

			Add(val,mintime,maxtime);
		}
		for (const Link* link=from->Next(); link!=from; link=link->Next())	{
			RecursiveAdd(intree,inchrono,link->Out(),withmean);
		}
	}

	void RecursiveAdd(BranchVarTree<UnitReal>* intree, Chronogram* inchrono, const Link* from, bool withmean)	{
		// if (from != intree->GetTree()->GetLCA("Procavia","Pan"))	{
		if (! from->isRoot())	{
			double tmp = intree->GetBranchVal(from->GetBranch())->val();
			double time = inchrono->GetBranchVal(from->GetBranch())->val();
			if (withmean)	{
				tmp /= time;
			}
			// double val = tmp;
			double val = log(tmp / (1 - tmp));
			double mintime = inchrono->GetNodeVal(from->GetNode())->val();
			double maxtime = inchrono->GetNodeVal(from->Out()->GetNode())->val();
			if (mintime > maxtime)	{
				cerr << mintime << '\t' << maxtime << '\n';
				exit(1);
			}

			Add(val,mintime,maxtime);
		}
		for (const Link* link=from->Next(); link!=from; link=link->Next())	{
			RecursiveAdd(intree,inchrono,link->Out(),withmean);
		}
	}

	vector<double> mean;
	vector<double> var;
	vector<int> count;
	int size;
	double scale;
	int scalecount;

};

#endif

