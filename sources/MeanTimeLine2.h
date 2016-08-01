

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
			os << - (1 - ((double) i) / size) * scale << '\t' << mean[i] << '\t' << sqrt(var[i]) << '\n';
		}
	}

	void Add(double val, double mintime, double maxtime)	{
		int minindex = (int) (size * (1 - maxtime));
		int maxindex = (int) (size * (1 - mintime));
		if ((minindex <0) || (minindex >=size))	{
			cerr << "overflow\n";
			cerr << minindex << '\t' << size << '\n';
			exit(1);
		}
		if (maxindex == size)	{
			maxindex --;
		}
		if ((maxindex <0) || (maxindex >=size))	{
			cerr << "overflow\n";
			cerr << maxindex << '\t' << size << '\n';
			exit(1);
		}
		for (int i=minindex; i<maxindex; i++)	{
			mean[i] += val;
			var[i] += val*val;
			count[i] ++;
		}
	}

	void Add(double up, double down, double mintime, double maxtime)	{
		int minindex = (int) (size * (1 - maxtime));
		int maxindex = (int) (size * (1 - mintime));
		if ((minindex <0) || (minindex >=size))	{
			cerr << "overflow\n";
			cerr << minindex << '\t' << size << '\t' << mintime << '\n';
			exit(1);
		}
		if (maxindex == size)	{
			maxindex --;
		}
		if ((maxindex <0) || (maxindex >=size))	{
			cerr << "overflow\n";
			cerr << maxindex << '\t' << size << '\t' << maxtime << '\n';
			exit(1);
		}
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

	/*
	void Add(BranchVarTree<PosReal>* intree, Chronogram* inchrono, bool withmean)	{
		// RecursiveAdd(intree, inchrono, intree->GetTree()->GetLCA("Procavia","Pan"), withmean);
		RecursiveAdd(intree, inchrono, intree->GetRoot(), withmean);
	}

	void Add(MultiVariateTreeProcess* intree, Chronogram* inchrono, int index)	{
		// RecursiveAdd(intree, inchrono, intree->GetTree()->GetLCA("Procavia","Pan"), withmean, index);
		RecursiveAdd(intree, inchrono, intree->GetRoot(), index);
	}
	*/

	private:

	double GetRelativeTime(Chronogram* inchrono, const Node* node)	{
		double time = inchrono->GetNodeVal(node)->val();
		CoalCalibratedChronogram* tmp = dynamic_cast<CoalCalibratedChronogram*>(inchrono);
		if (tmp)	{
			// time /= tmp->GetScale()->val();
			// cerr << time << '\t' << tmp->GetScale()->val() << '\t' << inchrono->GetNodeVal(inchrono->GetRoot()->GetNode())->val() << '\n';
			time /= inchrono->GetNodeVal(inchrono->GetRoot()->GetNode())->val();
		}
		return time;
	}

	void RecursiveAdd(MultiVariateTreeProcess* intree, Chronogram* inchrono, const Link* from, int index)	{
		// if (from != intree->GetTree()->GetLCA("Procavia","Pan"))	{
		if (! from->isRoot())	{
			double tmp1 = (*intree->GetMultiNormal(from))[index];
			double tmp2 = (*intree->GetMultiNormal(from->Out()))[index];
			// double val = (tmp1 + tmp2) / 2;
			double mintime = GetRelativeTime(inchrono,from->GetNode());
			double maxtime = GetRelativeTime(inchrono,from->Out()->GetNode());
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
			double mintime = GetRelativeTime(inchrono,from->GetNode());
			double maxtime = GetRelativeTime(inchrono,from->Out()->GetNode());
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
