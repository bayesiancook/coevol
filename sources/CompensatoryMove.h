

class MultiVariateRootCompensatoryMove : public MCUpdate, public Mnode {

	MultiVariateTreeProcess* process;
	Chronogram* chronogram;
	PhyloProcess* phyloprocess;
	double tuning;
	int index;

	public:

	MultiVariateRootCompensatoryMove(Chronogram* inchronogram, MultiVariateTreeProcess* inprocess, PhyloProcess* inphyloprocess, double intuning, int inindex){
		process = inprocess;
		chronogram = inchronogram;
		phyloprocess = inphyloprocess;
		tuning = intuning;
		index = inindex;
		process->RecursiveRegister(this,process->GetRoot());
		chronogram->RecursiveRegister(this,chronogram->GetRoot());
		phyloprocess->RegisterRootChildren(this);
	}

	double Move(double tuning_modulator){


		Corrupt(true);

		// make all the complicated calculus

		Link* root = chronogram->GetRoot();
		if (root->Next()->Next()->Next() != root)	{
			cerr << "error in Compensatory Move: not a binary tree\n";
			exit(1);
		}

		Link* rootleft = root->Next()->Out();
		Link* rootright = root->Next()->Next()->Out();

		double tl = chronogram->GetBranchTimeLength(rootleft);
		double tr = chronogram->GetBranchTimeLength(rootright);

		if (tl > 1)	{
			cerr << "time not correct : " << tl << '\n';
			exit(1);
		}
		if (tr > 1)	{
			cerr << "time not correct : " << tr << '\n';
			exit(1);
		}

		double r0 = process->GetExpVal(root,index);
		double rl = process->GetExpVal(rootleft,index);
		double rr = process->GetExpVal(rootright,index);

		double ll = tl * 0.5 * (r0 + rl);
		double lr = tr * 0.5 * (r0 + rr);

		double tot = ll + lr;
		double u = tuning * tuning_modulator * (Random::Uniform() - 0.5);

		double ml = ll + u;
		while ((ml < 0) || (ml > tot))	{
			if (ml < 0)	{
				ml = -ml;
			}
			if (ml > tot)	{
				ml = 2*tot - ml;
			}
		}
		double mr = tot - ml;

		double deltalength = ml - ll;
		if (deltalength > 0)	{
			phyloprocess->ShiftRootMappings(rootleft,rootright,root,deltalength);
		}
		else	{
			phyloprocess->ShiftRootMappings(rootright,rootleft,root,-deltalength);
		}

		double al = rl * (1 - tl);
		double ar = rr * (1 - tr);

		double bl = 2*ml + r0 + al;
		double br = 2*mr + r0 + ar;

		double dl = bl*bl - 8*r0*ml;
		double dr = br*br - 8*r0*mr;

		if (dl < 0)	{
			cerr << "negative discriminant: " << dl << '\n';
			exit(1);
		}
		if (dr < 0)	{
			cerr << "negative discriminant: " << dr << '\n';
			exit(1);
		}

		double ul = (bl - sqrt(dl))/2/r0;
		double ur = (br - sqrt(dr))/2/r0;

		if ((ul < 0) || (ul > 1+ 1e-6))	{
			cerr << "incorrect time : " << ul << '\n';
			exit(1);
		}
		if ((ur < 0) || (ur > 1 + 1e-6))	{
			cerr << "incorrect time : " << ur << '\n';
			exit(1);
		}

		double sl = (2*ml/ul - r0);
		double sr = (2*mr/ur - r0);

		double logHastings = log(fabs((1 - bl/sqrt(dl)) * (1 - br/sqrt(dr)) * (2*ll/tl/tl - r0) * (2*lr/tr/tr - r0)) / 4 / r0 / r0);

		logHastings += log(tl/ul) + log(tr/ur);

		int nl = chronogram->GetTree()->GetSize(rootleft);
		int nr = chronogram->GetTree()->GetSize(rootright);

		if (nl > 1)	{
			logHastings += (nl-1) * log(rl/sl);
		}
		else	{
			cerr << "in root compensation move: should deal with the case where root's daughter is leaf node\n";
			exit(1);
			// logHastings += log(sl/rl);
		}
		if (nr > 1)	{
			logHastings += (nr-1) * log(rr/sr);
		}
		else	{
			cerr << "in root compensation move: should deal with the case where root's daughter is leaf node\n";
			exit(1);
			// logHastings += log(sr/rr);
		}


		chronogram->MultiplyTimes(rootleft,rl/sl);
		process->RecursivePiecewiseTranslation(rootleft,log(sl/rl),index,1);

		chronogram->MultiplyTimes(rootright,rr/sr);
		process->RecursivePiecewiseTranslation(rootright,log(sr/rr),index,1);

		if (std::isinf(logHastings))	{
			cerr << "in root comp move: inf\n";
			exit(1);
		}

		if (std::isnan(logHastings))	{
			cerr << "in root comp move: nan\n";
			exit(1);
		}

		double logratio = Update() + logHastings;

		bool accepted = (log(Random::Uniform()) < logratio);
		if (! accepted)	{
			Corrupt(false);
			if (deltalength > 0)	{
				phyloprocess->ShiftRootMappings(rootright,rootleft,root,deltalength);
			}
			else	{
				phyloprocess->ShiftRootMappings(rootleft,rootright,root,-deltalength);
			}
			Restore();
		}
		return (double) accepted;
	}
};
