
#include "PhyloProcess.h"

int PhyloProcess::GetMaxNstate()	{
	if (!MaxNstate)	{
		MaxNstate = 0;
		for (int i=0; i<GetNsite(); i++)	{
			if (MaxNstate < GetNstate(i))	{
				MaxNstate = GetNstate(i);
			}
		}
	}
	return MaxNstate;
}

void PhyloProcess::SetName(string inname)	{
	RecursiveSetName(GetRoot(), inname);
}

void PhyloProcess::RecursiveSetName(const Link* from, string inname)	{
	string name = (from->isRoot()) ? " root " : " non root ";
	for (int i=0; i<GetNsite(); i++)	{
		if (! isMissing(from,i))	{
			GetPath(0,i)->SetName(inname + name);
		}
	}
	for (const Link* link=from->Next(); link!=from; link=link->Next())	{
		RecursiveSetName(link->Out(),inname);
	}
}

void PhyloProcess::CreateMissingMap()	{
	RecursiveCreateMissingMap(GetRoot());
	for (int i=0; i<GetNsite(); i++)	{
		FillMissingMap(GetRoot(), i);
	}
	ComputeTotalMissingPerSite(GetRoot());
}

void PhyloProcess::RecursiveCreateMissingMap(const Link* from)	{
	bool* array = new bool[GetNsite()];
	missingmap[from->GetNode()] = array;
	for (const Link* link=from->Next(); link!=from; link=link->Next())	{
		RecursiveCreateMissingMap(link->Out());
	}
}

void PhyloProcess::ComputeTotalMissingPerSite(const Link* from)	{
	int tot = 0;
	for (int i=0; i<GetNsite(); i++)	{
		if (isMissing(from,i))	{
			tot++;
		}
	}
	totmissingmap[from->GetNode()] = tot;
	for (const Link* link=from->Next(); link!=from; link=link->Next())	{
		ComputeTotalMissingPerSite(link->Out());
	}
}

bool PhyloProcess::FillMissingMap(const Link* from, int i)	{
	if (from->isLeaf())	{
		missingmap[from->GetNode()][i] = data->isMissing(from->GetNode()->GetIndex(),i);
		return missingmap[from->GetNode()][i];
	}
	bool allmiss = true;
	for (const Link* link=from->Next(); link!=from; link=link->Next())	{
		allmiss &= FillMissingMap(link->Out(),i);
	}
	missingmap[from->GetNode()][i] = allmiss;
	return allmiss;
}

void PhyloProcess::SetStateSpace()	{
	RecursiveSetStateSpace(GetRoot());
}

void PhyloProcess::RecursiveSetStateSpace(const Link* from)	{
	for (int i=0; i<GetNsite(); i++)	{
		if (! isMissing(from,i))	{
			GetPath(from->GetBranch(),i)->SetStateSpace(GetStateSpace());
		}
	}
	for (const Link* link=from->Next(); link!=from; link=link->Next())	{
		RecursiveSetStateSpace(link->Out());
	}
}

void PhyloProcess::RecursiveCreate(const Link* from)	{
	int* state = new int[GetNsite()];
	RandomBranchSitePath** path = CreateRandomBranchSitePath(from);
	pathmap[from->GetBranch()] = path;
	statemap[from->GetNode()] = state;
	for (const Link* link=from->Next(); link!=from; link=link->Next())	{
		RecursiveCreate(link->Out());
	}
}

RandomBranchSitePath** PhyloProcess::CreateRandomBranchSitePath(const Link* link)	{
	RandomBranchSitePath** array = new RandomBranchSitePath*[GetNsite()];
	for (int i=0; i<GetNsite(); i++)	{
		if (! isMissing(link,i))	{
			array[i] = CreateRandomBranchSitePath(link,i);
		}
		else	{
			array[i] = 0;
		}
	}
	return array;
}

void PhyloProcess::RecursiveDelete(const Link* from)	{
	for (const Link* link=from->Next(); link!=from; link=link->Next())	{
		RecursiveDelete(link->Out());
	}
	DeletePath(from);
	delete[] statemap[from->GetNode()];
}

void PhyloProcess::DeletePath(const Link* link)	{
	RandomBranchSitePath** path = pathmap[link->GetBranch()];
	for (int i=0; i<GetNsite(); i++)	{
		delete path[i];
	}
	delete[] path;
}

PhyloProcess::PhyloProcess(LengthTree* intree, SequenceAlignment* indata, bool inbranchmap)	{
	branchmap = inbranchmap;
	MaxNstate = 0;
	tree = intree;
	data = indata;
	maxtrial = DEFAULTMAXTRIAL;
}

void PhyloProcess::SetData(SequenceAlignment* indata)	{
	data = indata;
}

void PhyloProcess::Unfold()	{

	sitearray = new int[GetNsite()];
	sitelnL = new double[GetNsite()];
	for (int i=0; i<GetNsite(); i++)	{
		sitearray[i] = 1;
		//sitearray[i] = -1;
	}
	CreateMissingMap();
	ResetFlagMap(GetRoot(),true);
	RecursiveCreate(GetRoot());
	RecursiveCreateTBL(GetRoot(),GetMaxNstate());
	sitemapping = new PhyloProcessSiteMapping*[GetNsite()];
	for (int i=0; i<GetNsite(); i++)	{
		sitemapping[i] = new PhyloProcessSiteMapping(this,i);
	}
	ClampData();
	SetStateSpace();
}

void PhyloProcess::Cleanup()	{
	RecursiveDeleteTBL(GetRoot());
	RecursiveDelete(GetRoot());
	delete[] sitearray;
	delete[] sitelnL;
}

void PhyloProcess::RecursiveCreateTBL(const Link* from, int innstate)	{
	condlmap[from] = new double[innstate+1];
	for(Link* link=from->Next(); link!=from; link=link->Next())	{
		condlmap[link] = new double[innstate+1];
		RecursiveCreateTBL(link->Out(), innstate);
	}
}

void PhyloProcess::RecursiveDeleteTBL(const Link* from)	{
	for(Link* link=from->Next(); link!=from; link=link->Next())	{
		RecursiveDeleteTBL(link->Out());
		delete[] condlmap[link];
	}
	delete[] condlmap[from];
}

PhyloProcess::~PhyloProcess()	{
	Cleanup();
}

double PhyloProcess::SiteLogLikelihood(int site)	{
	if (isMissing(GetRoot(),site))	{
		return 0;
	}
	Pruning(GetRoot(),site);
	double ret = 0;
	double* t = GetCondLikelihood(GetRoot());
	const double* stat = GetPath(0,site)->GetStationary();
	// const double* stat = GetBranchSiteSubstitutionProcess(0,site)->GetStationary();
	for (int k=0; k<GetNstate(site); k++)	{
		ret += t[k] * stat[k];
	}
	if (ret == 0)	{
		cerr << "pruning : 0 \n";
		for (int k=0; k<GetNstate(site); k++)	{
			cerr << t[k] << '\t' << stat[k] << '\n';
		}
		exit(1);
	}
	return log(ret) + t[GetNstate(site)];
}

double PhyloProcess::FastSiteLogLikelihood(int site)	{
	if (isMissing(GetRoot(),site))	{
		return 0;
	}
	// Pruning(GetRoot(),site);
	double ret = 0;
	double* t = GetCondLikelihood(GetRoot());
	const double* stat = GetPath(0,site)->GetStationary();
	// const double* stat = GetBranchSiteSubstitutionProcess(0,site)->GetStationary();
	for (int k=0; k<GetNstate(site); k++)	{
		ret += t[k] * stat[k];
	}
	if (ret == 0)	{
		cerr << "pruning : 0 \n";
		for (int k=0; k<GetNstate(site); k++)	{
			cerr << t[k] << '\t' << stat[k] << '\n';
		}
		exit(1);
	}
	sitelnL[site] = log(ret) + t[GetNstate(site)];
	return sitelnL[site];
}

double PhyloProcess::GetFastLogProb()	{
	double total = 0;
	for (int i=0; i<GetNsite(); i++)	{
		/*
		if (sitelnL[i] == -1)	{
			Pruning(GetRoot(),i);
			FastSiteLogLikelihood(i);
		}
		*/
		total += sitelnL[i];
	}
	return total;
}

double PhyloProcess::GetLogProb()	{
	// return GetPathLogProb();
	double total = 0;
	for (int i=0; i<GetNsite(); i++)	{
		/*
		if (sitelnL[i] == -1)	{
			Pruning(GetRoot(),i);
			FastSiteLogLikelihood(i);
		}
		total += sitelnL[i];
		*/
		total += SiteLogLikelihood(i);
	}
	return total;
}

double PhyloProcess::GetPathLogProb()	{
	double total = 0;
	for (int i=0; i<GetNsite(); i++)	{
		if (! isMissing(GetRoot()->GetNode(),i))	{
			total += GetPathLogProb(GetRoot(),i);
		}
	}
	return total;
}

double PhyloProcess::GetPathLogProb(const Link* from, int site)	{
	double total = GetPath(from->GetBranch(),site)->GetFastLogProb();
	/*
	double total2 = GetPath(from->GetBranch(),site)->GetLogProb();
	if (fabs(total - total2) > 1e-6)	{
		cerr << "error in phyloprocess: get path log prob\n";
		cerr << total - total2 << '\n';
		exit(1);
	}
	*/
	for(const Link* link=from->Next(); link!=from; link=link->Next())	{
		if (! isMissing(link->Out()->GetNode(),site))	{
			total += GetPathLogProb(link->Out(),site);
		}
	}
	return total;
}

void PhyloProcess::Pruning(const Link* from, int site)	{

	double* t = GetCondLikelihood(from);
	if (! flagmap[from->GetNode()])	{
		for (int k=0; k<GetNstate(site); k++)	{
			t[k] = 0;
		}
		t[GetNstate(site)] = 0;
		t[GetState(from->GetNode(),site)] = 1;
	}
	else	{
	if (from->isLeaf())	{
		int totcomp = 0;
		for (int k=0; k<GetNstate(site); k++)	{
			if (isDataCompatible(from->GetNode()->GetIndex(),site,k))	{
				t[k] = 1;
				totcomp ++;
			}
			else	{
				t[k] = 0;
			}
		}
		if (! totcomp)	{
			cerr << "error : no compatibility\n";
			cerr << GetData(from->GetNode()->GetIndex(),site) << '\n';
			exit(1);
		}

		t[GetNstate(site)] = 0;
	}
	else	{
		for (int k=0; k<GetNstate(site); k++)	{
			t[k] = 1.0;
		}
		t[GetNstate(site)] = 0;
		for(const Link* link=from->Next(); link!=from; link=link->Next())	{
			if (! isMissing(link->GetBranch(),site))	{
				double* tbl = GetCondLikelihood(link);
				Pruning(link->Out(),site);
				GetPath(link->GetBranch(),site)->BackwardPropagate(GetCondLikelihood(link->Out()),GetCondLikelihood(link));
				for (int k=0; k<GetNstate(site); k++)	{
					t[k] *= tbl[k];
				}
				t[GetNstate(site)] += tbl[GetNstate(site)];
			}
		}
		double max = 0;
		for (int k=0; k<GetNstate(site); k++)	{
			if (t[k] <0)	{
				/*
				cerr << "error in pruning: negative prob : " << t[k] << "\n";
				exit(1);
				*/
				t[k] = 0;
			}
			if (max < t[k])	{
				max = t[k];
			}
		}
		if (max == 0)	{
			cerr << "max = 0\n";
			cerr << "error in pruning: null likelihood\n";
			if (from->isRoot())	{
				cerr << "is root\n";
			}
			// GetBranchSiteSubstitutionProcess(from->GetBranch(),site)->GetSubMatrix()->ToStream(cerr);
			// cerr << GetBranchSiteSubstitutionProcess(from->GetBranch(),site)->GetTime() << '\n';
			// cerr << GetBranchSiteSubstitutionProcess(from->GetBranch(),site)->GetRate() << '\n';
			cerr << '\n';
			/*
			for (int k=0; k<GetNstate(site); k++)	{
				cerr << t[k] << '\t';
			}
			cerr << '\n';
			*/
			exit(1);
			max = 1e-20;
		}
		for (int k=0; k<GetNstate(site); k++)	{
			t[k] /= max;
		}
		t[GetNstate(site)] += log(max);
	}
	}
}

void PhyloProcess::BackupNodeStates(const Link* from, int site)	{

	bkstatemap[from->GetNode()] = statemap[from->GetNode()][site];
	for (const Link* link=from->Next(); link!=from; link=link->Next())	{
		BackupNodeStates(link->Out(),site);
	}
}

void PhyloProcess::RestoreNodeStates(const Link* from, int site)	{

	statemap[from->GetNode()][site] = bkstatemap[from->GetNode()];
	for (const Link* link=from->Next(); link!=from; link=link->Next())	{
		RestoreNodeStates(link->Out(),site);
	}
}

void PhyloProcess::PruningAncestral(const Link* from, int site)	{
	if (! flagmap[from->GetNode()])	{
		cerr << "error in pruning ancestral: flag map not set up correctly\n";
		exit(1);
	}
	int& state = GetState(from->GetNode(),site);
	if (from->isRoot())	{
		double* aux = new double[GetNstate(site)];
		double* cumulaux = new double[GetNstate(site)];
		try	{
			double* tbl = GetCondLikelihood(from);
			const double* stat = GetPath(0,site)->GetStationary();
			// const double* stat = GetBranchSiteSubstitutionProcess(0,site)->GetStationary();
			double tot = 0;
			for (int k=0; k<GetNstate(site); k++)	{
				aux[k] = stat[k] * tbl[k];
				tot += aux[k];
				cumulaux[k] = tot;
			}
			double u = tot * Random::Uniform();
			int s = 0;
			while ((s < GetNstate(site)) && (cumulaux[s] < u))	{
				s ++;
			}
			if (s == GetNstate(site))	{
				cerr << "error in pruning ancestral: overflow\n";
				exit(1);
			}
			state = s;
		}
		catch(...)	{
			cerr << "in root::PruningAncestral\n";
			for (int k=0; k<GetNstate(site); k++)	{
				cerr << aux[k] << '\n';
			}
			exit(1);
			throw;
		}
		delete[] aux;
		delete[] cumulaux;
	}
	for(const Link* link=from->Next(); link!=from; link=link->Next())	{
		if ((flagmap[link->Out()->GetNode()]) && (!isMissing(link->Out()->GetNode(),site)))	{
			double* aux = new double[GetNstate(site)];
			double* cumulaux = new double[GetNstate(site)];
			try	{
				for (int k=0; k<GetNstate(site); k++)	{
					aux[k] = 1;
				}
				GetPath(link->GetBranch(),site)->GetFiniteTimeTransitionProb(state,aux);
				// GetBranchSiteSubstitutionProcess(link->GetBranch(),site)->GetFiniteTimeTransitionProb(state,aux);
				double* tbl = GetCondLikelihood(link->Out());
				for (int k=0; k<GetNstate(site); k++)	{
					aux[k] *= tbl[k];
				}

				// dealing with numerical problems:
				double max = 0;
				for (int k=0; k<GetNstate(site); k++)	{
					if (aux[k] <0)	{
						aux[k] = 0;
					}
					if (max < aux[k])	{
						max = aux[k];
					}
				}
				if (max == 0)	{
					const double* stat = GetPath(0,site)->GetStationary();
					// const double* stat = GetBranchSiteSubstitutionProcess(0,site)->GetStationary();
					for (int k=0; k<GetNstate(site); k++)	{
						aux[k] = stat[k];
					}
				}
				// end of dealing with dirty numerical problems
				double tot = 0;
				for (int k=0; k<GetNstate(site); k++)	{
					tot += aux[k];
					cumulaux[k] = tot;
				}
				double u = tot * Random::Uniform();
				int s = 0;
				while ((s < GetNstate(site)) && (cumulaux[s] < u))	{
					s ++;
				}
				if (s == GetNstate(site))	{
					cerr << "error in pruning ancestral: overflow\n";
					exit(1);
				}
				int& nodestate = GetState(link->Out()->GetNode(),site);
				nodestate = s;
			}
			catch(...)	{
				cerr << "in internal leave::PruningAncestral\n";
				for (int k=0; k<GetNstate(site); k++)	{
					cerr << aux[k] << '\n';
				}
				exit(1);
				throw;
			}
			delete[] aux;
			delete[] cumulaux;
			PruningAncestral(link->Out(),site);
		}
	}
}

double PhyloProcess::RecordPruningAncestralLogProb(const Link* from, int site)	{
	double total = 0;
	if (! flagmap[from->GetNode()])	{
		cerr << "error in pruning ancestral: flag map not set up correctly\n";
		exit(1);
	}
	int& state = GetState(from->GetNode(),site);
	if (from->isRoot())	{
		double* aux = new double[GetNstate(site)];
		try	{
			double* tbl = GetCondLikelihood(from);
			const double* stat = GetPath(0,site)->GetStationary();
			// const double* stat = GetBranchSiteSubstitutionProcess(0,site)->GetStationary();
			double tot = 0;
			for (int k=0; k<GetNstate(site); k++)	{
				aux[k] = stat[k] * tbl[k];
				tot += aux[k];
			}
			total += log(aux[state] / tot);
		}
		catch(...)	{
			cerr << "in root::PruningAncestral\n";
			for (int k=0; k<GetNstate(site); k++)	{
				cerr << aux[k] << '\n';
			}
			exit(1);
			throw;
		}
		delete[] aux;
	}
	for(const Link* link=from->Next(); link!=from; link=link->Next())	{
		if ((flagmap[link->Out()->GetNode()]) && (!isMissing(link->Out()->GetNode(),site)))	{
			double* aux = new double[GetNstate(site)];
			try	{
				for (int k=0; k<GetNstate(site); k++)	{
					aux[k] = 1;
				}
				GetPath(link->GetBranch(),site)->GetFiniteTimeTransitionProb(state,aux);
				// GetBranchSiteSubstitutionProcess(link->GetBranch(),site)->GetFiniteTimeTransitionProb(state,aux);
				double* tbl = GetCondLikelihood(link->Out());
				for (int k=0; k<GetNstate(site); k++)	{
					aux[k] *= tbl[k];
				}

				// dealing with numerical problems:
				double max = 0;
				for (int k=0; k<GetNstate(site); k++)	{
					if (aux[k] <0)	{
						aux[k] = 0;
					}
					if (max < aux[k])	{
						max = aux[k];
					}
				}
				if (max == 0)	{
					const double* stat = GetPath(0,site)->GetStationary();
					// const double* stat = GetBranchSiteSubstitutionProcess(0,site)->GetStationary();
					for (int k=0; k<GetNstate(site); k++)	{
						aux[k] = stat[k];
					}
				}
				// end of dealing with dirty numerical problems
				double tot = 0;
				for (int k=0; k<GetNstate(site); k++)	{
					tot += aux[k];
				}
				int& nodestate = GetState(link->Out()->GetNode(),site);
				total += log(aux[nodestate] / tot);
			}
			catch(...)	{
				cerr << "in internal leave::PruningAncestral\n";
				for (int k=0; k<GetNstate(site); k++)	{
					cerr << aux[k] << '\n';
				}
				exit(1);
				throw;
			}
			delete[] aux;
			total += RecordPruningAncestralLogProb(link->Out(),site);
		}
	}
	return total;
}

void PhyloProcess::drawSample()	{
	if (clampdata) {
		cerr << "mapping substitutions\n";
		ResampleSub();
		cerr << "ok\n";
		cerr << '\n';
	}
	else	{
		cerr << "mapping substitutions\n";
		cerr << '\n';
		PriorSample();
	}

}

void PhyloProcess::RootPosteriorDraw(int site)	{
	double* aux = new double[GetNstate(site)];
	double* tbl = GetCondLikelihood(GetRoot());
	const double* stat = GetPath(0,site)->GetStationary();
	// const double* stat = GetBranchSiteSubstitutionProcess(0,site)->GetStationary();
	for (int k=0; k<GetNstate(site); k++)	{
		aux[k] = stat[k] * tbl[k];
	}
	GetState(GetRoot()->GetNode(),site) = Random::DrawFromDiscreteDistribution(aux,GetNstate(site));
	delete[] aux;
}

void PhyloProcess::PriorSample(const Link* from, int site, bool rootprior)	{
	int& state = GetState(from->GetNode(),site);
	if (from->isRoot())	{
		if (rootprior)	{
			state = GetPath(0,site)->DrawStationary();
			// state = GetBranchSiteSubstitutionProcess(0,site)->DrawStationary();
		}
		else	{
			RootPosteriorDraw(site);
		}
	}
	for(const Link* link=from->Next(); link!=from; link=link->Next())	{
		if (!isMissing(link->GetBranch(),site))	{
			GetState(link->Out()->GetNode(),site) = GetPath(link->GetBranch(),site)->DrawFiniteTime(state);
			// GetState(link->Out()->GetNode(),site) = GetBranchSiteSubstitutionProcess(link->GetBranch(),site)->DrawFiniteTime(state);
			PriorSample(link->Out(),site,rootprior);
		}
	}
}

void PhyloProcess::ResampleState()	{
	for (int i=0; i<GetNsite(); i++)	{
		ResampleState(i);
	}
}

void PhyloProcess::ResampleState(int site)	{
	if (! isMissing(GetRoot()->GetNode(),site))	{
		Pruning(GetRoot(),site);
		// FastSiteLogLikelihood(site);
		PruningAncestral(GetRoot(),site);
	}
}

void PhyloProcess::ResampleSub()	{

	pruningchrono.Start();
	for (int i=0; i<GetNsite(); i++)	{
		if (sitearray[i])	{
			if (! isMissing(GetRoot()->GetNode(),i))	{
				ResampleState(i);
			}
		}
	}
	pruningchrono.Stop();

	resamplechrono.Start();
	for (int i=0; i<GetNsite(); i++)	{
		if (sitearray[i])	{
			if (! isMissing(GetRoot()->GetNode(),i))	{
				ResampleSub(GetRoot(),i);
			}
		}
	}
	resamplechrono.Stop();
}

void PhyloProcess::ResampleSub(int site)	{
	if (! isMissing(GetRoot()->GetNode(),site))	{
		ResampleState(site);
		ResampleSub(GetRoot(),site);
	}
}

void PhyloProcess::ResampleSub(const Link* from, int site)	{
	if (from->isRoot())	{
		GetPath(0,site)->ResampleRoot(GetState(from->GetNode(),site));
	}
	for (const Link* link=from->Next(); link!=from; link=link->Next())	{
		if (! isMissing(link->Out(),site))	{
			GetPath(link->GetBranch(),site)->Resample(GetState(link->GetNode(),site),GetState(link->Out()->GetNode(),site));
			/*
			if (branchmap)	{
				GetPath(link->GetBranch(),site)->Resample(GetState(link->GetNode(),site),GetState(link->Out()->GetNode(),site));
			}
			else	{
				GetPath(link->GetBranch(),site)->ResampleEnds(GetState(link->GetNode(),site),GetState(link->Out()->GetNode(),site));
			}
			*/
			ResampleSub(link->Out(),site);
		}
	}
}

/*********
/ MH moves
********/

void PhyloProcess::SetProposalMatrices()	{
	RecursiveSetProposalMatrices(GetRoot());
}

void PhyloProcess::RecursiveSetProposalMatrices(const Link* from)	{

	for (int i=0; i<GetNsite(); i++)	{
		if (! isMissing(from,i))	{
			GetPath(from->GetBranch(),i)->SetProposalMatrix(GetProposalMatrix(from->GetBranch(),i));
		}
	}
	for (const Link* link=from->Next(); link!=from; link=link->Next())	{
		RecursiveSetProposalMatrices(link->Out());
	}
}

void PhyloProcess::RecursiveRegister(const Link* from, int site, Mnode* mnode)	{

	if (! isMissing(from->GetBranch(),site))	{
		GetPath(from->GetBranch(),site)->Register(mnode);
		if (flagmap[from->GetNode()])	{
			for (const Link* link=from->Next(); link!=from; link=link->Next())	{
				RecursiveRegister(link->Out(),site,mnode);
			}
		}
	}
}

void PhyloProcess::SwapMatrices(const Link* from, int site)	{

	if (! isMissing(from->GetBranch(),site))	{
		GetPath(from->GetBranch(),site)->SwapMatrices();
		if (flagmap[from->GetNode()])	{
			for (const Link* link=from->Next(); link!=from; link=link->Next())	{
				SwapMatrices(link->Out(),site);
			}
		}
	}
}

double PhyloProcess::ProposeResampleSub(const Link* from, int site)	{

	if (! isMissing(from->GetBranch(),site))	{
		double tot = 0;
		if (from->isRoot())	{
			tot = GetPath(from->GetBranch(),site)->ProposeMoveRoot(GetState(from->GetNode(),site));
		}
		if (flagmap[from->GetNode()])	{
			for (const Link* link=from->Next(); link!=from; link=link->Next())	{
				if (! isMissing(link->Out(),site))	{
					tot += GetPath(link->GetBranch(),site)->ProposeMove(GetState(link->GetNode(),site),GetState(link->Out()->GetNode(),site));
					tot += ProposeResampleSub(link->Out(),site);
				}
			}
		}
		return tot;
	}
	return 0;
}

double PhyloProcess::MHMove(int nrep, double s01, double s10)	{

	int nacc = 0;
	for (int rep=0; rep<nrep; rep++)	{
		// reset flags
		// choose a subset of nodes rootlink
		int site = (int) (GetNsite() * Random::Uniform());
		while (isMissing(GetRoot()->GetNode(),site))	{
			site = (int) (GetNsite() * Random::Uniform());
		}

		const Link* rootlink;
		if (s01 == -1)	{
			ResetFlagMap(GetRoot(),true);
			rootlink = GetRoot();
		}
		else	{
			ResetFlagMap(GetRoot(),false);
			int sw = 0;
			rootlink = ChooseNodeSet(GetRoot(),s01,s10,sw);
		}

		// choose a site
		Mnode* mnode = new Mnode;
		RecursiveRegister(rootlink,site,mnode);

		mnode->Corrupt(true);

		BackupNodeStates(rootlink,site);

		SwapMatrices(rootlink,site);

		Pruning(rootlink,site);
		// double loghastings = 0;

		// loghastings += RecordPruningAncestralLogProb(rootlink,site);

		PruningAncestral(rootlink,site);

		// loghastings -= RecordPruningAncestralLogProb(rootlink,site);

		double loghastings = ProposeResampleSub(rootlink,site);

		SwapMatrices(rootlink,site);

		double logratio = mnode->Update() + loghastings;

		bool accepted = (log(Random::Uniform()) < logratio);
		if (accepted)	{
			nacc++;
		}
		else	{
			RestoreNodeStates(rootlink,site);
			mnode->Corrupt(false);
			mnode->Restore();
		}
		delete mnode;
	}
	return ((double) nacc) / nrep;
}

void PhyloProcess::ResetFlagMap(const Link* from, bool in)	{
	flagmap[from->GetNode()] = in;
	for (const Link* link=from->Next(); link!=from; link=link->Next())	{
		ResetFlagMap(link->Out(),in);
	}
}

const Link* PhyloProcess::ChooseNodeSet(const Link* from, double s01, double s10, int sw)	{

	const Link* ret = 0;
	if (sw == 0)	{
		if (Random::Uniform() < s01)	{
			sw = 1;
			ret = from;
			flagmap[from->GetNode()] = true;
		}
	}
	else if (sw == 1)	{
		if (Random::Uniform() < s10)	{
			sw = 2;
		}
	}
	if (sw < 2)	{
		for (const Link* link=from->Next(); link!=from; link=link->Next())	{
			const Link* tmp = ChooseNodeSet(link->Out(),s01,s10,sw);
			if (! ret && tmp)	{
				ret = tmp;
			}
		}
	}
	return ret;
}


void PhyloProcess::PriorSample()	{

	for (int i=0; i<GetNsite(); i++)	{
		if (! isMissing(GetRoot()->GetNode(),i))	{
			PriorSample(GetRoot(),i,true);
		}
	}
}

void PhyloProcess::PostPredSample(bool rootprior)	{
	for (int i=0; i<GetNsite(); i++)	{
		PostPredSample(i,rootprior);
	}
}

void PhyloProcess::PostPredSample(int site, bool rootprior)	{
	// why pruning?
	if (! isMissing(GetRoot()->GetNode(),site))	{
		Pruning(GetRoot(),site);
		PriorSample(GetRoot(),site,rootprior);
	}
}

double PhyloProcess::Move(double tuning)	{
	for (int i=0; i<GetNsite(); i++)	{
		sitearray[i] = (Random::Uniform() < tuning);
	}
	ResampleSub();
	return 1;
}

void PhyloProcess::GetLeafData(SequenceAlignment* data)	{
	RecursiveGetLeafData(GetRoot(),data);
}

void PhyloProcess::RecursiveGetLeafData(const Link* from, SequenceAlignment* data)	{

	if (from->isLeaf())	{
		for (int site=0; site<GetNsite(); site++)	{
			int state = GetState(from->GetNode(),site);
			int obsstate = GetData(from->GetNode()->GetIndex(),site);
			if (obsstate != unknown)	{
				data->SetState(from->GetNode()->GetIndex(),site,state);
			}
			else	{
				data->SetState(from->GetNode()->GetIndex(),site,unknown);
			}
		}
	}
	for (const Link* link=from->Next(); link!=from; link=link->Next())	{
		RecursiveGetLeafData(link->Out(),data);
	}
}


void PhyloProcess::GetLeafFreqs(const Link* from, double** taxfreq)	{

	if (from->isLeaf())	{
		for (int site=0; site<GetNsite(); site++)	{
			int state = GetState(from->GetNode(),site);
			int obsstate = GetData(from->GetNode()->GetIndex(),site);
			if (obsstate != unknown)	{
				taxfreq[from->GetNode()->GetIndex()][state]++;
			}
		}
	}
	for (const Link* link=from->Next(); link!=from; link=link->Next())	{
		GetLeafFreqs(link->Out(),taxfreq);
	}
}

double PhyloProcess::CompositionalHeterogeneityIndex(ostream& os)	{

	double** taxfreq = new double*[GetNtaxa()];
	for (int j=0; j<GetNtaxa(); j++)	{
		taxfreq[j] = new double[GetMaxNstate()];
		for (int k=0; k<GetMaxNstate(); k++)	{
			taxfreq[j][k] = 0;
		}
	}

	// recursive filling of the array of frequencies
	GetLeafFreqs(GetRoot(),taxfreq);

	// make global freqs out of tax-specific freqs
	double* globalfreq = new double[GetMaxNstate()];
	for (int k=0; k<GetMaxNstate(); k++)	{
		globalfreq[k] = 0;
		for (int j=0; j<GetNtaxa(); j++)	{
			globalfreq[k] += taxfreq[j][k];
		}
	}

	// normalise
	double total = 0;
	for (int k=0; k<GetMaxNstate(); k++)	{
		total += globalfreq[k];
	}
	for (int k=0; k<GetMaxNstate(); k++)	{
		globalfreq[k] /= total;
	}
	for (int j=0; j<GetNtaxa(); j++)	{
		double total = 0;
		for (int k=0; k<GetMaxNstate(); k++)	{
			total += taxfreq[j][k];
		}
		for (int k=0; k<GetMaxNstate(); k++)	{
			taxfreq[j][k] /= total;
			os << taxfreq[j][k] << '\t';
		}
		os << '\n';
	}
	os << '\n';

	// compute max distance
	double maxdist = 0;
	for (int j=0; j<GetNtaxa(); j++)	{
		double dist = 0;
		for (int k=0; k<GetMaxNstate(); k++)	{
			double tmp = (taxfreq[j][k] - globalfreq[k]);
			dist += tmp * tmp;
		}
		if (maxdist < dist)	{
			maxdist = dist;
		}
	}
	return maxdist;
}

