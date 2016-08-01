#ifndef PROFILEINFINITEMIXTUREPHYLOPROCESS_H
#define PROFILEINFINITEMIXTUREPHYLOPROCESS_H


class ProfileInfiniteMixturePhyloProcess : public MatrixInfiniteMixturePhyloProcess<Profile>	{

	public:

	ProfileInfiniteMixturePhyloProcess(LengthTree* intree, MatrixInfiniteMixture<Profile>* inmatmix,  CodonSequenceAlignment* indata) : MatrixInfiniteMixturePhyloProcess<Profile>(intree,inmatmix,indata)	{data=indata;}

	int BranchSiteNonsynSubCount(Link* from, int site)      {
		int count = 0;
		for (Link* blink=from->Next(); blink!=from; blink=blink->Next())        {
			Plink* plink = GetPath(blink->GetBranch(), site)->Init();
			while (plink != GetPath(blink->GetBranch(), site)->Last())      {
				if (! data->GetCodonStateSpace()->Synonymous(plink->GetState(), plink->Next()->GetState()))     {
					count++;
				}
				plink=plink->Next();
			}
			count += BranchSiteNonsynSubCount(blink->Out(), site);
		}
		return count;
        }

	double NonsynSubMean()    {
		int count = 0;
		for (int site=0; site<GetNsite(); site++)       {
			count += BranchSiteNonsynSubCount(GetRoot(), site);
		}
		return ( (double)(count)/GetNsite() );
	}

	double NonsynSubVariance()	{
		int count = 0;
		int sum = 0;
		double var = 0;
		for (int site=0; site<GetNsite(); site++)	{
			count = BranchSiteNonsynSubCount(GetRoot(), site);
			var += (count * count);
			sum += count;
		}
		double mean = (double)(sum)/GetNsite();
		var /= GetNsite();
		var -= (mean * mean);
		return var;
	}

	private:
	CodonSequenceAlignment* data;

};


#endif
