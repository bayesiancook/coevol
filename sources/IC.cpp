
#include "Tree.h"
#include "ContinuousData.h"
#include "dagostino.h"
#include <vector>

double ComputeIndependentContrasts(const Link* from, ContinuousData* data, int index, vector<double>& ic, Tree* tree, double& hm)	{

	if (from->isLeaf())	{
		int tax = data->GetTaxonSet()->GetTaxonIndex(from->GetNode()->GetName());
		if (tax == -1)	{
			cerr << "error : " << from->GetNode()->GetName() << '\n';
			exit(1);
		}
		double tmp = data->GetState(tax, index);
		hm = 0;
		return log(tmp);
	}
	else	{
		double hm1 = 0;
		double hm2 = 0;
		double tmp1 = ComputeIndependentContrasts(from->Next()->Out(),data,index,ic,tree,hm1);
		double tmp2 = ComputeIndependentContrasts(from->Next()->Next()->Out(),data,index,ic,tree,hm2);

		string length1 = from->Next()->GetBranch()->GetName();
		string length2 = from->Next()->Next()->GetBranch()->GetName();
		double l1 = atof(length1.c_str()) + hm1;
		double l2 = atof(length2.c_str()) + hm2;

		cout << (tmp1 - tmp2) / sqrt(l1 + l2) << '\n';
		// cout << tree->GetLeftMost(from) << '\t' << tree->GetRightMost(from) << '\t' <<  l1 << '\t' << l2 << '\t' <<(tmp2 - tmp1) / sqrt(l1 + l2) << '\n';
		ic.push_back((tmp1 - tmp2) / sqrt(l1 + l2));
		double tot = (tmp1 * l2 + tmp2 * l1) / (l1 + l2);
		hm = l1 * l2 / (l1 + l2);
		return tot;
	}

		/*
		double tot = 0;
		double totweight = 0;
		for (const Link* link=from->Next(); link!=from; link=link->Next())	{
			double tmp = ComputeIndependentContrast(link->Out(),index);
			string length = link->GetBranch()->GetName();
			double l = atof(length.c_str());
			cerr << length << '\t' << l << '\n';
			double weight = 1.0 / l;
			tot += weight * tmp;
			totweight += weight;
		}
		return tot / totweight;
		*/

	return 0;
}

int main(int argc, char* argv[])	{

	string datafile = argv[1];
	string treefile = argv[2];

	Tree* tree = new Tree(treefile);
	ContinuousData* contdata = new FileContinuousData(datafile);

	const TaxonSet* taxset = contdata->GetTaxonSet();
	tree->RegisterWith(taxset);
	tree->ToStream(cerr);

	vector<double> ic1;
	double hm = 0;
	ComputeIndependentContrasts(tree->GetRoot(),contdata,0,ic1,tree,hm);
	vector<double> ic2;
	ComputeIndependentContrasts(tree->GetRoot(),contdata,1,ic2,tree,hm);

	cout << "ic\n";
	for (int i=0; i<ic1.size(); i++)	{
		cout << ic1[i] << '\t' << ic2[i] << '\n';
	}

	cerr << "vector size: " << ic1.size() << '\n';
	cerr << "dago\n\n";
	cerr << DAgostinosKandZ(ic1) << '\n';
}

