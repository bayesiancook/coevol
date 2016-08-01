
class BranchSumStat {

	public:

	BranchSumStat(StateSpace* instatespace)	: statspace(instatespace) {}

	// accessors

	protected:
	map<int,int> rootcount;
	map< pair<int,int>, int> paircount;

};

class NodeSumStat	{

	NodeSumStat(StateSpace* instatespace)	: statspace(instatespace) {}

	// accessors

	Tree* tree;
	StateSpace* statespace;

	map<const Node*, map<int,int> > rootcount;

};

class BranchSumStatTree	{

	BranchSumStatTree(const Tree* intree, PhyloProcess* inprocess)	{

	}

	map<const Node*, NodeSumStat> nodemap;
};

class NodeSumStatTree	{


	map<const Branch* BranchSumStat> branchmap;
};


class CodonBranchSumStat : public BranchSumStat	{

};

class CodonNodeSumStat : public NodeSumStat	{

};

class CodonBranchSumStatTree	{


};


