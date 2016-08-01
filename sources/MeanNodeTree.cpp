
#include "FileDrawTree.h"

#include <fstream>
#include <sstream>

class MeanNodeTree : public Tree	{

	public:

	MeanNodeTree(string filename, string list) : Tree(filename), flag(false)	{
		ifstream is(list.c_str());
		is >> size;
		for (int i=0; i<size; i++)	{
			string name;
			is >> name;
			cerr << name << '\n';
			FileTree tree(name);
			RecursiveAdd(&tree, tree.GetRoot());
		}

		RecursiveNormalise(GetRoot());
		flag = true;

	}

	void RecursiveNormalise(const Link* from)	{
		nodeval[from->GetNode()] /= size;
		for (const Link* link = from->Next(); link != from; link=link->Next())	{
			RecursiveNormalise(link->Out());
		}
	}

	void RecursiveAdd(FileTree* tree, const Link* from)	{
		const Link* leaf1 = tree->GetLeftMostLink(from);
		const Link* leaf2 = tree->GetRightMostLink(from);
		string name1 = tree->GetLeafNodeName(leaf1);
		string name2 = tree->GetLeafNodeName(leaf2);
		const Link* link = GetLCA(name1,name2);
		if (! link)	{
			cerr << "error in mean node tree: did not find common ancestor of " << name1 << " and " << name2 << '\n';
			ToStream(cerr);
			exit(1);
		}
		nodeval[link->GetNode()] += tree->GetNodeVal(from);
		for (const Link* link = from->Next(); link != from; link=link->Next())	{
			RecursiveAdd(tree,link->Out());
		}
	}

	string GetNodeName(const Link* link)	const {
		if (! flag)	{
			return Tree::GetNodeName(link);
		}
		ostringstream s;
		if (link->isLeaf())	{
			s << Tree::GetNodeName(link) << "_" << nodeval.find(link->GetNode())->second;
		}
		else	{
			s << nodeval.find(link->GetNode())->second;
		}
		return s.str();
	}

	string GetBranchName(const Link* link)	{
		return Tree::GetBranchName(link);
	}

	protected:

	bool flag;
	map<const Node*, double> nodeval;
	int size;
};


int main(int argc, char* argv[])	{

	string ref = argv[1];
	string list = argv[2];
	ofstream os(argv[3]);

	MeanNodeTree* nodetree = new MeanNodeTree(ref,list);
	nodetree->ToStream(os);
	os.close();
}

