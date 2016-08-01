
// #include "CompoTree.h"
#include "FileDrawTree.h"

class DrawCompTree : public BubbleTree {

	public:

	DrawCompTree(Tree* intree1, Tree* intree2)	{
		tree1 = intree1;
		tree2 = intree2;
		cerr << "recursive set val\n";
		RecursiveSetVal(tree1->GetRoot(),tree2->GetRoot());
	}

	const Link* GetRoot() const {
		return tree1->GetRoot();
	}

	void RecursiveSetVal(const Link* from1, const Link* from2)	{

		string name1 = from1->GetNode()->GetName();
		string name2 = from2->GetNode()->GetName();
		string leafname1 = "";
		string leafname2 = "";
		if (from1->isLeaf())	{
			leafname1 = GetValFromLeafName(name1);
			leafname2 = GetValFromLeafName(name2);
		}
		else	{
			leafname1 = name1;
			leafname2 = name2;
		}
		// cerr << leafname1 << '\t' << leafname2 << '\n';
		val1[from1->GetNode()] = atof(leafname1.c_str());
		val2[from1->GetNode()] = atof(leafname2.c_str());

		const Link* link1 = from1->Next();
		const Link* link2 = from2->Next();
		while (link1 != from1)	{
			RecursiveSetVal(link1->Out(),link2->Out());
			link1 = link1->Next();
			link2 = link2->Next();
		}
	}

	double GetNodeVal1(const Link* from)	{
		return val1[from->GetNode()];
	}

	double GetNodeVal2(const Link* from)	{
		return val2[from->GetNode()];
	}
	/*
	double GetNodeVal1(const Link* from)	{
		cerr << from << '\n';
		string tmp = from->GetNode()->GetName();
		cerr << tmp << '\n';
		double temp = atof(tmp.c_str());
		if (temp <= 0)	{
			cerr << "error: value should be positive\n";
			cerr << tmp << '\t' << temp << '\n';
			exit(1);
		}
		return temp;
	}

	double GetNodeVal1(const Link* from)	{

		string name = from->GetNode()->GetName();
		string s = from->isLeaf() ? ".+_" : "";
		regex re;
		s += "(\\d*(\\.\\d+)?(e\\+\\d+)?)";
		re.assign(s,regex_constants::icase);
		cmatch c;
		if (!regex_match(name.c_str(),c,re))	{
			cerr << "error in filetree\n";
			exit(1);
		}
		string rets(c[1].first,c[1].second);
		double ret = atof(rets.c_str());
		if (ret <= 0)	{
			cerr << "error : " << ret << '\n';
			exit(1);
		}
		return ret;

	}

	double GetNodeVal2(const Link* from)	{

		string name1 = tree1->GetLeftMost(from);
		string name2 = tree1->GetRightMost(from);
		string leafname1 = ParseLeafName(name1);
		string leafname2 = ParseLeafName(name2);

		cerr << name1 << '\t' << name2 << '\t' << leafname1 << '\t' << leafname2 << '\n';

		const Link* from2 = tree2->GetLCA(leafname1,leafname2);
		string name = from2->GetNode()->GetName();
		cerr << name << '\n';

		string s = from->isLeaf() ? ".+_" : "";
		regex re;
		s += "(\\d*(\\.\\d+)?(e\\+\\d+)?)";
		re.assign(s,regex_constants::icase);
		cmatch c;
		if (!regex_match(name.c_str(),c,re))	{
			cerr << "error in filetree\n";
			exit(1);
		}
		string rets(c[1].first,c[1].second);
		double ret = atof(rets.c_str());
		if (ret <= 0)	{
			cerr << "error : " << ret << '\n';
			exit(1);
		}
		return ret;

	}

	double GetNodeVal2(const Link* from)	{
		cerr << from << '\n';
		cerr << tree1->GetLeftMost(from) << '\t' << tree1->GetRightMost(from) << '\n';
		const Link* from2 = tree2->GetLCA(tree1->GetLeftMost(from),tree1->GetRightMost(from));
		string tmp = from2->GetNode()->GetName();
		double temp = atof(tmp.c_str());
		if (temp <= 0)	{
			cerr << "error: value should be positive\n";
			cerr << tmp << '\t' << temp << '\n';
			exit(1);
		}
		return temp;
	}
	*/

	double GetNodeMin(const Link* from)	{
		double tmp1 = GetNodeVal1(from);
		double tmp2 = GetNodeVal2(from);
		return (tmp1 < tmp2) ? tmp1 : tmp2;
	}

	double GetNodeMax(const Link* from)	{
		double tmp1 = GetNodeVal1(from);
		double tmp2 = GetNodeVal2(from);
		return (tmp1 > tmp2) ? tmp1 : tmp2;
	}

	double GetNodeVal(const Link* from)	{
		double tmp1 = GetNodeVal1(from);
		double tmp2 = GetNodeVal2(from);
		return 0.5 * (tmp1 + tmp2);
	}

	double GetLength(const Link* from) const {
		if (! from->GetBranch())	{
			return 0;
		}
		string tmp = from->GetBranch()->GetName();
		double temp = atof(tmp.c_str());
		if (temp <= 0)	{
			cerr << "error : bl is not positive\n";
			cerr << tmp << '\t' << temp << '\n';
			exit(1);
		}
		return temp;
	}


	string GetLeafNodeName(const Link* from)	const {
		if (! from->isLeaf())	{
			cerr << "error: get leaf node name called on non leaf link\n";
			exit(1);
		}
		string name = from->GetNode()->GetName();
		return ParseLeafName(name);
	}

	string GetValFromLeafName(string name) const {
		string s = ".+";
		s += "_(\\d*(\\.\\d+)?(e\\+\\d+)?)";
		regex re;
		re.assign(s,regex_constants::icase);
		cmatch c;
		if (!regex_match(name.c_str(),c,re))	{
			cerr << "error in filetree\n";
			cerr << "get leaf node name\n";
			cerr << name << '\n';
			exit(1);
		}
		string rets(c[1].first,c[1].second);
		return StringReplace('_'," ",rets);
	}

	string ParseLeafName(string name) const {
		string s = "(.+)";
		s += "_(\\d*(\\.\\d+)?(e\\+\\d+)?)";
		regex re;
		re.assign(s,regex_constants::icase);
		cmatch c;
		if (!regex_match(name.c_str(),c,re))	{
			cerr << "error in filetree\n";
			cerr << "get leaf node name\n";
			cerr << name << '\n';
			exit(1);
		}
		string rets(c[1].first,c[1].second);
		return StringReplace('_'," ",rets);
	}

	string GetNodeName(const Link* from)	const	{
		if (from->isLeaf())	{
			return GetLeafNodeName(from);
		}
		return "";
	}

	protected:

	Tree* tree1;
	Tree* tree2;
	map<const Node*,double> val1;
	map<const Node*,double> val2;
};


int main(int argc, char* argv[])	{

	string file1 = "";
	string file2 = "";
	string out = "";
	double nodescale = 5;
	double nodepower = 1;
	double x = 6;
	double y = 8;
	// double x = 10;
	// double y = 17;
	bool sizeset = false;
	bool minmax = false;
	int withnodetext = 0;
	int prec = 2;
	bool withheader = true;
	bool unibranch = false;

	// double fontsize = 8;
	double fontsize = 6;
	// string style = "article";
	string style = "beamer";
	double xoffset = 0;
	bool fontset = false;

	int Nstate = Nnuc;
	double thickness = 0.06;

	try	{

		if (argc == 1)	{
			throw(0);
		}

		int i = 1;
		while (i < argc)	{
			string s = argv[i];

			if (s == "-ns")	{
				i++;
				nodescale = atof(argv[i]);
			}
			else if (s == "-np")	{
				i++;
				nodepower= atof(argv[i]);
			}
			else if (s == "-fs")	{
				i++;
				fontsize = atof(argv[i]);
				fontset = true;
			}
			else if (s == "-aa")	{
				Nstate = Naa;
			}
			else if (s == "-branch")	{
				unibranch = true;
			}
			else if (s == "-th")	{
				i++;
				thickness = atof(argv[i]);
			}
			else if (s == "-xoff")	{
				i++;
				xoffset = atof(argv[i]);
			}
			else if (s == "-x")	{
				i++;
				x = atof(argv[i]);
				sizeset = true;
			}
			else if (s == "-y")	{
				i++;
				y = atof(argv[i]);
				sizeset = true;
			}
			else if (s == "-beamer")	{
				style = "beamer";
				if (! sizeset)	{
					x = 6;
					y = 8;
				}
				if (! fontset)	{
					fontsize = 4;
				}
			}
			else if (s == "-article")	{
				style = "article";
				if (! sizeset)	{
					x = 15;
					y = 25;
				}
				if (! fontset)	{
					fontsize = 8;
				}
			}
			else if (s == "-o")	{
				i++;
				out = argv[i];
			}
			else if (s == "-h")	{
				withheader = false;
			}
			else if (s == "-mm")	{
				minmax = true;
			}
			else if (s == "-nt")	{
				withnodetext = 1;
			}
			else if (s == "-lt")	{
				withnodetext = 2;
			}
			else if (s == "-p")	{
				i++;
				prec = atoi(argv[i]);
			}
			else	{
				if (i != (argc -2))	{
					throw(0);
				}
				file1 = argv[i];
				i++;
				file2 = argv[i];
			}
			i++;
		}
		if (file1 == "")	{
			throw(0);
		}
		if (file2 == "")	{
			throw(0);
		}
	}
	catch(...)	{
		cerr << file1 << '\t' << file2 << '\n';
		cerr << "drawtree [-x <xsize> -y <ysize> -o <outfile> -ns <node_val_scale> -mm] <inputfile> <inputfile2> \n";
		cerr << '\n';
		exit(1);
	}

	if (out == "")	{
		out = "out";
	}

	Tree* tree1 = new Tree(file1);
	Tree* tree2 = new Tree(file2);
	tree1->ToStream(cerr);
	tree2->ToStream(cerr);
	DrawCompTree* tree = new DrawCompTree(tree1,tree2);
	tree->SetDocumentStyle(style);
	tree->SetScale(x,y);
	tree->SetFontSize(fontsize);
	tree->SetNodeScale(nodescale);
	tree->SetNodePower(nodepower);
	cerr << "draw\n";
	tree->Draw(out);
	cerr << "ok\n";
}

