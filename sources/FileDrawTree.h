
#ifndef FILEDRAWTREE_H
#define FILEDRAWTREE_H

#include "DrawTree.h"

#include "StringStreamUtils.h"
#include "BiologicalSequences.h"
#include <regex>


static const string number = "\\-?\\d*(\\.\\d+)?(e\\+\\d+)?";

void FileTree::Check(const Link* from)	{

	string name = from->GetNode()->GetName();
	regex re;
	int minmax = 2;
	string s = from->isLeaf() ? ".+_" : "";
	s += "(" + number + ")_(" + number + ")";
	re.assign(s,regex_constants::icase);
	if (!regex_match(name,re))	{
		minmax = 1;
		string s = from->isLeaf() ? ".+_" : "";
		s += "(" + number + ")";
		re.assign(s,regex_constants::icase);
		if (!regex_match(name,re))	{
			minmax = 0;
		}
	}
	if (from->isLeaf())	{
		if (leafval == -1)	{
			leafval = minmax;
		}
		else if (leafval != minmax)	{
			cerr << "error : not all leaves have same format\n";
			cerr << leafval << '\t' << minmax << '\n';
			cerr << name << '\n';
			exit(1);
		}
		else	{
			leafval = minmax;
		}
	}
	else	{
		if (internalval == -1)	{
			internalval = minmax;
		}
		else if (internalval != minmax)	{
			cerr << "error : not all internal nodes have same format\n";
			cerr << internalval << '\t' << minmax << '\n';
			cerr << "node name : " << name << '\n';
			exit(1);
		}
		else	{
			internalval = minmax;
		}
	}
	for (const Link* link=from->Next(); link!=from; link=link->Next())	{
		Check(link->Out());
	}
}

string FileTree::GetLeafNodeName(const Link* from){

	if (!from->isLeaf())	{
		cerr << "error :get leaf node name\n";
		exit(1);
	}
	string name = from->GetNode()->GetName();
	string s = "(.+)";
	regex re;
	if (HasLeafMinMax())	{
		s += "(" + number + ")_(" + number + ")";
	}
	else if (HasLeafVal())	{
		s += "_(" + number + ")";
	}
	else	{
	}
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

double FileTree::GetNodeVal(const Link* from)	{

	string name = from->GetNode()->GetName();
	string s = from->isLeaf() ? ".+_" : "";
	regex re;
	s += "(" + number + ")";
	re.assign(s,regex_constants::icase);
	cmatch c;
	if (!regex_match(name.c_str(),c,re))	{
		cerr << "error in filetree\n";
		exit(1);
	}
	string rets(c[1].first,c[1].second);
	double ret = atof(rets.c_str());
	if (rescale)	{
		ret /= rescale;
	}
	return ret;
}


double FileTree::GetNodeMin(const Link* from)	{

	string name = from->GetNode()->GetName();
	string s = from->isLeaf() ? ".+_" : "";
	regex re;
	s += "(" + number + ")_(" + number + ")";
	re.assign(s,regex_constants::icase);
	cmatch c;
	if (!regex_match(name.c_str(),c,re))	{
		cerr << "error in filetree\n";
		exit(1);
	}
	string rets(c[1].first,c[1].second);
	double ret = atof(rets.c_str());
	if (rescale)	{
		ret /= rescale;
	}
	return ret;
}

double FileTree::GetNodeMax(const Link* from)	{

	string name = from->GetNode()->GetName();
	string s = from->isLeaf() ? ".+_" : "";
	regex re;
	s += "(" + number + ")_(" + number + ")";
	re.assign(s,regex_constants::icase);
	cmatch c;
	if (!regex_match(name.c_str(),c,re))	{
		cerr << "error in filetree\n";
		exit(1);
	}
	string rets(c[4].first,c[4].second);
	double ret = atof(rets.c_str());
	if (rescale)	{
		ret /= rescale;
	}
	return ret;
}

#endif
