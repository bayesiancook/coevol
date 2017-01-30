
#include <list>
#include <sstream>
#include <fstream>
#include "Tree.h"

bool NewickTree::simplify = false;

void NewickTree::ToStream(ostream& os) const {
	if (simplify)	{
		ToStreamSimplified(os,GetRoot());
	}
	else	{
		ToStream(os,GetRoot());
	}
	os << ";\n";
}

double NewickTree::ToStreamSimplified(ostream& os, const Link* from) const {

	if (!from->isLeaf())	{
		if (from->Next()->Next() == from)	{
			double tot = ToStreamSimplified(os,from->Next()->Out());
			tot += atof(GetBranchName(from).c_str());
			return tot;
		}
		else	{
			os << '(';
			for (const Link* link=from->Next(); link!=from; link=link->Next())	{
				double tmp = ToStreamSimplified(os,link->Out());
				os << ':' << tmp;
				if (link->Next() != from)	{
					os << ',';
				}
			}
			os << ')';
		}
	}
	else	{
	}
	os << GetNodeName(from);
	/*
	if (!from->isRoot())	{
		string brval = GetBranchName(from);
		if (brval != "")	{
			os << ':' << brval;
		}
	}
	*/
	if (from->isRoot())	{
		return 0;
	}
	return atof(GetBranchName(from).c_str());
}

void NewickTree::ToStream(ostream& os, const Link* from) const {

	if (!from->isLeaf())	{
		os << '(';
		for (const Link* link=from->Next(); link!=from; link=link->Next())	{
			ToStream(os,link->Out());
			if (link->Next() != from)	{
				os << ',';
			}
		}
		os << ')';
	}
	// if (from->isLeaf())	{
		os << GetNodeName(from);
	// }
	if (!from->isRoot())	{
		string brval = GetBranchName(from);
		if (brval != "")	{
			os << ':' << brval;
		}
	}
}

Tree::Tree()	{
	root = 0;
	taxset = 0;
}

Tree::Tree(const Tree* from)	{
	taxset = 0;
	root = new Link(from->root);
	root->InsertOut(root);
	RecursiveClone(from->root,root);
}

void Tree::RecursiveClone(const Link* from, Link* to)	{
	Node* node = new Node(from->GetNode());
	to->SetNode(node);
	const Link* linkfrom = from->Next();
	Link* linkto = to;
	while (linkfrom != from)	{
		Link* newnext = new Link(linkfrom); // newnext points to same node and branch as linkfrom
		newnext->SetNode(node);
		linkto->Insert(newnext);
		Link* newout = new Link(linkfrom->Out()); // idem, same node and branch as linkfrom->Out()
		newout->InsertOut(newnext);
		Branch* branch = new Branch(linkfrom->GetBranch());
		newnext->SetBranch(branch);
		newout->SetBranch(branch);
		RecursiveClone(linkfrom->Out(),newout);
		linkfrom = linkfrom->Next();
		linkto = linkto->Next();
	}
}


void Tree::RecursiveDelete(Link* from)	{
	Link* link = from->Next();
	while (link != from)	{
		delete link->Out()->GetNode();
		delete link->GetBranch();
		RecursiveDelete(link->Out());
		Link* keep = link->Next();
		delete link;
		link = keep;
	}
	delete link;
}

/*
void Tree::RecursiveCreateTBL(Link* from, int Nstate)	{
	for(Link* link=from->Next(); link!=from; link=link->Next())	{
		delete[] link->tbl;
		link->tbl = new double[Nstate+1];
		RecursiveCreateTBL(link->Out(), Nstate);
	}
	delete[] from->tbl;
	from->tbl = new double[Nstate+1];
}

void Tree::RecursiveDeleteTBL(Link* from)	{
	for(Link* link=from->Next(); link!=from; link=link->Next())	{
		RecursiveDeleteTBL(link->Out());
		delete[] link->tbl;
	}
	delete[] from->tbl;
}
*/

Tree::~Tree()	{
	if (root)	{
		RecursiveDelete(root);
		root = 0;
	}
}


void Tree::DeleteNextLeaf(Link* previous){
	Link* link = previous->Next();
	if(!link->Out()->isLeaf()){
		cout << "Bad call of DeleteNextLeaf, it must be call on a link pointing on a leaf\n";
		exit(1);
	}
	previous->Next()->Next()->AppendTo(previous);
	delete link->Out()->GetNode();
	delete link->Out();
	delete link->GetBranch();
	delete link;
}

void Tree::DeleteUnaryNode(Link* from){
	if(!from->isUnary()){
		cout << "Bad call of DeleteUnaryNode, node is not unary\n";
		exit(1);
	}
	if(from->isRoot()){
		Link* newroot = from->Next()->Out();
		newroot->SetBranch(from->GetBranch());
		delete from->Next()->GetBranch();
		delete from->GetNode();
		delete from->Next();
		delete from;
		root = newroot;
		root->InsertOut(root);
	}
	else{
		ostringstream sum;
		sum << atof((from->GetBranch()->GetName()).c_str()) + atof((from->Next()->GetBranch()->GetName()).c_str());
		from->GetBranch()->SetName(sum.str());
		from->Next()->Out()->SetBranch(from->GetBranch());
		from->Out()->InsertOut(from->Next()->Out());
		delete from->Next()->GetBranch();
		delete from->GetNode();
		delete from->Next();
		delete from;
	}
}


void Tree::RootAt(Link* from)	{
	Link* link = from->Next();
	root->AppendTo(from);
	link->AppendTo(root);
	root->SetNode(from->GetNode());
}

void Tree::EraseInternalNodeName()	{
	EraseInternalNodeName(GetRoot());
}

void Tree::EraseInternalNodeName(Link * from)	{
	if (! from->isLeaf())	{
		from->GetNode()->SetName("");
	}
	for(Link* link=from->Next(); link!=from; link=link->Next())	{
		EraseInternalNodeName(link->Out());
	}
}

void Tree::RegisterWith(const TaxonSet* taxset)	{
	int tot = 0;
	if(!RegisterWith(taxset,GetRoot(),tot)){
		cout << "There is no match between the tree and the sequences.\n";
		exit(1);
	}
	if (tot != taxset->GetNtaxa())	{
		cerr << "error : non matching number of taxa : " << tot << '\t' << taxset->GetNtaxa() << '\n';
		cerr << "some taxa in the dataset are not present in the tree\n";
		exit(1);
	}
}

bool Tree::RegisterWith(const TaxonSet* taxset, Link* from, int& tot)	{
	if (from->isLeaf())	{
		int i = taxset->GetTaxonIndex(from->GetNode()->GetName());
		if (i != -1)	{
			from->GetNode()->SetIndex(i);
			tot++;
		}
		return(i != -1);
	}
	else{
		Link* previous = from;
		while(previous->Next() != from)	{
			if(RegisterWith(taxset,previous->Next()->Out(),tot)){
				previous = previous->Next();
			}
			else{
				// cout << "delete !!\n";
				DeleteNextLeaf(previous);
			}
		}
		if (! from)	{
			cerr << "from is null\n";
			exit(1);
		}
		if(from->isUnary()){
			// cerr << "DELETE UNARY NODE\n";
			DeleteUnaryNode(from);
		}
		return (!from->isLeaf());
	}
}


typedef list<string>::const_iterator csit;

Tree::Tree(string filename)	{

	ifstream is(filename.c_str());
	if (!is)	{
		cerr << "cannot find file : " << filename << '\n';
		exit(1);
	}
	ReadFromStream(is);
}

void Tree::ReadFromStream(istream& is)	{

	string expr = "";
	int cont = 1;
	while (cont)	{
		string s;
		is >> s;
		unsigned int k = 0;
		while ((k<s.length()) && (s[k] != ';'))	k++;
		expr += s.substr(0,k);
		cont = (! is.eof()) && (k==s.length());
	}
	SetRoot(ParseGroup(expr,0));
}

Link* Tree::ParseList(string input, Node* node)	{

	try	{

		// parse input as a list of strings separated by ','
		list<string> lst;
		int n = input.size();
		int k = 0;
		int brack = 0;
		int b = 0;
		while (k<n)	{
			char c = input[k];
			if (c == '(') brack++;
			if (c == ')') brack--;
			if ((!brack) && (c == ','))	{
				lst.push_back((string) (input.substr(b,k-b)));
				b = k+1;
			}

			if (brack<0)	{
				cout << "in parse list : too many )\n";
				cout << input.substr(0,k) << '\n';
				cout << input << '\n';
				throw;
			}
			k++;
		}
		if (brack)	{
			cout << "in parse list : too many (\n";
			cout << input << '\n';
			throw;
		}
		lst.push_back(input.substr(b,k-b));

		// make a circular single link chain around the node
		// with one link for each term of the list
		// and call parse group on each term
		Link* firstlink = new Link;
		Link* prevlink = firstlink;
		firstlink->SetNode(node);
		for (csit i=lst.begin(); i!=lst.end(); i++)	{
			Link* link = new Link;
			link->SetNode(node);
			link->AppendTo(prevlink);
			ParseGroup(*i,link);
			prevlink = link;
		}
		firstlink->AppendTo(prevlink);
		return firstlink;
	}
	catch(...)	{
		cout << "exit in parse list\n";
		exit(1);
	}
}

Link* Tree::ParseGroup(string input, Link* from)	{

	try	{

		// parse input as (body)nodeval:branchval

		string body = "";
		unsigned int k = 0;
		if (input[0] == '(')	{
			int brack = 1;
			k = 1;
			while ((k<input.length()) && brack)	{
				char c = input[k];
				if (c == '(') brack++;
				if (c == ')') brack--;
				k++;
			}
			if (brack)	{
				cout << "in parse group: too many (\n";
				cout << input << '\n';
				throw;
			}
			body = input.substr(1,k-2);
		}

		int b = k;
		while ((k<input.length()) && (input[k]!=':')) k++;
		string nodeval = input.substr(b,k-b);

		string branchval = "";
		if (k<input.length())	{
			branchval = input.substr(k+1, input.length() - k);
		}

		// make a new node and a new branch
		Node* node = new Node(nodeval);

		// call parse body
		Link* link = 0;
		if (body != "")	{
			link = ParseList(body,node);
		}
		else	{
			link = new Link;
			link->SetNode(node);
		}
		if (from)	{
			Branch* branch = new Branch(branchval);
			link->SetBranch(branch);
			from->SetBranch(branch);
			link->InsertOut(from);
		}
		return link;
	}
	catch(...)	{
		cout << "exit in parse group\n";
		exit(1);
	}
}

void Tree::Subdivide(Link* from, int Ninterpol)	{

	for (Link* link=from->Next(); link!=from; link=link->Next())	{
		Subdivide(link->Out(),Ninterpol);
	}
	// if ((! from->isLeaf()) && (! from->isRoot()))	{
	if (! from->isRoot())	{
		double l = atof(from->GetBranch()->GetName().c_str());
		if (l <= 0)	{
			cerr << "warning : non strictly positive branch length : " << l << '\n';
			cerr << "correcting and setting to 0.001\n";
			l = 0.001;
		}

		ostringstream s;
		s << l / Ninterpol;

		delete from->GetBranch();

		Link* current = from;
		Link* final = from->Out();
		int i = 0;
		while (i < Ninterpol-1)	{
			Link* link1 = new Link;
			Link* link2 = new Link;
			Branch* newbranch = new Branch(s.str());
			Node* newnode = new Node();
			current->SetBranch(newbranch);
			link1->SetNext(link2);
			link2->SetNext(link1);
			link1->SetBranch(newbranch);
			link1->SetNode(newnode);
			link2->SetNode(newnode);
			current->SetOut(link1);
			link1->SetOut(current);
			current = link2;
			i++;
		}
		current->SetOut(final);
		final->SetOut(current);
		Branch* newbranch = new Branch(s.str());
		final->SetBranch(newbranch);
		current->SetBranch(newbranch);
	}
}


int Tree::CountInternalNodes(const Link* from)	{
	int total = 0;
	if (! from->isLeaf())	{
	// if ((! from->isLeaf()) && (! from->isRoot()))	{
		total = 1;
		for (const Link* link=from->Next(); link!=from; link=link->Next())	{
			total += CountInternalNodes(link->Out());
		}
	}
	return total;
}

const Link* Tree::ChooseInternalNode(const Link* from, const Link*& fromup, int& n)	{
	if (from->isLeaf())	{
		return 0;
	}
	const Link* ret = 0;
	if (!n)	{
		ret = from;
	}
	else	{
		n--;
		for (const Link* link=from->Next(); link!=from; link=link->Next())	{
			if (!ret)	{
				const Link* tmp = ChooseInternalNode(link->Out(), fromup, n);
				if (tmp)	{
					ret = tmp;
				}
			}
		}
		if (ret && (! fromup))	{
			fromup = from;
		}
	}
	return ret;
}

int Tree::CountNodes(const Link* from)	{
	int total = 1;
	for (const Link* link=from->Next(); link!=from; link=link->Next())	{
		total += CountNodes(link->Out());
	}
	return total;
}

const Link* Tree::ChooseNode(const Link* from, const Link*& fromup, int& n)	{
	const Link* ret = 0;
	if (!n)	{
		ret = from;
	}
	else	{
		n--;
		for (const Link* link=from->Next(); (!ret && (link!=from)); link=link->Next())	{
			const Link* tmp = ChooseNode(link->Out(), fromup, n);
			if (tmp)	{
				ret = tmp;
			}
		}
		if (ret && (! fromup))	{
			fromup = from;
		}
	}
	return ret;
}
/*
void Tree::Print(ostream& os, const Link* from)	const {

	if (!from->isLeaf())	{
		os << '(';
		for (const Link* link=from->Next(); link!=from; link=link->Next())	{
			Print(os,link->Out());
			if (link->Next() != from)	{
				os << ',';
			}
		}
		os << ')';
	}
	os << GetNodeName(from);
	if (!from->isRoot())	{
		string brval = GetBranchName(from);
		if (brval != "")	{
			os << ':' << brval;
		}
	}
}

void Tree::Print(ostream& os)	const {
	Print(os,GetRoot());
	os << ";\n";
}
*/
