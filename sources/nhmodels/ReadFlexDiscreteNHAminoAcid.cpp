
#include "Sample.h"
#include "FlexDiscreteNHLGAminoAcidModel.h"
#include "MeanValTree.h"
#include "BranchMeanTree.h"

class MeanAllocTree  {

	public:

	MeanAllocTree(Tree* intree) {
		tree = intree;
		Reset();
	}

	Tree* GetTree() {return tree;}

	void Reset()	{
		RecursiveReset(GetTree()->GetRoot());
		size = 0;
	}

	void Normalise()	{
		RecursiveNormalise(GetTree()->GetRoot());
	}

	void Add(BinaryAllocationTree* sample, LengthTree* chronogram)	{
		RecursiveAdd(sample, chronogram, GetTree()->GetRoot());
		size++;
	}

	void ToStream(ostream& os)	{
		RecursiveSetName(os,GetTree()->GetRoot());
		tree->Print(os);
	}

	private:

	void RecursiveAdd(BinaryAllocationTree* sample, LengthTree* chronogram, Link* from)	{
		for(const Link* link=from->Next(); link!=from; link=link->Next())	{
			RecursiveAdd(sample, chronogram, link->Out());
			branchval[link->GetBranch()] += sample->GetBranchVal(link->GetBranch())->val();
			double time = chronogram->GetBranchVal(link->GetBranch())->val();
			meantime[link->GetBranch()] += time;
			vartime[link->GetBranch()] += time * time;
		}
	}

	void RecursiveReset(Link* from)	{
		for(const Link* link=from->Next(); link!=from; link=link->Next())	{
			RecursiveReset(link->Out());
			branchval[link->GetBranch()] = 0;
			meantime[link->GetBranch()] = 0 ;
			vartime[link->GetBranch()] = 0;
		}
	}

	void RecursiveNormalise(Link* from)	{
		for(const Link* link=from->Next(); link!=from; link=link->Next())	{
			RecursiveNormalise(link->Out());
			branchval[link->GetBranch()] /= size;
			meantime[link->GetBranch()] /= size;
			vartime[link->GetBranch()] /= size;
		}
	}

	void RecursiveSetName(ostream& os, Link* from)	{
		for(const Link* link=from->Next(); link!=from; link=link->Next())	{
			ostringstream s;
			if (link->Out()->isLeaf())	{
				string t = link->Out()->GetNode()->GetName();
				unsigned int i = 0;
				char c = ' ';
				while ((i<t.size()) && (c != '@'))	{
					c = t[i];
					if (c != '@')	{
						s << c;
					}
					i++;
				}
				s << '@';
			}
			s << Translate(branchval[link->GetBranch()]);
			link->Out()->GetNode()->SetName(s.str());
			RecursiveSetName(os, link->Out());
			ostringstream s2;
			s2 << meantime[link->GetBranch()];
			link->GetBranch()->SetName(s2.str());
		}
	}

	string Translate(double from)	{
		ostringstream s;
		s << ((int) (100 * from));
		return s.str();
	}

	Tree* tree;
	int size;


	map<const Branch*,double> meantime;
	map<const Branch*,double> vartime;
	map<const Branch*, double> branchval;

};

class FlexDiscreteNHAminoAcidSample : public Sample	{

	private:
	string modeltype;
	string datafile;
	string treefile;
	string contdatafile;
	bool withsep;
	GeneticCodeType type;

	public:

	string GetModelType() {return modeltype;}

	// FlexDiscreteNHAminoAcidModel* GetModel() {return (FlexDiscreteNHAminoAcidModel*) model;}
	FlexDiscreteNHAminoAcidModel* GetModel() {return (FlexDiscreteNHAminoAcidModel*) model;}

	FlexDiscreteNHAminoAcidSample(string filename, int inburnin, int inevery, int inuntil) : Sample(filename,inburnin,inevery,inuntil)	{
		Open();
	}

	void Open()	{

		// open <name>.param
		ifstream is((name + ".param").c_str());

		// check that file exists
		if (!is)	{
			cerr << "error : cannot find file : " << name << ".param\n";
			exit(1);
		}

		// read model type, and other standard fields
		is >> modeltype;
		is >> type;
		is >> datafile >> treefile >> contdatafile;
		is >> withsep;
		is >> chainevery >> chainuntil >> chainsize;
		cerr << modeltype << '\n';
		cerr << type << '\n';
		cerr << datafile << '\t' << treefile << '\t' << contdatafile << '\n';
		cerr << withsep << '\n';
		cerr << every << '\t' << until << '\t' << size << '\n';

		// the chain's saving frequency, upper limit and current size
		// not to be confused with the sample's subsampling frequency, upper limit and size

		// make a new model depending on the type obtained from the file
		if (modeltype == "FLEXDISCRETENHAMINOACID")	{
			model = new FlexDiscreteNHAminoAcidModel(datafile,treefile,contdatafile,withsep,true,type);
		}
		else	{
			cerr << "error when opening file "  << name << " : does not recognise model type : " << modeltype << '\n';
			exit(1);
		}

		model->FromStream(is);
		// model->Update();

		OpenChainFile();
	}

	void Read()	{

		int Ncont = GetModel()->Ncont;

		MeanAllocTree** meanalloctree = new MeanAllocTree*[Ncont];
		for (int k=0; k<Ncont; k++)	{
			meanalloctree[k] = new MeanAllocTree(GetModel()->GetTree());
			meanalloctree[k]->Reset();
		}

		cerr << "size : " << size << '\n';

		for (int i=0; i<size; i++)	{
			cerr << '.';

			GetNextPoint();

			for (int k=0; k<Ncont; k++)	{
				meanalloctree[k]->Add(GetModel()->GetAllocationTree(k), GetModel()->GetGamTree());
			}
		}
		cerr << '\n';

		for (int k=0; k<Ncont; k++)	{
			meanalloctree[k]->Normalise();
			ostringstream s;
			s << GetName() << ".postmeanalloc" << k << ".tre";
			ofstream os(s.str().c_str());
			meanalloctree[k]->ToStream(os);
		}

		cerr << "done\n";
	}
};


int main(int argc, char* argv[])	{

	int burnin = 0;
	int every = 1;
	int until = -1;
	string name;

	int ppred = 0;
	int site = 0;

	try	{

		if (argc == 1)	{
			throw(0);
		}

		int i = 1;
		while (i < argc)	{
			string s = argv[i];

			if (s == "-pp")	{
				ppred = 1;
				i++;
				site = atoi(argv[i]);
			}
			else if ( (s == "-x") || (s == "-extract") )	{
				i++;
				if (i == argc) throw(0);
				s = argv[i];
				burnin = atoi(argv[i]);
				i++;
				if (i == argc) throw(0);
				s = argv[i];
				every = atoi(argv[i]);
				i++;
				if (i == argc) throw(0);
				s = argv[i];
				until = atoi(argv[i]);
			}
			else	{
					if (i != (argc -1))	{
					throw(0);
				}
				name = argv[i];
			}
			i++;
		}
		if (name == "")	{
			throw(0);
		}
	}
	catch(...)	{
		cerr << "readpb [-x <burnin> <every> <until>] <chainname> \n";
		cerr << '\n';
		exit(1);
	}

	FlexDiscreteNHAminoAcidSample sample(name,burnin,every,until);
	sample.Read();

}

