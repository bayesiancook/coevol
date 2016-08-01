
#include "Sample.h"
#include "MixBGCModel.h"
#include "MeanValTree.h"
#include "MeanICTree.h"
#include "MeanCovMatrix.h"
#include "StringStreamUtils.h"
#include "MeanChronogram.h"
#include "MeanChronoBubbleTree.h"

#include "TexTab.h"


class MeanVarTree : public NewickTree {

	public:

	MeanVarTree(Tree* intree, bool inwithroot = false) : tree(intree), withRoot(inwithroot)	{
		Reset();
	}

	const Tree* GetTree() const {return tree;}

	const Link* GetRoot() const {return GetTree()->GetRoot();}

	void Reset()	{
		RecursiveReset(GetTree()->GetRoot());
		size = 0;
	}

	void Add(BGCModel* sample)	{
		RecursiveAdd(sample,GetRoot());
		size++;
	}

	void ResetCorrel()	{
		Rec1ResetCorrel(GetRoot());
	}

	void AddCorrel(BGCModel* sample)	{
		RecursiveAdd(sample,GetRoot());
		Rec1AddCorrel(sample,GetRoot());
		size++;
	}

	void NormaliseCorrel()	{
		RecursiveNormalise(GetRoot());
		Rec1NormaliseCorrel(GetRoot());
	}

	void TabulateCorrel(ostream& os)	{
		Rec1TabulateCorrel(os,GetRoot());
	}

	void Normalise()	{
		RecursiveNormalise(GetTree()->GetRoot());
	}

	void Tabulate(ostream& os)	{
		RecursiveTabulate(os,GetTree()->GetRoot());
	}

	void TabulatePart2(ostream& os)	{
		RecursiveTabulatePart2(os,GetTree()->GetRoot());
	}

	void TabulatePart(ostream& os)	{
		RecursiveTabulatePart(os,GetTree()->GetRoot());
	}

	double GetMeanTime(const Branch* branch) const	{
		map<const Branch*, double>::const_iterator i = meantime.find(branch);
		return i->second;
	}

	double GetMeanBGC(const Branch* branch) const	{
		map<const Branch*, double>::const_iterator i = meanbgc.find(branch);
		return i->second;
	}

	double GetVarBGC(const Branch* branch) const	{
		map<const Branch*, double>::const_iterator i = varbgc.find(branch);
		return i->second;
	}

	string GetNodeName(const Link* link) const {
		return link->GetNode()->GetName();
	}

	string GetBranchName(const Link* link) const {
		ostringstream s;
		s << GetVarBGC(link->GetBranch());
		return s.str();
	}

	protected:

	void Rec1AddCorrel(BGCModel* sample, const Link* from)	{
		Rec2AddCorrel(sample,from,GetRoot());
		for(const Link* link=from->Next(); link!=from; link=link->Next())	{
			Rec1AddCorrel(sample, link->Out());
		}
	}

	void Rec2AddCorrel(BGCModel* sample, const Link* from1, const Link* from2)	{
		if ((!from1->isRoot()) && (!from2->isRoot()))	{
			double time = 0;
			double r = sample->GetBranchCorrelation(from1,from2,time);
			pair<const Branch*, const Branch*> p(from1->GetBranch(),from2->GetBranch());
			covarbgc[p] += r;
			deltatime[p] += time;
		}
		for (const Link* link = from2->Next(); link!=from2; link=link->Next())	{
			Rec2AddCorrel(sample,from1,link->Out());
		}

	}

	void Rec1NormaliseCorrel(const Link* from)	{
		Rec2NormaliseCorrel(from,GetRoot());
		for(const Link* link=from->Next(); link!=from; link=link->Next())	{
			Rec1NormaliseCorrel(link->Out());
		}
	}

	void Rec2NormaliseCorrel(const Link* from1, const Link* from2)	{
		if ((!from1->isRoot()) && (!from2->isRoot()))	{
			pair<const Branch*, const Branch*> p(from1->GetBranch(),from2->GetBranch());
			covarbgc[p] /= size;
			deltatime[p] /= size;
		}
		for (const Link* link = from2->Next(); link!=from2; link=link->Next())	{
			Rec2NormaliseCorrel(from1,link->Out());
		}

	}

	void Rec1ResetCorrel(const Link* from)	{
		Rec2ResetCorrel(from,GetRoot());
		for(const Link* link=from->Next(); link!=from; link=link->Next())	{
			Rec1ResetCorrel(link->Out());
		}
	}

	void Rec2ResetCorrel(const Link* from1, const Link* from2)	{
		if ((!from1->isRoot()) && (!from2->isRoot()))	{
		// if ((from1 != from2) && (!from1->isRoot()) && (!from2->isRoot()))	{
			pair<const Branch*, const Branch*> p(from1->GetBranch(),from2->GetBranch());
			covarbgc[p] = 0;
			deltatime[p] = 0;
		}
		for (const Link* link = from2->Next(); link!=from2; link=link->Next())	{
			Rec2ResetCorrel(from1,link->Out());
		}
	}

	void Rec1TabulateCorrel(ostream& os, const Link* from)	{
		Rec2TabulateCorrel(os,from,GetRoot());
		for(const Link* link=from->Next(); link!=from; link=link->Next())	{
			Rec1TabulateCorrel(os, link->Out());
		}
	}

	void Rec2TabulateCorrel(ostream& os, const Link* from1, const Link* from2)	{
		if ((!from1->isRoot()) && (!from2->isRoot()))	{
		// if ((from1 != from2) && (!from1->isRoot()) && (!from2->isRoot()))	{
			cerr << meantime[from1->GetBranch()] << '\t' << meantime[from2->GetBranch()] << '\n';
			if ((meantime[from1->GetBranch()] > 0.3) || (meantime[from2->GetBranch()] > 0.3))	{
				pair<const Branch*, const Branch*> p(from1->GetBranch(),from2->GetBranch());
				if (deltatime[p] == 0)	{
					cerr << GetTree()->GetLeftMost(from1) << '\t' << GetTree()->GetRightMost(from1) << '\t';
					cerr << '\t' << GetTree()->GetLeftMost(from2) << '\t' << GetTree()->GetRightMost(from2) << '\n';
				}
				os << deltatime[p] << '\t' << covarbgc[p] << '\n';
			}
		}
		for (const Link* link = from2->Next(); link!=from2; link=link->Next())	{
			Rec2TabulateCorrel(os,from1,link->Out());
		}
	}

	void RecursiveAdd(BGCModel* sample, const Link* from)	{
		if (withRoot || (!from->isRoot()))	{
			double time = sample->GetLengthTree()->GetBranchVal(from->GetBranch())->val();
			meantime[from->GetBranch()] += time;
			double bgc = sample->GetBGCProcess()->GetBranchVal(from->GetBranch())->val();
			meanbgc[from->GetBranch()] += bgc;
			double vbgc = sample->GetBranchVarCoeff(from);
			varbgc[from->GetBranch()] += vbgc;
			meangreater1[from->GetBranch()] += sample->GetPropGreaterThan(from,1.0);
			meangreater5[from->GetBranch()] += sample->GetPropGreaterThan(from,5.0);
			meangreater10[from->GetBranch()] += sample->GetPropGreaterThan(from,10.0);
		}
		for(const Link* link=from->Next(); link!=from; link=link->Next())	{
			RecursiveAdd(sample, link->Out());
		}
	}

	void RecursiveNormalise(const Link* from)	{
		if (withRoot || (!from->isRoot()))	{
			meantime[from->GetBranch()] /= size;
			meangreater1[from->GetBranch()] /= size;
			meangreater5[from->GetBranch()] /= size;
			meangreater10[from->GetBranch()] /= size;
			meanbgc[from->GetBranch()] /= size;
			varbgc[from->GetBranch()] /= size;
		}
		for(const Link* link=from->Next(); link!=from; link=link->Next())	{
			RecursiveNormalise(link->Out());
		}
	}

	void RecursiveReset(const Link* from)	{
		if (withRoot || (!from->isRoot()))	{
			meantime[from->GetBranch()] = 0;
			meangreater1[from->GetBranch()] = 0;
			meangreater5[from->GetBranch()] = 0;
			meangreater10[from->GetBranch()] = 0;
			meanbgc[from->GetBranch()] = 0;
			varbgc[from->GetBranch()] = 0;
		}
		for(const Link* link=from->Next(); link!=from; link=link->Next())	{
			RecursiveReset(link->Out());
		}
	}

	void RecursiveTabulatePart(ostream& os, const Link* from)	{
		if (withRoot || (!from->isRoot()))	{
			os << GetTree()->GetLeftMost(from) << '\t' << GetTree()->GetRightMost(from) << '\t';
			os << '\t' << meantime[from->GetBranch()];
			os << '\t' << meanbgc[from->GetBranch()];
			os << '\t' << varbgc[from->GetBranch()];
			os << '\t' << sqrt(varbgc[from->GetBranch()]);
			os << '\t' << meangreater1[from->GetBranch()];
			os << '\t' << meangreater5[from->GetBranch()];
			os << '\t' << meangreater10[from->GetBranch()];
			os << '\n';
		}
		for(const Link* link=from->Next(); link!=from; link=link->Next())	{
			RecursiveTabulatePart(os,link->Out());
		}
	}

	void RecursiveTabulatePart2(ostream& os, const Link* from)	{
		if (withRoot || (!from->isRoot()))	{
			os << GetTree()->GetLeftMost(from) << '\t' << GetTree()->GetRightMost(from) << '\t';
			int tmp = 0;
			tmp = (int) (100 * meangreater1[from->GetBranch()]);
			if (tmp)	{
				os << tmp;
			}
			else	{
				os << '-';
			}
			os << '/';

			/*
			tmp = (int) (100 * meangreater5[from->GetBranch()]);
			if (tmp)	{
				os << tmp;
			}
			else	{
				os << '-';
			}
			os << '/';
			*/

			tmp = (int) (100 * meangreater10[from->GetBranch()]);
			if (tmp)	{
				os << tmp;
			}
			else	{
				os << '-';
			}
			os << '\n';
		}
		for(const Link* link=from->Next(); link!=from; link=link->Next())	{
			RecursiveTabulatePart2(os,link->Out());
		}
	}

	void RecursiveTabulate(ostream& os, const Link* from)	{
		if (withRoot || (!from->isRoot()))	{
			os << GetTree()->GetLeftMost(from) << '\t' << GetTree()->GetRightMost(from);
			os << '\t' << meantime[from->GetBranch()];
			os << '\t' << meanbgc[from->GetBranch()];
			os << '\t' << varbgc[from->GetBranch()];
			os << '\t' << sqrt(varbgc[from->GetBranch()]);
			/*
			int tmp = 0;
			tmp = (int) (100 * meangreater1[from->GetBranch()]);
			if (tmp)	{
				os << tmp;
			}
			else	{
				os << '-';
			}
			os << '/';

			tmp = (int) (100 * meangreater5[from->GetBranch()]);
			if (tmp)	{
				os << tmp;
			}
			else	{
				os << '-';
			}
			os << '/';

			tmp = (int) (100 * meangreater10[from->GetBranch()]);
			if (tmp)	{
				os << tmp;
			}
			else	{
				os << '-';
			}
			*/
			/*
			os << '\t' << meangreater1[from->GetBranch()];
			os << '\t' << meangreater5[from->GetBranch()];
			os << '\t' << meangreater10[from->GetBranch()];
			*/
			os << '\n';
		}
		for(const Link* link=from->Next(); link!=from; link=link->Next())	{
			RecursiveTabulate(os,link->Out());
		}
	}


	map<const Branch*,double> meantime;
	map<const Branch*,double> meangreater1;
	map<const Branch*,double> meangreater5;
	map<const Branch*,double> meangreater10;
	map<const Branch*,double> meanbgc;
	map<const Branch*,double> varbgc;

	map<pair<const Branch*, const Branch*>, double> covarbgc;
	map<pair<const Branch*, const Branch*>, double> deltatime;

	Tree* tree;
	bool withRoot;
	int size;

};


class BGCSample : public Sample	{

	private:
	string modeltype;
	string datafile;
	string treefile;
	string contdatafile;
	string calibfile;
	string suffstatfile;
	string rootfile;

	double rootage;
	double rootstdev;

	int chronoprior;
	double meanchi;
	double meanchi2;

	bool clampdiag;
	bool autoregressive;
	bool clamptree;
	bool meanexp;

	double priorsigma;
	double fixalpha;

	int conjpath;
	int contdatatype;

	bool normalise;
	int nrep;

	int df;

	int clampbgcoffset;
	int flexrho;

	int nsplit;

	double lambda;
	double at2cg;
	double at2gc;
	double gc2cg;

	double cg2at;
	double cg2gc;
	double cg2ta;
	int discgam;
	int discn;
	int triplet;
	int clampreg;

	public:

	string GetModelType() {return modeltype;}

	BGCModel* GetModel() {return (BGCModel*) model;}

	BGCSample(string filename, int inburnin, int inevery, int inuntil) : Sample(filename,inburnin,inevery,inuntil)	{
		suffstatfile = "None";
		rootfile = "None";
		autoregressive = false;
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

		priorsigma = 1;
		fixalpha = 0;

		nsplit = 1;

		clampbgcoffset = 0;
		flexrho = 0;

		// read model type, and other standard fields
		is >> modeltype;
		is >> datafile >> treefile >> contdatafile;
		is >> calibfile >> rootage >> rootstdev;
		is >> chronoprior >> meanchi >> meanchi2;
		is >> clampdiag;
		is >> conjpath;
		is >> contdatatype;
		is >> meanexp;
		is >> normalise >> nrep;
		is >> df;
		is >> priorsigma;
		is >> nsplit;
		is >> clamptree;
		is >> suffstatfile;
		is >> autoregressive;
		is >> rootfile;

		int check;
		is >> check;
		if (check)	{
			is >> fixalpha;
			is >> check;
			if (check)	{
				is >> clampbgcoffset;
				is >> flexrho;
				is >> check;
				if (check)	{
					is >> lambda;
					is >> at2cg;
					is >> at2gc;
					is >> gc2cg;
					is >> check;
					if (check)	{
						is >> discgam;
						is >> cg2at;
						is >> cg2gc;
						is >> cg2ta;
						is >> discn;
						is >> check;
						if (check)	{
							is >> triplet;
							is >> check;
							if (check)	{
								is >> clampreg;
								is >> check;
								if (check)	{
									cerr << "error when reading model\n";
									exit(1);
								}
							}
						}
					}
				}
			}
		}

		is >> chainevery >> chainuntil >> chainsize;
		// the chain's saving frequency, upper limit and current size
		// not to be confused with the sample's subsampling frequency, upper limit and size

		// make a new model depending on the type obtained from the file
		if (modeltype == "BGC")	{
			model = new BGCModel(datafile,treefile,contdatafile,calibfile,rootage,rootstdev,chronoprior,meanchi,meanchi2,priorsigma,df,clampdiag,autoregressive,conjpath,contdatatype,clamptree,meanexp,normalise,nrep,nsplit,suffstatfile,rootfile,fixalpha,clampbgcoffset,flexrho,lambda,at2cg,at2gc,gc2cg,cg2at,cg2gc,cg2ta,discgam,discn,triplet,clampreg,false);
		}
		else	{
			cerr << "error when opening file "  << name << " : does not recognise model type : " << modeltype << '\n';
			exit(1);
		}

		// read model (i.e. chain's last point) from <name>.param
		model->FromStream(is);
		// model->Update();

		// open <name>.chain, and prepare stream and stream iterator
		OpenChainFile();
		// now, size is defined (it is the total number of points with which this Sample object will make all its various posterior averages)
		// all these points can be accessed to (only once) by repeated calls to GetNextPoint()

		cerr << "number of points to read : " << size << '\n';
		cerr << '\n';
	}

	void ReadVariance(string name)	{

		list<double> varcoeff;
		double meanvarcoeff = 0;
		double varvarcoeff = 0;

		MeanVarTree* vartree = new MeanVarTree(GetModel()->GetFineGrainedTree(),false);

		// cycle over the sample
		for (int i=0; i<size; i++)	{
			cerr << '.';

			// get next point -> will be stored into "model", and thus, will be accessible through GetModel()

			GetNextPoint();
			GetModel()->UpdateBGC();
			double tmp = GetModel()->GetBranchVarCoeff();
			varcoeff.push_front(tmp);
			meanvarcoeff += tmp;
			varvarcoeff += tmp * tmp;
			vartree->Add(GetModel());
		}
		cerr << '\n';
		cerr << "normalise\n";

		vartree->Normalise();
		ofstream os((GetName() + ".var").c_str());
		vartree->Tabulate(os);

		ofstream os2((GetName() + ".branchval").c_str());
		vartree->TabulatePart(os2);

		ofstream os3((GetName() + ".branchval2").c_str());
		vartree->TabulatePart2(os3);
		os3 << "END END END\n";

		ofstream bos((GetName() + ".var.tre").c_str());
		vartree->ToStream(bos);

		meanvarcoeff /= size;
		varvarcoeff /= size;
		varvarcoeff -= meanvarcoeff * meanvarcoeff;

		/*
		double c = 0.975;
		varcoeff.sort();
		int n = ((int) (((double) varcoeff.size()) * (1-c)));
		list<double>::const_iterator i = varcoeff.begin();
		for (int j=0; j<n; j++)	{
			i++;
		}
		double min = *i;
		n = ((int) (((double) varcoeff.size()) * c));
		i = varcoeff.begin();
		for (int j=0; j<n; j++)	{
			i++;
		}
		double max = *i;
		*/

		ofstream vos((GetName() + ".varcoeff").c_str());
		// vos << name << '\t' << '&' << '\t' << meanvarcoeff << "\\\\\n";
		// vos << meanvarcoeff << '\t' << min << '\t' << max << '\n';
		// vos << meanvarcoeff << '\t' << sqrt(meanvarcoeff) << '\n';
		vos << meanvarcoeff << '\n';

		cerr << '\n';
		cerr << "tabulated variances in " << GetName() << ".var\n";
		cerr << '\n';
	}

	void ReadBranchCorrel(string name)	{

		MeanVarTree* vartree = new MeanVarTree(GetModel()->GetFineGrainedTree(),false);
		vartree->ResetCorrel();

		for (int i=0; i<size; i++)	{
			cerr << '.';

			// get next point -> will be stored into "model", and thus, will be accessible through GetModel()

			GetNextPoint();
			GetModel()->UpdateBGC();
			GetModel()->GetCalibratedChronogram()->specialUpdate();
			vartree->AddCorrel(GetModel());
		}
		cerr << '\n';
		cerr << "normalise\n";

		vartree->NormaliseCorrel();
		ofstream os((GetName() + ".correl").c_str());
		vartree->TabulateCorrel(os);

		cerr << '\n';
		cerr << "tabulated correlations in " << GetName() << ".correl\n";
		cerr << '\n';
	}

	void ReadSuccessiveBranchesCorrel()	{

		double mean = 0;
		double above10 = 0;
		double above5 = 0;
		double joint10 = 0;
		double joint5 = 0;
		double two10 = 0;
		double two5 = 0;
		double one10 = 0;
		double one5 = 0;
		double zero5 = 0;

		for (int i=0; i<size; i++)	{
			cerr << '.';

			// get next point -> will be stored into "model", and thus, will be accessible through GetModel()

			GetNextPoint();
			GetModel()->UpdateBGC();
			mean += GetModel()->GetSuccessiveBranchesCorrelation();
			double j10 = 0;
			double j5 = 0;
			double t10 =0;
			double t5 =0;
			double o10 = 0;
			double o5 = 0;
			double z10 = 0;
			double z5 = 0;
			above10 += GetModel()->GetMeanFractionAbove(10.0,j10,t10,o10,z10);
			above5 += GetModel()->GetMeanFractionAbove(5.0,j5,t5,o5,z5);
			joint10 += j10;
			joint5 += j5;
			two5 += t5;
			one5 += o5;
			zero5 += z5;
		}
		cerr << '\n';

		mean /= size;
		above10 /= size;
		above5 /= size;
		joint10 /= size;
		joint5 /= size;
		two10 /= size;
		two5 /= size;
		one10 /= size;
		one5 /= size;
		zero5 /= size;
		cerr << "mean correlation coeff between successive branches : " << mean << '\n';
		cerr << "mean number of branches with gene-spec bgc > 10x mean for that gene : " << above10 << '\n';
		cerr << "mean number of branches with gene-spec bgc > 5x mean for that gene : " << above5 << '\n';
		cerr << "fraction of genes with 0 such branch   : " << zero5 << '\n';
		cerr << "fraction of genes with 1 such branch   : " << one5 << '\n';
		cerr << "fraction of genes with 2 such branches : " << two5 << '\n';
		cerr << "fraction among them for which these 2 branches are successive : " << joint5  / two5 << '\n';

		cerr << '\n';
	}

	void ReadTexTable1(string modelname)	{

		int Ncont = GetModel()->GetContMatrix()->GetDim();
		list<double>* reg = new list<double>[Ncont];
		list<double> lambda;
		list<double> branchvar;
		list<double> alpha;

		double* correl = new double[Ncont];
		for (int i=0; i<Ncont; i++)	{
			correl[i] = 0;
		}

		double meanalpha = 0;

		ofstream fullos((GetName() + ".slopedist").c_str());

		for (int i=0; i<size; i++)	{
			cerr << '.';

			GetNextPoint();

			lambda.push_front(GetModel()->GetLambda());

			GetModel()->GetBGCProcess()->specialUpdate();
			GetModel()->UpdateGeneProcesses();
			branchvar.push_front(GetModel()->GetGeneVarCoeff());

			CovMatrix& m = *(GetModel()->GetContMatrix());

			for (int k=0; k<Ncont; k++)	{
				double tmp = GetModel()->GetRegCoef(k);
				fullos << tmp << '\t';
				reg[k].push_front(tmp);

				double cov = 0;
				for (int j=0; j<Ncont; j++)	{
					cov += GetModel()->GetRegCoef(j) * m[k][j];
				}
				double varb = GetModel()->GetBGCVar();
				for (int j1=0; j1<Ncont; j1++)	{
					for (int j2=0; j2<Ncont; j2++)	{
						varb += GetModel()->GetRegCoef(j1) * GetModel()->GetRegCoef(j2) * m[j1][j2];
					}
				}
				double cor = cov / sqrt(varb * m[k][k]);
				correl[k] += cor;
			}
			fullos << '\n';

			alpha.push_front(GetModel()->GetAlpha());
			meanalpha += GetModel()->GetAlpha();

		}
		cerr << '\n';

		for (int k=0; k<Ncont; k++)	{
			correl[k] /= size;
		}

		meanalpha /= size;

		ofstream cout((GetName() + "_table1.tex").c_str());

		if (modelname != "")	{
			cout << modelname;
		}
		for (int k=0; k<Ncont; k++)	{
			cout << "&\t";
			cout << textabentry(reg[k],true,true,true,2);
			cout << "&\t" << ((double) ((int) (100 * correl[k]))) / 100;
		}
		/*
		cout << textabentry(lambda,true,false,true,2);
		cout <<  "&\t";
		if (meanalpha)	{
			cout << textabentry(alpha,true,false,true,2);
		}
		else	{
			cout << "\\multicolumn{6}{c}{-}";
		}
		*/
		/*
		cout << textabentry(branchvar,true,false,false,2);
		*/
		cout << "\\\\\n";
		cout.close();

		ofstream os((GetName() + "_lambda.tex").c_str());
		os << textabentry(lambda,true,false,true,2);
		os.close();
	}

	void ReadNuc(string modelname)	{
		list<double> meanat2gc;
		list<double> meangc2at;
		list<double> meangc2ta;
		list<double> meanat2cg;
		list<double> meangc2cg;

		for (int i=0; i<size; i++)	{
			cerr << '.';

			GetNextPoint();

			if (!GetModel()->Reversible())	{
				meanat2gc.push_front(GetModel()->at2gc->val());
				meangc2at.push_front(GetModel()->cg2ta->val());
				meangc2ta.push_front(GetModel()->cg2at->val());
				meanat2cg.push_front(GetModel()->at2cg->val());
				meangc2cg.push_front(GetModel()->cg2gc->val());
			}
			else	{
				double tmp = GetModel()->at2ta->val() * GetModel()->lambda->val() / (1.0 + GetModel()->lambda->val());
				meanat2gc.push_front(GetModel()->at2gc->val() * 1.0 / (1.0 + GetModel()->lambda->val()) / tmp);
				meangc2at.push_front(GetModel()->at2gc->val() * GetModel()->lambda->val() / (1.0 + GetModel()->lambda->val()) / tmp);
				meangc2ta.push_front(GetModel()->at2cg->val() * GetModel()->lambda->val() / (1.0 + GetModel()->lambda->val()) / tmp);
				meanat2cg.push_front(GetModel()->at2cg->val() * 1.0 / (1.0 + GetModel()->lambda->val()) / tmp);
				meangc2cg.push_front(GetModel()->gc2cg->val() * 1.0 / (1.0 + GetModel()->lambda->val()) / tmp);
			}
		}

		cerr << '\n';

		ofstream cout((GetName() + "_table2.tex").c_str());
		cout << modelname << "&\t";
		cout << textabentry(meanat2gc,true,false,true,2);
		cout <<  "&\t";
		cout << textabentry(meangc2at,true,false,true,2);
		cout <<  "&\t";
		cout << textabentry(meangc2ta,true,false,true,2);
		cout <<  "&\t";
		cout << textabentry(meanat2cg,true,false,true,2);
		cout <<  "&\t";
		cout << textabentry(meangc2cg,true,false,true,2);
		cout << "\\\\\n";
		cout.close();

	}

	/*
	void ReadNuc()	{
		double meanat2gc = 0;
		double meangc2at = 0;
		double meanat2ta = 0;
		double meangc2ta = 0;
		double meanat2cg = 0;
		double meangc2cg = 0;

		for (int i=0; i<size; i++)	{
			cerr << '.';

			GetNextPoint();

			double tmp = GetModel()->at2ta->val() * GetModel()->lambda->val() / (1.0 + GetModel()->lambda->val());

			meanat2gc += GetModel()->at2gc->val() * 1.0 / (1.0 + GetModel()->lambda->val()) / tmp;
			meangc2at += GetModel()->at2gc->val() * GetModel()->lambda->val() / (1.0 + GetModel()->lambda->val()) / tmp;
			meanat2ta += GetModel()->at2ta->val() * GetModel()->lambda->val() / (1.0 + GetModel()->lambda->val()) / tmp;
			meangc2ta += GetModel()->at2cg->val() * GetModel()->lambda->val() / (1.0 + GetModel()->lambda->val()) / tmp;
			meanat2cg += GetModel()->at2cg->val() * 1.0 / (1.0 + GetModel()->lambda->val()) / tmp;
			meangc2cg += GetModel()->gc2cg->val() * 1.0 / (1.0 + GetModel()->lambda->val()) / tmp;
		}

		cerr << '\n';

		meanat2gc /= size;
		meangc2at /= size;
		meanat2ta /= size;
		meangc2ta /= size;
		meanat2cg /= size;
		meangc2cg /= size;

		cout << meanat2gc << '\t' << 4.29 / 1.29 << '\n';
		cout << meangc2at << '\t' << 9.61 / 1.29 << '\n';
		cout << meanat2ta << '\t' << 1 << '\n';
		cout << meangc2ta << '\t' << 2.58 / 1.29 << '\n';
		cout << meanat2cg << '\t' << 1.52 / 1.29 << '\n';
		cout << meangc2cg << '\t' << 2.95 / 1.29 << '\n';

	}
	*/

	void ReadIC()	{

		MeanICTree* ictree = new MeanICTree(GetModel()->GetContProcess());
		// MeanICTree* ppictree = new MeanICTree(GetModel()->GetContProcess());
		ictree->Reset();
		// ppictree->Reset();

		MeanICTree* bictree = new MeanICTree(GetModel()->GetLogBGCProcess(),GetModel()->GetLengthTree());
		// MeanICTree* ppbictree = new MeanICTree(GetModel()->GetLogBGCProcess(),GetModel()->GetLengthTree());
		bictree->Reset();
		// ppbictree->Reset();

		/*
		int dim = GetModel()->GetNcont();
		double* pp = new double[dim];
		for (int k=0; k<dim; k++)	{
			pp[k] = 0;
		}

		int bdim = GetModel()->GetNreg();
		double* bpp = new double[bdim];
		for (int k=0; k<bdim; k++)	{
			bpp[k] = 0;
		}
		*/

		for (int i=0; i<size; i++)	{
			cerr << '.';
			GetNextPoint();
			ictree->Add();
			bictree->Add();
			// GetModel()->Update();
			/*
			GetModel()->SampleProcess();
			ppictree->Add();
			ppbictree->Add();
			for (int k=0; k<dim; k++)	{
				if (ppictree->GetGmean(k) > ictree->GetGmean(k))	{
					pp[k]++;
				}
			}
			for (int k=0; k<bdim; k++)	{
				if (ppbictree->GetGmean(k) > bictree->GetGmean(k))	{
					bpp[k]++;
				}
			}
			*/
		}
		cerr << '\n';
		ictree->Normalise();
		// ppictree->Normalise();
		ofstream os((name + ".postmeanic").c_str());
		ictree->Tabulate(os,false);
		ofstream los((name + ".postmeanleafic").c_str());
		ictree->Tabulate(los,true);
		cerr << "mean independent contrasts in " << name  << ".postmeanic\n";
		cerr << "for terminal branches only in " << name  << ".postmeanleafic\n";
		/*
		cerr << "pp positive trend : \n";
		for (int k=0; k<dim; k++)	{
			cerr << k << '\t' << pp[k] / size << '\n';
		}
		cerr << '\n';
		*/
		ictree->DAgostinosKandZ();

		bictree->Normalise();
		// ppbictree->Normalise();
		ofstream bos((name + ".postmeanbic").c_str());
		bictree->Tabulate(bos,false);
		ofstream lbos((name + ".postmeanleafbic").c_str());
		bictree->Tabulate(lbos,true);
		cerr << "mean independent contrasts for BGC in " << name  << ".postmeanbic\n";
		cerr << "for terminal branches only in " << name  << ".postmeanleafbic\n";
		/*
		cerr << "pp positive trend : \n";
		for (int k=0; k<bdim; k++)	{
			cerr << k << '\t' << bpp[k] / size << '\n';
		}
		cerr << '\n';
		*/
		bictree->DAgostinosKandZ();
	}

	void ReadGeneBGC(bool printlog, bool printmean, bool printci, bool printstdev, bool withleaf, bool withinternal)	{

		int Ngene = GetModel()->GetNgene();
		MeanExpNormTree** genetree = new MeanExpNormTree*[Ngene];
		for (int gene=0; gene<Ngene; gene++)	{
			genetree[gene] = new MeanExpNormTree(GetModel()->GetFineGrainedTree(),false,printlog,printmean,printci,printstdev,withleaf,withinternal);
		}

		for (int i=0; i<size; i++)	{
			cerr << '.';

			// get next point -> will be stored into "model", and thus, will be accessible through GetModel()

			GetNextPoint();

			GetModel()->UpdateBGC();

			for (int gene=0; gene<Ngene; gene++)	{
				GetModel()->UpdateGeneProcess(gene);
				genetree[gene]->Add((BranchVarTree<PosReal>*) GetModel()->GetGeneProcess(gene),GetModel()->GetLengthTree());
			}
		}
		for (int gene=0; gene<Ngene; gene++)	{
			genetree[gene]->Normalise();
			ostringstream s;
			s << GetName() << ".genebgc" << gene << ".tre";
			ofstream os(s.str().c_str());
			genetree[gene]->ToStream(os);
		}
	}

	void Read(bool printlog, bool printmean, bool printci, bool printstdev, bool withleaf, bool withinternal, bool tex, double xscale, double yscale, double nodescale, double nodepower, double barwidth, int fontsize, bool bubbletext, bool withheader, bool postdist)	{

		// cerr << "df : " << GetModel()->df << '\n';
		int Ncont = GetModel()->GetContMatrix()->GetDim();
		MeanCovMatrix*  mat = new MeanCovMatrix(Ncont);
		MeanExpNormTree* meanbgc = new MeanExpNormTree(GetModel()->GetFineGrainedTree(),false,printlog,printmean,printci,printstdev,withleaf,withinternal,0,0);
		MeanExpNormTree* meansynrate = new MeanExpNormTree(GetModel()->GetFineGrainedTree(),false,printlog,printmean,printci,printstdev,withleaf,withinternal,0,0);
		MeanExpNormTree** tree = new MeanExpNormTree*[Ncont];
		for (int k=0; k<Ncont; k++)	{
			tree[k] = new MeanExpNormTree(GetModel()->GetFineGrainedTree(),false,printlog,printmean,printci,printstdev,withleaf,withinternal);
		}

		MeanChronogram* meanchrono = new MeanChronogram(GetModel()->GetTree());

		double* correl = new double[Ncont];
		double* meanreg = new double[Ncont];
		double* varreg = new double[Ncont];
		double* ppreg = new double[Ncont];
		for (int i=0; i<Ncont; i++)	{
			correl[i] = 0;
			meanreg[i] = 0;
			varreg[i] = 0;
			ppreg[i] = 0;
		}

		MeanExpNormTree* meanrho = 0;
		MeanBranchTree* meanbranchrho = 0;
		if (GetModel()->FlexRho())	{
			meanrho = new MeanExpNormTree(GetModel()->GetFineGrainedTree(),false,printlog,printmean,printci,printstdev,withleaf,withinternal,0,0);
			meanbranchrho = new MeanBranchTree(GetModel()->GetFineGrainedTree(),false);
		}

		double mnrho = 0;
		double varrho = 0;
		double meanstatvar = 0;
		double varstatvar = 0;


		// cycle over the sample
		for (int i=0; i<size; i++)	{
			cerr << '.';

			// get next point -> will be stored into "model", and thus, will be accessible through GetModel()

			GetNextPoint();

			GetModel()->UpdateBGC();

			meanchrono->Add(GetModel()->GetChronogram());

			double tmp1 = GetModel()->GetRho();
			double tmp2 = GetModel()->GetStatVar();
			mnrho += tmp1;
			varrho += tmp1 * tmp1;
			meanstatvar += tmp2;
			varstatvar += tmp2 * tmp2;
			GetModel()->GetBGCProcess()->specialUpdate();
			if (GetModel()->Split())	{
				GetModel()->GetSplitLengthTree()->specialUpdate();
			}
			meanbgc->Add(GetModel()->GetLogBGCProcess(), GetModel()->GetLengthTree(), 0, GetModel()->GetLogBGCOffset());
			// meanbgc->Add(GetModel()->GetBGCProcess(), GetModel()->GetLengthTree(), false);
			meansynrate->Add((NodeVarTree<Real>*) (GetModel()->GetSynRateTree()), GetModel()->GetLengthTree());

			if (GetModel()->FlexRho())	{
				meanrho->Add(GetModel()->GetRhoProcess(), GetModel()->GetLengthTree(), true);
				meanbranchrho->Add(GetModel()->GetRhoProcess());
			}

			CovMatrix& m = *(GetModel()->GetContMatrix());
			mat->Add(&m);

			for (int k=0; k<Ncont; k++)	{
				double tmp = GetModel()->GetRegCoef(k);
				meanreg[k] += tmp;
				varreg[k] += tmp * tmp;
				if (tmp > 0)	{
					ppreg[k]++;
				}
				tree[k]->Add(GetModel()->GetContProcess(), GetModel()->GetLengthTree(), k);

				double cov = 0;
				for (int j=0; j<Ncont; j++)	{
					cov += GetModel()->GetRegCoef(j) * m[k][j];
				}
				double varb = GetModel()->GetBGCVar();
				for (int j1=0; j1<Ncont; j1++)	{
					for (int j2=0; j2<Ncont; j2++)	{
						varb += GetModel()->GetRegCoef(j1) * GetModel()->GetRegCoef(j2) * m[j1][j2];
					}
				}
				double cor = cov / sqrt(varb * m[k][k]);
				correl[k] += cor;
			}

		}
		cerr << '\n';
		cerr << "normalise\n";

		mnrho /= size;
		varrho /= size;
		varrho -= mnrho * mnrho;
		meanstatvar /= size;
		varstatvar /= size;
		varstatvar -= meanstatvar * meanstatvar;
		cout << '\n';
		cout << "rho      : " << mnrho << '\t' << sqrt(varrho) << '\n';
		cout << "stat var : " << meanstatvar << '\t' << sqrt(varstatvar) << '\n';
		cout << '\n';

		meanchrono->Normalise();
		ofstream chos((GetName() + ".postmeandates.tre").c_str());
		meanchrono->ToStream(chos);

		ofstream cchos((GetName() + ".postmeandates.tab").c_str());
		meanchrono->Tabulate(cchos);

		for (int k=0; k<Ncont; k++)	{
			tree[k]->Normalise();
			ostringstream s;
			s << GetName() << ".postmean" << k+1 << ".tre";
			ofstream os(s.str().c_str());
			if (GetModel()->Split())	{
				tree[k]->Simplify();
			}
			tree[k]->ToStream(os);
			cerr << "reconstructed variation in continuous character # " << k+1 << " in "  << name << ".postmean" << k+1 << ".tre\n";
			cerr << "pp of mean leaf values > root value : " << tree[k]->GetPPLeafRoot() << '\n';
			if (tex)	{
				MeanChronoBubbleTree* textree = new MeanChronoBubbleTree(meanchrono,tree[k],xscale,yscale,nodescale,nodepower,barwidth,fontsize,bubbletext,withheader);
				textree->Draw((s.str() + ".tex").c_str());
			}
			if (postdist)	{
				ostringstream s;
				s << GetName() << ".postdist" << k+1 << ".tab";
				ofstream os(s.str().c_str());
				tree[k]->TabulateDistribution(os);
			}
		}

		if (GetModel()->FlexRho())	{
			meanrho->Normalise();
			ofstream rhoos((GetName() + ".postmeanrho.tre").c_str());
			meanrho->ToStream(rhoos);
			meanbranchrho->Normalise();
			ofstream branchrhoos((GetName() + ".postmeanbranchrho.tre").c_str());
			meanbranchrho->ToStream(branchrhoos);
		}

		meanbgc->Normalise();
		ofstream bgcos((GetName() + ".postmeanbgc.tre").c_str());
		if (GetModel()->Split())	{
			meanbgc->Simplify();
		}
		meanbgc->ToStream(bgcos);
		ofstream bbgcos((GetName() + ".postmeanbgc.tab").c_str());
		meanbgc->Tabulate(bbgcos);

		cerr << "BGC: pp of mean leaf values > root value : " << meanbgc->GetPPLeafRoot() << '\n';

		meansynrate->Normalise();
		ofstream sos((GetName() + ".postmeansynrate.tre").c_str());
		meansynrate->ToStream(sos);
		ofstream ssos((GetName() + ".postmeansynrate.tab").c_str());
		meansynrate->Tabulate(ssos);

		ofstream sumos((GetName() + ".slopes").c_str());
		for (int k=0; k<Ncont; k++)	{
			meanreg[k] /= size;
			varreg[k] /= size;
			ppreg[k] /= size;
			varreg[k] -= meanreg[k] * meanreg[k];
			correl[k] /= size;
			sumos << "mean reg: " << k << '\t' << meanreg[k] << " +/- " << sqrt(varreg[k]) << '\t' << ppreg[k] << '\t' << correl[k] << '\t' << correl[k] * correl[k] <<  '\n';
		}

		mat->Normalize();
		ofstream cout((GetName() + ".cov").c_str());
		cout << *mat;
		cerr << "covariance matrix in " << name << ".cov\n";
		cerr << '\n';

	}
};


int main(int argc, char* argv[])	{

	int burnin = 0;
	int every = 1;
	int until = -1;
	string name;

	string truefile = "";

	bool printlog = false;
	bool printmean = false;
	bool printci = true;
	bool printstdev = false;
	bool withleaf = true;
	bool withinternal = true;

	bool withheader = true;
	bool postdist = false;

	bool tex = false;
	double nodescale = 5;
	double nodepower = 1;
	double barwidth = 0.04;
	double x = 6;
	double y = 8;
	bool withnodebubbles = true;
	bool withtimeintervals = true;
	int fontsize = 4;

	bool bubbletext = false;

	bool variance = false;
	bool correl = false;
	bool ic = false;

	bool genebgc = false;

	int table = 0;
	string modelname = "modelname";

	try	{

		if (argc == 1)	{
			throw(0);
		}

		int i = 1;
		while (i < argc)	{
			string s = argv[i];

			if (s == "-tex")	{
				tex = true;
			}
			else if (s == "-gene")	{
				genebgc = true;
			}
			else if (s == "-table1")	{
				table = 1;
			}
			else if ((s == "-nuc") || (s == "-table2"))	{
				table = 2;
			}
			else if (s == "-name")	{
				i++;
				modelname = argv[i];
			}
			else if (s == "-h")	{
				withheader = false;
			}
			else if (s == "-postdist")	{
				postdist = true;
			}
			else if (s == "-ns")	{
				i++;
				nodescale = atof(argv[i]);
			}
			else if (s == "-np")	{
				i++;
				nodepower = atof(argv[i]);
			}
			else if (s == "-bw")	{
				i++;
				barwidth = atof(argv[i]);
			}
			else if (s == "+bt")	{
				bubbletext = true;
			}
			else if (s == "-bt")	{
				bubbletext = false;
			}
			else if (s == "-fs")	{
				i++;
				fontsize = atoi(argv[i]);
			}
			else if (s == "-xs")	{
				i++;
				x = atof(argv[i]);
			}
			else if (s == "-ys")	{
				i++;
				y = atof(argv[i]);
			}
			else if (s == "+nb")	{
				withnodebubbles = true;
			}
			else if (s == "-nb")	{
				withnodebubbles = false;
			}
			else if (s == "+ti")	{
				withtimeintervals = true;
			}
			else if (s == "-ti")	{
				withtimeintervals = false;
			}
			else if (s == "+log")	{
				printlog = true;
			}
			else if (s == "-log")	{
				printlog = false;
			}
			else if (s == "+mean")	{
				printmean = true;
			}
			else if (s == "-mean")	{
				printmean = false;
			}
			else if (s == "+stdev")	{
				printstdev = true;
			}
			else if (s == "-stdev")	{
				printstdev = false;
			}
			else if (s == "+ci")	{
				printci = true;
			}
			else if (s == "-ci")	{
				printci = false;
			}
			else if (s == "+leaf")	{
				withleaf = true;
			}
			else if (s == "-leaf")	{
				withleaf = false;
			}
			else if (s == "+internal")	{
				withinternal = true;
			}
			else if (s == "-internal")	{
				withinternal = false;
			}
			else if (s == "-var")	{
				variance = true;
			}
			else if (s == "-correl")	{
				correl = true;
			}
			else if (s == "-ic")	{
				ic = true;
			}
			else if ( (s == "-x") || (s == "-extract") )	{
				i++;
				if (i == argc) throw(0);
				s = argv[i];
				if (! IsInt(s))	{
					throw(0);
				}
				burnin = atoi(argv[i]);
				i++;
				if (i == argc) throw(0);
				s = argv[i];
				if (IsInt(s))	{
					every = atoi(argv[i]);
					i++;
					if (i == argc) throw(0);
					s = argv[i];
					if (IsInt(s))	{
						until = atoi(argv[i]);
					}
					else	{
						i--;
					}
				}
				else	{
					i--;
				}
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
		cerr << "readcoevol [-x <burnin> <every> <until>] <chainname> \n";
		cerr << '\n';
		exit(1);
	}

	BGCSample sample(name,burnin,every,until);

	if (table == 2)	{
		sample.ReadNuc(modelname);
	}
	else if (table == 1)	{
		sample.ReadTexTable1(modelname);
	}
	else if (genebgc)	{
		sample.ReadGeneBGC(printlog,printmean,printci,printstdev,withleaf,withinternal);
	}
	else if (correl)	{
		sample.ReadSuccessiveBranchesCorrel();
		// sample.ReadBranchCorrel(modelname);
	}
	else if (variance)	{
		sample.ReadVariance(modelname);
	}
	else if (ic)	{
		sample.ReadIC();
	}
	else	{
		sample.Read(printlog,printmean,printci,printstdev,withleaf,withinternal,tex,x,y,nodescale,nodepower,barwidth,fontsize,bubbletext,withheader,postdist);
	}
}



