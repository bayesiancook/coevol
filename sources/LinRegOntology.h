
#ifndef LINREGONTO_H
#define LINREGONTO_H

#include "LinRegOmega.h"
#include "Ontology.h"

class GOMean : public Dvar<Real>	{

	public:

	GOMean(Ontology* inontology, BidimIIDNormal* inbeta, int ingene, int incont, BidimIIDBernouilli* intoggle=0)	{
		ontology = inontology;
		beta = inbeta;
		gene = ingene;
		cont = incont;
		toggle = intoggle;
		beta->RegisterChildHorizontal(this,cont);
		if (toggle)	{
			toggle->RegisterChildHorizontal(this,cont);
		}
		specialUpdate();
	}

	GOMean(DAGnode* inroot)	{
		beta = 0;
		Register(inroot);
		specialUpdate();
	}

	Ontology* GetOntology()	{
		return ontology;
	}

	protected:

	void specialUpdate()	{
		if (beta)	{
			double total = 0;
			if (toggle)	{
				for (int l=0; l<beta->GetN(); l++)	{
					total += ontology->GetZ(gene,l) * beta->GetCell(l,cont)->val() * toggle->GetCell(l,cont)->val();
				}
			}
			else	{
				for (int l=0; l<beta->GetN(); l++)	{
					total += ontology->GetZ(gene,l) * beta->GetCell(l,cont)->val();
				}
			}
			setval(total);
		}
		else	{
			setval(0);
		}
	}

	Ontology* ontology;
	BidimIIDNormal* beta;
	BidimIIDBernouilli* toggle;
	int gene;
	int cont;
};

class BidimGOMean	{

	public:

	BidimGOMean(Ontology* inontology, BidimIIDNormal* inbeta, BidimIIDBernouilli* intoggle = 0, int inextra = 0, DAGnode* inroot = 0)	{
		ontology = inontology;
		beta = inbeta;
		toggle = intoggle;
		extra = inextra;
		root = inroot;
		if (ontology->GetNconcept() != beta->GetN())	{
			cerr << "error in bidim go mean\n";
			cerr << ontology->GetNconcept() << '\t' << beta->GetN() << '\n';
			exit(1);
		}
		Create();
	}

	~BidimGOMean()	{
		Delete();
	}

	int GetNgene() {return ontology->GetNgene();}
	int GetNconcept() {return ontology->GetNconcept();}
	int GetNcont() {return beta->GetP();}
	int GetP() {return GetNcont() + extra;}

	Ontology* GetOntology() {return ontology;}

	Dvar<Real>* GetMean(int gene, int cont)	{
		return mean[gene][cont];
	}

	protected:

	void Create()	{
		mean = new GOMean**[GetNgene()];
		for (int gene=0; gene<GetNgene(); gene++)	{
			mean[gene] = new GOMean*[GetNcont() + extra];
			for (int k=0; k<GetNcont(); k++)	{
				mean[gene][k] = new GOMean(ontology,beta,gene,k,toggle);
			}
			for (int k=GetNcont(); k<GetNcont()+extra; k++)	{
				mean[gene][k] = new GOMean(root);
			}
		}
	}

	void Delete()	{
		for (int gene=0; gene<GetNgene(); gene++)	{
			for (int k=0; k<GetNcont() + extra; k++)	{
				delete mean[gene][k];
			}
			delete[] mean[gene];
		}
		delete[] mean;
	}

	Ontology* ontology;
	BidimIIDNormal* beta;
	BidimIIDBernouilli* toggle;
	GOMean*** mean;
	int extra;
	DAGnode* root;
};

class BidimGONormal : public BidimArray<Real>	{

	public:

	BidimGONormal(BidimGOMean* inmean, Var<PosReal>** invar) :
			BidimArray<Real>(inmean->GetNgene(),inmean->GetP())	{

		mean = inmean;
		var = invar;
		Create();
	}

	~BidimGONormal()	{
		Delete();
	}

	Ontology* GetOntology() {return mean->GetOntology();}
	int GetNgene() {return mean->GetNgene();}
	int GetNcont() {return mean->GetP();}
	int GetNconcept() {return mean->GetNconcept();}

	protected:

	Rvar<Real>* CreateCell(int i, int j)	{
		return new Normal(mean->GetMean(i,j),var[j]);
	}

	BidimGOMean* mean;
	Var<PosReal>** var;
};

class BidimNormal : public BidimArray<Real>	{

	public:

	BidimNormal(BidimArray<Real>* inmean, BidimArray<PosReal>* invar) : BidimArray<Real>(inmean->GetN(), inmean->GetP())	{
		mean = inmean;
		var = invar;
		Create();
	}

	~BidimNormal()	{
		Delete();
	}

	protected:

	Rvar<Real>* CreateCell(int i, int j)	{
		return new Normal(mean->GetCell(i,j),var->GetCell(i,j));
	}

	BidimArray<Real>* mean;
	BidimArray<PosReal>* var;
};

class BidimIIDGamma : public BidimArray<PosReal>	{

	public:

	BidimIIDGamma(int inN, int inP, Var<PosReal>** inalpha, Var<PosReal>** inbeta) : BidimArray<PosReal>(inN,inP)	{
		alpha = inalpha;
		beta = inbeta;
		singlealpha = 0;
		singlebeta = 0;
		Create();
	}

	BidimIIDGamma(int inN, int inP, Var<PosReal>* inalpha, Var<PosReal>* inbeta) : BidimArray<PosReal>(inN,inP)	{
		alpha = 0;
		beta = 0;
		singlealpha = inalpha;
		singlebeta = inbeta;
		Create();
	}

	~BidimIIDGamma()	{
		Delete();
	}

	protected:

	Rvar<PosReal>* CreateCell(int i, int j)	{
		if (alpha)	{
			return new Gamma(alpha[j],beta[j]);
		}
		return new Gamma(singlealpha,singlebeta);
	}

	Var<PosReal>* singlealpha;
	Var<PosReal>* singlebeta;
	Var<PosReal>** alpha;
	Var<PosReal>** beta;
};

#endif

