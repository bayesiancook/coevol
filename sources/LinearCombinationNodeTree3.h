#ifndef NODELINEARCOMBINATIONTREE2_H
#define NODELINEARCOMBINATIONTREE2_H

#include "ValTree.h"


#include <cmath>
#include <sstream>

class SynRateLinearCombination : public Dvar<Real> {
	
	public :
	
	SynRateLinearCombination(Var<RealVector>* inx, Var<PosReal>* inrootage, int inL) {
		x = inx;
		rootage = inrootage;
        L = inL;
		Register(x);
		Register(rootage);
		specialUpdate();
	}
			
		
	void specialUpdate() {
        // dS = u / gentime per tree depth = u / # days per generation * # days per tree depth
        // rootage is assumed to be in Myr
        double a = (*x)[L] - (*x)[L+2] + log(rootage->val() * 365 * 1000 * 1000);
		setval(a);	
	}
		
					
	
	private :			
	
	Var<RealVector>* x;
	Var<PosReal>* rootage;
    int L;
};		

class SynRateLinearCombinationNodeTree : public NodeValPtrTree<Dvar<Real> > {
	
	public :
	
	SynRateLinearCombinationNodeTree(NodeVarTree<RealVector>* inprocess, Var<PosReal>* inrootage, int inL)  {
		process	= inprocess;
        rootage = inrootage;
        L = inL;
		RecursiveCreate(GetRoot());
	}
	
	~SynRateLinearCombinationNodeTree() {
		RecursiveDelete(GetRoot());
	}
	
	void specialUpdate()	{
		specialUpdate(GetRoot());
	}
	
	Tree* GetTree() {
		return process->GetTree();
	}	
		
	private :
	
	void specialUpdate(Link* from)	{
        GetNodeVal(from->GetNode())->specialUpdate();
		for(Link* link=from->Next(); link!=from; link=link->Next())	{
			specialUpdate(link->Out());
		}
	}

	Dvar<Real>* CreateNodeVal(const Link* link){
        return new SynRateLinearCombination(process->GetNodeVal(link->GetNode()), rootage, L);
	}
		
	
	NodeVarTree<RealVector>* process;
	Var<PosReal>* rootage;
    int L;
};	

class OmegaLinearCombination : public Dvar<Real> {
	
	public :
	
	
	OmegaLinearCombination(Var<RealVector>* inx, Var<Real>* ingamma, Var<Real>* inbeta, int inL) {
		x = inx;
		gamma = ingamma;
        beta = inbeta;
        L = inL;
        if ((L != 0) && (L != 1))   {
            cerr << "error in omega linear combination\n";
            exit(1);
        }
		Register(x);
		Register(gamma);
		Register(beta);
		specialUpdate();
	}
	
	void specialUpdate() {
        double a = -gamma->val() * (*x)[L+1] + beta->val();
        if (L == 1)  {
            double b = exp((*x)[0]);
            setval(log(a+b));
        }
        else    {
            setval(a);
        }
	}		
	
	private :
	
	Var<RealVector>* x;
	Var<Real>* gamma;
	Var<Real>* beta;
    int L;
};		

class OmegaLinearCombinationNodeTree : public NodeValPtrTree<Dvar<Real> > {
	
	public :
	
	OmegaLinearCombinationNodeTree(NodeVarTree<RealVector>* inprocess, Var<Real>* ingamma, Var<Real>* inbeta, int inL)  {
		process	= inprocess;
		gamma = ingamma;
        beta = inbeta;
        L = inL;
		RecursiveCreate(GetRoot());
	}
	
	~OmegaLinearCombinationNodeTree() {
		RecursiveDelete(GetRoot());
	}
	
	void specialUpdate()	{
		specialUpdate(GetRoot());
	}
	
	Tree* GetTree() {
		return process->GetTree();
	}	
		
	private :
	
	void specialUpdate(Link* from)	{
        GetNodeVal(from->GetNode())->specialUpdate();
		for(Link* link=from->Next(); link!=from; link=link->Next())	{
			specialUpdate(link->Out());
		}
	}

	Dvar<Real>* CreateNodeVal(const Link* link){
        return new OmegaLinearCombination(process->GetNodeVal(link->GetNode()), gamma, beta, L);
	}
		
	
	NodeVarTree<RealVector>* process;
	Var<Real>* gamma;
	Var<Real>* beta;
    int L;
};	

class LogNe: public Rvar<Real> {
	
	public :
	
	LogNe(Var<RealVector>* inx, Var<PosReal>* invar, int inL)  {
		x = inx;
        var = invar;
        L = inL;
		Register(x);
        Register(var);
        Sample();
	}
			
    void drawSample()   {
        /*
        double v = (*x)[L+1] + sqrt(var->val()) * Random::sNormal();
        setval(v);
        */
        setval((*x)[L+1]);
    }

    double logProb()    {
        return - 0.5 * (val() - (*x)[L+1]) * (val() - (*x)[L+1]) / var->val() - 0.5 * log(var->val());
	}
		
	private :			
	
	Var<RealVector>* x;
    Var<PosReal>* var;
    int L;
};		

class SynNumber : public Rvar<Int>  {

    public:

    SynNumber(Var<RealVector>* inx, Var<Real>* inlogNe, int inL, int inSNPnumber)    {
        x = inx;
        logNe = inlogNe;
        L = inL;
        SNPnumber = inSNPnumber;
        Register(x);
        Register(logNe);
        // Sample();
    }

    double getrate()    {

            double logpis = (*x)[L] + logNe->val() + log(4.0);
            double rate = SNPnumber * exp(logpis);
            return rate;
    }

    void drawSample()   {

        double rate = getrate();
        double p = exp(rate) * Random::Uniform();
        int n = 0;
        double ratio = 1;
        double total = 1;
        while (total<p)	{
            n++;
            ratio *=  rate / n;
            total += ratio;
        }
        setval(n);
    }

    double logProb()    {

        double rate = getrate();
        double ret = -rate + *this * log(rate);
        /*
        for (int i=2; i<=*this; i++)	{
            ret -= log((double) i);
        }
        */
        return ret;
    }

    private:

    Var<RealVector>* x;
    Var<Real>* logNe;
    int L;
    int SNPnumber;
};

class NonSynNumber : public Rvar<Int>  {

    public:

	NonSynNumber(Var<RealVector>* inx, Var<Real>* inlogNe, Var<Real>* ingamma, Var<Real>* inbeta, Var<Real>* indelta, int inL, int inSNPnumber) {
		x = inx;
        logNe = inlogNe;
		gamma = ingamma;
        beta = inbeta;
        delta = indelta;
        L = inL;
        SNPnumber = inSNPnumber;
		Register(x);
        Register(logNe);
		Register(gamma);
		Register(beta);
        Register(delta);
        // Sample();
	}
	
    double getrate()    {
        // piN/piS = b N_e ** -gamma
        // where b = beta * NI if model is constrained on same sequence, and beta otherwise
        double b = beta->val();
        if (delta) {
            b += log(delta->val());
        }
        double logpin = (*x)[L] + log(4.0) + (1-gamma->val()) * logNe->val() + b;
        double rate = SNPnumber * exp(logpin);
        return rate;
    }

    void drawSample()   {

        double rate = getrate();
        double p = exp(rate) * Random::Uniform();
        int n = 0;
        double ratio = 1;
        double total = 1;
        while (total<p)	{
            n++;
            ratio *=  rate / n;
            total += ratio;
        }
        setval(n);
    }

    double logProb()    {

        double rate = getrate();
        double ret = -rate + *this * log(rate);
        /*
        for (int i=2; i<=*this; i++)	{
            ret -= log((double) i);
        }
        */
        return ret;
    }

    private:

	Var<RealVector>* x;
    Var<Real>* logNe;
	Var<Real>* gamma;
	Var<Real>* beta;
    Var<Real>* delta;
    int L;
    int SNPnumber;
};

/*

class LogPiS: public Rvar<Real> {
	
	public :
	
	LogPiS(Var<RealVector>* inx, Var<PosReal>* invar, int inL)  {
		x = inx;
        var = invar;
        L = inL;
		Register(x);
        Register(var);
        Sample();
	}
			
    double getmean()    {
        double a = (*x)[L] + (*x)[L+1] + log(4.0);
        return a;
    }

    void drawSample()   {
        double x = getmean() + sqrt(var->val()) * Random::sNormal();
        setval(x);
    }

    double logProb()    {
        double a = getmean();
        return - 0.5 * (val() - a) * (val() - a) / var->val() - 0.5 * log(var->val());
	}
		
	private :			
	
	Var<RealVector>* x;
    Var<PosReal>* var;
    int L;
};		

class LogPiNPiS: public Rvar<Real> {
	
	public :
	
	LogPiNPiS(Var<RealVector>* inx, Var<Real>* ingamma, Var<Real>* inbeta, Var<PosReal>* inNI, Var<PosReal>* invar, int inL) {
		x = inx;
		gamma = ingamma;
        beta = inbeta;
        NI = inNI;
        var = invar;
        L = inL;
		Register(x);
		Register(gamma);
		Register(beta);
        Register(NI);
        Register(var);
        Sample();
	}
	
    double getmean()    {
        // piN/piS = b N_e ** -gamma
        // where b = beta * NI if model is constrained on same sequence, and beta otherwise
        double b = beta->val();
        if (NI) {
            b += log(NI->val());
        }
        double a = -gamma->val() * (*x)[L+1] + b;
        return a;
    }

    void drawSample()   {
        double x = getmean() + sqrt(var->val()) * Random::sNorm();
        setval(x);
    }

    double logProb()    {
        double a = getmean();
        return - 0.5 * (val() - a) * (val() - a) / var->val() - 0.5 * log(var->val());
	}
		
	private :
	
	Var<RealVector>* x;
	Var<Real>* gamma;
	Var<Real>* beta;
    Var<PosReal>* NI;
    Var<PosReal>* var;
    int L;
};		
class PiSLinearCombination : public Dvar<Real> {
	
	public :
	
	PiSLinearCombination(Var<RealVector>* inx, int inL)  {
		x = inx;
        L = inL;
		Register(x);
		specialUpdate();
	}
			
		
	void specialUpdate() {
        // piS = 4 Ne u
        double a = (*x)[L] + (*x)[L+1] + log(4.0);
		setval(a);	
	}
		
	private :			
	
	Var<RealVector>* x;
    int L;
};		

class PiSLinearCombinationNodeTree : public NodeValPtrTree<Dvar<Real> > {
	
	public :
	
	PiSLinearCombinationNodeTree(NodeVarTree<RealVector>* inprocess, int inL)  {
		process	= inprocess;
        L = inL;
		RecursiveCreate(GetRoot());
	}
	
	~PiSLinearCombinationNodeTree() {
		RecursiveDelete(GetRoot());
	}
	
	void specialUpdate()	{
		specialUpdate(GetRoot());
	}
	
	Tree* GetTree() {
		return process->GetTree();
	}	
		
	private :
	
	void specialUpdate(Link* from)	{
        GetNodeVal(from->GetNode())->specialUpdate();
		for(Link* link=from->Next(); link!=from; link=link->Next())	{
			specialUpdate(link->Out());
		}
	}

	Dvar<Real>* CreateNodeVal(const Link* link){
        return new SynRateLinearCombination(process->GetNodeVal(link->GetNode()), L);
	}
		
	
	NodeVarTree<RealVector>* process;
    int L;
};	

class PiNPiSLinearCombination : public Dvar<Real> {
	
	public :
	
	
	PiNPiSLinearCombination(Var<RealVector>* inx, Var<Real>* ingamma, Var<Real>* inbeta, Var<PosReal>* inNI, int inL) {
		x = inx;
		gamma = ingamma;
        beta = inbeta;
        NI = inNI;
        L = inL;
		Register(x);
		Register(gamma);
		Register(beta);
        Register(NI);
		specialUpdate();
	}
	
	void specialUpdate() {
        // piN/piS = b N_e ** -gamma
        // where b = beta * NI if model is constrained on same sequence, and beta otherwise
        double b = inbeta->val();
        if (NI) {
            b += log(NI->val());
        }
        double a = -gamma->val() * (*x)[L+1] + b;
        setval(a);
	}		
	
	private :
	
	Var<RealVector>* x;
	Var<Real>* gamma;
	Var<Real>* beta;
    Var<PosReal>* NI;
    int L;
};		

class PiNPiSLinearCombinationNodeTree : public NodeValPtrTree<Dvar<Real> > {
	
	public :
	
	PiNPiSLinearCombinationNodeTree(NodeVarTree<RealVector>* inprocess, Var<Real>* ingamma, Var<Real>* inbeta, Var<PosReal>* inNI, int inL)  {
		process	= inprocess;
		gamma = ingamma;
        beta = inbeta;
        NI = inNI;
        L = inL;
		RecursiveCreate(GetRoot());
	}
	
	~PiNPiSLinearCombinationNodeTree() {
		RecursiveDelete(GetRoot());
	}
	
	void specialUpdate()	{
		specialUpdate(GetRoot());
	}
	
	Tree* GetTree() {
		return process->GetTree();
	}	
		
	private :
	
	void specialUpdate(Link* from)	{
        GetNodeVal(from->GetNode())->specialUpdate();
		for(Link* link=from->Next(); link!=from; link=link->Next())	{
			specialUpdate(link->Out());
		}
	}

	Dvar<Real>* CreateNodeVal(const Link* link){
        return new PiNPiSLinearCombination(process->GetNodeVal(link->GetNode()), gamma, beta, NI, L);
	}
		
	
	NodeVarTree<RealVector>* process;
	Var<Real>* gamma;
	Var<Real>* beta;
    Var<PosReal>* NI;
    int L;
};	
*/

#endif	
	
	
	
	
	
