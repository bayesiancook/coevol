
#pragma once

class ShiftedNodeAge : public Dvar<PosReal> {

    public:

    ShiftedNodeAge(Var<PosReal>* inage, Var<Real>* inlogNe, Var<RealVector>* inx, Var<PosReal>* inrootage, int inidxgentime)   {

        age = inage;
        logNe = inlogNe;
        x = inx;
        rootage = inrootage;
        idxgentime = inidxgentime;
        Register(age);
        Register(logNe);
        Register(x);
        specialUpdate();
    }

    // age is relative to root
    // shift node age by Ne * tau / rootage / 365.10^6
    void specialUpdate()    {
        if (age->val() > 1.0)   {
            cerr << "error in shifted node age: original age is not in (0,1)\n";
            cerr << age->val() << '\n';
            exit(1);
        }
        setval(age->val() + 2 * exp(logNe->val() + (*x)[idxgentime] - 6*log(10.0)) / rootage->val() / 365);
    }

    protected:

    Var<PosReal>* age;
    Var<Real>* logNe;
    Var<RealVector>* x;
    Var<PosReal>* rootage;
    int idxgentime;
};

class ShiftedChronogram : public NodeValPtrTree<Dvar<PosReal> >  {

    // given relative time of a node
    // and absolute root age
    // -> absolute 
    // at a given node:
    // generation time + Ne
    // -> 2 N_e * tau / rootage / 365.10^6: shift
    // however, shifted node should be older than daughter nodes

    public:

    ShiftedChronogram(NodeVarTree<PosReal>* inchronogram, NodeVarTree<Real>* inlogNe, NodeVarTree<RealVector>* inprocess, Var<PosReal>* inrootage, int inidxgentime)    {
        chronogram = inchronogram;
        logNe = inlogNe;
        process = inprocess;
        rootage = inrootage;
        idxgentime = inidxgentime;
		RecursiveCreate(GetRoot());
    }

	Tree* GetTree(){
		return chronogram->GetTree();
	}

    void specialUpdate()    {
        specialUpdate(GetRoot());
    }

	protected:

    void specialUpdate(const Link* from)    {
        GetNodeVal(from->GetNode())->specialUpdate();
        for (Link* link=from->Next(); link!=from; link=link->Next())    {
            specialUpdate(link->Out());
        }
    }

	Dvar<PosReal>* CreateNodeVal (const Link* link){
        return new ShiftedNodeAge(
                chronogram->GetNodeVal(link->GetNode()),
                logNe->GetNodeVal(link->GetNode()),
                process->GetNodeVal(link->GetNode()),
                rootage,
                idxgentime);
	}

    NodeVarTree<PosReal>* chronogram;
    NodeVarTree<Real>* logNe;
    NodeVarTree<RealVector>* process;
    Var<PosReal>* rootage;
    int idxgentime;
};

