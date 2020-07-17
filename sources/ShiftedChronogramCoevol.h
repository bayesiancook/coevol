
#pragma once

class ShiftedNodeAge : public Dvar<PosReal> {

    public:

    ShiftedNodeAge(Var<PosReal>* inage, Var<RealVector>* inx, int inidxpiS, int inidxdS)   {

        age = inage;
        x = inx;
        idxpiS = inidxpiS;
        idxdS = inidxdS;
        Register(age);
        Register(x);
        specialUpdate();
    }

    // segregation time in ancestral population:
    // delta = 2 * N_e * tau
    //
    // piS = 4 N_e u    ->  N_e = piS / 4 / u
    // u = dS * tau     ->  N_e = piS / dS / 4 / tau
    //                  ->  delta = 2 * N_e * tau = piS / dS / 2
    //
    void specialUpdate()    {
        if (age->val() > 1.0)   {
            cerr << "error in shifted node age: original age is not in (0,1)\n";
            cerr << age->val() << '\n';
            exit(1);
        }
        if (age->val() < 1e-4)  {
            setval(age->val());
        }
        else    {
            setval(age->val() + 0.5 * exp((*x)[idxpiS] - (*x)[idxdS]));
        }
    }

    protected:

    Var<PosReal>* age;
    Var<RealVector>* x;
    int idxpiS;
    int idxdS;
};

class ShiftedChronogram : public NodeValPtrTree<Dvar<PosReal> >  {

    public:

    ShiftedChronogram(NodeVarTree<PosReal>* inchronogram, NodeVarTree<RealVector>* inprocess, int inidxpiS, int inidxdS)    {
        chronogram = inchronogram;
        process = inprocess;
        idxpiS = inidxpiS;
        idxdS = inidxdS;
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
                process->GetNodeVal(link->GetNode()),
                idxpiS, idxdS);
	}

    NodeVarTree<PosReal>* chronogram;
    NodeVarTree<RealVector>* process;
    int idxpiS;
    int idxdS;
};

