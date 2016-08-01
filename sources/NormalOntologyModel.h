
#ifndef NORMALONTOLOGY_H
#define NORMALONTOLOGY_H

#include "RandomTypes.h"
#include "Jeffreys.h"
#include "IID.h"
#include "ProbModel.h"
#include "LinRegOmega.h"
#include "LinRegOntology.h"
#include "MatrixAlgebra.h"

// const double minjeff = 0.0000001;
// const double maxjeff = 10000000;
const double minjeff = 1e-6;
const double maxjeff = 1e-2;
const double minuni = -1000;
const double maxuni = 1000;

class NormalOntologyModel : public ProbModel, public virtual MatrixAlgebra {

	public:

	string ontologyfile;
	string datafile;

	int Ngene;
	int Ncont;
	int mingene;
	int minconcept;

	string* genename;
	double** genemean;
	double** genevar;

	Ontology* ontology;

	Const<Real>* Zero;
	Const<PosReal>* One;

	JeffreysIIDArray* jeffkappa;
	Var<PosReal>** kappa;
	Var<Real>** betamean;
	BidimIIDNormal* beta;

	Beta** theta;
	BidimIIDBernouilli* toggle;

	BidimGOMean* alphamean;

	JeffreysIIDArray* jeffalphavar;
	Var<PosReal>** alphavar;
	BidimGONormal* alpha;

	BidimNormal* GeneMean;
	BidimIIDGamma* GeneVar;

	double** Lambda;
	double** alphaz;
	double** betabar;
	double** alphahat;
	double** alphabar;
	CovMatrix* CovQ;
	double** Q;
	double** InvQ;
	double** InvMZ;
	CovMatrix* CovM;
	double** M;
	double** InvM;
	double* conjbeta;
	double* conjalpha;

	int nactive;
	int* index;
	CovMatrix* CovredM;
	double** redM;
	double** InvredM;
	double** redLambda;

	bool fixvar;

	int GetNgene()	{
		return Ngene;
	}

	int GetNcont()	{
		return Ncont;
	}

	int GetNconcept()	{
		return ontology->GetNconcept();
	}

	Ontology* GetOntology()	{
		return ontology;
	}

	NormalOntologyModel(string indatafile, string inontologyfile, int inmingene, int inminconcept, int withtoggle, double inlogvar)	{

		cerr << "create model\n";

		fixvar = (inlogvar != -100);

		datafile = indatafile;
		ontologyfile = inontologyfile;
		mingene = inmingene;
		minconcept = inminconcept;

		cerr << "datafile : " << datafile << '\n';
		ifstream is(datafile.c_str());
		int tmpNgene;
		is >> tmpNgene >> Ncont;
		cerr << "Ngene : " << tmpNgene << '\n';
		cerr << "Ncont : " << Ncont << '\n';
		string* tmpgenename = new string[tmpNgene];
		double** tmpgenemean = new double*[tmpNgene];
		double** tmpgenevar = new double*[tmpNgene];
		for (int gene=0; gene<tmpNgene; gene++)	{
			is >> tmpgenename[gene];
			tmpgenemean[gene] = new double[GetNcont()];
			tmpgenevar[gene] = new double[GetNcont()];
			double tmp1, tmp2;
			for (int cont=0; cont<GetNcont(); cont++)	{
				is >> tmp1 >> tmp2;
				tmpgenemean[gene][cont] = tmp1;
				tmpgenevar[gene][cont] = tmp2 * tmp2;
			}
		}

		int* include = new int[tmpNgene];

		Ontology* tmpontology = new Ontology(ontologyfile);
		ontology = new Ontology(tmpontology,tmpNgene,tmpgenename,mingene,minconcept,include);
		Ngene = 0;
		for (int g=0; g<tmpNgene; g++)	{
			Ngene += include[g];
		}
		if (Ngene != ontology->GetNgene())	{
			cerr << "error : non matching reduxed gene sets\n";
			exit(1);
		}
		genename = new string[Ngene];
		genemean = new double*[Ngene];
		genevar = new double*[Ngene];
		int gene = 0;
		for (int g=0; g<tmpNgene; g++)	{
			if (include[g])	{
				genename[gene] = tmpgenename[g];
				genemean[gene] = tmpgenemean[g];
				genevar[gene] = tmpgenevar[g];
				gene++;
			}
			else	{
				delete[] tmpgenemean[g];
				delete[] tmpgenevar[g];
			}
		}
		delete[] tmpgenename;
		delete[] tmpgenemean;
		delete[] tmpgenevar;

		delete[] include;

		cerr << "after reduction: " << ontology->GetNgene() << '\t' << ontology->GetNconcept() << '\n';
		cerr << GetNgene() << '\t' << GetNconcept() << '\t' << GetNcont() << '\n';
		for (int gene=0; gene<GetNgene(); gene++)	{
			cerr << genename[gene];
			for (int k=0; k<GetNcont(); k++)	{
				cerr << '\t' << genemean[gene][k] << '\t' << genevar[gene][k];
			}
			cerr << '\n';
		}
		cerr << '\n';

		Zero = new Const<Real>(0);
		One = new Const<PosReal>(1);

		jeffkappa = new JeffreysIIDArray(GetNcont(),minjeff,maxjeff,Zero);
		kappa = new Var<PosReal>*[GetNcont()];
		betamean = new Var<Real>*[GetNcont()];
		for (int cont=0; cont<GetNcont(); cont++)	{
			kappa[cont] = jeffkappa->GetVal(cont);
			if (fixvar)	{
				jeffkappa->GetVal(cont)->setval(exp(inlogvar* log(10.0)));
				jeffkappa->GetVal(cont)->Clamp();
			}
			else	{
				jeffkappa->GetVal(cont)->setval(0.000002);
			}
			betamean[cont] = Zero;
		}

		beta = new BidimIIDNormal(GetNconcept(),GetNcont(),betamean,kappa);

		theta = 0;
		toggle = 0;
		if (withtoggle)	{
			theta = new Beta*[GetNcont()];
			for (int cont=0; cont<GetNcont(); cont++)	{
				theta[cont] = new Beta(One,One);
			}
			toggle = new BidimIIDBernouilli(GetNconcept(),GetNcont(),theta);
			for (int cont=0; cont<GetNcont(); cont++)	{
				toggle->SetAt(1.0,cont);
			}
		}

		alphamean = new BidimGOMean(ontology,beta,toggle);

		jeffalphavar = new JeffreysIIDArray(GetNcont(),minjeff,maxjeff,Zero);
		alphavar = new Var<PosReal>*[GetNcont()];
		for (int cont=0; cont<GetNcont(); cont++)	{
			alphavar[cont] = jeffalphavar->GetVal(cont);
			jeffalphavar->GetVal(cont)->setval(0.001);
		}

		alpha = new BidimGONormal(alphamean,alphavar);
		/*
		for (int k=0; k<GetNcont(); k++)	{
			for (int i=0; i<GetNgene(); i++)	{
				alpha->GetCell(i,k)->setval(genemean[i][k]);
				alpha->GetCell(i,k)->Clamp();
			}
		}
		*/

		GeneVar = new BidimIIDGamma(GetNgene(),GetNcont(),One,One);
		GeneVar->ClampAt(genevar);
		GeneMean = new BidimNormal(alpha,GeneVar);
		GeneMean->ClampAt(genemean);

		RootRegister(Zero);
		RootRegister(One);
		Register();

		Update();

		MakeScheduler();

		CreateConjugateMatrices();

		cerr << "model created\n";
	}

	~NormalOntologyModel()	{}

	void MakeScheduler();
	/*
	double Move(double tuning = 1)	{
		scheduler.Cycle(1,1,true,false);
		return 1;
	}
	*/

	bool withToggle()	{
		return toggle;
	}

	BidimIIDNormal* GetBeta() {return beta;}
	BidimGONormal* GetAlpha() {return alpha;}
	BidimIIDBernouilli* GetToggle() {return toggle;}

	Beta* GetTheta(int cont) {return theta[cont];}

	double GetBeta(int concept, int cont)	{
		if (toggle)	{
			return beta->GetCell(concept,cont)->val() * toggle->GetCell(concept,cont)->val();
		}
		return beta->GetCell(concept,cont)->val();
	}

	double GetMeanBeta(int cont)	{
		double mean = 0;
		for (int k=0; k<GetNconcept(); k++)	{
			mean += GetBeta(k,cont);
		}
		mean /= GetNconcept();
		return mean;
	}

	double GetVarBeta(int cont)	{
		double mean = 0;
		double var = 0;
		int tot = 0;
		for (int k=0; k<GetNconcept(); k++)	{
			double tmp = GetBeta(k,cont);
			if ((!toggle) || (toggle->GetCell(k,cont)->val()))	{
				tot++;
			}
			mean += tmp;
			var += tmp * tmp;
		}
		mean /= GetNconcept();
		var /= GetNconcept();
		var -= mean * mean;
		if (!tot)	{
			return 0;
		}
		return var;
	}

	double GetAlpha(int gene, int cont)	{
		return alpha->GetCell(gene,cont)->val();
	}

	void drawSample()	{
		cerr << "in normal ontology draw sample\n";
		exit(1);
	}

	double GetLogProb()	{
		double total = 0;
		total += jeffkappa->GetLogProb();
		if (toggle)	{
			for (int cont=0; cont<GetNcont(); cont++)	{
				for (int k=0; k<GetNconcept(); k++)	{
					if (GetBeta(k,cont))	{
						total += beta->GetCell(k,cont)->GetLogProb();
					}
				}
			}
		}
		else	{
			total += beta->GetLogProb();
		}
		if (toggle)	{
			for (int cont=0; cont<GetNcont(); cont++)	{
				total += theta[cont]->GetLogProb();
			}
			total += toggle->GetLogProb();
		}
		total += jeffalphavar->GetLogProb();
		total += alpha->GetLogProb();
		total += GeneMean->GetLogProb();
		return total;
	}

	void TraceHeader(ostream& os)	{
		os << "logl";
		for (int cont=0; cont<GetNcont(); cont++)	{
			os << "\tbetamean" << cont;
			os << "\tbetavar" << cont;
			os << "\tjeffbeta" << cont;
			if (toggle)	{
				os << "\ttoggle" << cont;
			}
			os << "\talphavar" << cont;

		}
		os << '\n';
	}

	void Trace(ostream& os)	{
		os << GetLogProb();
		for (int cont=0; cont<GetNcont(); cont++)	{
			os << '\t' << GetMeanBeta(cont) << '\t' << GetVarBeta(cont);
			// os << '\t' << beta->GetMean(cont) << '\t' << beta->GetVar(cont);
			os << '\t' << log(jeffkappa->GetVal(cont)->val());
			if (toggle)	{
				os << '\t' << toggle->GetMean(cont);
			}
			os << '\t' << log(jeffalphavar->GetVal(cont)->val());
		}
		os << '\n';
	}

	void ToStream(ostream& os)	{
		jeffkappa->ToStream(os);
		os << '\n';
		jeffalphavar->ToStream(os);
		os << '\n';
		if (toggle)	{
			for (int cont=0; cont<GetNcont(); cont++)	{
				os << *theta[cont] << '\t';
			}
			os << '\n';
			toggle->ToStream(os);
		}
		beta->ToStream(os);
		os << '\n';
		alpha->ToStream(os);
		os << '\n';
	}

	void FromStream(istream& is)	{
		jeffkappa->FromStream(is);
		jeffalphavar->FromStream(is);
		if (toggle)	{
			for (int cont=0; cont<GetNcont(); cont++)	{
				is >> *theta[cont];
			}
			toggle->FromStream(is);
		}
		beta->FromStream(is);
		alpha->FromStream(is);
	}

	void CreateConjugateMatrices()	{
		Lambda = ontology->GetLambda();
		alphaz = MatrixCreate(GetNconcept(),1);
		betabar = MatrixCreate(GetNconcept(),1);
		CovM = new CovMatrix(GetNconcept());
		M = CovM->GetMatrix();
		MatrixSetIdentity(M,GetNconcept());
		CovM->CorruptDiag();
		CovM->Diagonalise();
		InvM = CovM->GetInvMatrix();
		conjbeta = new double[GetNconcept()];

		redLambda = MatrixCreate(GetNconcept(),GetNconcept());
		index = new int[GetNconcept()];
		CovredM = 0;

		/*
		alphahat = MatrixCreate(GetNgene(),1);
		alphabar = MatrixCreate(GetNgene(),1);
		CovQ = new CovMatrix(GetNgene());
		Q = CovQ->GetMatrix();
		MatrixSetIdentity(Q,GetNgene());
		CovQ->CorruptDiag();
		CovQ->Diagonalise();
		InvQ = CovQ->GetInvMatrix();
		conjalpha = new double[GetNgene()];

		InvMZ = MatrixCreate(GetNconcept(),GetNgene());
		*/

	}

	void ResampleAlpha()	{
		for (int i=0; i<GetNgene(); i++)	{
			for (int j=0; j<GetNcont(); j++)	{
				double tau0 = 1.0 / alphavar[j]->val();
				double tau1 = 1.0 / genevar[i][j];
				double tau2 = tau0 + tau1;
				double m = (tau0 * alphamean->GetMean(i,j)->val() + tau1 * genemean[i][j]) / (tau0 + tau1);
				double v = 1.0 / tau2;
				double postalpha = m + Random::sNormal() * sqrt(v);
				alpha->GetCell(i,j)->setval(postalpha);
			}
		}
	}

	int ComputeToggleM(int cont)	{

		// count number of entries that are currently active
		nactive = 0;
		for (int k=0; k<GetNconcept(); k++)	{
			if (toggle->GetCell(k,cont)->val())	{
				index[nactive] = k;
				nactive++;
			}
		}

		for (int k=0; k<nactive; k++)	{
			for (int l=0; l<nactive; l++)	{
				redLambda[k][l] = Lambda[index[k]][index[l]];
			}
		}

		// M = tau Lambda + kappa I (Nconcept x Nconcept)
		delete CovredM;
		CovredM = new CovMatrix(nactive);
		redM = CovredM->GetMatrix();
		MatrixSetIdentity(redM,nactive);
		MatrixScalarProduct(redM,alphavar[cont]->val()/kappa[cont]->val(),nactive,nactive);
		MatrixAdd(redM,redLambda,nactive,nactive);
		MatrixScalarProduct(redM,1.0/alphavar[cont]->val(),nactive,nactive);
		CovredM->CorruptDiag();
		int failed = CovredM->Diagonalise();
		InvredM = CovredM->GetInvMatrix();
		return failed;
	}

	void ComputeToggleBetaBar(int cont)	{

		for (int k=0; k<nactive; k++)	{
			alphaz[k][0] = 0;
			for (int i=0; i<GetNgene(); i++)	{
				alphaz[k][0] += GetAlpha(i,cont) * ontology->GetZ(i,index[k]);
			}
		}
		MatrixScalarProduct(alphaz,1.0/alphavar[cont]->val(),nactive,1);
		MatrixProduct(InvredM,alphaz,betabar,nactive,nactive,1);
	}

	void ComputeM(int cont)	{

		// M = tau Lambda + kappa I (Nconcept x Nconcept)
		MatrixSetIdentity(M,GetNconcept());
		MatrixScalarProduct(M,alphavar[cont]->val()/kappa[cont]->val(),GetNconcept(),GetNconcept());
		MatrixAdd(M,Lambda,GetNconcept(),GetNconcept());
		MatrixScalarProduct(M,1.0/alphavar[cont]->val(),GetNconcept(),GetNconcept());
		CovM->CorruptDiag();
		CovM->Diagonalise();
	}

	void ComputeBetaBar(int cont)	{

		// betabar = M^-1 alpha Z
		for (int k=0; k<GetNconcept(); k++)	{
			alphaz[k][0] = 0;
			for (int i=0; i<GetNgene(); i++)	{
				alphaz[k][0] += GetAlpha(i,cont) * ontology->GetZ(i,k);
			}
		}
		MatrixScalarProduct(alphaz,1.0/alphavar[cont]->val(),GetNconcept(),1);

		MatrixProduct(InvM,alphaz,betabar,GetNconcept(),GetNconcept(),1);
	}

	/*
	void ComputeQ(int cont)	{

		// Q = D + tau^2  Z' M-1 Z + tau I (Ngene x Ngene)
		// where D is gene var diagonal matrix for cont

		double** Z = ontology->GetZmatrix();
		double** transZ = ontology->GetTransZmatrix();

		double tau = 1.0 / jeffalphavar->GetVal(cont)->val();

		MatrixProduct(InvM,Z,InvMZ,GetNconcept(),GetNconcept(),GetNgene());
		MatrixProduct(transZ,InvMZ,Q,GetNgene(),GetNconcept(),GetNgene());
		MatrixScalarProduct(Q,tau*tau,GetNgene(),GetNgene());

		for (int i=0; i<GetNgene(); i++)	{
			Q[i][i] += tau + 1.0 / genevar[i][cont];
		}

		CovQ->CorruptDiag();
		CovQ->Diagonalise();
	}

	void ComputeAlphaBar(int cont)	{

		// alphabar = R^-1 D alphahat

		for (int i=0; i<GetNgene(); i++)	{
			alphahat[i][0] = genemean[i][cont] / genevar[i][cont];
		}

		MatrixProduct(InvQ,alphahat,alphabar,GetNgene(),GetNgene(),1);
	}
	*/

	double CheckMatrixInverse(double** M, double** InvM, int N)	{
		double max = 0;
		for (int i=0; i<N; i++)	{
			for (int j=0; j<N; j++)	{
				double tot = 0;
				for (int k=0; k<N; k++)	{
					tot += M[i][k] * InvM[k][j];
				}
				if (i==j)	{
					tot -= 1;
				}
				if (max < fabs(tot))	{
					max = fabs(tot);
				}
			}
		}
		return max;
	}

	/*
	double MarginalLogProb(int cont)	{

		ComputeM(cont);
		ComputeQ(cont);
		ComputeAlphaBar(cont);

		// alphabar = gene mean
		// S = alphabar' R alphabar
		double S = 0;
		for (int i=0; i<GetNgene(); i++)	{
			for (int j=0; j<GetNgene(); j++)	{
				S += alphabar[i][0] * Q[i][j] * alphabar[j][0];
			}
		}

		double tau = 1.0 / jeffalphavar->GetVal(cont)->val();
		double kappa = 1.0 / jeffkappa->GetVal(cont)->val();

		double total = 0;
		total += (GetNgene() - 1) * log(tau);
		total += (GetNconcept() - 1) * log(kappa);
		total -= CovM->GetLogDeterminant();
		total -= CovQ->GetLogDeterminant();
		total -= S;

		return 0.5 * total;
	}

	double MHResampleAlphaVar(double tuning, int cont)	{

		double bk = jeffalphavar->GetVal(cont)->val();
		double m = tuning * (Random::Uniform() - 0.5);
		double e = exp(m);
		double logratio = -MarginalLogProb(cont);
		double newval = bk*e;
		jeffalphavar->GetVal(cont)->setval(newval);
		logratio += MarginalLogProb(cont);
		logratio += m;
		int acc = ((newval > minjeff) && (newval < maxjeff) && (log(Random::Uniform()) < logratio));
		if (! acc)	{
			jeffalphavar->GetVal(cont)->setval(bk);
		}
		return acc;
	}

	double MHResampleBetaVar(double tuning, int cont)	{

		double bk = jeffkappa->GetVal(cont)->val();
		double m = tuning * (Random::Uniform() - 0.5);
		double e = exp(m);
		double logratio = -MarginalLogProb(cont);
		double newval = bk*e;
		jeffkappa->GetVal(cont)->setval(newval);
		logratio += MarginalLogProb(cont);
		logratio += m;
		int acc = ((newval > minjeff) && (newval < maxjeff) && (log(Random::Uniform()) < logratio));
		if (! acc)	{
			jeffkappa->GetVal(cont)->setval(bk);
		}
		return acc;
	}

	double ResampleAlphaBeta(double tuning = 1, int nrep = 0)	{

		double acc = 0;
		double tot = 0;
		// model is in fact split into independent components
		for (int cont=0; cont<GetNcont(); cont++)	{

			for (int rep=0; rep<nrep; rep++)	{
				acc += MHResampleAlphaVar(tuning,cont);
				tot++;
				acc += MHResampleAlphaVar(0.1 * tuning,cont);
				tot++;
				acc += MHResampleBetaVar(tuning,cont);
				tot++;
				acc += MHResampleBetaVar(0.1 * tuning,cont);
				tot++;
			}

			ComputeM(cont);
			ComputeQ(cont);

			ComputeAlphaBar(cont);

			CovQ->drawValInv(conjalpha);
			for (int i=0; i<GetNgene(); i++)	{
				alpha->GetCell(i,cont)->setval(conjalpha[i] + alphabar[i][0]);
			}

			ComputeBetaBar(cont);

			CovM->drawValInv(conjbeta);
			for (int k=0; k<GetNconcept(); k++)	{
				beta->GetCell(k,cont)->setval(conjbeta[k] + betabar[k][0]);
			}
		}
		return acc / tot;
	}
	*/

	void ResampleBeta(int cont)	{

		ComputeM(cont);
		ComputeBetaBar(cont);

		CovM->drawValInv(conjbeta);

		for (int k=0; k<GetNconcept(); k++)	{
			beta->GetCell(k,cont)->setval(conjbeta[k] + betabar[k][0]);
		}
	}

	double LogProb(int cont)	{

		double total = 0;
		total += nactive * log(jeffkappa->GetVal(cont)->val());
		total += CovredM->GetLogDeterminant();
		for (int k=0; k<nactive; k++)	{
			for (int l=0; l<nactive; l++)	{
				total -= betabar[k][0] * redM[k][l] * betabar[l][0];
			}
		}
		return -0.5 * total;
	}

	int ResampleToggleBeta(int cont, int nrep)	{

		double logprob = 0;
		int failed =0;

		if (nrep)	{
			failed = ComputeToggleM(cont);
			ComputeToggleBetaBar(cont);

			double logprob = LogProb(cont);
		}

		for (int rep=0; rep<nrep; rep++)	{
			for (int k=0; k<GetNconcept(); k++)	{
				int bkstate = toggle->GetCell(k,cont)->val();
				double bklogprob = logprob;
				double bknactive = nactive;

				toggle->GetCell(k,cont)->setval(1 - bkstate);

				ComputeToggleM(cont);
				ComputeToggleBetaBar(cont);

				logprob = LogProb(cont);

				double delta = logprob - bklogprob;

				/*
				if (! bkstate)	{
					delta += log(theta[cont]->val()) - log(1 - theta[cont]->val());
				}
				else	{
					delta -= log(theta[cont]->val()) - log(1 - theta[cont]->val());
				}
				*/
				if (! bkstate)	{
					delta += log((double) (bknactive+1)) - log((double) (GetNconcept() - bknactive));
				}
				else	{
					delta += log((double) (GetNconcept() - bknactive + 1)) - log((double) bknactive);
				}
				int acc = (log(Random::Uniform()) < delta);
				if (! acc)	{
					toggle->GetCell(k,cont)->setval(bkstate);
					logprob = bklogprob;
					nactive = bknactive;
				}
			}
		}

		failed = ComputeToggleM(cont);
		if (! failed)	{
		ComputeToggleBetaBar(cont);
		CovredM->drawValInv(conjbeta);

		for (int k=0; k<GetNconcept(); k++)	{
			// beta->GetCell(k,cont)->setval(0);
			beta->GetCell(k,cont)->Sample();
		}

		for (int k=0; k<nactive; k++)	{
			beta->GetCell(index[k],cont)->setval(conjbeta[k] + betabar[k][0]);
		}
		}
		else	{
			return 0;
		}

		double x = Random::sGamma(1.0 + nactive);
		double y = Random::sGamma(1.0 + GetNconcept() - nactive);
		theta[cont]->setval(x / (x+y));
		return 1;
	}

	void ResampleAlphaVar()	{

		for (int cont=0; cont<GetNcont(); cont++)	{
			jeffalphavar->GetVal(cont)->Corrupt(false);
			double s2 = 0;
			for (int i=0; i<GetNgene(); i++)	{
				double tot = 0;
				for (int k=0; k<GetNconcept(); k++)	{
					tot += ontology->GetZ(i,k) * GetBeta(k,cont);
				}
				double tmp = alphamean->GetMean(i,cont)->val() - alpha->GetCell(i,cont)->val();
				s2 += tmp*tmp;
			}
			double tmp = Random::Gamma(0.5 * GetNgene(), 0.5 * s2);
			double var = 1.0 / tmp;
			if ((var > minjeff) && (var < maxjeff))	{
				jeffalphavar->GetVal(cont)->setval(var);
			}
			jeffalphavar->GetVal(cont)->Update();
		}
	}

	void ResampleBetaVar()	{

		for (int cont = 0; cont<GetNcont(); cont++)	{
			jeffkappa->GetVal(cont)->Corrupt(false);
			double s2 = 0;
			for (int k=0; k<GetNconcept(); k++)	{
				double tmp = beta->GetCell(k,cont)->val();
				s2 += tmp * tmp;
			}
			double tmp = Random::Gamma(0.5 * GetNconcept(), 0.5 * s2);
			double var = 1.0 / tmp;
			if ((var > minjeff) && (var < maxjeff))	{
				jeffkappa->GetVal(cont)->setval(var);
			}
			jeffkappa->GetVal(cont)->Update();
		}
	}

	void ResampleToggleBetaVar()	{

		for (int cont = 0; cont<GetNcont(); cont++)	{
			jeffkappa->GetVal(cont)->Corrupt(false);
			double s2 = 0;
			int nactive = 0;
			for (int k=0; k<GetNconcept(); k++)	{
				if (GetBeta(k,cont))	{
					double tmp = beta->GetCell(k,cont)->val();
					s2 += tmp * tmp;
					nactive++;
				}
			}
			if (nactive)	{
				double tmp = Random::Gamma(0.5 * nactive, 0.5 * s2);
				double var = 1.0 / tmp;
				if ((var > minjeff) && (var < maxjeff))	{
					jeffkappa->GetVal(cont)->setval(var);
				}
			}
			else	{
				jeffkappa->GetVal(cont)->Sample();
			}
			for (int k=0; k<GetNconcept(); k++)	{
				if (!GetBeta(k,cont))	{
					beta->GetCell(k,cont)->Sample();
				}
			}
			jeffkappa->GetVal(cont)->Update();
		}
	}
};

class ConjugateAlphaVarMove : public MCUpdate	{

	NormalOntologyModel* model;

	public:

	ConjugateAlphaVarMove(NormalOntologyModel* inmodel)	{
		model = inmodel;
	}

	double Move(double tuning=1)	{
		model->ResampleAlphaVar();
		return 1;
	}
};

class ConjugateToggleBetaVarMove : public MCUpdate	{

	NormalOntologyModel* model;

	public:

	ConjugateToggleBetaVarMove(NormalOntologyModel* inmodel)	{
		model = inmodel;
	}

	double Move(double tuning=1)	{
		model->ResampleToggleBetaVar();
		return 1;
	}
};

class ConjugateBetaVarMove : public MCUpdate	{

	NormalOntologyModel* model;

	public:

	ConjugateBetaVarMove(NormalOntologyModel* inmodel)	{
		model = inmodel;
	}

	double Move(double tuning=1)	{
		model->ResampleBetaVar();
		return 1;
	}
};

class ConjugateAlphaMove : public MCUpdate, public Mnode	{

	private:

	NormalOntologyModel* model;

	public:

	ConjugateAlphaMove(NormalOntologyModel* inmodel)	{
		model = inmodel;
		model->GetAlpha()->RegisterArray(this);
	}

	double Move(double tuning = 1)	{

		Corrupt(false);
		model->ResampleAlpha();
		Update();
		return 1;
	}
};

class ConjugateBetaMove : public MCUpdate, public Mnode	{

	private:

	NormalOntologyModel* model;
	int cont;

	public:

	ConjugateBetaMove(NormalOntologyModel* inmodel, int incont)	{
		model = inmodel;
		cont = incont;
		model->GetBeta()->RegisterArray(this,cont);
	}

	double Move(double tuning = 1)	{

		Corrupt(false);
		model->ResampleBeta(cont);
		Update();
		return 1;
	}
};

class ConjugateToggleBetaMove : public MCUpdate, public Mnode	{

	private:

	NormalOntologyModel* model;
	int cont;
	int nrep;

	public:

	ConjugateToggleBetaMove(NormalOntologyModel* inmodel, int incont, int innrep)	{
		model = inmodel;
		cont = incont;
		nrep = innrep;
		model->GetBeta()->RegisterArray(this,cont);
		model->GetToggle()->RegisterArray(this,cont);
		model->GetTheta(cont)->Register(this);
	}

	double Move(double tuning = 1)	{

		Corrupt(false);
		double ret = model->ResampleToggleBeta(cont,nrep);
		Update();
		return ret;
	}
};

/*
class ConjugateAlphaBetaMove : public MCUpdate	{

	private:

	NormalOntologyModel* model;
	double tuning;
	int nrep;

	public:

	ConjugateAlphaBetaMove(NormalOntologyModel* inmodel, double intuning, int innrep)	{
		model = inmodel;
		tuning = intuning;
		nrep = innrep;
	}

	double Move(double tuning = 1)	{

		model->ResampleAlphaBeta(tuning,nrep);
		model->Update();
		return 1;
	}
};
*/

void NormalOntologyModel::MakeScheduler()	{

	bool conjugate = true;

	if (conjugate)	{
	// conjugate
	int nrep = 1;
	if (toggle)	{
		for (int rep=0; rep<nrep; rep++)	{
			scheduler.Register(new ConjugateAlphaMove(this),1,"conj alpha");
			for (int i=0; i<1; i++)	{
				scheduler.Register(new SimpleMove(toggle,1),1,"toggle");
				for (int cont=0; cont<GetNcont(); cont++)	{
					scheduler.Register(new ConjugateToggleBetaMove(this,cont,0),1,"conj beta");
					// scheduler.Register(new SimpleMove(theta[cont],0.1),10,"theta");
				}
			}
			/*
			for (int cont=0; cont<GetNcont(); cont++)	{
				scheduler.Register(new ConjugateToggleBetaMove(this,cont,1),1,"conj beta");
			}
			*/
			/*
			for (int cont=0; cont<GetNcont(); cont++)	{
				scheduler.Register(new SimpleMove(theta[cont],0.1),100,"theta");
			}
			*/
			/*
			scheduler.Register(new ConjugateAlphaVarMove(this),1,"conjugate alpha var");
			// scheduler.Register(new ConjugateToggleBetaVarMove(this),1,"conjugate beta var");
			scheduler.Register(new SimpleMove(jeffkappa,1),10,"kappa");
			scheduler.Register(new SimpleMove(jeffkappa,0.1),10,"kappa");
			*/
		}
	}
	else	{
		for (int rep=0; rep<nrep; rep++)	{
			scheduler.Register(new ConjugateAlphaMove(this),1,"conj alpha");
			for (int cont=0; cont<GetNcont(); cont++)	{
				scheduler.Register(new ConjugateBetaMove(this,cont),1,"conj beta");
			}
			scheduler.Register(new ConjugateAlphaVarMove(this),1,"conjugate alpha var");
			scheduler.Register(new ConjugateBetaVarMove(this),1,"conjugate beta var");
		}
	}

	}
	else	{
	// MH
	int nrep = 1;
	for (int rep=0; rep<nrep; rep++)	{
		scheduler.Register(new SimpleMove(alpha,1),1,"alpha");
		scheduler.Register(new SimpleMove(beta,1),1,"beta");
		scheduler.Register(new SimpleMove(jeffalphavar,1),100,"alpha var");
		scheduler.Register(new SimpleMove(jeffalphavar,0.1),100,"alpha var");
		scheduler.Register(new SimpleMove(jeffkappa,1),100,"kappa");
		scheduler.Register(new SimpleMove(jeffkappa,0.1),100,"kappa");
		if (toggle)	{
			scheduler.Register(new SimpleMove(toggle,1),1,"toggle");
			for (int cont=0; cont<GetNcont(); cont++)	{
				scheduler.Register(new SimpleMove(theta[cont],0.1),100,"theta");
			}
		}
	}
	}
}


#endif

