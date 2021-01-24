
#ifndef KALMANMULTIVARIATETREEPROCESS
#define KALMANMULTIVARIATETREEPROCESS

#include "ConjugateMultiVariateTreeProcess.h"

class KalmanMultiVariateTreeProcess : public virtual MultiVariateTreeProcess, public MatrixAlgebra {

	public:

	KalmanMultiVariateTreeProcess() : MultiVariateTreeProcess()	{
	}


	KalmanMultiVariateTreeProcess(Var<CovMatrix>* insigma, LengthTree* intree, LengthTree* inscaletree = 0, Var<RealVector>* indrift = 0, Var<RealVector>* inrootmean = 0, Var<PosRealVector>* inrootvar = 0) : MultiVariateTreeProcess(insigma,intree,inscaletree,indrift,inrootmean,inrootvar) {
	}

	void KalmanMove(int k, int n)	{
		// cerr << "create\n";
		KalmanCreate(k,n);
		// cerr << "backward\n";
		KalmanBackward(GetRoot(),k,n);
		// cerr << "forward\n";
		KalmanForward(GetRoot(),k,n);
		// cerr << "delete\n";
		KalmanDelete(k,n);
		// cerr << "quit\n";
		// cerr << '\n';
		// exit(1);
	}


	protected:

	void KalmanCreate(int k, int n)	{
		sigma->CorruptDiag();
		sigma->Diagonalise();
		// sigma->Invert();
		Omega = sigma->GetInvMatrix();
		Omegaxy = MatrixCreate(k,n);
		MatrixSet(Omegaxy,Omega,k,n,0,k);
		RecursiveKalmanCreate(GetRoot(),k,n);
	}

	void RecursiveKalmanCreate(const Link* from, int k, int n)	{
		KLambda[from] = new CovMatrix(k);
		KM[from] = new CovMatrix(k);
		KD[from] = new CovMatrix(k);
		KF[from] = MatrixCreate(k,1);
		Kgamma[from] = MatrixCreate(k,1);
		Kmu[from] = MatrixCreate(k,1);

		for (const Link* link=from->Next(); link!=from; link=link->Next())	{
			RecursiveKalmanCreate(link->Out(),k,n);
		}
	}

	void KalmanDelete(int k, int n)	{
		RecursiveKalmanDelete(GetRoot(),k,n);
		MatrixDelete(Omegaxy,k);
	}

	void RecursiveKalmanDelete(const Link* from, int k, int n)	{
		for (const Link* link=from->Next(); link!=from; link=link->Next())	{
			RecursiveKalmanDelete(link->Out(),k,n);
		}

		MatrixDelete(Kmu[from],k);
		MatrixDelete(Kgamma[from],k);
		MatrixDelete(KF[from],k);
		delete KD[from];
		delete KM[from];
		delete KLambda[from];
	}

	void KalmanBackward(const Link* from, int k, int n)	{

		CovMatrix* covM = KM[from];
		double** M = covM->GetMatrix();
		double** Omega = sigma->GetInvMatrix();

		if (from->isLeaf())	{
			bool cl = false;
			for (int i=0; i<k; i++)	{
				if (! cl)	{
					if (GetMultiNormal(from)->ClampVector[i])	{
						cl = true;
					}
				}
				else	{
					if (!GetMultiNormal(from)->ClampVector[i])	{
						cerr << "error in Kalman Backward: works only for all clamped or all unclamped leaf data\n";
						exit(1);
					}
				}
			}
			if (cl)	{
				double nu = 1.0 / GetBranchLength(from->GetBranch());
				MatrixSet(M,Omega,k,k);
				MatrixScalarProduct(M,nu,k,k);
				double** x = MatrixCreate(k,1);
				for (int i=0; i<k; i++)	{
					x[i][0] = (*GetMultiNormal(from))[i];
				}
				double** gamma = Kgamma[from];
				MatrixProduct(M,x,gamma,k,k,1);
				// cerr << "leaf M gamma : " << M[0][0] << '\t' << gamma[0][0] << '\t' << x[0][0] << '\n';
				double** E = MatrixCreate(k,n);
				double** y = MatrixCreate(n,1);
				for (int i=0; i<n; i++)	{
					y[i][0] = (*GetMultiNormal(from))[i+k] - (*GetMultiNormal(from->Out()))[i+k];
				}
				MatrixProduct(Omegaxy,y,E,k,n,1);
				MatrixScalarProduct(E,nu,k,1);
				MatrixAdd(gamma,E,k,1);
				MatrixDelete(x,k);
				MatrixDelete(y,n);
				MatrixDelete(E,k);
			}
			else	{
				MatrixSet(M,0,k,k);
				double** gamma = Kgamma[from];
				MatrixSet(gamma,0,k,1);
			}
		}
		else	{
			// here D is K in the manuscript
			// and gamma is M*gamma in the manuscript
			// compute D and mu from the Ms and gammas of children
			CovMatrix* covD = KD[from];
			double** D = covD->GetMatrix();
			MatrixSet(D,0,k,k);
			double** tmp = MatrixCreate(k,1);
			MatrixSet(tmp,0,k,1);

			for (const Link* link=from->Next(); link!=from; link=link->Next())	{

				KalmanBackward(link->Out(),k,n);

				MatrixAdd(D,KM[link->Out()]->GetMatrix(),k,k);
				MatrixAdd(tmp,Kgamma[link->Out()],k,1);
			}

			// covD->Invert();
			double** mu = Kmu[from];
			if (covD->GetMax() < 1e-8)	{
				MatrixSet(mu,0,k,1);
			}
			else	{
				covD->CorruptDiag();
				covD->Diagonalise();
				MatrixProduct(covD->GetInvMatrix(),tmp,mu,k,k,1);
			}
			// cerr << "internal D mu : " << D[0][0] << '\t' << mu[0][0] << '\n';
			MatrixDelete(tmp,k);

			// compute own M and gamma
			if (! from->isRoot())	{
			double nu = 1.0 / GetBranchLength(from->GetBranch());
			double** Lambda = KLambda[from]->GetMatrix();
			MatrixSet(Lambda,Omega,k,k);
			MatrixScalarProduct(Lambda,nu,k,k);
			MatrixAdd(Lambda,D,k,k);
			// KLambda[from]->Invert();
			KLambda[from]->CorruptDiag();
			KLambda[from]->Diagonalise();
			double** InvLambda = KLambda[from]->GetInvMatrix();

			double** E = MatrixCreate(k,n);
			double** y = MatrixCreate(n,1);
			for (int i=0; i<n; i++)	{
				// y[0][i] = (*GetMultiNormal(from))[i+k];
				y[i][0] = (*GetMultiNormal(from))[i+k] - (*GetMultiNormal(from->Out()))[i+k];
			}
			MatrixProduct(Omegaxy,y,E,k,n,1);
			MatrixScalarProduct(E,-nu,k,1);

			double** F = KF[from];
			MatrixProduct(D,mu,F,k,k,1);
			MatrixAdd(F,E,k,1);

			double** G = MatrixCreate(k,k);
			double** H = MatrixCreate(k,k);
			MatrixProduct(Omega,InvLambda,G,k,k,k);
			MatrixScalarProduct(G,-nu,k,k);
			MatrixSetIdentity(H,k);
			MatrixAdd(H,G,k,k);
			MatrixScalarProduct(G,-1,k,k);
			MatrixScalarProduct(H,nu,k,k);

			double** M = KM[from]->GetMatrix();
			double** gamma = Kgamma[from];
			MatrixProduct(H,Omega,M,k,k,k);
			MatrixProduct(G,F,gamma,k,k,1);
			MatrixScalarProduct(E,-1,k,1);
			MatrixAdd(gamma,E,k,1);

			MatrixDelete(E,k);
			MatrixDelete(y,n);
			MatrixDelete(G,k);
			MatrixDelete(H,k);
			// cerr << "internal M gamma : " << M[0][0] << '\t' << gamma[0][0] << '\n';
			}
		}
	}

	void KalmanForward(const Link* from, int k, int n, const Link* parent = 0)	{
		if (from->isRoot())	{
			// draw N(mu,D-1)
			KD[from]->CorruptDiag();
			KD[from]->Diagonalise();
			GetMultiNormal(from)->DrawNormalFromPrecision(Kmu[from],KD[from]);
		}
		else	{

			if (from->isLeaf())	{
				bool cl = false;
				for (int i=0; i<k; i++)	{
					if (! cl)	{
						if (GetMultiNormal(from)->ClampVector[i])	{
							cl = true;
						}
					}
					else	{
						if (!GetMultiNormal(from)->ClampVector[i])	{
							cerr << "error in Kalman Forward: works only for all clamped or all unclamped leaf data\n";
							exit(1);
						}
					}
				}
				if (!cl)	{
					double** alpha = MatrixCreate(k,1);
					double** tmp = MatrixCreate(k,n);
					double** y = MatrixCreate(n,1);

					double nu = 1.0 / GetBranchLength(from->GetBranch());

					double** Lambda = KLambda[from]->GetMatrix();
					MatrixSet(Lambda,Omega,k,k);
					MatrixScalarProduct(Lambda,nu,k,k);
					KLambda[from]->CorruptDiag();
					KLambda[from]->Diagonalise();
					double** InvLambda = KLambda[from]->GetInvMatrix();

					for (int i=0; i<n; i++)	{
						// y[0][i] = (*GetMultiNormal(from))[i+k];
						y[i][0] = (*GetMultiNormal(from))[i+k] - (*GetMultiNormal(from->Out()))[i+k];
					}

					MatrixProduct(Omegaxy,y,tmp,k,n,1);
					MatrixProduct(InvLambda,tmp,alpha,k,k,1);
					MatrixScalarProduct(alpha,-nu,k,1);
					for (int i=0; i<k; i++)	{
						alpha[i][0] += (*GetMultiNormal(parent))[i];
					}
					GetMultiNormal(from)->DrawNormalFromPrecision(alpha,KLambda[from]);

					MatrixDelete(alpha,k);
					MatrixDelete(tmp,k);
					MatrixDelete(y,k);
				}
			}
			else	{

				// grep value of parent x
				double** xup = MatrixCreate(k,1);
				for (int i=0; i<k; i++)	{
					xup[i][0] = (*GetMultiNormal(parent))[i];
				}

				double nu = 1.0 / GetBranchLength(from->GetBranch());

				double** tmp = MatrixCreate(k,1);
				double** alpha = MatrixCreate(k,1);

				MatrixProduct(Omega,xup,tmp,k,k,1);
				MatrixScalarProduct(tmp,nu,k,1);
				MatrixAdd(tmp,KF[from],k,1);
				KLambda[from]->CorruptDiag();
				KLambda[from]->Diagonalise();
				double** InvLambda = KLambda[from]->GetInvMatrix();
				MatrixProduct(InvLambda,tmp,alpha,k,k,1);

				// draw N(alpha,Lambda-1)
				// cerr << "from lambda : " << KLambda[from]->GetMatrix()[0][0] << '\n';
				GetMultiNormal(from)->DrawNormalFromPrecision(alpha,KLambda[from]);

				MatrixDelete(tmp,k);
				MatrixDelete(alpha,k);
				MatrixDelete(xup,k);

			}
		}

		// recursive call
		for (const Link* link=from->Next(); link!=from; link=link->Next())	{
			KalmanForward(link->Out(),k,n,from);
		}
	}

	double** Omega;
	double** Omegaxy;

	map<const Link*, CovMatrix*> KLambda;
	map<const Link*, double**> KF;
	map<const Link*, CovMatrix*> KM;
	map<const Link*, double**> Kgamma;
	map<const Link*, CovMatrix*> KD;
	map<const Link*, double**> Kmu;

};

class ExternalKalmanMultiVariateTreeProcess : public virtual MultiVariateTreeProcess, public MatrixAlgebra {

	public:

	ExternalKalmanMultiVariateTreeProcess() : MultiVariateTreeProcess()	{
	}


	ExternalKalmanMultiVariateTreeProcess(Var<CovMatrix>* insigma, LengthTree* intree, LengthTree* inscaletree = 0, Var<RealVector>* indrift = 0, Var<RealVector>* inrootmean = 0, Var<PosRealVector>* inrootvar = 0) : MultiVariateTreeProcess(insigma,intree,inscaletree,indrift,inrootmean,inrootvar) {
	}

	void KalmanMove()	{
		// cerr << "create\n";
		// KalmanCreate();
		cerr << "backward\n";
		KalmanBackward(GetRoot());
		cerr << "forward\n";
		KalmanForward(GetRoot());
		// cerr << "delete\n";
		// KalmanDelete();
		// cerr << "quit\n";
		// cerr << '\n';
		// exit(1);
		cerr << "AFTER KALMAN\n";
		for (int i=0; i<GetDim(); i++)	{
			cerr << GetMean(i) << '\t';
			if (std::isinf(GetMean(i)) || std::isnan(GetMean(i)))	{
				cerr << "error : inf or nan\n";
				cerr << (*this) << '\n';
				exit(1);
			}
		}
		cerr << '\n';
	}

	double** GetW()	{
		return W;
	}

	double** GetKy(const Link* from)	{
		return Ky[from];
	}

	void KalmanCreate()	{

		k = GetDim();
		CovOmegaW = new CovMatrix(k);
		OmegaW = CovOmegaW->GetMatrix();
		// OmegaW = MatrixCreate(k,k);
		W = MatrixCreate(k,k);

		RecursiveKalmanCreate(GetRoot());
	}

	void KalmanDelete()	{
		RecursiveKalmanDelete(GetRoot());
		delete CovOmegaW;
		// MatrixDelete(OmegaW,k);
		MatrixDelete(W,k);
	}

	void KalmanInit()	{
		sigma->CorruptDiag();
		sigma->Diagonalise();
		Omega = sigma->GetInvMatrix();
		MatrixSet(OmegaW,Omega,k,k);
		MatrixAdd(OmegaW,W,k,k);
		CovOmegaW->CorruptDiag();
		CovOmegaW->Diagonalise();
	}

	protected:

	void RecursiveKalmanCreate(const Link* from)	{
		KLambda[from] = new CovMatrix(k);
		KM[from] = new CovMatrix(k);
		KD[from] = new CovMatrix(k);
		KF[from] = MatrixCreate(k,1);
		Ky[from] = MatrixCreate(k,1);
		Kgamma[from] = MatrixCreate(k,1);
		Kmu[from] = MatrixCreate(k,1);

		for (const Link* link=from->Next(); link!=from; link=link->Next())	{
			RecursiveKalmanCreate(link->Out());
		}

	}

	void RecursiveKalmanDelete(const Link* from)	{
		for (const Link* link=from->Next(); link!=from; link=link->Next())	{
			RecursiveKalmanDelete(link->Out());
		}

		MatrixDelete(Kmu[from],k);
		MatrixDelete(Kgamma[from],k);
		MatrixDelete(KF[from],k);
		MatrixDelete(Ky[from],k);
		delete KD[from];
		delete KM[from];
		delete KLambda[from];
	}

	void KalmanBackward(const Link* from)	{

		CovMatrix* covM = KM[from];
		double** M = covM->GetMatrix();

		if (from->isLeaf())	{
			int clamp = 0;
			for (int i=0; i<k; i++)	{
				if (GetMultiNormal(from)->ClampVector[i])	{
					clamp++;
				}
			}
			if (clamp == k)	{
				cerr << "clamp k\n";
				exit(1);
				double nu = 1.0 / GetBranchLength(from->GetBranch());
				MatrixSet(M,OmegaW,k,k);
				MatrixScalarProduct(M,nu,k,k);

				double** x = MatrixCreate(k,1);
				for (int i=0; i<k; i++)	{
					x[i][0] = (*GetMultiNormal(from))[i];
				}
				double** gamma = Kgamma[from];
				MatrixProduct(M,x,gamma,k,k,1);
				double** E = MatrixCreate(k,1);
				double** y = Ky[from];
				MatrixSet(E,y,k,1);
				MatrixScalarProduct(E,-nu,k,1);
				MatrixAdd(gamma,E,k,1);
				MatrixDelete(E,k);
				MatrixDelete(x,k);
			}
			else if (clamp == 0)	{
				cerr << "clamp 0\n";
				exit(1);
				MatrixSet(M,0,k,k);
				double** gamma = Kgamma[from];
				MatrixSet(gamma,0,k,1);
			}
			else	{
				int K = 0;
				for (int i=0; i<GetDim(); i++)	{
					if (! GetMultiNormal(from)->ClampVector[i])	{
						K++;
					}
				}
				int L = GetDim() - K;
				// not clamped
				int* index1 = new int[K];
				// clamped
				int* index2 = new int[L];
				int kk = 0;
				int ll = 0;

				for (int i=0; i<GetDim(); i++)	{
					if (GetMultiNormal(from)->ClampVector[i])	{
						index2[ll] = i;
						ll++;
					}
					else	{
						index1[kk] = i;
						kk++;
					}
				}

				double nu = 1.0 / GetBranchLength(from->GetBranch());
				MatrixSet(M,OmegaW,k,k);
				MatrixScalarProduct(M,nu,k,k);

				double** Lambda22 = MatrixCreate(L,L);
				double** Lambda12 = MatrixCreate(K,L);
				double** Lambda21 = MatrixCreate(L,K);
				CovMatrix* CovLambda11 = new CovMatrix(K);
				double** Lambda11 = CovLambda11->GetMatrix();

				for (int i=0; i<K; i++)	{
					for (int j=0; j<K; j++)	{
						Lambda11[i][j] = M[index1[i]][index1[j]];
					}
				}

				for (int i=0; i<K; i++)	{
					for (int j=0; j<L; j++)	{
						Lambda12[i][j] = M[index1[i]][index2[j]];
					}
				}

				for (int i=0; i<L; i++)	{
					for (int j=0; j<K; j++)	{
						Lambda21[i][j] = M[index2[i]][index1[j]];
					}
				}

				for (int i=0; i<L; i++)	{
					for (int j=0; j<L; j++)	{
						Lambda22[i][j] = M[index2[i]][index2[j]];
					}
				}

				CovLambda11->CorruptDiag();
				CovLambda11->Diagonalise();
				double** InvLambda11 = CovLambda11->GetInvMatrix();

				double** tmp12 = MatrixCreate(K,L);
				double** tmp22 = MatrixCreate(L,L);
				MatrixProduct(InvLambda11,Lambda12,tmp12,K,K,L);
				MatrixProduct(Lambda21,tmp12,tmp22,L,K,L);
				MatrixScalarProduct(tmp22,-1,L,L);
				MatrixAdd(Lambda22,tmp22,L,L);

				for (int i=0; i<k; i++)	{
					for (int j=0; j<k; j++)	{
						M[i][j] = 0;
					}
				}
				for (int i=0; i<L; i++)	{
					for (int j=0; j<L; j++)	{
						M[index2[i]][index2[j]] = Lambda22[i][j];
					}
				}

				double** x = MatrixCreate(k,1);
				for (int i=0; i<k; i++)	{
					x[i][0] = (*GetMultiNormal(from))[i];
				}
				for (int i=0; i<K; i++)	{
					x[index1[i]][0] = 0;
				}

				double** InvOmegaW = CovOmegaW->GetInvMatrix();
				double** E = MatrixCreate(k,1);
				double** y = Ky[from];
				MatrixProduct(InvOmegaW,y,E,k,k,1);
				MatrixScalarProduct(E,-1,k,1);
				MatrixAdd(x,E,k,1);

				double** gamma = Kgamma[from];
				MatrixProduct(M,E,gamma,k,k,1);

				MatrixDelete(E,k);
				MatrixDelete(x,k);
				MatrixDelete(Lambda12,K);
				MatrixDelete(Lambda21,L);
				MatrixDelete(Lambda22,L);
				delete CovLambda11;
				MatrixDelete(tmp12,K);
				MatrixDelete(tmp22,L);
			}
		}
		else	{
			// here D is K in the manuscript
			// and gamma is M*beta in the manuscript
			// compute D and mu from the Ms and gammas of children
			CovMatrix* covD = KD[from];
			double** D = covD->GetMatrix();
			MatrixSet(D,0,k,k);
			double** tmp = MatrixCreate(k,1);
			MatrixSet(tmp,0,k,1);

			for (const Link* link=from->Next(); link!=from; link=link->Next())	{

				KalmanBackward(link->Out());

				MatrixAdd(D,KM[link->Out()]->GetMatrix(),k,k);
				MatrixAdd(tmp,Kgamma[link->Out()],k,1);
			}

			// covD->Invert();
			double** mu = Kmu[from];
			if (covD->GetMax() < 1e-8)	{
				MatrixSet(mu,0,k,1);
			}
			else	{
				covD->CorruptDiag();
				covD->Diagonalise();
				MatrixProduct(covD->GetInvMatrix(),tmp,mu,k,k,1);
			}
			// cerr << "internal D mu : " << D[0][0] << '\t' << mu[0][0] << '\n';
			MatrixDelete(tmp,k);

			// now compute own M and gamma
			if (! from->isRoot())	{

			// inverse branch length
			double nu = 1.0 / GetBranchLength(from->GetBranch());

			// compute Lambda and its inverse
			double** Lambda = KLambda[from]->GetMatrix();
			MatrixSet(Lambda,OmegaW,k,k);
			MatrixScalarProduct(Lambda,nu,k,k);
			MatrixAdd(Lambda,D,k,k);
			// KLambda[from]->Invert();
			KLambda[from]->CorruptDiag();
			KLambda[from]->Diagonalise();
			double** InvLambda = KLambda[from]->GetInvMatrix();

			// M <- v (Omega + W) - v^2 (Omega + W) Lambda^-1 (Omega + W)
			double** G = MatrixCreate(k,k);
			double** H = MatrixCreate(k,k);
			MatrixProduct(OmegaW,InvLambda,G,k,k,k);
			MatrixScalarProduct(G,-nu,k,k);
			MatrixSetIdentity(H,k);
			MatrixAdd(H,G,k,k);
			MatrixScalarProduct(H,nu,k,k);
			double** M = KM[from]->GetMatrix();
			MatrixProduct(H,OmegaW,M,k,k,k);
			MatrixScalarProduct(G,-1,k,k);
			// now, G = nu (Omega + W) Lambda^-1

			double** gamma = Kgamma[from];

			double** E = MatrixCreate(k,1);
			double** y = Ky[from];
			MatrixSet(E,y,k,1);
			MatrixScalarProduct(E,nu,k,1);

			double** F = KF[from];
			MatrixProduct(D,mu,F,k,k,1);
			MatrixAdd(F,E,k,1);

			MatrixProduct(G,F,gamma,k,k,1);
			MatrixScalarProduct(E,-1,k,1);
			MatrixAdd(gamma,E,k,1);

			MatrixDelete(E,k);
			// MatrixDelete(y,n);
			MatrixDelete(G,k);
			MatrixDelete(H,k);
			// cerr << "internal M gamma : " << M[0][0] << '\t' << gamma[0][0] << '\n';
			}
		}
	}

	void KalmanForward(const Link* from, const Link* parent = 0)	{
		if (from->isRoot())	{

			int clamp = 0;
			for (int i=0; i<k; i++)	{
				if (GetMultiNormal(from)->ClampVector[i])	{
					clamp++;
				}
			}
			if (clamp < k)	{
			// draw N(mu,D-1)
			KD[from]->CorruptDiag();
			KD[from]->Diagonalise();
			double** mu = Kmu[from];
			for (int i=0; i<k; i++)	{
				if (GetMultiNormal(from)->ClampVector[i])	{
					mu[i][0] = (*GetMultiNormal(from))[i];
				}
			}
			GetMultiNormal(from)->DrawNormalFromPrecisionWithConstraint(Kmu[from],KD[from]);
			}
		}
		else	{

			if (from->isLeaf())	{
				int clamp = 0;
				for (int i=0; i<k; i++)	{
					if (GetMultiNormal(from)->ClampVector[i])	{
						clamp++;
					}
				}
				if (clamp < k)	{
					double** alpha = MatrixCreate(k,1);

					double nu = 1.0 / GetBranchLength(from->GetBranch());

					double** Lambda = KLambda[from]->GetMatrix();
					MatrixSet(Lambda,OmegaW,k,k);
					MatrixScalarProduct(Lambda,nu,k,k);
					KLambda[from]->CorruptDiag();
					KLambda[from]->Diagonalise();
					double** InvLambda = KLambda[from]->GetInvMatrix();

					double** y = Ky[from];
					MatrixProduct(InvLambda,y,alpha,k,k,1);
					MatrixScalarProduct(alpha,nu,k,1);
					for (int i=0; i<k; i++)	{
						alpha[i][0] += (*GetMultiNormal(parent))[i];
					}
					if (! clamp)	{
						GetMultiNormal(from)->DrawNormalFromPrecision(alpha,KLambda[from]);
					}
					else	{
						GetMultiNormal(from)->DrawNormalFromPrecisionWithConstraint(alpha,KLambda[from]);
					}

					MatrixDelete(alpha,k);
				}
			}
			else	{

				// grep value of parent x
				double** xup = MatrixCreate(k,1);
				for (int i=0; i<k; i++)	{
					xup[i][0] = (*GetMultiNormal(parent))[i];
				}

				double nu = 1.0 / GetBranchLength(from->GetBranch());

				double** tmp = MatrixCreate(k,1);
				double** alpha = MatrixCreate(k,1);

				MatrixProduct(OmegaW,xup,tmp,k,k,1);
				MatrixScalarProduct(tmp,nu,k,1);
				MatrixAdd(tmp,KF[from],k,1);

				KLambda[from]->CorruptDiag();
				KLambda[from]->Diagonalise();
				double** InvLambda = KLambda[from]->GetInvMatrix();
				MatrixProduct(InvLambda,tmp,alpha,k,k,1);

				// draw N(alpha,Lambda-1)
				// cerr << "from lambda : " << KLambda[from]->GetMatrix()[0][0] << '\n';
				GetMultiNormal(from)->DrawNormalFromPrecision(alpha,KLambda[from]);

				MatrixDelete(tmp,k);
				MatrixDelete(alpha,k);
				MatrixDelete(xup,k);

			}
		}

		// recursive call
		for (const Link* link=from->Next(); link!=from; link=link->Next())	{
			KalmanForward(link->Out(),from);
		}
	}

	int k;
	double** Omega;
	double** OmegaW;
	CovMatrix* CovOmegaW;
	double** W;
	map<const Link*, CovMatrix*> KLambda;
	map<const Link*, double**> KF;
	map<const Link*, double**> Ky;
	map<const Link*, CovMatrix*> KM;
	map<const Link*, double**> Kgamma;
	map<const Link*, CovMatrix*> KD;
	map<const Link*, double**> Kmu;

};

class ConditionalExternalKalmanMultiVariateTreeProcess : public virtual MultiVariateTreeProcess, public MatrixAlgebra {

	public:

	ConditionalExternalKalmanMultiVariateTreeProcess() : MultiVariateTreeProcess()	{
	}


	ConditionalExternalKalmanMultiVariateTreeProcess(Var<CovMatrix>* insigma, LengthTree* intree, LengthTree* inscaletree = 0, Var<RealVector>* indrift = 0, Var<RealVector>* inrootmean = 0, Var<PosRealVector>* inrootvar = 0) : MultiVariateTreeProcess(insigma,intree,inscaletree,indrift,inrootmean,inrootvar) {
	}

	void KalmanMove()	{
		// cerr << "create\n";
		// KalmanCreate();
		KalmanBackward(GetRoot());
		KalmanForward(GetRoot());
		// cerr << "delete\n";
		// KalmanDelete();
		// cerr << "quit\n";
		// cerr << '\n';
		// exit(1);
		/*
		cerr << "AFTER KALMAN\n";
		for (int i=0; i<GetDim(); i++)	{
			cerr << GetMean(i) << '\t';
			if (isinf(GetMean(i)) || isnan(GetMean(i)))	{
				cerr << "error : inf or nan\n";
				cerr << (*this) << '\n';
				exit(1);
			}
		}
		cerr << '\n';
		*/
	}

	double** GetW()	{
		return W;
	}

	double** GetKa(const Link* from)	{
		return Ka[from];
	}

	void KalmanCreate(int* profile)	{

		index = new int[GetDim()];
		k = 0;
		n = 0;
		for (int i=0; i<GetDim(); i++)	{
			if (profile[i] == 1)	{
				k++;
				index[i] = 1;
			}
			else	{
				n++;
				index[i] = 0;
			}
		}
		// not clamped
		index1 = new int[k];
		// clamped
		index2 = new int[n];
		int ii = 0;
		int jj = 0;
		for (int i=0; i<GetDim(); i++)	{
			if (index[i] == 1)	{
				index1[ii] = i;
				ii++;
			}
			else	{
				index2[jj] = i;
				jj++;
			}
		}

		CovQ = new CovMatrix(GetDim());
		Qyy = MatrixCreate(k,k);
		Qyz = MatrixCreate(k,n);
		W = MatrixCreate(GetDim(),GetDim());

		RecursiveKalmanCreate(GetRoot());
	}

	void KalmanDelete()	{
		RecursiveKalmanDelete(GetRoot());
		MatrixDelete(W,GetDim());
		MatrixDelete(Qyz,k);
		MatrixDelete(Qyy,k);
		delete CovQ;
		delete[] index;
		delete[] index1;
		delete[] index2;
	}

	void KalmanInit()	{
		sigma->CorruptDiag();
		sigma->Diagonalise();
		Omega = sigma->GetInvMatrix();
		Q = CovQ->GetMatrix();
		MatrixSet(Q,Omega,GetDim(),GetDim());
		MatrixAdd(Q,W,GetDim(),GetDim());
		CovQ->CorruptDiag();
		CovQ->Diagonalise();
		InvQ = CovQ->GetInvMatrix();
		for (int i=0; i<k; i++)	{
			for (int j=0; j<k; j++)	{
				Qyy[i][j] = Q[index1[i]][index1[j]];
			}
		}
		for (int i=0; i<k; i++)	{
			for (int j=0; j<n; j++)	{
				Qyz[i][j] = Q[index1[i]][index2[j]];
			}
		}
	}

	protected:

	void RecursiveKalmanCreate(const Link* from)	{
		KLambda[from] = new CovMatrix(k);
		Ka[from] = MatrixCreate(GetDim(),1);
		KM[from] = new CovMatrix(k);
		KD[from] = new CovMatrix(k);
		KF[from] = MatrixCreate(k,1);
		Kgamma[from] = MatrixCreate(k,1);
		Kmu[from] = MatrixCreate(k,1);

		for (const Link* link=from->Next(); link!=from; link=link->Next())	{
			RecursiveKalmanCreate(link->Out());
		}

	}

	void RecursiveKalmanDelete(const Link* from)	{
		for (const Link* link=from->Next(); link!=from; link=link->Next())	{
			RecursiveKalmanDelete(link->Out());
		}

		MatrixDelete(Kmu[from],k);
		MatrixDelete(Kgamma[from],k);
		MatrixDelete(KF[from],k);
		MatrixDelete(Ka[from],k);
		delete KD[from];
		delete KM[from];
		delete KLambda[from];
	}

	void KalmanBackward(const Link* from)	{

		// inverse branch length
		double nu = 1.0 / GetBranchLength(from->GetBranch());

		// compute a, b and split b into y and z part
		double** a = Ka[from];
		double** b = MatrixCreate(GetDim(),1);
		MatrixProduct(InvQ,a,b,GetDim(),GetDim(),1);
		double** by = MatrixCreate(k,1);
		double** bz = MatrixCreate(n,1);
		for (int i=0; i<k; i++)	{
			by[i][0] = b[index1[i]][0];
		}
		for (int i=0; i<n; i++)	{
			bz[i][0] = b[index2[i]][0];
		}

		// compute c = vu[Qyy by + Qyz (Deltaz - bz)]
		double** c = MatrixCreate(k,1);
		MatrixScalarProduct(bz,-1,n,1);
		for (int i=0; i<n; i++)	{
			bz[i][0] += (*GetMultiNormal(from))[index2[i]] - (*GetMultiNormal(from->Out()))[index2[i]];
		}
		MatrixProduct(Qyz,bz,c,k,n,1);
		double** Qyyby = MatrixCreate(k,1);
		MatrixProduct(Qyy,by,Qyyby,k,k,1);
		MatrixAdd(c,Qyyby,k,1);
		MatrixScalarProduct(c,nu,k,1);

		// take pointers
		CovMatrix* covM = KM[from];
		double** M = covM->GetMatrix();
		double** gamma = Kgamma[from];

		if (from->isLeaf())	{

			// check clamp status
			int clamp = 0;
			for (int i=0; i<k; i++)	{
				if (GetMultiNormal(from)->ClampVector[index1[i]])	{
					clamp++;
				}
			}
			if ((clamp > 0) && (clamp < k))	{
				cerr << "error in kalman backward: clamping profile\n";
				exit(1);
			}

			if (clamp == k)	{

				// F and Lambda : don't care

				MatrixSet(M,Qyy,k,k);
				MatrixScalarProduct(M,nu,k,k);

				double** y = MatrixCreate(k,1);
				for (int i=0; i<k; i++)	{
					y[i][0] = (*GetMultiNormal(from))[index1[i]];
				}
				MatrixProduct(M,y,gamma,k,k,1);
				MatrixScalarProduct(c,-1,k,1);
				MatrixAdd(gamma,c,k,1);
				MatrixDelete(y,k);

			}
			else {

				// compute F = c
				double** F = KF[from];
				MatrixSet(F,c,k,1);

				// compute Lambda and its inverse
				double** Lambda = KLambda[from]->GetMatrix();
				MatrixSet(Lambda,Qyy,k,k);
				MatrixScalarProduct(Lambda,nu,k,k);
				KLambda[from]->CorruptDiag();
				KLambda[from]->Diagonalise();

				// set M = 0 and gamma = 0
				MatrixSet(M,0,k,k);
				MatrixSet(gamma,0,k,1);
			}
		}
		else	{

			// here D is K in the manuscript
			// compute D and mu from the Ms and gammas of children
			CovMatrix* covD = KD[from];
			double** D = covD->GetMatrix();
			MatrixSet(D,0,k,k);
			double** tmp = MatrixCreate(k,1);
			MatrixSet(tmp,0,k,1);

			for (const Link* link=from->Next(); link!=from; link=link->Next())	{

				KalmanBackward(link->Out());

				MatrixAdd(D,KM[link->Out()]->GetMatrix(),k,k);
				MatrixAdd(tmp,Kgamma[link->Out()],k,1);
			}

			double** mu = Kmu[from];
			if (covD->GetMax() < 1e-8)	{
				MatrixSet(mu,0,k,1);
			}
			else	{
				covD->CorruptDiag();
				covD->Diagonalise();
				MatrixProduct(covD->GetInvMatrix(),tmp,mu,k,k,1);
			}
			MatrixDelete(tmp,k);

			// now compute own M and gamma
			if (! from->isRoot())	{

			// compute F = D mu + c
			double** F = KF[from];
			MatrixProduct(D,mu,F,k,k,1);
			MatrixAdd(F,c,k,1);

			// compute Lambda and its inverse
			double** Lambda = KLambda[from]->GetMatrix();
			MatrixSet(Lambda,Qyy,k,k);
			MatrixScalarProduct(Lambda,nu,k,k);
			MatrixAdd(Lambda,D,k,k);
			// KLambda[from]->Invert();
			KLambda[from]->CorruptDiag();
			KLambda[from]->Diagonalise();

			double** InvLambda = KLambda[from]->GetInvMatrix();

			// M = v [I - v Qyy Lambda^-1] Qyy
			double** G = MatrixCreate(k,k);
			double** H = MatrixCreate(k,k);
			MatrixProduct(Qyy,InvLambda,G,k,k,k);
			MatrixScalarProduct(G,-nu,k,k);
			MatrixSetIdentity(H,k);
			MatrixAdd(H,G,k,k);
			MatrixScalarProduct(H,nu,k,k);
			MatrixProduct(H,Qyy,M,k,k,k);
			MatrixDelete(H,k);

			// gamma = v Qyy F - c
			double** gamma = Kgamma[from];
			MatrixScalarProduct(G,-1,k,k);
			MatrixProduct(G,F,gamma,k,k,1);
			MatrixScalarProduct(c,-1,k,1);
			MatrixAdd(gamma,c,k,1);
			MatrixDelete(G,k);

			}
		}

		MatrixDelete(c,k);
		MatrixDelete(Qyyby,k);
		MatrixDelete(bz,n);
		MatrixDelete(by,k);
		MatrixDelete(b,GetDim());
	}

	void KalmanForward(const Link* from, const Link* parent = 0)	{
		if (from->isRoot())	{

			int clamp = 0;
			for (int i=0; i<k; i++)	{
				if (GetMultiNormal(from)->ClampVector[index1[i]])	{
					clamp++;
				}
			}
			if ((clamp > 0) && (clamp < k))	{
				cerr << "error in kalman forward: root clamping profile\n";
				exit(1);
			}
			if (!clamp)	{
				// draw N(mu,D-1)
				KD[from]->CorruptDiag();
				KD[from]->Diagonalise();
				double** mu = Kmu[from];
				GetMultiNormal(from)->ConditionalDrawNormalFromPrecision(Kmu[from],KD[from],index1,k);
			}
		}
		else	{

			int clamp = 0;
			if (from->isLeaf())	{
				for (int i=0; i<k; i++)	{
					if (GetMultiNormal(from)->ClampVector[index1[i]])	{
						clamp++;
					}
				}
				if ((clamp > 0) && (clamp < k))	{
					cerr << "error in kalman forward: root clamping profile\n";
					exit(1);
				}
			}

			if ((! from->isLeaf()) || (! clamp))	{

				// grep value of parent x
				double** yup = MatrixCreate(k,1);
				for (int i=0; i<k; i++)	{
					yup[i][0] = (*GetMultiNormal(parent))[index1[i]];
				}

				double nu = 1.0 / GetBranchLength(from->GetBranch());

				// calculate alpha
				double** alpha = MatrixCreate(k,1);
				double** tmp = MatrixCreate(k,1);

				MatrixProduct(Qyy,yup,tmp,k,k,1);
				MatrixScalarProduct(tmp,nu,k,1);
				MatrixAdd(tmp,KF[from],k,1);

				double** InvLambda = KLambda[from]->GetInvMatrix();
				MatrixProduct(InvLambda,tmp,alpha,k,k,1);

				// draw y ~ N(alpha,Lambda^-1)
				GetMultiNormal(from)->ConditionalDrawNormalFromPrecision(alpha,KLambda[from],index1,k);

				MatrixDelete(tmp,k);
				MatrixDelete(alpha,k);
				MatrixDelete(yup,k);

			}
		}

		// recursive call
		for (const Link* link=from->Next(); link!=from; link=link->Next())	{
			KalmanForward(link->Out(),from);
		}
	}

	int k;
	int n;
	int* index;
	int* index1;
	int* index2;

	double** Omega;
	CovMatrix* CovQ;
	double** Q;
	double** InvQ;
	double** Qyy;
	double** Qyz;
	double** W;

	map<const Link*, CovMatrix*> KLambda;
	map<const Link*, double**> KF;
	map<const Link*, double**> Ka;
	map<const Link*, CovMatrix*> KM;
	map<const Link*, double**> Kgamma;
	map<const Link*, CovMatrix*> KD;
	map<const Link*, double**> Kmu;

};

/*
class ConditionalExternalKalmanMultiVariateTreeProcess : public virtual MultiVariateTreeProcess, public MatrixAlgebra {

	public:

	ConditionalExternalKalmanMultiVariateTreeProcess() : MultiVariateTreeProcess()	{
	}


	ConditionalExternalKalmanMultiVariateTreeProcess(Var<CovMatrix>* insigma, LengthTree* intree, LengthTree* inscaletree = 0, Var<RealVector>* indrift = 0, Var<RealVector>* inrootmean = 0, Var<PosRealVector>* inrootvar = 0) : MultiVariateTreeProcess(insigma,intree,inscaletree,indrift,inrootmean,inrootvar) {
	}

	void KalmanMove()	{
		// cerr << "create\n";
		// KalmanCreate();
		cerr << "backward\n";
		KalmanBackward(GetRoot());
		cerr << "forward\n";
		KalmanForward(GetRoot());
		// cerr << "delete\n";
		// KalmanDelete();
		// cerr << "quit\n";
		// cerr << '\n';
		// exit(1);
		cerr << "AFTER KALMAN\n";
		for (int i=0; i<GetDim(); i++)	{
			cerr << GetMean(i) << '\t';
			if (isinf(GetMean(i)) || isnan(GetMean(i)))	{
				cerr << "error : inf or nan\n";
				cerr << (*this) << '\n';
				exit(1);
			}
		}
		cerr << '\n';
	}

	double** GetW()	{
		return WW;
	}

	double** GetKy(const Link* from)	{
		return Ky[from];
	}

	void KalmanCreate(int* profile)	{

		cerr << "kalman create\n";
		index = new int[GetDim()];
		k = 0;
		n = 0;
		for (int i=0; i<GetDim(); i++)	{
			cerr << profile[i] << '\t';
			if (profile[i] == 1)	{
				k++;
				index[i] = 1;
			}
			else	{
				n++;
				index[i] = 0;
			}
		}
		cerr << '\n';
		cerr << k << '\t' << n << '\n';
		// not clamped
		index1 = new int[k];
		// clamped
		index2 = new int[n];
		int ii = 0;
		int jj = 0;
		for (int i=0; i<GetDim(); i++)	{
			if (index[i] == 1)	{
				index1[ii] = i;
				ii++;
			}
			else	{
				index2[jj] = i;
				jj++;
			}
		}
		Omegaxy = MatrixCreate(k,n);
		WW = MatrixCreate(GetDim(),GetDim());
		W = MatrixCreate(k,k);
		CovOmegaW = new CovMatrix(k);

		RecursiveKalmanCreate(GetRoot());
		cerr << "create ok\n";
	}

	void KalmanDelete()	{
		RecursiveKalmanDelete(GetRoot());
		delete CovOmegaW;
		// MatrixDelete(OmegaW,k);
		MatrixDelete(W,k);
		delete[] index;
		delete[] index1;
		delete[] index2;
		MatrixDelete(Omegaxy,k);
	}

	void KalmanInit()	{
		cerr << "init\n";
		sigma->CorruptDiag();
		sigma->Diagonalise();
		Omega = sigma->GetInvMatrix();
		for (int i=0; i<k; i++)	{
			for (int j=0; j<n; j++)	{
				Omegaxy[i][j] = Omega[index1[i]][index2[j]];
			}
		}

		OmegaW = CovOmegaW->GetMatrix();
		for (int i=0; i<k; i++)	{
			for (int j=0; j<k; j++)	{
				W[i][j] = WW[index1[i]][index1[j]];
				OmegaW[i][j] += Omega[index1[i]][index1[j]] + WW[index1[i]][index1[j]];
			}
		}
		CovOmegaW->CorruptDiag();
		CovOmegaW->Diagonalise();
		cerr << "init ok\n";
	}

	protected:

	void RecursiveKalmanCreate(const Link* from)	{
		KLambda[from] = new CovMatrix(k);
		KM[from] = new CovMatrix(k);
		KD[from] = new CovMatrix(k);
		KF[from] = MatrixCreate(k,1);
		Ky[from] = MatrixCreate(GetDim(),1);
		Kgamma[from] = MatrixCreate(k,1);
		Kmu[from] = MatrixCreate(k,1);

		for (const Link* link=from->Next(); link!=from; link=link->Next())	{
			RecursiveKalmanCreate(link->Out());
		}

	}

	void RecursiveKalmanDelete(const Link* from)	{
		for (const Link* link=from->Next(); link!=from; link=link->Next())	{
			RecursiveKalmanDelete(link->Out());
		}

		MatrixDelete(Kmu[from],k);
		MatrixDelete(Kgamma[from],k);
		MatrixDelete(KF[from],k);
		MatrixDelete(Ky[from],k);
		delete KD[from];
		delete KM[from];
		delete KLambda[from];
	}

	void KalmanBackward(const Link* from)	{

		CovMatrix* covM = KM[from];
		double** M = covM->GetMatrix();

		if (from->isLeaf())	{
			int clamp = 0;
			for (int i=0; i<k; i++)	{
				if (GetMultiNormal(from)->ClampVector[index1[i]])	{
					clamp++;
				}
			}
			if ((clamp > 0) && (clamp < k))	{
				cerr << "error in kalman backward: clamping profile\n";
				exit(1);
			}
			if (clamp == k)	{

				double nu = 1.0 / GetBranchLength(from->GetBranch());
				MatrixSet(M,OmegaW,k,k);
				MatrixScalarProduct(M,nu,k,k);

				double** x = MatrixCreate(k,1);
				for (int i=0; i<k; i++)	{
					x[i][0] = (*GetMultiNormal(from))[index1[i]];
				}
				double** gamma = Kgamma[from];
				MatrixProduct(M,x,gamma,k,k,1);

				double** z = MatrixCreate(n,1);
				for (int i=0; i<n; i++)	{
					z[i][0] = - (*GetMultiNormal(from))[index2[i]] + (*GetMultiNormal(from->Out()))[index2[i]];
				}
				double** E = MatrixCreate(k,1);
				MatrixProduct(Omegaxy,z,E,k,n,1);

				double** yy = Ky[from];
				double** y = MatrixCreate(k,1);
				for (int i=0; i<k; i++)	{
					y[i][0] = yy[index1[i]][0];
				}
				MatrixAdd(E,y,k,1);
				MatrixScalarProduct(E,-nu,k,1);
				MatrixAdd(gamma,E,k,1);

				MatrixDelete(E,k);
				MatrixDelete(z,n);
				MatrixDelete(x,k);
				MatrixDelete(y,k);
			}
			else {
				MatrixSet(M,0,k,k);
				double** gamma = Kgamma[from];
				MatrixSet(gamma,0,k,1);
			}
		}
		else	{
			// here D is K in the manuscript
			// and gamma is M*beta in the manuscript
			// compute D and mu from the Ms and gammas of children
			CovMatrix* covD = KD[from];
			double** D = covD->GetMatrix();
			MatrixSet(D,0,k,k);
			double** tmp = MatrixCreate(k,1);
			MatrixSet(tmp,0,k,1);

			for (const Link* link=from->Next(); link!=from; link=link->Next())	{

				KalmanBackward(link->Out());

				MatrixAdd(D,KM[link->Out()]->GetMatrix(),k,k);
				MatrixAdd(tmp,Kgamma[link->Out()],k,1);
			}

			// covD->Invert();
			double** mu = Kmu[from];
			if (covD->GetMax() < 1e-8)	{
				MatrixSet(mu,0,k,1);
			}
			else	{
				covD->CorruptDiag();
				covD->Diagonalise();
				MatrixProduct(covD->GetInvMatrix(),tmp,mu,k,k,1);
			}
			// cerr << "internal D mu : " << D[0][0] << '\t' << mu[0][0] << '\n';
			MatrixDelete(tmp,k);

			// now compute own M and gamma
			if (! from->isRoot())	{

			// inverse branch length
			double nu = 1.0 / GetBranchLength(from->GetBranch());

			// compute Lambda and its inverse
			double** Lambda = KLambda[from]->GetMatrix();
			MatrixSet(Lambda,OmegaW,k,k);
			MatrixScalarProduct(Lambda,nu,k,k);
			MatrixAdd(Lambda,D,k,k);
			// KLambda[from]->Invert();
			KLambda[from]->CorruptDiag();
			KLambda[from]->Diagonalise();
			double** InvLambda = KLambda[from]->GetInvMatrix();

			// M <- v (Omega + W) - v^2 (Omega + W) Lambda^-1 (Omega + W)
			double** G = MatrixCreate(k,k);
			double** H = MatrixCreate(k,k);
			MatrixProduct(OmegaW,InvLambda,G,k,k,k);
			MatrixScalarProduct(G,-nu,k,k);
			MatrixSetIdentity(H,k);
			MatrixAdd(H,G,k,k);
			MatrixScalarProduct(H,nu,k,k);
			double** M = KM[from]->GetMatrix();
			MatrixProduct(H,OmegaW,M,k,k,k);
			MatrixScalarProduct(G,-1,k,k);
			// now, G = nu (Omega + W) Lambda^-1

			double** gamma = Kgamma[from];

			double** z = MatrixCreate(n,1);
			for (int i=0; i<n; i++)	{
				z[i][0] = - (*GetMultiNormal(from))[index2[i]] + (*GetMultiNormal(from->Out()))[index2[i]];
			}
			double** E = MatrixCreate(k,1);
			MatrixProduct(Omegaxy,z,E,k,n,1);

			double** yy = Ky[from];
			double** y = MatrixCreate(k,1);
			for (int i=0; i<k; i++)	{
				y[i][0] = yy[index1[i]][0];
			}
			MatrixAdd(E,y,k,1);
			MatrixScalarProduct(E,nu,k,1);

			double** F = KF[from];
			MatrixProduct(D,mu,F,k,k,1);
			MatrixAdd(F,E,k,1);

			MatrixProduct(G,F,gamma,k,k,1);
			MatrixScalarProduct(E,-1,k,1);
			MatrixAdd(gamma,E,k,1);

			MatrixDelete(E,k);
			MatrixDelete(z,n);
			MatrixDelete(y,k);
			MatrixDelete(G,k);
			MatrixDelete(H,k);
			// cerr << "internal M gamma : " << M[0][0] << '\t' << gamma[0][0] << '\n';
			}
		}
	}

	void KalmanForward(const Link* from, const Link* parent = 0)	{
		if (from->isRoot())	{

			int clamp = 0;
			for (int i=0; i<k; i++)	{
				if (GetMultiNormal(from)->ClampVector[index1[i]])	{
					clamp++;
				}
			}
			if ((clamp > 0) && (clamp < k))	{
				cerr << "error in kalman forward: root clamping profile\n";
				exit(1);
			}
			if (!clamp)	{
				// draw N(mu,D-1)
				KD[from]->CorruptDiag();
				KD[from]->Diagonalise();
				double** mu = Kmu[from];
				for (int i=0; i<k; i++)	{
					if (GetMultiNormal(from)->ClampVector[i])	{
						mu[i][0] = (*GetMultiNormal(from))[i];
					}
				}
				GetMultiNormal(from)->ConditionalDrawNormalFromPrecision(Kmu[from],KD[from],index1,k);
			}
		}
		else	{

			if (from->isLeaf())	{
				int clamp = 0;
				for (int i=0; i<k; i++)	{
					if (GetMultiNormal(from)->ClampVector[index1[i]])	{
						clamp++;
					}
				}
				if ((clamp > 0) && (clamp < k))	{
					cerr << "error in kalman forward: root clamping profile\n";
					exit(1);
				}
				if (! clamp)	{
					double** alpha = MatrixCreate(k,1);

					double nu = 1.0 / GetBranchLength(from->GetBranch());

					double** Lambda = KLambda[from]->GetMatrix();
					MatrixSet(Lambda,OmegaW,k,k);
					MatrixScalarProduct(Lambda,nu,k,k);
					KLambda[from]->CorruptDiag();
					KLambda[from]->Diagonalise();
					double** InvLambda = KLambda[from]->GetInvMatrix();

					double** yy = Ky[from];
					double** y = MatrixCreate(k,1);
					for (int i=0; i<k; i++)	{
						y[i][0] = yy[index1[i]][0];
					}
					MatrixProduct(InvLambda,y,alpha,k,k,1);

					double** z = MatrixCreate(n,1);
					for (int i=0; i<n; i++)	{
						z[i][0] = - (*GetMultiNormal(from))[index2[i]] + (*GetMultiNormal(from->Out()))[index2[i]];
					}
					double** E = MatrixCreate(k,1);
					MatrixProduct(Omegaxy,z,E,k,n,1);
					MatrixAdd(alpha,E,k,1);

					MatrixScalarProduct(alpha,nu,k,1);
					for (int i=0; i<k; i++)	{
						alpha[i][0] += (*GetMultiNormal(parent))[i];
					}

					GetMultiNormal(from)->ConditionalDrawNormalFromPrecision(alpha,KLambda[from],index1,k);

					MatrixDelete(alpha,k);
					MatrixDelete(E,k);
					MatrixDelete(y,k);
					MatrixDelete(z,n);
				}
			}
			else	{

				// grep value of parent x
				double** xup = MatrixCreate(k,1);
				for (int i=0; i<k; i++)	{
					xup[i][0] = (*GetMultiNormal(parent))[i];
				}

				double nu = 1.0 / GetBranchLength(from->GetBranch());

				double** tmp = MatrixCreate(k,1);
				double** alpha = MatrixCreate(k,1);

				MatrixProduct(OmegaW,xup,tmp,k,k,1);
				MatrixScalarProduct(tmp,nu,k,1);
				MatrixAdd(tmp,KF[from],k,1);

				KLambda[from]->CorruptDiag();
				KLambda[from]->Diagonalise();
				double** InvLambda = KLambda[from]->GetInvMatrix();
				MatrixProduct(InvLambda,tmp,alpha,k,k,1);

				// draw N(alpha,Lambda-1)
				// cerr << "from lambda : " << KLambda[from]->GetMatrix()[0][0] << '\n';
				GetMultiNormal(from)->ConditionalDrawNormalFromPrecision(alpha,KLambda[from],index1,k);

				MatrixDelete(tmp,k);
				MatrixDelete(alpha,k);
				MatrixDelete(xup,k);

			}
		}

		// recursive call
		for (const Link* link=from->Next(); link!=from; link=link->Next())	{
			KalmanForward(link->Out(),from);
		}
	}

	int k;
	int n;
	int* index;
	int* index1;
	int* index2;

	double** Omega;
	double** Omegaxy;
	double** OmegaW;
	CovMatrix* CovOmegaW;
	double** W;
	double** WW;
	map<const Link*, CovMatrix*> KLambda;
	map<const Link*, double**> KF;
	map<const Link*, double**> Ky;
	map<const Link*, CovMatrix*> KM;
	map<const Link*, double**> Kgamma;
	map<const Link*, CovMatrix*> KD;
	map<const Link*, double**> Kmu;

};
*/

class MultiVariateKalmanMove : public MCUpdate, public Mnode {

	KalmanMultiVariateTreeProcess* tree;
	int k;
	int n;

	public:

	MultiVariateKalmanMove(KalmanMultiVariateTreeProcess* intree, int ink, int inn){
		tree = intree;
		k = ink;
		n = inn;
		tree->RecursiveRegister(this,tree->GetRoot());
	}

	double Move(double tuning_modulator){
		Corrupt(false);
		tree->KalmanMove(k,n);
		Update();
		return 1;
	}
};

class ConjugateKalmanMultiVariateTreeProcess : public KalmanMultiVariateTreeProcess, public ConjugateMultiVariateTreeProcess	{

	public:

	ConjugateKalmanMultiVariateTreeProcess() {}

	ConjugateKalmanMultiVariateTreeProcess(ConjugateInverseWishart* insigma, LengthTree* intree, LengthTree* inscaletree = 0, Var<RealVector>* indrift = 0, Var<RealVector>* inrootmean = 0, Var<PosRealVector>* inrootvar = 0) : ConjugateMultiVariateTreeProcess(insigma,intree,inscaletree,indrift,inrootmean,inrootvar) {}

};

class ConjugateExternalKalmanMultiVariateTreeProcess : public ExternalKalmanMultiVariateTreeProcess, public ConjugateMultiVariateTreeProcess	{

	public:

	ConjugateExternalKalmanMultiVariateTreeProcess() {}

	ConjugateExternalKalmanMultiVariateTreeProcess(ConjugateInverseWishart* insigma, LengthTree* intree, LengthTree* inscaletree = 0, Var<RealVector>* indrift = 0, Var<RealVector>* inrootmean = 0, Var<PosRealVector>* inrootvar = 0) : ConjugateMultiVariateTreeProcess(insigma,intree,inscaletree,indrift,inrootmean,inrootvar) {}

};

class ConjugateConditionalExternalKalmanMultiVariateTreeProcess : public ConditionalExternalKalmanMultiVariateTreeProcess, public ConjugateMultiVariateTreeProcess	{

	public:

	ConjugateConditionalExternalKalmanMultiVariateTreeProcess() {}

	ConjugateConditionalExternalKalmanMultiVariateTreeProcess(ConjugateInverseWishart* insigma, LengthTree* intree, LengthTree* inscaletree = 0, Var<RealVector>* indrift = 0, Var<RealVector>* inrootmean = 0, Var<PosRealVector>* inrootvar = 0) : ConjugateMultiVariateTreeProcess(insigma,intree,inscaletree,indrift,inrootmean,inrootvar) {}

};

#endif

