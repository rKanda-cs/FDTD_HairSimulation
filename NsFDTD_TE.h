#ifndef _NS_FDTD_TE_
#define _NS_FDTD_TE_
#include "FDTD_TE.h"

class NsFDTD_TE: public FDTD_TE{
	typedef FDTD_TE super;
private:
	double R_P, R_M;
public:
	NsFDTD_TE();
	~NsFDTD_TE();
	bool calc();
	void field();
	void PMLfield();
private: 
	complex<double> Dx2_n(complex<double> *p, int i, int j, int t){
		return (p[index(i,j+1, t)] + p[index(i, j-1, t)] - p[index(i-1,j+1, t)] - p[index(i-1, j-1, t)])/2.0;
	};

	complex<double> Dy2_n(complex<double> *p, int i, int j, int t){
		return (p[index(i+1,j, t)] + p[index(i-1, j, t)] - p[index(i+1,j-1, t)] - p[index(i-1, j-1, t)])/2.0;
	};

	void CalcE(){
		int Npml = mField->getNpml();
	#ifdef _OPENMP
	#pragma omp parallel
	#endif
		{
#ifdef _OPENMP
#pragma omp for
#endif
			//ìdäEÇÃåvéZEx
			for (int i = Npml; i < mField->getNpx() - Npml; i++) {
				for (int j = Npml; j < mField->getNpy() - Npml; j++) {
					EX(i, j, +1) = CEX(i, j)*EX(i, j, 0)
						+ CEXLY(i, j)*(HZ(i, j+1, 0) - HZ(i, j, 0));
				}
			}

#ifdef _OPENMP
#pragma omp for
#endif
			//ìdäEÇÃåvéZEy
			for (int i = Npml; i < mField->getNpx() - Npml; i++) {
				for (int j = Npml; j < mField->getNpy() - Npml; j++) {
					EY(i, j, +1) = CEY(i, j)*EY(i, j, 0)
						- CEYLX(i, j)*(HZ(i+1, j, 0) - HZ(i, j, 0));
				}
			}

		}
	}

	//é•äEÇÃåvéZ Hz(i+1/2, j+1/2) -> Hz[i,j]
	void CalcH(){
		int Npml = mField->getNpml();
#ifdef _OPENMP
#pragma omp parallel for
#endif
		for (int i = Npml; i < mField->getNpx() - Npml; i++) {
			for (int j = Npml; j < mField->getNpy() - Npml; j++) {
				HZ(i, j, +1) = HZ(i, j, 0)
					- CHZLH(i, j)*(R_P *(EY(i, j, +1) - EY(i-1, j, +1) - (EX(i, j, +1) - EX(i, j-1, +1)))
						+ R_M * (Dx2_n(Ey, i, j, +1) - Dy2_n(Ex, i, j, +1)));
			}
		}
	}

	//Ns_PML
	void CalcE_PML() {
		int Npml = mField->getNpml();
#ifdef _OPENMP
#pragma omp parallel
#endif
		{
#ifdef _OPENMP
#pragma omp for
#endif
			//ìdäEÇÃåvéZEx
			for (int i = 1; i < mField->getNpx() - 1; i++) {
				for (int j = 0; j < mField->getNpy() - 1; j++) {
					if (i < Npml || i > mField->getNpx() - Npml - 1 || j < Npml || j > mField->getNpy() - Npml - 1)
						EX(i, j, +1) = BEYP(i, j)*BEYM(i, j)*EX(i, j, 0)
						+ BEYP(i, j)*CEXLY(i, j)*(HZ(i, j+1, 0) - HZ(i, j, 0));
				}
			}

#ifdef _OPENMP
#pragma omp for
#endif
			//ìdäEÇÃåvéZEy
			for (int i = 0; i < mField->getNpx() - 1; i++) {
				for (int j = 1; j < mField->getNpy() - 1; j++) {
					if (i < Npml || i > mField->getNpx() - Npml - 1 || j < Npml || j > mField->getNpy() - Npml - 1)
						EY(i, j, +1) = BEXP(i, j)*BEXM(i, j)*EY(i, j, 0)
						- BEXP(i, j)*CEYLX(i, j)*(HZ(i+1, j, 0) - HZ(i, j, 0));
				}
			}

		}
	}

	//é•äEÇÃåvéZ Hz(i+1/2, j+1/2) -> Hz[i,j]
	void CalcH_PML() {
		int Npml = mField->getNpml();
#ifdef _OPENMP
#pragma omp parallel
#endif
		{
#ifdef _OPENMP
#pragma omp for
#endif
			for (int i = 1; i < mField->getNpx() - 1; i++) {
				for (int j = 1; j < mField->getNpy() - 1; j++) {
					if (i < Npml || i > mField->getNpx() - Npml - 1 || j < Npml || j > mField->getNpy() - Npml - 1)
						HZX(i, j, +1) = BHZXP(i, j)*BHZXM(i, j)*HZX(i, j, 0)
						- BHZXP(i, j)*CHZLH(i, j)*(R_P * (EY(i, j, +1) - EY(i-1, j, +1))
							+ R_M * Dx2_n(Ey, i, j, +1));
				}
			}
#ifdef _OPENMP
#pragma omp for
#endif
			for (int i = 1; i < mField->getNpx() - 1; i++) {
				for (int j = 1; j < mField->getNpy() - 1; j++) {
					if (i < Npml || i > mField->getNpx() - Npml - 1 || j < Npml || j > mField->getNpy() - Npml - 1)
						HZY(i, j, +1) = BHZYP(i, j)*BHZYM(i, j)*HZY(i, j, 0)
						+ BHZYP(i, j)*CHZLH(i, j)*(R_P *(EX(i, j, +1) - EX(i, j-1, +1))
							+ R_M * Dy2_n(Ex, i, j, +1));
				}
			}
#ifdef _OPENMP
#pragma omp for
#endif
			for (int i = 1; i < mField->getNpx() - 1; i++)
				for (int j = 1; j < mField->getNpy() - 1; j++)
					if (i < Npml || i > mField->getNpx() - Npml - 1 || j < Npml || j > mField->getNpy() - Npml - 1)
						HZ(i, j, +1) = HZX(i, j, +1) + HZY(i, j, +1);
		}
	}

	void absorbing();
	bool EndTask();

	void ReStart(){
		super::Initialize();
		field();
	}
};
#endif //_NS_FDTD_TE	

			/*
	//todo i+1 - i Ç≈ÇÕÇ»Ç≠ i - (i-1) jÇ‡ìØól HÇÃï˚Ç™+0.5êiÇÒÇ≈Ç¢ÇÈÇ©ÇÁ
#ifdef _OPENMP
#pragma omp parallel
#endif
{
    #ifdef _OPENMP
    #pragma omp for
    #endif
	//ìdäEÇÃåvéZEx
	for(int i=2; i<mField->getNx()-1; i++)
		for(int j=2; j<mField->getNy(); j++)
			EX(i,j, +1) = CEX(i,j)*EX(i,j, 0)
								+ CEXLY(i,j)*( R_P*(HZ(i,j, 0) - HZ(i,j-1, 0))
													+ R_M*Dy2_n(Hz,i,j-1, 0)
													);
    #ifdef _OPENMP
    #pragma omp for
    #endif
	//ìdäEÇÃåvéZEy
	for(int i=2; i<mField->getNx(); i++)
		for(int j=2; j<mField->getNy()-1; j++)
			EY(i,j, +1) = CEY(i,j)*EY(i,j, 0) 
								- CEYLX(i,j)*( R_P*(HZ(i,j, 0) - HZ(i-1,j, 0))
													+ R_M*Dx2_n(Hz,i-1,j, 0)
													);
}
	}

	//é•äEÇÃåvéZ Hz(i+1/2, j+1/2) -> Hz[i,j]
	//todo i - (i-1)Ç≈ÇÕÇ»Ç≠ i+1 - i ? jÇ‡ìØól HÇÃï˚Ç™+0.5êiÇÒÇ≈Ç¢ÇÈÇ©ÇÁ
	void CalcH(){
#ifdef _OPENMP
#pragma omp parallel for
#endif
	for(int i=1; i<mField->getNx()-1; i++)
		for(int j=1; j<mField->getNy()-1; j++)
			HZ(i,j, +1) = HZ(i,j, 0)
								-CHZLH(i,j)*( EY(i+1,j  , +1) - EY(i,j, +1) )
								+CHZLH(i,j)*( EX(i  ,j+1, +1) - EX(i,j, +1) );
	}
	*/
