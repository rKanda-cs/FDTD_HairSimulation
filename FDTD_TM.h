#ifndef _FDTD_TM_H
#define _FDTD_TM_H
#define _USE_MATH_DEFINES
#include "Solver.h"
#include "Object.h"
#include<math.h>

enum NTF_FLAG{
	ALLSTR = true,
	ONLYTOTAL = false
};
const int ANGLE = 180;

class FDTD_TM: public Solver{
	typedef Solver super;

protected:
	Complex *Ez, *Hx, *Hy;
	Complex *Ezx, *Ezy; //pml用
	double *C_EZ, *C_EZLH, *C_HXLY, *C_HYLX;
	double *C_EZX, *C_EZY, *C_EZXLX, *C_EZYLY;
	double *C_HX, *C_HY;
	double *EPS_EZ, *EPS_HX, *EPS_HY;
	double *B_EZXp, *B_EZXm, *B_EZYp, *B_EZYm; //pml付きNsFDTD用
	double *B_HXp, *B_HXm, *B_HYp, *B_HYm;

public:
	FDTD_TM();
	virtual ~FDTD_TM();
	virtual bool calc()=0;
	void draw();
	void field();

protected:
	virtual void Initialize();

	void NTFFindexform(string label, NTFF::output f = NTFF::REFLEC);

	void PML();

	void NsScatteredWave(int angle);	//散乱波
	void IncidentWave(int angle);
	void IncidentWaveH(int angle);

	//ゲッター(インライン化することでオーバーヘッド削減, なるか知らんけど)
	Complex& EZ(const int &i, const int &j, const int &t){
		return Ez[pmlIndex(i,j,t)];
	};

	Complex& HX(const int &i, const int &j, const int &t){
		return Hx[pmlIndex(i,j,t)];
	};

	Complex& HY(const int &i, const int &j, const int &t){
		return Hy[pmlIndex(i,j,t)];
	};

	Complex& EZX(const int &i, const int &j, const int &t){
		return Ezx[pmlIndex(i,j,t)];
	}
	
	Complex& EZY(const int &i, const int &j, const int &t){
		return Ezy[pmlIndex(i,j,t)];
	}

	Complex& EZ(const int &i, const int &j) {
		return Ez[pmlIndex(i, j)];
	};

	Complex& HX(const int &i, const int &j) {
		return Hx[pmlIndex(i, j)];
	};

	Complex& HY(const int &i, const int &j) {
		return Hy[pmlIndex(i, j)];
	};

	Complex& EZX(const int &i, const int &j) {
		return Ezx[pmlIndex(i, j)];
	}

	Complex& EZY(const int &i, const int &j) {
		return Ezy[pmlIndex(i, j)];
	}

	double& CEZX(const int &i, const int &j){
		return C_EZX[pmlIndex(i,j)];
	}

	double& CEZY(const int &i, const int &j){
		return C_EZY[pmlIndex(i,j)];
	}

	double& CEZXLX(const int &i, const int &j){
		return C_EZXLX[pmlIndex(i,j)];
	}

	double& CEZYLY(const int &i, const int &j){
		return C_EZYLY[pmlIndex(i,j)];
	}

	double& CEZ(const int &i, const int &j){
		return C_EZ[pmlIndex(i,j)];
	};

	double& CEZLH(const int &i, const int &j){
		return C_EZLH[pmlIndex(i,j)];
	};

	double& CHXLY(const int &i, const int &j){
		return C_HXLY[pmlIndex(i,j)];
	};

	double& CHYLX(const int &i, const int &j){
		return C_HYLX[pmlIndex(i,j)];
	};

	double& CHX(const int &i, const int &j){
		return C_HX[pmlIndex(i,j)];
	};

	double& CHY(const int &i, const int &j){
		return C_HY[pmlIndex(i,j)];
	};
	double& EPSEZ(const int &i, const int &j){
		return EPS_EZ[pmlIndex(i,j)];
	};

	double& EPSHX(const int &i, const int &j){
		return EPS_HX[pmlIndex(i,j)];
	};

	double& EPSHY(const int &i, const int &j){
		return EPS_HY[pmlIndex(i,j)];
	};

	double& BEZXP(const int &i, const int &j) {
		return B_EZXp[pmlIndex(i, j)];
	};

	double& BEZXM(const int &i, const int &j) {
		return B_EZXm[pmlIndex(i, j)];
	};

	double& BEZYP(const int &i, const int &j) {
		return B_EZYp[pmlIndex(i, j)];
	};

	double& BEZYM(const int &i, const int &j) {
		return B_EZYm[pmlIndex(i, j)];
	};

	double& BHXP(const int &i, const int &j) {
		return B_HXp[pmlIndex(i, j)];
	};

	double& BHXM(const int &i, const int &j) {
		return B_HXm[pmlIndex(i, j)];
	};

	double& BHYP(const int &i, const int &j) {
		return B_HYp[pmlIndex(i, j)];
	};

	double& BHYM(const int &i, const int &j) {
		return B_HYm[pmlIndex(i, j)];
	};

	complex<double> EZ_NTF(const int &i, const int &j){
	//	return EZ(i,j,0);	//todo ?

		double _sin  = 4*pow(sin( w_s/ sqrt(EPSEZ(i,j)) / 2.0), 2.0);
		double _sin2 = 2*sin(w_s/ sqrt(EPSEZ(i,j)));

		double ez_r = -(EZ(i,j,+1).real() + EZ(i,j,-1).real() - EZ(i,j, 0).real() )/_sin;
		double ez_i = -(EZ(i,j,+1).real() - EZ(i,j,-1).real() ) / _sin2;
		//return complex<double>(ez_r, ez_i);
		
		//_sin = 1.0;
		return -( EZ(i,j, +1) + EZ(i,j, -1) - EZ(i,j, 0))/_sin;

		//complex<double> ez_i = -( EZ(i,j,+1) - EZ(i,j, -1) ) / (2*sin(w_s/EPSEZ(i,j)*DT_S)/DT_S);

		//return complex<double>(ez_r.real(), ez_i.imag());
	}

	complex<double> HX_NTF(const int &i, const int &j){
	//	return HX(i,j,0);	//todo ?

		double _sin = 4*pow(sin(w_s/sqrt(EPSHX(i,j))/2.0), 2.0);
		double _sin2 = 2*sin(w_s/ sqrt(EPSHX(i,j)));
		double hx_r = -( HX(i,j, +1).real() + HX(i,j, -1).real() - HX(i,j, 0).real())/_sin;
		double hx_i = -( HX(i,j, +1).real() - HX(i,j, -1).real() ) / _sin2;
		
		//_sin = 1.0;
		//return complex<double>(hx_r, hx_i);
		return -(HX(i,j, +1) + HX(i,j, -1) - HX(i,j, 0))/_sin;
	}

	complex<double> HY_NTF(const int &i, const int &j){
		return HY(i,j,0);

		double _sin = 4*pow(sin(w_s/sqrt(EPSHY(i,j)/EPSILON_0_S)*DT_S/2.0),2.0)/DT_S/DT_S;
		double _sin2 = 2*sin(w_s/ sqrt(EPSHY(i,j)));
		double hy_r = -( HY(i,j, +1).real() + HY(i,j, -1).real() - HY(i,j, 0).real())/_sin;
		double hy_i = -( HY(i,j, +1).real() - HY(i,j, -1).real() ) / _sin2;
		
		//_sin = 1.0;
		//return complex<double>(hy_r, hy_i);
		return -(HY(i,j, +1) + HY(i,j, -1) - HY(i,j, 0))/_sin;
	}

	complex<double> EZi_NTF(const int &i, const int &j){
		double _sin = 4*pow(sin(w_s/sqrt(EPSEZ(i,j)/EPSILON_0_S)*DT_S/2.0),2.0)/DT_S/DT_S;
		double rad = wave_angle/180.0*M_PI;
		/*
		complex<double> ez_i = ( polar(1.0, k_s*(i*cos(rad) + j*sin(rad)) - w_s*(time+DT_S))
							   - polar(1.0, k_s*(i*cos(rad) + j*sin(rad)) - w_s*(time-DT_S))
							   ) / (2*sin(w_s/EPSEZ(i,j)*DT_S)/DT_S);
							   */
		complex<double> ez_i = -( EZ(i,j,+1) - EZ(i,j, -1) ) / (2*sin(w_s/EPSEZ(i,j)*DT_S)/DT_S);
		return -(EZ(i,j, +1) + EZ(i,j, -1) - EZ(i,j, 0))/_sin + ez_i;
	}

	complex<double> HXi_NTF(const int &i, const int &j){
		double _sin = 4*pow(sin(w_s/sqrt(EPSHX(i,j)/EPSILON_0_S)*DT_S/2.0),2.0)/DT_S/DT_S;

		//_sin = 1.0;
		return -(HX(i,j, +1) + HX(i,j, -1) - HX(i,j, 0))/_sin;
	}

	complex<double> HYi_NTF(const int &i, const int &j){
		double _sin = 4*pow(sin(w_s/sqrt(EPSHY(i,j)/EPSILON_0_S)*DT_S/2.0),2.0)/DT_S/DT_S;

		//_sin = 1.0;
		return -(HY(i,j, +1) + HY(i,j, -1) - HY(i,j, 0))/_sin;
	}

	void SaveData(string prefix = "");
	void OpenData(string prefix = "");
};

#endif //_FDTD_TM_H

				