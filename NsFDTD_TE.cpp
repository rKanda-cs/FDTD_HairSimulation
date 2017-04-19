#include"NsFDTD_TE.h"
#include <omp.h>

NsFDTD_TE::NsFDTD_TE()
:FDTD_TE()
{
	cout << "NsFDTD_TE Constructor" << endl;
};

NsFDTD_TE::~NsFDTD_TE(){
	cout << "NsFDTD_TE Destructor" << endl;
};

bool NsFDTD_TE::calc(){	

	CalcE();	//電界の計算
	CalcE_PML();

	NsScatteredWave(wave_angle);
	//absorbing();

	CalcH();		//磁界の計算 Hz(i+1/2, j+1/2) -> Hz[i,j]
	CalcH_PML();

	//pointLightSource(Hz);

	if (time > maxStep) {
		MiePrint(Ey, "time" + to_s(maxStep) + "_PML" + to_s(mField->getNpml()) + "_NsTE_");
		capture(to_s(time));
		//capture(to_s(Inv_Nano_S(lambda_s)));
		return EndTask();
	}
	return true;
};


bool NsFDTD_TE::EndTask(){
	cout << "End Task" << endl;
	string label = "";

	NTFFindexform(label, NTFF::NTFFDATA | NTFF::TOTAL);

	//終了条件の確認
	if( !Terminate())
		return false;

	ReStart();

	return true;
}

void NsFDTD_TE::field(){		
	super::field();	//誘電率の設定
	setWorkingDirPass(MakeDir("Ns"));

	R_M = NsCoef();		//差分演算用の計算定数の設定
	R_P = 1.0-R_M;

	//σ=0を用いて最適化した定数の計算
	double mu;
	double sig = 0;	//真空の誘電率, 導電率, 透磁率
	for(int i=0; i<mField->getNpx(); i++){
		for(int j=0; j<mField->getNpy(); j++){
			mu = MU_0_S;
			double u_hz = sin(w_s / sqrt(EPSHZ(i,j) / EPSILON_0_S) * DT_S / 2) / sin(k_s*DT_S/2);
			double u_ex = sin(w_s / sqrt(EPSEX(i,j) / EPSILON_0_S) * DT_S / 2) / sin(k_s*DT_S/2);
			double u_ey = sin(w_s / sqrt(EPSEY(i,j) / EPSILON_0_S) * DT_S / 2) / sin(k_s*DT_S/2);
			CEX(i,j)   = 1.0;
			CEXLY(i,j) = u_ex*sqrt(mu/EPSEX(i,j));

			CEY(i,j)   = 1.0;
			CEYLX(i,j) = u_ey*sqrt(mu/EPSEY(i,j));	

			CHZLH(i,j) = u_hz*sqrt(EPSHZ(i,j)/mu);;

		}
	}
	cout << "Ns_TE_field" << endl;

	PMLfield();
}

void NsFDTD_TE::PMLfield() {
	double mu = MU_0_S;

	for (int i = 0; i < mField->getNpx(); i++) {
		for (int j = 0; j < mField->getNpy(); j++) {
			double sig_x = mField->sigmaX(i, j);			//σx, σx*, σy, σy* 　　<- B-PMLの係数
			double sig_xx = mu / EPSILON_0_S * sig_x;
			double sig_y = mField->sigmaY(i, j);
			double sig_yy = mu / EPSILON_0_S * sig_y;

			if (SIG(i, j) != 0) {
				sig_x = SIG(i, j);
				sig_y = SIG(i, j);
				sig_xx = mu / EPSILON_0_S * sig_x;
				sig_yy = mu / EPSILON_0_S * sig_y;
			}

			double ax = sig_x * DT_S / (2 * EPSEX(i, j));
			double ay = sig_y * DT_S / (2 * EPSEY(i, j));
			double axx = sig_xx * DT_S / (2 * mu);
			double ayy = sig_yy * DT_S / (2 * mu);
			
			/*
			BEZXP(i, j) = 1 + tanh(ax);
			BEZXM(i, j) = 1 - tanh(ax);
			BEZYP(i, j) = 1 + tanh(ay);
			BEZYM(i, j) = 1 - tanh(ay);
			BHXP(i, j) = 1 + tanh(axx);
			BHXM(i, j) = 1 - tanh(axx);
			BHYP(i, j) = 1 + tanh(ayy);
			BHYM(i, j) = 1 - tanh(ayy);
			*/
			BEXP(i, j) = 1 + (tanh(ax) / (1 + tanh(ax)*tanh(axx)));
			BEXM(i, j) = 1 - (tanh(ax) / (1 + tanh(ax)*tanh(axx)));
			BEYP(i, j) = 1 + (tanh(ay) / (1 + tanh(ay)*tanh(ayy)));
			BEYM(i, j) = 1 - (tanh(ay) / (1 + tanh(ay)*tanh(ayy)));
			BHZXP(i, j) = 1 + (tanh(axx) / (1 + tanh(ax)*tanh(axx)));
			BHZXM(i, j) = 1 - (tanh(axx) / (1 + tanh(ax)*tanh(axx)));
			BHZYP(i, j) = 1 + (tanh(ayy) / (1 + tanh(ay)*tanh(ayy)));
			BHZYM(i, j) = 1 - (tanh(ayy) / (1 + tanh(ay)*tanh(ayy)));

			BEXP(i, j) = 1 / BEXP(i, j);
			BEYP(i, j) = 1 / BEYP(i, j);
			BHZXP(i, j) = 1 / BHZXP(i, j);
			BHZYP(i, j) = 1 / BHZYP(i, j);
		}
	}
	cout << "PML_field" << endl;
}

void NsFDTD_TE::absorbing(){
	//H(i,j)だったやつ（今までちゃんと動いていたやつ）
	
	absorbing_nsRL(Ey, 0,	 LEFT);
	absorbing_nsRL(Ey, mField->getNpx()-2, RIGHT);
	absorbing_nsTB(Ex, 0,    BOTTOM);
	absorbing_nsTB(Ex, mField->getNpy()-2, TOP);
	/*
	absorbing_nsRL(Ey, 1,	 LEFT);
	absorbing_nsRL(Ey, mField->getNx()-1, RIGHT);
	absorbing_nsTB(Ex, 1,    BOTTOM);
	absorbing_nsTB(Ex, mField->getNy()-1, TOP);
	*/
}



//double a = 0; //α = σ/(2ε) より, 今はσ=0としているのでαも0
//double u = sqrt( (_pow(sin(sqrt(w_s*w_s - a*a)*DT_S/2 ),2)+_pow( sinh(a*DT_S/2),2) )/ (_pow(sin(k_s*DT_S/2),2)*cosh(a*DT_S))  );
// a = 0, tanh(0) = sinh(0) = 0, cosh(0) = 1　を用いて最適化する


