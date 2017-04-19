#include "Object.h"
#include "StFDTD_TM.h"

using namespace std;

StFDTD_TM::StFDTD_TM()
:FDTD_TM()
{
	cout << "StFDTD_TM Constructor" << endl;
}

StFDTD_TM::~StFDTD_TM(){
	cout << "StFDTD_TM Destructor" << endl;
}

void StFDTD_TM::field(){
	super::field();
	setWorkingDirPass(MakeDir("St"));
	setWorkingDirPass(MakeDir(to_s(wave_angle) + "deg_lambda" + to_s((int)Inv_Nano_S(lambda_s))));

	double sig_x, sig_y, sig_xx, sig_yy;//σx, σx*, σy, σy* 　　<- B-PMLの係数
	for(int i=0; i<mField->getNpx(); i++){
		for(int j=0; j<mField->getNpy(); j++){
			sig_x = mField->sigmaX(i,j);
			sig_xx = MU_0_S/EPSILON_0_S * sig_x;
			sig_y = mField->sigmaY(i,j);
			sig_yy = MU_0_S/EPSILON_0_S * sig_y;

			if (SIG(i, j) != 0) {
				sig_x = SIG(i, j);
				sig_y = SIG(i, j);
				sig_xx = MU_0_S / EPSILON_0_S * sig_x;
				sig_yy = MU_0_S / EPSILON_0_S * sig_y;
			}

			//Δt = 1, μ(i,j) = μ0 
 			CEZX(i,j)   = MaxwellCoef(EPSEZ(i,j), sig_x);		// 1- σ/ε			  1/ε
			CEZXLX(i,j) = MaxwellCoef2(EPSEZ(i,j), sig_x);		// --------			--------
																// 1+ σ/ε			1+ σ/ε
			CEZY(i,j)   = MaxwellCoef(EPSEZ(i,j), sig_y);
			CEZYLY(i,j) = MaxwellCoef2(EPSEZ(i,j), sig_y);

			CHX(i,j)    = MaxwellCoef(MU_0_S, sig_yy);			// 1- σ/μ			  1/μ
			CHXLY(i,j)  = MaxwellCoef2(MU_0_S, sig_yy);			// --------			--------
																// 1+ σ/μ			1+ σ/μ
			CHY(i,j)    = MaxwellCoef(MU_0_S, sig_xx);
			CHYLX(i,j)  = MaxwellCoef2(MU_0_S, sig_xx);
		}
	}
}

bool StFDTD_TM::calc(){

	CalcE();
	//EZX(mField->getNpx()/2, mField->getNpy()/2) += 0.5*ray_coef*polar(1.0, w_s*time);
	//EZY(mField->getNpx()/2, mField->getNpy()/2) += 0.5*ray_coef*polar(1.0, w_s*time);
	NsScatteredWave(wave_angle);
	//scatteredWave(Ezx, EPS_EZ);
	//scatteredWave(Ezy, EPS_EZ);
	//scatteredWave(Ez, EPS_EZ);

	//absorbing();
	CalcH();

	/*
	if(file){
		file << Ez[index(9*mField->getNx()/10, mField->getNy()/2, +1)] << endl;
	}
	else{		
		file.open("../../Fourie/Fourie/TM_data2.txt");
	}
	*/
	
	if (time > maxStep) {
		MiePrint(Ez, "time" + to_s(maxStep) + "_PML" + to_s(mField->getNpml()) + "_StTM_");
		capture(to_s(time));
		//capture(to_s((int)Inv_Nano_S(lambda_s)));
		return EndTask();
	}
	return true;
}

//----------------終了時の仕事------------------------//
bool StFDTD_TM::EndTask() {
	cout << "End Task" << endl;

	string label = "";
	NTFFindexform("", NTFF::NTFFDATA | NTFF::TOTAL);	// label -> "" にしたけど動くか確認してない.

														//終了条件の確認
	if (!Terminate())
		return false;

	ReStart();
	return true;
}

//----------------吸収境界-----------------------//

void StFDTD_TM::absorbing(){
	absorbing_stRL(Ez,0,	LEFT);
	absorbing_stRL(Ez,mField->getNx()-1,	RIGHT);
	absorbing_stTB(Ez,0,	BOTTOM);
	absorbing_stTB(Ez,mField->getNy()-1,	TOP);
}

//----------------周期境界-----------------------//
void StFDTD_TM::cycle(){
	cycle_stRL(Ez,0   ,LEFT);
	cycle_stRL(Ez,mField->getNx()-1,RIGHT);
	//i=0 にi=Nx-2をコピー, i=Nx-1にi=1をコピー
	for(int i=0; i<mField->getNx(); i++){
		EZ(i,0 ,+1) = EZ(i, mField->getNy()-2, +1);
		EZ(i,mField->getNy()-1, +1) = EZ(i,1, +1);
	}
}

/*
	ButtonFactory::setButton("epsx  ", EPSHX(mField->getNpx() - mField->getNpml()-8 ,mField->getNpy()/2));
	ButtonFactory::setButton("epsy  ", EPSHY(mField->getNpx() - mField->getNpml()-8 ,mField->getNpy()/2));
	ButtonFactory::setButton("CHYLX  ", CHYLX(mField->getNpx() - mField->getNpml()-8 ,mField->getNpy()/2));
	ButtonFactory::setButton("CHXLY  ", CHXLY(mField->getNpx() - mField->getNpml()-8 ,mField->getNpy()/2));
	ButtonFactory::setButton("CEZXLX  ", CEZXLX(mField->getNpx() - mField->getNpml()-8 ,mField->getNpy()/2));
	ButtonFactory::setButton("CEZYLY  ", CEZYLY(mField->getNpx() - mField->getNpml()-8 ,mField->getNpy()/2));
	ButtonFactory::setButton("CHY  ", CHY(mField->getNpx() - mField->getNpml()-8 ,mField->getNpy()/2));
	ButtonFactory::setButton("CHX  ", CHX(mField->getNpx() - mField->getNpml()-8 ,mField->getNpy()/2));
	ButtonFactory::setButton("CEZX  ", CEZX(mField->getNpx() - mField->getNpml()-8 ,mField->getNpy()/2));
	ButtonFactory::setButton("CEZY  ", CEZY(mField->getNpx() - mField->getNpml()-8 ,mField->getNpy()/2));
	ButtonFactory::setButton("Hy(right)", HY(mField->getNpx() - mField->getNpml()-8 ,mField->getNpy()/2).real());
	ButtonFactory::setButton("Hy(top)", HY(mField->getNpx()/2 ,mField->getNpy()- mField->getNpml()).real());
*/