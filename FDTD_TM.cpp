#define _USE_MATH_DEFINES
#define _CRTDBG_MAP_ALLOC
#include <stdlib.h>
#include <crtdbg.h>
#include "Object.h"
#include "FDTD_TM.h"
#include "Vec.h"
#include <math.h>
#include <cstring>
#define new ::new(_NORMAL_BLOCK, __FILE__, __LINE__)

using namespace std;
const int n_ang = 180;
FDTD_TM::FDTD_TM()
:Solver()
{
	//領域確保
	Ez  = new complex<double>[3*mField->getNcel()];		//Ez(i,j)      → Ez(i,j)
	Hx  = new complex<double>[3*mField->getNcel()];		//Hx(i, j+0.5) → Hx(i,j)
	Hy  = new complex<double>[3*mField->getNcel()];		//Hy(i+0.5, j) → Hy(i,j) を意味する
	Ezx = new complex<double>[3*mField->getNcel()];
	Ezy = new complex<double>[3*mField->getNcel()];

	//計算用定数配列
	C_EZ     = new double[mField->getNcel()];	//Cez(i, j)       → CEZ(i,j)
	C_EZLH	 = new double[mField->getNcel()];	//Cezlh(i, j)     → CEZLH(i,j)
	C_HXLY   = new double[mField->getNcel()];	//Chxly(i, j+0.5) → CHXLY(i,j)
	C_HYLX   = new double[mField->getNcel()];	//Chylx(i+0.5, j) → CHYLX(i,j) を意味する

	C_EZX   = new double[mField->getNcel()]; 
	C_EZY   = new double[mField->getNcel()];
	C_EZXLX = new double[mField->getNcel()];
	C_EZYLY = new double[mField->getNcel()];
	C_HX    = new double[mField->getNcel()];
	C_HY    = new double[mField->getNcel()];
	C_HXLY  = new double[mField->getNcel()];	//Chxly(i, j+0.5) → CHXLY(i,j)
	C_HYLX  = new double[mField->getNcel()];	//Chylx(i+0.5, j) → CHYLX(i,j) を意味する

	EPS_EZ = new double[mField->getNcel()];
	EPS_HX = new double[mField->getNcel()];
	EPS_HY = new double[mField->getNcel()];

	B_EZXp = new double[mField->getNcel()];
	B_EZXm = new double[mField->getNcel()];
	B_EZYp = new double[mField->getNcel()];
	B_EZYm = new double[mField->getNcel()];
	B_HXp = new double[mField->getNcel()];
	B_HXm = new double[mField->getNcel()];
	B_HYp = new double[mField->getNcel()];
	B_HYm = new double[mField->getNcel()];

	//領域初期化
	for(int i=0; i<3*mField->getNcel(); i++)
			Ez[i] = Hx[i] = Hy[i] = Ezx[i] = Ezy[i] = 0;

	cout << "FDTD_TM Constructor" << endl;
}

FDTD_TM::~FDTD_TM(){
	delete[] Ez;
	delete[] Hx;
	delete[] Hy;
	delete[] Ezx;
	delete[] Ezy;
	delete[] C_HX;
	delete[] C_HY;
	delete[] C_EZX;
	delete[] C_EZY;
	delete[] C_EZXLX;
	delete[] C_EZYLY;
	delete[] C_HXLY;
	delete[] C_HYLX;
	
	delete[] C_EZ;
	delete[] C_EZLH;

	delete[] EPS_EZ;
	delete[] EPS_HX;
	delete[] EPS_HY;

	delete[] B_EZXp;
	delete[] B_EZXm;
	delete[] B_EZYp;
	delete[] B_EZYm;
	delete[] B_HXp;
	delete[] B_HXm;
	delete[] B_HYp;
	delete[] B_HYm;

	cout << "FDTD_TM Destructor" << endl;
}

void FDTD_TM::Initialize(){
	super::Initialize();
	//領域初期化
	for(int i=0; i<3*mField->getNcel(); i++)
			Ez[i] = Hx[i] = Hy[i] = Ezx[i] = Ezy[i] = 0;
}

void FDTD_TM::OpenData(string prefix){
	string str = prefix+ "data/" + getWaveData();	//データオープン
	cout << str << endl;

	open_data(Ez, str +"Ez");
	open_data(Hx, str +"Hx");
	open_data(Hy, str +"Hy");
}

void FDTD_TM::SaveData(string prefix){
	string str = MakeDir(prefix+ "data");
	cout << str << endl;
	save_data(Ez, str + "/" + getWaveData() + "Ez");
	save_data(Hx, str + "/" + getWaveData() + "Hx");
	save_data(Hy, str + "/" + getWaveData() + "Hy");
}

void FDTD_TM::field(){
	setWorkingDirPass(mModel->mkdir(getDataDirPass()));	//データを保存するディレクトリの設定
	setWorkingDirPass(MakeDir("TM"));
	for(int i=0; i<mField->getNpx(); i++){
		for(int j=0; j<mField->getNpy(); j++){
			EPSEZ(i,j) = mModel->calcEPS(i,j);
			EPSHX(i,j) = mModel->calcEPS(i,j+0.5);
			EPSHY(i,j) = mModel->calcEPS(i+0.5,j);
			N_S(i,j)   = sqrt( mModel->calcEPS(i+0.5,j+0.5) / EPSILON_0_S );	//屈折率 n = √(ε/ ε0) (ここでは μ = μ0としているので)

			if(mField->sig == true)
				SIG(i, j) = mModel->calcSIG(i, j, lambda_s);	// Hair/incidenceモデルのみ(メラニン色素考慮有)
		}
	}
}

//散乱波
void FDTD_TM::NsScatteredWave(int ang){
	double rad = ang*M_PI/180;	//ラジアン変換
	double _cos = cos(rad)*k_s;	//毎回計算すると時間かかりそうだから,代入しておく
	double _sin = sin(rad)*k_s;	
	//double a = (1-exp(-_pow(0.01*time,2)));		//不連続に入射されるのを防ぐ為の係数
	for(int i=mField->getNpml(); i<mField->getNpx()-mField->getNpml(); i++){
		for(int j=mField->getNpml(); j<mField->getNpy()-mField->getNpml(); j++){
			if(EPSEZ(i,j) == EPSILON_0_S) continue;
			double n  = sqrt(EPSEZ(i,j)/EPSILON_0_S);				//todo σ = 0が前提の話
			double u0 = sin(w_s/n*H_S/2.0)  / sin(  k_s*H_S/2.0);
			double u1 = sin(w_s/n*DT_S/2.0) / sin(n*k_s*H_S/2.0);
			double _n = u0/u1;
			double ikx = i*_cos + j*_sin;		//k_s*(i*cos + j*sin)
			EZX(i,j,+1) += 0.5*ray_coef*(1.0/(_n*n) - 1)*(polar(1.0, ikx - w_s*(time+DT_S)) - polar(1.0, ikx - w_s*time));
			EZY(i,j,+1) += 0.5*ray_coef*(1.0/(_n*n) - 1)*(polar(1.0, ikx - w_s*(time+DT_S)) - polar(1.0, ikx - w_s*time)); 
			EZ(i, j,+1) += ray_coef*(1.0/(_n*n) - 1)*(polar(1.0, ikx - w_s*(time+DT_S)) - polar(1.0, ikx - w_s*time));
		}
	}
}

//入射波
void FDTD_TM::IncidentWave(int ang){
	double rad = ang*M_PI/180;	//ラジアン変換
	double _cos = cos(rad), _sin = sin(rad);	//毎回計算すると時間かかりそうだから,代入しておく
	//double a = (1-exp(-_pow(0.01*time,2)));		//不連続に入射されるのを防ぐ為の係数
	double a = ray_coef;
	for(int i=0; i<mField->getNx()-1; i++){
		for(int j=0; j<mField->getNy()-1; j++){
			double ikx = k_s*(i*_cos + j*_sin);
			EZX(i, j, +1) += 0.5*polar(1.0, ikx - w_s*(time + DT_S));
			EZY(i, j, +1) += 0.5*polar(1.0, ikx - w_s*(time + DT_S));
			EZ(i, j, +1) += polar(1.0, ikx - w_s*(time+DT_S));
		}
	}
}

void FDTD_TM::IncidentWaveH(int ang){
	double rad = ang*M_PI/180;	//ラジアン変換
	double _cos = cos(rad), _sin = sin(rad);	//毎回計算すると時間かかりそうだから,代入しておく
	//double a = (1-exp(-_pow(0.01*time,2)));		//不連続に入射されるのを防ぐ為の係数
		double a = ray_coef;
	for(int i=0; i<mField->getNx()-1; i++){
		for(int j=0; j<mField->getNy()-1; j++){
			double ikx = k_s*(i*_cos + j*_sin);
			HY(i,j, +1) += - k_s/(w_s*MU_0_S) * polar(1.0, ikx - w_s*(time+DT_S) );
			HX(i,j, +1) += + k_s/(w_s*MU_0_S) * polar(1.0, ikx - w_s*(time+DT_S) );
		}
	}
}

void FDTD_TM::draw(){
	//super::draw(Ezx, Ezy);	
	super::draw(Ez);	
}

bool FDTD_TM::calc(){
	return true;
}

void FDTD_TM::PML(){

}

//-------------------Near_To_Far変換---------------------//
void FDTD_TM::NTFFindexform(string label, NTFF::output flag){
	cout << "NTF start" << endl;

	int cx = mField->getNpx()/2;	//座標変換後の原点
	int cy = mField->getNpy()/2;
	
	double r0 = 1.0e6;
    complex<double> iu( 0.0, 1.0 );		//単位虚数ベクトル
	complex<double> Coef = sqrt( iu*k_s/(8*M_PI*r0) ) * exp( iu*k_s*r0 );	//√(jw/8πcr)e^(jkr)
    int offset = 5;		// closed line offset

	int lt,rt, tp, bm;		//閉曲面の場所	//計算領域内
	tp = mField->getNpy() - mField->getNpml() - offset;			//上から-5
	bm = mField->getNpml() + offset;			//下から5
	rt = mField->getNpx() - mField->getNpml() - offset;		//右から-5
	lt = mField->getNpml() + offset;			//左から5

	double sum = 0;			//総和
	double strength[360];
	int max_angle = 360;
	for(int ang=0; ang<max_angle; ang++){
		//double rad = ang*M_PI/180.0;
		double rad = (ang-90)*M_PI / 180.0;
		Vec2<double> r(cos(rad), sin(rad));		//遠方の方向ベクトル.
		Vec2<double> r2;						//中心からセルまでの距離

		complex<double> Nz( 0, 0 );
        complex<double> Lx( 0, 0 );
        complex<double> Ly( 0, 0 );
        complex<double> C_EZ, C_HX, C_HY;
		
        // (left,bottom) -> (right,bottom) n=(0, -1)
        for ( int i=lt; i<rt; i++ ) {
            r2   = Vec2<double>( i-cx, bm-cy);	//中心からセルまでの距離
			C_EZ =  EZ(i,bm, +1);
			C_HX = 0.5 * ( HX(i,bm, +1) + HX(i,bm-1, +1) ); 

            Nz  += C_HX * exp( iu * k_s * In_prod(r,r2) );
            Lx  += C_EZ * exp( iu * k_s * In_prod(r,r2) );
        }

        // (right,bottom) -> (right,top) n=(1,0)
        for ( int j=bm; j<tp; j++ ) {
            r2   = Vec2<double>( rt-cx, j-cy);
			C_EZ = EZ(rt,j, +1);
			C_HY = 0.5 * ( HY(rt,j, +1) + HY(rt-1,j, +1) );

            Nz  += C_HY * exp( iu * k_s * In_prod(r,r2) );
            Ly  += C_EZ * exp( iu * k_s * In_prod(r,r2) );	
        }

        // (right,top) -> (left,top)  n=(0,1)
        for ( int i=lt; i<rt; i++ ) {
            r2    = Vec2<double>( i-cx, tp-cy );
			C_EZ = EZ(i,tp, +1);
			C_HX = 0.5 * ( HX(i,tp, +1) + HX(i,tp-1, +1) );
            Nz   -= C_HX * exp( iu * k_s * In_prod(r,r2) );		
            Lx   -= C_EZ * exp( iu * k_s * In_prod(r,r2) );
        }

        // (left,top) -> (left,bottom)
        for ( int j=bm; j<tp; j++ ) {
            r2    = Vec2<double>( lt-cx, j-cy);
            C_EZ  = EZ(lt,j, +1);
            C_HY  = 0.5 * ( HY(lt,j, +1) + HY(lt-1,j, +1));
            Nz   -= C_HY * exp( iu * k_s * In_prod(r,r2) );		
            Ly   -= C_EZ * exp( iu * k_s * In_prod(r,r2) );
        }
 
        // Get Ephi
        complex<double> EZZ = Coef * ( Z_0_S*Nz - Lx*sin(rad) + Ly*cos(rad) ); //岡田先輩の実装では Lx, Lyに-がつく
        
		strength[int(ang)] = norm(EZZ);
		sum += norm(EZZ);
	}

	//NTFFの結果を出力
	if( (flag & NTFF::NTFFDATA) == NTFF::NTFFDATA ){
		ofstream fp = WriteOpen(MakeDir("NTFF") + label + getWaveData());
		for(int i=0; i<max_angle; i++)
			fp << strength[i] << endl;
	}

	//NTFFの結果の総和を出力
	if( (flag & NTFF::TOTAL) == NTFF::TOTAL){
		ofstream fp= WriteOpen(MakeDir("NTFF") + label + "WaveAngleStrength", DATAFILE::ADD); //falseは追記モード
		fp << "(" << to_s((int)Inv_Nano_S(lambda_s)) << "," << wave_angle << ")  " << sum << endl;		//波長ごとの総和を出力
	}

	//反射率を出力
	if( (flag & NTFF::REFLEC) == NTFF::REFLEC){
		ofstream fp = WriteOpen(MakeDir("Reflection") + label + getWaveData());
		for(int ang = 0; ang < 360; ang++)
			fp << strength[ang] / sum << endl;
	}

	cout << "NTF finish" << endl;
}
