#ifndef _SOLVER_H
#define _SOLVER_H

#include "Object.h"
#include "MenuWindow.h"
#include "Model.h"
#include <direct.h>
#include<stack>
#include<math.h>
#include<iomanip>
#include "PhysicalConstant.h"
#include "Field.h"

#define _USE_MATH_DEFINES

enum DIRECT{
	LEFT = +1,
	RIGHT = -1,
	TOP  = -1,
	BOTTOM = +1
};

//NTTFする時の書き出しオプション
namespace NTFF{
	typedef unsigned int output;
	const output NTFFDATA = 1;			//変換したデータ（エネルギー）を書き出す
	const output TOTAL    = 1 << 1;		//変換したデータの合計を最後の行に書き出す
	const output REFLEC   = 1 << 2;		//変換したデータを合計で割った反射率を書き出す
}

namespace DATAFILE{
	typedef bool filemode;
	const filemode    ADD = false;
	const filemode REMOVE = true;
}

template <class T> class Range
{
	T Min;
	T Max;
	T interval;
public:
	Range(T mi, T ma)
	{
		Min = min(mi,ma);
		Max = max(mi,ma);
	}
	
	Range(T mi, T ma, T in)
	{
		Min = min(mi,ma);
		Max = max(mi,ma);
		interval = in;
	}

	Range(T m):Min(m), Max(m)
	{
	}

	Range()
	{
	}

	T MAX(){ return Max;}

	T MIN(){ return Min; }

	T INTERVAL(){ return interval; }
};

class Solver{
private:
	string DataDir;			//データを置くディレクトリへのパス
	string WorkingDir;		//ルートからワーキングディレクトリへのパス
	stack<double> LamList;
	stack<int>	  WaveAngleList;
protected:
	const double H_S, DT_S;	//1セルの大きさのシミュレーション値, 物理量, 1ステップあたりの時間
	double time;	//時間
	double lambda_s, w_s, k_s, T_s;		//シミュレーション用の値
	int		wave_angle;	//波の角度
	double	*n_s;		//屈折率, 誘電率
	double  *Sig_hair;	//メラニン色素の吸収係数
	double ray_coef;
	int maxStep;
	Range<double>	LambdaRange;
	Range<int>		WaveAngleRange;
	FazzyModel	*mModel;
	Field		*mField;	//フィールド
public:
	Solver();
	virtual ~Solver();
	virtual bool calc() = 0;
	virtual void draw() = 0;
	virtual void field() = 0;
	double getTime(){ return time; }
	void nextTime(){
		time += DT_S;							//時間の更新
		ray_coef = 1-exp(-0.0001*time*time);	//波が不連続に入射されるのを防ぐための係数
		if (((int)time) % 100 == 0)		cout << time << endl;
/*
		if ((int)time == 1)
			modelCheck();		//散乱体モデルをcaptureして終了
*/
/*
		if ((int)time == 1)
			modelCheck2();		//途中段階をcapture
*/
/*
		if (((int)time) % 100 == 1) {
			ostringstream oss;
			oss << std::setfill('0') << std::right << std::setw(5) << time - 1;
			string name = oss.str();
			capture(name);		//100stepごとに出力を記録　avi動画作成用
		}
*/
	};

protected:
	void draw(Complex *p);
	void draw(Complex *p, Complex *q);
	void draw_model();
	void modelCheck();
	void modelCheck2();

	//---------------------------------------------------------//
	//--------------------入射波-------------------------------//
	//---------------------------------------------------------//
	void linearLightSource(complex<double> *p);	//横方向の線光源
	void pointLightSource(complex<double> *p);	//中心に点光源
	void scatteredWave(complex<double> *p, double *eps);

	virtual void Initialize(){
		time = 0;	//時間を0に
		DataDir		=  "../DataSet/";
	}

	void SetWaveParameter(double lam){
		lambda_s = lam;				//波長設定
		k_s      = 2*PI/lambda_s;	//波数
		w_s      = LIGHT_SPEED_S*k_s;			//角周波数
		T_s      = 2*M_PI/w_s;		//周期
	}

	double bilinear_interpolation(complex<double> *, double, double);//バイライナーインタポレーション補間

	Color color(double);
	bool neighber(int i, int j);

	bool nextLambda(){
		SetWaveParameter(lambda_s + LambdaRange.INTERVAL());
		if(lambda_s > LambdaRange.MAX())
			return false;
		return true;
	}

	bool nextWaveAngle(){
		wave_angle += WaveAngleRange.INTERVAL();
		if(wave_angle > WaveAngleRange.MAX())
			return false;

		return true;
	}

	bool Terminate(){
//		return false; //今は常に終了

		if( !nextLambda()){
			if(!nextWaveAngle()){
				if(!mModel->update(10)){
					return false;
				}
				else{
					wave_angle = WaveAngleRange.MIN();
					SetWaveParameter(LambdaRange.MIN());
				}
			}
			else{
				SetWaveParameter(LambdaRange.MIN());
			}
		}
		cout << to_s(Inv_Nano_S(lambda_s)) + "nm" << endl;
		cout << to_s(wave_angle) + "deg" << endl;
		return true;
	}

	//吸収境界
	void absorbing_stRL(complex<double> *p, int _X, enum DIRECT a);	//適用配列と, 壁の位置, 左右上下どの壁か
	void absorbing_stTB(complex<double> *p, int _Y, enum DIRECT a);
	void absorbing_nsRL(complex<double> *p, int _X, enum DIRECT a);
	void absorbing_nsTB(complex<double> *p, int _Y, enum DIRECT a);

	//周期境界
	void cycle_stRL(complex<double> *p, int _X, enum DIRECT a);


	void MiePrint(complex<double>*, string);	//Mie散乱解の出力

	void save_data(complex<double> *data, string name);
	void open_data(complex<double> *data, string name);
	void capture(string name);

	// 中心差分の二回微分 x方向	
	complex<double> Dx2(complex<double>*p, int i, int j, int t){
		return p[index(i+1,j, t)] + p[index(i-1,j, t)] - 2.0*p[index(i,j, t)];
	};

	// 中心差分の二回微分 y方向 
	complex<double> Dy2(complex<double>*p, int i, int j, int t){
		return p[index(i,j+1, t)] + p[index(i,j-1, t)] - 2.0*p[index(i,j, t)];
	};

	//中心差分の二回微分, t方向
	complex<double> Dt2(complex<double>* p, int i, int j){
		return p[index(i,j, +1)] + p[index(i,j, -1)] - 2.0*p[index(i,j, 0)];
	};

	//DxDy2 = Dx2Dy
	complex<double> DxDy2(complex<double> *p, int i, int j, int t) {
		return     p[index(i + 1, j - 1, t)] + p[index(i + 1, j + 1, t)] + p[index(i - 1, j - 1, t)] + p[index(i - 1, j + 1, t)]
			- 2.0*(p[index(i + 1, j, t)] + p[index(i - 1, j, t)] + p[index(i, j + 1, t)] + p[index(i, j - 1, t)])
			+ 4.0*p[index(i, j, t)];
	};

	int index(const int& i, const int& j, const int& t){
		int k = ( t + (int)time + 3) % 3;
		return k*mField->getNcel() + index(i,j);
		//return index(i,j);
	};

	int index(const int& i, const int& j){
		return mField->index(i,j);
	};

	//pml用の配列番号取得
	int pmlIndex(const int &i, const int &j, const int &t){
		int k = ( t + (int)time + 3) % 3;
		return k*mField->getNcel() + pmlIndex(i,j);
		//return pmlIndex(i,j); //todo
	}

	int pmlIndex(const int &i, const int &j){
		return mField->pmlIndex(i,j);
	}


	//ゲッター
	double& N_S(const int& i, const int& j){
		return n_s[index(i,j)];
	}

	double& SIG(const int& i, const int& j) {
		return Sig_hair[index(i, j)];
	}

	ofstream WriteOpen(string name, DATAFILE::filemode f=DATAFILE::REMOVE){
		ofstream ofp;
		if(f){
			ofp.open(DataDir + name + ".txt");
		}
		else{
			ofp.open(DataDir + name + ".txt",ios::out | ios::app);
		}
		ofp.setf(ios::fixed, ios::floatfield);	//固定小数点, 20桁制度で指定
		ofp.precision(30);
		return ofp;
	}

	ifstream ReadOpen(string name){
		ifstream ifp(DataDir + name + ".txt");
		if(!ifp){ throw "cannot open file"; }
		return ifp;
	}


	//作成するファイルの名前を波長と、入射角度で返す
	string getDataName(){
		return WorkingDir + getWaveData();
	}

	//波の情報を返す
	string getWaveData(){
		string str = "(WL=" + to_s(Inv_Nano_S(lambda_s)) + "nm,AOI=" + to_s(wave_angle) + "deg)";
		return str;
	}

	string getWorkingDirPass(){
		return WorkingDir;
	}

	void setWorkingDirPass(string wk){
		//WorkingDir = wk;
		DataDir += wk;
		cout << DataDir + WorkingDir << endl;
	}

	string getDataDirPass(){
		return DataDir;
	}

	string MakeDir(string name){
		if(_mkdir((DataDir + name).c_str()) == 0) cout << "makeNEWdir" <<endl;
		return name + "/"; 
	}

	//nano prefix in simulation
	double Nano_S(double a){
		return mField->nanoToCell(a);
	}

	//inverse Nano_S
	double Inv_Nano_S(double a){
		return mField->cellToNano(a);
	}

	//(1 - γ0)/2 を返す
	double NsCoef(){
		//計算用定数の設定
		double kx_s = k_s *          pow(2.0, -0.25);	//2の-4分の1乗
		double ky_s = k_s * sqrt(1 - pow(2.0, -0.5 ));
		double sin2_kx = pow(sin(kx_s*H_S/2), 2);
		double sin2_ky = pow(sin(ky_s*H_S/2), 2);
		double sin2_k  = pow(sin(k_s *H_S/2), 2);
		return (sin2_kx + sin2_ky - sin2_k)/(4*sin2_kx*sin2_ky);		//(1-γ0)/2 = R_M, ( R_P = 1とするとおかしくなる todo?)
	}


	//Maxwell方程式用の係数のひな形 Δt = 1
	//ep_mu εかμ(Eの係数->ε, Hの係数-> μ
	//sig  σ
	double MaxwellCoef(double ep_mu, double sig){
		return (1.0 - sig/ep_mu)/ (1.0 + sig/ep_mu);
	}

	double MaxwellCoef2(double ep_mu, double sig){
		return 1.0/ep_mu / (1.0 + sig/ep_mu);
	}
};

double fazzy_circle(double r, int x0, int y0, int x, int y);


#endif //_SOLVER_H