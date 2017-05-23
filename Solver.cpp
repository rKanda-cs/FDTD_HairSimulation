#include "Solver.h"
#define _CRTDBG_MAP_ALLOC
#include <stdlib.h>
#include <crtdbg.h>
#include <math.h>
#include <stdio.h> 
#include <string>
#include <algorithm>
#include <glut.h>
#include <fstream>
#include <iostream>
#include <omp.h>

#define new ::new(_NORMAL_BLOCK, __FILE__, __LINE__)
#define _USE_MATH_DEFINES

Solver::Solver()
	:H_S(1.0), DT_S(1.0)
{
	mField = new Field(32000, 128000, 50, 10); //width, height, Δh, Npml
	LambdaRange    = Range<double>(Nano_S(380), Nano_S(700), Nano_S(10));
	WaveAngleRange = Range<int>   (135, 135, 30);

	SetWaveParameter( LambdaRange.MIN() );
	wave_angle  = WaveAngleRange.MIN();

	time = 0;
	maxStep  = 8000;
	mField->sig = false;		//吸収係数σの有無　有：true / 無：false (FazzyHair_incidenceModelのみ選択、 その他の場合false)

	n_s     = new double[mField->getNcel()];	//屈折率
	Sig_hair = new double[mField->getNcel()];	//吸光率

	//mModel	= new FazzySlabModel(mField);
	//mModel	= new FazzyMieModel(mField, lambda_s);
	mModel	= new FazzyHair_incidenceModel(mField);
	//mModel	= new FazzyHair_normalModel(mField);
	//mModel	= new FazzyHair_NONcuticleModel(mField);

	DataDir		=  "../DataSet/";
	WorkingDir  =  "";

	cout << "Solver Constructor" << endl;

#ifdef _OPENMP
	cout << "OpenMP used" << endl;
	cout << "The number of processors is " << omp_get_num_procs() << endl;
	cout << "OpenMP : Enabled (Max # of threads = " << omp_get_max_threads() << ")" << endl;
#else
	cout << "OpenMP not used" << endl;
#endif
}

Solver::~Solver(){
	delete[] n_s;
	delete[] Sig_hair;
	delete mModel;
	delete mField;
	cout << "Solver Destructor" << endl;
}

//Bilinear Interpolation補間  実数値の配列番号が引数(四隅からの比で値を出す)
double Solver::bilinear_interpolation(complex<double> *p, double x, double y){
	int i = floor(x);
	int j = floor(y);
	double dx = (x - 1.0*i);	//小数点以下の値
	double dy = (y - 1.0*j);	

	return     norm(p[index(i,    j)]) * (1.0-dx)*(1.0-dy)
	         + norm(p[index(i+1,  j)]) * dx*(1.0-dy)
	         + norm(p[index(i,  j+1)]) * (1.0-dx)*dy
	         + norm(p[index(i+1,j+1)]) * dx*dy;
}

//カラーマップ
Color Solver::color(double phi){
	double range = 2.0;
	Color c;
	double ab_phi = abs(phi);
/*
//	double a = ab_phi < 1 ? (ab_phi <  0.34 ? min(1.0, max(0.0, 3*ab_phi)) : (-1.5*ab_phi+2) ) : 0.5;
	double a = ab_phi < range ? (ab_phi <  range/3.0 ? 3.0/range*ab_phi : (-3.0/4.0/range*ab_phi+1.25) ) : 0.5;
	c.red  = phi > 0 ? a:0;
	c.blue = phi < 0 ? a:0;
	c.green = min(1.0, max(0.0, -3*ab_phi+2));
*/
	//Analyzerと値を合わせる
	double nv = max(0.0, min(1.0, ab_phi));
	// Get color
	if (phi >= 0.75) { c.red = 1.0; c.green = 4.0*(1.0 - nv); c.blue = 0.0; }
	else if (phi >= 0.5) { c.red = 4.0*(nv - 0.5); c.green = 1.0; c.blue = 0.0; }
	else if (phi >= 0.25) { c.red = 0.0; c.green = 1.0; c.blue = 4.0*(0.5 - nv); }
	else { c.red = 0.0; c.green = 4.0*nv; c.blue = 1.0; }
	
	return c;
}

void Solver::MiePrint(complex<double>* p, string name){
	printf("Mie print\n");
	ofstream ofs = WriteOpen(name+"Mie");

	if(ofs) {
		for(int i=0; i<=180; i++){
			double _x = 1.2*lambda_s*cos(i*PI/180) + mField->getNpx()/2;
			double _y = 1.2*lambda_s*sin(i*PI/180) + mField->getNpy()/2;
			double _val = bilinear_interpolation(p,_x,_y);
			ofs << _val << endl;
		}
		printf("output finished\n");
	}
	else{
		printf("No file");
	}
}


//---------------------------------------//
//				 入射波				 	 //
//--------------------------------------//
//点光源の入射
void Solver::pointLightSource(complex<double> *p){
	p[index(mField->getNpx()/2, mField->getNpy()/2, +1)] += ray_coef*polar(1.0, w_s*time);
}

//線光源の入射
void Solver::linearLightSource(complex<double> *p){
	for(int i=1; i<mField->getNy()-1; i++)
		p[index(5,i, +1)]  +=ray_coef*polar(1.0, w_s*time);	//todo 境界上にも入れていいのか?
}

//散乱波の計算
void Solver::scatteredWave(complex<double> *p, double *eps){
	double rad = wave_angle*M_PI/180;	//ラジアン変換
	double _cos = cos(rad), _sin = sin(rad);	//毎回計算すると時間かかりそうだから,代入しておく
	std::complex<double> MinusOne(-1, 0), I;
	I = sqrt(MinusOne);

	for(int i=mField->getNpml(); i<mField->getNpx(); i++){
		for(int j=mField->getNpml(); j<mField->getNpy(); j++){
			if( N_S(i,j) == 1.0 ) continue;		//屈折率が1なら散乱は起きない
			double ikx = k_s*(i*_cos + j*_sin);
			p[index(i,j, +1)] += ray_coef*(1/_pow(N_S(i,j), 2)-1)
				                    *(polar(1.0, ikx-w_s*(time+DT_S))+polar(1.0, ikx-w_s*(time-DT_S))-2.0*polar(1.0, ikx-w_s*time)); 
/*		
			p[index(i, j)] += ray_coef*(EPSILON_0_S / eps[index(i, j)] - 1)*(
				cos(ikx - w_s*(time + 0.5)) + I*sin(ikx - w_s*(time + 0.5))
				- cos(ikx - w_s*(time - 0.5)) - I*sin(ikx - w_s*(time - 0.5))
				);
*/		}
	}
}


//4近傍を調べる
bool Solver::neighber(int _x, int _y){
	for(int i=-2; i<=1; i++)
		if(n_s[index(_min(mField->getNx()-1, _max(0, _x+i%2)), _min(mField->getNy()-1, _max(0,_y+(i+1)%2)) )] != n_s[index(_x,_y)]) return false;

	return true;
}

//-------------------------------------------//
//------------St吸収境界--------------------//
//------------------------------------------//
void Solver::absorbing_stRL(complex<double> *p, int X, enum DIRECT offset){
	double u;
	for(int j=1; j<mField->getNpy()-1; j++){
		u = LIGHT_SPEED_S*DT_S/n_s[index(X,j)];
		if(j == 1 || j == mField->getNpy()-2)		// 四隅の横は一次元吸収境界
			p[index(X,j, +1)] = p[index(X+offset,j, 0)] + (1- u)/(1+u)*(p[index(X,j, 0)] - p[index(X+offset,j, +1)]);

		else						//それ以外は二次元吸収境界
			p[index(X,j, +1)] = - p[index(X+offset,j, -1)] 
								     - (1-u)/(1+u)*(p[index(X,j, -1)] + p[index(X+offset,j, +1)]) 
								     +     2/(1+u)*(p[index(X,j,  0)] + p[index(X+offset,j,  0)]) 
								     + u*u/(1+u)/2*( Dy2(p, X,j, 0)   +  Dy2(p, X+offset,j, 0)	);
												        //  dy^2 φn     +   dy^2 φb
	}
}

void Solver::absorbing_stTB(complex<double> *p, int Y, enum DIRECT offset){
	double u;
	for(int i=1; i<mField->getNpx()-1; i++){
		u = LIGHT_SPEED_S*DT_S/n_s[index(i,Y)];

		if(i==1 || i==mField->getNpx()-2)	//四隅の横は一次元吸収境界
			p[index(i,Y, +1)]    = p[index(i,Y+offset, 0)]    + (1- u)/(1+u)*(p[index(i,Y, 0)]    - p[index(i,Y+offset, +1)]);

		else				//二次元吸収境界
			p[index(i,Y, +1)]  = -p[index(i,Y+offset, -1)]    
								  - (1-u)/(1+u)*(p[index(i,Y, -1)] + p[index(i,Y+offset, +1)]) 
							      +     2/(1+u)*(p[index(i,Y, 0)]  + p[index(i,Y+offset, 0)])     
								  + u*u/(1+u)/2*( Dx2(p, i,Y, 0)   + Dx2(p, i,Y+offset, 0) 	);
												  //  dx^2 φn     +   dx^2 φb
	}		
}

//-----------------------------------------------------------//
//-------------------------Ns吸収境界------------------------//
//----------------------------------------------------------//
/**左右の壁のNS吸収境界
** 適用配列
** 適用するx座標
** 右か左か
*/
void Solver::absorbing_nsRL(complex<double> *p, int X, enum DIRECT offset){
	double kx_s = 1/sqrt(sqrt(2.0)) * k_s;
	double ky_s = sqrt(1 - 1/sqrt(2.0) ) * k_s;
	double u1, u2;
	double w_b  = w_s*DT_S/2;
	double k_b  = k_s*H_S/2;
	double kx_b = kx_s*H_S/2;
	double ky_b = ky_s*H_S/2;

	for (int j = 1; j < mField->getNpy() - 1; j++) {
		u1 = tan(w_b / n_s[index(X, j)]) / tan(k_b);
		u2 = 2 * _pow(sin(w_b / N_S(X, j)), 2) / _pow(sin(ky_b), 2) * (1 - tan(kx_b) / tan(k_b));

		if(j == 1 || j == mField->getNpy()-2)		// 四隅の横は一次元吸収境界
			p[index(X,j, +1)] = p[index(X+offset,j, 0)] + (1- u1)/(1+u1)*(p[index(X,j, 0)] - p[index(X+offset,j, +1)]);

/*		if (j == 1) {
			complex <double> D2_1 = p[index(j + 1, X, 0)] + p[index(mField->getNpy() - 2, X, 0)] - 2.0*p[index(j, X, 0)];
			complex <double> D2_2 = p[index(j + 1, X + offset, 0)] + p[index(mField->getNpy() - 2, X + offset, 0)] - 2.0*p[index(j, X + offset, 0)];
			p[index(j, X, +1)] = -p[index(j, X + offset, -1)]
				- (1 - u1) / (1 + u1)*(p[index(j, X, -1)] + p[index(j, X + offset, +1)])
				+ 2 / (1 + u1)*(p[index(j, X, 0)] + p[index(j, X + offset, 0)])
				+ u2*u2 / (1 + u1) / 2 * (D2_1 + D2_2 + Dt2(p, X, j) + Dt2(p, X + offset, j));
		}
		else if (j == mField->getNpy() - 2) {
			complex <double> D2_1 = p[index(1, X, 0)] + p[index(j - 1, X, 0)] - 2.0*p[index(j, X, 0)];
			complex <double> D2_2 = p[index(1, X + offset, 0)] + p[index(j - 1, X + offset, 0)] - 2.0*p[index(j, X + offset, 0)];
			p[index(j, X, +1)] = -p[index(j, X + offset, -1)]
				- (1 - u1) / (1 + u1)*(p[index(j, X, -1)] + p[index(j, X + offset, +1)])
				+ 2 / (1 + u1)*(p[index(j, X, 0)] + p[index(j, X + offset, 0)])
				+ u2*u2 / (1 + u1) / 2 * (D2_1 + D2_2 + Dt2(p, X, j) + Dt2(p, X + offset, j));
		}
*/		else {						//それ以外は二次元吸収境界
			p[index(X, j, +1)] = -p[index(X + offset, j, -1)]
				- (1 - u1) / (1 + u1)*(p[index(X, j, -1)] + p[index(X + offset, j, +1)])
				+ 2 / (1 + u1)*(p[index(X, j, 0)] + p[index(X + offset, j, 0)])
				+ u2*u2 / (1 + u1) / 2 * (Dy2(p, X, j, 0) + Dy2(p, X + offset, j, 0) + Dt2(p, X, j) + Dt2(p, X + offset, j));
											//  dy^2 φn     +   dy^2 φb
		}
	}
}

/**上下の壁のNS吸収境界
** 適用配列
** 適用するy座標
** 上か下か
*/
void Solver::absorbing_nsTB(complex<double> *p, int Y, enum DIRECT offset){
	double kx_s = 1/sqrt(sqrt(2.0)) * k_s;
	double ky_s = sqrt(1 - 1/sqrt(2.0) ) * k_s;
	double u1, u2;
	double w_b  = w_s*DT_S/2;
	double k_b  = k_s*H_S/2;
	double kx_b = kx_s*H_S/2;
	double ky_b = ky_s*H_S/2;

	for(int i=1; i<mField->getNpx()-1; i++){
		u1 = tan(w_b/n_s[index(i,Y)]) / tan(k_b);
		u2 = 2 * _pow(sin(w_b/n_s[index(i,Y)]), 2) / _pow(sin(ky_b),2) * (1 - tan(kx_b)/tan(k_b));

		if(i==1 || i==mField->getNpx()-2)	//四隅の横は一次元吸収境界
			p[index(i,Y, +1)] = p[index(i,Y+offset, 0)] + (1- u1)/(1+u1)*(p[index(i,Y, 0)] - p[index(i,Y+offset, +1)]);

/*		if (i == 1) {
			complex <double> D2_1 = p[index(i + 1, Y, 0)] + p[index(mField->getNpx() - 2, Y, 0)] - 2.0*p[index(i, Y, 0)];
			complex <double> D2_2 = p[index(i + 1, Y + offset, 0)] + p[index(mField->getNpx() - 2, Y + offset, 0)] - 2.0*p[index(i, Y + offset, 0)];
			p[index(i, Y, +1)] = -p[index(i, Y + offset, -1)]
				- (1 - u1) / (1 + u1)*(p[index(i, Y, -1)] + p[index(i, Y + offset, +1)])
				+ 2 / (1 + u1)*(p[index(i, Y, 0)] + p[index(i, Y + offset, 0)])
				+ u2*u2 / (1 + u1) / 2 * (D2_1 + D2_2 + Dt2(p, i, Y) + Dt2(p, i, Y + offset));
		}
		else if (i == mField->getNpx() - 2) {
			complex <double> D2_1 = p[index(1, Y, 0)] + p[index(i - 1, Y, 0)] - 2.0*p[index(i, Y, 0)];
			complex <double> D2_2 = p[index(1, Y + offset, 0)] + p[index(i - 1, Y + offset, 0)] - 2.0*p[index(i, Y + offset, 0)];
			p[index(i, Y, +1)] = -p[index(i, Y + offset, -1)]
				- (1 - u1) / (1 + u1)*(p[index(i, Y, -1)] + p[index(i, Y + offset, +1)])
				+ 2 / (1 + u1)*(p[index(i, Y, 0)] + p[index(i, Y + offset, 0)])
				+ u2*u2 / (1 + u1) / 2 * (D2_1 + D2_2 + Dt2(p, i, Y) + Dt2(p, i, Y + offset));
		}
*/		else {				//二次元吸収境界
			p[index(i, Y, +1)] = -p[index(i, Y + offset, -1)]
				- (1 - u1) / (1 + u1)*(p[index(i, Y, -1)] + p[index(i, Y + offset, +1)])
				+ 2 / (1 + u1)*(p[index(i, Y, 0)] + p[index(i, Y + offset, 0)])
				+ u2*u2 / (1 + u1) / 2 * (Dx2(p, i, Y, 0) + Dx2(p, i, Y + offset, 0) + Dt2(p, i, Y) + Dt2(p, i, Y + offset));
											//  dx^2 φn     +   dx^2 φb
		}
	}		
}

//----------------------------------------------//
//----------------周期境界-----------------------//
//-----------------------------------------------//

void Solver::cycle_stRL(complex<double> *p, int X, enum DIRECT offset){
	double u;
	for(int i=1; i<mField->getNy()-1; i++){
		u = LIGHT_SPEED_S*DT_S/n_s[index(X,i)];
			//二次元吸収境界
		p[index(X,i, +1)] = - p[index(X+offset,i, -1)]     
							- (1-u)/(1+u)*(p[index(X,i, -1)] + p[index(X+offset,i, +1)])    
							+     2/(1+u)*(p[index(X,i,  0)] + p[index(X+offset,i,  0)]) 
							+ u*u/(1+u)/2*(  Dy2(p,X,i, 0)   +   Dy2(p,X+offset,i, 0) );
												//  dy^2 φn     +   dy^2 φb
	}		
}

//---------------------------描画------------------------------//
//左下が(0,0), 右上が(Nx,Ny)
void Solver::draw(Complex *p, Complex *q){
	double N = max(mField->getNx(),mField->getNy());
	double ws = 2.0*MAIN_WINDOW_X/WINDOW_W/N;
	double hs = 2.0*MAIN_WINDOW_H/WINDOW_H/N;
	for (int i = mField->getNpml(); i < mField->getNpx()-mField->getNpml(); i++){
		for (int j = mField->getNpml(); j < mField->getNpy()-mField->getNpml(); j++){
			int x = i-mField->getNpml();
			int y = j-mField->getNpml();
			Color c = color( norm(p[index(i,j)] + q[index(i,j)]) );
			//Color c = color(30.0*(p[index(i,j)].real() + q[index(i,j)].real()));
			glColor3d(c.red, c.green, c.blue);
			if(j==mField->getNpy()/2 || i==mField->getNpx()/2) glColor3d(1,1,1);
			
			glRectd(x*ws-1, y*hs-1, (x+1.0)*ws-1, (y+1.0)*hs-1);
				
		}
	}

	draw_model();
}

//---------------------------描画------------------------------//
//左下が(0,0), 右上が(Nx,Ny)
void Solver::draw(Complex *p){
	double N = max(mField->getNx(),mField->getNy());
	double ws = 2.0*MAIN_WINDOW_X/WINDOW_W/N;
	double hs = 2.0*MAIN_WINDOW_H/WINDOW_H/N;
	for (int i = mField->getNpml(); i < mField->getNpx()-mField->getNpml(); i++){
		for (int j = mField->getNpml(); j < mField->getNpy()-mField->getNpml(); j++){
			int x = i-mField->getNpml();
			int y = j-mField->getNpml();
			Color c = color( norm(p[index(i,j)]) );
			//Color c = color(30.0*p[index(i,j)].real());
			glColor3d(c.red, c.green, c.blue);
			if(j==mField->getNpy()/2 || i==mField->getNpx()/2) glColor3d(1,1,1);
			
			glRectd(x*ws-1, y*hs-1, (x+1.0)*ws-1, (y+1.0)*hs-1);
				
		}
	}

	draw_model();
}

//散乱体の描画
void Solver::draw_model(){
	double N = max(mField->getNx(),mField->getNy());
	double ws = 2.0*MAIN_WINDOW_X/WINDOW_W/N;
	double hs = 2.0*MAIN_WINDOW_H/WINDOW_H/N;
	for (int i = mField->getNpml(); i < mField->getNpx()-mField->getNpml(); i++){
		for (int j = mField->getNpml(); j < mField->getNpy()-mField->getNpml(); j++){
			int x = i-mField->getNpml();
			int y = j-mField->getNpml();
			//媒質境界
			const double n = N_S(i, j);	//ここで,屈折率を書き換えてはいけない
			const double s = SIG(i, j);

			if(n == 1.0) continue;	//屈折率が1ならとばす
			glEnable(GL_BLEND);
			glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
			glDepthMask(GL_FALSE);
			glColor4d(0.7/(n+s), 0.7/(n+s), 0.7/(n+s), 0.6);
			glDepthMask(GL_TRUE);
		glRectd(x*ws-1, y*hs-1, (x+1.0)*ws-1, (y+1.0)*hs-1);	
		}
	}

}

void Solver::modelCheck() {
	capture("../../model");
	exit(-1);
}

//-----------------データの保存----------------------//
void Solver::save_data(complex<double> *data, string name){
	ofstream out = WriteOpen(name);

	for(int k=-1; k<2; k++)
		for(int i=0; i<mField->getNx(); i++)
			for(int j=0; j<mField->getNy(); j++)
				out << data[index(i,j, k)] << endl;

}

void Solver::open_data(complex<double> *data, string name){
	ifstream in = ReadOpen(name);

	for(int k=-1; k<2; k++)
		for(int i=0; i<mField->getNx(); i++)
			for(int j=0; j<mField->getNy(); j++){
				if(!in) {throw  "cannot open file" ;}
				in >> data[index(i,j, k)];
			}
}


//OpenCV関係
#include <opencv2/opencv.hpp> // インクルードファイル指定
#include <opencv2/highgui/highgui.hpp>
//#include <opencv2/opencv_lib.hpp> // 静的リンクライブラリの指定
#include <opencv2/core/core.hpp>
#pragma comment(lib,"opencv_highgui2411.lib")
#pragma comment(lib,"opencv_core2411.lib")

// @brief 現在の画面の状態をキャプチャしてpngに保存する //

void Solver::capture(string name)
{
	int width = WINDOW_W, height = WINDOW_H;

	cv::Mat cvmtx(cv::Size(width, height), CV_8UC4, cv::Scalar(0, 0, 0));//黒で初期化
																		 // 画像のキャプチャ
	glReadBuffer(GL_FRONT);// フロントを読み込む様に設定する
	glPixelStorei(GL_UNPACK_ALIGNMENT, 1); // 初期値は4
	glReadPixels(0, 0, width, height, GL_BGRA_EXT, GL_UNSIGNED_BYTE, (void*)cvmtx.data);
	//上下逆にする
	cv::flip(cvmtx, cvmtx, 0);
	// 画像の書き出し 
	cv::imwrite(DataDir + name + ".jpg", cvmtx);
}

//---------------------------------------------//
//--------------PML----------------------------//
//--------------------------------------------//

/*
Color Solver::color2(double a){
	Color c;
	c.red  = max(0.0, min(1.0, a-1.0));
	c.blue = max(0.0, min(-a+1, 1.0));
	c.green = max(0.0, min(1.0, -2*abs(a-1)+1));
	return c;
}

*/