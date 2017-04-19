#ifndef _OBJECT_H
#define _OBJECT_H
#include <stdlib.h>
#include <math.h>
#include <stdio.h> 
#include <string>
#include <algorithm>
#include <fstream>
#include <iostream>
#include <sstream>
#include <glut.h>
#include <vector>
#include <map>
#include <complex>

using namespace std;

#define _USE_MATH_DEFINES


typedef complex<double> Complex;
const int REGION_X = 512;
const int REGION_Y = 512;
const int WINDOW_W = REGION_X + 200;
const int WINDOW_H = REGION_Y;

//メインウィンドウの大きさと位置
const int MAIN_WINDOW_W = 200;
const int MAIN_WINDOW_H = WINDOW_H;
const int MAIN_WINDOW_X = WINDOW_W - MAIN_WINDOW_W;
const int MAIN_WINDOW_Y = 0;

const bool CHECK   = true;		//ボタンに表示した内容をファイルにとっておくかどうか
const bool UNCHECK = false;


const double nano  = 0.000000001;
const double micro = 1000*nano;
const double miri  = 1000000*nano;
const double PI    = 3.14159265358979323846;

const string Fourie = "../../Fourie/Fouie";

class Object{
public:
	virtual int calc() = 0;
	virtual void draw() = 0;
};

struct Color{
	double red;
	double green;
	double blue;
	double alpha;
}typedef color;

template<class T> T _max(T a, T b){
	return a < b ? b:a;
}

template<class T> T _min(T a, T b){
	return b < a ? b:a;
}

template<class T> double _pow(T a, int b){
	if(b==0) return 1.0;				//0乗は1
	if(b<0)  return 1.0/_pow(a, -b);	//a^(-b) = 1/a^b 
	T c =  _pow(a, b/2);
	if(b%2 == 0)		
		return c*c;
	else             
		return c*c*a;
}

template<class T> string to_s(T n){
	string rt;
	stringstream ss;
	ss << n;
	ss >> rt;
	return rt;
}

double s_to_d(const std::string& str);

void save_val(string name, double val);

color lamda_to_color(double lamda);

template<class T> bool between(T &a, T &mi, T &ma){
	return (mi<=a && a<=ma);
}

#endif //_OBJECT_H