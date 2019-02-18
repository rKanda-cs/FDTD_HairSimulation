#ifndef _MODEL_H
#define _MODEL_H
#include"PhysicalConstant.h"
#include<string>
#include<direct.h>
#include"Object.h"
#include"Field.h"

using namespace std;
enum INTEG{
	D_XY = 0,
	D_X,
	D_Y
};

class FazzyModel{
protected:
	Field *mField;
	
public:
	FazzyModel(Field *field){
		mField = field;
	};
	virtual string mkdir(string root) = 0;
	virtual double calcEPS(const double&, const double&, enum INTEG = D_XY) = 0;
	virtual double calcSIG(const double&, const double&, const double lam, enum INTEG = D_XY) = 0;
	virtual bool update(int) = 0;
	virtual void Initialize()=0;
};

class FazzySlabModel :public FazzyModel{
	const double ep1, ep2;
	const int width1, width2, layer;
	bool pain;
	double lay[7];	// pain = trueのとき、剥がれの角度  pain = falseのとき、層の有無(0or1)
	int random[7];
public:
	FazzySlabModel(Field*);
	string mkdir(string root);
	double calcEPS(const double&, const double&, enum INTEG = D_XY);
	double calcSIG(const double&, const double&, const double lam, enum INTEG = D_XY) {
		return 0;
	}
	bool update(int){
		//return true;
		return false;
	}
	void Initialize()
	{
	}
};

class FazzyMieModel :public FazzyModel{
	double ep;
	double r;
public:
	FazzyMieModel(Field* f, double _r);
	string mkdir(string root);	//ディレクトリ作成

	//誘電率計算
	double calcEPS(const double& x, const double& y, enum INTEG f);
	double calcSIG(const double&, const double&, const double lam, enum INTEG = D_XY) {
		return 0;
	}
	bool update(int a){
		  return true;
	}
	void Initialize()
	{
	}
};

class FazzyHair_incidenceModel :public FazzyModel {
	const double ep1, ep2;
	const int alpha;
	double alphaR, length, ln, lx, ly, rn, cwidth, r;
public:
	FazzyHair_incidenceModel(Field*);
	string mkdir(string root);
	double calcEPS(const double&, const double&, enum INTEG = D_XY);
	double calcSIG(const double&, const double&, const double lam, enum INTEG = D_XY);
	bool update(int) {
		//return true;
		return false;
	}
	void Initialize()
	{
	}
};

class FazzyHair_incidenceLayerModel :public FazzyModel {
	const double ep1, ep2, ep3;
	const int alpha;
	double alphaR, r, rn, cwidth, cn, cmc, mn, length, ln, lx, ly;
	/*
	alpha:傾き		alphaR:傾き(ラジアン)
	r:毛皮質範囲の半径(μm)							rn:毛皮質範囲の半径(nmシミュレーション値)
	cwidth:キューティクル厚さ(μm)					cn:キューティクル厚さ(nmシミュレーション値)
	cmc:CMC幅(μm)									mn:CMC幅(nmシミュレーション値)
	length:キューティクル長さ(μm)					ln:キューティクル長さ(nmシミュレーション値)
	ly:キューティクル範囲(nmシミュレーション値)		lx:x方向長さ(nmシミュレーション値)
	*/
public:
	FazzyHair_incidenceLayerModel(Field*);
	string mkdir(string root);
	double calcEPS(const double&, const double&, enum INTEG = D_XY);
	double calcSIG(const double&, const double&, const double lam, enum INTEG = D_XY);
	bool update(int) {
		//return true;
		return false;
	}
	void Initialize()
	{
	}
};

class FazzyHair_normalModel :public FazzyModel {
	const double ep1, ep2, e;
	const int r;
	double rn, ax, by;
public:
	FazzyHair_normalModel(Field*);
	string mkdir(string root);
	double calcEPS(const double&, const double&, enum INTEG f);
	double calcSIG(const double&, const double&, const double lam, enum INTEG = D_XY) {
		return 0;
	}
	bool update(int) {
		//return true;
		return false;
	}
	void Initialize()
	{
	}
};

class FazzyHair_NONcuticleModel :public FazzyModel {
	const double ep1, ep2;
	const int r;
	double rn;
public:
	FazzyHair_NONcuticleModel(Field*);
	string mkdir(string root);
	double calcEPS(const double&, const double&, enum INTEG = D_XY);
	double calcSIG(const double&, const double&, const double lam, enum INTEG = D_XY);
	bool update(int) {
		//return true;
		return false;
	}
	void Initialize()
	{
	}
};

class FazzyMorphoModel :public FazzyModel {
	const double ep1, ep2;
	const int X = 100;
	const int Y = 400;
	bool p[100][400];

public:
	FazzyMorphoModel(Field*);
	string mkdir(string root);
	double calcEPS(const double&, const double&, enum INTEG = D_XY);
	double calcSIG(const double&, const double&, const double lam, enum INTEG = D_XY) {
		return 0;
	}
	bool update(int) {
		//return true;
		return false;
	}
	void Initialize()
	{
	}
};



#endif //_MODEL_H