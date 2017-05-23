#include"Model.h"
#include"Field.h"

/*---------------------------------------------*/
/*--------------円柱Mie散乱--------------------*/
/*---------------------------------------------*/
FazzyMieModel::FazzyMieModel(Field *f, double _r):
FazzyModel(f),r(_r)
{
	ep = 1.6*1.6*EPSILON_0_S;			//誘電率 = (屈折率)^2
		cout << "r=" + to_s((int)mField->cellToNano(r)) << endl;
}

string FazzyMieModel::mkdir(string root){
	_mkdir((root + "Mie").c_str());

	string name = "Mie/" + to_s((int)(mField->cellToNano(r))) +"nm,"+ mField->getStringCellInfo();
	_mkdir((root + name).c_str());	//ディレクトリの作成
	return name + "/";
}

double FazzyMieModel::calcEPS(const double& x, const double& y, enum INTEG f){

	double mx = x - mField->getNpml(); //計算領域内へ写像
	double my = y - mField->getNpml();
	if (mx < 0 || my < 0 || mx >= mField->getNx() || my >= mField->getNy()) return EPSILON_0_S;

	double _x = mx - 0.5*mField->getNx();//(N_X/2, N_Y/2)を原点にシフト
	double _y = my - 0.5*mField->getNy();

	//中心との距離が半径+√2/2セル以上なら, 完全に媒質の外(念のため, 半径+1 以上か調べている)
	if(_x*_x + _y*_y >= pow(r+1, 2.0))

		return EPSILON_0_S;

	//中心との距離が, 半径-√2/2セル以内なら, 完全に媒質の外
	if(_x*_x + _y*_y <= pow(r-1, 2.0))
		return ep;

	double s=0;

	double a = 1.0,b=1.0;
	if(f == D_X) b = 0;
	if(f == D_Y) a = 0;
	for(double i=-16+0.5; i<16; i+=1)
		for(double j=-16+0.5; j<16; j+=1)
			if(pow(_x+a*i/32.0, 2.0) + pow(_y+b*j/32.0, 2.0) <= r*r)
				s+=1; 
	s /= 32.0*32.0;
	return ep*s + EPSILON_0_S*(1-s);
}

/*---------------------------------------------*/
/*--------------多層膜-------------------------*/
/*---------------------------------------------*/
FazzySlabModel::FazzySlabModel(Field* f):
FazzyModel(f), ep1(2.0*2.0*EPSILON_0_S), ep2(EPSILON_0_S), width1(250), width2(50)
{
}

double FazzySlabModel::calcEPS(const double& x, const double& y, enum INTEG f){
//左100nmから,250nm間隔で50nmのスラブを入れていく  **左250nmから(L70.71)10nmスラブに変更(L73)
//多層膜
	
	double mx = x - mField->getNpml(); //計算領域内へ写像
	double my = y - mField->getNpml();

	if(mx < 0 || my < 0 || mx >= mField->getNx() || my >= mField->getNy() ) return EPSILON_0_S;

	int k    = (int)(mField->cellToNano(mx) - 250)%250;
	double l =      (mField->cellToNano(mx) - 250)/250;

	if( k > 0 && k <=10 && l < 5)
		return ep1;
	else
		return ep2;

}

string FazzySlabModel::mkdir(string root){
	_mkdir((root + "SlabModel").c_str());

	string name = "SlabModel/" + mField->getStringCellInfo();
	_mkdir((root + name).c_str());	//ディレクトリの作成
	return name + "/";
}

/*---------------------------------------------*/
/*---------------毛髪--------------------------*/
/*---------------------------------------------*/

/*----------------------------*/
/*-----------縦断面-----------*/
/*----------------------------*/
FazzyHair_incidenceModel::FazzyHair_incidenceModel(Field* f) :
	FazzyModel(f), ep1(1.55*1.55*EPSILON_0_S), ep2(EPSILON_0_S), alpha(3), cwidth(1), r(32+(8-cwidth))
	//alpha:キューティクルの角度(deg)  length:キューティクルの幅(μm)  r:毛の半径(μm)(半径+キューティクルが重なる領域)
{
	alphaR = alpha * PI / 180;
	length = cwidth / sin(alphaR);
	cout << "キューティクルの角度 : " + to_s(alpha) + "deg" << endl;
	cout << "キューティクル幅 : " + to_s(cwidth) + "micro" << endl;
	cout << "キューティクル1枚の露出幅 : " + to_s(length) + "micro" << endl;
}

double FazzyHair_incidenceModel::calcEPS(const double& x, const double& y, enum INTEG f) {
	alphaR = alpha * PI / 180;
	ln = mField->nanoToCell(length * 1000);
	lx = ln * cos(alphaR);
	rn = mField->nanoToCell(r * 1000);
	
	double mx = x - mField->getNpml(); //計算領域内へ写像
	double my = y - mField->getNpml();
	double cx = mField->getNx() / 2;
	double cy = mField->getNy() / 2;
	
	double h = mField->nanoToCell(0*1000);

	if (mx < 0 || my < 0 || mx >= mField->getNx() || my >= mField->getNy()) return ep2;	//PML層
	if (my < cy)	my = 2 * cy - my;		//x軸に対して線対称


	/***************************************************/
	//上面だけキューティクルなし
//	if (y - mField->getNpml() >= mField->nanoToCell(32 * 1000) + cy)	return ep2;
	/***************************************************/


	int c = mField->getNx() / lx + 1;		//計算範囲内のキューティクルの数
	for (int i = 0; i < c; i++) {
		if (mx > i * lx + h && mx < (i + 1) * lx + h && mx < mField->getNx() - h) {
			//			if (my > tan(alphaR) * (mx - lx*i) + cy + rn)	return ep2;
			//			else return ep1;		//Fuzzyなし(Staircaseモデル)

			double dy1 = my - (tan(alphaR) * (mx - lx*i - h) + cy + rn);
			double dy2 = my - (tan(alphaR) * ((mx - lx*i - h) + 1) + cy + rn);
			double s;
			if (dy1 > 0 && dy2 > 0) return ep2;		//キューティクル直線の外側 (1)
			if (fabs(dy1) > 1 && fabs(dy2) > 1) return ep1;		//キューティクル直線の内側 (2)

			if (dy1 <= 0 && dy2 <= 0) {
				if (fabs(dy1) <= 1 && fabs(dy2) <= 1) {
					s = (fabs(dy1) + fabs(dy2)) * 1.0 / 2.0;
					return ep1 * s + ep2 * (1 - s);		// (3)
				}
				if (fabs(dy1) < 1 && fabs(dy2) > 1) {
					s = (1 - fabs(dy1)) * ((my - cy - rn) / tan(alphaR) - (mx - lx*i - h)) / 2;
					return ep2 * s + ep1 * (1 - s);		// (4)
				}
			}
			if (dy1 > 0 && dy2 < 0) {
				s = fabs(dy2) * (((mx - lx*i - h) + 1) - (my - cy - rn) / tan(alphaR)) / 2;
				return ep1 * s + ep2 * (1 - s);		// (5)
			}
		}
		else
			continue;
			//break;
	}

	return ep2;
}

double FazzyHair_incidenceModel::calcSIG(const double& x, const double& y, const double lam, enum INTEG f) {
	rn = mField->nanoToCell(32 * 1000);

	double mx = x - mField->getNpml(); //計算領域内へ写像
	double my = y - mField->getNpml();
	double cx = mField->getNx() / 2;
	double cy = mField->getNy() / 2;

	double h = mField->nanoToCell(0 * 1000);

	if (mx < 0 || my < 0 || mx >= mField->getNx() || my >= mField->getNy()) return 0;	//PML層
	if (my < cy)	my = 2 * cy - my;		//x軸に対して線対称

	if (my <= rn + cy) {
		int k = (int)(mField->cellToNano(my) - mField->cellToNano(cy) - 1500) % 4000;
		double l = (mField->cellToNano(my) - mField->cellToNano(cy) - 1500) / 4000;

		if (k > 0 && k <= 1000 && l < 8)	return 1.0;		//黒色の色素
		else	return 0;
	}
	else  return 0;
}

string FazzyHair_incidenceModel::mkdir(string root) {
	_mkdir((root + "HairModel").c_str());
	string name;

	if (mField->sig == false) {
		_mkdir((root + "HairModel/incidenceplane").c_str());				//吸収係数なしの場合
		name = "HairModel/incidenceplane/" + mField->getStringCellInfo();
	}
	else if (mField->sig == true) {
		_mkdir((root + "HairModel/incidenceplane_withSig").c_str());		//吸収係数ありの場合
		name = "HairModel/incidenceplane_withSig/" + mField->getStringCellInfo();
	}

	_mkdir((root + name).c_str());	//ディレクトリの作成
	
	return name + "/";
}

/*------------------------------------------------*/
/*-----------縦断面(キューティクルなし)-----------*/
/*------------------------------------------------*/
FazzyHair_NONcuticleModel::FazzyHair_NONcuticleModel(Field* f) :
	FazzyModel(f), ep1(1.55*1.55*EPSILON_0_S), ep2(EPSILON_0_S), r(32)
	//r:毛の半径(μm)
{
	cout << "キューティクル : なし"<< endl;
}

double FazzyHair_NONcuticleModel::calcEPS(const double& x, const double& y, enum INTEG f) {
	rn = mField->nanoToCell(r * 1000);

	double mx = x - mField->getNpml(); //計算領域内へ写像
	double my = y - mField->getNpml();
	double cx = mField->getNx() / 2;
	double cy = mField->getNy() / 2;

	double h = mField->nanoToCell(0 * 1000);

	if (mx < 0 || my < 0 || mx >= mField->getNx() || my >= mField->getNy()) return ep2;	//PML層
	if (my < cy)	my = 2 * cy - my;		//x軸に対して線対称

	if (my <= rn + cy)	return ep1;
	else  return ep2;
}

double FazzyHair_NONcuticleModel::calcSIG(const double& x, const double& y, const double lam, enum INTEG f) {
	rn = mField->nanoToCell(r * 1000);

	double mx = x - mField->getNpml(); //計算領域内へ写像
	double my = y - mField->getNpml();
	double cx = mField->getNx() / 2;
	double cy = mField->getNy() / 2;

	double h = mField->nanoToCell(0 * 1000);

	if (mx < 0 || my < 0 || mx >= mField->getNx() || my >= mField->getNy()) return 0;	//PML層
	if (my < cy)	my = 2 * cy - my;		//x軸に対して線対称

	if (my <= rn + cy) {
		int k = (int)(mField->cellToNano(my) - mField->cellToNano(cy) - 1500) % 4000;
		double l = (mField->cellToNano(my) - mField->cellToNano(cy) - 1500) / 4000;

		if (k > 0 && k <= 1000 && l < 8)	return 1.0;		//黒色の色素
		else	return 0;
	}
	else  return 0;
}

string FazzyHair_NONcuticleModel::mkdir(string root) {
	_mkdir((root + "HairModel").c_str());
	string name;

	if (mField->sig == false) {
		_mkdir((root + "HairModel/NONcuticle").c_str());				//吸収係数なしの場合
		name = "HairModel/NONcuticle/" + mField->getStringCellInfo();
	}
	else if (mField->sig == true) {
		_mkdir((root + "HairModel/NONcuticle_withSig").c_str());		//吸収係数ありの場合
		name = "HairModel/NONcuticle_withSig/" + mField->getStringCellInfo();
	}

	_mkdir((root + name).c_str());	//ディレクトリの作成
	return name + "/";
}

/*----------------------------*/
/*-----------横断面-----------*/
/*----------------------------*/
FazzyHair_normalModel::FazzyHair_normalModel(Field* f) :
	FazzyModel(f), ep1(1.55*1.55*EPSILON_0_S), ep2(EPSILON_0_S), e(0.6), r(32)
	//a:離心率  r:毛の半径(μm)
{
	cout << "楕円の離心率 = " + to_s((double)e) << endl;
}

double FazzyHair_normalModel::calcEPS(const double& x, const double& y, enum INTEG f) {
	rn = mField->nanoToCell(r * 1000);
	ax = rn;
	by = ax * sqrt(1 - e*e);

	double mx = x - mField->getNpml(); //計算領域内へ写像
	double my = y - mField->getNpml();
	double cx = mField->getNx() / 2;
	double cy = mField->getNy() / 2;

	if (mx < 0 || my < 0 || mx >= mField->getNx() || my >= mField->getNy()) return ep2;	//PML層

	double _x = mx - 0.5*mField->getNx();//(N_X/2, N_Y/2)を原点にシフト
	double _y = my - 0.5*mField->getNy();

	//中心との距離が半径+√2/2セル以上なら, 完全に媒質の外(念のため, 半径+1 以上か調べている)
	double _ax = ax+1, _by = by+1;
	if ((_x*_x) / (_ax*_ax) + (_y*_y) / (_by*_by) >= 1)
		return ep2;

	_ax = ax - 1;
	_by = by - 1;
	//中心との距離が, 半径-√2/2セル以内なら, 完全に媒質の外
	if ((_x*_x) / (_ax*_ax) + (_y*_y) / (_by*_by) <= 1)
		return ep1;

	double s = 0;

	double a = 1.0, b = 1.0;
	if (f == D_X) b = 0;
	if (f == D_Y) a = 0;
	for (double i = -16 + 0.5; i<16; i += 1)
		for (double j = -16 + 0.5; j<16; j += 1)
			if (pow(_x + a*i / 32.0, 2.0) / (ax*ax) + pow(_y + b*j / 32.0, 2.0) / (by*by) <= 1)
				s += 1;
	s /= 32.0*32.0;
	return ep1*s + ep2*(1 - s);
}

string FazzyHair_normalModel::mkdir(string root) {
	_mkdir((root + "HairModel").c_str());
	_mkdir((root + "HairModel/normalplane").c_str());
	
	string name = "HairModel/normalplane/e=" + to_s((double)e);
	_mkdir((root + name).c_str());	//ディレクトリの作成
	
	name = "HairModel/normalplane/e=" + to_s((double)e) + "/" + mField->getStringCellInfo();
	_mkdir((root + name).c_str());	//ディレクトリの作成
	return name + "/";
}