#include"Model.h"
#include"Field.h"
#include<time.h>
#include<fstream>
#include<stdlib.h>

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

	//中心との距離が, 半径-√2/2セル以内なら, 完全に媒質の中
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
FazzyModel(f), ep1(1.55*1.55*EPSILON_0_S), ep2(EPSILON_0_S), width1(530), width2(500), layer(7)
//[width1]nm間隔で[width2]nmのスラブを[layer]枚入れる  //隙間=width1-width2
{
//	pain = true;
	pain = false;

/*	cout << "random number" << endl;		// 層の乱数調整
	for (int i = 0; i < layer; i++) {
		random[i] = ((rand() % 5) - 4) * 20;
		cout << random[i] << endl;
	}
*/
	if (pain == true) {
		cout << "cuticle scratch" << endl;

		//--- 数値指定の場合 ---//
		lay[6] = 12; lay[5] = -1; lay[4] = -1; lay[3] = 6; lay[2] = -1; lay[1] = -1; lay[0] = 0;
		for (int i = layer; i > 0; i--) {
			cout << "slab" << i << " " << lay[i-1] << "deg" << endl;
		}
		//----------------------//
/*
		//--- ランダム指定の場合 ---//
		srand((unsigned int)time(NULL));

		lay[6] = 25 + (int)(rand()*(30 - 25 + 1.0) / (1.0 + RAND_MAX));
		cout << "slab7 " << lay[6] << "deg" << endl;
		for (int i = 5; i > 0; i--) {
			int t = 1 + (int)(rand()*(10 - 1 + 1.0) / (1.0 + RAND_MAX));
			lay[i] = lay[i + 1] - t;
			if (lay[i] < 0)	lay[i] = 0;
			cout << "slab" << i+1 << " " << lay[i] << "deg" << endl;
		}
		lay[0] = 0;
		cout << "slab1 " << lay[0] << "deg" << endl;
		//--------------------------//
*/

		for (int i = 0; i < 7; i++) {
			lay[i] = lay[i] * PI / 180;
			cout << lay[i] << endl;
		}
	}

	else if (pain == false) {
		lay[0] = 1; lay[1] = 1; lay[2] = 1; lay[3] = 1; lay[4] = 1; lay[5] = 1; lay[6] = 1;		// 層あり = 1  層なし = 0
	}
}

double FazzySlabModel::calcEPS(const double& x, const double& y, enum INTEG f){
//多層膜
	double h = 0;  //余白[nm]

	double mx = x - mField->getNpml(); //計算領域内へ写像
	double my = y - mField->getNpml();

	if(mx < 0 || my < 0 || mx >= mField->getNx() || my >= mField->getNy() ) return EPSILON_0_S;		// 計算範囲外
//	if (mx < mField->nanoToCell(h) || mx >= mField->getNx() - mField->nanoToCell(h))	return ep2;		// 左右の余白

//-----------------------------------------
/*
	// 隙間なし(対照実験1)
	if (my < mField->nanoToCell(h) || my >= mField->getNx() - mField->nanoToCell(h))	return ep2;
	return ep1;
*/

/*
	// CMC空間を同じ誘電率で埋める場合(対照実験2)
	if (my < mField->nanoToCell((width2 * layer) + (width1 - width2)*(layer - 1)))	return ep1;
	return ep2;
*/
//-----------------------------------------

	if (pain == true) {		// 剥がれによる傾きあり
		for (int i = 0; i < layer; i++) {		// iの最大値 = スラブの枚数
/*			// Fuzzyなし(Staircaseモデル)
			if (mField->cellToNano(my) > (mField->cellToNano(mx) * tan(lay[i]) + (h + i*width1)) && mField->cellToNano(my) < (mField->cellToNano(mx) * tan(lay[i]) + (h + i*width1 + width2)))
				return ep1;
*/
			// ver2 剥がれとめくれの両立
			if (mField->cellToNano(my) > (mField->cellToNano(mx) * tan(lay[i]) + (h + i*width1)) && mField->cellToNano(my) < (mField->cellToNano(mx) * tan(lay[i]) + (h + i*width1 + width2))) {
				if (lay[i] < 0)	// 例外処理で抜く層を指定
					return ep2;
				else
					return ep1;
			}
			

/*			// Fuzzyモデル
			double dy11 = my - (tan(lay[i]) * mx +			(mField->nanoToCell(h) + i*mField->nanoToCell(width1)));
			double dy12 = my - (tan(lay[i]) * (mx + 1) +	(mField->nanoToCell(h) + i*mField->nanoToCell(width1)));
			double dy21 = my - (tan(lay[i]) * mx +			(mField->nanoToCell(h) + i*mField->nanoToCell(width1) + mField->nanoToCell(width2)));
			double dy22 = my - (tan(lay[i]) * (mx + 1) +	(mField->nanoToCell(h) + i*mField->nanoToCell(width1) + mField->nanoToCell(width2)));
			double s;

			if ((dy21 < -1 && dy22 < -1) && (dy11 > 0 && dy12 > 0)) return ep1;		//キューティクル直線の内側 (2)

			if (dy11 <= 0 && dy12 <= 0) {
				if (fabs(dy11) <= 1 && fabs(dy12) <= 1) {
					s = (fabs(dy11) + fabs(dy12)) * 1.0 / 2.0;
					return ep2 * s + ep1 * (1 - s);		// (3)
				}
				if (fabs(dy11) < 1 && fabs(dy12) > 1) {
					s = (1 - fabs(dy11)) * (my / tan(lay[i]) - mx) / 2.0;
					return ep1 * s + ep2 * (1 - s);		// (4)
				}
			}
			if (dy11 > 0 && dy12 < 0) {
				s = fabs(dy12) * ((mx + 1) - my / tan(lay[i])) / 2.0;
				return ep2 * s + ep1 * (1 - s);		// (5)
			}
/* なぜかエラーになる　以下要改善
			if (dy21 <= 0 && dy22 <= 0) {
				if (fabs(dy21) <= 1 && fabs(dy22) <= 1) {
					s = (fabs(dy21) + fabs(dy22)) * 1.0 / 2.0;
					return ep1 * s + ep2 * (1 - s);		// (3)
				}
				if (fabs(dy21) < 1 && fabs(dy22) > 1) {
					s = (1 - fabs(dy21)) * (my / tan(lay[i]) - mx) / 2.0;
					return ep2 * s + ep1 * (1 - s);		// (4)
				}
			}
			if (dy21 > 0 && dy22 < 0) {
				s = fabs(dy22) * ((mx + 1) - my / tan(lay[i])) / 2.0;
				return ep1 * s + ep2 * (1 - s);		// (5)
			}
*/		}
		return ep2;
	}

	else if(pain == false) {	// 正常
		int k = (int)(mField->cellToNano(my) - h) % width1;
		double l = (mField->cellToNano(my) - h) / width1;

/*		// ver1 -----------------------//
		if (k > 0 && k <= width2 && l < layer)		// l = スラブの枚数
			return ep1;
		else
			return ep2;
		// ----------------------------//
*/
		// ver2 抜け層あり ------------//
		if (k > 0 && k <= width2 && l < layer) {		// l = スラブの枚数
			for (int i = 0; i < layer; i++) {
				if (l > i && l <= i + 1) {
					if (lay[i] == 1)	return ep1;
					else if (lay[i] == 0)	return ep2;
				}
			}
		}
		else
			return ep2;
		// ----------------------------//

/*		// ver3 乱数処理 --------------//
		for (int i = 0; i < layer; i++) {
			if (l > i && l <= i + 1) {
				if (k > 0 && k <= width2 + random[i] && l < layer) {		// l = スラブの枚数
					if (lay[i] == 1)	return ep1;
					else if (lay[i] == 0)	return ep2;
				}
				else
					return ep2;
			}
		}
		return ep2;
*/		// ----------------------------//
	}

}

string FazzySlabModel::mkdir(string root){
	_mkdir((root + "SlabModel").c_str());
	string name;

	if (pain == false) {
		name = "SlabModel/" + mField->getStringCellInfo();		//傷みなしの場合
	}
	else if (pain == true) {
		_mkdir((root + "SlabModel/scratch").c_str());			//傷みありの場合
		name = "SlabModel/scratch/" + mField->getStringCellInfo();
	}

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
	FazzyModel(f), ep1(1.55*1.55*EPSILON_0_S), ep2(EPSILON_0_S), alpha(5), cwidth(0.5), r(32+(4.36-cwidth))
	//alpha:キューティクルの角度(deg)  cwidth:キューティクルの幅(μm)  r:毛の半径(μm)(半径+キューティクルが重なる領域)
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

/*--------------------------------------------------*/
/*-----------縦断面(多層膜キューティクル)-----------*/
/*--------------------------------------------------*/
FazzyHair_incidenceLayerModel::FazzyHair_incidenceLayerModel(Field* f) :
	FazzyModel(f), ep1(1.55*1.55*EPSILON_0_S), ep2(EPSILON_0_S), ep3(mField->eps3),alpha(5), cwidth(0.5), length(50), r(32), cmc(0.06)
	//alpha:キューティクルの角度(deg)  cwidth:キューティクルの厚さ(μm)  length:キューティクルの長さ(μm)	r:毛皮質範囲の半径(μm)  cmc:CMC範囲(μm)
	//ep3:付着物質		屈折率: 水=1.33
{
	alphaR = alpha * PI / 180;
	int n = length * sin(alphaR) / (cmc + cwidth);

	cout << "キューティクルの角度 : " + to_s(alpha) + "deg" << endl;
	cout << "キューティクル厚さ : " + to_s(cwidth) + "micro" << endl;
	cout << "キューティクル長さ : " + to_s(length) + "micro" << endl;
	cout << "キューティクル1枚の露出幅 : " + to_s((cmc+cwidth)/sin(alphaR)) + "micro" << endl;
	cout << "キューティクル範囲幅 : " + to_s(length*sin(alphaR)) + "micro" << endl;
	cout << "キューティクル重なり枚数 : " + to_s(n) + "枚" << endl;
}

double FazzyHair_incidenceLayerModel::calcEPS(const double& x, const double& y, enum INTEG f) {
	alphaR = alpha * PI / 180;
	ln = mField->nanoToCell(length * 1000);
	lx = ln * cos(alphaR);
	ly = ln * sin(alphaR);
	rn = mField->nanoToCell(r * 1000);
	cn = mField->nanoToCell(cwidth * 1000);
	mn = mField->nanoToCell(cmc * 1000);

	double mx = x - mField->getNpml(); //計算領域内へ写像
	double my = y - mField->getNpml();
	double cx = mField->getNx() / 2;
	double cy = mField->getNy() / 2;

	if (mx < 0 || my < 0 || mx >= mField->getNx() || my >= mField->getNy()) return ep2;	//PML層

	/**** 毛髪全体 ****/
	if (my < cy)	my = 2 * cy - my;		//x軸に対して線対称
	if (my <= rn + cy)		return ep1;		//毛皮質部分
	my = my - rn - cy;		//軸移動

	/****************** キューティクル部分のみ ******************/
	/* Field推奨サイズ                                          */
	/* Field(8000, 8000, 5, 10) Field(16000, 8000, 10, 10) など */
	/************************************************************/
//	if (my <= 100)		return ep1;		//毛皮質部分
//	my = my - 100;		//軸移動


	if (my > ly)	return ep2;
	for (int i = -(ly / (cn + mn)); (i*mn + i*cn) / tan(alphaR) < mField->getNx(); i++) {
		if (my <= ly - cn) {
			if (my >= tan(alphaR) * (mx - (i*mn + (i + 1)*cn) / tan(alphaR)) + 1 && my <= tan(alphaR) * (mx - (i*mn + i*cn) / tan(alphaR)) - 1)
				return ep1;

			else if (my >= tan(alphaR) * (mx - ((i + 1)*mn + (i + 1)*cn) / tan(alphaR)) + 1 && my <= tan(alphaR) * (mx - (i*mn + (i + 1)*cn) / tan(alphaR)) - 1)
				return ep3;

			else if ((my > tan(alphaR) * (mx - (i*mn + i*cn) / tan(alphaR)) - 1 && my < tan(alphaR) * (mx - (i*mn + i*cn) / tan(alphaR)) + 1)
				|| (my > tan(alphaR) * (mx - (i*mn + (i + 1)*cn) / tan(alphaR)) - 1 && my < tan(alphaR) * (mx - (i*mn + (i + 1)*cn) / tan(alphaR)) + 1)) {
				
				double s = 0;
				for (double a = -16; a < 16; a += 1)
					for (double b = -16; b < 16; b += 1)
						if (my + b / 32.0 >= tan(alphaR) * (mx + a / 32.0 - (i*mn + (i + 1)*cn) / tan(alphaR)) && my + b / 32.0 <= tan(alphaR) * (mx + a / 32.0 - (i*mn + i*cn) / tan(alphaR)) && mx <= (ly + i*mn + i*cn) / tan(alphaR))
							s += 1;
				s /= 32.0*32.0;
				return ep1*s + ep3*(1 - s);
				
			}
			
		}
		else if (my > ly - cn) {
			if (mx >= (my + i*mn + i*cn + 1) / tan(alphaR) && mx <= (ly + i*mn + i*cn) / tan(alphaR))
				return ep1;

			else if (mx < (my + i*mn + i*cn + 1) / tan(alphaR) && mx >(my + i*mn + i*cn - 1) / tan(alphaR)) {

				double s = 0;
				for (double a = -16; a < 16; a += 1)
					for (double b = -16; b < 16; b += 1)
						if (mx + a / 32.0 >= (my + b / 32.0 + i*mn + i*cn) / tan(alphaR) && mx <= (ly + i*mn + i*cn) / tan(alphaR))
							s += 1;
				s /= 32.0*32.0;
				return ep1*s + ep3*(1 - s);
			}	
		}
	}

	return ep3;
}

double FazzyHair_incidenceLayerModel::calcSIG(const double& x, const double& y, const double lam, enum INTEG f) {
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

string FazzyHair_incidenceLayerModel::mkdir(string root) {
	_mkdir((root + "HairModel").c_str());
	string name;

	if (mField->sig == false) {
		if (mField->ep == false) {
			_mkdir((root + "HairModel/incidenceLayer").c_str());				//吸収係数なし&付着なしの場合
			name = "HairModel/incidenceLayer/" + mField->getStringCellInfo();
		}
		else if (mField->ep == true) {
			_mkdir((root + "HairModel/incidenceLayerWet").c_str());				//吸収係数なし&付着ありの場合
			name = "HairModel/incidenceLayerWet/" + mField->getStringCellInfo();
		}
	}
	else if (mField->sig == true) {
		if (mField->ep == false) {
			_mkdir((root + "HairModel/incidenceLayer_withSig").c_str());		//吸収係数あり&付着なしの場合
			name = "HairModel/incidenceLayer_withSig/" + mField->getStringCellInfo();
		}
		else if (mField->ep == true) {
			_mkdir((root + "HairModel/incidenceLayerWet_withSig").c_str());		//吸収係数あり&付着ありの場合
			name = "HairModel/incidenceLayerWet_withSig/" + mField->getStringCellInfo();
		}
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
	FazzyModel(f), ep1(1.55*1.55*EPSILON_0_S), ep2(EPSILON_0_S), e(0), r(40)
	//e:離心率  r:毛の半径(μm)
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


/*---------------------------------------------*/
/*--------------morphoModel-------------------------*/
/*---------------------------------------------*/
FazzyMorphoModel::FazzyMorphoModel(Field* f) :
	FazzyModel(f), ep1(1.55*1.55*EPSILON_0_S), ep2(EPSILON_0_S)
{
	ifstream ifp("morphoModel_cell.txt");
//	ifstream ifp("morphoModel_cell2.csv");

	if (!ifp) {
		cout << "not found input File" << endl;
		exit(1);
	}
	
	for (int y = 0; y < Y; y++) {
		for (int x = 0; x < X; x++) {
			ifp >> p[x][y];
		}
	}
}

double FazzyMorphoModel::calcEPS(const double& x, const double& y, enum INTEG f) {
	//多層膜
	double H = 1500;  //余白[nm]
	double h = mField->nanoToCell(H);

	double mx = x - mField->getNpml(); //計算領域内へ写像
	double my = y - mField->getNpml();

	if (mx < 0 || my < 0 || mx >= mField->getNx() || my >= mField->getNy()) return EPSILON_0_S;
	if (mx < h || mx >= mField->getNx() - h)	return ep2;

	mx -= h;
	if (p[(int)mx][(int)my] == 0)
		return ep2;
	else
		return ep1;

}

string FazzyMorphoModel::mkdir(string root) {
	_mkdir((root + "MorphoModel").c_str());

	string name = "MorphoModel/" + mField->getStringCellInfo();
	_mkdir((root + name).c_str());	//ディレクトリの作成
	return name + "/";
}