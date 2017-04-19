#define _USE_MATH_DEFINES
#define _CRTDBG_MAP_ALLOC
#include <stdlib.h>
#include <crtdbg.h>
#include "Simulator.h"
#include "EventState.h"
#include "Solver.h"

using namespace std;
Simulator::Simulator() {
	//Field *mField = new Field(2000, 2000, 10, 20); //width, height, Δh, Npml
	//Model *mModel	 = new FazzyMieModel(mField, lambda_s);
	int mode;

	cout << "1: StFDTD_TM" << endl;
	cout << "2: StFDTD_TE" << endl;
	cout << "3: NsFDTD_TM" << endl;
	cout << "4: NsFDTD_TE" << endl;
	cin >> mode;
	if(mode == 1)	solv = new StFDTD_TM();
	else if(mode == 2)	solv = new StFDTD_TE();
	else if(mode == 3)	solv = new NsFDTD_TM();
	else if(mode == 4)	solv = new NsFDTD_TE();
	else exit(-1);

//	solv = new StFDTD_TM();
//	solv = new StFDTD_TE();
//	solv = new NsFDTD_TM();
//	solv = new NsFDTD_TE();

	solv->field();
}

Simulator::~Simulator() {
	cout << "Simulator Destructor" << endl;
	ButtonFactory::deleteAllButton();	//ボタンの削除, ウィンドウを閉じて終了した際に,すでに解放されたボタンを解放する恐れがある.(この順番だと大丈夫かも)
	delete solv;						//solverの削除
};

int Simulator::calc()
{
	try {
		solv->nextTime();
		return solv->calc();
	}
	catch (char *err) {
		cout << err << endl;
		return 0;
	}
	return 1;
}

void Simulator::draw()
{
	//	return;
	if (((int)solv->getTime()) % 20 != 0) return;
	solv->draw();					//シミュレーション状況描画
	ButtonFactory::draw();		//ボタンを描画
}
