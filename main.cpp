#define _CRTDBG_MAP_ALLOC
#include <stdlib.h>
#include <crtdbg.h>
#include <math.h>
#include <stdio.h>
#include <string>
#include <glut.h>
#include "Simulator.h"
#include "EventState.h"
#include<Windows.h>
#define new ::new(_NORMAL_BLOCK, __FILE__, __LINE__)

//----------------------------------------------------
// 変数の宣言
//----------------------------------------------------
const int WindowPositionX = 200;  //生成するウィンドウ位置のX座標
const int WindowPositionY = 200;  //生成するウィンドウ位置のY座標
const char WindowTitle[] = "シミュレーション";  //ウィンドウのタイトル

Simulator *wave_simulator = new Simulator();

//----------------------------------------------------
// 関数プロトタイプ（後に呼び出す関数名と引数の宣言）
//----------------------------------------------------
void Initialize(void);   //初期設定時に呼び出す関数
void Idle(void);         //アイドル時に呼び出す関数
void Display(void);      //画面描画時に呼び出す関数
void Mouse(int _button , int _state , int _x , int _y);
void Motion(int _x,int _y);
void KeyBoard(unsigned char key, int x, int y);
void SpecialKeyBoard(int key, int x, int y);

//----------------------------------------------------
// メイン関数
//----------------------------------------------------

int main(int argc, char *argv[]){
//	cout.setf(ios::fixed, ios::floatfield);	//固定小数点, 18桁制度で指定
//	cout.precision(18);
	
	_CrtSetDbgFlag ( _CRTDBG_ALLOC_MEM_DF | _CRTDBG_LEAK_CHECK_DF );	//メモリリーク検出用
	glutInit(&argc, argv);                                     //環境の初期化
	glutInitWindowPosition(WindowPositionX, WindowPositionY);  //ウィンドウの位置の指定
	glutInitWindowSize(WINDOW_W, WINDOW_H);					//ウィンドウサイズの指定
	glutInitDisplayMode(GLUT_RGBA | GLUT_DEPTH | GLUT_DOUBLE); //ディスプレイモードの指定
	glutCreateWindow(WindowTitle);                             //ウィンドウの作成
	glutIdleFunc(Idle);                                        //プログラムアイドル状態時に呼び出される関数
	glutDisplayFunc(Display);                                  //描画時に呼び出される関数を指定する（関数名：Display）
	glutMouseFunc(Mouse);
	glutKeyboardFunc(KeyBoard);
	glutSpecialFunc(SpecialKeyBoard);
	glutMotionFunc(Motion);
	Initialize();                                              //初期設定の関数を呼び出す
	glutMainLoop();  
	return 0;

}
//----------------------------------------------------
// 初期設定の関数
//----------------------------------------------------
void Initialize(void){
	
}
//----------------------------------------------------
// アイドル時に呼び出される関数
//----------------------------------------------------
void Idle(){

	//if( getSpecialKeyState() == GLUT_KEY_UP){
		if(!wave_simulator->calc()){
			delete wave_simulator;
			exit(0);
		}
	//}

	if( getKeyState() == _KEY_ESC){
		delete wave_simulator;
		exit(0);
	}
	
	Sleep(1);
  setKeyState(-1);			//キーボード初期化
  setSpecialKeyState(-1);
  glutPostRedisplay(); //glutDisplayFunc()を１回実行する
}

void Mouse(int _button , int _state , int _x , int _y){
	setMouseState( _button,  _state, _x, _y);
}

void Motion(int _x,int _y){
	setPos(_x,_y);
}

void KeyBoard(unsigned char key, int x, int y){
	setKeyState(key);
}

void SpecialKeyBoard(int key, int x, int y){
	setSpecialKeyState(key);
}
//----------------------------------------------------
// 描画の関数
//----------------------------------------------------
void Display(void) {
	glClear(GL_COLOR_BUFFER_BIT);

	wave_simulator->draw();

	glutSwapBuffers(); //glutInitDisplayMode(GLUT_DOUBLE)でダブルバッファリングを利用可
	glFlush();
}
