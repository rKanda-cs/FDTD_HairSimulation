#include "EventState.h"
#include "stdio.h"
#include <glut.h>

int MOUSE_BUTTON = 0;
int MOUSE_STATE = 0;
int MOUSE_X = 0;
int MOUSE_Y = 0;
int MOTION = 0;
int MOUSE_ST_X = 0;
int MOUSE_ST_Y = 0;
unsigned char KEY_STATE = -1;
int SPECIALKEY_STATE = -1;

void setMouseState(int _button, int _state, int _x, int _y ){
	MOUSE_BUTTON = _button;
	MOUSE_STATE  = _state;
	MOUSE_X      = _x;
	MOUSE_Y      = _y;
	if(MOUSE_STATE == GLUT_DOWN) {MOTION = 1; MOUSE_ST_X = _x; MOUSE_ST_Y = _y;}
	if(MOUSE_STATE == GLUT_UP  ) MOTION = 0; 
}

void setPos(int _x,int _y){
	MOUSE_X = _x;
	MOUSE_Y = _y;
}
int getMotion(){
	return MOTION;
}

int getMouseButton(){ 
	return MOUSE_BUTTON;
}

int getMouseState(){
	return MOUSE_STATE;
}

int getMouseX(){
	return MOUSE_X;
}

int getMouseY(){
	return MOUSE_Y;
}

int getMouseStX(){
	return MOUSE_ST_X;
}

int getMouseStY(){
	return MOUSE_ST_Y;
}


void setKeyState(unsigned char _key){
	KEY_STATE = _key;
}

unsigned char getKeyState(){
	return KEY_STATE;
}

void setSpecialKeyState(int _key){
	SPECIALKEY_STATE = _key;
}

int getSpecialKeyState(){
	return SPECIALKEY_STATE;
}