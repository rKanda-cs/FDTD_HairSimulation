#ifndef _EVENT_STATE_H_
#define _EVENT_STATE_H_
const int _KEY_DEL = 127;
const int _KEY_0   = '0';
const int _KEY_9   = '9';
const int _KEY_ESC = '\033';
int getMouseButton();
int getMouseState();
int getMouseX();
int getMouseY();
int getMouseStX();
int getMouseStY();
unsigned char getKeyState();
void setSpecialKeyState(int);
void setMouseState(int, int, int, int);
void setKeyState(unsigned char);
int getSpecialKeyState();
int getMotion();
void setPos(int,int);

#endif