#ifndef _SIMULATOR_H
#define _SIMULATOR_H

#include"MenuWindow.h"
#include "NsFDTD_TM.h"
#include "NsFDTD_TE.h"
#include "StFDTD_TM.h"
#include "StFDTD_TE.h"
#include "Object.h"

using namespace std;

class Simulator:public Object{
private:
	Solver *solv;	//Solver

public:
	Simulator();
	~Simulator();

	virtual int calc();
	virtual void draw();
};


#endif //_SIMULATOR_H