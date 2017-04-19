#include "Object.h"

double s_to_d(const std::string& str){
	double rt;
	stringstream ss;
	ss << str;
	ss >> rt;
	return rt;
}
/*
var rgbTable:Array = [
    { L:380 ,R:0.06076 ,G:0.00000 ,B:0.11058 }
   ,{ L:390 ,R:0.08700 ,G:0.00000 ,B:0.16790 }
   ,{ L:400 ,R:0.13772 ,G:0.00000 ,B:0.26354 }
   ,{ L:410 ,R:0.20707 ,G:0.00000 ,B:0.39852 }
   ,{ L:420 ,R:0.31129 ,G:0.00000 ,B:0.60684 }
   ,{ L:430 ,R:0.39930 ,G:0.00000 ,B:0.80505 }
   ,{ L:440 ,R:0.40542 ,G:0.00000 ,B:0.87684 }
   ,{ L:450 ,R:0.34444 ,G:0.00000 ,B:0.88080 }
   ,{ L:460 ,R:0.11139 ,G:0.00000 ,B:0.86037 }
   ,{ L:470 ,R:0.00000 ,G:0.15233 ,B:0.77928 }
   ,{ L:480 ,R:0.00000 ,G:0.38550 ,B:0.65217 }
   ,{ L:490 ,R:0.00000 ,G:0.49412 ,B:0.51919 }
   ,{ L:500 ,R:0.00000 ,G:0.59271 ,B:0.40008 }
   ,{ L:510 ,R:0.00000 ,G:0.69549 ,B:0.25749 }
   ,{ L:520 ,R:0.00000 ,G:0.77773 ,B:0.00000 }
   ,{ L:530 ,R:0.00000 ,G:0.81692 ,B:0.00000 }
   ,{ L:540 ,R:0.00000 ,G:0.82625 ,B:0.00000 }
   ,{ L:550 ,R:0.00000 ,G:0.81204 ,B:0.00000 }
   ,{ L:560 ,R:0.47369 ,G:0.77626 ,B:0.00000 }
   ,{ L:570 ,R:0.70174 ,G:0.71523 ,B:0.00000 }
   ,{ L:580 ,R:0.84922 ,G:0.62468 ,B:0.00000 }
   ,{ L:590 ,R:0.94726 ,G:0.49713 ,B:0.00000 }
   ,{ L:600 ,R:0.99803 ,G:0.31072 ,B:0.00000 }
   ,{ L:610 ,R:1.00000 ,G:0.00000 ,B:0.00000 }
   ,{ L:620 ,R:0.95520 ,G:0.00000 ,B:0.00000 }
   ,{ L:630 ,R:0.86620 ,G:0.00000 ,B:0.00000 }
   ,{ L:640 ,R:0.76170 ,G:0.00000 ,B:0.00000 }
   ,{ L:650 ,R:0.64495 ,G:0.00000 ,B:0.00000 }
   ,{ L:660 ,R:0.52857 ,G:0.00000 ,B:0.00000 }
   ,{ L:670 ,R:0.41817 ,G:0.00000 ,B:0.00000 }
   ,{ L:680 ,R:0.33202 ,G:0.00000 ,B:0.00000 }
   ,{ L:690 ,R:0.25409 ,G:0.00000 ,B:0.00000 }
   ,{ L:700 ,R:0.19695 ,G:0.00000 ,B:0.00000 }
   ,{ L:710 ,R:0.15326 ,G:0.00000 ,B:0.00000 }
   ,{ L:720 ,R:0.11902 ,G:0.00000 ,B:0.00000 }
   ,{ L:730 ,R:0.09063 ,G:0.00000 ,B:0.00000 }
   ,{ L:740 ,R:0.06898 ,G:0.00000 ,B:0.00000 }
   ,{ L:750 ,R:0.05150 ,G:0.00000 ,B:0.00000 }
   ,{ L:760 ,R:0.04264 ,G:0.00000 ,B:0.00000 }
   ,{ L:770 ,R:0.03666 ,G:0.00000 ,B:0.00794 }
   ,{ L:780 ,R:0.00000 ,G:0.00000 ,B:0.00000 }
   ];
var wlenStep:Number = 10;
var wlenMin:Number  = 380;
var wlenMax:Number  = 780;

static map<double, struct Color> RGBTable;

color lambda_to_color(double lamda){ 
   wlen               = Math.max(wlenMin,wlen);
   wlen               = Math.min(wlenMax,wlen);
   var widx:int       = (wlen-wlenMin)/wlenStep;
   var wlenSub:Number = (wlen-wlenMin)-(widx*wlenStep);
   var ret            = new Object();
   ret.R              = rgbTable[widx].R;
   ret.G              = rgbTable[widx].G;
   ret.B              = rgbTable[widx].B;
   if( wlenSub==0 ) return ret;
   var rate:Number    = wlenSub/wlenStep;
   ret.R  += (rgbTable[widx+1].R-rgbTable[widx].R)*rate;
   ret.G  += (rgbTable[widx+1].G-rgbTable[widx].G)*rate;
   ret.B  += (rgbTable[widx+1].B-rgbTable[widx].B)*rate;
   return ret;
   }
}
*/