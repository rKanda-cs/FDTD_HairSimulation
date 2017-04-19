#ifndef _PHYSICAL_H
#define _PHYSICAL_H

#include <math.h>
/*-----------------------------*/
/*物理定数のシミュレーション値 */
/*-----------------------------*/
const double LIGHT_SPEED_S  = 0.7;	//光速
const double EPSILON_0_S    = 1.0;	//真空の誘電率
const double MU_0_S	        = 1.0/0.7/0.7;	//真空の透磁率　//Debugモードだと,0.7をLight_SPEED_Sと書くと,0で初期化されてしまう.(Release　ならok)
const double Z_0_S	        = sqrt(1.0/0.7/0.7/1.0);	//波動インピーダンス z = √(μ/ε)
#endif