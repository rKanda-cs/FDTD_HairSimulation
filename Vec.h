#ifndef _VEC_H
#define _VEC_H
#include"Object.h"

template<class T> 
class Vec2{
public:
	T x,y;
	Vec2(T _x, T _y):x(_x),y(_y){};
	Vec2(){};
};

template<class T> double In_prod(const Vec2<T> &a, const Vec2<T> &b){
	return a.x*b.x + a.y*b.y;
}

template<class T> double Ou_prod(const Vec2<T> &a, const Vec2<T> &b){
	return a.x*b.y - a.y*b.x;
}

#endif