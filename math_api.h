#pragma once
#ifndef MATH_API_H
#define MATH_API_H

#endif // MATH_API_H
#include <math.h>
#include "QTextStream"
#include "QDebug"
inline int ind(int i,int j){
    return (i+1)*3+(j+1);
}

class Vector2d{
public:
    double x;
    double y;
    Vector2d();
    Vector2d(double a);
    Vector2d(double a,double b);

    double mag();
    Vector2d  operator % (Vector2d a);

    friend Vector2d operator+(Vector2d a,Vector2d b);

    friend Vector2d operator*(double a,Vector2d b);

    friend Vector2d operator*(Vector2d a,double b);

    friend double operator*(Vector2d a,Vector2d b);

    friend Vector2d operator>=(Vector2d a,Vector2d b);

    friend Vector2d operator<=(Vector2d a,Vector2d b);

};
QDebug operator<< (QDebug out,const Vector2d &r);
QTextStream &operator <<(QTextStream &s, const Vector2d &r);
std::ostream& operator <<(std::ostream& cout_,const Vector2d r);
