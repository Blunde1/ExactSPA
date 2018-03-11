// Copyright (C) 2016-2017 Berent Lunde

/**
* Makes the following available for the TMB user:
* The "cType<T>" complex AD data T.
* The arithmetic follows standard arithmetic for complex variables.
* It is defined to work together with the standard TMB data T "T".
* Constructor: cType<T> z; defaults to 0+0*i.
*              cType<T> z((T)Re,(T)Im) = Re + i*Im.
* Compound assignements (+=, -=, *=, /=).
* Unary operators (+, -).
* Arithmetic (+, -, *, /).
* Relational and comparison operators (==, !=).
* Standard functions for complex variables (abs, arg, conj, real, imag).
* Exponential functions (exp, log).
* Power functions (pow, sqrt).
* Trigonometric functions (sin, cos, tan, asin, acos, atan, sinh, cosh, tanh).
*
*/

#ifndef __COMPLEX_HPP_INCLUDED__
#define __COMPLEX_HPP_INCLUDED__

template<typename T>
struct cType{

    // MEMBER FUNCTIONS
    T r, i;

    // Constructor
    cType(void) {r=0;i=0;}
    cType(T r_, T i_=0) : r(r_), i(i_) {}
    cType(const cType& c) : r(c.r), i(c.i) {}

    // Compound assignements
    cType& operator =(const T& t){
        r = t;
        i = 0;
        return *this;
    }
    cType& operator =(const cType& c){
        r = c.r;
        i = c.i;
        return *this;
    }
    cType& operator +=(const T& t){
        r = r + t;
        return *this;
    }
    cType& operator +=(const cType& c){
        r = r + c.r;
        i = i + c.i;
        return *this;
    }
    cType& operator -=(const T& t){
        r = r - t;
        return *this;
    }
    cType& operator -=(const cType& c){
        r = r - c.r;
        i = i - c.i;
        return *this;
    }
    cType& operator *=(const T& t){
        r = r*t;
        i = i*t;
        return *this;
    }
    cType& operator *=(const cType& c){
        T tmp_r, tmp_i;
        tmp_r = r*c.r - i*c.i;
        tmp_i = r*c.i + i*c.r;
        r = tmp_r, i=tmp_i;
        return *this;
    }
    cType& operator /=(const T& t){
        r = r / t;
        i = i / t;
        return *this;
    }
    cType& operator /=(const cType& c){
        T div = c.r*c.r + c.i*c.i, tmp_r, tmp_i;
        tmp_r = (r*c.r + i*c.i)/div;
        tmp_i = (i*c.r - r*c.i)/div;
        r = tmp_r, i = tmp_i;
        return *this;
    }
};

// NON-MEMBER FUNCTIONS
// Unary operators
template<typename T>
cType<T> operator+ (const cType<T>& z){
    cType<T> c = z;
    return c;
}
template<typename T>
cType<T> operator- (const cType<T>& z){
    cType<T> c = z, zero;
    return zero -= c;
}

// Arithmetic
template<typename T>
cType<T> operator +(const cType<T>& c, const T& t){
    cType<T> res = c;
    return res+=t;
}
template<typename T>
cType<T> operator +(const T& t, const cType<T>& c){
    cType<T> res = c;
    return res+=t;
}
template<typename T>
cType<T> operator +(const cType<T>& c_1, const cType<T>& c_2){
    cType<T> res = c_1;
    return res += c_2;
}
template<typename T>
cType<T> operator -(const cType<T>& c, const T& t){
    cType<T> res = c;
    return res-=t;
}
template<typename T>
cType<T> operator -(const T& t, const cType<T>& c){
    cType<T> res(t,0);
    return res -= c;
}
template<typename T>
cType<T> operator -(const cType<T>& c_1, const cType<T>& c_2){
    cType<T> res = c_1;
    return res -= c_2;
}
template<typename T>
cType<T> operator *(const cType<T>& c, const T& t){
    cType<T> res = c;
    return res *= t;
}
template<typename T>
cType<T> operator *(const T& t, const cType<T>& c){
    cType<T> res(t,0);
    return res *= c;
}
template<typename T>
cType<T> operator *(const cType<T>& c_1, const cType<T>& c_2){
    cType<T> c1 = c_1, c2=c_2;
    return c1 *= c2;
}
template<typename T>
cType<T> operator /(const cType<T>& c, const T& t){
    cType<T> res = c;
    return res /= t;
}
template<typename T>
cType<T> operator /(const T& t, const cType<T>& c){
    cType<T> res(t,0);
    return res /= c;
}
template<typename T>
cType<T> operator /(const cType<T>& c_1, const cType<T>& c_2){
    cType<T> res = c_1;
    return res /= c_2;
}

// Relational and comparison operators
template<typename T>
bool operator ==(const cType<T>& lhs, const cType<T>& rhs){
    cType<T> c_1=lhs, c_2 = rhs;
    if(c_1.r == c_2.r && c_1.i == c_2.i) {return true;}
    else{return false;}
}
template<typename T>
bool operator ==(const cType<T>& lhs, const T& rhs){
    cType<T> c_1=lhs, c_2(rhs,0);
    if(c_1.r == c_2.r && c_1.i == c_2.i) {return true;}
    else{return false;}
}
template<typename T>
bool operator ==(const T& lhs, const cType<T>& rhs){
    cType<T> c_1(lhs,0), c_2=rhs;
    if(c_1.r == c_2.r && c_1.i == c_2.i) {return true;}
    else{return false;}
}
template<typename T>
bool operator !=(const cType<T>& lhs, const cType<T>& rhs){
    cType<T> c_1=lhs, c_2 = rhs;
    return !(c_1==c_2);
}
template<typename T>
bool operator !=(const cType<T>& lhs, const T& rhs){
    cType<T> c_1=lhs, c_2(rhs,0);
    return !(c_1==c_2);
}
template<typename T>
bool operator !=(const T& lhs, const cType<T>& rhs){
    cType<T> c_1(lhs,0), c_2=rhs;
    return !(c_1==c_2);
}

// Standard functions
template<typename T>
T abs(const cType<T>& z){
    cType<T> c = z;
    T abs = sqrt(c.r*c.r + c.i*c.i);
    return abs;
}
// atan2 function
template<typename T>
T atan2(T y, T x){
    T res;
    if(x>0){res = atan(y/x); }
    else if(y >= 0 && x<0){ res = atan(y/x) + M_PI;}
    else if(y<0 && x<0){res = atan(y/x) - M_PI;}
    else if(y>0 && x==0){res=M_PI/2.0;}
    else if(y<0 && x==0){res=-M_PI/2.0;}
    return res;
}

template<typename T>
T arg(const cType<T>& z){ // Returns Arg(z)
    cType<T> c = z;
    T arg;
    if(c.i==0 && c.r>0){arg = 0; }
    else if(c.i==0 && c.r<0){arg = (T)M_PI; }
    else{ arg = atan2(c.i,c.r);}
    return arg;
}
template<typename T>
cType<T> conj(const cType<T>& z) {
    cType<T> c = z;
    c.i = - c.i;
    return c;
}
template<typename T>
T real(const cType<T>& z){
    cType<T> c = z;
    return c.r;
}
template<typename T>
T imag(const cType<T>& z){
    cType<T> c = z;
    return c.i;
}

// Exponential functions
template<typename T>
cType<T> exp(const cType<T>& z) {
    cType<T> c = z;
    T temp = exp(c.r)*cos(c.i);
    c.i = exp(c.r)*sin(c.i);
    c.r = temp;
    return c;
}
template<typename T>
cType<T> log(const cType<T>& z){
    cType<T> c = z;
    T r = abs(z), theta = arg(z);
    c.r = log(r);
    c.i = theta;
    return c;
}

// Power functions
template<typename T>
cType<T> pow(const cType<T>& z1, const cType<T>& z2){
    cType<T> c1=z1, c2=z2, c3;
    T a=c1.r,b=c1.i, c=c2.r,d=c2.i;
    c3.r = pow((a*a+b*b),(c/2))*exp(-d*arg(c1))*
        (cos(c*arg(c1)+0.5*d*log(a*a+b*b)));
    c3.i = pow((a*a+b*b),(c/2))*exp(-d*arg(c1))*
        (sin(c*arg(c1)+0.5*d*log(a*a+b*b)));
    return c3;
}
template<typename T>
cType<T> pow(const cType<T>& z, const T& t){
    cType<T> c1=z, c2(t,0);
    c1 = pow(c1,c2);
    return c1;
}
template<typename T>
cType<T> pow(const T& t, const cType<T>& z){
    cType<T> c1=z, c2(t,0);
    c1 = pow(c2,c1);
    return c1;
}
template<typename T>
cType<T> sqrt(const cType<T>& z){
    cType<T> c1=z, c2(0.5,0);
    c1 = pow(c1,c2);
    return c1;
}

// Trigonometric functions
template<typename T>
cType<T> sin(const cType<T>& z){
    cType<T> c = z, i(0,1);
    c = (exp(i*c) - exp((-(T)1)*i*c)) / ((T)2 * i);
    return c;
}
template<typename T>
cType<T> cos(const cType<T>& z){
    cType<T> c = z, i(0,1);
    c = (exp(i*c) + exp(c/i)) / ((T)2);
    return c;
}
template<typename T>
cType<T> tan(const cType<T>& z){
    cType<T> c = z;
    c = sin(c) / cos(c);
    return c;
}
template<typename T>
cType<T> asin(const cType<T>& z){
    cType<T> c = z, i(0,1);
    c = (-(T)1)*i*log( i*c + sqrt((T)1 - (c*c)) );
    return c;
}
template<typename T>
cType<T> acos(const cType<T>& z){
    cType<T> c = z;
    c = asin((-(T)1 ) * c) + (T)(M_PI/2);
    return c;
}
template<typename T>
cType<T> atan(const cType<T>& z){
    cType<T> c = z, i(0,1);
    c = (-(T)(0.5))*i*(log((T)1 - (i*c)) - log((T)1 + (i*c) ) );
    return c;
}


template<typename T>
cType<T> sinh(const cType<T>& z){
    cType<T> c = z, i(0,1);
    c = sin((-(T)1)*i*c) * i;
    return c;
}
template<typename T>
cType<T> cosh(const cType<T>& z){
    cType<T> c = z, i(0,1);
    c = cos(c/i);
    return c;
}
template<typename T>
cType<T> tanh(const cType<T>& z){
    cType<T> c = z;
    c = sinh(c) / cosh(c);
    return c;
}

#endif // __COMPLEX_HPP_INCLUDED__
