#ifndef _Z_VECTOR2_H_
#define _Z_VECTOR2_H_

#include <math.h>

class vector2 {
public:
    double x, y;

    // constructor
    vector2(void);
    vector2(const double&, const double&);
    vector2(const int&, const int&);

    // operators
    vector2& operator  = (const vector2&);
    vector2& operator += (const vector2&);
    vector2& operator -= (const vector2&);
    vector2& operator *= (const double&);
    vector2& operator *= (const int&);
    vector2& operator /= (const double&);
    vector2& operator /= (const int&);

    friend vector2 operator + (const vector2&);
    friend vector2 operator - (const vector2&);
    friend vector2 operator + (const vector2&, const vector2&);
    friend vector2 operator - (const vector2&, const vector2&);

    // friend operators
    friend vector2 operator * (const double&, const vector2&);
    friend vector2 operator * (const int&, const vector2&);
    friend vector2 operator * (const vector2&, const double&);
    friend vector2 operator * (const vector2&, const int&);
    friend vector2 operator / (const vector2&, const double&);
    friend vector2 operator / (const vector2&, const int&);

    // member functions
    void set(const double&, const double&);
    void set(const int&, const int&);
    void set(void);
    void reset(void);

    double abs(void) const;
    double abs2(void) const;

    // friend functions
    friend double dot(const vector2&, const vector2&);
    friend double abs(const vector2&);

};

// ---- constructor ----
inline vector2::vector2(void)
    : x(0.0), y(0.0) {}

inline vector2::vector2(const double& _x, const double& _y)
    : x(_x), y(_y) {}

inline vector2::vector2(const int& _x, const int& _y)
    : x(static_cast<double>(_x)), y(static_cast<double>(_y)) {}

// ----  operators ------
inline vector2& vector2::operator = (const vector2& v) {
    x = v.x;
    y = v.y;
    return *this;
}

inline vector2& vector2::operator += (const vector2& v) {
    x += v.x;
    y += v.y;
    return *this;
}

inline vector2& vector2::operator -= (const vector2& v) {
    x -= v.x;
    y -= v.y;
    return *this;
}

inline vector2& vector2::operator *= (const double& d) {
    x *= d;
    y *= d;
    return *this;
}

inline vector2& vector2::operator *= (const int& i) {
    double d = static_cast<double>(i);
    x *= d;
    y *= d;
    return *this;
}

inline vector2& vector2::operator /= (const double& d) {
    x /= d;
    y /= d;
    return *this;
}

inline vector2& vector2::operator /= (const int& i) {
    double d = static_cast<double>(i);
    x /= d;
    y /= d;
    return *this;
}

inline vector2 operator + (const vector2& v) {
    return v;
}

inline vector2 operator - (const vector2& v) {
    return vector2(-v.x, -v.y);
}

inline vector2 operator + (const vector2& a, const vector2& b) {
    return vector2(a.x + b.x, a.y + b.y);
}

inline vector2 operator - (const vector2& a, const vector2& b) {
    return vector2(a.x - b.x, a.y - b.y);
}

inline vector2 operator * (const vector2& v, const double& d) {
    return vector2(d * v.x, d * v.y);
}

inline vector2 operator / (const vector2& v, const double& d) {
    return vector2(v.x / d, v.y / d);
}

// ---- member functions ----
inline void vector2::set(const double& _x, const double& _y) {
    x = _x;
    y = _y;
}

inline void vector2::set(const int& _x, const int& _y) {
    x = static_cast<double>(_x);
    y = static_cast<double>(_y);
}

inline void vector2::set(void) {
    x = 0.0;
    y = 0.0;
}

inline void vector2::reset(void) {
    x = 0.0;
    y = 0.0;
}

inline double vector2::abs2(void) const {
    return (x * x + y * y);
}

inline double vector2::abs(void) const {
    return sqrt(x * x + y * y);
}

// ---- useful functions ----
inline double dot(const vector2& a, const vector2& b) {
    return (a.x * b.x + a.y * b.y);
}

inline double abs(const vector2& v) {
    return v.abs();
}

#endif  // _Z_VECTOR2_H_
