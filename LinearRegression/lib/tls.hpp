#ifndef TLS_HPP_INCLUDE
#define TLS_HPP_INCLUDE

#include <cmath>
#include <iostream>

class Tangent {
public:
  mutable double v;
  mutable double t;

  Tangent(const double &x);
  Tangent();
  Tangent &operator=(const Tangent &x);
  Tangent &operator+=(const Tangent &x);
  Tangent &operator-=(const Tangent &x);
  Tangent &operator*=(const Tangent &x);
  Tangent &operator/=(const Tangent &x);
};

Tangent operator*(const Tangent &x1, const Tangent &x2);
Tangent operator*(const Tangent &x1, double x2);
Tangent operator*(double x1, const Tangent &x2);

Tangent operator/(const Tangent &x1, const Tangent &x2);
Tangent operator/(const Tangent &x1, double x2);
Tangent operator/(double x1, const Tangent &x2);

Tangent operator+(const Tangent &x1, const Tangent &x2);
Tangent operator+(const Tangent &x1, double x2);
Tangent operator+(double x1, const Tangent &x2);

Tangent operator-(const Tangent &x1, const Tangent &x2);
Tangent operator-(const Tangent &x1, double x2);
Tangent operator-(double x1, const Tangent &x2);

Tangent operator-(const Tangent &x);

Tangent sin(const Tangent &x);
Tangent cos(const Tangent &x);
Tangent exp(const Tangent &x);
Tangent sqrt(const Tangent &x);
Tangent log(const Tangent &x);

std::istream &operator>>(std::istream &in, const Tangent &x);

double abs(const Tangent &x);

bool operator==(const Tangent &, const Tangent &);
bool operator==(const Tangent &, double);
bool operator==(double, const Tangent &);

bool operator!=(const Tangent &, const Tangent &);
bool operator!=(const Tangent &, double);
bool operator!=(double, const Tangent &);

bool operator<(const Tangent &, const Tangent &);
bool operator<(const Tangent &, double);
bool operator<(double, const Tangent &);

bool operator<=(const Tangent &, const Tangent &);
bool operator<=(const Tangent &, double);
bool operator<=(double, const Tangent &);

bool operator>(const Tangent &, const Tangent &);
bool operator>(const Tangent &, double);
bool operator>(double, const Tangent &);

bool operator>=(const Tangent &, const Tangent &);
bool operator>=(const Tangent &, double);
bool operator>=(double, const Tangent &);

#endif
