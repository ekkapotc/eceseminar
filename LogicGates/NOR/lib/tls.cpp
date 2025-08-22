#include "tls.hpp"
#include <cmath>
#include <iostream>

using namespace std;

Tangent::Tangent(const double &x) : v(x), t(0){};
Tangent::Tangent() : v(0), t(0){};

Tangent &Tangent::operator=(const Tangent &x) {
  if (this == &x)
    return *this;
  v = x.v;
  t = x.t;
  return *this;
}

Tangent operator*(const Tangent &x1, const Tangent &x2) {
  Tangent tmp;
  tmp.v = x1.v * x2.v;
  tmp.t = x1.t * x2.v + x1.v * x2.t;
  return tmp;
}

Tangent operator*(const Tangent &x1, double x2) {
  Tangent tmp;
  tmp.v = x1.v * x2;
  tmp.t = x1.t * x2;
  return tmp;
}

Tangent operator*(double x1, const Tangent &x2) {
  Tangent tmp;
  tmp.v = x1 * x2.v;
  tmp.t = x2.t * x1;
  return tmp;
}

Tangent operator/(const Tangent &x1, const Tangent &x2) {
  Tangent tmp;
  tmp.v = x1.v / x2.v;
  tmp.t = (1.0 / x2.v) * x1.t - (x1.v / (x2.v * x2.v)) * x2.t;
  return tmp;
}

Tangent operator/(const Tangent &x1, double x2) {
  Tangent tmp;
  tmp.v = x1.v / x2;
  tmp.t = (1.0 / x2) * x1.t;
  return tmp;
}

Tangent operator/(double x1, const Tangent &x2) {
  Tangent tmp;
  tmp.v = x1 / x2.v;
  tmp.t = -(x1 / (x2.v * x2.v)) * x2.t;
  return tmp;
}

Tangent operator+(const Tangent &x1, const Tangent &x2) {
  Tangent tmp;
  tmp.v = x1.v + x2.v;
  tmp.t = x1.t + x2.t;
  return tmp;
}

Tangent operator+(const Tangent &x1, double x2) {
  Tangent tmp;
  tmp.v = x1.v + x2;
  tmp.t = x1.t;
  return tmp;
}

Tangent operator+(double x1, const Tangent &x2) {
  Tangent tmp;
  tmp.v = x1 + x2.v;
  tmp.t = x2.t;
  return tmp;
}

Tangent operator-(const Tangent &x1, const Tangent &x2) {
  Tangent tmp;
  tmp.v = x1.v - x2.v;
  tmp.t = x1.t - x2.t;
  return tmp;
}

Tangent operator-(const Tangent &x1, double x2) {
  Tangent tmp;
  tmp.v = x1.v - x2;
  tmp.t = x1.t;
  return tmp;
}

Tangent operator-(double x1, const Tangent &x2) {
  Tangent tmp;
  tmp.v = x1 - x2.v;
  tmp.t = -x2.t;
  return tmp;
}

Tangent operator-(const Tangent &x) {
  Tangent tmp;
  tmp.v = -x.v;
  tmp.t = -x.t;
  return tmp;
}

Tangent sin(const Tangent &x) {
  Tangent tmp;
  tmp.v = sin(x.v);
  tmp.t = cos(x.v) * x.t;
  return tmp;
}

Tangent cos(const Tangent &x) {
  Tangent tmp;
  tmp.v = cos(x.v);
  tmp.t = -sin(x.v) * x.t;
  return tmp;
}

Tangent exp(const Tangent &x) {
  Tangent tmp;
  tmp.v = exp(x.v);
  tmp.t = tmp.v * x.t;
  return tmp;
}

Tangent sqrt(const Tangent &x) {
  Tangent tmp;
  tmp.v = sqrt(x.v);
  tmp.t = (0.5 / sqrt(x.v)) * x.t;
  return tmp;
}

Tangent log(const Tangent &x) {
  Tangent tmp;
  tmp.v = log(x.v);
  tmp.t = (1.0 / x.v) * x.t;
  return tmp;
}

istream &operator>>(istream &in, const Tangent &x) {
  in >> x.v;
  return in;
}

double abs(const Tangent &x) { return std::abs(x.v); }

bool operator==(const Tangent &x1, const Tangent &x2) { return x1.v == x2.v; }

bool operator==(const Tangent &x1, double x2) { return x1.v == x2; }

bool operator==(double x1, const Tangent &x2) { return x1 == x2.v; }

bool operator!=(const Tangent &x1, const Tangent &x2) { return x1.v != x2.v; }

bool operator!=(const Tangent &x1, double x2) { return x1.v != x2; }

bool operator!=(double x1, const Tangent &x2) { return x1 != x2.v; }

bool operator>(const Tangent &x1, const Tangent &x2) { return x1.v > x2.v; }

bool operator>(const Tangent &x1, double x2) { return x1.v > x2; }

bool operator>(double x1, const Tangent &x2) { return x1 > x2.v; }

bool operator>=(const Tangent &x1, const Tangent &x2) { return x1.v >= x2.v; }

bool operator>=(const Tangent &x1, double x2) { return x1.v >= x2; }

bool operator>=(double x1, const Tangent &x2) { return x1 >= x2.v; }

bool operator<(const Tangent &x1, const Tangent &x2) { return x1.v < x2.v; }

bool operator<(const Tangent &x1, double x2) { return x1.v < x2; }

bool operator<(double x1, const Tangent &x2) { return x1 < x2.v; }

bool operator<=(const Tangent &x1, const Tangent &x2) { return x1.v <= x2.v; }

bool operator<=(const Tangent &x1, double x2) { return x1.v <= x2; }

bool operator<=(double x1, const Tangent &x2) { return x1 <= x2.v; }
