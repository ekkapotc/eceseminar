#include <cmath>
#include <iostream>
#include <vector>

#include "adj.hpp"

int currPtr = 0;

ADJTapeEntry * adjTape{nullptr};

std::vector<int> indepPtr;
std::vector<int> depPtr;

Adjoint::Adjoint():va{0}, v{0.0}{}

Adjoint::Adjoint(const double & x) {
  v = x;
  adjTape[currPtr].oc = OpCode::CONST;
  adjTape[currPtr].v = x;
  va = currPtr++;

}

Adjoint::Adjoint(const Adjoint &x) {
  if (this != &x) {
    v = x.v;
    adjTape[currPtr].oc = OpCode::ASG;
    adjTape[currPtr].v = v;
    adjTape[currPtr].arg1 = x.va;
    va = currPtr++;
  }

#ifdef DCO_A1S_DEBUG
  std::cout << "Adjoint( const Adjoint & )" << std::endl;
#endif
}

Adjoint &Adjoint::operator=(const Adjoint &x) {
  if (this == &x)
    return *this;

  v = x.v;
  adjTape[currPtr].oc = OpCode::ASG;
  adjTape[currPtr].v = v;
  adjTape[currPtr].arg1 = x.va;
  va = currPtr++;

#ifdef DCO_A1S_DEBUG
  std::cout << "operator=( const Adjoint & )" << std::endl;
#endif

  return *this;
}

/*Adjoint & Adjoint::operator=( double x )
{
  v = x;
  adjTape[currPtr].oc = OpCode::CONST;
  adjTape[currPtr].v = x;
  va = currPtr++;
  std::cout << "operator=( double )" << std::endl;
  return *this;
}*/

Adjoint &Adjoint::operator+=(double x) {
  v += x;
  adjTape[currPtr].oc = OpCode::ASG;
  adjTape[currPtr].v = v;
  adjTape[currPtr].arg1 = va;
  va = currPtr++;

#ifdef DCO_A1S_DEBUG
  std::cout << "operator+=( double )" std::endl;
#endif

  return *this;
}

Adjoint &Adjoint::operator+=(const Adjoint &x) {
  v += x.v;
  adjTape[currPtr].oc = OpCode::AA_ADD;
  adjTape[currPtr].v = v;
  adjTape[currPtr].arg1 = va;
  adjTape[currPtr].arg2 = x.va;
  va = currPtr++;
  return *this;
}

Adjoint &Adjoint::operator-=(double x) {
  v -= x;
  adjTape[currPtr].oc = OpCode::ASG;
  adjTape[currPtr].v = v;
  adjTape[currPtr].arg1 = va;
  va = currPtr++;
  return *this;
}

Adjoint &Adjoint::operator-=(const Adjoint &x) {
  v -= x.v;
  adjTape[currPtr].oc = OpCode::AA_SUB;
  adjTape[currPtr].v = v;
  adjTape[currPtr].arg1 = va;
  adjTape[currPtr].arg2 = x.va;
  va = currPtr++;
  return *this;
}

Adjoint operator*(const Adjoint &x1, const Adjoint &x2) {
  Adjoint tmp;
  tmp.v = x1.v * x2.v;
  adjTape[currPtr].oc = OpCode::AA_MUL;
  adjTape[currPtr].arg1 = x1.va;
  adjTape[currPtr].arg2 = x2.va;
  adjTape[currPtr].v = tmp.v;
  tmp.va = currPtr++;

#ifdef DCO_A1S_DEBUG
  std::cout << "operator*( const Adjoint & , const Adjoint & )"
            << std::endl;
#endif

  return tmp;
}

Adjoint operator*(const Adjoint &x1, double x2) {
  Adjoint tmp;
  tmp.v = x1.v * x2;
  adjTape[currPtr].oc = OpCode::AP_MUL;
  adjTape[currPtr].arg1 = x1.va;
  adjTape[currPtr].v = tmp.v;
  adjTape[currPtr].val2 = x2;
  tmp.va = currPtr++;
  return tmp;
}

Adjoint operator*(double x1, const Adjoint &x2) {
  Adjoint tmp;
  tmp.v = x1 * x2.v;
  adjTape[currPtr].oc = OpCode::PA_MUL;
  adjTape[currPtr].arg2 = x2.va;
  adjTape[currPtr].v = tmp.v;
  adjTape[currPtr].val1 = x1;
  tmp.va = currPtr++;
  return tmp;
}

Adjoint operator/(const Adjoint &x1, const Adjoint &x2) {
  Adjoint tmp;
  tmp.v = x1.v / x2.v;
  adjTape[currPtr].oc = OpCode::AA_DIV;
  adjTape[currPtr].arg1 = x1.va;
  adjTape[currPtr].arg2 = x2.va;
  adjTape[currPtr].v = tmp.v;
  tmp.va = currPtr++;
  return tmp;
}

Adjoint operator/(const Adjoint &x1, double x2) {
  Adjoint tmp;
  tmp.v = x1.v / x2;
  adjTape[currPtr].oc = OpCode::AP_DIV;
  adjTape[currPtr].arg1 = x1.va;
  adjTape[currPtr].v = tmp.v;
  adjTape[currPtr].val2 =
      x2; // x2 is needed for the calculation of the partial derivative
  tmp.va = currPtr++;
  return tmp;
}

Adjoint operator/(double x1, const Adjoint &x2) {
  Adjoint tmp;
  tmp.v = x1 / x2.v;
  adjTape[currPtr].oc = OpCode::PA_DIV;
  adjTape[currPtr].arg2 = x2.va;
  adjTape[currPtr].v = tmp.v;
  adjTape[currPtr].val1 = x1;
  adjTape[currPtr].val2 = x2.v;
  tmp.va = currPtr++;
  return tmp;
}

Adjoint operator+(const Adjoint &x1, const Adjoint &x2) {
  Adjoint tmp;
  tmp.v = x1.v + x2.v;
  adjTape[currPtr].oc = OpCode::AA_ADD;
  adjTape[currPtr].arg1 = x1.va;
  adjTape[currPtr].arg2 = x2.va;
  adjTape[currPtr].v = tmp.v;
  tmp.va = currPtr++;

#ifdef DCO_A1S_DEBUG
  std::cout << "operator+( const Adjoint & , const Adjoint & )"
            << std::endl;
#endif

  return tmp;
}

Adjoint operator+(const Adjoint &x1, double x2) {
  Adjoint tmp;
  tmp.v = x1.v + x2;
  adjTape[currPtr].oc = OpCode::AP_ADD;
  adjTape[currPtr].arg1 = x1.va;
  adjTape[currPtr].v = tmp.v;
  tmp.va = currPtr++;

#ifdef DCO_A1S_DEBUG
  std::cout << "operator+( const Adjoint & , double )" << std::endl;
#endif

  return tmp;
}

Adjoint operator+(double x1, const Adjoint &x2) {
  Adjoint tmp;
  tmp.v = x1 + x2.v;
  adjTape[currPtr].oc = OpCode::PA_ADD;
  adjTape[currPtr].arg2 = x2.va;
  adjTape[currPtr].v = tmp.v;
  tmp.va = currPtr++;

#ifdef DCO_A1S_DEBUG
  std::cout << "operator+( double , const Adjoint & )" << std::endl;
#endif

  return tmp;
}

Adjoint operator-(const Adjoint &x1, const Adjoint &x2) {
  Adjoint tmp;
  tmp.v = x1.v - x2.v;
  adjTape[currPtr].oc = OpCode::AA_SUB;
  adjTape[currPtr].arg1 = x1.va;
  adjTape[currPtr].arg2 = x2.va;
  adjTape[currPtr].v = tmp.v;
  tmp.va = currPtr++;

#ifdef DCO_A1S_DEBUG
  std::cout << "operator-( const Adjoint & , const Adjoint & )"
            << std::endl;
#endif

  return tmp;
}

Adjoint operator-(const Adjoint &x1, double x2) {
  Adjoint tmp;
  tmp.v = x1.v - x2;
  adjTape[currPtr].oc = OpCode::AP_SUB;
  adjTape[currPtr].arg1 = x1.va;
  adjTape[currPtr].v = tmp.v;
  tmp.va = currPtr++;

#ifdef DCO_A1S_DEBUG
  std::cout << "operator-( const Adjoint & , double )" << std::endl;
#endif

  return tmp;
}

Adjoint operator-(double x1, const Adjoint &x2) {
  Adjoint tmp;
  tmp.v = x1 - x2.v;
  adjTape[currPtr].oc = OpCode::PA_SUB;
  adjTape[currPtr].arg2 = x2.va;
  adjTape[currPtr].v = tmp.v;
  tmp.va = currPtr++;

#ifdef DCO_A1S_DEBUG
  std::cout << "operator-( double , const Adjoint & )" << std::endl;
#endif

  return tmp;
}

Adjoint operator-(const Adjoint &x) {
  Adjoint tmp;
  tmp.v = -x.v;
  adjTape[currPtr].oc = OpCode::MINUS;
  adjTape[currPtr].arg1 = x.va;
  adjTape[currPtr].v = tmp.v;
  tmp.va = currPtr++;

#ifdef DCO_A1S_DEBUG
  std::cout << "operator-( const Adjoint & )" << std::endl;
#endif

  return tmp;
}

Adjoint sin(const Adjoint &x) {
  Adjoint tmp;
  tmp.v = std::sin(x.v);
  adjTape[currPtr].oc = OpCode::SIN;
  adjTape[currPtr].arg1 = x.va;
  adjTape[currPtr].v = tmp.v;
  tmp.va = currPtr++;

#ifdef DCO_A1S_DEBUG
  std::cout << "sin( const Adjoint & )" << std::endl;
#endif

  return tmp;
}

Adjoint cos(const Adjoint &x) {
  Adjoint tmp;
  tmp.v = std::cos(x.v);
  adjTape[currPtr].oc = OpCode::COS;
  adjTape[currPtr].arg1 = x.va;
  adjTape[currPtr].v = tmp.v;
  tmp.va = currPtr++;

#ifdef DCO_A1S_DEBUG
  std::cout << "cos( const Adjoint & )" << std::endl;
#endif

  return tmp;
}

Adjoint exp(const Adjoint &x) {
  Adjoint tmp;
  tmp.v = std::exp(x.v);
  adjTape[currPtr].oc = OpCode::EXP;
  adjTape[currPtr].arg1 = x.va;
  adjTape[currPtr].v = tmp.v;
  tmp.va = currPtr++;

#ifdef DCO_A1S_DEBUG
  std::cout << "exp( const Adjoint & )" << std::endl;
#endif

  return tmp;
}

Adjoint log(const Adjoint &x) {
  Adjoint tmp;
  tmp.v = std::log(x.v);
  adjTape[currPtr].oc = OpCode::LOG;
  adjTape[currPtr].arg1 = x.va;
  adjTape[currPtr].v = tmp.v;
  tmp.va = currPtr++;

#ifdef DCO_A1S_DEBUG
  std::cout << "log( const Adjoint & )" << std::endl;
#endif

  return tmp;
}

Adjoint log10(const Adjoint &x) {
  Adjoint tmp;
  tmp.v = std::log10(x.v);
  adjTape[currPtr].oc = OpCode::LOG10;
  adjTape[currPtr].arg1 = x.va;
  adjTape[currPtr].v = tmp.v;
  tmp.va = currPtr++;

#ifdef DCO_A1S_DEBUG
  std::cout << "log10( const Adjoint & )" << std::endl;
#endif

  return tmp;
}

Adjoint tan(const Adjoint &x) {
  Adjoint tmp;
  tmp.v = std::tan(x.v);
  adjTape[currPtr].oc = OpCode::TAN;
  adjTape[currPtr].arg1 = x.va;
  adjTape[currPtr].v = tmp.v;
  tmp.va = currPtr++;

#ifdef DCO_A1S_DEBUG
  std::cout << "tan( const Adjoint & )" << std::endl;
#endif

  return tmp;
}

Adjoint sqrt(const Adjoint &x) {
  Adjoint tmp;
  tmp.v = std::sqrt(x.v);
  adjTape[currPtr].oc = OpCode::SQRT;
  adjTape[currPtr].arg1 = x.va;
  adjTape[currPtr].v = tmp.v;
  tmp.va = currPtr++;

#ifdef DCO_A1S_DEBUG
  std::cout << "sqrt( const Adjoint & )" << std::endl;
#endif

  return tmp;
}

double abs(const Adjoint &x) { return std::abs(x.v); }

bool operator==(const Adjoint &x1, const Adjoint &x2) {
  return x1.v == x2.v;
}

bool operator==(const Adjoint &x1, double x2) { return x1.v < x2; }

bool operator==(double x1, const Adjoint &x2) { return x1 < x2.v; }

bool operator!=(const Adjoint &x1, const Adjoint &x2) {
  return x1.v != x2.v;
}

bool operator!=(const Adjoint &x1, double x2) { return x1.v != x2; }

bool operator!=(double x1, const Adjoint &x2) { return x1 != x2.v; }

bool operator<(const Adjoint &x1, const Adjoint &x2) {
  return x1.v < x2.v;
}

bool operator<(const Adjoint &x1, double x2) { return x1.v < x2; }

bool operator<(double x1, const Adjoint &x2) { return x1 < x2.v; }

bool operator<=(const Adjoint &x1, const Adjoint &x2) {
  return x1.v <= x2.v;
}

bool operator<=(const Adjoint &x1, double x2) { return x1.v <= x2; }

bool operator<=(double x1, const Adjoint &x2) { return x1 <= x2.v; }

bool operator>(const Adjoint &x1, const Adjoint &x2) {
  return x1.v > x2.v;
}

bool operator>(const Adjoint &x1, double x2) { return x1.v > x2; }

bool operator>(double x1, const Adjoint &x2) { return x1 > x2.v; }

bool operator>=(const Adjoint &x1, const Adjoint &x2) {
  return x1.v >= x2.v;
}

bool operator>=(const Adjoint &x1, double x2) { return x1.v >= x2; }

bool operator>=(double x1, const Adjoint &x2) { return x1 >= x2.v; }

std::istream &operator>>(std::istream &in, const Adjoint &x) {
  in >> x.v;
  return in;
}

void createTape(size_t size) {
  adjTape = new ADJTapeEntry[size];
}

void destroyTape() { delete[] adjTape; }

void registerIndependent(const Adjoint &x) {
  indepPtr.push_back(x.va);
}

void registerDependent(const Adjoint &x) {
  depPtr.push_back(x.va);
}

void printTape() {
  std::cout << "TAPE:" << std::endl;

  for (int i = 0; i < currPtr; i++) {
    std::cout << i << ": [ " << static_cast<int>(adjTape[i].oc) << " , "
         << adjTape[i].arg1 << " , " << adjTape[i].arg2 << " , "
         << adjTape[i].v << " , " << adjTape[i].a << " ] " << std::endl;
  }
}

void resetTape() {
  for (int i = 0; i < currPtr; i++) {
    adjTape[i].a = 0;
  }

  currPtr = 0;

  indepPtr.clear();
  depPtr.clear();
}

void interpretTape() {
  for (int i = currPtr - 1; i >= 0; i--) {
    switch (adjTape[i].oc) {
    case OpCode::CONST: {
      adjTape[i].a = 0.0;
      break;
    }

    case OpCode::ASG: {
      adjTape[adjTape[i].arg1].a += adjTape[i].a;
      break;
    }

    case OpCode::AA_ADD: {
      adjTape[adjTape[i].arg1].a += adjTape[i].a;
      adjTape[adjTape[i].arg2].a += adjTape[i].a;
      break;
    }

    case OpCode::AP_ADD: {
      adjTape[adjTape[i].arg1].a += adjTape[i].a;
      break;
    }

    case OpCode::PA_ADD: {
      adjTape[adjTape[i].arg2].a += adjTape[i].a;
      break;
    }

    case OpCode::AA_SUB: {
      adjTape[adjTape[i].arg1].a += adjTape[i].a;
      adjTape[adjTape[i].arg2].a -= adjTape[i].a;
      break;
    }

    case OpCode::AP_SUB: {
      adjTape[adjTape[i].arg1].a += adjTape[i].a;
      break;
    }

    case OpCode::PA_SUB: {
      adjTape[adjTape[i].arg2].a -= adjTape[i].a;
      break;
    }

    case OpCode::MINUS: {
      adjTape[adjTape[i].arg1].a -= adjTape[i].a;
    }

    case OpCode::AA_MUL: {
      adjTape[adjTape[i].arg1].a +=
          adjTape[adjTape[i].arg2].v * adjTape[i].a;
      adjTape[adjTape[i].arg2].a +=
          adjTape[adjTape[i].arg1].v * adjTape[i].a;
      break;
    }

    case OpCode::AP_MUL: {
      adjTape[adjTape[i].arg1].a +=
          (adjTape[i].val2) * adjTape[i].a;
      // adjTape[adjTape[i].arg1].a+=(adjTape[i].v/adjTape[adjTape[i].arg1].v)*adjTape[i].a;
      break;
    }

    case OpCode::PA_MUL: {
      adjTape[adjTape[i].arg2].a +=
          (adjTape[i].val1) * adjTape[i].a;
      // adjTape[adjTape[i].arg2].a+=(adjTape[i].v/adjTape[adjTape[i].arg2].v)*adjTape[i].a;
      break;
    }

    case OpCode::SIN: {
      adjTape[adjTape[i].arg1].a +=
          cos(adjTape[adjTape[i].arg1].v) * adjTape[i].a;
      break;
    }

    case OpCode::COS: {
      adjTape[adjTape[i].arg1].a -=
          sin(adjTape[adjTape[i].arg1].v) * adjTape[i].a;
      break;
    }

    case OpCode::EXP: {
      adjTape[adjTape[i].arg1].a +=
          exp(adjTape[adjTape[i].arg1].v) * adjTape[i].a;
      break;
    }

    case OpCode::LOG: {
      adjTape[adjTape[i].arg1].a +=
          (1.0 / adjTape[adjTape[i].arg1].v) * adjTape[i].a;
      break;
    }

    case OpCode::LOG10: {
      adjTape[adjTape[i].arg1].a +=
          (1.0 / (adjTape[adjTape[i].arg1].v * std::log(10))) *
          adjTape[i].a;
      break;
    }

    case OpCode::TAN: {
      adjTape[adjTape[i].arg1].a +=
          (1 + (adjTape[adjTape[i].arg1].v) *
                   (adjTape[adjTape[i].arg1].v)) *
          adjTape[i].a;
      break;
    }

    case OpCode::AA_DIV: {
      adjTape[adjTape[i].arg1].a +=
          (1.0 / adjTape[adjTape[i].arg2].v) * adjTape[i].a;
      adjTape[adjTape[i].arg2].a -=
          (adjTape[adjTape[i].arg1].v /
           std::pow(adjTape[adjTape[i].arg2].v, 2)) *
          adjTape[i].a;
      break;
    }

    case OpCode::AP_DIV: {
      adjTape[adjTape[i].arg1].a +=
          (1.0 / adjTape[i].val2) * adjTape[i].a;
      break;
    }

    case OpCode::PA_DIV: {
      adjTape[adjTape[i].arg2].a -=
          (adjTape[i].val1 / std::pow(adjTape[i].val2, 2.0)) *
          adjTape[i].a;
      break;
    }

    case OpCode::SQRT: {
      if (adjTape[adjTape[i].arg1].v > 0.0)
        adjTape[adjTape[i].arg1].a +=
            0.5 * (1.0 / std::sqrt(adjTape[adjTape[i].arg1].v)) *
            adjTape[i].a;
      else
        adjTape[adjTape[i].arg1].a +=
            0.5 *
            (1.0 / std::sqrt(adjTape[adjTape[i].arg1].v + 1.0e-12)) *
            adjTape[i].a;
      break;
    }
    }
  }
}

double value(const Adjoint &x) { return x.v; }

double derivative(const int i) {
  return adjTape[indepPtr[i]].a;
}

