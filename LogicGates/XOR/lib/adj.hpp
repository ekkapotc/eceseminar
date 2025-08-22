#ifndef ADJ_HPP_INCLUDED_
#define ADJ_HPP_INCLUDED_

#include <cmath>
#include <vector>
#include <iostream>

constexpr int TAPE_SIZE{100'000'000};

extern int  currPtr;
extern std::vector<int> indepPtr;
extern std::vector<int> depPtr;

enum class OpCode {
  UNDEF = -1,
  CONST = 0,
  ASG,

  AA_ADD,
  AP_ADD,
  PA_ADD,

  AA_SUB,
  AP_SUB,
  PA_SUB,

  AA_MUL,
  AP_MUL,
  PA_MUL,

  SIN,
  COS,
  EXP,
  MINUS,

  P_INC,
  A_INC,

  P_DEC,
  A_DEC,

  LOG,
  LOG10,
  TAN,

  AA_DIV,
  AP_DIV,
  PA_DIV,

  SQRT,
  AP_POW
};

struct ADJTapeEntry{
  OpCode oc{OpCode::UNDEF};

  int arg1{0};
  int arg2{0};

  double v{0.0};
  double a{0.0};

  double val1{0.0};
  double val2{0.0};
};

struct Adjoint{

  mutable int va;
  mutable double v;

  Adjoint();

  Adjoint(const double &);

  Adjoint(const Adjoint &);

  Adjoint &operator=(const Adjoint &);

  // Adjoint & operator=( double );

  Adjoint &operator+=(double);

  Adjoint &operator+=(const Adjoint &);

  Adjoint &operator-=(double);

  Adjoint &operator-=(const Adjoint &);


  double getValue();

  double getAdjoint();

  void setValue(double );
 
  void setAdjoint(double ); 
};

Adjoint operator*(const Adjoint &, const Adjoint &);
Adjoint operator*(const Adjoint &, double);
Adjoint operator*(double, const Adjoint &);

Adjoint operator/(const Adjoint &, const Adjoint &);
Adjoint operator/(const Adjoint &, double);
Adjoint operator/(double, const Adjoint &);

Adjoint operator+(const Adjoint &, const Adjoint &);
Adjoint operator+(const Adjoint &, double);
Adjoint operator+(double, const Adjoint &);

Adjoint operator-(const Adjoint &, const Adjoint &);
Adjoint operator-(const Adjoint &, double);
Adjoint operator-(double, const Adjoint &);
Adjoint operator-(const Adjoint &);

Adjoint sin(const Adjoint &);
Adjoint cos(const Adjoint &);
Adjoint exp(const Adjoint &);
Adjoint log(const Adjoint &);
Adjoint log10(const Adjoint &);
Adjoint tan(const Adjoint &);
Adjoint sqrt(const Adjoint &);

double abs(const Adjoint &);

bool operator==(const Adjoint &, const Adjoint &);
bool operator==(const Adjoint &, double);
bool operator==(double, const Adjoint &);

bool operator!=(const Adjoint &, const Adjoint &);
bool operator!=(const Adjoint &, double);
bool operator!=(double, const Adjoint &);

bool operator<(const Adjoint &, const Adjoint &);
bool operator<(const Adjoint &, double);
bool operator<(double, const Adjoint &);

bool operator<=(const Adjoint &, const Adjoint &);
bool operator<=(const Adjoint &, double);
bool operator<=(double, const Adjoint &);

bool operator>(const Adjoint &, const Adjoint &);
bool operator>(const Adjoint &, double);
bool operator>(double, const Adjoint &);

bool operator>=(const Adjoint &, const Adjoint &);
bool operator>=(const Adjoint &, double);
bool operator>=(double, const Adjoint &);

std::istream &operator>>(std::istream &in, const Adjoint &x);

void createTape(size_t size);
void destroyTape();
void registerIndependent(const Adjoint &);
void registerDependent(const Adjoint &);
void printTape();
void interpretTape();
void resetTape();


#endif
