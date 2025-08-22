#include <cmath>
#include <cstdlib>
#include <ctime>
#include <iostream>
#include <tuple>
#include <vector>

#include "adj.hpp"

/*
 * In this problem,
 * we will determine the two parameters m and c of a given linear model 'y = mx
 * + c'
 */

#define value(x) (x).v

extern ADJTapeEntry *adjTape;

template <typename T> using Tuple = std::tuple<T, T>;

template <typename T> using Vector = std::vector<T>;

template <typename T, typename S> T J(T m, T c, Vector<Tuple<S>> &data) {

  T acc{0.0};
  for (auto &e : data) {
    S x = std::get<0>(e);
    S y_exp = std::get<1>(e);
    T y_act = m * x + c;
    acc += (y_exp - y_act) * (y_exp - y_act);
  }

  return acc / data.size();
}

int main(int argc, char *argv[]) {

  constexpr size_t MAX_EPOCHS{50};

  Vector<Tuple<double>> training_data = {
      {0.0, 1.0}, {1.0, 3.0}, {2.0, 5.0}, {3.0, 7.0}};

  createTape(100000);

  auto no_samples = training_data.size();

  srand(time(nullptr));

  Adjoint m{static_cast<double>(rand()) / RAND_MAX};
  Adjoint c{static_cast<double>(rand()) / RAND_MAX};

  std::cout << "Epoch " << 0 << std::endl;
  std::cout << "\t m = " << value(m) << " , c = " << value(c)
            << " , cost = " << J(value(m), value(c), training_data)
            << std::endl;
  double step{1e-2};

  for (auto epoch{0}; epoch < MAX_EPOCHS; epoch++) {

    registerIndependent(m);
    registerIndependent(c);

    Adjoint cost = J(m, c, training_data);

    registerDependent(cost);

    adjTape[cost.va].a = 1.0;

    interpretTape();

    double dJdm = adjTape[indepPtr[0]].a;
    double dJdc = adjTape[indepPtr[1]].a;

    m -= step * dJdm;
    c -= step * dJdc;

    std::cout << "Epoch " << epoch + 1 << std::endl;
    std::cout << "\t m = " << value(m) << " , c = " << value(c)
              << " , cost = " << J(value(m), value(c), training_data)
              << std::endl;
    resetTape();
  }

  destroyTape();

  return 0;
}
