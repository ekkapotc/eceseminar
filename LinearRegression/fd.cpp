#include <cmath>
#include <cstdlib>
#include <ctime>
#include <iostream>
#include <tuple>
#include <vector>

/*
 * In this problem,
 * we will determine the two parameters m and c of a given linear model 'y = mx
 * + c'
 */

template <typename T> using Tuple = std::tuple<T, T>;

template <typename T> using Vector = std::vector<T>;

double J(double m, double c, Vector<Tuple<double>> &data) {

  double acc{0.0};

  for (auto &e : data) {
    auto x = std::get<0>(e);
    auto y_exp = std::get<1>(e);
    auto y_act = m * x + c;
    acc += (y_exp - y_act) * (y_exp - y_act);
  }

  return acc / data.size();
}

int main(int argc, char *argv[]) {

  Vector<Tuple<double>> training_data = {
      {0.0, 1.0}, {1.0, 3.0}, {2.0, 5.0}, {3.0, 7.0}};

  auto no_samples = training_data.size();

  srand(time(nullptr));

  double m{static_cast<double>(rand()) / RAND_MAX};
  double c{static_cast<double>(rand()) / RAND_MAX};

  double cost{J(m, c, training_data)};

  std::cout << "Iteration " << 0 << std::endl;
  std::cout << "\t m = " << m << " , c = " << c
            << " , cost = " << J(m, c, training_data) << std::endl;

  double epsilon{1e-6};
  double step{1e-2};

  for (auto i{0}; i < 5000; i++) {

    double dJdm = (J(m + epsilon, c, training_data) - cost) / epsilon;
    double dJdc = (J(m, c + epsilon, training_data) - cost) / epsilon;

    m -= step * dJdm;
    c -= step * dJdc;

    std::cout << "Iteration " << i + 1 << std::endl;

    cost = J(m, c, training_data);

    std::cout << "\t m = " << m << " , c = " << c << " , cost = " << cost
              << std::endl;
  }

  return 0;
}
