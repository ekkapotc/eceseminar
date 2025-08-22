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

double dJ_dm(double m, double c, Vector<Tuple<double>> &data) {
  double deriv{0.0};
  for (auto &e : data) {
    auto x = std::get<0>(e);
    auto y_exp = std::get<1>(e);
    auto y_act = m * x + c;
    deriv += (y_exp - y_act) * (-x);
  }
  deriv *= 2.0;
  deriv /= data.size();

  return deriv;
}

double dJ_dc(double m, double c, Vector<Tuple<double>> &data) {
  double deriv{0.0};
  for (auto &e : data) {
    auto x = std::get<0>(e);
    auto y_exp = std::get<1>(e);
    auto y_act = m * x + c;
    deriv += (y_exp - y_act) * (-1.0);
  }
  deriv *= 2.0;
  deriv /= data.size();

  return deriv;
}

int main(int argc, char *argv[]) {

  size_t constexpr MAX_EPOCHS{50};

  Vector<Tuple<double>> training_data = {
      {0.0, 1.0}, {1.0, 3.0}, {2.0, 5.0}, {3.0, 7.0}};

  auto no_samples = training_data.size();

  srand(time(nullptr));

  double m{static_cast<double>(rand()) / RAND_MAX};
  double c{static_cast<double>(rand()) / RAND_MAX};

  double cost{J(m, c, training_data)};

  std::cout << "Epoch " << 0 << std::endl;
  std::cout << "\t m = " << m << " , c = " << c
            << " , cost = " << J(m, c, training_data) << std::endl;

  double step{1e-2};

  for (auto i{0}; i < MAX_EPOCHS; i++) {

    double dJdm = dJ_dm(m, c, training_data);
    double dJdc = dJ_dc(m, c, training_data);

    m -= step * dJdm;
    c -= step * dJdc;

    std::cout << "Epoch " << i + 1 << std::endl;
    std::cout << "\t m = " << m << " , c = " << c
              << " , cost = " << J(m, c, training_data) << std::endl;
  }

  return 0;
}
