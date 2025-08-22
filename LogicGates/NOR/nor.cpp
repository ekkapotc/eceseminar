#include <iostream>
#include <vector>
#include <tuple>
#include <ctime>
#include <cstdlib>
#include <cmath>

#include "tls.hpp"

/*
* In this problem, 
* we will determine the two parameters m and c of a given linear model 'y = mx + c'
*/

#define value(x) 	(x).v
#define derivative(x) 	(x).t

template<typename T>
using Tuple = std::tuple<T,T,T>;

template<typename T>
using Vector = std::vector<T>;
   
template <typename T>
T sigmoid( const T & x ){
    return 1.0/(1.0+exp(-x));
}

template<typename T>
T neuralNet(T x1 , T x2 , T w11, T w12, T w21 , T w22 , T b1 , T b2 , T w3 , T w4 , T b3){
   T z1 = sigmoid(w11*x1+w21*x2+b1);
   T z2 = sigmoid(w12*x1+w22*x2+b2);
   T z3 = sigmoid(w3*z1+w4*z2+b3);
   return z3;
}

template<typename T>
T cost( T w11 , T w12 , T w21 , T w22 , 
          T b1  , T b2 ,
          T w3  , T w4 , 
          T b3  , 
          Vector<Tuple<double>> & data ){

          T acc {0.0};

	  for( auto & e : data ){
	     T x1 = std::get<0>(e);
	     T x2 = std::get<1>(e);
	     T y_actual  = std::get<2>(e);
	     T y_expected = neuralNet(x1,x2,w11,w12,w21,w22,b1,b2,w3,w4,b3);
	     acc = acc + (y_expected-y_actual)*(y_expected-y_actual);
	  }

          acc = acc / data.size();
          return acc;
}

int main(int argc,char * argv[]){
  
  Vector<Tuple<double>> training_data = {
	{0.0 , 0.0 , 1.0},
	{1.0 , 0.0 , 0.0},
	{0.0 , 1.0 , 0.0},
	{1.0 , 1.0 , 0.0}
  };

  auto no_samples = training_data.size();

  srand(time(nullptr));
  
  Tangent w11 = static_cast<double>(rand()) /RAND_MAX;
  Tangent w12 = static_cast<double>(rand()) /RAND_MAX;
  Tangent w21 = static_cast<double>(rand()) /RAND_MAX;
  Tangent w22 = static_cast<double>(rand()) /RAND_MAX;
  Tangent b1 = static_cast<double>(rand()) /RAND_MAX;
  Tangent b2 = static_cast<double>(rand()) /RAND_MAX;
  Tangent w3 = static_cast<double>(rand()) /RAND_MAX;
  Tangent w4 = static_cast<double>(rand()) /RAND_MAX;
  Tangent b3 = static_cast<double>(rand()) /RAND_MAX;
 
  std::cout << "Pre-trained Model: " << std::endl;
  //use pre-trained model
  for ( auto & e  : training_data ){
    Tangent x1 = std::get<0>(e);
    Tangent x2 = std::get<1>(e);
    Tangent z = neuralNet(x1,x2,w11,w12,w21,w22,b1,b2,w3,w4,b3); 
    std::cout << value(x1) << " NOR " << value(x2) << " = " << value(z) << std::endl; 
  }
 
  double step{1e-2};

  for( auto i{0} ; i<500000 ; i++ ){
      derivative(w11) = 1.0;
      double dJdw11 = derivative(cost(w11,w12,w21,w22,b1,b2,w3,w4,b3,training_data));
      derivative(w11) = 0.0;
      
      derivative(w12) = 1.0;
      double dJdw12 = derivative(cost(w11,w12,w21,w22,b1,b2,w3,w4,b3,training_data));
      derivative(w12) = 0.0;
      
      derivative(w21) = 1.0;
      double dJdw21 = derivative(cost(w11,w12,w21,w22,b1,b2,w3,w4,b3,training_data));
      derivative(w21) = 0.0;

      derivative(w22) = 1.0;
      double dJdw22 = derivative(cost(w11,w12,w21,w22,b1,b2,w3,w4,b3,training_data));
      derivative(w22) = 0.0;

      derivative(b1) = 1.0;
      double dJdb1 = derivative(cost(w11,w12,w21,w22,b1,b2,w3,w4,b3,training_data));
      derivative(b1) = 0.0;
      
      derivative(b2) = 1.0;
      double dJdb2 = derivative(cost(w11,w12,w21,w22,b1,b2,w3,w4,b3,training_data));
      derivative(b2) = 0.0;

      derivative(w3) = 1.0;
      double dJdw3 = derivative(cost(w11,w12,w21,w22,b1,b2,w3,w4,b3,training_data));
      derivative(w3) = 0.0;

      derivative(w4) = 1.0;
      double dJdw4 = derivative(cost(w11,w12,w21,w22,b1,b2,w3,w4,b3,training_data));
      derivative(w4) = 0.0;

      derivative(b3) = 1.0;
      double dJdb3 = derivative(cost(w11,w12,w21,w22,b1,b2,w3,w4,b3,training_data));
      derivative(b3) = 0.0;


      //update parameters
      w11 = w11 - step*dJdw11;
      w12 = w12 - step*dJdw12;
      w21 = w21 - step*dJdw21;
      w22 = w22 - step*dJdw22;
      b1  = b1  - step*dJdb1;
      b2  = b2  - step*dJdb2;
      w3  = w3  - step*dJdw3;
      w4  = w4  - step*dJdw4;
      b3  = b3  - step*dJdb3;   
  }

  std::cout << "Post-trained Model: " << std::endl;
  //use trained model
  for ( auto & e  : training_data ){
    Tangent x1 = std::get<0>(e);
    Tangent x2 = std::get<1>(e);
    Tangent z = neuralNet(x1,x2,w11,w12,w21,w22,b1,b2,w3,w4,b3); 
    std::cout << value(x1) << " NOR " << value(x2) << " = " << value(z) << std::endl; 
 }

  return 0;
}
