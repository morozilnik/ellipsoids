#include <chrono>
#include <iostream>
#include <Eigen/Eigen>
#include "Petunin.hpp"


int main()
{
   const int testSize = 100000;
   const int dimensions = 5;
   Eigen::MatrixXd inputData = Eigen::MatrixXd::Random(dimensions, testSize) * 1000;
   Diploma::PetuninEllipse EllipseGen1(inputData);
   std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
   EllipseGen1.calculateEllises(inputData);
   std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
   std::cout << "Time difference = " << std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() << std::endl;
   
   Diploma::HarPelet EllipseGen2(inputData);
   begin = std::chrono::steady_clock::now();
   EllipseGen2.calculateEllises(inputData);
   end = std::chrono::steady_clock::now();
   std::cout << "Time difference = " << std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() << std::endl;


   return 0;
}