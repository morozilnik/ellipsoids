#include <chrono>
#include <iostream>
#include <Eigen/Eigen>
#include "Petunin.hpp"
#include "HarPelet.hpp"
#include "Malandain.hpp"

int main()
{
   const int testSize = 1000;
   const int dimensions = 2;
   srand(time(NULL));
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
 
   // Third method: Malandain
   Diploma::Malandain EllipseGen3(inputData);
   begin = std::chrono::steady_clock::now();
   EllipseGen3.calculateEllises(inputData);
   end = std::chrono::steady_clock::now();
   std::cout << "Time difference = " << std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() << std::endl;


   return 0;
}