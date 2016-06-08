#include <algorithm>
#include <chrono>
#include <iostream>
#include <Eigen/Sparse>

#include "Petunin.hpp"

using namespace Eigen;
using namespace std;

namespace Diploma
{


   PetuninEllipse::PetuninEllipse(const PointsVec& input)
   {
   }

   void PetuninEllipse::calculateEllises(const PointsVec& input)
   {
      // 1. Find Diameter
      std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
      double diam = findDiameter(input, D1, D2);
      std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
      std::cout << "Diameter time = " << std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() << std::endl;
      
      cout << "Calculated diameter = " << diam
         << "\nWith D1 = " << D1 << "\nD2 = " << D2 << std::endl;
      HyperCube out;
      findRectangle(input, out);
   }

   const std::pair<Point, Point> PetuninEllipse::getCenter(int i) const
   {
      if (i < F1.size() && F1.size() == F2.size())
      {
         return std::make_pair(F1.col(i), F2.col(i));
      }
      else
      {
         return std::pair<Point, Point>();
      }
   }

   double PetuninEllipse::findDiameter(const PointsVec& input, Point& A, Point& B)
   {
      std::cout << "Inside PetuninEllipse:findDiameter\n";
      double maxDistance = 0;
      for (auto i = 0; i < input.cols(); i++)
      {
         for (auto j = i + 1; j < input.cols(); j++)
         {
            double distance = dist(input.col(i), input.col(j));
            if (distance > maxDistance)
            {
               maxDistance = distance;
               A = input.col(i);
               B = input.col(j);
            }

         }
      }
      return maxDistance;
   }

   void PetuninEllipse::findRectangle(const PointsVec& input, HyperCube output)
   {
      double distance = dist(D1, D2);
      //Change coordinate system, so that D1 D2 is on OX1 axis
      // 1. Translation
      PointsVec points = input.colwise() - D1;
      Point P = D2 - D1;

      // 2. Rotation. Create rotation mat as composition of one-dimentional rotations.
      VectorXd Cos = P / P.norm();

      int dimensions = input.rows();
      MatrixXd rotation = MatrixXd::Identity(dimensions, dimensions);
      
      for (auto i = 1; i < Cos.size(); i++)
      {
         MatrixXd axRot = MatrixXd::Zero(dimensions, dimensions);
         double cosine = Cos(i);
         double sine = sqrt(1 - cosine * cosine);
         if (D2(0) > 0 != D2(i) > 0)
         {
            sine = -sine;
         }

         axRot(0, 0) = cosine;
         axRot(0, i) = sine;
         axRot(i, 0) = -sine;
         axRot(i, i) = cosine;
         rotation *= axRot;
      }

      // 3. Rotate all points on given matrix
      points = rotation * points;

      // 4. Find parallelepiped containing all points
      auto boundingBox = findBoundingBox(points, std::vector<int>());
      // 5. Shrink to square
      points.colwise() -= boundingBox.pointMin;
      Point Zero = VectorXd(dimensions);
      for (int i = 0; i < dimensions; i++)
      {
         double div = (boundingBox.pointMax(i) - boundingBox.pointMin(i));
         points.row(i) /= div;
         Zero(i) = 0.;
      }
      std::vector<double> distances(input.cols());
      
      for (int i = 0; i < distances.size(); i++)
      {
         distances[i] = dist(Zero, input.col(i));
      }
      std::sort(distances.begin(), distances.end());

   }



} //namespace Diploma
