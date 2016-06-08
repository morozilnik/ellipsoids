#include <algorithm>
#include "helpers.hpp"

namespace Diploma
{
   double dist(const Point& A, const Point& B)
   {
      return (A - B).norm();
   }

   BoundingBox findBoundingBox(const PointsVec& input, std::vector<int>& indices)
   {
      BoundingBox R;
      if (indices.empty())
      {
         R.pointMin = input.rowwise().minCoeff();
         R.pointMax = input.rowwise().maxCoeff();
      }
      else
      {
         R.pointMin = R.pointMax = input.col(indices.front());
      }

      for (int i : indices)
      {
         for (int j = 0; j < input.rows(); j++)
         {
            R.pointMin(j) = std::min(R.pointMin(j), input(j, i));
            R.pointMax(j) = std::max(R.pointMax(j), input(j, i));
         }
      }

      Point diff = R.pointMax - R.pointMin;
      R.maxSide = diff.maxCoeff(&R.maxDimension);
      R.center = R.pointMin + (diff / 2);
      R.radius = diff.norm() / 2;
      return R;
   }

} //namespace Diploma
