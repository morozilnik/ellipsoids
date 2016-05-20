#pragma once

#include <utility>
#include "helpers.hpp"

namespace Diploma
{
   /*Base class that implements Petunin ellipses generation with diameter search O(n^2)*/
   class PetuninEllipse
   {
   public:
      PetuninEllipse(const PointsVec& input);
      void calculateEllises(const PointsVec& input);
      const std::pair<Point, Point> getCenter(int i) const;

   protected:
      /** Return diameter
      * param[out] A, B - diametrial points
      */
      virtual double findDiameter(const PointsVec& input, Point& A, Point& B);

	   void findRectangle(const PointsVec& input, HyperCube output);

	   Point D1;
	   Point D2;
      PointsVec F1;
      PointsVec F2;
   };

} //namespace Dimploma
