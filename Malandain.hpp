#include "Petunin.hpp"

namespace Diploma
{
   typedef std::pair<int, int> pii;
   class Malandain : public PetuninEllipse
   {
   public:
      Malandain(const PointsVec& input);
   protected:
      /** Return diameter
      * param[out] A, B - diametrial points
      */
      virtual double findDiameter(const PointsVec& input, Point& A, Point& B);
   private:
      // Finds pair of points which are furthest from one another
      // Starts with first active index
      double findDoubleNormal(const PointsVec& input, pii& DN);
      double findMostDistantPoint(const PointsVec& input, int from, int& res);
      int getOuterPoint(const PointsVec& input, pii pq);
      std::vector <bool> activeIndices;
   };
} //namespace Diploma