#include "Malandain.hpp"

namespace Diploma
{

   Malandain::Malandain(const PointsVec& input)
      : PetuninEllipse(input)
      , activeIndices(input.cols(), true)
   {

   }

   double Malandain::findDiameter(const PointsVec& input, Point& A, Point& B)
   {
      double diam = 0;
      bool stop = false;
      pii pq(-1, -1);
      pii res;
      do
      {
         double curDiam = findDoubleNormal(input, pq);
         if (diam < curDiam)
         {
            diam = curDiam;
            res = pq;
            pq.first = getOuterPoint(input, pq);
            if (pq.first < 0)
            {
               stop = true;
            }
         }
         else
         {
            stop = true;
         }
      } while (!stop);


      return diam;
   }

   // Finds pair of points which are furthest from one another
   // Starts with first active index
   double Malandain::findDoubleNormal(const PointsVec& input, pii& DN)
   {
      int p = DN.first;
      if (p < 0)
      {
         p = *std::find(begin(activeIndices), end(activeIndices), true);
      }

      int q;
      double diam = 0;
      double prevDiam;
      do
      {
         prevDiam = diam;
         activeIndices[p] = false;
         double d = findMostDistantPoint(input, p, q);
         if (d > diam)
         {
            diam = d;
            DN = std::make_pair(p, q);
            p = q;
         }
      } while (diam > prevDiam);

      return diam;
   }

   double Malandain::findMostDistantPoint(const PointsVec& input, int from, int& res)
   {
      int q = from + 1;
      if (q == input.cols())
      {
         return -1.;
      }

      double maxDist = -1.;
      for (int i = q; i < input.cols(); i++)
      {
         if (activeIndices[i])
         {
            double d = dist(input.col(from), input.col(i));
            if (d > maxDist)
            {
               maxDist = d;
               q = i;
            }
         }
      }
      res = q;
      return maxDist;
   }

   int Malandain::getOuterPoint(const PointsVec& input, pii pq)
   {
      Point P = input.col(pq.first);
      Point Q = input.col(pq.second);
      Point C = (P + Q) / 2;
      double d = dist(C, P);
      int res = -1;
      for (int i = 0; i < input.cols(); i++)
      {
         if (activeIndices[i] && i != pq.first && i != pq.second
            && dist(C, input.col(i)) > d)
         {
            d = dist(C, input.col(i));
            res = i;
         }
      }
      return res;
   }

} // namespace Diploma