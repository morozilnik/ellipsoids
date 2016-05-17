#pragma once

#include <utility>
#include "types.hpp"

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

   struct BoundingBox;
   struct FairSplitTreeNode;
   struct NodePair;

   /*Class that implemets diameter search with HarPelet algorithm with epsilon accuracy*/
   class HarPelet : public PetuninEllipse
   {
   public:
      explicit HarPelet(const PointsVec& input);
      ~HarPelet();

   protected:
      virtual double findDiameter(const PointsVec& input, Point& A, Point&B);

   private:
      void splitTreeNode(const PointsVec& input, FairSplitTreeNode* parent);
      void pushNewPair(NodePair& uv);

      FairSplitTreeNode* root;
      Point pDiam;
      Point qDiam;
      double currDiam;
      // The correct value of diameter is searched within epsilon
      const double epsilon;

      // Heap with nodes pairs
      std::vector<NodePair> Pcurr;
      //This keeps refernces to nodes, so we can delete them later
      std::vector<FairSplitTreeNode *> FairSplitTree;
      const PointsVec& inputData;


   };

} //namespace Dimploma
