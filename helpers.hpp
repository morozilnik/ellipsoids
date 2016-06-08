#pragma once
/** This is geometrical types for Moroz Anton Diploma */


#include <vector>
#include <cmath>

#include <Eigen/Eigen>

namespace Diploma
{
	typedef Eigen::VectorXd Point;
	typedef Eigen::MatrixXd PointsVec;
	typedef PointsVec HyperCube;


   struct BoundingBox
   {
      double maxSide;
      int maxDimension;
      Eigen::VectorXd pointMin;
      Eigen::VectorXd pointMax;
      Eigen::VectorXd center;
      double radius;
   };

   struct FairSplitTreeNode{
   public: // methods
      bool empty()
      {
         return indices.empty();
      }

      std::size_t size()
      {
         return indices.size();
      }
   public: // fields
      // Enumeration of nodes in tree
      std::vector<int> indices;
      FairSplitTreeNode* parent;
      FairSplitTreeNode* left;
      FairSplitTreeNode* right;
      BoundingBox boundingBox;
   };

   struct NodePair
   {
      NodePair(FairSplitTreeNode* a, FairSplitTreeNode* b)
      : i1(a), i2(b)
      {}
      FairSplitTreeNode* i1;
      FairSplitTreeNode* i2;
      double M;
      const bool operator< (const NodePair& r) const
      {
         return M < r.M;
      }
   };

   double dist(const Point& A, const Point& B);
   BoundingBox findBoundingBox(const PointsVec& input, std::vector<int>& indices);


} //namespace Diploma
