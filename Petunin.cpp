#include <algorithm>
#include <numeric>
#include <Eigen/Sparse>

#include "Petunin.hpp"

using namespace Eigen;
using namespace std;

namespace Dimploma
{
   PetuninEllipse::PetuninEllipse(const PointsVec& input)
   {
   }

   void PetuninEllipse::calculateEllises(const PointsVec& input)
   {
      // 1. Find Diameter
      double diam = findDiameter(input, D1, D2);


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
      double maxDistance = 0;
      for (auto i = 0; i < input.cols(); i++)
      {
         for (auto j = i + 1; j < input.cols(); j++)
         {
            if (i != j)
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
      SparseMatrix<double> rotation(points.cols(), points.rows());
      
      for (auto i = 1; i < Cos.size(); i++)
      {

         std::vector<Triplet<double> > mat(4);
         double sin = std::sqrt(1 - Cos(i) * Cos(i));
         mat[0] = Triplet<double>(0, 0, Cos(i));
         mat[0] = Triplet<double>(0, i, Cos(i));
      }
   }

   struct BoundingBox
   {
      double maxSide;
      int maxDimension;
      VectorXd pointMin;
      VectorXd pointMax;
      VectorXd center;
      double radius;
   };

   struct NodePair
   {
      NodePair(int a, int b)
      : i1(a), i2(b)
      {}
      int i1;
      int i2;
      double M;
      const bool operator< (const NodePair& r) const
      {
         return M < r.M;
      }
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
            R.pointMin(j) = std::min(R.pointMin(j), input(i, j));
            R.pointMax(j) = std::max(R.pointMax(j), input(i, j));
         }
      }

      auto diff = R.pointMax - R.pointMin;

      R.maxSide = diff.minCoeff(&R.maxDimension);
      R.center = R.pointMin + diff / 2;
      R.radius = (diff / 2).norm();
      return R;
   }

   void splitTreeNode(
      const PointsVec& input,
      FairSplitTreeNode* parent)
   {
      FairSplitTreeNode* left = parent->left;
      FairSplitTreeNode* right = parent->right;
      left->parent = parent;
      right->parent = parent;
      const BoundingBox& bigBox = parent->boundingBox;
      
      // Dividing in halves according to biggest dimension of bounding box
      left->indices.reserve(parent->indices.size() / 2);
      right->indices.reserve(parent->indices.size() / 2);
      int j = bigBox.maxDimension;
      double threshold = bigBox.center(j);

      for (auto i : parent->indices)
      {
         if (input(i, j) > threshold)
         {
            right->indices.push_back(i);
         }
         else
         {
            left->indices.push_back(i);
         }
      }
      left->boundingBox = findBoundingBox(input, left->indices);
      right->boundingBox = findBoundingBox(input, right->indices);
      
   }

   double M_Measure(FairSplitTreeNode* u, FairSplitTreeNode* v)
   {
      return dist(u->boundingBox.center, v->boundingBox.center) +
         u->boundingBox.radius + v->boundingBox.radius;
   }

   double findDiameterHarPelet(const PointsVec& input, Point&A, Point&B)
   {
      // We initialize algorithm with single node having all tree
      // PCurr = (root(T), root(T)

      FairSplitTreeNode root;
      root.parent = 0;
      root.boundingBox = findBoundingBox(input, root.indices);
      root.indices.resize(input.cols());
      std::iota(root.indices.begin(), root.indices.end(), 0);
      std::vector<FairSplitTreeNode *> FairSplitTree(1);
      FairSplitTree[0] = &root;

      
      // Starting approximation as max side of bounding box (sqrt(d) approximation)
      VectorXd pDiam, qDiam;
      double currDiam = root.boundingBox.maxSide;
      // Heap with nodes pairs
      std::vector<NodePair> Pcurr;
      NodePair rootPair(0, 0);
      //Upper bound for diameter
      rootPair.M = root.boundingBox.radius * 2;
      Pcurr.push_back(rootPair);
      
      // The correct value of diameter is searched within epsilon
      double epsilon = 0.1;
      double prevDiam;

      //Start iterative process
      do
      {
         prevDiam = currDiam;
         NodePair& currentPair = Pcurr.front();

         // We can throw away a pair if P(v) or P(u) is empty or 
         // M(u,v) <= (1 + eps) * currDiam
         if (FairSplitTree[currentPair.i1]->empty() || 
             FairSplitTree[currentPair.i2]->empty() ||
             currentPair.M < (1 + epsilon) * currDiam)
         {
            pop_heap(begin(Pcurr), end(Pcurr));
            Pcurr.pop_back();
            continue;
         }
         else //When decided not to trow away we split parts
         {
            FairSplitTreeNode* nodeToSplit =
               FairSplitTree[currentPair.i1]->size() >= FairSplitTree[currentPair.i2]->size() ?
               FairSplitTree[currentPair.i1] : FairSplitTree[currentPair.i2];
            if (nodeToSplit->size() > 1)
            {
               if (!nodeToSplit->left || !nodeToSplit->right)
               {
                  nodeToSplit->left = new FairSplitTreeNode();
                  nodeToSplit->right = new FairSplitTreeNode();
               }
               splitTreeNode(input, nodeToSplit);
            }
            
         }
         

         
      } while (currDiam - prevDiam > epsilon && !Pcurr.empty());
            
      
      return currDiam;
   }



} //namespace Dimploma
