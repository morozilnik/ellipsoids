#include <algorithm>
#include <iostream>
#include <numeric>
#include "HarPelet.hpp"

namespace Diploma
{
   namespace { //Anonimous

      double M_Measure(FairSplitTreeNode* u, FairSplitTreeNode* v)
      {
         return dist(u->boundingBox.center, v->boundingBox.center) +
            u->boundingBox.radius + v->boundingBox.radius;
      }

   } // namespace anon

   void HarPelet::splitTreeNode(
      const PointsVec& input,
      FairSplitTreeNode* parent)
   {
      FairSplitTreeNode* left = new FairSplitTreeNode();
      FairSplitTreeNode* right = new FairSplitTreeNode();
      parent->left = left;
      parent->right = right;
      left->parent = parent;
      right->parent = parent;
      left->left = 0;
      left->right = 0;
      right->left = 0;
      right->right = 0;
      const BoundingBox& bigBox = parent->boundingBox;

      // Dividing in halves according to biggest dimension of bounding box
      left->indices.reserve(parent->indices.size() / 2);
      right->indices.reserve(parent->indices.size() / 2);
      int j = bigBox.maxDimension;
      double threshold = bigBox.pointMin(j) + bigBox.maxSide / 2;

      for (auto i : parent->indices)
      {
         if (input(j, i) > threshold)
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

   HarPelet::HarPelet(const PointsVec& input)
      : epsilon(0.005)
      , PetuninEllipse(input)
      , inputData(input)
   {
      // We initialize algorithm with single node having all tree
      // PCurr = (root(T), root(T)
      root = new FairSplitTreeNode();
      root->parent = 0;
      root->left = 0;
      root->right = 0;
      root->boundingBox = findBoundingBox(input, root->indices);
      root->indices.resize(input.cols());
      std::iota(root->indices.begin(), root->indices.end(), 0);
      currDiam = root->boundingBox.maxSide;

      NodePair rootPair(root, root);
      //Upper bound for diameter
      rootPair.M = root->boundingBox.radius * 2;
      Pcurr.push_back(rootPair);
   }

   HarPelet::~HarPelet()
   {
      //Deleting all tree nodes
      delete root;
      for (auto pt : FairSplitTree)
      {
         delete pt;
      }
   }

   void HarPelet::pushNewPair(NodePair& uv)
   {
      if (uv.M > (1 + epsilon) * currDiam)
      {
         Pcurr.push_back(uv);
         std::push_heap(begin(Pcurr), end(Pcurr));
         int index1 = rand() % uv.i1->size();
         int index2 = rand() % uv.i2->size();
         double newDistance =
            dist(inputData.col(uv.i1->indices[index1]),
            inputData.col(uv.i2->indices[index2]));
         if (newDistance > currDiam)
         {
            currDiam = newDistance;
            pDiam = inputData.col(uv.i1->indices[index1]);
            qDiam = inputData.col(uv.i2->indices[index2]);
         }
      }
   }

   double HarPelet::findDiameter(const PointsVec& input, Point&A, Point&B)
   {
      std::cout << "Inside HarPelet:findDiameter\n";
      double prevDiam;

      int iteration = 0;
      //Start iterative process
      do
      {
         prevDiam = currDiam;
         NodePair currentPair = Pcurr.front();

         pop_heap(begin(Pcurr), end(Pcurr));
         Pcurr.pop_back();
         // We can throw away a pair if P(v) or P(u) is empty or
         // M(u,v) <= (1 + eps) * currDiam
         if (currentPair.i1->size() < 2 &&
            currentPair.i2->size() < 2 ||
            currentPair.M < (1 + epsilon) * currDiam)
         {
            continue;
         }
         else //When decided not to trow away we split parts
         {
            FairSplitTreeNode* nodeToSplit =
               currentPair.i1->boundingBox.maxSide >= currentPair.i2->boundingBox.maxSide ?
               currentPair.i1 : currentPair.i2;
            if (nodeToSplit->size() > 1)
            {
               if (!nodeToSplit->left || !nodeToSplit->right)
               {
                  splitTreeNode(input, nodeToSplit);
                  FairSplitTree.push_back(nodeToSplit->left);
                  FairSplitTree.push_back(nodeToSplit->right);
               }

               // Create new pairs. Different cases if i1 == i2 and else.
               if (currentPair.i1 == currentPair.i2)
               {
                  NodePair ll(nodeToSplit->left, nodeToSplit->left);
                  ll.M = M_Measure(ll.i1, ll.i2);
                  pushNewPair(ll);
                  NodePair rr(nodeToSplit->right, nodeToSplit->right);
                  rr.M = M_Measure(rr.i1, rr.i2);
                  pushNewPair(rr);
                  NodePair lr(nodeToSplit->left, nodeToSplit->right);
                  lr.M = M_Measure(lr.i1, lr.i2);
                  pushNewPair(lr);
               }
               else
               {
                  FairSplitTreeNode* otherNode =
                     nodeToSplit == currentPair.i1 ?
                     currentPair.i2 : currentPair.i1;
                  NodePair ll(nodeToSplit->left, otherNode);
                  ll.M = M_Measure(ll.i1, ll.i2);
                  pushNewPair(ll);
                  NodePair rr(nodeToSplit->right, otherNode);
                  rr.M = M_Measure(rr.i1, rr.i2);
                  pushNewPair(rr);
               } //endif (i1 == i2)

            } // endif (nodeToSplit.size > 1)
         } // endif (process pair)
      } while (/*currDiam - prevDiam > epsilon &&*/ !Pcurr.empty());

      A = pDiam;
      B = qDiam;
      return currDiam;
   }
} // namespace Diploma