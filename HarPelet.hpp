#include "Petunin.hpp"

namespace Diploma
{

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
} // namespace Diploma