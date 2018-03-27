/**
 *  @file   larpandoracontent/LArObjects/LArAdaBoostDecisionTree.h
 *
 *  @brief  Header file for the lar adaptive boosted decision tree class.
 *
 *  $Log: $
 */
#ifndef LAR_ADABOOST_DECISION_TREE_H
#define LAR_ADABOOST_DECISION_TREE_H 1

#include "Pandora/AlgorithmTool.h"
#include "Pandora/StatusCodes.h"

#include "larpandoracontent/LArHelpers/LArMvaHelper.h"

#include "larpandoracontent/LArObjects/LArMvaInterface.h"

#include <functional>
#include <map>
#include <vector>

namespace lar_content
{

/**
 *  @brief  AdaBoostDecisionTree class
 */
class AdaBoostDecisionTree : public MvaInterface
{
public:
    /**
     *  @brief  Constructor.
     */
    AdaBoostDecisionTree();

    /**
     *  @brief  Copy constructor
     * 
     *  @param  rhs the AdaBoostDecisionTree to copy
     */
    AdaBoostDecisionTree(const AdaBoostDecisionTree &rhs);

    /**
     *  @brief  Assignment operator
     * 
     *  @param  rhs the AdaBoostDecisionTree to assign
     */
    AdaBoostDecisionTree &operator=(const AdaBoostDecisionTree &rhs);

    /**
     *  @brief  Destructor
     */
    ~AdaBoostDecisionTree();

    /**
     *  @brief  Initialize the bdt model 
     *
     *  @param  parameterLocation the location of the model
     *  @param  bdtName the name of the model
     *
     *  @return success
     */
    pandora::StatusCode Initialize(const std::string &parameterLocation, const std::string &bdtName);

    /**
     *  @brief  Classify the set of input features based on the trained model
     *
     *  @param  features the input features
     *
     *  @return the classification 
     */
    bool Classify(const LArMvaHelper::MvaFeatureVector &features) const;

    /**
     *  @brief  Calculate the classification score for a set of input features, based on the trained model
     *
     *  @param  features the input features
     *
     *  @return the classification score
     */
    double CalculateClassificationScore(const LArMvaHelper::MvaFeatureVector &features) const;

    /**
     *  @brief  Calculate the classification probability for a set of input features, based on the trained model
     *
     *  @param  features the input features
     *
     *  @return the classification probability
     */
    double CalculateProbability(const LArMvaHelper::MvaFeatureVector &features) const;

private:
    /**
     *  @brief Node class used for representing a decision tree
     */
    class Node
    {
    public:
        /**
         *  @brief  Constructor
         *  
         *  @param  nodeId id of node
         *  @param  parentNodeId id of parent node
         *  @param  leftChildNodeId id of left child node
         *  @param  rightChildNodeId id of right child node
         *  @param  isLeaf is node a leaf
         *  @param  threshold threshold for cut
         *  @param  variableId index
         *  @param  outcome, only written if leaf node 
         */
        Node(const int nodeId, const int parentNodeId, const int leftChildNodeId, const int rightChildNodeId, const bool isLeaf, const double threshold, const int variableId, const bool outcome = false);

        /**
         *  @brief  Constructor using xml handle to set member variables
         *
         *  @param  pXmlHandle xml handle to use when setting member variables 
         */
        Node(const pandora::TiXmlHandle *const pXmlHandle);

        /**
         *  @brief  Copy constructor
         * 
         *  @param  rhs the node to copy
         */
        Node(const Node &rhs);

        /**
         *  @brief  Assignment operator
         * 
         *  @param  rhs the node to assign
         */
        Node &operator=(const Node &rhs);

        /**
         *  @brief  Destructor
         */
        ~Node();

        /**
         *  @brief  Return node id
         *
         *  @return node id
         */
        int GetNodeId() const;

        /**
         *  @brief  Return parent node id
         *
         *  @return parent node id
         */
        int GetParentNodeId() const;

        /**
         *  @brief  Return left child node id
         *
         *  @return left child node id
         */
        int GetLeftChildNodeId() const;

        /**
         *  @brief  Return right child node id
         *
         *  @return right child node id
         */
        int GetRightChildNodeId() const;

        /**
         *  @brief  Return is the node a leaf
         *
         *  @return is node a leaf
         */
        bool IsLeaf() const;

        /**
         *  @brief  Return node threshold
         *
         *  @return threshold cut
         */
        double GetThreshold() const;

        /**
         *  @brief  Return cut variable
         *
         *  @return variable cut on
         */
        int GetVariableId() const;

        /**
         *  @brief  Return outcome
         *
         *  @return outcome of cut
         */
        bool GetOutcome() const;

    private:
        int       m_nodeId;               ///< Node id
        int       m_parentNodeId;         ///< Parent node id
        int       m_leftChildNodeId;      ///< Left child node id
        int       m_rightChildNodeId;     ///< Right child node id
        bool      m_isLeaf;               ///< Is node a leaf
        double    m_threshold;            ///< Threshold used for decision if decision node
        int       m_variableId;           ///< Variable cut on for decision if decision node
        bool      m_outcome;              ///< Outcome if leaf node
    };

    typedef std::map<int, const Node*> IdToNodeMap;

    /**
     *  @brief  WeakClassifier class containing a decision tree and a weight
     */
    class WeakClassifier
    {
    public: 
        /**
         *  @brief  Constructor, set hierarchy for nodes 
         *
         *  @param  idToNodeMap map of node id to node
         *  @param  weight applied to decision tree
         *  @param  treeId index of tree in decision tree forest 
         */
        WeakClassifier(const IdToNodeMap &idToNodeMap, const double weight, const int treeId);

        /**
         *  @brief  Constructor using xml handle to set member variables
         *
         *  @param  pXmlHandle xml handle to use when setting member variables 
         */
        WeakClassifier(const pandora::TiXmlHandle *const pXmlHandle);

        /**
         *  @brief  Copy constructor
         * 
         *  @param  rhs the weak classifier to copy
         */
        WeakClassifier(const WeakClassifier &rhs);

        /**
         *  @brief  Assignment operator
         * 
         *  @param  rhs the weak classifier to assign
         */
        WeakClassifier &operator=(const WeakClassifier &rhs);

        /**
         *  @brief  Destructor
         */
        ~WeakClassifier();

        /**
         *  @brief  Predict signal or background based on trained data
         *
         *  @param  features the input features 
         *
         *  @return is signal or background
         */
        bool Predict(const LArMvaHelper::MvaFeatureVector &features) const;

        /**
         *  @brief  Evalute node and return outcome
         *
         *  @param  nodeId current node id 
         *  @param  features the input features 
         *
         *  @return is signal or background node
         */
        bool EvaluateNode(const int nodeId, const LArMvaHelper::MvaFeatureVector &features) const;

        /**
         *  @brief  Get boost weight for weak classifier
         *
         *  @return weight for decision tree
         */
        double GetWeight() const;

        /**
         *  @brief  Get tree id for weak classifier
         *
         *  @return tree id 
         */
        int GetTreeId() const;

    private:
        IdToNodeMap     m_idToNodeMap; ///< Decision tree nodes
        double          m_weight;      ///< Boost weight 
        int             m_treeId;      ///< Decision tree id        
    };

    typedef std::vector<const WeakClassifier*> WeakClassifiers;

    /**
     *  @brief  StrongClassifier class used in application of adaptive boost decision tree
     */
    class StrongClassifier
    {
    public:
        /**
         *  @brief  Constructor
         */
        StrongClassifier(const WeakClassifiers &weakClassifiers);

        /**
         *  @brief  Constructor using xml handle to set member variables
         *
         *  @param  pXmlHandle xml handle to use when setting member variables 
         */
        StrongClassifier(const pandora::TiXmlHandle *const pXmlHandle);

        /**
         *  @brief  Copy constructor
         * 
         *  @param  rhs the strong classifier to copy
         */
        StrongClassifier(const StrongClassifier &rhs);

        /**
         *  @brief  Assignment operator
         * 
         *  @param  rhs the strong classifier to assign
         */
        StrongClassifier &operator=(const StrongClassifier &rhs);

        /**
         *  @brief  Destructor
         */
        ~StrongClassifier();

        /**
         *  @brief  Predict signal or background based on trained data
         *
         *  @param  features the input features 
         *
         *  @return return score produced from trained model
         */
        double Predict(const LArMvaHelper::MvaFeatureVector &features) const;

    private: 
        /**
         *  @brief  Read xml element and if weak classifier add to member variables
         */
        pandora::StatusCode ReadComponent(pandora::TiXmlElement *pCurrentXmlElement);

        WeakClassifiers     m_weakClassifiers;     ///< Vector of weak classifers 
    };

    /**
     *  @brief  Calculate score for input features using strong classifier
     *
     *  @param  features the input features
     *
     *  @return score
     */
    double CalculateScore(const LArMvaHelper::MvaFeatureVector &features) const;

    StrongClassifier     *m_pStrongClassifier;           ///< Strong adaptive boost tree classifier 
};

//------------------------------------------------------------------------------------------------------------------------------------------

inline int AdaBoostDecisionTree::Node::GetNodeId() const 
{
    return m_nodeId;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline int AdaBoostDecisionTree::Node::GetParentNodeId() const 
{
    return m_parentNodeId;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline int AdaBoostDecisionTree::Node::GetLeftChildNodeId() const 
{
    return m_leftChildNodeId;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline int AdaBoostDecisionTree::Node::GetRightChildNodeId() const 
{
    return m_rightChildNodeId;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline bool AdaBoostDecisionTree::Node::IsLeaf() const
{
    return m_isLeaf;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline double AdaBoostDecisionTree::Node::GetThreshold() const
{
    return m_threshold;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline int AdaBoostDecisionTree::Node::GetVariableId() const
{
    return m_variableId;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline bool AdaBoostDecisionTree::Node::GetOutcome() const
{
    return m_outcome;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline double AdaBoostDecisionTree::WeakClassifier::GetWeight() const 
{
    return m_weight;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline int AdaBoostDecisionTree::WeakClassifier::GetTreeId() const
{
    return m_treeId;
}

} // namespace lar_content

#endif // #ifndef LAR_ADABOOST_DECISION_TREE_H
