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
#include "larpandoracontent/LArObjects/LArMultivariateAnalyisBaseClass.h"

#include <functional>
#include <map>
#include <vector>

namespace lar_content
{

/**
 *  @brief  AdaBoostDecisionTree class
 */
class AdaBoostDecisionTree : public MultivariateAnalyisBaseClass
{
public:
    /**
     *  @brief  Constructor.
     */
    AdaBoostDecisionTree();

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
    bool Classify(const DoubleVector &features) const;

    /**
     *  @brief  Calculate the classification score for a set of input features, based on the trained model
     *
     *  @param  features the input features
     *
     *  @return the classification score
     */
    double CalculateClassificationScore(const DoubleVector &features) const;

    /**
     *  @brief  Calculate the classification probability for a set of input features, based on the trained model
     *
     *  @param  features the input features
     *
     *  @return the classification probability
     */
    double CalculateProbability(const DoubleVector &features) const;

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
         *  @brief  Return node id
         */
        int GetNodeId() const;

        /**
         *  @brief  Return parent node id
         */
        int GetParentNodeId() const;

        /**
         *  @brief  Return left child node id
         */
        int GetLeftChildNodeId() const;

        /**
         *  @brief  Return right child node id
         */
        int GetRightChildNodeId() const;

        /**
         *  @brief  Return is the node a leaf
         */
        bool IsLeaf() const;

        /**
         *  @brief  Return node threshold
         */
        double GetThreshold() const;

        /**
         *  @brief  Return cut variable
         */
        int GetVariableId() const;

        /**
         *  @brief  Return outcome
         */
        bool GetOutcome() const;

    private:
        const int       m_nodeId;               ///< Node id
        const int       m_parentNodeId;         ///< Parent node id
        const int       m_leftChildNodeId;      ///< Left child node id
        const int       m_rightChildNodeId;     ///< Right child node id
        const bool      m_isLeaf;               ///< Is node a leaf
        const double    m_threshold;            ///< Threshold used for decision if decision node
        const int       m_variableId;           ///< Variable cut on for decision if decision node
        const bool      m_outcome;              ///< Outcome if leaf node
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
         */
        WeakClassifier(const IdToNodeMap &idToNodeMap, const double weight, const int treeId);

        /**
         *  @brief  Predict signal or background based on trained data
         *
         *  @param  features the input features 
         */
        bool Predict(const DoubleVector &features) const;

        /**
         *  @brief  Evalute node and return outcome
         *
         *  @param  nodeId current node id 
         *  @param  features the input features 
         */
        bool EvaluateNode(const int nodeId, const DoubleVector &features) const;

        /**
         *  @brief  Get boost weight for weak classifier
         */
        double GetWeight() const;

        /**
         *  @brief  Get tree id for weak classifier
         */
        int GetTreeId() const;

    private:
        IdToNodeMap     m_idToNodeMap; ///< Decision tree nodes
        const double    m_weight;      ///< Boost weight 
        const int       m_treeId;      ///< Decision tree id        
    };

    typedef std::vector<WeakClassifier> WeakClassifiers;

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
         *  @brief  Predict signal or background based on trained data
         *
         *  @param  features the input features 
         */
        double Predict(const DoubleVector &features) const;

    private: 
        WeakClassifiers     m_weakClassifiers;     ///< Vector of weak classifers 
    };

    /**
     *  @brief  Read the bdt parameters from an xml file
     *
     *  @param  bdtFileName the sml file name
     *  @param  bdtName the name of the bdt
     *  @param  weakClassifiers vector of weak classifiers contained in xml
     *
     *  @return success
     */
    pandora::StatusCode ReadXmlFile(const std::string &bdtFileName, const std::string &bdtName, WeakClassifiers &weakClassifiers) const;

    /**
     *  @brief  Check whether the strong classifier has been correctly initialized
     */
    void CheckInitialization() const;

    /**
     *  @brief  Read the component at the current xml element
     *
     *  @param  pCurrentXmlElement address of the current xml element
     *  @param  weakClassifiers vector of weak classifiers contained in xml
     *
     *  @return success
     */
    pandora::StatusCode ReadComponent(pandora::TiXmlElement *pCurrentXmlElement, WeakClassifiers &weakClassifiers) const;

    /**
     *  @brief  Read the decision tree component at the current xml handle
     *
     *  @param  currentHandle the current xml handle
     *  @param  weakClassifiers vector of weak classifiers contained in xml
     *
     *  @return success
     */
    pandora::StatusCode ReadDecisionTree(pandora::TiXmlHandle &currentHandle, WeakClassifiers &weakClassifiers) const;

    /**
     *  @brief  Read the node component at the current xml handle
     *
     *  @param  currentHandle the current xml handle
     *  @param  idToNodeMap map of id to node
     *
     *  @return success
     */
    pandora::StatusCode ReadNode(pandora::TiXmlHandle &currentHandle, IdToNodeMap &idToNodeMap) const;

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
