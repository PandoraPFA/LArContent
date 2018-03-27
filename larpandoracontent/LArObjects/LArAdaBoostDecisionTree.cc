/**
 *  @file   larpandoracontent/LArObjects/LArAdaBoostDecisionTree.cc
 *
 *  @brief  Implementation of the lar adaptive boost decision tree class.
 *
 *  $Log: $
 */

#include "Helpers/XmlHelper.h"

#include "larpandoracontent/LArObjects/LArAdaBoostDecisionTree.h"

namespace lar_content
{

AdaBoostDecisionTree::AdaBoostDecisionTree() :
    m_pStrongClassifier(nullptr)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

AdaBoostDecisionTree::AdaBoostDecisionTree(const AdaBoostDecisionTree &rhs)
{
    m_pStrongClassifier = new StrongClassifier(*rhs.m_pStrongClassifier);
}

//------------------------------------------------------------------------------------------------------------------------------------------

AdaBoostDecisionTree &AdaBoostDecisionTree::operator=(const AdaBoostDecisionTree &rhs)
{
    if (this != &rhs)
    {
        m_pStrongClassifier = new StrongClassifier(*rhs.m_pStrongClassifier);
    }

    return *this;
}

//------------------------------------------------------------------------------------------------------------------------------------------

AdaBoostDecisionTree::~AdaBoostDecisionTree()
{
    delete m_pStrongClassifier;
}

//------------------------------------------------------------------------------------------------------------------------------------------

pandora::StatusCode AdaBoostDecisionTree::Initialize(const std::string &bdtXmlFileName, const std::string &bdtName)
{
    if (m_pStrongClassifier)
    {
        std::cout << "AdaBoostDecisionTree: AdaBoostDecisionTree was already initialized" << std::endl;
        return pandora::STATUS_CODE_FAILURE;
    }

    pandora::TiXmlDocument xmlDocument(bdtXmlFileName);

    if (!xmlDocument.LoadFile())
    {
        std::cout << "AdaBoostDecisionTree::Initialize - Invalid xml file." << std::endl;
        return pandora::STATUS_CODE_FAILURE;
    }

    const pandora::TiXmlHandle xmlDocumentHandle(&xmlDocument);
    pandora::TiXmlNode *pContainerXmlNode(pandora::TiXmlHandle(xmlDocumentHandle).FirstChildElement().Element());

    while (pContainerXmlNode)
    {
        if (pContainerXmlNode->ValueStr() != "AdaBoostDecisionTree")
            return pandora::STATUS_CODE_FAILURE;

        const pandora::TiXmlHandle currentHandle(pContainerXmlNode);

        std::string currentName;
        PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, pandora::XmlHelper::ReadValue(currentHandle, "Name", currentName));

        if (currentName.empty() || (currentName.size() > 1000))
        {
            std::cout << "AdaBoostDecisionTree::Initialize - Implausible AdaBoostDecisionTree name extracted from xml." << std::endl;
            return pandora::STATUS_CODE_INVALID_PARAMETER;
        }

        if (currentName == bdtName)
            break;

        pContainerXmlNode = pContainerXmlNode->NextSibling();
    }

    if (!pContainerXmlNode)
    {
        std::cout << "AdaBoostDecisionTree: Could not find an AdaBoostDecisionTree of name " << bdtName << std::endl;
        return pandora::STATUS_CODE_NOT_FOUND;
    }

    pandora::TiXmlHandle xmlHandle(pContainerXmlNode);

    try
    { 
        m_pStrongClassifier = new StrongClassifier(&xmlHandle);
    }
    catch (pandora::StatusCodeException &statusCodeException)
    {
        if (pandora::STATUS_CODE_INVALID_PARAMETER == statusCodeException.GetStatusCode())
        {
            std::cout << "AdaBoostDecisionTree: Initialization failure, unknown component in xml file." << std::endl;
            delete m_pStrongClassifier;            
        }
        else if (pandora::STATUS_CODE_FAILURE == statusCodeException.GetStatusCode())
        {
            std::cout << "AdaBoostDecisionTree: Node definition does not contain expected leaf or branch variables." << std::endl;
        }
    }

    return pandora::STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool AdaBoostDecisionTree::Classify(const LArMvaHelper::MvaFeatureVector &features) const
{
    return (this->CalculateScore(features) > 0.0 ? true : false);
}

//------------------------------------------------------------------------------------------------------------------------------------------

double AdaBoostDecisionTree::CalculateClassificationScore(const LArMvaHelper::MvaFeatureVector &features) const
{
    return this->CalculateScore(features);
}

//------------------------------------------------------------------------------------------------------------------------------------------

double AdaBoostDecisionTree::CalculateProbability(const LArMvaHelper::MvaFeatureVector &features) const
{
    // ATTN: BDT score, once normalised by total weight, is confined to the range -1 to +1.  This linear mapping places the score in the 
    // range 0 to 1 so that it may be interpreted as a probability. 
    return (this->CalculateScore(features) + 1.0) * 0.5;
}

//------------------------------------------------------------------------------------------------------------------------------------------

double AdaBoostDecisionTree::CalculateScore(const LArMvaHelper::MvaFeatureVector &features) const
{
    if (!m_pStrongClassifier)
    {
        std::cout << "AdaBoostDecisionTree: Attempting to use an uninitialized bdt" << std::endl;
        throw pandora::StatusCodeException(pandora::STATUS_CODE_NOT_INITIALIZED);
    }

    double score(0.0);

    try
    {
        // ATTN: Add consistency check for number of features, bearing in min not all features in a bdt may be used
        score = m_pStrongClassifier->Predict(features);
    }
    catch (pandora::StatusCodeException &statusCodeException)
    {
        if (pandora::STATUS_CODE_NOT_FOUND == statusCodeException.GetStatusCode())
        {
            std::cout << "AdaBoostDecisionTree: Caught exception thrown when trying to cut on an unknown variable." << std::endl;
        }
        else if (pandora::STATUS_CODE_INVALID_PARAMETER == statusCodeException.GetStatusCode())
        {
            std::cout << "AdaBoostDecisionTree: Caught exception thrown when classifier weights sum to zero indicating defunct classifier." << std::endl;
        }
        else if (pandora::STATUS_CODE_OUT_OF_RANGE == statusCodeException.GetStatusCode())
        {
            std::cout << "AdaBoostDecisionTree: Caught exception thrown when heirarchy in decision tree is incomplete." << std::endl;
        }
        else
        {
            std::cout << "AdaBoostDecisionTree: Unexpected exception thrown." << std::endl;
        }
        throw statusCodeException;
    }

    return score;
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

AdaBoostDecisionTree::Node::Node(const int nodeId, const int parentNodeId, const int leftChildNodeId, const int rightChildNodeId, const bool isLeaf, const double threshold, const int variableId, const bool outcome) : 
    m_nodeId(nodeId),
    m_parentNodeId(parentNodeId),
    m_leftChildNodeId(leftChildNodeId),
    m_rightChildNodeId(rightChildNodeId),
    m_isLeaf(isLeaf),
    m_threshold(threshold),
    m_variableId(variableId),
    m_outcome(outcome)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

AdaBoostDecisionTree::Node::Node(const pandora::TiXmlHandle *const pXmlHandle) 
{
    PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, pandora::XmlHelper::ReadValue(*pXmlHandle, "NodeId", m_nodeId));
    PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, pandora::XmlHelper::ReadValue(*pXmlHandle, "ParentNodeId", m_parentNodeId));

    pandora::StatusCode leftChildNodeIdStatusCode(pandora::XmlHelper::ReadValue(*pXmlHandle, "LeftChildNodeId", m_leftChildNodeId));
    pandora::StatusCode rightChildNodeIdStatusCode(pandora::XmlHelper::ReadValue(*pXmlHandle, "RightChildNodeId", m_rightChildNodeId));
    pandora::StatusCode thresholdStatusCode(pandora::XmlHelper::ReadValue(*pXmlHandle, "Threshold", m_threshold));
    pandora::StatusCode variableIdStatusCode(pandora::XmlHelper::ReadValue(*pXmlHandle, "VariableId", m_variableId));
    pandora::StatusCode outcomeStatusCode(pandora::XmlHelper::ReadValue(*pXmlHandle, "Outcome", m_outcome));

    if (pandora::STATUS_CODE_SUCCESS == leftChildNodeIdStatusCode || pandora::STATUS_CODE_SUCCESS == rightChildNodeIdStatusCode ||
        pandora::STATUS_CODE_SUCCESS == thresholdStatusCode || pandora::STATUS_CODE_SUCCESS == variableIdStatusCode)
    {
        m_isLeaf = false; 
        m_outcome = false;
    }
    else if (outcomeStatusCode == pandora::STATUS_CODE_SUCCESS)
    {
        m_isLeaf = true; 
        m_leftChildNodeId = std::numeric_limits<int>::max();
        m_rightChildNodeId = std::numeric_limits<int>::max();
        m_threshold = std::numeric_limits<double>::max();
        m_variableId = std::numeric_limits<int>::max();
    }
    else
    {
        throw pandora::StatusCodeException(pandora::STATUS_CODE_FAILURE);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

AdaBoostDecisionTree::Node::Node(const Node &rhs) : 
    m_nodeId(rhs.m_nodeId),
    m_parentNodeId(rhs.m_parentNodeId),
    m_leftChildNodeId(rhs.m_leftChildNodeId),
    m_rightChildNodeId(rhs.m_rightChildNodeId),
    m_isLeaf(rhs.m_isLeaf),
    m_threshold(rhs.m_threshold),
    m_variableId(rhs.m_variableId),
    m_outcome(rhs.m_outcome)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

AdaBoostDecisionTree::Node &AdaBoostDecisionTree::Node::operator=(const Node &rhs)
{
    if (this != &rhs)
    {
        m_nodeId = rhs.m_nodeId;
        m_parentNodeId = rhs.m_parentNodeId;
        m_leftChildNodeId = rhs.m_leftChildNodeId;
        m_rightChildNodeId = rhs.m_rightChildNodeId;
        m_isLeaf = rhs.m_isLeaf;
        m_threshold = rhs.m_threshold;
        m_variableId = rhs.m_variableId;
        m_outcome = rhs.m_outcome;
    }

    return *this;
}

//------------------------------------------------------------------------------------------------------------------------------------------

AdaBoostDecisionTree::Node::~Node()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

AdaBoostDecisionTree::WeakClassifier::WeakClassifier(const IdToNodeMap &idToNodeMap, const double weight, const int treeId) : 
    m_idToNodeMap(idToNodeMap),
    m_weight(weight),
    m_treeId(treeId)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

AdaBoostDecisionTree::WeakClassifier::WeakClassifier(const pandora::TiXmlHandle *const pXmlHandle)
{
    for (pandora::TiXmlElement *pHeadTiXmlElement = pXmlHandle->FirstChildElement().ToElement(); pHeadTiXmlElement != NULL; pHeadTiXmlElement = pHeadTiXmlElement->NextSiblingElement())
    {
        if ("TreeIndex" == pHeadTiXmlElement->ValueStr())
        {
            PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, pandora::XmlHelper::ReadValue(*pXmlHandle, "TreeIndex", m_treeId));
        }
        else if ("TreeWeight" == pHeadTiXmlElement->ValueStr())
        {
            PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, pandora::XmlHelper::ReadValue(*pXmlHandle, "TreeWeight", m_weight));
        }
        else if ("Node" == pHeadTiXmlElement->ValueStr())
        {
            pandora::TiXmlHandle nodeHandle(pHeadTiXmlElement);
            const Node *pNode = new Node(&nodeHandle);
            m_idToNodeMap.insert(IdToNodeMap::value_type(pNode->GetNodeId(), pNode));
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

AdaBoostDecisionTree::WeakClassifier::WeakClassifier(const WeakClassifier &rhs) : 
    m_weight(rhs.m_weight),
    m_treeId(rhs.m_treeId)
{
    for (const auto &iter : rhs.m_idToNodeMap)
    {
        const Node *pNode = new Node(*iter.second);
        m_idToNodeMap.insert(IdToNodeMap::value_type(pNode->GetNodeId(), pNode));
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

AdaBoostDecisionTree::WeakClassifier &AdaBoostDecisionTree::WeakClassifier::operator=(const WeakClassifier &rhs)
{
    if (this != &rhs)
    {
        for (const auto &iter : rhs.m_idToNodeMap)
        {
            const Node *pNode = new Node(*iter.second);
            m_idToNodeMap.insert(IdToNodeMap::value_type(pNode->GetNodeId(), pNode));
        }

        m_weight = rhs.m_weight;
        m_treeId = rhs.m_treeId;
    }

    return *this;
}

//------------------------------------------------------------------------------------------------------------------------------------------

AdaBoostDecisionTree::WeakClassifier::~WeakClassifier()
{
    for (const auto &iter : m_idToNodeMap)
    {
        const Node *pNode(iter.second);
        delete pNode;
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool AdaBoostDecisionTree::WeakClassifier::Predict(const LArMvaHelper::MvaFeatureVector &features) const 
{   
    return this->EvaluateNode(0, features); 
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool AdaBoostDecisionTree::WeakClassifier::EvaluateNode(const int nodeId, const LArMvaHelper::MvaFeatureVector &features) const 
{
    const Node *pActiveNode(nullptr);

    if (m_idToNodeMap.find(nodeId) != m_idToNodeMap.end())
        pActiveNode = m_idToNodeMap.at(nodeId);
    else 
        throw pandora::StatusCodeException(pandora::STATUS_CODE_OUT_OF_RANGE);

    if (pActiveNode->IsLeaf())
        return pActiveNode->GetOutcome();

    if (features.size() <= pActiveNode->GetVariableId())
        throw pandora::StatusCodeException(pandora::STATUS_CODE_NOT_FOUND);

    if (features.at(pActiveNode->GetVariableId()) <= pActiveNode->GetThreshold())
        return this->EvaluateNode(pActiveNode->GetLeftChildNodeId(), features);
    else 
        return this->EvaluateNode(pActiveNode->GetRightChildNodeId(), features);
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

AdaBoostDecisionTree::StrongClassifier::StrongClassifier(const WeakClassifiers &weakClassifiers) : 
    m_weakClassifiers(weakClassifiers)    
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

AdaBoostDecisionTree::StrongClassifier::StrongClassifier(const pandora::TiXmlHandle *const pXmlHandle) 
{
    pandora::TiXmlElement *pCurrentXmlElement = pXmlHandle->FirstChild().Element();

    while (pCurrentXmlElement)
    {
        if (pandora::STATUS_CODE_SUCCESS != this->ReadComponent(pCurrentXmlElement))
            throw pandora::StatusCodeException(pandora::STATUS_CODE_INVALID_PARAMETER); 

        pCurrentXmlElement = pCurrentXmlElement->NextSiblingElement();
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

AdaBoostDecisionTree::StrongClassifier::StrongClassifier(const StrongClassifier &rhs)
{
    for (const WeakClassifier *pWeakClassifier : rhs.m_weakClassifiers)
    {
        const WeakClassifier *pCopyWeakClassifier = new WeakClassifier(*pWeakClassifier);
        m_weakClassifiers.push_back(pCopyWeakClassifier);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

AdaBoostDecisionTree::StrongClassifier &AdaBoostDecisionTree::StrongClassifier::operator=(const StrongClassifier &rhs)
{
    if (this != &rhs)
    {
        for (const WeakClassifier *pWeakClassifier : rhs.m_weakClassifiers)
        {
            const WeakClassifier *pCopyWeakClassifier = new WeakClassifier(*pWeakClassifier);
            m_weakClassifiers.push_back(pCopyWeakClassifier);
        }
    }

    return *this;
}

//------------------------------------------------------------------------------------------------------------------------------------------

AdaBoostDecisionTree::StrongClassifier::~StrongClassifier()
{
    for (const WeakClassifier *pWeakClassifier : m_weakClassifiers)
        delete pWeakClassifier;
}

//------------------------------------------------------------------------------------------------------------------------------------------

double AdaBoostDecisionTree::StrongClassifier::Predict(const LArMvaHelper::MvaFeatureVector &features) const
{
    double score(0.0), weights(0.0);

    for (const WeakClassifier *pWeakClassifier : m_weakClassifiers)
    {
        weights += pWeakClassifier->GetWeight();

        if (pWeakClassifier->Predict(features))
            score += pWeakClassifier->GetWeight();
        else
            score -= pWeakClassifier->GetWeight();
    }

    if (weights > std::numeric_limits<double>::min())
        score /= weights;
    else
        throw pandora::StatusCodeException(pandora::STATUS_CODE_INVALID_PARAMETER);

    return score;
}

//------------------------------------------------------------------------------------------------------------------------------------------

pandora::StatusCode AdaBoostDecisionTree::StrongClassifier::ReadComponent(pandora::TiXmlElement *pCurrentXmlElement)
{
    const std::string componentName(pCurrentXmlElement->ValueStr());
    pandora::TiXmlHandle currentHandle(pCurrentXmlElement);

    if ((std::string("Name") == componentName) || (std::string("Timestamp") == componentName))
        return pandora::STATUS_CODE_SUCCESS;

    if (std::string("DecisionTree") == componentName)
    {
        const WeakClassifier *pWeakClassifier = new WeakClassifier(&currentHandle);
        m_weakClassifiers.push_back(pWeakClassifier);
        return pandora::STATUS_CODE_SUCCESS;
    }

    return pandora::STATUS_CODE_INVALID_PARAMETER;
}

} // namespace lar_content
