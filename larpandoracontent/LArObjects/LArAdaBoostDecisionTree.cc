/**
 *  @file   larpandoracontent/LArObjects/LArAdaBoostDecisionTree.cc
 *
 *  @brief  Implementation of the lar adaptive boost decision tree class.
 *
 *  $Log: $
 */

#include "Helpers/XmlHelper.h"

#include "larpandoracontent/LArObjects/LArAdaBoostDecisionTree.h"

using namespace pandora;

namespace lar_content
{

AdaBoostDecisionTree::AdaBoostDecisionTree() : m_pStrongClassifier(nullptr)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

AdaBoostDecisionTree::AdaBoostDecisionTree(const AdaBoostDecisionTree &rhs)
{
    m_pStrongClassifier = new StrongClassifier(*(rhs.m_pStrongClassifier));
}

//------------------------------------------------------------------------------------------------------------------------------------------

AdaBoostDecisionTree &AdaBoostDecisionTree::operator=(const AdaBoostDecisionTree &rhs)
{
    if (this != &rhs)
        m_pStrongClassifier = new StrongClassifier(*(rhs.m_pStrongClassifier));

    return *this;
}

//------------------------------------------------------------------------------------------------------------------------------------------

AdaBoostDecisionTree::~AdaBoostDecisionTree()
{
    delete m_pStrongClassifier;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode AdaBoostDecisionTree::Initialize(const std::string &bdtXmlFileName, const std::string &bdtName)
{
    if (m_pStrongClassifier)
    {
        std::cout << "AdaBoostDecisionTree: AdaBoostDecisionTree was already initialized" << std::endl;
        return STATUS_CODE_ALREADY_INITIALIZED;
    }

    TiXmlDocument xmlDocument(bdtXmlFileName);

    if (!xmlDocument.LoadFile())
    {
        std::cout << "AdaBoostDecisionTree::Initialize - Invalid xml file." << std::endl;
        return STATUS_CODE_INVALID_PARAMETER;
    }

    const TiXmlHandle xmlDocumentHandle(&xmlDocument);
    TiXmlNode *pContainerXmlNode(TiXmlHandle(xmlDocumentHandle).FirstChildElement().Element());

    while (pContainerXmlNode)
    {
        if (pContainerXmlNode->ValueStr() != "AdaBoostDecisionTree")
            return STATUS_CODE_FAILURE;

        const TiXmlHandle currentHandle(pContainerXmlNode);

        std::string currentName;
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(currentHandle, "Name", currentName));

        if (currentName.empty() || (currentName.size() > 1000))
        {
            std::cout << "AdaBoostDecisionTree::Initialize - Implausible AdaBoostDecisionTree name extracted from xml." << std::endl;
            return STATUS_CODE_INVALID_PARAMETER;
        }

        if (currentName == bdtName)
            break;

        pContainerXmlNode = pContainerXmlNode->NextSibling();
    }

    if (!pContainerXmlNode)
    {
        std::cout << "AdaBoostDecisionTree: Could not find an AdaBoostDecisionTree of name " << bdtName << std::endl;
        return STATUS_CODE_NOT_FOUND;
    }

    const TiXmlHandle xmlHandle(pContainerXmlNode);

    try
    {
        m_pStrongClassifier = new StrongClassifier(&xmlHandle);
    }
    catch (StatusCodeException &statusCodeException)
    {
        delete m_pStrongClassifier;

        if (STATUS_CODE_INVALID_PARAMETER == statusCodeException.GetStatusCode())
            std::cout << "AdaBoostDecisionTree: Initialization failure, unknown component in xml file." << std::endl;

        if (STATUS_CODE_FAILURE == statusCodeException.GetStatusCode())
            std::cout << "AdaBoostDecisionTree: Node definition does not contain expected leaf or branch variables." << std::endl;

        return statusCodeException.GetStatusCode();
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool AdaBoostDecisionTree::Classify(const LArMvaHelper::MvaFeatureVector &features) const
{
    return ((this->CalculateScore(features) > 0.) ? true : false);
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
    return (this->CalculateScore(features) + 1.) * 0.5;
}

//------------------------------------------------------------------------------------------------------------------------------------------

double AdaBoostDecisionTree::CalculateScore(const LArMvaHelper::MvaFeatureVector &features) const
{
    if (!m_pStrongClassifier)
    {
        std::cout << "AdaBoostDecisionTree: Attempting to use an uninitialized bdt" << std::endl;
        throw StatusCodeException(STATUS_CODE_NOT_INITIALIZED);
    }

    try
    {
        // TODO: Add consistency check for number of features, bearing in mind not all features in a bdt may be used
        return m_pStrongClassifier->Predict(features);
    }
    catch (StatusCodeException &statusCodeException)
    {
        if (STATUS_CODE_NOT_FOUND == statusCodeException.GetStatusCode())
        {
            std::cout << "AdaBoostDecisionTree: Caught exception thrown when trying to cut on an unknown variable." << std::endl;
        }
        else if (STATUS_CODE_INVALID_PARAMETER == statusCodeException.GetStatusCode())
        {
            std::cout << "AdaBoostDecisionTree: Caught exception thrown when classifier weights sum to zero indicating defunct classifier."
                      << std::endl;
        }
        else if (STATUS_CODE_OUT_OF_RANGE == statusCodeException.GetStatusCode())
        {
            std::cout << "AdaBoostDecisionTree: Caught exception thrown when heirarchy in decision tree is incomplete." << std::endl;
        }
        else
        {
            std::cout << "AdaBoostDecisionTree: Unexpected exception thrown." << std::endl;
        }

        throw statusCodeException;
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

AdaBoostDecisionTree::Node::Node(const TiXmlHandle *const pXmlHandle) :
    m_nodeId(0), m_parentNodeId(0), m_leftChildNodeId(0), m_rightChildNodeId(0), m_isLeaf(false), m_threshold(0.), m_variableId(0), m_outcome(false)
{
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(*pXmlHandle, "NodeId", m_nodeId));
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(*pXmlHandle, "ParentNodeId", m_parentNodeId));

    const StatusCode leftChildNodeIdStatusCode(XmlHelper::ReadValue(*pXmlHandle, "LeftChildNodeId", m_leftChildNodeId));
    const StatusCode rightChildNodeIdStatusCode(XmlHelper::ReadValue(*pXmlHandle, "RightChildNodeId", m_rightChildNodeId));
    const StatusCode thresholdStatusCode(XmlHelper::ReadValue(*pXmlHandle, "Threshold", m_threshold));
    const StatusCode variableIdStatusCode(XmlHelper::ReadValue(*pXmlHandle, "VariableId", m_variableId));
    const StatusCode outcomeStatusCode(XmlHelper::ReadValue(*pXmlHandle, "Outcome", m_outcome));

    if (STATUS_CODE_SUCCESS == leftChildNodeIdStatusCode || STATUS_CODE_SUCCESS == rightChildNodeIdStatusCode ||
        STATUS_CODE_SUCCESS == thresholdStatusCode || STATUS_CODE_SUCCESS == variableIdStatusCode)
    {
        m_isLeaf = false;
        m_outcome = false;
    }
    else if (outcomeStatusCode == STATUS_CODE_SUCCESS)
    {
        m_isLeaf = true;
        m_leftChildNodeId = std::numeric_limits<int>::max();
        m_rightChildNodeId = std::numeric_limits<int>::max();
        m_threshold = std::numeric_limits<double>::max();
        m_variableId = std::numeric_limits<int>::max();
    }
    else
    {
        throw StatusCodeException(STATUS_CODE_FAILURE);
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

AdaBoostDecisionTree::WeakClassifier::WeakClassifier(const TiXmlHandle *const pXmlHandle) : m_weight(0.), m_treeId(0)
{
    for (TiXmlElement *pHeadTiXmlElement = pXmlHandle->FirstChildElement().ToElement(); pHeadTiXmlElement != NULL;
         pHeadTiXmlElement = pHeadTiXmlElement->NextSiblingElement())
    {
        if ("TreeIndex" == pHeadTiXmlElement->ValueStr())
        {
            PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(*pXmlHandle, "TreeIndex", m_treeId));
        }
        else if ("TreeWeight" == pHeadTiXmlElement->ValueStr())
        {
            PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(*pXmlHandle, "TreeWeight", m_weight));
        }
        else if ("Node" == pHeadTiXmlElement->ValueStr())
        {
            const TiXmlHandle nodeHandle(pHeadTiXmlElement);
            const Node *pNode = new Node(&nodeHandle);
            m_idToNodeMap.insert(IdToNodeMap::value_type(pNode->GetNodeId(), pNode));
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

AdaBoostDecisionTree::WeakClassifier::WeakClassifier(const WeakClassifier &rhs) : m_weight(rhs.m_weight), m_treeId(rhs.m_treeId)
{
    for (const auto &mapEntry : rhs.m_idToNodeMap)
    {
        const Node *pNode = new Node(*(mapEntry.second));
        m_idToNodeMap.insert(IdToNodeMap::value_type(pNode->GetNodeId(), pNode));
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

AdaBoostDecisionTree::WeakClassifier &AdaBoostDecisionTree::WeakClassifier::operator=(const WeakClassifier &rhs)
{
    if (this != &rhs)
    {
        for (const auto &mapEntry : rhs.m_idToNodeMap)
        {
            const Node *pNode = new Node(*(mapEntry.second));
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
    for (const auto &mapEntry : m_idToNodeMap)
        delete mapEntry.second;
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
    {
        pActiveNode = m_idToNodeMap.at(nodeId);
    }
    else
    {
        throw StatusCodeException(STATUS_CODE_OUT_OF_RANGE);
    }

    if (pActiveNode->IsLeaf())
        return pActiveNode->GetOutcome();

    if (static_cast<int>(features.size()) <= pActiveNode->GetVariableId())
        throw StatusCodeException(STATUS_CODE_NOT_FOUND);

    if (features.at(pActiveNode->GetVariableId()).Get() <= pActiveNode->GetThreshold())
    {
        return this->EvaluateNode(pActiveNode->GetLeftChildNodeId(), features);
    }
    else
    {
        return this->EvaluateNode(pActiveNode->GetRightChildNodeId(), features);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

AdaBoostDecisionTree::StrongClassifier::StrongClassifier(const TiXmlHandle *const pXmlHandle)
{
    TiXmlElement *pCurrentXmlElement = pXmlHandle->FirstChild().Element();

    while (pCurrentXmlElement)
    {
        if (STATUS_CODE_SUCCESS != this->ReadComponent(pCurrentXmlElement))
            throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);

        pCurrentXmlElement = pCurrentXmlElement->NextSiblingElement();
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

AdaBoostDecisionTree::StrongClassifier::StrongClassifier(const StrongClassifier &rhs)
{
    for (const WeakClassifier *const pWeakClassifier : rhs.m_weakClassifiers)
        m_weakClassifiers.emplace_back(new WeakClassifier(*pWeakClassifier));
}

//------------------------------------------------------------------------------------------------------------------------------------------

AdaBoostDecisionTree::StrongClassifier &AdaBoostDecisionTree::StrongClassifier::operator=(const StrongClassifier &rhs)
{
    if (this != &rhs)
    {
        for (const WeakClassifier *const pWeakClassifier : rhs.m_weakClassifiers)
            m_weakClassifiers.emplace_back(new WeakClassifier(*pWeakClassifier));
    }

    return *this;
}

//------------------------------------------------------------------------------------------------------------------------------------------

AdaBoostDecisionTree::StrongClassifier::~StrongClassifier()
{
    for (const WeakClassifier *const pWeakClassifier : m_weakClassifiers)
        delete pWeakClassifier;
}

//------------------------------------------------------------------------------------------------------------------------------------------

double AdaBoostDecisionTree::StrongClassifier::Predict(const LArMvaHelper::MvaFeatureVector &features) const
{
    double score(0.), weights(0.);

    for (const WeakClassifier *const pWeakClassifier : m_weakClassifiers)
    {
        weights += pWeakClassifier->GetWeight();

        if (pWeakClassifier->Predict(features))
        {
            score += pWeakClassifier->GetWeight();
        }
        else
        {
            score -= pWeakClassifier->GetWeight();
        }
    }

    if (weights > std::numeric_limits<double>::epsilon())
    {
        score /= weights;
    }
    else
    {
        throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);
    }

    return score;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode AdaBoostDecisionTree::StrongClassifier::ReadComponent(TiXmlElement *pCurrentXmlElement)
{
    const std::string componentName(pCurrentXmlElement->ValueStr());
    TiXmlHandle currentHandle(pCurrentXmlElement);

    if ((std::string("Name") == componentName) || (std::string("Timestamp") == componentName))
        return STATUS_CODE_SUCCESS;

    if (std::string("DecisionTree") == componentName)
    {
        m_weakClassifiers.emplace_back(new WeakClassifier(&currentHandle));
        return STATUS_CODE_SUCCESS;
    }

    return STATUS_CODE_INVALID_PARAMETER;
}

} // namespace lar_content
