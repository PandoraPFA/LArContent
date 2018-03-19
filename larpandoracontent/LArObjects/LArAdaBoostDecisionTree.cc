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

pandora::StatusCode AdaBoostDecisionTree::Initialize(const std::string &parameterLocation, const std::string &bdtName)
{
    if (!m_pStrongClassifier)
    {
        std::cout << "AdaBoostDecisionTree: bdt was already initialized" << std::endl;
        return pandora::STATUS_CODE_FAILURE;
    }

    WeakClassifiers weakClassifiers;
    const pandora::StatusCode statusCode(this->ReadXmlFile(parameterLocation, bdtName, weakClassifiers));

    if (pandora::STATUS_CODE_SUCCESS != statusCode)
        throw pandora::StatusCodeException(statusCode);

    // ATTN: Add consistency check for number of features, bearing in min not all features in a bdt may be used
    m_pStrongClassifier = new StrongClassifier(weakClassifiers);

    return pandora::STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool AdaBoostDecisionTree::Classify(const DoubleVector &features) const
{
    this->CheckInitialization();
    return (m_pStrongClassifier->Predict(features) > 0.0 ? true : false);
}

//------------------------------------------------------------------------------------------------------------------------------------------

double AdaBoostDecisionTree::CalculateClassificationScore(const DoubleVector &features) const
{
    this->CheckInitialization();
    return m_pStrongClassifier->Predict(features);
}

//------------------------------------------------------------------------------------------------------------------------------------------

double AdaBoostDecisionTree::CalculateProbability(const DoubleVector &features) const
{
    this->CheckInitialization();
    return (m_pStrongClassifier->Predict(features) + 1.0) * 0.5;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void AdaBoostDecisionTree::CheckInitialization() const
{
    if (m_pStrongClassifier == nullptr)
    {
        std::cout << "AdaBoostDecisionTree: Attempting to use an uninitialized bdt" << std::endl;
        throw pandora::StatusCodeException(pandora::STATUS_CODE_NOT_INITIALIZED);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

pandora::StatusCode AdaBoostDecisionTree::ReadXmlFile(const std::string &bdtFileName, const std::string &bdtName, WeakClassifiers &weakClassifiers) const
{
    pandora::TiXmlDocument xmlDocument(bdtFileName);

    if (!xmlDocument.LoadFile())
    {
        std::cout << "AdaBoostDecisionTree::Initialize - Invalid xml file." << std::endl;
        return pandora::STATUS_CODE_FAILURE;
    }

    const pandora::TiXmlHandle xmlDocumentHandle(&xmlDocument);
    pandora::TiXmlNode *pContainerXmlNode(pandora::TiXmlHandle(xmlDocumentHandle).FirstChildElement().Element());

    // Find the xml container for the bdt with the required name
    while (pContainerXmlNode)
    {
        if (pContainerXmlNode->ValueStr() != "AdaBoostDecisionTree")
            return pandora::STATUS_CODE_FAILURE;

        const pandora::TiXmlHandle currentHandle(pContainerXmlNode);

        std::string currentName;
        PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, pandora::XmlHelper::ReadValue(currentHandle, "Name", currentName));

        if (currentName.empty() || (currentName.size() > 1000))
        {
            std::cout << "AdaBoostDecisionTree::Initialize - Implausible bdt name extracted from xml." << std::endl;
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

    // Read the components of this bdt container
    pandora::TiXmlHandle localHandle(pContainerXmlNode);
    pandora::TiXmlElement *pCurrentXmlElement = localHandle.FirstChild().Element();

    while (pCurrentXmlElement)
    {
        if (pandora::STATUS_CODE_SUCCESS != this->ReadComponent(pCurrentXmlElement, weakClassifiers))
        {
            std::cout << "AdaBoostDecisionTree: Unknown component in xml file" << std::endl;
            return pandora::STATUS_CODE_FAILURE;
        }

        pCurrentXmlElement = pCurrentXmlElement->NextSiblingElement();
    }
    return pandora::STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

pandora::StatusCode AdaBoostDecisionTree::ReadComponent(pandora::TiXmlElement *pCurrentXmlElement, WeakClassifiers &weakClassifiers) const 
{
    const std::string componentName(pCurrentXmlElement->ValueStr());
    pandora::TiXmlHandle currentHandle(pCurrentXmlElement);

    if ((std::string("Name") == componentName) || (std::string("Timestamp") == componentName))
        return pandora::STATUS_CODE_SUCCESS;

    if (std::string("DecisionTree") == componentName)
        return this->ReadDecisionTree(currentHandle, weakClassifiers);

    return pandora::STATUS_CODE_INVALID_PARAMETER;
}

//------------------------------------------------------------------------------------------------------------------------------------------

pandora::StatusCode AdaBoostDecisionTree::ReadDecisionTree(pandora::TiXmlHandle &currentHandle, WeakClassifiers &weakClassifiers) const
{
    IdToNodeMap idToNodeMap;
    double boostWeight(0.);
    int treeIndex(0);

    for (pandora::TiXmlElement *pHeadTiXmlElement = currentHandle.FirstChildElement().ToElement(); pHeadTiXmlElement != NULL; pHeadTiXmlElement = pHeadTiXmlElement->NextSiblingElement())
    {
        if ("TreeIndex" == pHeadTiXmlElement->ValueStr())
        {
            PANDORA_RETURN_RESULT_IF_AND_IF(pandora::STATUS_CODE_SUCCESS, pandora::STATUS_CODE_NOT_FOUND, !=, pandora::XmlHelper::ReadValue(currentHandle,
                "TreeIndex", treeIndex));
        }
        else if ("BoostWeight" == pHeadTiXmlElement->ValueStr())
        {
            PANDORA_RETURN_RESULT_IF_AND_IF(pandora::STATUS_CODE_SUCCESS, pandora::STATUS_CODE_NOT_FOUND, !=, pandora::XmlHelper::ReadValue(currentHandle,
                "BoostWeight", boostWeight));
        }
        else if ("Node" == pHeadTiXmlElement->ValueStr())
        {
            pandora::TiXmlHandle nodeHandle(pHeadTiXmlElement);
            this->ReadNode(nodeHandle, idToNodeMap);
        }
    }

    weakClassifiers.emplace_back(idToNodeMap, boostWeight, treeIndex);

    return pandora::STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

pandora::StatusCode AdaBoostDecisionTree::ReadNode(pandora::TiXmlHandle &currentHandle, IdToNodeMap &idToNodeMap) const
{
    int nodeId(-1);
    PANDORA_RETURN_RESULT_IF_AND_IF(pandora::STATUS_CODE_SUCCESS, pandora::STATUS_CODE_NOT_FOUND, !=, pandora::XmlHelper::ReadValue(currentHandle,
        "NodeId", nodeId));

    int variableId(-1);
    PANDORA_RETURN_RESULT_IF_AND_IF(pandora::STATUS_CODE_SUCCESS, pandora::STATUS_CODE_NOT_FOUND, !=, pandora::XmlHelper::ReadValue(currentHandle,
        "VariableId", variableId));

    double threshold(0.);
    PANDORA_RETURN_RESULT_IF_AND_IF(pandora::STATUS_CODE_SUCCESS, pandora::STATUS_CODE_NOT_FOUND, !=, pandora::XmlHelper::ReadValue(currentHandle,
        "Threshold", threshold));

    int leftDaughterId(-1);
    PANDORA_RETURN_RESULT_IF_AND_IF(pandora::STATUS_CODE_SUCCESS, pandora::STATUS_CODE_NOT_FOUND, !=, pandora::XmlHelper::ReadValue(currentHandle,
        "LeftDaughterId", leftDaughterId));

    int rightDaughterId(-1);
    PANDORA_RETURN_RESULT_IF_AND_IF(pandora::STATUS_CODE_SUCCESS, pandora::STATUS_CODE_NOT_FOUND, !=, pandora::XmlHelper::ReadValue(currentHandle,
        "RightDaughterId", rightDaughterId));

    int parentId(-1);
    PANDORA_RETURN_RESULT_IF_AND_IF(pandora::STATUS_CODE_SUCCESS, pandora::STATUS_CODE_NOT_FOUND, !=, pandora::XmlHelper::ReadValue(currentHandle,
        "ParentId", parentId));

    bool outcome(false);
    PANDORA_RETURN_RESULT_IF_AND_IF(pandora::STATUS_CODE_SUCCESS, pandora::STATUS_CODE_NOT_FOUND, !=, pandora::XmlHelper::ReadValue(currentHandle,
        "Outcome", outcome));

    bool isLeaf(leftDaughterId == -1 ? true : false);

    const Node *pNode = new Node(nodeId, parentId, leftDaughterId, rightDaughterId, isLeaf, threshold, variableId, outcome);
    idToNodeMap.insert(IdToNodeMap::value_type(nodeId, pNode));

    return pandora::STATUS_CODE_SUCCESS;
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
//------------------------------------------------------------------------------------------------------------------------------------------

AdaBoostDecisionTree::WeakClassifier::WeakClassifier(const IdToNodeMap &idToNodeMap, const double weight, const int treeId) : 
    m_idToNodeMap(idToNodeMap),
    m_weight(weight),
    m_treeId(treeId)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool AdaBoostDecisionTree::WeakClassifier::Predict(const DoubleVector &features) const 
{   
    return this->EvaluateNode(0, features); 
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool AdaBoostDecisionTree::WeakClassifier::EvaluateNode(const int nodeId, const DoubleVector &features) const 
{
    const Node *pActiveNode(nullptr);

    if (m_idToNodeMap.find(nodeId) != m_idToNodeMap.end())
    {
        pActiveNode = m_idToNodeMap.at(nodeId);
    }

    if (pActiveNode->IsLeaf())
        return pActiveNode->GetOutcome();

    if (std::find(features.begin(), features.end(), pActiveNode->GetVariableId()) == features.end())
    {
        std::cout << "AdaBoostDecisionTree: Attempting to cut on unknown variable" << std::endl;
        throw pandora::StatusCodeException(pandora::STATUS_CODE_NOT_FOUND);
    }

    if (features.at(pActiveNode->GetVariableId()) <= pActiveNode->GetThreshold())
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

AdaBoostDecisionTree::StrongClassifier::StrongClassifier(const WeakClassifiers &weakClassifiers) : 
    m_weakClassifiers(weakClassifiers)    
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

double AdaBoostDecisionTree::StrongClassifier::Predict(const DoubleVector &features) const
{
    double score(0.0);
    double weights(0.0);

    for (const WeakClassifier &weakClassifier : m_weakClassifiers)
    {
        weights += weakClassifier.GetWeight();

        if (weakClassifier.Predict(features))
	{
            score += weakClassifier.GetWeight();
        }
        else
        {
            score -= weakClassifier.GetWeight();
        }
    }

    if (weights > std::numeric_limits<double>::min())
    {
        score /= weights;
    }
    else
    {
        std::cout << "AdaBoostDecisionTree: Classifier weights sum to zero indicating defunct classifier" << std::endl;
        throw pandora::StatusCodeException(pandora::STATUS_CODE_INVALID_PARAMETER);
    }

    return score;
}

} // namespace lar_content
