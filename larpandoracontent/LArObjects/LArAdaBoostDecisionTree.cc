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
    IDToNodeMap idToNodeMap;
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

pandora::StatusCode AdaBoostDecisionTree::ReadNode(pandora::TiXmlHandle &currentHandle, IDToNodeMap &idToNodeMap) const
{
    int nodeID(-1);
    PANDORA_RETURN_RESULT_IF_AND_IF(pandora::STATUS_CODE_SUCCESS, pandora::STATUS_CODE_NOT_FOUND, !=, pandora::XmlHelper::ReadValue(currentHandle,
        "NodeID", nodeID));

    int variableID(-1);
    PANDORA_RETURN_RESULT_IF_AND_IF(pandora::STATUS_CODE_SUCCESS, pandora::STATUS_CODE_NOT_FOUND, !=, pandora::XmlHelper::ReadValue(currentHandle,
        "VariableID", variableID));

    double threshold(0.);
    PANDORA_RETURN_RESULT_IF_AND_IF(pandora::STATUS_CODE_SUCCESS, pandora::STATUS_CODE_NOT_FOUND, !=, pandora::XmlHelper::ReadValue(currentHandle,
        "Threshold", threshold));

    int leftDaughterID(-1);
    PANDORA_RETURN_RESULT_IF_AND_IF(pandora::STATUS_CODE_SUCCESS, pandora::STATUS_CODE_NOT_FOUND, !=, pandora::XmlHelper::ReadValue(currentHandle,
        "LeftDaughterID", leftDaughterID));

    int rightDaughterID(-1);
    PANDORA_RETURN_RESULT_IF_AND_IF(pandora::STATUS_CODE_SUCCESS, pandora::STATUS_CODE_NOT_FOUND, !=, pandora::XmlHelper::ReadValue(currentHandle,
        "RightDaughterID", rightDaughterID));

    int parentID(-1);
    PANDORA_RETURN_RESULT_IF_AND_IF(pandora::STATUS_CODE_SUCCESS, pandora::STATUS_CODE_NOT_FOUND, !=, pandora::XmlHelper::ReadValue(currentHandle,
        "ParentID", parentID));

    bool outcome(false);
    PANDORA_RETURN_RESULT_IF_AND_IF(pandora::STATUS_CODE_SUCCESS, pandora::STATUS_CODE_NOT_FOUND, !=, pandora::XmlHelper::ReadValue(currentHandle,
        "Outcome", outcome));

    bool isLeaf(leftDaughterID == -1 ? true : false);

    const Node *pNode = new Node(nodeID, parentID, leftDaughterID, rightDaughterID, isLeaf, threshold, variableID, outcome);
    idToNodeMap.insert(IDToNodeMap::value_type(nodeID, pNode));

    return pandora::STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

AdaBoostDecisionTree::Node::Node(const int &nodeID, const int &parentNodeID, const int &leftChildNodeID, const int &rightChildNodeID, const bool &isLeaf, const double &threshold, const int &variableID, const bool outcome) : 
    m_nodeID(nodeID),
    m_parentNodeID(parentNodeID),
    m_leftChildNodeID(leftChildNodeID),
    m_rightChildNodeID(rightChildNodeID),
    m_isLeaf(isLeaf),
    m_threshold(threshold),
    m_variableID(variableID),
    m_outcome(outcome)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

AdaBoostDecisionTree::WeakClassifier::WeakClassifier(IDToNodeMap &idToNodeMap, const double &weight, const int &treeID) : 
    m_idToNodeMap(idToNodeMap),
    m_weight(weight),
    m_treeID(treeID)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool AdaBoostDecisionTree::WeakClassifier::Predict(const DoubleVector &features) const 
{   
    return this->EvaluateNode(0, features); 
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool AdaBoostDecisionTree::WeakClassifier::EvaluateNode(const int &nodeID, const DoubleVector &features) const 
{
    const Node *pActiveNode(nullptr);

    if (m_idToNodeMap.find(nodeID) != m_idToNodeMap.end())
    {
        pActiveNode = m_idToNodeMap.at(nodeID);
    }

    if (pActiveNode->IsLeaf())
        return pActiveNode->GetOutcome();

    if (features.at(pActiveNode->GetVariableID()) <= pActiveNode->GetThreshold())
    {
        return this->EvaluateNode(pActiveNode->GetLeftChildNodeID(), features);
    }
    else 
    {
        return this->EvaluateNode(pActiveNode->GetRightChildNodeID(), features);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

AdaBoostDecisionTree::StrongClassifier::StrongClassifier(WeakClassifiers &weakClassifiers) : 
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

    score /= weights;

    return score;
}

} // namespace lar_content
