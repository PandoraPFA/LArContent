/**
 *  @file   larpandoracontent/LArObjects/LArSupportVectorMachine.cc
 *
 *  @brief  Implementation of the lar support vector machine class.
 *
 *  $Log: $
 */

#include "Helpers/XmlHelper.h"

#include "larpandoracontent/LArObjects/LArSupportVectorMachine.h"

using namespace pandora;

namespace lar_content
{

SupportVectorMachine::SupportVectorMachine() :
    m_isInitialized(false),
    m_enableProbability(false),
    m_probAParameter(0.),
    m_probBParameter(0.),
    m_standardizeFeatures(true),
    m_nFeatures(0),
    m_bias(0.),
    m_scaleFactor(1.),
    m_kernelType(QUADRATIC),
    m_kernelFunction(QuadraticKernel),
    m_kernelMap{{LINEAR, LinearKernel}, {QUADRATIC, QuadraticKernel}, {CUBIC, CubicKernel}, {GAUSSIAN_RBF, GaussianRbfKernel}}
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode SupportVectorMachine::Initialize(const std::string &parameterLocation, const std::string &svmName)
{
    if (m_isInitialized)
    {
        std::cout << "SupportVectorMachine: svm was already initialized" << std::endl;
        return STATUS_CODE_ALREADY_INITIALIZED;
    }

    this->ReadXmlFile(parameterLocation, svmName);

    // Check the sizes of sigma and scale factor if they are to be used as divisors
    if (m_standardizeFeatures)
    {
        for (const FeatureInfo &featureInfo : m_featureInfoList)
        {
            if (featureInfo.m_sigmaValue < std::numeric_limits<double>::epsilon())
            {
                std::cout << "SupportVectorMachine: could not standardize parameters because sigma value was too small" << std::endl;
                throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);
            }
        }
    }

    // Check the number of features is consistent.
    m_nFeatures = m_featureInfoList.size();

    for (const SupportVectorInfo &svInfo : m_svInfoList)
    {
        if (svInfo.m_supportVector.size() != m_nFeatures)
        {
            std::cout << "SupportVectorMachine: the number of features in the xml file was inconsistent" << std::endl;
            throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);
        }
    }

    // There's the possibility of a user-defined kernel that doesn't use this as a divisor but let's be safe
    if (m_scaleFactor < std::numeric_limits<double>::epsilon())
    {
        std::cout << "SupportVectorMachine: could not evaluate kernel because scale factor was too small" << std::endl;
        throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);
    }

    m_isInitialized = true;
    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void SupportVectorMachine::ReadXmlFile(const std::string &svmFileName, const std::string &svmName)
{
    TiXmlDocument xmlDocument(svmFileName);

    if (!xmlDocument.LoadFile())
    {
        std::cout << "SupportVectorMachine::Initialize - Invalid xml file." << std::endl;
        throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);
    }

    const TiXmlHandle xmlDocumentHandle(&xmlDocument);
    TiXmlNode *pContainerXmlNode(TiXmlHandle(xmlDocumentHandle).FirstChildElement().Element());

    // Try to find the svm container with the required name
    while (pContainerXmlNode)
    {
        if (pContainerXmlNode->ValueStr() != "SupportVectorMachine")
            throw StatusCodeException(STATUS_CODE_FAILURE);

        const TiXmlHandle currentHandle(pContainerXmlNode);

        std::string currentName;
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(currentHandle, "Name", currentName));

        if (currentName.empty() || (currentName.size() > 1000))
        {
            std::cout << "SupportVectorMachine::Initialize - Implausible svm name extracted from xml." << std::endl;
            throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);
        }

        if (currentName == svmName)
            break;

        pContainerXmlNode = pContainerXmlNode->NextSibling();
    }

    if (!pContainerXmlNode)
    {
        std::cout << "SupportVectorMachine: Could not find an svm by the name " << svmName << std::endl;
        throw StatusCodeException(STATUS_CODE_NOT_FOUND);
    }

    // Read the components of this svm container
    TiXmlHandle localHandle(pContainerXmlNode);
    TiXmlElement *pCurrentXmlElement = localHandle.FirstChild().Element();

    while (pCurrentXmlElement)
    {
        if (STATUS_CODE_SUCCESS != this->ReadComponent(pCurrentXmlElement))
        {
            std::cout << "SupportVectorMachine: Unknown component in xml file" << std::endl;
            throw StatusCodeException(STATUS_CODE_FAILURE);
        }

        pCurrentXmlElement = pCurrentXmlElement->NextSiblingElement();
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode SupportVectorMachine::ReadComponent(TiXmlElement *pCurrentXmlElement)
{
    const std::string componentName(pCurrentXmlElement->ValueStr());
    const TiXmlHandle currentHandle(pCurrentXmlElement);

    if ((std::string("Name") == componentName) || (std::string("Timestamp") == componentName))
        return STATUS_CODE_SUCCESS;

    if (std::string("Machine") == componentName)
        return this->ReadMachine(currentHandle);

    if (std::string("Features") == componentName)
        return this->ReadFeatures(currentHandle);

    if (std::string("SupportVector") == componentName)
        return this->ReadSupportVector(currentHandle);

    return STATUS_CODE_INVALID_PARAMETER;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode SupportVectorMachine::ReadMachine(const TiXmlHandle &currentHandle)
{
    int kernelType(0);
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(currentHandle, "KernelType", kernelType));

    double bias(0.);
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(currentHandle, "Bias", bias));

    double scaleFactor(0.);
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(currentHandle, "ScaleFactor", scaleFactor));

    bool standardize(true);
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(currentHandle, "Standardize", standardize));

    bool enableProbability(false);
    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(currentHandle, "EnableProbability", enableProbability));

    double probAParameter(0.);
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(currentHandle, "ProbAParameter", probAParameter));

    double probBParameter(0.);
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(currentHandle, "ProbBParameter", probBParameter));

    m_kernelType = static_cast<KernelType>(kernelType);
    m_bias = bias;
    m_scaleFactor = scaleFactor;
    m_enableProbability = enableProbability;
    m_probAParameter = probAParameter;
    m_probBParameter = probBParameter;

    if (kernelType != USER_DEFINED) // if user-defined, leave it so it alone can be set before/after initialization
        m_kernelFunction = m_kernelMap.at(m_kernelType);

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode SupportVectorMachine::ReadFeatures(const TiXmlHandle &currentHandle)
{
    std::vector<double> muValues;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadVectorOfValues(currentHandle, "MuValues", muValues));

    std::vector<double> sigmaValues;
    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadVectorOfValues(currentHandle, "SigmaValues", sigmaValues));

    if (muValues.size() != sigmaValues.size())
    {
        std::cout << "SupportVectorMachine: could not add feature info because the size of mu (" << muValues.size()
                  << ") did not match "
                     "the size of sigma ("
                  << sigmaValues.size() << ")" << std::endl;
        return STATUS_CODE_INVALID_PARAMETER;
    }

    m_featureInfoList.reserve(muValues.size());

    for (std::size_t i = 0; i < muValues.size(); ++i)
        m_featureInfoList.emplace_back(muValues.at(i), sigmaValues.at(i));

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode SupportVectorMachine::ReadSupportVector(const TiXmlHandle &currentHandle)
{
    double yAlpha(0.0);
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(currentHandle, "AlphaY", yAlpha));

    std::vector<double> values;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadVectorOfValues(currentHandle, "Values", values));

    LArMvaHelper::MvaFeatureVector valuesFeatureVector;
    for (const double &value : values)
        valuesFeatureVector.emplace_back(value);

    m_svInfoList.emplace_back(yAlpha, valuesFeatureVector);
    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

double SupportVectorMachine::CalculateClassificationScoreImpl(const LArMvaHelper::MvaFeatureVector &features) const
{
    if (!m_isInitialized)
    {
        std::cout << "SupportVectorMachine: could not perform classification because the svm was uninitialized" << std::endl;
        throw StatusCodeException(STATUS_CODE_NOT_INITIALIZED);
    }

    if (m_svInfoList.empty())
    {
        std::cout << "SupportVectorMachine: could not perform classification because the initialized svm had no support vectors in the model"
                  << std::endl;
        throw StatusCodeException(STATUS_CODE_NOT_INITIALIZED);
    }

    LArMvaHelper::MvaFeatureVector standardizedFeatures;
    standardizedFeatures.reserve(m_nFeatures);

    if (m_standardizeFeatures)
    {
        for (std::size_t i = 0; i < m_nFeatures; ++i)
            standardizedFeatures.push_back(m_featureInfoList.at(i).StandardizeParameter(features.at(i).Get()));
    }

    double classScore(0.);
    for (const SupportVectorInfo &supportVectorInfo : m_svInfoList)
    {
        classScore += supportVectorInfo.m_yAlpha *
            m_kernelFunction(supportVectorInfo.m_supportVector, (m_standardizeFeatures ? standardizedFeatures : features), m_scaleFactor);
    }

    return classScore + m_bias;
}

} // namespace lar_content
