/**
 *  @file   larpandoracontent/LArObjects/LArSupportVectorMachine.h
 * 
 *  @brief  Header file for the lar support vector machine class.
 * 
 *  $Log: $
 */
#ifndef LAR_SUPPORT_VECTOR_MACHINE_H
#define LAR_SUPPORT_VECTOR_MACHINE_H 1

#include "Pandora/AlgorithmTool.h"

#include <functional>
#include <typeindex>

namespace lar_content
{

template<typename T, typename TALG, typename ...Ts>
class SVMFeatureTool;
    
/**
 *  @brief  SVMFeatureToolBase base class template
 */
template<typename TALG, typename ...Ts>
class SVMFeatureToolBase : public pandora::AlgorithmTool
{
public:
    /**
     *  @brief  Default virtual destructor
     */
    virtual ~SVMFeatureToolBase() = default;
    
    using FeatureToolMap = std::multimap<std::type_index, SVMFeatureToolBase<TALG, Ts...> *>; ///< Alias template for a map between feature types and feature tools
    
    template <typename T>
    using FeatureTool = SVMFeatureTool<T, TALG, Ts...>;
};

//------------------------------------------------------------------------------------------------------------------------------------------

template <typename TALG, typename ...Ts>
using SVMFeatureToolMap = std::multimap<std::type_index, SVMFeatureToolBase<TALG, Ts...> *>; ///< Alias template for a map between feature types and feature tools

//------------------------------------------------------------------------------------------------------------------------------------------

/**
 *  @brief  SVMFeatureTool class template
 */
template<typename T, typename TALG, typename ...Ts>
class SVMFeatureTool : public SVMFeatureToolBase<TALG, Ts...>
{
public: 
    typedef typename std::decay<T>::type ReturnType;
    typedef std::vector<ReturnType> FeatureList;
    
    SVMFeatureTool();
    
    /**
     *  @brief  Run the algorithm tool
     *
     *  @param  pAlg address of the calling algorithm 
     *  @param  args arguments to pass to the tool
     *
     *  @return the feature
     */
    virtual ReturnType Run(const TALG *const pAlg, Ts ... args) = 0;
};

//------------------------------------------------------------------------------------------------------------------------------------------
    
/**
 *  @brief  SupportVectorMachine base class
 */
class SupportVectorMachineBase
{
public:
    typedef std::vector<double> DoubleVector;
    typedef std::vector<float>  FloatVector;
    typedef std::vector<int>    IntVector;
    typedef std::vector<bool>   BoolVector;

    /**
     *  @brief  KernelType enum
     */
    enum KernelType
    {
        USER_DEFINED = 0,
        LINEAR       = 1,
        QUADRATIC    = 2,
        CUBIC        = 3,
        GAUSSIAN_RBF = 4
    };

    //--------------------------------------------------------------------------------------------------------------------------------------

    /**
     *  @brief  Default virtual destructor
     */
    virtual ~SupportVectorMachineBase() = default;
};

//------------------------------------------------------------------------------------------------------------------------------------------

/**
 *  @brief  SupportVectorMachine class template
 */
template <std::size_t NFEATURES>
class SupportVectorMachine : SupportVectorMachineBase
{
public:
    typedef std::array<double, NFEATURES> FeatureArray;    
    typedef std::function<double(const FeatureArray &, const FeatureArray &, const double)> KernelFunction;
   
    /**
     *  @brief  Default constructor.
     */
    SupportVectorMachine();
    
    /**
     *  @brief  Initialize the SVM using a serialized model.
     * 
     *  @param  parameterLocation the location of the model
     *  @param  svmName the name of the model
     * 
     *  @return success
     */
    pandora::StatusCode Initialize(const std::string &parameterLocation, const std::string &svmName);
    
    //--------------------------------------------------------------------------------------------------------------------------------------
    
    /**
     *  @brief  Make a classification for a set of input features, based on the trained model
     *  
     *  @param  features the input features
     * 
     *  @return the predicted boolean class
     */
    bool Classify(const FeatureArray &features) const;
    
    /**
     *  @brief  Calculate the classification score for a set of input features, based on the trained model 
     *  
     *  @param  features the input features
     * 
     *  @return the classification score
     */
    double CalculateClassificationScore(const FeatureArray &features) const;
    
    /**
     *  @brief  Get whether this SVM is initialized
     *  
     *  @return whether the model is initialized
     */
    bool IsInitialized() const;
    
    /**
     *  @brief  Get the number of features
     *  
     *  @return the number of features
     */
    std::size_t GetNumFeatures() const;
    
    /**
     *  @brief  Set the kernel function to use
     *  
     *  @param kernelFunction the kernel function
     */
    void SetKernelFunction(KernelFunction kernelFunction);

private:
    /**
     *  @brief  SupportVectorInfo class
     */
    class SupportVectorInfo
    {
    public:
        /**
         *  @brief  Constructor
         *  
         *  @param  yAlpha the alpha value multiplied by the y-value for the support vector
         *  @param  supportVector the support vector
         */
        SupportVectorInfo(const double yAlpha, FeatureArray supportVector);

        double       m_yAlpha;        ///< The alpha-value multiplied by the y-value for the support vector
        FeatureArray m_supportVector; ///< The support vector
    };
    
    //--------------------------------------------------------------------------------------------------------------------------------------
    
    /**
     *  @brief  FeatureInfo class.
     */
    class FeatureInfo
    {
    public:
        /**
         *  @brief  FeatureInfo class.
         * 
         *  @param  muValue the average value of this feature
         *  @param  sigmaValue the stddev of this feature
         */
        FeatureInfo(const double muValue, const double sigmaValue);
        
        /**
         *  @brief  Default constructor to allow default-construction of (uninitialized) SVMs.
         */
        FeatureInfo() : m_muValue(0.0), m_sigmaValue(0.0) {}
    
        /**
         *  @brief  Standardize a parameter corresponding to this feature
         * 
         *  @param  parameter the parameter
         * 
         *  @return the standardized parameter
         */
        double StandardizeParameter(const double parameter) const;
    
        double    m_muValue;       ///< The average value of this feature
        double    m_sigmaValue;    ///< The stddev of this feature
    };
    
    //--------------------------------------------------------------------------------------------------------------------------------------

    typedef std::vector<SupportVectorInfo>     SVInfoList;
    typedef std::array<FeatureInfo, NFEATURES> FeatureInfoArray;
    
    typedef std::map<KernelType, KernelFunction> KernelMap;
    
    bool              m_isInitialized;       ///< Whether this SVM has been initialized
    
    bool              m_standardizeFeatures; ///< Whether to standardize the features
    double            m_bias;                ///< The bias term
    double            m_scaleFactor;         ///< The kernel scale factor
    
    SVInfoList        m_svInfoList;          ///< The list of SupportVectorInfo objects
    FeatureInfoArray  m_featureInfoList;     ///< The list of FeatureInfo objects
    
    KernelType        m_kernelType;          ///< The kernel type
    KernelFunction    m_kernelFunction;      ///< The kernel function
    
    const KernelMap   m_kernelMap{{LINEAR,       LinearKernel},
                                  {QUADRATIC,    QuadraticKernel},
                                  {CUBIC,        CubicKernel},
                                  {GAUSSIAN_RBF, GaussianRbfKernel}};    ///< Map from the kernel types to the kernel functions
    
    //--------------------------------------------------------------------------------------------------------------------------------------
    
    /**
     *  @brief  Read the SVM parameters from an xml file
     * 
     *  @param  svmFileName the sml file name
     *  @param  svmName the name of the SVM
     */
    void ReadXmlFile(const std::string &svmFileName, const std::string &svmName);
    
    /**
     *  @brief  Read the component at the current xml element
     * 
     *  @param  pCurrentXmlElement address of the current xml element
     * 
     *  @return success
     */
    pandora::StatusCode ReadComponent(pandora::TiXmlElement *pCurrentXmlElement);
    
    /**
     *  @brief  Read the machine component at the current xml handle
     * 
     *  @param  currentHandle the current xml handle
     * 
     *  @return success
     */
    pandora::StatusCode ReadMachine(const pandora::TiXmlHandle &currentHandle);
    
    /**
     *  @brief  Read the feature component at the current xml handle
     * 
     *  @param  currentHandle the current xml handle
     * 
     *  @return success
     */
    pandora::StatusCode ReadFeatures(const pandora::TiXmlHandle &currentHandle);
    
    /**
     *  @brief  Read the support vector component at the current xml handle
     * 
     *  @param  currentHandle the current xml handle
     * 
     *  @return success
     */
    pandora::StatusCode ReadSupportVector(const pandora::TiXmlHandle &currentHandle);
    
    /**
     *  @brief  Set the machine parameters
     *  
     *  @param  kernelType the type of kernel to use
     *  @param  bias the bias term
     *  @param  scaleFactor the optional kernel scale factor
     * 
     *  @return success
     */
    pandora::StatusCode SetMachineParameters(KernelType kernelType, const double bias, const double scaleFactor);
    
    /**
     *  @brief  Set the feature parameters
     *  
     *  @param  muValues the average feature values
     *  @param  sigmaValues the stdev of the feature values
     * 
     *  @return success
     */
    pandora::StatusCode SetFeatureParameters(const DoubleVector &muValues, const DoubleVector &sigmaValues);
    
    /**
     *  @brief  Add a support vector
     *  
     *  @param  yAlpha the alpha value multiplied by the y value for this vector
     *  @param  values the components of the vector
     * 
     *  @return success
     */
    pandora::StatusCode AddSupportVector(const double yAlpha, const DoubleVector &values);
    
    /**
     *  @brief  Implementation method for calculating the classification score using the trained model.
     * 
     *  @param  features the array of features
     * 
     *  @return the classification score
     */
    double CalculateClassificationScoreImpl(const FeatureArray &features) const;
    
    //--------------------------------------------------------------------------------------------------------------------------------------
    
    /**
     *  @brief  An inhomogeneous quadratic kernel
     *  
     *  @param  supportVector the support vector
     *  @param  features the features
     *  @param  scale factor the scale factor
     * 
     *  @return result of the kernel operation
     */
    static double QuadraticKernel(const FeatureArray &supportVector, const FeatureArray &features, const double scaleFactor = 1.f);
    
    /**
     *  @brief  An inhomogeneous cubic kernel
     *  
     *  @param  supportVector the support vector
     *  @param  features the features
     *  @param  scale factor the scale factor
     * 
     *  @return result of the kernel operation
     */
    static double CubicKernel(const FeatureArray &supportVector, const FeatureArray &features, const double scaleFactor = 1.f);
    
    /**
     *  @brief  A linear kernel
     *  
     *  @param  supportVector the support vector
     *  @param  features the features
     *  @param  scale factor the scale factor
     * 
     *  @return result of the kernel operation
     */
    static double LinearKernel(const FeatureArray &supportVector, const FeatureArray &features, const double scaleFactor = 1.f);
    
    /**
     *  @brief  A gaussian RBF kernel
     *  
     *  @param  supportVector the support vector
     *  @param  features the features
     *  @param  scale factor the scale factor
     * 
     *  @return result of the kernel operation
     */
    static double GaussianRbfKernel(const FeatureArray &supportVector, const FeatureArray &features, const double scaleFactor = 1.f);
};

//------------------------------------------------------------------------------------------------------------------------------------------

template <std::size_t NFEATURES>
SupportVectorMachine<NFEATURES>::SupportVectorMachine() :
    m_isInitialized(false),
    m_standardizeFeatures(true),
    m_bias(0.0),
    m_scaleFactor(1.0),
    m_kernelType(QUADRATIC),
    m_kernelFunction(QuadraticKernel)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

template <std::size_t NFEATURES>
pandora::StatusCode SupportVectorMachine<NFEATURES>::Initialize(const std::string &parameterLocation, const std::string &svmName)
{
    if (m_isInitialized)
    {
        std::cout << "SupportVectorMachine: SVM was already initialized" << std::endl;
        return pandora::STATUS_CODE_FAILURE;
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
                throw pandora::StatusCodeException(pandora::STATUS_CODE_INVALID_PARAMETER); 
            }
        }
    }
    
    // There's the possibility of a user-defined kernel that doesn't use this as a divisor but let's be safe
    if (m_scaleFactor < std::numeric_limits<double>::epsilon()) 
    {
        std::cout << "SupportVectorMachine: could not evaluate kernel because scale factor was too small" << std::endl;
        throw pandora::StatusCodeException(pandora::STATUS_CODE_INVALID_PARAMETER); 
    }
        
    m_isInitialized = true;
    return pandora::STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

template <std::size_t NFEATURES>
inline bool SupportVectorMachine<NFEATURES>::Classify(const FeatureArray &features) const
{
    return (this->CalculateClassificationScoreImpl(features) > 0.0);
}
    
//------------------------------------------------------------------------------------------------------------------------------------------
    
template <std::size_t NFEATURES>
inline double SupportVectorMachine<NFEATURES>::CalculateClassificationScore(const FeatureArray &features) const
{
    return this->CalculateClassificationScoreImpl(features);
}

//------------------------------------------------------------------------------------------------------------------------------------------
    
template <std::size_t NFEATURES>
inline bool SupportVectorMachine<NFEATURES>::IsInitialized() const
{
    return m_isInitialized;
}

//------------------------------------------------------------------------------------------------------------------------------------------
    
template <std::size_t NFEATURES>
inline std::size_t SupportVectorMachine<NFEATURES>::GetNumFeatures() const
{
    return NFEATURES;
}

//------------------------------------------------------------------------------------------------------------------------------------------

template <std::size_t NFEATURES>
inline void SupportVectorMachine<NFEATURES>::SetKernelFunction(KernelFunction kernelFunction)
{
    m_kernelFunction = std::move(kernelFunction);
}
    
//------------------------------------------------------------------------------------------------------------------------------------------

template <std::size_t NFEATURES>
void SupportVectorMachine<NFEATURES>::ReadXmlFile(const std::string &svmFileName, const std::string &svmName)
{
    pandora::TiXmlDocument xmlDocument(svmFileName);

    if (!xmlDocument.LoadFile())
    {
        std::cout << "SupportVectorMachine::Initialize - Invalid xml file." << std::endl;
        throw pandora::StatusCodeException(pandora::STATUS_CODE_FAILURE);
    }

    const pandora::TiXmlHandle xmlDocumentHandle(&xmlDocument);
    pandora::TiXmlNode *pContainerXmlNode(pandora::TiXmlHandle(xmlDocumentHandle).FirstChildElement().Element());
    
    // Try to find the SVM container with the required name
    while (pContainerXmlNode)
    {
        if (pContainerXmlNode->ValueStr() != "SupportVectorMachine")
            throw pandora::StatusCodeException(pandora::STATUS_CODE_FAILURE);
            
        const pandora::TiXmlHandle currentHandle(pContainerXmlNode);
        
        std::string currentName;
        PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, pandora::XmlHelper::ReadValue(currentHandle, "Name", currentName));
            
        if (currentName == svmName)
            break;
        
        pContainerXmlNode = pContainerXmlNode->NextSibling();
    }
    
    if (!pContainerXmlNode)
    {
        std::cout << "SupportVectorMachine: Could not find an SVM by the name " << svmName << std::endl;
        throw pandora::StatusCodeException(pandora::STATUS_CODE_NOT_FOUND);
    }
    
    // Read the components of this SVM container
    pandora::TiXmlHandle localHandle(pContainerXmlNode);
    pandora::TiXmlElement *pCurrentXmlElement = localHandle.FirstChild().Element();
    
    while (pCurrentXmlElement)
    {
        if (pandora::STATUS_CODE_SUCCESS != this->ReadComponent(pCurrentXmlElement))
        {
            std::cout << "SupportVectorMachine: Unknown component in xml file" << std::endl;
            throw pandora::StatusCodeException(pandora::STATUS_CODE_FAILURE);
        }
        
        pCurrentXmlElement = pCurrentXmlElement->NextSiblingElement();
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

template <std::size_t NFEATURES>
pandora::StatusCode SupportVectorMachine<NFEATURES>::ReadComponent(pandora::TiXmlElement *pCurrentXmlElement)
{
    const std::string componentName(pCurrentXmlElement->ValueStr());
    const pandora::TiXmlHandle currentHandle(pCurrentXmlElement);

    if ((std::string("Name") == componentName) || (std::string("Timestamp") == componentName))
        return pandora::STATUS_CODE_SUCCESS;
    
    if (std::string("Machine") == componentName)
        return this->ReadMachine(currentHandle);
    
    if (std::string("Features") == componentName)
        return this->ReadFeatures(currentHandle);
        
    if (std::string("SupportVector") == componentName)
        return this->ReadSupportVector(currentHandle);
        
    return pandora::STATUS_CODE_INVALID_PARAMETER;
}

//------------------------------------------------------------------------------------------------------------------------------------------

template <std::size_t NFEATURES>
pandora::StatusCode SupportVectorMachine<NFEATURES>::ReadMachine(const pandora::TiXmlHandle &currentHandle)
{
    int kernelType(0);
    PANDORA_RETURN_RESULT_IF_AND_IF(pandora::STATUS_CODE_SUCCESS, pandora::STATUS_CODE_NOT_FOUND, !=, pandora::XmlHelper::ReadValue(currentHandle,
            "KernelType", kernelType));
    
    double bias(0.0);
    PANDORA_RETURN_RESULT_IF_AND_IF(pandora::STATUS_CODE_SUCCESS, pandora::STATUS_CODE_NOT_FOUND, !=, pandora::XmlHelper::ReadValue(currentHandle,
            "Bias", bias));
            
    double scaleFactor(0.0);    
    PANDORA_RETURN_RESULT_IF_AND_IF(pandora::STATUS_CODE_SUCCESS, pandora::STATUS_CODE_NOT_FOUND, !=, pandora::XmlHelper::ReadValue(currentHandle,
            "ScaleFactor", scaleFactor));
            
    bool standardize(true);
    PANDORA_RETURN_RESULT_IF_AND_IF(pandora::STATUS_CODE_SUCCESS, pandora::STATUS_CODE_NOT_FOUND, !=, pandora::XmlHelper::ReadValue(currentHandle,
            "Standardize", standardize));
    
    return this->SetMachineParameters(static_cast<SupportVectorMachineBase::KernelType>(kernelType), bias, scaleFactor);
}    

//------------------------------------------------------------------------------------------------------------------------------------------

template <std::size_t NFEATURES>
pandora::StatusCode SupportVectorMachine<NFEATURES>::ReadFeatures(const pandora::TiXmlHandle &currentHandle)
{
    DoubleVector muValues;    
    PANDORA_RETURN_RESULT_IF_AND_IF(pandora::STATUS_CODE_SUCCESS, pandora::STATUS_CODE_NOT_FOUND, !=, pandora::XmlHelper::ReadVectorOfValues(currentHandle,
            "MuValues", muValues));
    
    DoubleVector sigmaValues;
    PANDORA_RETURN_RESULT_IF_AND_IF(pandora::STATUS_CODE_SUCCESS, pandora::STATUS_CODE_NOT_FOUND, !=, pandora::XmlHelper::ReadVectorOfValues(currentHandle,
            "SigmaValues", sigmaValues));
    
    return this->SetFeatureParameters(muValues, sigmaValues);
}

//------------------------------------------------------------------------------------------------------------------------------------------

template <std::size_t NFEATURES>
pandora::StatusCode SupportVectorMachine<NFEATURES>::ReadSupportVector(const pandora::TiXmlHandle &currentHandle)
{
    double yAlpha(0.0);
    PANDORA_RETURN_RESULT_IF_AND_IF(pandora::STATUS_CODE_SUCCESS, pandora::STATUS_CODE_NOT_FOUND, !=, pandora::XmlHelper::ReadValue(currentHandle,
            "AlphaY", yAlpha));
    
    DoubleVector values;
    PANDORA_RETURN_RESULT_IF_AND_IF(pandora::STATUS_CODE_SUCCESS, pandora::STATUS_CODE_NOT_FOUND, !=, pandora::XmlHelper::ReadVectorOfValues(currentHandle,
            "Values", values));
    
    return this->AddSupportVector(yAlpha, values);
}

//------------------------------------------------------------------------------------------------------------------------------------------

template <std::size_t NFEATURES>
pandora::StatusCode SupportVectorMachine<NFEATURES>::SetMachineParameters(const KernelType kernelType, const double bias, const double scaleFactor)
{
    m_kernelType = kernelType;
    m_bias = bias;
    m_scaleFactor = scaleFactor;
    
    if (kernelType != USER_DEFINED) // if user-defined, leave it so it alone can be set before/after initialization
        m_kernelFunction = m_kernelMap.at(m_kernelType);
        
    return pandora::STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

template <std::size_t NFEATURES>
pandora::StatusCode SupportVectorMachine<NFEATURES>::SetFeatureParameters(const DoubleVector &muValues, const DoubleVector &sigmaValues)
{
    if ((muValues.size() != NFEATURES) || (sigmaValues.size() != NFEATURES)) // do all bounds-checking during initialization
    {
        std::cout << "SupportVectorMachine: could not add feature info because the size of mu (" << muValues.size() << ") and/or the size of sigma (" 
                  << sigmaValues.size() << ") did not match the expected number of features (" << NFEATURES << ")" << std::endl;
        return pandora::STATUS_CODE_INVALID_PARAMETER;
    }
    
    for (std::size_t i = 0; i < NFEATURES; ++i)
        m_featureInfoList[i] = {muValues[i], sigmaValues[i]};
        
    return pandora::STATUS_CODE_SUCCESS;
}
    
//------------------------------------------------------------------------------------------------------------------------------------------
  
template <std::size_t NFEATURES>  
pandora::StatusCode SupportVectorMachine<NFEATURES>::AddSupportVector(const double yAlpha, const DoubleVector &values)
{
    if (values.size() != NFEATURES) // do all bounds-checking during initialization
    {
        std::cout << "SupportVectorMachine: could not add a support vector because the number of features in the vector (" << values.size() 
                  << ") did not match the expected number of features (" << NFEATURES << ")" << std::endl;
        return pandora::STATUS_CODE_INVALID_PARAMETER;
    }
    
    FeatureArray valueArray;
    std::copy_n(values.begin(), NFEATURES, valueArray.begin());
    m_svInfoList.emplace_back(yAlpha, valueArray);
    
    return pandora::STATUS_CODE_SUCCESS;
}
    
//------------------------------------------------------------------------------------------------------------------------------------------

template <std::size_t NFEATURES>
double SupportVectorMachine<NFEATURES>::CalculateClassificationScoreImpl(const FeatureArray &features) const
{    
    if (!m_isInitialized)
    {
        std::cout << "SupportVectorMachine: could not perform classification because the SVM was uninitialized" << std::endl;
        throw pandora::StatusCodeException(pandora::STATUS_CODE_NOT_INITIALIZED); 
    }
    
    if (m_svInfoList.empty())
    {
        std::cout << "SupportVectorMachine: could not perform classification because the initialized SVM had no support vectors in the model" << std::endl;
        throw pandora::StatusCodeException(pandora::STATUS_CODE_NOT_INITIALIZED); 
    }
    
    // No need for bounds checking by construction, std::array::[] operator faster than std::array::at
    FeatureArray standardizedFeatures;
    if (m_standardizeFeatures)
    {
        for (std::size_t i = 0; i < NFEATURES; ++i)
            standardizedFeatures[i] = m_featureInfoList[i].StandardizeParameter(features[i]);
    }

    double classScore(0.0);
    for (const SupportVectorInfo &supportVectorInfo : m_svInfoList)
    {
        classScore += supportVectorInfo.m_yAlpha * 
            m_kernelFunction(supportVectorInfo.m_supportVector, (m_standardizeFeatures ? standardizedFeatures : features), m_scaleFactor);
    }

    return classScore + m_bias;
}

//------------------------------------------------------------------------------------------------------------------------------------------

template <std::size_t NFEATURES>
inline double SupportVectorMachine<NFEATURES>::LinearKernel(const FeatureArray &supportVector, const FeatureArray &features, const double scaleFactor)
{
    // No need for bounds checking by construction, std::array::[] operator faster than std::array::at
    double total(0.0);
    for (std::size_t i = 0; i < NFEATURES; ++i)
        total += supportVector[i] * features[i];
        
    return total / (scaleFactor * scaleFactor);
}

//------------------------------------------------------------------------------------------------------------------------------------------

template <std::size_t NFEATURES>
inline double SupportVectorMachine<NFEATURES>::QuadraticKernel(const FeatureArray &supportVector, const FeatureArray &features, const double scaleFactor)
{
    // No need for bounds checking by construction, std::array::[] operator faster than std::array::at
    double total(0.0);
    for (std::size_t i = 0; i < NFEATURES; ++i)
        total += supportVector[i] * features[i];
        
    total = total / (scaleFactor * scaleFactor) + 1.0;    
    return total * total;
}

//------------------------------------------------------------------------------------------------------------------------------------------

template <std::size_t NFEATURES>
inline double SupportVectorMachine<NFEATURES>::CubicKernel(const FeatureArray &supportVector, const FeatureArray &features, const double scaleFactor)
{
    // No need for bounds checking by construction, std::array::[] operator faster than std::array::at
    double total(0.0);
    for (std::size_t i = 0; i < NFEATURES; ++i)
        total += supportVector[i] * features[i];
        
    total = total / (scaleFactor * scaleFactor) + 1.0;
    return total * total * total;
}

//------------------------------------------------------------------------------------------------------------------------------------------

template <std::size_t NFEATURES>
inline double SupportVectorMachine<NFEATURES>::GaussianRbfKernel(const FeatureArray &supportVector, const FeatureArray &features, const double scaleFactor)
{
    // No need for bounds checking by construction, std::array::[] operator faster than std::array::at
    double total(0.0);
    for (std::size_t i = 0; i < NFEATURES; ++i)
        total += (supportVector[i] - features[i]) * (supportVector[i] - features[i]);
        
    return std::exp(-scaleFactor * total);
}

//------------------------------------------------------------------------------------------------------------------------------------------

template <std::size_t NFEATURES>
inline SupportVectorMachine<NFEATURES>::SupportVectorInfo::SupportVectorInfo(const double yAlpha, FeatureArray supportVector) : 
    m_yAlpha(yAlpha), 
    m_supportVector(std::move(supportVector)) 
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

template <std::size_t NFEATURES>
inline SupportVectorMachine<NFEATURES>::FeatureInfo::FeatureInfo(const double muValue, const double sigmaValue) : 
    m_muValue(muValue),
    m_sigmaValue(sigmaValue) 
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

template <std::size_t NFEATURES>
inline double SupportVectorMachine<NFEATURES>::FeatureInfo::StandardizeParameter(const double parameter) const
{
    return (parameter - m_muValue) / m_sigmaValue; // sigma size already checked
}

//------------------------------------------------------------------------------------------------------------------------------------------

template <typename T, typename TALG, typename ...Ts>
SVMFeatureTool<T, TALG, Ts...>::SVMFeatureTool() 
{ 
    static_assert((std::is_same<ReturnType, double>::value || std::is_same<ReturnType, float>::value || std::is_same<ReturnType, int>::value || 
        std::is_same<ReturnType, bool>::value), 
        "SVMFeatureTool must produce a double, float, int or bool output"); 
}
} // namespace lar_content

#endif // #ifndef LAR_SUPPORT_VECTOR_MACHINE_H
