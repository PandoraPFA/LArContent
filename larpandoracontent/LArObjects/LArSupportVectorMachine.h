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
#include "Pandora/StatusCodes.h"

#include <functional>
#include <map>
#include <vector>

namespace lar_content
{

/**
 *  @brief  SVMFeatureTool class template
 */
template <typename ...Ts>
class SVMFeatureTool : public pandora::AlgorithmTool
{
public:
    typedef std::vector<double> DoubleVector;
    typedef std::vector<SVMFeatureTool<Ts...> *> FeatureToolVector;

    /**
     *  @brief  Default constructor.
     */
    SVMFeatureTool() = default;

    /**
     *  @brief  Run the algorithm tool
     *
     *  @param  pAlg address of the calling algorithm
     *  @param  args arguments to pass to the tool
     *
     *  @return the feature
     */
    virtual void Run(DoubleVector &featureVector, Ts... args) = 0;
};

template <typename ...Ts>
using SVMFeatureToolVector = std::vector<SVMFeatureTool<Ts...> *>;

//------------------------------------------------------------------------------------------------------------------------------------------

/**
 *  @brief  SupportVectorMachine class
 */
class SupportVectorMachine
{
public:
    typedef std::vector<double> DoubleVector;
    typedef std::function<double(const DoubleVector &, const DoubleVector &, const double)> KernelFunction;

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

    /**
     *  @brief  Make a classification for a set of input features, based on the trained model
     *
     *  @param  features the input features
     *
     *  @return the predicted boolean class
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
     *  @brief  Query whether this SVM is initialized
     *
     *  @return whether the model is initialized
     */
    bool IsInitialized() const;

    /**
     *  @brief  Get the number of features
     *
     *  @return the number of features
     */
    unsigned int GetNFeatures() const;

    /**
     *  @brief  Set the kernel function to use
     *
     *  @param  kernelFunction the kernel function
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
         *  @param  supportVector the support vector, passed by value then uses move semantics for efficiency
         */
        SupportVectorInfo(const double yAlpha, DoubleVector supportVector);

        double       m_yAlpha;        ///< The alpha-value multiplied by the y-value for the support vector
        DoubleVector m_supportVector; ///< The support vector
    };

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
        FeatureInfo();

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

    typedef std::vector<SupportVectorInfo> SVInfoList;
    typedef std::vector<FeatureInfo>       FeatureInfoVector;

    typedef std::map<KernelType, KernelFunction> KernelMap;

    bool              m_isInitialized;       ///< Whether this SVM has been initialized

    bool              m_standardizeFeatures; ///< Whether to standardize the features
    unsigned int      m_nFeatures;           ///< The number of features
    double            m_bias;                ///< The bias term
    double            m_scaleFactor;         ///< The kernel scale factor

    SVInfoList        m_svInfoList;          ///< The list of SupportVectorInfo objects
    FeatureInfoVector m_featureInfoList;     ///< The list of FeatureInfo objects

    KernelType        m_kernelType;          ///< The kernel type
    KernelFunction    m_kernelFunction;      ///< The kernel function
    KernelMap         m_kernelMap;           ///< Map from the kernel types to the kernel functions

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
     *  @brief  Implementation method for calculating the classification score using the trained model.
     *
     *  @param  features the vector of features
     *
     *  @return the classification score
     */
    double CalculateClassificationScoreImpl(const DoubleVector &features) const;

    /**
     *  @brief  An inhomogeneous quadratic kernel
     *
     *  @param  supportVector the support vector
     *  @param  features the features
     *  @param  scale factor the scale factor
     *
     *  @return result of the kernel operation
     */
    static double QuadraticKernel(const DoubleVector &supportVector, const DoubleVector &features, const double scaleFactor = 1.);

    /**
     *  @brief  An inhomogeneous cubic kernel
     *
     *  @param  supportVector the support vector
     *  @param  features the features
     *  @param  scale factor the scale factor
     *
     *  @return result of the kernel operation
     */
    static double CubicKernel(const DoubleVector &supportVector, const DoubleVector &features, const double scaleFactor = 1.);

    /**
     *  @brief  A linear kernel
     *
     *  @param  supportVector the support vector
     *  @param  features the features
     *  @param  scale factor the scale factor
     *
     *  @return result of the kernel operation
     */
    static double LinearKernel(const DoubleVector &supportVector, const DoubleVector &features, const double scaleFactor = 1.);

    /**
     *  @brief  A gaussian RBF kernel
     *
     *  @param  supportVector the support vector
     *  @param  features the features
     *  @param  scale factor the scale factor
     *
     *  @return result of the kernel operation
     */
    static double GaussianRbfKernel(const DoubleVector &supportVector, const DoubleVector &features, const double scaleFactor = 1.);
};

//------------------------------------------------------------------------------------------------------------------------------------------

inline bool SupportVectorMachine::Classify(const DoubleVector &features) const
{
    return (this->CalculateClassificationScoreImpl(features) > 0.);
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline double SupportVectorMachine::CalculateClassificationScore(const DoubleVector &features) const
{
    return this->CalculateClassificationScoreImpl(features);
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline bool SupportVectorMachine::IsInitialized() const
{
    return m_isInitialized;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline unsigned int SupportVectorMachine::GetNFeatures() const
{
    return m_nFeatures;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline void SupportVectorMachine::SetKernelFunction(KernelFunction kernelFunction)
{
    m_kernelFunction = std::move(kernelFunction);
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline double SupportVectorMachine::LinearKernel(const DoubleVector &supportVector, const DoubleVector &features, const double scaleFactor)
{
    const double denominator(scaleFactor * scaleFactor);
    if (denominator < std::numeric_limits<double>::epsilon())
        throw pandora::StatusCodeException(pandora::STATUS_CODE_INVALID_PARAMETER);

    double total(0.);
    for (unsigned int i = 0; i < features.size(); ++i)
        total += supportVector.at(i) * features.at(i);

    return total / denominator;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline double SupportVectorMachine::QuadraticKernel(const DoubleVector &supportVector, const DoubleVector &features, const double scaleFactor)
{
    const double denominator(scaleFactor * scaleFactor);
    if (denominator < std::numeric_limits<double>::epsilon())
        throw pandora::StatusCodeException(pandora::STATUS_CODE_INVALID_PARAMETER);

    double total(0.);
    for (unsigned int i = 0; i < features.size(); ++i)
        total += supportVector.at(i) * features.at(i);

    total = total / denominator + 1.;
    return total * total;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline double SupportVectorMachine::CubicKernel(const DoubleVector &supportVector, const DoubleVector &features, const double scaleFactor)
{
    const double denominator(scaleFactor * scaleFactor);
    if (denominator < std::numeric_limits<double>::epsilon())
        throw pandora::StatusCodeException(pandora::STATUS_CODE_INVALID_PARAMETER);

    double total(0.);
    for (unsigned int i = 0; i < features.size(); ++i)
        total += supportVector.at(i) * features.at(i);

    total = total / denominator + 1.;
    return total * total * total;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline double SupportVectorMachine::GaussianRbfKernel(const DoubleVector &supportVector, const DoubleVector &features, const double scaleFactor)
{
    double total(0.);
    for (unsigned int i = 0; i < features.size(); ++i)
        total += (supportVector.at(i) - features.at(i)) * (supportVector.at(i) - features.at(i));

    return std::exp(-scaleFactor * total);
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline SupportVectorMachine::SupportVectorInfo::SupportVectorInfo(const double yAlpha, DoubleVector supportVector) :
    m_yAlpha(yAlpha),
    m_supportVector(std::move(supportVector))
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline SupportVectorMachine::FeatureInfo::FeatureInfo(const double muValue, const double sigmaValue) :
    m_muValue(muValue),
    m_sigmaValue(sigmaValue)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline SupportVectorMachine::FeatureInfo::FeatureInfo() :
    m_muValue(0.),
    m_sigmaValue(0.)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline double SupportVectorMachine::FeatureInfo::StandardizeParameter(const double parameter) const
{
    if (m_sigmaValue < std::numeric_limits<double>::epsilon())
        throw pandora::StatusCodeException(pandora::STATUS_CODE_INVALID_PARAMETER);

    return (parameter - m_muValue) / m_sigmaValue;
}

} // namespace lar_content

#endif // #ifndef LAR_SUPPORT_VECTOR_MACHINE_H