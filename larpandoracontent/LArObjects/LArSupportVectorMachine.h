/**
 *  @file   larpandoracontent/LArObjects/LArSupportVectorMachine.h
 *
 *  @brief  Header file for the lar support vector machine class.
 *
 *  $Log: $
 */
#ifndef LAR_SUPPORT_VECTOR_MACHINE_H
#define LAR_SUPPORT_VECTOR_MACHINE_H 1

#include "larpandoracontent/LArHelpers/LArMvaHelper.h"

#include "larpandoracontent/LArObjects/LArMvaInterface.h"

#include "Helpers/XmlHelper.h"
#include "Pandora/StatusCodes.h"

#include <functional>
#include <map>
#include <vector>

//------------------------------------------------------------------------------------------------------------------------------------------

namespace lar_content
{

/**
 *  @brief  SupportVectorMachine class
 */
class SupportVectorMachine : public MvaInterface
{
public:
    typedef std::function<double(const LArMvaHelper::MvaFeatureVector &, const LArMvaHelper::MvaFeatureVector &, const double)> KernelFunction;

    /**
     *  @brief  KernelType enum
     */
    enum KernelType
    {
        USER_DEFINED = 0,
        LINEAR = 1,
        QUADRATIC = 2,
        CUBIC = 3,
        GAUSSIAN_RBF = 4
    };

    /**
     *  @brief  Default constructor.
     */
    SupportVectorMachine();

    /**
     *  @brief  Initialize the svm using a serialized model.
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

    /**
     *  @brief  Query whether this svm is initialized
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
        SupportVectorInfo(const double yAlpha, LArMvaHelper::MvaFeatureVector supportVector);

        double m_yAlpha;                                ///< The alpha-value multiplied by the y-value for the support vector
        LArMvaHelper::MvaFeatureVector m_supportVector; ///< The support vector
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
         *  @brief  Default constructor to allow default-construction of (uninitialized) svms.
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

        double m_muValue;    ///< The average value of this feature
        double m_sigmaValue; ///< The stddev of this feature
    };

    typedef std::vector<SupportVectorInfo> SVInfoList;
    typedef std::vector<FeatureInfo> FeatureInfoVector;

    typedef std::map<KernelType, KernelFunction> KernelMap;

    bool m_isInitialized; ///< Whether this svm has been initialized

    bool m_enableProbability; ///< Whether to enable probability calculations
    double m_probAParameter;  ///< The first-order score coefficient for mapping to a probability using the logistic function
    double m_probBParameter;  ///< The score offset parameter for mapping to a probability using the logistic function

    bool m_standardizeFeatures; ///< Whether to standardize the features
    unsigned int m_nFeatures;   ///< The number of features
    double m_bias;              ///< The bias term
    double m_scaleFactor;       ///< The kernel scale factor

    SVInfoList m_svInfoList;             ///< The list of SupportVectorInfo objects
    FeatureInfoVector m_featureInfoList; ///< The list of FeatureInfo objects

    KernelType m_kernelType;         ///< The kernel type
    KernelFunction m_kernelFunction; ///< The kernel function
    KernelMap m_kernelMap;           ///< Map from the kernel types to the kernel functions

    /**
     *  @brief  Read the svm parameters from an xml file
     *
     *  @param  svmFileName the sml file name
     *  @param  svmName the name of the svm
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
    double CalculateClassificationScoreImpl(const LArMvaHelper::MvaFeatureVector &features) const;

    /**
     *  @brief  An inhomogeneous quadratic kernel
     *
     *  @param  supportVector the support vector
     *  @param  features the features
     *  @param  scale factor the scale factor
     *
     *  @return result of the kernel operation
     */
    static double QuadraticKernel(
        const LArMvaHelper::MvaFeatureVector &supportVector, const LArMvaHelper::MvaFeatureVector &features, const double scaleFactor = 1.);

    /**
     *  @brief  An inhomogeneous cubic kernel
     *
     *  @param  supportVector the support vector
     *  @param  features the features
     *  @param  scale factor the scale factor
     *
     *  @return result of the kernel operation
     */
    static double CubicKernel(
        const LArMvaHelper::MvaFeatureVector &supportVector, const LArMvaHelper::MvaFeatureVector &features, const double scaleFactor = 1.);

    /**
     *  @brief  A linear kernel
     *
     *  @param  supportVector the support vector
     *  @param  features the features
     *  @param  scale factor the scale factor
     *
     *  @return result of the kernel operation
     */
    static double LinearKernel(
        const LArMvaHelper::MvaFeatureVector &supportVector, const LArMvaHelper::MvaFeatureVector &features, const double scaleFactor = 1.);

    /**
     *  @brief  A gaussian RBF kernel
     *
     *  @param  supportVector the support vector
     *  @param  features the features
     *  @param  scale factor the scale factor
     *
     *  @return result of the kernel operation
     */
    static double GaussianRbfKernel(
        const LArMvaHelper::MvaFeatureVector &supportVector, const LArMvaHelper::MvaFeatureVector &features, const double scaleFactor = 1.);
};

//------------------------------------------------------------------------------------------------------------------------------------------

inline bool SupportVectorMachine::Classify(const LArMvaHelper::MvaFeatureVector &features) const
{
    return (this->CalculateClassificationScoreImpl(features) > 0.);
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline double SupportVectorMachine::CalculateClassificationScore(const LArMvaHelper::MvaFeatureVector &features) const
{
    return this->CalculateClassificationScoreImpl(features);
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline double SupportVectorMachine::CalculateProbability(const LArMvaHelper::MvaFeatureVector &features) const
{
    if (!m_enableProbability)
    {
        std::cout << "LArSupportVectorMachine: cannot calculate probabilities for this SVM" << std::endl;
        throw pandora::STATUS_CODE_NOT_INITIALIZED;
    }

    // Use the logistic function to map the linearly-transformed score on the interval (-inf,inf) to a probability on [0,1] - the two free
    // parameters in the linear transformation are trained such that the logistic map produces an accurate probability
    const double scaledScore = m_probAParameter * this->CalculateClassificationScoreImpl(features) + m_probBParameter;

    return 1. / (1. + std::exp(scaledScore));
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

inline double SupportVectorMachine::LinearKernel(
    const LArMvaHelper::MvaFeatureVector &supportVector, const LArMvaHelper::MvaFeatureVector &features, const double scaleFactor)
{
    const double denominator(scaleFactor * scaleFactor);
    if (denominator < std::numeric_limits<double>::epsilon())
        throw pandora::StatusCodeException(pandora::STATUS_CODE_INVALID_PARAMETER);

    double total(0.);
    for (unsigned int i = 0; i < features.size(); ++i)
        total += supportVector.at(i).Get() * features.at(i).Get();

    return total / denominator;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline double SupportVectorMachine::QuadraticKernel(
    const LArMvaHelper::MvaFeatureVector &supportVector, const LArMvaHelper::MvaFeatureVector &features, const double scaleFactor)
{
    const double denominator(scaleFactor * scaleFactor);
    if (denominator < std::numeric_limits<double>::epsilon())
        throw pandora::StatusCodeException(pandora::STATUS_CODE_INVALID_PARAMETER);

    double total(0.);
    for (unsigned int i = 0; i < features.size(); ++i)
        total += supportVector.at(i).Get() * features.at(i).Get();

    total = total / denominator + 1.;
    return total * total;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline double SupportVectorMachine::CubicKernel(
    const LArMvaHelper::MvaFeatureVector &supportVector, const LArMvaHelper::MvaFeatureVector &features, const double scaleFactor)
{
    const double denominator(scaleFactor * scaleFactor);
    if (denominator < std::numeric_limits<double>::epsilon())
        throw pandora::StatusCodeException(pandora::STATUS_CODE_INVALID_PARAMETER);

    double total(0.);
    for (unsigned int i = 0; i < features.size(); ++i)
        total += supportVector.at(i).Get() * features.at(i).Get();

    total = total / denominator + 1.;
    return total * total * total;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline double SupportVectorMachine::GaussianRbfKernel(
    const LArMvaHelper::MvaFeatureVector &supportVector, const LArMvaHelper::MvaFeatureVector &features, const double scaleFactor)
{
    double total(0.);
    for (unsigned int i = 0; i < features.size(); ++i)
        total += (supportVector.at(i).Get() - features.at(i).Get()) * (supportVector.at(i).Get() - features.at(i).Get());

    return std::exp(-scaleFactor * total);
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline SupportVectorMachine::SupportVectorInfo::SupportVectorInfo(const double yAlpha, LArMvaHelper::MvaFeatureVector supportVector) :
    m_yAlpha(yAlpha), m_supportVector(std::move(supportVector))
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline SupportVectorMachine::FeatureInfo::FeatureInfo(const double muValue, const double sigmaValue) :
    m_muValue(muValue), m_sigmaValue(sigmaValue)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline SupportVectorMachine::FeatureInfo::FeatureInfo() : m_muValue(0.), m_sigmaValue(0.)
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
