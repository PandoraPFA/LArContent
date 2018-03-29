/**
 *  @file   larpandoracontent/LArObjects/LArMvaInterface.h
 *
 *  @brief  Header file for the lar multivariate analysis interface class.
 *
 *  $Log: $
 */
#ifndef LAR_MVA_INTERFACE_H
#define LAR_MVA_INTERFACE_H 1

#include "Pandora/PandoraInputTypes.h"

#include <vector>

namespace lar_content
{

/**
 *  @brief  MvaTypes class
 */
class MvaTypes
{
public:
    typedef pandora::PandoraInputType<double> MvaFeature;
    typedef std::vector<MvaFeature> MvaFeatureVector;
};

//------------------------------------------------------------------------------------------------------------------------------------------

/**
 *  @brief  MvaInterface class
 */
class MvaInterface
{
public:
    /**
     *  @brief  Classify the set of input features based on the trained model
     *
     *  @param  features the input features
     *
     *  @return the classification
     */
    virtual bool Classify(const MvaTypes::MvaFeatureVector &features) const = 0;

    /**
     *  @brief  Calculate the classification score for a set of input features, based on the trained model
     *
     *  @param  features the input features
     *
     *  @return the classification score
     */
    virtual double CalculateClassificationScore(const MvaTypes::MvaFeatureVector &features) const = 0;

    /**
     *  @brief  Calculate the classification probability for a set of input features, based on the trained model
     *
     *  @param  features the input features
     *
     *  @return the classification probability
     */
    virtual double CalculateProbability(const MvaTypes::MvaFeatureVector &features) const = 0;

    /**
     *  @brief  Destructor
     */
    virtual ~MvaInterface() = default;
};

} // namespace lar_content

#endif // #ifndef LAR_MVA_INTERFACE_H
