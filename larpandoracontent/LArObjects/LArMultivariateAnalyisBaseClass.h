/**
 *  @file   larpandoracontent/LArObjects/MultivariateAnalyisBaseClass.h
 *
 *  @brief  Header file for the lar multivariate analysis base class.
 *
 *  $Log: $
 */
#ifndef LAR_MULTIVARIATE_ANALYSIS_BASE_CLASS_H
#define LAR_MULTIVARIATE_ANALYSIS_BASE_CLASS_H 1

namespace lar_content
{

/**
 *  @brief  MultivariateAnalyisBaseClass class
 */
class MultivariateAnalyisBaseClass
{
public:
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
};

} // namespace lar_content

#endif // #ifndef LAR_MULTIVARIATE_ANALYSIS_BASE_CLASS_H
