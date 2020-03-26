/**
 *  @file   larpandoracontent/LArHelpers/DiscreteCumulativeDistributionHelper.h
 *
 *  @brief  Header file for the discrete cumulative distribution helper helper class.
 *
 *  $Log: $
 */
#ifndef LAR_DISCRETE_CUMULATIVE_DISTRIBUTION_HELPER_H
#define LAR_DISCRETE_CUMULATIVE_DISTRIBUTION_HELPER_H 1

#include "larpandoracontent/LArObjects/LArDiscreteCumulativeDistribution.h"
#include "larpandoracontent/LArHelpers/LArClusterHelper.h"
#include "Objects/CaloHit.h"

/**
 *  @brief  DiscreteCumulativeDistributionHelper class
 */

namespace lar_content
{

class LArDiscreteCumulativeDistributionHelper
{
public:
    /**
     *  @brief  Get the KS test statistic between two cumulative distributions
     *
     *  @param  distributionA the first discrete cumulatie distribution
     *  @param  distributionB the second discrete cumulatie distribution
     *
     *  @return the KS test statistic
     */
    static float CalculateKSTestStatistic(const DiscreteCumulativeDistribution &distributionA, 
            const DiscreteCumulativeDistribution &distributionB);

    /**
     *  @brief  Find the Y value for a particular X value
     *
     *  @param  distribution the discrete cumulatie distribution
     *  @param  x the X value
     *
     *  @return the Y value
     */
    static float FindY(const DiscreteCumulativeDistribution &distribution, const float &x);

    /**
     *  @brief  Fill and create a cumulatie distribution using a CaloHitList
     *
     *  @param  caloHitList the caloHitList
     *  @param  distribution the empty discrete cumulative distribution
     *
     *  @return
     */
    static void CreateDistributionFromCaloHits(pandora::CaloHitList caloHitList, 
            DiscreteCumulativeDistribution &distribution);

    /**
     *  @brief  Full computation of PValue(trueKs>measuredKs|distA==distB)
     *
     *  @param  distributionA the first discrete cumulatie distribution
     *  @param  distributionB the second discrete cumulatie distribution
     *
     *  @return PValue(trueKs>measuredKs|distA==distB)
     */
    static float CalculatePValueWithKSTestStatistic(const DiscreteCumulativeDistribution &distributionA, 
	    const DiscreteCumulativeDistribution &distributionB);

    /**
     *  @brief  Computation of one term of the sum in PValue(trueKs>measuredKs|distA==distB)
     *                                                                                                                                                                                                   
     *  @param  indice of the sum term
     *  @param  measuredKs
     *  @param  distributionA the first discrete cumulatie distribution                                                                                                                                       *  @param  distributionB the second discrete cumulatie distribution
     *
     *  @return sum term of PValue(trueKs>measuredKs|distA==distB)
     */
    static float CalculatePValueSumTerm(const int i, const float &ks, 
	    const DiscreteCumulativeDistribution &distributionA, const DiscreteCumulativeDistribution &distributionB);

};

} // namespace lar_content
#endif // #ifndef LAR_DISCRETE_CUMULATIVE_DISTRIBUTION_HELPER_H
