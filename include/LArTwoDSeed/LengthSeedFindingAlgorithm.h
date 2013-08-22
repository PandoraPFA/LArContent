/**
 *  @file   LArContent/include/LArTwoDSeed/LengthSeedFindingAlgorithm.h
 * 
 *  @brief  Header file for the length seed finding algorithm class.
 * 
 *  $Log: $
 */
#ifndef LAR_LENGTH_SEED_FINDING_ALGORITHM_H
#define LAR_LENGTH_SEED_FINDING_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

#include "SeedFindingBaseAlgorithm.h"

namespace lar
{

/**
 *  @brief  LengthSeedFindingAlgorithm class
 */
class LengthSeedFindingAlgorithm : public SeedFindingBaseAlgorithm
{
public:
    /**
     *  @brief  Factory class for instantiating algorithm
     */
    class Factory : public pandora::AlgorithmFactory
    {
    public:
        pandora::Algorithm *CreateAlgorithm() const;
    };

private:
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);
    void GetSeedClusterList(const pandora::ClusterVector &candidateClusters, pandora::ClusterList &seedClusterList) const;

    int                 m_initialLengthCut;         ///< 
    int                 m_finalLengthCut;           ///< 
    int                 m_initialChangeIter;        ///< 
    int                 m_finalChangeIter;          ///< 
};

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::Algorithm *LengthSeedFindingAlgorithm::Factory::CreateAlgorithm() const
{
    return new LengthSeedFindingAlgorithm();
}

} // namespace lar

#endif // #ifndef LAR_LENGTH_SEED_FINDING_ALGORITHM_H
