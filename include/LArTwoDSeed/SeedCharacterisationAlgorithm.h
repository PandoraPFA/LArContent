/**
 *  @file   LArContent/include/LArTwoDSeed/SeedCharacterisationAlgorithm.h
 * 
 *  @brief  Header file for the seed characterisation algorithm class.
 * 
 *  $Log: $
 */
#ifndef LAR_SEED_CHARACTERISATION_ALGORITHM_H
#define LAR_SEED_CHARACTERISATION_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

#include "SeedBranchGrowingAlgorithm.h"

namespace lar
{

/**
 *  @brief  SeedCharacterisationAlgorithm class
 */
class SeedCharacterisationAlgorithm : public SeedBranchGrowingAlgorithm
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
    pandora::StatusCode Run();
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    /**
     *  @brief  Get the seed association list
     * 
     *  @param  particleSeedVector the list of all particle seeds
     *  @param  seedAssociationList to receive the populated seed association list
     */
    void GetSeedAssociationList(const pandora::ClusterVector &particleSeedVector, SeedAssociationList &seedAssociationList) const;

    std::string         m_trackSeedClusterListName;     ///< The track seed cluster list name
    std::string         m_showerSeedClusterListName;    ///< The shower seed cluster list name
};

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::Algorithm *SeedCharacterisationAlgorithm::Factory::CreateAlgorithm() const
{
    return new SeedCharacterisationAlgorithm();
}

} // namespace lar

#endif // #ifndef LAR_SEED_CHARACTERISATION_ALGORITHM_H
