/**
 *  @file   LArContent/include/LArTwoDReco/LArSeedFinding/CheatingClusterCharacterisationAlgorithm.h
 * 
 *  @brief  Header file for the cheating cluster characterisation algorithm class.
 * 
 *  $Log: $
 */
#ifndef LAR_CHEATING_CLUSTER_CHARACTERISATION_ALGORITHM_H
#define LAR_CHEATING_CLUSTER_CHARACTERISATION_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

namespace lar_content
{

/**
 *  @brief  CheatingClusterCharacterisationAlgorithm class
 */
class CheatingClusterCharacterisationAlgorithm : public pandora::Algorithm
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

    /**
     *  @brief  Whether cluster is identified as a clear track
     *
     *  @param  pCluster address of the relevant cluster
     * 
     *  @return boolean
     */
    bool IsClearTrack(const pandora::Cluster *const pCluster) const;

    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    pandora::StringVector   m_inputClusterListNames;        ///< The names of the input cluster lists
};

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::Algorithm *CheatingClusterCharacterisationAlgorithm::Factory::CreateAlgorithm() const
{
    return new CheatingClusterCharacterisationAlgorithm();
}

} // namespace lar_content

#endif // #ifndef LAR_CHEATING_CLUSTER_CHARACTERISATION_ALGORITHM_H
