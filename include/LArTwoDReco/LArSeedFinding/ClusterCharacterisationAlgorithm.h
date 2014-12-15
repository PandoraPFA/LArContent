/**
 *  @file   LArContent/include/LArTwoDReco/LArSeedFinding/ClusterCharacterisationAlgorithm.h
 * 
 *  @brief  Header file for the cluster characterisation algorithm class.
 * 
 *  $Log: $
 */
#ifndef LAR_CLUSTER_CHARACTERISATION_ALGORITHM_H
#define LAR_CLUSTER_CHARACTERISATION_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

namespace lar_content
{

/**
 *  @brief  ClusterCharacterisationAlgorithm class
 */
class ClusterCharacterisationAlgorithm : public pandora::Algorithm
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

    /**
     *  @brief  Default constructor
     */
    ClusterCharacterisationAlgorithm();

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

    unsigned int            m_slidingFitWindow;             ///< The layer window for the sliding linear fits
    unsigned int            m_minHitsInCluster;             ///< The minimum number of hits in a clear track candidate
    float                   m_maxLayerGapFraction;          ///< The maximum (sliding fit) layer gap fraction for a clear track candidate
    float                   m_maxWidthPerUnitLength;        ///< The maximum width per unit length for a clear track candidate
};

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::Algorithm *ClusterCharacterisationAlgorithm::Factory::CreateAlgorithm() const
{
    return new ClusterCharacterisationAlgorithm();
}

} // namespace lar_content

#endif // #ifndef LAR_CLUSTER_CHARACTERISATION_ALGORITHM_H
