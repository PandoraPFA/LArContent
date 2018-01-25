/**
 *  @file   larpandoracontent/LArTrackShowerId/PfoCharacterisationBaseAlgorithm.h
 * 
 *  @brief  Header file for the pfo characterisation base algorithm class.
 * 
 *  $Log: $
 */
#ifndef LAR_PFO_CHARACTERISATION_BASE_ALGORITHM_H
#define LAR_PFO_CHARACTERISATION_BASE_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

namespace lar_content
{

/**
 *  @brief  PfoCharacterisationBaseAlgorithm class
 */
class PfoCharacterisationBaseAlgorithm : public pandora::Algorithm
{
public:
    /**
     *  @brief  Default constructor
     */
    PfoCharacterisationBaseAlgorithm();

    /**
     *  @brief  Destructor
     */
    virtual ~PfoCharacterisationBaseAlgorithm();

protected:
    pandora::StatusCode Run();

    /**
     *  @brief  Whether pfo is identified as a clear track using its three clusters
     *
     *  @param  pPfo address of the relevant pfo
     *
     *  @return boolean
     */
    virtual bool IsClearTrack3x2D(const pandora::ParticleFlowObject *const pPfo) const;

    /**
     *  @brief  Whether pfo is identified as a clear track
     *
     *  @param  pPfo address of the relevant pfo
     * 
     *  @return boolean
     */
    virtual bool IsClearTrack(const pandora::ParticleFlowObject *const pPfo) const = 0;

    /**
     *  @brief  Whether cluster is identified as a clear track
     *
     *  @param  pCluster address of the relevant cluster
     * 
     *  @return boolean
     */
    virtual bool IsClearTrack(const pandora::Cluster *const pCluster) const = 0;

    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    std::string             m_trackPfoListName;             ///< The track pfo list name
    std::string             m_showerPfoListName;            ///< The shower pfo list name
    pandora::StringVector   m_inputPfoListNames;            ///< The names of the input pfo lists

    bool                    m_updateClusterIds;             ///< Whether to update daughter cluster particle id labels to match pfo id
    bool                    m_postBranchAddition;           ///< Whether to use configuration for shower clusters post branch addition
    bool                    m_useThreeDInformation;         ///< Whether to use PFO and 3D information or clusters for characterisation
    unsigned int            m_minTrackLikeViews;            ///< The minimum number of track-like views to declare a pfo as track-like
};

} // namespace lar_content

#endif // #ifndef LAR_PFO_CHARACTERISATION_BASE_ALGORITHM_H
