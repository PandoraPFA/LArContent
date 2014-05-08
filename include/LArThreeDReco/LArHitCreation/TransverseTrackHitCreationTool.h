/**
 *  @file   LArContent/include/LArThreeDReco/LArHitCreation/TransverseTrackHitCreationTool.h
 * 
 *  @brief  Header file for the transverse track hit creation tool.
 * 
 *  $Log: $
 */
#ifndef TRANSVERSE_TRACK_HIT_CREATION_TOOL_H
#define TRANSVERSE_TRACK_HIT_CREATION_TOOL_H 1

#include "LArThreeDReco/LArHitCreation/ThreeDHitCreationAlgorithm.h"

namespace lar
{

/**
 *  @brief  TransverseTrackHitCreationTool class
 */
class TransverseTrackHitCreationTool : public HitCreationTool
{
public:
    /**
     *  @brief  Factory class for instantiating algorithm tool
     */
    class Factory : public pandora::AlgorithmToolFactory
    {
    public:
        pandora::AlgorithmTool *CreateAlgorithmTool() const;
    };

    void Run(ThreeDHitCreationAlgorithm *pAlgorithm, const pandora::ParticleFlowObject *const pPfo, const pandora::CaloHitList &inputTwoDHits,
        pandora::CaloHitList &newThreeDHits, pandora::CaloHitList &omittedTwoDHits);

private:
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    /**
     *  @brief  Get the primary u, v and w clusters in a provided pfo
     * 
     *  @param  pPfo address of the pfo
     *  @param  pClusterU to receive the address of cluster u
     *  @param  pClusterV to receive the address of cluster v
     *  @param  pClusterW to receive the address of cluster w
     */
    void GetClusters(const pandora::ParticleFlowObject *const pPfo, pandora::Cluster *&pClusterU, pandora::Cluster *&pClusterV,
        pandora::Cluster *&pClusterW) const;

    /**
     *  @brief  Create three dimensional hits, using an input list of two dimensional hits and two associated sliding fit results
     * 
     *  @param  pAlgorithm address of the calling algorithm
     *  @param  inputTwoDHits the list of input two dimensional hits
     *  @param  fitResult1 the first sliding fit result
     *  @param  fitResult2 the second sliding fit result
     *  @param  newThreeDHits to receive the new three dimensional hits
     *  @param  omittedTwoDHits to receive the two dimensional hits for which no three dimensional hits could be created
     */
    void CreateThreeDHits(ThreeDHitCreationAlgorithm *pAlgorithm, const pandora::CaloHitList &inputTwoDHits, const TwoDSlidingFitResult &fitResult1,
        const TwoDSlidingFitResult &fitResult2, pandora::CaloHitList &newThreeDHits, pandora::CaloHitList &omittedTwoDHits) const;

    /**
     *  @brief  Get the three dimensional position using a provided two dimensional calo hit and sliding linear fits in the other two views
     * 
     *  @param  pCaloHit2D address of the two dimensional calo hit
     *  @param  fitResult1 the first sliding fit result
     *  @param  fitResult2 the second sliding fit result
     *  @param  position3D to receive the three dimensional position
     *  @param  chiSquared to receive the chi squared value
     */
    void GetPosition3D(const pandora::CaloHit *const pCaloHit2D, const TwoDSlidingFitResult &fitResult1, const TwoDSlidingFitResult &fitResult2,
        pandora::CartesianVector &position3D, float &chiSquared) const;

    unsigned int    m_slidingFitWindow;         ///< The layer window for the sliding linear fits
    bool            m_useChiSquaredApproach;    ///< Whether to obtain y, z positions via chi2 approach, or projected position approach
    float           m_sigmaFitMultiplier;       ///< The multiplier from sigma hit (i.e. sigma uvw from transformation calculator) to sigma fit
    float           m_chiSquaredCut;            ///< The chi squared cut (accept only values below the cut value)
};

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::AlgorithmTool *TransverseTrackHitCreationTool::Factory::CreateAlgorithmTool() const
{
    return new TransverseTrackHitCreationTool();
}

} // namespace lar

#endif // #ifndef TRANSVERSE_TRACK_HIT_CREATION_TOOL_H
