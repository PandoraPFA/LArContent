/**
 *  @file   LArContent/include/LArThreeDReco/LArHitCreation/ShowerHitCreationTool.h
 * 
 *  @brief  Header file for the shower hit creation tool.
 * 
 *  $Log: $
 */
#ifndef SHOWER_HIT_CREATION_TOOL_H
#define SHOWER_HIT_CREATION_TOOL_H 1

#include "LArThreeDReco/LArHitCreation/ThreeDHitCreationAlgorithm.h"

namespace lar
{

/**
 *  @brief  ShowerHitCreationTool class
 */
class ShowerHitCreationTool : public HitCreationTool
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
        pandora::CaloHitList &newThreeDHits);

private:
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
     *  @brief  Create three dimensional hits, using a list of input two dimensional hits and the hits (contained in the same particle)
     *          from the other two views
     * 
     *  @param  pAlgorithm address of the calling algorithm
     *  @param  inputTwoDHits the list of input two dimensional hits
     *  @param  caloHitList1 the first 
     *  @param  caloHitList2 the second 
     *  @param  newThreeDHits to receive the new three dimensional hits
     */
    void CreateThreeDHits(ThreeDHitCreationAlgorithm *pAlgorithm, const pandora::CaloHitList &inputTwoDHits, const pandora::CaloHitList &caloHitList1,
        const pandora::CaloHitList &caloHitList2, pandora::CaloHitList &newThreeDHits) const;

    /**
     *  @brief  Filter a list of calo hits to find those within a specified tolerance of a give x position
     * 
     *  @param  x the x position
     *  @param  xTolerance the x tolerance
     *  @param  inputCaloHitList the input calo hit list
     *  @param  outputCaloHitList to receive the output calo hit list
     */
    void FilterCaloHits(const float x, const float xTolerance, const pandora::CaloHitList &inputCaloHitList, pandora::CaloHitList &outputCaloHitList) const;

    /**
     *  @brief  Get the three dimensional position for to a two dimensional calo hit, using the hit and a list of candidate matched
     *          hits in the other two views
     * 
     *  @param  pCaloHit2D address of the two dimensional calo hit
     *  @param  caloHitList1 the list of candidate hits in view 1
     *  @param  caloHitList2 the list of candidate hits in view 2
     *  @param  position3D to receive the three dimensional position
     *  @param  chiSquared to receive the chi squared value
     */
    void GetPosition3D(const pandora::CaloHit *const pCaloHit2D, const pandora::CaloHitList &caloHitList1, const pandora::CaloHitList &caloHitList2,
        pandora::CartesianVector &position3D, float &chiSquared) const;

    /**
     *  @brief  Get the three dimensional position for to a two dimensional calo hit, using the hit and a pair of candidate matched
     *          hits in the other two views
     * 
     *  @param  pCaloHit2D address of the two dimensional calo hit
     *  @param  pCaloHit1 address of the candidate calo hit in view 1
     *  @param  pCaloHit2 address of the candidate calo hit in view 2
     *  @param  position3D to receive the three dimensional position
     *  @param  chiSquared to receive the chi squared value
     */
    void GetPosition3D(const pandora::CaloHit *const pCaloHit2D, const pandora::CaloHit *const pCaloHit1, const pandora::CaloHit *const pCaloHit2,
        pandora::CartesianVector &position3D, float &chiSquared) const;

    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    float           m_xTolerance;               ///< The x tolerance to use when looking for associated calo hits between views
    bool            m_useDeltaXCorrection;      ///< Whether to add a term to chi squared accounting for hit combination delta x values
    float           m_sigmaX;                   ///< Resolution in x dimension, used for delta x correction to chi squared
    bool            m_useChiSquaredApproach;    ///< Whether to obtain y, z positions via chi2 approach, or projected position approach
    float           m_chiSquaredCut;            ///< The chi squared cut (accept only values below the cut value)
};

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::AlgorithmTool *ShowerHitCreationTool::Factory::CreateAlgorithmTool() const
{
    return new ShowerHitCreationTool();
}

} // namespace lar

#endif // #ifndef SHOWER_HIT_CREATION_TOOL_H
