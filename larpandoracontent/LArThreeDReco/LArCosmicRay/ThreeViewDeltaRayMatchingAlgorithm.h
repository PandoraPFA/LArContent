/**
 *  @file   larpandoracontent/LArThreeDReco/LArCosmicRay/ThreeViewDeltaRayMatching.h
 *
 *  @brief  Header file for the three view delta ray matching class
 *
 *  $Log: $
 */
#ifndef LAR_THREE_VIEW_DELTA_RAY_MATCHING_ALGORITHM_H
#define LAR_THREE_VIEW_DELTA_RAY_MATCHING_ALGORITHM_H 1

#include "Pandora/Algorithm.h"
#include "Pandora/AlgorithmTool.h"

#include "larpandoracontent/LArObjects/LArTrackOverlapResult.h"

#include "larpandoracontent/LArThreeDReco/LArThreeDBase/NViewMatchingAlgorithm.h"
#include "larpandoracontent/LArThreeDReco/LArThreeDBase/ThreeViewMatchingControl.h"

namespace lar_content
{
    
class DeltaRayTensorTool;

//------------------------------------------------------------------------------------------------------------------------------------------

/**
 *  @brief  ThreeViewDeltaRayMatchingAlgorithm class
 */
class ThreeViewDeltaRayMatchingAlgorithm : public NViewMatchingAlgorithm<ThreeViewMatchingControl<TransverseOverlapResult> >
{
public:
    typedef NViewMatchingAlgorithm<ThreeViewMatchingControl<TransverseOverlapResult> > BaseAlgorithm;

    /**
     *  @brief  Default constructor
     */
    ThreeViewDeltaRayMatchingAlgorithm();

    void SelectInputClusters(const pandora::ClusterList *const pInputClusterList, pandora::ClusterList &selectedClusterList) const;

private:
    void CalculateOverlapResult(const pandora::Cluster *const pClusterU, const pandora::Cluster *const pClusterV, const pandora::Cluster *const pClusterW);

    /**
     *  @brief  Calculate the overlap result for given group of clusters
     *
     *  @param  pClusterU the cluster from the U view
     *  @param  pClusterV the cluster from the V view
     *  @param  pClusterW the cluster from the W view
     *  @param  overlapResult to receive the overlap result
     */
    pandora::StatusCode CalculateOverlapResult(const pandora::Cluster *const pClusterU, const pandora::Cluster *const pClusterV, const pandora::Cluster *const pClusterW,
        TransverseOverlapResult &overlapResult) const;

    void ExamineOverlapContainer();
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    typedef std::vector<DeltaRayTensorTool*> TensorToolVector;
    TensorToolVector                  m_algorithmToolVector;          ///< The algorithm tool vector
    
    unsigned int                      m_nMaxTensorToolRepeats;
    unsigned int                      m_minClusterCaloHits;
    float                             m_pseudoChi2Cut;              ///< Pseudo chi2 cut for three view matching 
    float                             m_xOverlapWindow;             ///< The maximum allowed displacement in x position
    float                             m_minMatchedFraction;
    unsigned int                      m_minMatchedPoints;
};

//------------------------------------------------------------------------------------------------------------------------------------------
    
/**
 *  @brief  DeltaRayTensorTool class
 */
class DeltaRayTensorTool : public pandora::AlgorithmTool
{
public:
    typedef ThreeViewDeltaRayMatchingAlgorithm::MatchingType::TensorType TensorType;
    typedef std::vector<TensorType::ElementList::const_iterator> IteratorList;

    /**
     *  @brief  Run the algorithm tool
     *
     *  @param  pAlgorithm address of the calling algorithm
     *  @param  overlapTensor the overlap tensor
     *
     *  @return whether changes have been made by the tool
     */
    virtual bool Run(ThreeViewDeltaRayMatchingAlgorithm *const pAlgorithm, TensorType &overlapTensor) = 0;
};
    
} // namespace lar_content

#endif // #ifndef LAR_THREE_VIEW_SHOWERS_ALGORITHM_H
