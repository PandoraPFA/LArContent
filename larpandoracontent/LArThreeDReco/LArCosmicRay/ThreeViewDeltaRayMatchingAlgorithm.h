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

#include "larpandoracontent/LArThreeDReco/LArCosmicRay/NViewDeltaRayMatchingAlgorithm.h"
#include "larpandoracontent/LArThreeDReco/LArThreeDBase/ThreeViewMatchingControl.h"

#include "larpandoracontent/LArUtility/KDTreeLinkerAlgoT.h"

namespace lar_content
{
    
class DeltaRayTensorTool;

//------------------------------------------------------------------------------------------------------------------------------------------

/**
 *  @brief  ThreeViewDeltaRayMatchingAlgorithm class
 */
class ThreeViewDeltaRayMatchingAlgorithm : public NViewDeltaRayMatchingAlgorithm<ThreeViewMatchingControl<DeltaRayOverlapResult> >
{
public:
    typedef NViewDeltaRayMatchingAlgorithm<ThreeViewMatchingControl<DeltaRayOverlapResult> > BaseAlgorithm;
    typedef ThreeViewDeltaRayMatchingAlgorithm::MatchingType::TensorType TensorType;

    /**
     *  @brief  Default constructor
     */
    ThreeViewDeltaRayMatchingAlgorithm();

    void GetConnectedElements(const pandora::Cluster *const pCluster1, const bool hasAssociatedMuon, TensorType::ElementList &elementList, pandora::ClusterSet &checkedClusters);

    void GetUnambiguousElements(const bool hasAssociatedMuon, TensorType::ElementList &elementList);

    bool DoesClusterPassTesorThreshold(const pandora::Cluster *const pCluster) const;    

    std::string GetClusteringAlgName() const;
    
private:
    void CalculateOverlapResult(const pandora::Cluster *const pClusterU, const pandora::Cluster *const pClusterV, const pandora::Cluster *const pClusterW);

    pandora::StatusCode CalculateOverlapResult(const pandora::Cluster *const pClusterU, const pandora::Cluster *const pClusterV, const pandora::Cluster *const pClusterW,
        DeltaRayOverlapResult &overlapResult) const;

    void FindCommonMuonParents(const pandora::Cluster *const pClusterU, const pandora::Cluster *const pClusterV, const pandora::Cluster *const pClusterW, pandora::PfoList &commonMuonPfoList) const;

    void ExamineOverlapContainer();
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    typedef std::vector<DeltaRayTensorTool*> TensorToolVector;
    TensorToolVector                  m_algorithmToolVector;          ///< The algorithm tool vector
    
    unsigned int                      m_nMaxTensorToolRepeats;
    unsigned int                      m_minClusterCaloHits;

    std::string  m_reclusteringAlgorithmName;
};

//------------------------------------------------------------------------------------------------------------------------------------------

inline std::string ThreeViewDeltaRayMatchingAlgorithm::GetClusteringAlgName() const
{
    return m_reclusteringAlgorithmName;
}
    
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

#endif // #ifndef LAR_THREE_VIEW_DELTA_RAY_MATCHING_ALGORITHM_H
