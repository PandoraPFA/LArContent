/**
 *  @file   larpandoracontent/LArThreeDReco/LArCosmicRay/TwoViewDeltaRayMatching.h
 *
 *  @brief  Header file for the two view delta ray matching class
 *
 *  $Log: $
 */
#ifndef LAR_TWO_VIEW_DELTA_RAY_MATCHING_ALGORITHM_H
#define LAR_TWO_VIEW_DELTA_RAY_MATCHING_ALGORITHM_H 1

#include "Pandora/Algorithm.h"
#include "Pandora/AlgorithmTool.h"

#include "larpandoracontent/LArObjects/LArTrackTwoViewOverlapResult.h"

#include "larpandoracontent/LArThreeDReco/LArCosmicRay/NViewDeltaRayMatchingAlgorithm.h"
#include "larpandoracontent/LArThreeDReco/LArThreeDBase/TwoViewMatchingControl.h"

#include "larpandoracontent/LArUtility/KDTreeLinkerAlgoT.h"

namespace lar_content
{
    
class DeltaRayMatrixTool;

//------------------------------------------------------------------------------------------------------------------------------------------

/**
 *  @brief  TwoViewDeltaRayMatchingAlgorithm class
 */
class TwoViewDeltaRayMatchingAlgorithm : public NViewDeltaRayMatchingAlgorithm<TwoViewMatchingControl<TwoViewDeltaRayOverlapResult> >
{
public:
    typedef NViewDeltaRayMatchingAlgorithm<TwoViewMatchingControl<TwoViewDeltaRayOverlapResult> > BaseAlgorithm;
    typedef TwoViewDeltaRayMatchingAlgorithm::MatchingType::MatrixType MatrixType;
    typedef std::vector<pandora::HitType> HitTypeVector;
    
    /**
     *  @brief  Default constructor
     */
    TwoViewDeltaRayMatchingAlgorithm();

    HitTypeVector GetHitTypeVector();

    const pandora::Cluster *GetCluster(const MatrixType::Element &element, const pandora::HitType &hitType);

    bool DoesClusterPassTesorThreshold(const pandora::Cluster *const pCluster) const;

    void RemoveThirdViewCluster(const pandora::Cluster *const pCluster);


    const std::string &GetThirdViewClusterListName() const;

    bool CreatePfo(const MatrixType::Element &element);

    std::string GetClusteringAlgName() const;
    
private:
    void CalculateOverlapResult(const pandora::Cluster *const pCluster1, const pandora::Cluster *const pCluster2, const pandora::Cluster *const pCluster3);

    pandora::StatusCode CalculateOverlapResult(const pandora::Cluster *const pCluster1, const pandora::Cluster *const pCluster2, const pandora::Cluster *const pCluster3,
        const pandora::CaloHitList &projectedHits, TwoViewDeltaRayOverlapResult &overlapResult) const;    
    
    pandora::StatusCode CalculateOverlapResult(const pandora::Cluster *const pCluster1, const pandora::Cluster *const pCluster2,
        TwoViewDeltaRayOverlapResult &overlapResult) const;

    void FindCommonMuonParents(const pandora::Cluster *const pCluster1, const pandora::Cluster *const pCluster2, pandora::PfoList &commonMuonPfoList) const;
    
    void CollectThirdViewClusters(const pandora::Cluster *const pCluster1, const pandora::Cluster *const pCluster2, const pandora::CartesianPointVector &projectedPositions,
        pandora::ClusterList &matchedClusters) const;

    void GetBestMatchedCluster(const pandora::Cluster *const pCluster1, const pandora::Cluster *const pCluster2, const pandora::PfoList &commonMuonPfoList,
        const pandora::ClusterList &matchedClusters, const pandora::Cluster *&pBestMatchedCluster, float &reducedChiSquared) const;
    
    void ExamineOverlapContainer();
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    void MergeThirdView(const MatrixType::Element &element, const pandora::Cluster *const pSeedCluster);

    void GetBestMatchedAvailableCluster(const pandora::ClusterList &matchedClusters, const pandora::Cluster *&pBestMatchedCluster) const;

    void GrowThirdView(const MatrixType::Element &element, ProtoParticle &protoParticle);

    std::string m_inputClusterListName;    

    typedef std::vector<DeltaRayMatrixTool*> MatrixToolVector;
    MatrixToolVector                  m_algorithmToolVector;          ///< The algorithm tool vector
    
    unsigned int                      m_nMaxMatrixToolRepeats;
    unsigned int                      m_minClusterCaloHits;
    float                             m_maxDistanceFromPrediction;
    float m_maxGoodMatchReducedChiSquared;

    std::string  m_reclusteringAlgorithmName;
};

//------------------------------------------------------------------------------------------------------------------------------------------

inline std::string TwoViewDeltaRayMatchingAlgorithm::GetClusteringAlgName() const
{
    return m_reclusteringAlgorithmName;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline const std::string &TwoViewDeltaRayMatchingAlgorithm::GetThirdViewClusterListName() const
{
    return m_inputClusterListName;
}
    
/**
 *  @brief  DeltaRayTensorTool class
 */
class DeltaRayMatrixTool : public pandora::AlgorithmTool
{
public:
    typedef TwoViewDeltaRayMatchingAlgorithm::MatchingType::MatrixType MatrixType;
    typedef std::vector<MatrixType::ElementList::const_iterator> IteratorList;

    /**
     *  @brief  Run the algorithm tool
     *
     *  @param  pAlgorithm address of the calling algorithm
     *  @param  overlapTensor the overlap tensor
     *
     *  @return whether changes have been made by the tool
     */
    virtual bool Run(TwoViewDeltaRayMatchingAlgorithm *const pAlgorithm, MatrixType &matrixTensor) = 0;
    };
    
} // namespace lar_content

#endif // #ifndef LAR_TWO_VIEW_DELTA_RAY_MATCHING_ALGORITHM_H
