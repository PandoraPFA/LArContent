/**
 *  @file   larpandoracontent/LArReclustering/ThreeDMultiReclusteringAlgorithm.h
 *
 *  @brief  Header file for the reclustering algorithm class.
 *
 *  $Log: $
 */

#ifndef LAR_THREE_D_MULTI_RECLUSTERING_ALGORITHM_H
#define LAR_THREE_D_MULTI_RECLUSTERING_ALGORITHM_H 1

#include "Pandora/Algorithm.h"
#include "Pandora/PandoraInternal.h"

#include "larpandoracontent/LArReclustering/ThreeDReclusteringFigureOfMeritBaseTool.h"

namespace lar_content
{

/**
  *  @brief  ThreeDMultiReclusteringAlgorithm class
  */
class ThreeDMultiReclusteringAlgorithm : public pandora::Algorithm
{
public:

    /**
     *  @brief  Default constructor
     */
    ThreeDMultiReclusteringAlgorithm();

    /**
    *  @brief  Destructor
    */
    ~ThreeDMultiReclusteringAlgorithm();

private:
    float GetRandomFom();

    pandora::StatusCode Run();

    pandora::StatusCode FreeClustersFromPfos(const pandora::PfoList &pfos,
                                             std::map<pandora::HitType, pandora::ClusterList> &viewToFreeClusters,
                                             std::map<const pandora::Pfo *const, pandora::ClusterList> &pfoToFreeClusters);

    pandora::StatusCode FreeCaloHitsFromClusters(const pandora::ClusterList &clusters,
                                                 const pandora::HitType &view,
                                                 pandora::CaloHitList &freeCaloHits,
                                                 pandora::CaloHitList &freeIsoCaloHits);

    pandora::StatusCode GetClusterListName(const pandora::HitType &view, std::string &listName);

    pandora::StatusCode AddClustersToPfos(std::map<const pandora::Pfo *const, pandora::ClusterList> &pfoToClusters);

    pandora::StatusCode BuildNewPfos(const pandora::ClusterList &clusters3D,
                                     std::map<pandora::HitType, pandora::CaloHitList> &viewToFreeCaloHits2D,
                                     std::map<pandora::HitType, pandora::CaloHitList> &viewToFreeIsoCaloHits2D);

    pandora::StatusCode Build2DClustersFrom3D(const pandora::Cluster *const pCluster3D,
                                              std::map<pandora::HitType, pandora::CaloHitList> &viewToFreeCaloHits2D,
                                              std::map<pandora::HitType, pandora::CaloHitList> &viewToFreeIsoCaloHits2D,
                                              pandora::ClusterList &newClusters2D);

    pandora::StatusCode MopUpCaloHits(const pandora::CaloHitList &caloHits, const pandora::ClusterList &clusters, bool addAsIso);

    pandora::StatusCode CreatePfoFromClusters(const pandora::ClusterList &clusters);

    pandora::StatusCode DebugCurrentLists();

    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    std::string                              m_pfoListName;       ///< Name of list of pfos to consider for reclustering
    std::string                              m_cluster3DListName; ///< Name of list of 3D clusters that comprise the pfos
    std::vector<std::string>                 m_clusteringAlgs;    ///< The ordered list of clustering algorithms to use
    ThreeDReclusteringFigureOfMeritBaseTool *m_pFomAlgTool;       ///< The address of the figure of merit algorithm tool to use
    std::string                              m_clusterUListName;  ///
    std::string                              m_clusterVListName;  ///
    std::string                              m_clusterWListName;  ///< Name of list of 2D U,V,W clusters that may need to get fragmented after reclustring of 3D clusters
};

//------------------------------------------------------------------------------------------------------------------------------------------

inline float ThreeDMultiReclusteringAlgorithm::GetRandomFom() { return static_cast<float>(rand()) / RAND_MAX; }

inline pandora::StatusCode ThreeDMultiReclusteringAlgorithm::GetClusterListName(const pandora::HitType &view, std::string &listName)
{ 
    switch (view)
    {
        case pandora::HitType::TPC_3D:
            listName = m_cluster3DListName;
            break;
        case pandora::HitType::TPC_VIEW_U:
            listName = m_clusterUListName;
            break;
        case pandora::HitType::TPC_VIEW_V:
            listName = m_clusterVListName;
            break;
        case pandora::HitType::TPC_VIEW_W:
            listName = m_clusterWListName;
            break;
        default:
            return pandora::StatusCode::STATUS_CODE_FAILURE;
    }
    return pandora::StatusCode::STATUS_CODE_SUCCESS;
}

} // namespace lar_content

#endif // #ifndef LAR_THREE_D_MULTI_RECLUSTERING_ALGORITHM_H
