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

    pandora::StatusCode Run();

    /**
     *  @brief Remove clusters from the pfos and store them them in maps
     *
     *  @param pfos pfo list
     *  @param viewToFreeClusters output map of hit type to list of removed clusters
     *  @param pfoToFreeClusters output map of original pfo to list of removed clusters
     */
    pandora::StatusCode FreeClustersFromPfos(const pandora::PfoList &pfos,
                                             std::map<pandora::HitType, pandora::ClusterList> &viewToFreeClusters,
                                             std::map<const pandora::Pfo *const, pandora::ClusterList> &pfoToFreeClusters);
    /**
     *  @brief Add clusters to existing pfos according to a map
     *
     *  @param pfoToClusters map from pfo to a list of clusters
     */
    pandora::StatusCode AddClustersToPfos(std::map<const pandora::Pfo *const, pandora::ClusterList> &pfoToClusters);

    /**
     *  @brief Delete clusters of a single hit type and store the associated hits in a list
     *
     *  @param clusters list of clusters to be deleted
     *  @param view hit type of the clusters
     *  @param freeCaloHits output hits of the deleted clusters
     */
    pandora::StatusCode FreeCaloHitsFromClusters(const pandora::ClusterList &clusters,
                                                 const pandora::HitType &view,
                                                 pandora::CaloHitList &freeCaloHits);
    /**
     *  @brief Create new pfos from the reclustered 3D clusters and original 2D hits. The original 2D hits are put into new 2D clusters
     *         that follow the new 3D clusters,
     *
     *  @param clusters3D list of reclustered 3D clusters
     *  @param viewToFreeCaloHits2D map of hit type to original 2D hits that need to be clustered
     */
    pandora::StatusCode BuildNewPfos(const pandora::ClusterList &clusters3D,
                                     std::map<pandora::HitType, pandora::CaloHitList> &viewToFreeCaloHits2D);

    /**
     *  @brief Create 2D clusters following a 3D cluster
     *
     *  @param pCluster3D a 3D cluster
     *  @param viewToFreeCaloHits2D map of hit type to original 2D hits that need to be clustered
     *  @param newClusters2D output list of newly created 2D clusters
     */
    pandora::StatusCode Build2DClustersFrom3D(const pandora::Cluster *const pCluster3D,
                                              std::map<pandora::HitType, pandora::CaloHitList> &viewToFreeCaloHits2D,
                                              pandora::ClusterList &newClusters2D);

    /**
     *  @brief Add hits to their nearest cluster
     *
     *  @param caloHits list of 2D hits of the same hit type to be added to nearest cluster
     *  @param clusters list of 2D clusters of the same hit type as the caloHits to have hits added to them
     *  @param addAsIso bool for adding the hits as isolated or not
     */
    pandora::StatusCode MopUpCaloHits(const pandora::CaloHitList &caloHits, const pandora::ClusterList &clusters, bool addAsIso);

    /**
     *  @brief Create a new pfo from a list of clusters
     *
     *  @param clusters list of clusters, expect a 3D and up to 3 2D clusters.
     */
    pandora::StatusCode CreatePfoFromClusters(const pandora::ClusterList &clusters);

    /**
     *  @brief Get the cluster list name associated with a hit type
     *
     *  @param view hit type
     *  @param listName output list name
     */
    pandora::StatusCode GetClusterListName(const pandora::HitType &view, std::string &listName);

    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    std::string                              m_pfoListName;       ///< Name of list of pfos to consider for reclustering
    std::string                              m_cluster3DListName; ///< Name of list of 3D clusters that comprise the pfos
    std::vector<std::string>                 m_clusteringAlgs;    ///< The ordered list of clustering algorithms to use
    ThreeDReclusteringFigureOfMeritBaseTool *m_pFomAlgTool;       ///< The address of the figure of merit algorithm tool to use
    std::string                              m_clusterUListName;  ///< Name of list of 2D U clusters that may need reclustering according to new 3D clusters
    std::string                              m_clusterVListName;  ///< Name of list of 2D V clusters that may need reclustering according to new 3D clusters
    std::string                              m_clusterWListName;  ///< Name of list of 2D W clusters that may need reclustering according to new 3D clusters
};

//------------------------------------------------------------------------------------------------------------------------------------------

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
