/**
 *  @file   larpandoracontent/LArReclustering/ThreeDReclusteringAlgorithm.h
 *
 *  @brief  Header file for the reclustering algorithm class.
 *
 *  $Log: $
 */

#ifndef LAR_THREE_D_RECLUSTERING_ALGORITHM_H
#define LAR_THREE_D_RECLUSTERING_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

namespace lar_content
{

 class ClusteringTool; 
 
//------------------------------------------------------------------------------------------------------------------------------------------

 /**
  *  @brief  RecursivePfoMopUpAlgorithm class
  */
 class ThreeDReclusteringAlgorithm : public pandora::Algorithm
 {
 public:
     /**
      *  @brief  Default constructor
      */
    ThreeDReclusteringAlgorithm();

   /**
    *  @brief  Destructor
    */
    ~ThreeDReclusteringAlgorithm();

private:

    typedef std::vector<ClusteringTool *> ClusteringToolVector;

    pandora::StatusCode Run();


    /**
     *  @brief Create new TwoD clusters and Pfos for each  new ThreeD cluster in newClustersList
     *
     *  @param pPfoToRebuild the address of the original pfo to rebuild 
     *  @param newClustersList a reference to the new list of clusters obtained via the reclustering process 
     */
    pandora::StatusCode RebuildPfo(const pandora::Pfo *pPfoToRebuild, pandora::ClusterList &newClustersList);

     /**
     *  @brief Create new TwoD clusters for each  new ThreeD cluster in newClustersList 
     *
     *  @param pPfoToRebuild the address of the original pfo to rebuild 
     *  @param newClustersList a reference to the new list of clusters obtained via the reclustering process 
     */
    pandora::StatusCode BuildNewTwoDClusters(const pandora::Pfo *pPfoToRebuild, pandora::ClusterList &newClustersList);
    /**
     *  @brief Create new Pfos for each  new ThreeD cluster in newClustersList  
     *
     *  @param pPfoToRebuild the address of the original pfo to rebuild 
     *  @param newClustersList a reference to the new list of clusters obtained via the reclustering process 
     */
    pandora::StatusCode BuildNewPfos(const pandora::Pfo *pPfoToRebuild, pandora::ClusterList &newClustersList);

    /**
     *  @brief Loop over all specified figure of merit names, calculate figures of merit for the CaloHitList under consideration, and return the smallest FOM
     *
     *  @param mergedClusterCaloHitList3D the CaloHitList under consideration 
     *
     *  @return The figure of merit
     */
    float GetFigureOfMerit(pandora::CaloHitList mergedClusterCaloHitList3D);
    
    /**
     *  @brief Calculate the specified figure of merit for the CaloHitList under consideration, under consideration, and return the smallest FOM
     *
     *  @param figureOfMeritName the name of the figure of merit
     *  @param mergedClusterCaloHitList3D the CaloHitList under consideration 
     *
     *  @return The figure of merit
     */
    float GetFigureOfMerit(std::string figureOfMeritName, pandora::CaloHitList mergedClusterCaloHitList3D);

    /**
     *  @brief Loop over all specified figure of merit names, calculate figures of merit for each CaloHitList in the provided vector, under consideration, and return the smallest FOM
     *
     *  @param newClustersCaloHitLists3D the vector of CaloHitLists under consideration 
     *
     *  @return The figure of merit
     */
    float GetFigureOfMerit(std::vector<pandora::CaloHitList*> newClustersCaloHitList3D);

    /**
     *  @brief Calculate the specified figure of merit for each CaloHitList in the provided vector, under consideration, and return the smallest FOM
     *
     *  @param figureOfMeritName the name of the figure of merit
     *  @param newClustersCaloHitLists3D the vector of CaloHitLists under consideration 
     *
     *  @return The figure of merit 
     */
   float GetFigureOfMerit(std::string figureOfMeritName, std::vector<pandora::CaloHitList*> newClustersCaloHitLists3D);

    /** 
     *  @brief Get cheated FOM as an impurity: the fraction of hits that are NOT contributed by the main MC particle. If clustering was perfect, cheated FOM would always be 0.
     *
     * @param List of calo hits for the pfo to recluster
     *
     * @return The figure of merit (purity)
     */
    float GetCheatedFigureOfMerit(pandora::CaloHitList mergedClusterCaloHitList3D);
    
    /**
     *  @brief Select pfos to be reclustered if it passes reclustering criteria
     *
     *  @param pPfo the address of the pfo 
     *
     *  @return bool corresponding to decision: recluster (1), do not recluster (0)
     */
    bool PassesCutsForReclustering(const pandora::ParticleFlowObject *const pPfo);

    /**
     *  @brief Get the MC index for the true particle that contributes most energy to the calo hit 
     *
     *  @param pCaloHit address of the calo hit
     *
     *  @return MCParticle Index
     */
    int GetMainMcParticleIndex(const pandora::CaloHit *const pCaloHit);

    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    ClusteringToolVector m_algorithmToolVector; ///< The reclustering algorithm tool vector
    std::string m_pfoListName; ///< The input pfo list name
    pandora::StringVector m_figureOfMeritNames; ///< The names of the figures of merit to use
    std::string m_newPfosListNameAllAfterReclustering; ///< Name of the list that contains new Pfos after reclustering
    int m_hitThresholdForNewPfo; ///< Minimum nr. of hits to form new 3Dclusters
    std::string m_mcParticleListName; ///< The mc particle list name 
    float m_fomThresholdForReclustering; ///< A threshold on the minimum figure of merit for reclustering
    std::map<int,const pandora::Cluster*> m_newClustersUMap, m_newClustersVMap,m_newClustersWMap; ///< Per-view maps associating new 3D clusters with new 2D clusters 

};


//------------------------------------------------------------------------------------------------------------------------------------------

/**
 *  @brief  ClusteringTool class
 */
class ClusteringTool: public pandora::AlgorithmTool {

public:
    ClusteringTool() = default;
    virtual ~ClusteringTool() = default;
    virtual bool Run(const pandora::Algorithm *const pAlgorithm, std::vector<pandora::CaloHitList*> &newCaloHitListsVector) = 0;
};

} // namespace lar_content

#endif // #ifndef LAR_THREE_D_RECLUSTERING_ALGORITHM_H
