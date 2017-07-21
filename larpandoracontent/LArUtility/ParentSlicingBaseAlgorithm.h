/**
 *  @file   larpandoracontent/LArUtility/ParentSlicingBaseAlgorithm.h
 *
 *  @brief  Header file for the parent slicing base algorithm class.
 *
 *  $Log: $
 */
#ifndef LAR_PARENT_SLICING_BASE_ALGORITHM_H
#define LAR_PARENT_SLICING_BASE_ALGORITHM_H 1

#include "larpandoracontent/LArUtility/ParentBaseAlgorithm.h"

namespace lar_content
{

class SlicingTool;

//------------------------------------------------------------------------------------------------------------------------------------------

/**
 *  @brief  ParentSlicingBaseAlgorithm class
 */
class ParentSlicingBaseAlgorithm : public ParentBaseAlgorithm
{
public:
    /**
     *  @brief  Default constructor
     */
    ParentSlicingBaseAlgorithm();

    /**
     *  @brief  Destructor
     */
    virtual ~ParentSlicingBaseAlgorithm();

    /**
     *  @brief  Slice class
     */
    class Slice
    {
    public:
        pandora::CaloHitList    m_caloHitListU;                     ///< The u calo hit list
        pandora::CaloHitList    m_caloHitListV;                     ///< The v calo hit list
        pandora::CaloHitList    m_caloHitListW;                     ///< The w calo hit list
    };

    typedef std::vector<Slice> SliceList;

    /**
     *  @brief  Copy all the input hits in an event into a single slice
     *
     *  @param  sliceList the slice list to receive the single new slice
     */
    void CopyAllHitsToSingleSlice(SliceList &sliceList) const;

protected:
    /**
     *  @brief  Perform first-pass 3D event reconstruction
     */
    virtual void FastReconstruction() const = 0;

    /**
     *  @brief  Use first-pass 3D event reconstruction to slice events into separate, distinct interactions for processing
     *
     *  @param  sliceList the slice list to receive the slice list
     */
    virtual void PerformSlicing(SliceList &sliceList) const;

    /**
     *  @brief  Run two dimensional clustering, for a given slice identifier, using hit list names provided via algorithm config
     *
     *  @param  sliceIndexString the slice index string/identifier
     *  @param  clusteringAlgName the clustering algorithm name
     *  @param  existingClusterList whether the intent is to add clusters to an existing output list, or fill this list for first time
     *  @param  additionalTwoDAlgorithms the names of any additional two dimensional algorithms to process each new cluster list
     */
    void RunTwoDClustering(const std::string &sliceIndexString, const std::string &clusteringAlgName, const bool existingClusterList,
        const pandora::StringVector &additionalTwoDAlgorithms = pandora::StringVector()) const;

    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    bool                        m_shouldPerformSlicing;             ///< Whether to slice events into separate, distinct interactions for processing
    SlicingTool                *m_pSlicingTool;                     ///< The address of the slicing tool
    std::string                 m_slicingListDeletionAlgorithm;     ///< The name of the slicing list deletion algorithm
};

//------------------------------------------------------------------------------------------------------------------------------------------

/**
 *  @brief  SlicingTool class
 */
class SlicingTool : public pandora::AlgorithmTool
{
public:
    /**
     *  @brief  Run the algorithm tool
     *
     *  @param  pAlgorithm address of the calling algorithm
     *  @param  caloHitListNames the hit type to calo hit list name map
     *  @param  clusterListNames the hit type to cluster list name map
     *  @param  sliceList to receive the populated slice list
     */
    virtual void Slice(const ParentSlicingBaseAlgorithm *const pAlgorithm, const ParentSlicingBaseAlgorithm::HitTypeToNameMap &caloHitListNames,
        const ParentSlicingBaseAlgorithm::HitTypeToNameMap &clusterListNames, ParentSlicingBaseAlgorithm::SliceList &sliceList) = 0;
};

} // namespace lar_content

#endif // #ifndef LAR_PARENT_SLICING_BASE_ALGORITHM_H
