/**
 *  @file   larpandoracontent/LArControlFlow/SlicingAlgorithm.h
 *
 *  @brief  Header file for the master algorithm class.
 *
 *  $Log: $
 */
#ifndef LAR_SLICING_ALGORITHM_H
#define LAR_SLICING_ALGORITHM_H 1

#include "Pandora/Algorithm.h"
#include "Pandora/AlgorithmTool.h"

namespace lar_content
{

class EventSlicingBaseTool;
class SliceMonitoringTool;

//------------------------------------------------------------------------------------------------------------------------------------------

/**
 *  @brief  SlicingAlgorithm class
 */
class SlicingAlgorithm : public pandora::Algorithm
{
public:
    /**
     *  @brief  Slice class
     */
    class Slice
    {
    public:
        pandora::CaloHitList m_caloHitListU; ///< The u calo hit list
        pandora::CaloHitList m_caloHitListV; ///< The v calo hit list
        pandora::CaloHitList m_caloHitListW; ///< The w calo hit list
    };

    typedef std::vector<Slice> SliceList;
    typedef std::map<pandora::HitType, std::string> HitTypeToNameMap;

    /**
     *  @brief  Default constructor
     */
    SlicingAlgorithm();

private:
    pandora::StatusCode Run();
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    EventSlicingBaseTool *m_pEventSlicingTool;  ///< The address of the event slicing tool
    std::string m_slicingListDeletionAlgorithm; ///< The name of the slicing list deletion algorithm

    HitTypeToNameMap m_caloHitListNames; ///< The hit type to calo hit list name map
    HitTypeToNameMap m_clusterListNames; ///< The hit type to cluster list name map

    std::string m_sliceClusterListName; ///< The name of the output slice cluster list
    std::string m_slicePfoListName;     ///< The name of the output slice pfo list

    SliceMonitoringTool *m_pSliceMonitoringTool; ///< The address of the slice monitoring tool
};

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

/**
 *  @brief  EventSlicingBaseTool class
 */
class EventSlicingBaseTool : public pandora::AlgorithmTool
{
public:
    /**
     *  @brief  Run the slicing tool
     *
     *  @param  pAlgorithm address of the calling algorithm
     *  @param  caloHitListNames the hit type to calo hit list name map
     *  @param  clusterListNames the hit type to cluster list name map
     *  @param  sliceList to receive the populated slice list
     */
    virtual void RunSlicing(const pandora::Algorithm *const pAlgorithm, const SlicingAlgorithm::HitTypeToNameMap &caloHitListNames,
        const SlicingAlgorithm::HitTypeToNameMap &clusterListNames, SlicingAlgorithm::SliceList &sliceList) = 0;
};

} // namespace lar_content

#endif // #ifndef LAR_SLICING_ALGORITHM_H
