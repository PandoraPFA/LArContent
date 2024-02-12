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

#include "larpandoracontent/LArControlFlow/EventSlicingBaseTool.h"
#include "larpandoracontent/LArObjects/LArSlice.h"

namespace lar_content
{

typedef std::map<pandora::HitType, std::string> HitTypeToNameMap;

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

} // namespace lar_content

#endif // #ifndef LAR_SLICING_ALGORITHM_H
