/**
 *  @file   larpandoracontent/LArMonitoring/SliceIdMonitoringTool.h
 *
 *  @brief  Header file for the cosmic-ray tagging monitoring tool class.
 *
 *  $Log: $
 */

#ifndef LAR_SLICE_ID_MONITORING_TOOL_H
#define LAR_SLICE_ID_MONITORING_TOOL_H 1

#include "larpandoracontent/LArControlFlow/MasterAlgorithm.h"

namespace lar_content
{

class SliceIdMonitoringTool : public SliceIdBaseTool
{
public:
    /**
     *  @brief  Default constructor
     */
    SliceIdMonitoringTool();
    virtual ~SliceIdMonitoringTool();

    void SelectOutputPfos(const pandora::Algorithm *const pAlgorithm, const SliceHypotheses &nuSliceHypotheses, const SliceHypotheses &crSliceHypotheses, pandora::PfoList &selectedPfos, const PfoToFloatMap &pfotoprobabilitymapb, const SliceVector &sliceVector);

private:
    /**
     *  @brief  Count the number of neutrino induced hits in a given list using MC information
     *
     *  @param  caloHitSet input list of calo hits
     *
     *  @return the number of neutrino induced hits in the input list
     */
    pandora::CaloHitList CountNeutrinoInducedHits(const pandora::CaloHitList &caloHitList) const;
    

    /**
     *  @brief  Read settings
     */
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    float                                       m_minPurity;                  ///< The minimum purity to consider a Pfo as "pure"
     float                                       m_minSignificance;            ///< The minimum significance to consider a Pfo as "significant"

};



}

#endif



