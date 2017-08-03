/**
 *  @file   larpandoracontent/LArThreeDReco/LArEventBuilding/NeutrinoIdTool.h
 *
 *  @brief  Header file for the neutrino id tool class.
 *
 *  $Log: $
 */
#ifndef LAR_NEUTRINO_ID_TOOL_H
#define LAR_NEUTRINO_ID_TOOL_H 1

#include "larpandoracontent/LArUtility/ParentAlgorithm.h"

namespace lar_content
{

/**
 *  @brief  NeutrinoIdTool class
 */
class NeutrinoIdTool : public NeutrinoIdBaseTool
{
public:
    /**
     *  @brief  Default constructor
     */
    NeutrinoIdTool();

    void FillNeutrinoProperties(const pandora::PfoList *const pPfoList, ParentAlgorithm::SliceProperties &sliceProperties) const;
    void FillCosmicRayProperties(const pandora::PfoList *const pPfoList, ParentAlgorithm::SliceProperties &sliceProperties) const;
    bool GetNeutrinoSliceIndex(const ParentAlgorithm::SliceIndexToPropertiesMap &sliceIndexToPropertiesMap, unsigned int &neutrinoSliceIndex) const;

private:
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);
};

} // namespace lar_content

#endif // #ifndef LAR_NEUTRINO_ID_TOOL_H
