/**
 *  @file   larpandoracontent/LArCheating/CheatingNeutrinoIdTool.h
 *
 *  @brief  Header file for the cheating neutrino id tool class.
 *
 *  $Log: $
 */
#ifndef LAR_CHEATING_NEUTRINO_ID_TOOL_H
#define LAR_CHEATING_NEUTRINO_ID_TOOL_H 1

#include "larpandoracontent/LArUtility/MasterAlgorithm.h"

namespace lar_content
{

/**
 *  @brief  CheatingNeutrinoIdTool class
 */
class CheatingNeutrinoIdTool : public NeutrinoIdBaseTool
{
public:
    /**
     *  @brief  Default constructor
     */
    CheatingNeutrinoIdTool();

    void FillNeutrinoProperties(const pandora::PfoList *const pPfoList, SliceProperties &sliceProperties) const;
    void FillCosmicRayProperties(const pandora::PfoList *const pPfoList, SliceProperties &sliceProperties) const;
    bool GetNeutrinoSliceIndex(const SliceIndexToPropertiesMap &sliceIndexToPropertiesMap, unsigned int &neutrinoSliceIndex) const;

private:
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);
};

} // namespace lar_content

#endif // #ifndef LAR_CHEATING_NEUTRINO_ID_TOOL_H
