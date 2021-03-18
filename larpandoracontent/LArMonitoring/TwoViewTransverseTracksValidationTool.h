/**
 *  @file   larpandoracontent/LArMonitoring/TwoViewTransverseTracksValidationTool.h
 *
 *  @brief  Header file for the two view transverse tracks validation tool
 *
 *  $Log: $
 */
#ifndef TWO_VIEW_TRANSVERSE_TRACKS_VALIDATION_TOOL_H
#define TWO_VIEW_TRANSVERSE_TRACKS_VALIDATION_TOOL_H 1

#include "Pandora/AlgorithmTool.h"
namespace lar_content
{

/**
 *  @brief  TwoViewTransverseTracksValidationTool class
 */
class TwoViewTransverseTracksValidationTool : public pandora::AlgorithmTool
{
public:
    /**
     *  @brief  Default constructor
     */
    TwoViewTransverseTracksValidationTool();

    bool Run(/*to be filled in*/);

private:
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    //int m_DataMember;             ///<data memeber
};

} // namespace lar_content

#endif // #ifndef TWO_VIEW_TRANSVERSE_TRACKS_VALIDATION_TOOL_H
