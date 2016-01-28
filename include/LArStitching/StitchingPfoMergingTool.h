/**
 *  @file   LArContent/include/LArStitching/StitchingPfoMergingTool.h
 * 
 *  @brief  Header file for the stitching pfo merging tool class.
 * 
 *  $Log: $
 */
#ifndef LAR_STITCHING_PFO_MERGING_TOOL_H
#define LAR_STITCHING_PFO_MERGING_TOOL_H 1

#include "LArStitching/MultiPandoraApi.h"
#include "LArStitching/StitchingAlgorithm.h"

namespace lar_content
{

/**
 *  @brief  StitchingPfoMergingTool class
 */
class StitchingPfoMergingTool : public StitchingTool
{
public:
    /**
     *  @brief  Factory class for instantiating algorithm tool
     */
    class Factory : public pandora::AlgorithmToolFactory
    {
    public:
        pandora::AlgorithmTool *CreateAlgorithmTool() const;
    };

    void Run(const StitchingAlgorithm *const pAlgorithm, StitchingAlgorithm::StitchingInfo &stitchingInfo);

private:
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);
};

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::AlgorithmTool *StitchingPfoMergingTool::Factory::CreateAlgorithmTool() const
{
    return new StitchingPfoMergingTool();
}

} // namespace lar_content

#endif // #ifndef LAR_STITCHING_PFO_MERGING_TOOL_H
