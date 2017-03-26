/**
 *  @file   larpandoracontent/LArStitching/StitchingObjectCreationTool.h
 *
 *  @brief  Header file for the stitching object creation tool class.
 *
 *  $Log: $
 */
#ifndef LAR_STITCHING_OBJECT_CREATION_TOOL_H
#define LAR_STITCHING_OBJECT_CREATION_TOOL_H 1

#include "larpandoracontent/LArStitching/MultiPandoraApi.h"
#include "larpandoracontent/LArStitching/StitchingAlgorithm.h"

namespace lar_content
{

/**
 *  @brief  StitchingObjectCreationTool class
 */
class StitchingObjectCreationTool : public StitchingTool
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

    /**
     *  @brief  Default constructor
     */
    StitchingObjectCreationTool();

    typedef StitchingAlgorithm::StitchingInfo StitchingInfo;

    void Run(const StitchingAlgorithm *const pAlgorithm, StitchingInfo &stitchingInfo);

private:
    /**
     *  @brief  Recreate all (3D aspects of) pfos associated with a provided pandora instance
     *
     *  @param  pAlgorithm the address of the calling algorithm
     *  @param  pPandora the address of the input pandora instance acting as a source for pfos
     *  @param  volumeInfo the volume information for the input pandora instance
     *  @param  stitchingInfo to receive any modifications to the stitching info
     */
    void Recreate3DContent(const StitchingAlgorithm *const pAlgorithm, const pandora::Pandora *const pPandora, const VolumeInfo &volumeInfo,
        StitchingInfo &stitchingInfo) const;

    /**
     *  @brief  Recreate (3D aspects of) a provided pfo
     *
     *  @param  pAlgorithm the address of the calling algorithm
     *  @param  pInputPfo the address of the input pfo
     *  @param  pNewParentPfo the address of the parent for the new pfo, nullptr if not relevant
     *  @param  pPandora the address of the pandora instance that owns the input pfo
     *  @param  volumeInfo the volume information for the input pandora instance
     *  @param  stitchingInfo to receive any modifications to the stitching info
     */
    void Recreate3DContent(const StitchingAlgorithm *const pAlgorithm, const pandora::ParticleFlowObject *const pInputPfo,
        const pandora::ParticleFlowObject *const pNewParentPfo, const pandora::Pandora *const pPandora, const VolumeInfo &volumeInfo,
        StitchingInfo &stitchingInfo) const;

    /**
     *  @brief  Add details to the stitching information block
     *
     *  @param  pNewPfo the address of a new pfo, recreated from an input pfo
     *  @param  pPandora the address of the pandora instance that owns the input pfo
     *  @param  volumeInfo the volume information for the input pandora instance
     *  @param  stitchingInfo to receive any modifications to the stitching info
     */
    void AddStitchingInfo(const pandora::ParticleFlowObject *const pNewPfo, const pandora::Pandora *const pPandora,
        const VolumeInfo &volumeInfo, StitchingInfo &stitchingInfo) const;

    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

};

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::AlgorithmTool *StitchingObjectCreationTool::Factory::CreateAlgorithmTool() const
{
    return new StitchingObjectCreationTool();
}

} // namespace lar_content

#endif // #ifndef LAR_STITCHING_OBJECT_CREATION_TOOL_H
