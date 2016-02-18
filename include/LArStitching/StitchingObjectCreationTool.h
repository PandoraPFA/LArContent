/**
 *  @file   LArContent/include/LArStitching/StitchingObjectCreationTool.h
 * 
 *  @brief  Header file for the stitching object creation tool class.
 * 
 *  $Log: $
 */
#ifndef LAR_STITCHING_OBJECT_CREATION_TOOL_H
#define LAR_STITCHING_OBJECT_CREATION_TOOL_H 1

#include "LArStitching/MultiPandoraApi.h"
#include "LArStitching/StitchingAlgorithm.h"

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
    void Recreate3DContent(const pandora::Algorithm *const pAlgorithm, const pandora::Pandora *const pPandora, const VolumeInfo &volumeInfo,
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
    void Recreate3DContent(const pandora::Algorithm *const pAlgorithm, const pandora::ParticleFlowObject *const pInputPfo,
        const pandora::ParticleFlowObject *const pNewParentPfo, const pandora::Pandora *const pPandora, const VolumeInfo &volumeInfo, StitchingInfo &stitchingInfo) const;

    /**
     *  @brief  Create a new calo hit in the current pandora instance, based upon the provided input calo hit
     * 
     *  @param  pAlgorithm the address of the calling algorithm
     *  @param  pInputCaloHit the address of the input calo hit
     *  @param  pInputPfo the address of the associated input pfo
     *  @param  volumeInfo the volume information for the input pandora instance
     * 
     *  @return the address of the new calo hit
     */
    const pandora::CaloHit *CreateCaloHit(const pandora::Algorithm *const pAlgorithm, const pandora::CaloHit *const pInputCaloHit,
        const pandora::ParticleFlowObject *const pInputPfo, const VolumeInfo &volumeInfo) const;

    /**
     *  @brief  Create a new cluster in the current pandora instance, based upon the provided input cluster
     * 
     *  @param  pAlgorithm the address of the calling algorithm
     *  @param  pInputCluster the address of the input cluster
     *  @param  newCaloHitList the list of calo hits for the new cluster
     * 
     *  @return the address of the new cluster
     */
    const pandora::Cluster *CreateCluster(const pandora::Algorithm *const pAlgorithm, const pandora::Cluster *const pInputCluster,
        const pandora::CaloHitList &newCaloHitList) const;

    /**
     *  @brief  Create a new vertex in the current pandora instance, based upon the provided input vertex
     * 
     *  @param  pAlgorithm the address of the calling algorithm
     *  @param  pInputVertex the address of the input vertex
     *  @param  pInputPfo the address of the associated input pfo
     *  @param  volumeInfo the volume information for the input pandora instance
     * 
     *  @return the address of the new vertex
     */
    const pandora::Vertex *CreateVertex(const pandora::Algorithm *const pAlgorithm, const pandora::Vertex *const pInputVertex,
        const pandora::ParticleFlowObject *const pInputPfo, const VolumeInfo &volumeInfo) const;

    /**
     *  @brief  Create a new pfo in the current pandora instance, based upon the provided input pfo
     * 
     *  @param  pAlgorithm the address of the calling algorithm
     *  @param  pInputPfo the address of the input pfo
     *  @param  newClusterList the list of clusters for the new pfo
     *  @param  newVertexList the list of vertices for the new pfo
     * 
     *  @return the address of the new pfo
     */
    const pandora::ParticleFlowObject *CreatePfo(const pandora::Algorithm *const pAlgorithm, const pandora::ParticleFlowObject *const pInputPfo,
        const pandora::ClusterList &newClusterList, const pandora::VertexList &newVertexList) const;

    /**
     *  @brief  Add details to the stitching information block
     * 
     *  @param  pNewPfo the address of a new pfo, recreated from an input pfo
     *  @param  pPandora the address of the pandora instance that owns the input pfo
     *  @param  volumeInfo the volume information for the input pandora instance
     *  @param  stitchingInfo to receive any modifications to the stitching info
     */
    void AddStitchingInfo(const pandora::ParticleFlowObject *const pNewPfo, const pandora::Pandora *const pPandora, const VolumeInfo &volumeInfo,
        StitchingInfo &stitchingInfo) const;

    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    std::string     m_newClusterListName;           ///< The new cluster list name
    std::string     m_newVertexListName;            ///< The new vertex list name
    std::string     m_newPfoListName;               ///< The new pfo list name
};

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::AlgorithmTool *StitchingObjectCreationTool::Factory::CreateAlgorithmTool() const
{
    return new StitchingObjectCreationTool();
}

} // namespace lar_content

#endif // #ifndef LAR_STITCHING_OBJECT_CREATION_TOOL_H
