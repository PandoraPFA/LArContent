/**
 *  @file   larpandoracontent/LArStitching/StitchingAlgorithm.h
 *
 *  @brief  Header file for the Stitching algorithm class.
 *
 *  $Log: $
 */
#ifndef PANDORA_STITCHING_ALGORITHM_H
#define PANDORA_STITCHING_ALGORITHM_H 1

#include "larpandoracontent/LArThreeDReco/LArPfoMopUp/PfoMopUpBaseAlgorithm.h"

#include "larpandoracontent/LArStitching/MultiPandoraApi.h"

#include <unordered_map>

namespace lar_content
{

class StitchingTool;

//------------------------------------------------------------------------------------------------------------------------------------------

/**
 *  @brief  StitchingAlgorithm class
 */
class StitchingAlgorithm : public PfoMopUpBaseAlgorithm
{
public:
    /**
     *  @brief  Factory class for instantiating algorithm
     */
    class Factory : public pandora::AlgorithmFactory
    {
    public:
        pandora::Algorithm *CreateAlgorithm() const;
    };

    typedef std::unordered_map<const pandora::ParticleFlowObject*, int> PfoToVolumeIdMap;

    /**
     *  @brief  StitchingInfo class
     */
    class StitchingInfo
    {
    public:
        PfoToVolumeIdMap    m_pfoToVolumeIdMap;         ///< Mapping between Pfos and Volume IDs
    };

    /**
     *  @brief  Get the new/recreated cluster list name
     *
     *  @return the new/recreated cluster list name
     */
    const std::string &GetNewClusterListName() const;

    /**
     *  @brief  Get the new/recreated vertex list name
     *
     *  @return the new/recreated vertex list name
     */
    const std::string &GetNewVertexListName() const;

    /**
     *  @brief  Get the new/recreated pfo list name
     *
     *  @return the new/recreated pfo list name
     */
    const std::string &GetNewPfoListName() const;

    /**
     *  @brief  Create a new calo hit in the current pandora instance, based upon the provided input calo hit
     *
     *  @param  pInputCaloHit the address of the input calo hit
     *  @param  volumeInfo the volume information block for this drift volume
     *  @param  x0 the x0 correction relative to the input hit
     *
     *  @return the address of the new calo hit
     */
    const pandora::CaloHit *CreateCaloHit(const pandora::CaloHit *const pInputCaloHit, const VolumeInfo &volumeInfo, const float x0) const;

    /**
     *  @brief  Create a new cluster in the current pandora instance, based upon the provided input cluster
     *
     *  @param  pInputCluster the address of the input cluster
     *  @param  newCaloHitList the list of calo hits for the new cluster
     *
     *  @return the address of the new cluster
     */
    const pandora::Cluster *CreateCluster(const pandora::Cluster *const pInputCluster, const pandora::CaloHitList &newCaloHitList) const;

    /**
     *  @brief  Create a new vertex in the current pandora instance, based upon the provided input vertex
     *
     *  @param  pInputVertex the address of the input vertex
     *  @param  volumeInfo the volume information block for this drift volume
     *  @param  x0 the x0 correction relative to the input vertex
     *
     *  @return the address of the new vertex
     */
    const pandora::Vertex *CreateVertex(const pandora::Vertex *const pInputVertex, const VolumeInfo &volumeInfo, const float x0) const;

    /**
     *  @brief  Create a new pfo in the current pandora instance, based upon the provided input pfo
     *
     *  @param  pInputPfo the address of the input pfo
     *  @param  newClusterList the list of clusters for the new pfo
     *  @param  newVertexList the list of vertices for the new pfo
     *
     *  @return the address of the new pfo
     */
    const pandora::ParticleFlowObject *CreatePfo(const pandora::ParticleFlowObject *const pInputPfo,
        const pandora::ClusterList &newClusterList, const pandora::VertexList &newVertexList) const;

    /**
     *  @brief  Shift a Pfo hierarchy by a specified x0 value
     *
     *  @param  pPfo  the address of the parent pfo
     *  @param  stitchingInfo  the source for additional, local, stitching information
     *  @param  x0  the x0 correction relative to the input pfo
     */
    void ShiftPfoHierarchy(const pandora::ParticleFlowObject *const pPfo, const StitchingAlgorithm::StitchingInfo &stitchingInfo,
        const float x0) const;

    /**
     *  @brief  Shift a Pfo by a specified x0 value
     *
     *  @param  pInputPfo  the address of the input pfo
     *  @param  volumeInfo  the volume information block for this drift volume
     *  @param  x0  the x0 correction relative to the input pfo
     */
    void ShiftPfo(const pandora::ParticleFlowObject *const pInputPfo, const VolumeInfo &volumeInfo, const float x0) const;

    /**
     *  @brief  Stitch together a pair of pfos
     *
     *  @param  pPfoToEnlarge the address of the pfo to enlarge
     *  @param  pPfoToDelete the address of the pfo to delete (will become a dangling pointer)
     *  @param  stitchingInfo the source for additional, local, stitching information
     */
    void StitchPfos(const pandora::ParticleFlowObject *const pPfoToEnlarge, const pandora::ParticleFlowObject *const pPfoToDelete,
        StitchingAlgorithm::StitchingInfo &stitchingInfo) const;

private:
    pandora::StatusCode Run();
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    typedef std::vector<StitchingTool*> StitchingToolVector;
    StitchingToolVector     m_algorithmToolVector;      ///< The algorithm tool vector

    std::string             m_newClusterListName;       ///< The new/recreated cluster list name
    std::string             m_newVertexListName;        ///< The new/recreated vertex list name
    std::string             m_newPfoListName;           ///< The new/recreated pfo list name
};

//------------------------------------------------------------------------------------------------------------------------------------------

/**
 *  @brief  StitchingTool class
 */
class StitchingTool : public pandora::AlgorithmTool
{
public:
    /**
     *  @brief  Run the algorithm tool
     *
     *  @param  pAlgorithm address of the calling algorithm
     *  @param  stitchingInfo the source for additional, local, stitching information
     */
    virtual void Run(const StitchingAlgorithm *const pAlgorithm, StitchingAlgorithm::StitchingInfo &stitchingInfo) = 0;
};

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::Algorithm *StitchingAlgorithm::Factory::CreateAlgorithm() const
{
    return new StitchingAlgorithm();
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

inline const std::string &StitchingAlgorithm::GetNewClusterListName() const
{
    return m_newClusterListName;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline const std::string &StitchingAlgorithm::GetNewVertexListName() const
{
    return m_newVertexListName;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline const std::string &StitchingAlgorithm::GetNewPfoListName() const
{
    return m_newPfoListName;
}

} // namespace lar_content

#endif // #ifndef PANDORA_STITCHING_ALGORITHM_H
