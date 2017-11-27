/**
 *  @file   larpandoracontent/LArControlFlow/BeamParticleIdTool.cc
 *
 *  @brief  Implementation of the beam particle id tool class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArControlFlow/BeamParticleIdTool.h"

using namespace pandora;

namespace lar_content
{

BeamParticleIdTool::BeamParticleIdTool() :
    m_selectAllBeamParticles(true),
    m_selectOnlyFirstSliceBeamParticles(false)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

void BeamParticleIdTool::SelectOutputPfos(const SliceHypotheses &beamSliceHypotheses, const SliceHypotheses &crSliceHypotheses, PfoList &selectedPfos)
{
    if (beamSliceHypotheses.size() != crSliceHypotheses.size())
        throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);

    for (unsigned int sliceIndex = 0, nSlices = beamSliceHypotheses.size(); sliceIndex < nSlices; ++sliceIndex)
    {
        const PfoList &sliceOutput((m_selectAllBeamParticles || (m_selectOnlyFirstSliceBeamParticles && (0 == sliceIndex))) ?
            beamSliceHypotheses.at(sliceIndex) : crSliceHypotheses.at(sliceIndex));

        selectedPfos.insert(selectedPfos.end(), sliceOutput.begin(), sliceOutput.end());
    }

/*
    for (unsigned int sliceIndex = 0, nSlices = beamSliceHypotheses.size(); sliceIndex < nSlices; ++sliceIndex)
    {
        const PfoList pfoListBeam(beamSliceHypotheses.at(sliceIndex));

        PfoList allConnectedPfoList;
        LArPfoHelper::GetAllConnectedPfos(pfoListBeam, allConnectedPfoList);

        CaloHitList caloHitList3D;
        for (const ParticleFlowObject *const pPfo : allConnectedPfoList)
        {
            // Get Clusters
            ClusterList clusterList;
            LArPfoHelper::GetThreeDClusterList(pPfo, clusterList);

            // Get Calo Hit List
            for (const Cluster *const pCluster : clusterList)
            {
                pCluster->GetOrderedCaloHitList().FillCaloHitList(caloHitList3D);
            }
        }

        CartesianVector centroid(0.f, 0.f, 0.f);
        LArPcaHelper::EigenVectors eigenVecs;
        LArPcaHelper::EigenValues eigenValues(0.f, 0.f, 0.f);
        LArPcaHelper::RunPca(caloHitList3D, centroid, eigenValues, eigenVecs);

        // Major axis of slice
        CartesianVector fitDirection(eigenVecs.at(0));







        const PfoList &sliceOutput((m_selectAllBeamParticles || (m_selectOnlyFirstSliceBeamParticles && (0 == sliceIndex))) ?
            beamSliceHypotheses.at(sliceIndex) : crSliceHypotheses.at(sliceIndex));

        selectedPfos.insert(selectedPfos.end(), sliceOutput.begin(), sliceOutput.end());
    }
*/
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode BeamParticleIdTool::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "SelectAllBeamParticles", m_selectAllBeamParticles));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "SelectOnlyFirstSliceBeamParticles", m_selectOnlyFirstSliceBeamParticles));

    if (m_selectAllBeamParticles == m_selectOnlyFirstSliceBeamParticles)
    {
        std::cout << "BeamParticleIdTool::ReadSettings - exactly one of SelectAllBeamParticles and SelectOnlyFirstSliceBeamParticles must be true" << std::endl;
        return STATUS_CODE_INVALID_PARAMETER;
    }

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
